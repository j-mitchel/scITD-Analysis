library(Seurat)
library(SeuratDisk)
library(sccore)

pbmc <- LoadH5Seurat("/home/jmitchel/data/cite_seq_data/pbmc_multimodal.h5seurat")

# looking at the meta data
head(pbmc@meta.data)
Idents(pbmc)[1:10]
length(unique(pbmc@meta.data$celltype.l3)) # 58 cell types at finest resolution
length(unique(pbmc@meta.data$donor)) # 8 donors

dim(pbmc@assays$SCT@counts) # gene x cells
dim(pbmc@assays$ADT@counts) # proteins x cells

# remove cells labeled as "Doublet"
pbmc <- subset(pbmc, idents = 'Doublet', invert=TRUE)

# getting UMAP plot
pbmc@reductions$umap
DimPlot(pbmc, reduction = "umap", group.by = 'celltype.l2', label = TRUE) + NoLegend()

pdf(file = "/home/jmitchel/figures/for_paper_v2/cite_seq_umap.pdf", useDingbats = FALSE,
    width = 6, height = 5.5)
DimPlot(pbmc, reduction = "umap", group.by = 'celltype.l3', label = TRUE, repel = TRUE,
        label.size = 2.85) + NoLegend()
dev.off()

## trying with my own collapse fn from scITD
cp <- get_sums(t(pbmc@assays$SCT@counts),as.factor(pbmc@meta.data$celltype.l3))
cp <- cp[2:nrow(cp),]
dim(cp)

# normalize and log transform the pseudo-cells
cp <- t(cp)
lib_sizes <- Matrix::colSums(cp)
cp <- sweep(cp,MARGIN=2,lib_sizes,FUN='/') * 10000
# log transform result
cp <- log1p(cp)

# now cluster and draw dendrogram
cp_cor <- cor(cp)
hc <- hclust(as.dist(1-cp_cor)) 
plot(hc, hang = -1, cex = 0.6)

tree_height_range <- seq(0,.45,.001)
n_comps <- 10
n_top_genes <- 2000
nclust_checked <- c()
all_ratios <- c()
tracked_height <- c()
for (height in tree_height_range) {
  clusts <- cutree(hc, h = height)
  n_clusts <- length(unique(clusts))
  if (!(n_clusts %in% nclust_checked)) {
    print(n_clusts)
    nclust_checked <- c(nclust_checked,n_clusts) #add num clusts to checked vec
    
    ## rename cell idents by new clusters
    pbmc@meta.data['tmp_clusts'] <- sapply(1:nrow(pbmc@meta.data), function(x){
      orig_ct <- pbmc@meta.data[x,'celltype.l3']
      new_clust <- as.character(clusts[orig_ct])
      names(new_clust) <- NULL
      return(new_clust)
    })
    
    ## compute intra-individual variability at this clustering
    av_self <- calc_mean_self_dists(pbmc,n_comps,n_top_genes)
    
    ## compute inter-individual variability at this clustering
    n_donors <- length(unique(pbmc@meta.data$donor))
    av_btwn <- calc_mean_btwn_dist(pbmc,n_comps*n_donors,n_top_genes)
    
    print('av self:')
    print(av_self)
    print('av btwn:')
    print(av_btwn)
    
    # store result
    all_ratios <- c(all_ratios,av_self/av_btwn)
    tracked_height <- c(tracked_height,height)
  }
}

# saveRDS(list(all_ratios,tracked_height),file='/home/jmitchel/data/cite_seq_data/ratios_heights.rds')
test <- readRDS(file='/home/jmitchel/data/cite_seq_data/ratios_heights.rds')


calc_mean_self_dists <- function(pbmc,n_comps,n_top_genes) {
  all_d <- unique(pbmc@meta.data$donor)
  all_self_dists <- c()
  for (d in all_d) {
    print(d)
    cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$donor==d]
    pbmc_sub <- subset(pbmc,cells = cells_keep)
    d_self_dists <- calc_self_dist(pbmc_sub,n_comps,n_top_genes)
    all_self_dists <- c(all_self_dists,d_self_dists)
    print('')
    print('')
  }
  mean_self_dists <- mean(all_self_dists)
  return(mean_self_dists)
}

calc_self_dist <- function(pbmc_sub,n_comps,n_top_genes) {
  counts <- pbmc_sub@assays$SCT@counts
  tmp_clusts <- unique(pbmc_sub@meta.data$tmp_clusts)
  all_dists <- c()
  for (cl in tmp_clusts) {
    print(cl)
    if (!(cl %in% unique(pbmc_sub@meta.data$tmp_clusts))) {
      next
    } else if (table(pbmc_sub@meta.data$tmp_clusts)[cl]<10) {
      next
    }
    cells_of_cl <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$tmp_clusts==cl]
    d_cl_dists <- mclapply(1:n_comps,function(k){
      cells_test <- sample(cells_of_cl,2)
      
      # get matrix of just the cells to compare
      tmp <- counts[,cells_test]
      
      # downsample cells to same lib sizes
      lib_sizes <- colSums(tmp)
      cell_to_downsamp <- order(lib_sizes,decreasing=TRUE)[1]
      n.molecules <- min(lib_sizes)
      profile <- tmp[,cell_to_downsamp]
      downsampled_counts <- rmultinom(1, n.molecules, profile/sum(profile))
      tmp[,cell_to_downsamp] <- downsampled_counts
      tmp <- t(tmp)
      
      # limit to top n expressed genes
      if (n_top_genes < ncol(tmp)) {
        tmp <- tmp[,rank(-colSums(tmp)) <= n_top_genes]
      }
      
      # normalize cells
      tmp <- as.matrix(t(tmp/pmax(1,rowSums(tmp))))
      
      # compute JS divergence
      dist.mat <- jsDist(tmp)

      # return results
      return(dist.mat[1,2])
    }, mc.cores=5)
    d_cl_dists <- unlist(d_cl_dists)
    all_dists <- c(all_dists,d_cl_dists)
  }
  return(all_dists)
}


calc_mean_btwn_dist <- function(pbmc,n_comps,n_top_genes) {
  counts <- pbmc@assays$SCT@counts
  tmp_clusts <- unique(pbmc@meta.data$tmp_clusts)
  all_dists <- c()
  for (cl in tmp_clusts) {
    print(cl)
    cells_of_cl <- rownames(pbmc@meta.data)[pbmc@meta.data$tmp_clusts==cl]
    d_cl_dists <- mclapply(1:n_comps,function(k){
      cells_test <- sample(cells_of_cl,2)
      
      # sample until get cells from different donors
      d_check <- pbmc@meta.data[cells_test,'donor']
      while (d_check[1]==d_check[2]) {
        cells_test <- sample(cells_of_cl,2)
        d_check <- pbmc@meta.data[cells_test,'donor']
      }
      
      # get matrix of just the cells to compare
      tmp <- counts[,cells_test]
      
      # downsample cells to same lib sizes
      lib_sizes <- colSums(tmp)
      cell_to_downsamp <- order(lib_sizes,decreasing=TRUE)[1]
      n.molecules <- min(lib_sizes)
      profile <- tmp[,cell_to_downsamp]
      downsampled_counts <- rmultinom(1, n.molecules, profile/sum(profile))
      tmp[,cell_to_downsamp] <- downsampled_counts
      tmp <- t(tmp)
      
      # limit to top n expressed genes
      if (n_top_genes < ncol(tmp)) {
        tmp <- tmp[,rank(-colSums(tmp)) <= n_top_genes]
      }
      
      # normalize cells
      tmp <- as.matrix(t(tmp/pmax(1,rowSums(tmp))))
      
      # compute JS divergence
      dist.mat <- jsDist(tmp)
      
      # return results
      return(dist.mat[1,2])
    }, mc.cores=15)
    d_cl_dists <- unlist(d_cl_dists)
    all_dists <- c(all_dists,d_cl_dists)
  }
  mean_dists <- mean(all_dists)
  return(mean_dists)
}

library(ggdendro)
library(ggplot2)
# library(tree)
## testing out plotting the dendrogram with a bar plot
dhc <- as.dendrogram(hc)
# Rectangular lines
ddata <- dendro_data(dhc, type = "rectangle")
p <- ggplot() + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), data=segment(ddata)) + 
  geom_text(aes(x = x, y = y-.001, label = label, angle = 0, hjust = 1), data=label(ddata), size=2.5) +
  scale_y_continuous(expand = c(0.15, 0)) +
  ylim(-.2,.45) +
  coord_flip() +
  theme_dendro()
p

# p2 <- matrix(ncol=2,nrow=2)
p2 <- cbind.data.frame(all_ratios,tracked_height)
colnames(p2) <- c('x','y')
# p2[,'x'] <- c(10,20)
# p2[,'y'] <- c(.1,.3)
p2 <- as.data.frame(p2)
p2_dat <- p2
p2 <- ggplot(p2_dat,aes(x=x,y=y)) +
  geom_point() +
  ylim(-.2,.45) +
  xlab('intra indv / inter indv dist') +
  ylab('dendrogram height') +
  coord_flip()
p2

pdf(file = "/home/jmitchel/figures/for_paper_v2/cite_seq_dendro.pdf", useDingbats = FALSE,
    width = 5, height = 7.75)
cowplot::plot_grid(p, p2, ncol = 1, align = "v",rel_heights = c(.7,.3))
dev.off()
