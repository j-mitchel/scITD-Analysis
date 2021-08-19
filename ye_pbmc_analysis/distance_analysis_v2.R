
library(Seurat)
library(devtools)
library(sccore)

pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')

# rename the donor meta column name
donor_ndx <- which(colnames(pbmc@meta.data)=='ind_cov_batch_cov')
colnames(pbmc@meta.data)[donor_ndx] <- 'donor'

# use the finer resolution clustering here - see what the column name is
unique(pbmc@meta.data$ct_cov)

DimPlot(pbmc, reduction = "umap", group.by = 'ct_cov', label = TRUE) + NoLegend()

# get pseudobulk cells
cp <- get_sums(t(pbmc@assays$RNA@counts),as.factor(pbmc@meta.data$ct_cov))
cp <- cp[2:nrow(cp),]
dim(cp)

# normalize and log transform the pseudo-cells
cp <- t(cp)
lib_sizes <- Matrix::colSums(cp)
cp <- sweep(cp,MARGIN=2,lib_sizes,FUN='/') * 100000
# log transform result
cp <- log1p(cp)


# now cluster and draw dendrogram
cp_cor <- cor(cp)
hc <- hclust(as.dist(1-cp_cor)) 
plot(hc, hang = -1, cex = 0.6)


tree_height_range <- seq(0,.12,.015)
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
      orig_ct <- pbmc@meta.data[x,'ct_cov']
      new_clust <- as.character(clusts[orig_ct])
      names(new_clust) <- NULL
      return(new_clust)
    })
    
    ## compute intra-individual variability per cluster
    dists_self <- calc_mean_self_dists(pbmc,n_comps,n_top_genes)
    
    ## compute inter-individual variability at this clustering
    n_donors <- length(unique(pbmc@meta.data$donor))
    av_btwn <- calc_mean_btwn_dist(pbmc,n_comps*10,n_top_genes)
    av_btwn <- unlist(av_btwn)
    names(av_btwn) <- as.character(names(av_btwn))
    
    ## calculate cluster-wise ratios
    all_cl_ratios <- dists_self/av_btwn[names(dists_self)]
    av_ratio <- mean(all_cl_ratios,na.rm=TRUE)
    
    print('all ratios:')
    print(all_cl_ratios)
    print('av_ratio:')
    print(av_ratio)
    
    # store summary results
    all_ratios <- c(all_ratios,av_ratio)
    tracked_height <- c(tracked_height,height)
  }
}


calc_mean_self_dists <- function(pbmc,n_comps,n_top_genes) {
  all_d <- unique(pbmc@meta.data$donor)
  all_d <- sample(all_d,20)
  my_clusts <- unique(pbmc@meta.data$tmp_clusts)
  myres <- matrix(NA,nrow=length(all_d),ncol=length(my_clusts)) # to store results
  rownames(myres) <- all_d
  colnames(myres) <- my_clusts
  
  all_self_dists <- lapply(1:length(all_d), function(i){
    d <- all_d[i]
    print(d)
    cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$donor==d]
    pbmc_sub <- subset(pbmc,cells = cells_keep)
    d_self_dists <- calc_self_dist(pbmc_sub,n_comps,n_top_genes)
    d_self_dists <- unlist(d_self_dists)
    names(d_self_dists) <- as.character(names(d_self_dists))
    return(d_self_dists)
  })
  
  names(all_self_dists) <- all_d
  
  ## loop through lists of results and populate the results matrix
  for (i in 1:length(all_self_dists)) {
    d <- names(all_self_dists)[i]
    d_res <- all_self_dists[[d]]
    myres[d,names(d_res)] <- d_res
  }
  ##
  myres_av <- colMeans(myres,na.rm=TRUE) # means per cluster
  return(myres_av)
}

calc_self_dist <- function(pbmc_sub,n_comps,n_top_genes) {
  counts <- pbmc_sub@assays$RNA@counts
  tmp_clusts <- unique(pbmc_sub@meta.data$tmp_clusts)
  all_dists <- list()
  for (cl in tmp_clusts) {
    print(cl)
    if (!(cl %in% unique(pbmc_sub@meta.data$tmp_clusts))) {
      next
    } else if (table(pbmc_sub@meta.data$tmp_clusts)[cl]<10) {
      next
    }
    cells_of_cl <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$tmp_clusts==cl]
    d_cl_dists <- lapply(1:n_comps,function(k){
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
    })
    d_cl_dists <- unlist(d_cl_dists)
    all_dists[[cl]] <- mean(d_cl_dists)
  }
  return(all_dists)
}


calc_mean_btwn_dist <- function(pbmc,n_comps,n_top_genes) {
  counts <- pbmc@assays$RNA@counts
  tmp_clusts <- unique(pbmc@meta.data$tmp_clusts)
  all_dists <- list()
  for (cl in tmp_clusts) {
    print(cl)
    cells_of_cl <- rownames(pbmc@meta.data)[pbmc@meta.data$tmp_clusts==cl]
    d_cl_dists <- lapply(1:n_comps,function(k){
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
    })
    d_cl_dists <- unlist(d_cl_dists)
    all_dists[[cl]] <- mean(d_cl_dists)
  }
  return(all_dists)
}


saveRDS(list(all_ratios,tracked_height),file='/home/jmitchel/data/lupus_data/ratios_heights.rds')



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
  ylim(-.1,.15) +
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
  ylim(-.1,.15) +
  xlab('intra indv / inter indv dist') +
  ylab('dendrogram height') +
  coord_flip()
p2

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_dendro.pdf", useDingbats = FALSE,
    width = 5, height = 7.75)
cowplot::plot_grid(p, p2, ncol = 1, align = "v",rel_heights = c(.7,.3))
dev.off()

