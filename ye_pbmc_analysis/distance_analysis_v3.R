library(Seurat)
library(SeuratDisk)
library(sccore)
library(gdata)

pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')

# rename the donor meta column name
donor_ndx <- which(colnames(pbmc@meta.data)=='ind_cov_batch_cov')
colnames(pbmc@meta.data)[donor_ndx] <- 'donor'

# use the finer resolution clustering here - see what the column name is
unique(pbmc@meta.data$ct_cov)

DimPlot(pbmc, reduction = "umap", group.by = 'ct_cov', label = TRUE) + NoLegend()
DimPlot(pbmc, reduction = "umap", group.by = 'donor', label = FALSE) + NoLegend()

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


tree_height_range <- seq(0,.126,.015)
n_comps <- 1
n_top_genes <- 100
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
    av_btwn <- calc_mean_btwn_dist(pbmc,n_comps,n_top_genes)
    av_btwn <- unlist(av_btwn)
    names(av_btwn) <- as.character(names(av_btwn))
    
    ## calculate cluster-wise ratios
    all_cl_ratios <- dists_self/av_btwn[names(dists_self)]
    all_cl_ratios <- all_cl_ratios[!is.na(all_cl_ratios)]
    
    # get cluster sizes
    c_counts <- table(pbmc@meta.data$tmp_clusts)
    c_counts <- c_counts[names(all_cl_ratios)]
    c_sizes <- c_counts/sum(c_counts)
    
    # get weighted av ratio across clusters
    av_ratio <- sum(all_cl_ratios*c_sizes)
    
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
      cells_test <- sample(cells_of_cl,min(200,length(cells_of_cl)))
      
      # get matrix of just the cells to compare
      tmp <- counts[,cells_test]
      
      # downsample cells to same lib sizes
      lib_sizes <- colSums(tmp)
      n.molecules <- min(lib_sizes)
      # loop through cells and downsample all to n.molecules
      for (cell_ndx in 1:ncol(tmp)) {
        profile <- tmp[,cell_ndx]
        downsampled_counts <- rmultinom(1, n.molecules, profile/sum(profile))
        tmp[,cell_ndx] <- downsampled_counts
      }
      
      tmp <- t(tmp)
      
      # limit to top n expressed genes
      if (n_top_genes < ncol(tmp)) {
        tmp <- tmp[,rank(-colSums(tmp)) <= n_top_genes]
      }
      
      # normalize cells
      tmp <- as.matrix(t(tmp/pmax(1,rowSums(tmp))))
      
      # ## trying log with cor dist
      # tmp <- log1p(tmp)
      # dist.mat <- cor(tmp)
      
      # compute JS divergence
      dist.mat <- jsDist(tmp)
      
      upper_tri <- upperTriangle(dist.mat)
      
      # return results
      return(upper_tri)
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
      cells_test <- sample(cells_of_cl,min(500,length(cells_of_cl)))
      
      # get matrix of just the cells to compare
      tmp <- counts[,cells_test]
      
      # downsample cells to same lib sizes
      lib_sizes <- colSums(tmp)
      n.molecules <- min(lib_sizes)
      # loop through cells and downsample all to n.molecules
      for (cell_ndx in 1:ncol(tmp)) {
        profile <- tmp[,cell_ndx]
        downsampled_counts <- rmultinom(1, n.molecules, profile/sum(profile))
        tmp[,cell_ndx] <- downsampled_counts
      }
      tmp <- t(tmp)
      
      # limit to top n expressed genes
      if (n_top_genes < ncol(tmp)) {
        tmp <- tmp[,rank(-colSums(tmp)) <= n_top_genes]
      }
      
      # normalize cells
      tmp <- as.matrix(t(tmp/pmax(1,rowSums(tmp))))
      
      # ## trying log with cor dist
      # tmp <- log1p(tmp)
      # dist.mat <- cor(tmp)
      
      # compute JS divergence
      dist.mat <- jsDist(tmp)
      
      # get donors each cell belongs to
      colnames(dist.mat) <- pbmc@meta.data[cells_test,'donor']
      rownames(dist.mat) <- pbmc@meta.data[cells_test,'donor']
      
      # make cell pairs of same donors to be NA
      dist.mat[outer(rownames(dist.mat), colnames(dist.mat), "==")] <- NA
      
      # get upper triangle
      upper_tri <- upperTriangle(dist.mat)
      upper_tri <- upper_tri[!is.na(upper_tri)]
      
      # return results
      return(upper_tri)
    })
    d_cl_dists <- unlist(d_cl_dists)
    all_dists[[cl]] <- mean(d_cl_dists)
  }
  return(all_dists)
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
  ylim(-.01,.13) +
  coord_flip() +
  xlab('') +
  ylab('')
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
  ylim(-.01,.13) +
  xlab('intra indv / inter indv dist') +
  ylab('dendrogram height') +
  coord_flip()
p2

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_dendro.pdf", useDingbats = FALSE,
    width = 5, height = 7.75)
cowplot::plot_grid(p, p2, ncol = 1, align = "v",rel_heights = c(.7,.3))
dev.off()



# saveRDS(list(all_ratios,tracked_height),file='/home/jmitchel/data/lupus_data/sle_ratios_heights.rds')




## looking at umap of a handful of donors to assess inter-individual variation
all_d <- unique(pbmc@meta.data$donor)
all_d <- sample(all_d,3)
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$donor%in%all_d]
pbmc_sub <- subset(pbmc,cells = cells_keep)
DimPlot(pbmc_sub, reduction = "umap", group.by = 'ct_cov', label = TRUE) + NoLegend()
DimPlot(pbmc_sub, reduction = "umap", group.by = 'donor', label = FALSE)


d <- 'IGTB508_IGTB508dmx_YE_8-2'
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$donor==d]
pbmc_sub <- subset(pbmc,cells = cells_keep)
d_self_dists <- calc_self_dist(pbmc_sub,n_comps,n_top_genes)
d_self_dists <- unlist(d_self_dists)
names(d_self_dists) <- as.character(names(d_self_dists))

# pbmc2 <- pbmc
d <- c('IGTB508_IGTB508dmx_YE_8-2','1165_1165dmx_YE_7-20')
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$donor%in%d]
pbmc <- subset(pbmc,cells = cells_keep)

## compute inter-individual variability at this clustering
av_btwn <- calc_mean_btwn_dist(pbmc,n_comps,n_top_genes)
av_btwn <- unlist(av_btwn)
names(av_btwn) <- as.character(names(av_btwn))

## calculate cluster-wise ratios
all_cl_ratios <- d_self_dists/av_btwn[names(d_self_dists)]
all_cl_ratios













### trying pca version
d <- 'IGTB508_IGTB508dmx_YE_8-2'
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$donor==d]
pbmc_sub <- subset(pbmc,cells = cells_keep)

d <- c('IGTB508_IGTB508dmx_YE_8-2','1165_1165dmx_YE_7-20')
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$donor%in%d]
pbmc_sub <- subset(pbmc,cells = cells_keep)

counts <- pbmc_sub@assays$RNA@counts

cells_of_cl <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$ct_cov=='cM']
cells_of_cl <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$ct_cov=='Bmem']
cells_of_cl <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$ct_cov=='Bnaive']
cells_of_cl <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$ct_cov=='T8em']
cells_of_cl <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$tmp_clusts==1]

cells_test <- sample(cells_of_cl,min(2000,length(cells_of_cl)))

#### other method
pbmc_sub = subset(pbmc_sub,cells=cells_test)
pbmc_sub <- FindVariableFeatures(pbmc_sub, selection.method = "vst", nfeatures = 2000)
pbmc_sub <- ScaleData(pbmc_sub)
pbmc_sub <- RunPCA(pbmc_sub, features = VariableFeatures(object = pbmc_sub))
exp_var <- (pbmc_sub@reductions[["pca"]]@stdev**2)/sum(pbmc_sub@reductions[["pca"]]@stdev**2)
exp_var[1:10]

lds1=pbmc_sub@reductions[["pca"]]@feature.loadings
lds2=pbmc_sub@reductions[["pca"]]@feature.loadings
g_int=intersect(rownames(lds1),rownames(lds2))
cor(lds1[g_int,1],lds2[g_int,1])
###

# get matrix of just the cells to compare
tmp <- counts[,cells_test]

# downsample cells to same lib sizes
lib_sizes <- colSums(tmp)
n.molecules <- min(lib_sizes)
# loop through cells and downsample all to n.molecules
for (cell_ndx in 1:ncol(tmp)) {
  profile <- tmp[,cell_ndx]
  downsampled_counts <- rmultinom(1, n.molecules, profile/sum(profile))
  tmp[,cell_ndx] <- downsampled_counts
}

tmp <- t(tmp)

# limit to top n expressed genes
if (n_top_genes < ncol(tmp)) {
  tmp <- tmp[,rank(-colSums(tmp)) <= n_top_genes]
}

# normalize cells
tmp <- as.matrix(t(tmp/pmax(1,rowSums(tmp))))

# log transform
tmp <- log1p(tmp)

# compute PCA
tmp <- t(tmp)
pres <- prcomp(tmp,center = T,scale. = T,retx = T)
pres[["sdev"]][1:10]

gp1<-pres[["x"]][,1]
gp1_2<-pres[["x"]][,2]
gp2<-pres[["x"]][,1]
gp2_2<-pres[["x"]][,2]

g_int <- intersect(names(gp1),names(gp2))
cor(gp1[g_int],gp2[g_int])

g_int <- intersect(names(gp1_2),names(gp2_2))
cor(gp1_2[g_int],gp2_2[g_int])

## seeing if PC1 is associated with the source donor
pc1 <- pres$x[,1,drop=F]
pc1 <- pbmc_sub@reductions[["pca"]]@cell.embeddings[,1,drop=F]
test <- cbind.data.frame(pc1,pbmc_sub@meta.data[rownames(pc1),'donor'])
colnames(test) <- c('pc','dnr')
lmres <- lm(pc~dnr,data=test)
summary(lmres)

exp_var <- (pres[["sdev"]]**2)/sum(pres[["sdev"]]**2)
exp_var[1:10]
