
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


tree_height_range <- seq(0,.12,.02)
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


calc_mean_self_dists <- function(pbmc,n_comps,n_top_genes) {
  all_d <- unique(pbmc@meta.data$donor)
  all_self_dists <- mclapply(1:length(all_d), function(i){
    d <- all_d[i]
    cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$donor==d]
    pbmc_sub <- subset(pbmc,cells = cells_keep)
    d_self_dists <- calc_self_dist(pbmc_sub,n_comps,n_top_genes)
    all_self_dists <- c(all_self_dists,d_self_dists)
  }, mc.cores=20)
  all_self_dists <- unlist(all_self_dists)
  mean_self_dists <- mean(all_self_dists)
  return(mean_self_dists)
}

calc_self_dist <- function(pbmc_sub,n_comps,n_top_genes) {
  counts <- pbmc_sub@assays$RNA@counts
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
    all_dists <- c(all_dists,d_cl_dists)
  }
  return(all_dists)
}


calc_mean_btwn_dist <- function(pbmc,n_comps,n_top_genes) {
  counts <- pbmc@assays$RNA@counts
  tmp_clusts <- unique(pbmc@meta.data$tmp_clusts)
  all_dists <- c()
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
    all_dists <- c(all_dists,d_cl_dists)
  }
  mean_dists <- mean(all_dists)
  return(mean_dists)
}





