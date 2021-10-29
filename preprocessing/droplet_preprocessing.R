

library(SeuratDisk)

### Batch: YE_7-19
### well: 1
Convert("/home/jmitchel/data/lupus_data/SLE_droplets/YE_7-19/YE_7-19-1/raw_gene_bc_matrices_h5.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc1 <- LoadH5Seurat("/home/jmitchel/data/lupus_data/SLE_droplets/YE_7-19/YE_7-19-1/raw_gene_bc_matrices_h5.h5seurat")
colnames(pbmc1@assays$RNA@counts) <- sapply(colnames(pbmc1@assays$RNA@counts), function(x){
  paste0(x,'-well1')
})
names(colnames(pbmc1@assays$RNA@counts)) <- NULL

### Batch: YE_7-19
### well: 2
Convert("/home/jmitchel/data/lupus_data/SLE_droplets/YE_7-19/YE_7-19-2/raw_gene_bc_matrices_h5.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc2 <- LoadH5Seurat("/home/jmitchel/data/lupus_data/SLE_droplets/YE_7-19/YE_7-19-2/raw_gene_bc_matrices_h5.h5seurat")
colnames(pbmc2@assays$RNA@counts) <- sapply(colnames(pbmc2@assays$RNA@counts), function(x){
  paste0(x,'-well2')
})
names(colnames(pbmc2@assays$RNA@counts)) <- NULL

### Batch: YE_7-19
### well: 3
Convert("/home/jmitchel/data/lupus_data/SLE_droplets/YE_7-19/YE_7-19-3/raw_gene_bc_matrices_h5.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc3 <- LoadH5Seurat("/home/jmitchel/data/lupus_data/SLE_droplets/YE_7-19/YE_7-19-3/raw_gene_bc_matrices_h5.h5seurat")
colnames(pbmc3@assays$RNA@counts) <- sapply(colnames(pbmc3@assays$RNA@counts), function(x){
  paste0(x,'-well3')
})
names(colnames(pbmc3@assays$RNA@counts)) <- NULL

### Batch: YE_7-19
### well: 4
Convert("/home/jmitchel/data/lupus_data/SLE_droplets/YE_7-19/YE_7-19-4/raw_gene_bc_matrices_h5.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc4 <- LoadH5Seurat("/home/jmitchel/data/lupus_data/SLE_droplets/YE_7-19/YE_7-19-4/raw_gene_bc_matrices_h5.h5seurat")
colnames(pbmc4@assays$RNA@counts) <- sapply(colnames(pbmc4@assays$RNA@counts), function(x){
  paste0(x,'-well4')
})
names(colnames(pbmc4@assays$RNA@counts)) <- NULL

## reduce data to just empty droplets (<10 UMIs)
batch_dat <- list(pbmc1,pbmc2,pbmc3,pbmc4)
batch_dat_sub <- list()
for (i in 1:length(batch_dat)) {
  b_cur <- batch_dat[[i]]
  b_cur_mat <- b_cur@assays$RNA@counts
  lib_sizes <- colSums(b_cur_mat)
  bcodes_keep <- colnames(b_cur_mat)[lib_sizes<10]
  b_cur_mat <- b_cur_mat[,bcodes_keep]
  batch_dat_sub[[i]] <- b_cur_mat
}


gene_totals <- rowSums(batch_dat_sub[[1]]) + rowSums(batch_dat_sub[[2]]) +
  rowSums(batch_dat_sub[[3]]) + rowSums(batch_dat_sub[[4]])
gene_fracs <- gene_totals / sum(gene_totals)

soupProf = data.frame(row.names = rownames(batch_dat_sub[[1]]), est = gene_fracs)
# saveRDS(soupProf,file="/home/jmitchel/data/lupus_data/SLE_droplets/YE_7-19/soupProf.rds")
















### Batch: YE_8-16
### well: 1
Convert("/home/jmitchel/data/lupus_data/SLE_droplets/YE_8-16/YE_8-16-1/raw_gene_bc_matrices_h5.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc1 <- LoadH5Seurat("/home/jmitchel/data/lupus_data/SLE_droplets/YE_8-16/YE_8-16-1/raw_gene_bc_matrices_h5.h5seurat")
colnames(pbmc1@assays$RNA@counts) <- sapply(colnames(pbmc1@assays$RNA@counts), function(x){
  paste0(x,'-well1')
})
names(colnames(pbmc1@assays$RNA@counts)) <- NULL

### Batch: YE_8-16
### well: 2
Convert("/home/jmitchel/data/lupus_data/SLE_droplets/YE_8-16/YE_8-16-2/raw_gene_bc_matrices_h5.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc2 <- LoadH5Seurat("/home/jmitchel/data/lupus_data/SLE_droplets/YE_8-16/YE_8-16-2/raw_gene_bc_matrices_h5.h5seurat")
colnames(pbmc2@assays$RNA@counts) <- sapply(colnames(pbmc2@assays$RNA@counts), function(x){
  paste0(x,'-well2')
})
names(colnames(pbmc2@assays$RNA@counts)) <- NULL

### Batch: YE_8-16
### well: 3
Convert("/home/jmitchel/data/lupus_data/SLE_droplets/YE_8-16/YE_8-16-3/raw_gene_bc_matrices_h5.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc3 <- LoadH5Seurat("/home/jmitchel/data/lupus_data/SLE_droplets/YE_8-16/YE_8-16-3/raw_gene_bc_matrices_h5.h5seurat")
colnames(pbmc3@assays$RNA@counts) <- sapply(colnames(pbmc3@assays$RNA@counts), function(x){
  paste0(x,'-well3')
})
names(colnames(pbmc3@assays$RNA@counts)) <- NULL

### Batch: YE_8-16
### well: 4
Convert("/home/jmitchel/data/lupus_data/SLE_droplets/YE_8-16/YE_8-16-4/raw_gene_bc_matrices_h5.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc4 <- LoadH5Seurat("/home/jmitchel/data/lupus_data/SLE_droplets/YE_8-16/YE_8-16-4/raw_gene_bc_matrices_h5.h5seurat")
colnames(pbmc4@assays$RNA@counts) <- sapply(colnames(pbmc4@assays$RNA@counts), function(x){
  paste0(x,'-well4')
})
names(colnames(pbmc4@assays$RNA@counts)) <- NULL

## reduce data to just empty droplets (<10 UMIs)
batch_dat <- list(pbmc1,pbmc2,pbmc3,pbmc4)
batch_dat_sub <- list()
for (i in 1:length(batch_dat)) {
  b_cur <- batch_dat[[i]]
  b_cur_mat <- b_cur@assays$RNA@counts
  lib_sizes <- colSums(b_cur_mat)
  bcodes_keep <- colnames(b_cur_mat)[lib_sizes<10]
  b_cur_mat <- b_cur_mat[,bcodes_keep]
  batch_dat_sub[[i]] <- b_cur_mat
}

gene_totals <- rowSums(batch_dat_sub[[1]]) + rowSums(batch_dat_sub[[2]]) +
  rowSums(batch_dat_sub[[3]]) + rowSums(batch_dat_sub[[4]])
gene_fracs <- gene_totals / sum(gene_totals)

soupProf = data.frame(row.names = rownames(batch_dat_sub[[1]]), est = gene_fracs)
# saveRDS(soupProf,file="/home/jmitchel/data/lupus_data/SLE_droplets/YE_8-16/soupProf.rds")





