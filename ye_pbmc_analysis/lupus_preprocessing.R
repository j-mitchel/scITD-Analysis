
library(Seurat)
library(SeuratDisk)

# convert h5ad to h5seurat
Convert("/home/jmitchel/data/lupus_data/Lupus_study_adjusted.h5ad", dest = "h5seurat", overwrite = TRUE)

# read the h5seurat file
pbmc <- LoadH5Seurat("/home/jmitchel/data/lupus_data/Lupus_study_adjusted.h5seurat")

# qc analyses
pbmc[['nCount_RNA']] <- colSums(pbmc[['RNA']]@counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nCount_RNA", "percent.mt"), ncol = 2)

# remove cells with high mitochondrial expression
pbmc <- subset(pbmc, subset = percent.mt < 10)

# normalize the data
pbmc <- NormalizeData(pbmc)

# limit donors to one sample each
cells_keep <- c()
for (d in unique(pbmc@meta.data$Genotype.ID)) {
  print(d)
  # subset meta data to a donor
  meta_sub <- pbmc@meta.data[pbmc@meta.data$Genotype.ID==d,]
  
  # get largest donor's batch
  batch_counts <- table((meta_sub$batch_cov))
  batch_use <- names(batch_counts)[which(batch_counts==max(batch_counts))]
  
  # get cell names corresponding to the donor, batch
  d_cells <- rownames(meta_sub)[meta_sub$batch_cov==batch_use]
  cells_keep <- c(cells_keep,d_cells)
}
pbmc <- subset(pbmc,cells = cells_keep)

# normalize using trimmed means edgeR method
samp_ndx <- sample(1:ncol(pbmc[['RNA']]@counts))
pbmc[['RNA']]@counts <- pbmc[['RNA']]@counts[,samp_ndx]
pbmc@meta.data <- pbmc@meta.data[samp_ndx,]
all_nf <- c()
from_cell <- 1
to_cell <- 50000
while (to_cell <= ncol(pbmc[['RNA']]@counts)) {
  print(from_cell)
  nf <- calcNormFactors(pbmc[['RNA']]@counts[,from_cell:to_cell])
  all_nf <- c(all_nf,nf)
  
  from_cell <- to_cell + 1
  if (to_cell== ncol(pbmc[['RNA']]@counts)) {
    to_cell <- ncol(pbmc[['RNA']]@counts) + 10
  } else {
    to_cell <- to_cell + 50000
    if (to_cell > ncol(pbmc[['RNA']]@counts)) {
      to_cell <- ncol(pbmc[['RNA']]@counts)
    }
  }
}

# divide by lib size and multiply by scale factor
lib_sizes <- Matrix::colSums(pbmc[['RNA']]@counts)
pbmc@assays$RNA <- sweep(pbmc[['RNA']]@counts,MARGIN=2,lib_sizes*all_nf,FUN='/') * 10000

# log transform result
pbmc@assays$RNA <- log1p(pbmc@assays$RNA)



# save the cleaned, normalized data
saveRDS()
