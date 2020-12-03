
library(Seurat)
library(SeuratDisk)

# convert h5ad to h5seurat
Convert("/home/jmitchel/data/lupus_data/Lupus_study_adjusted.h5ad", dest = "h5seurat", overwrite = TRUE)

# read the h5seurat file
pbmc <- LoadH5Seurat("/home/jmitchel/data/lupus_data/Lupus_study_adjusted.h5seurat")

# normalize the data
pbmc <- NormalizeData(pbmc)

# will use pbmc object directly in lupus_analysis script...



