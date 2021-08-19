# library(Seurat)
# library("readxl")
# library(Matrix)

# mat <- Matrix::readMM(file = "/home/jmitchel/data/covid_data/GSE158055_covid19_counts.mtx")
# feature.names = read.delim("/home/jmitchel/data/covid_data/GSE158055_covid19_features.tsv",
#                            header = FALSE,
#                            stringsAsFactors = FALSE)
# barcode.names = read.delim("/home/jmitchel/data/covid_data/GSE158055_covid19_barcodes.tsv",
#                            header = FALSE,
#                            stringsAsFactors = FALSE)
# 
# metadata <- as.data.frame(read_excel('/home/jmitchel/data/covid_data/GSE158055_sample_metadata.xlsx'))
# cell_annot <- read.csv(file='/home/jmitchel/data/covid_data/GSE158055_cell_annotation.csv')


# ## trying with the h5ad file
# library(Seurat)
# library(SeuratDisk)
# 
# # convert h5ad to h5seurat
# Convert("/home/jmitchel/data/covid_data/COVID19_subset.h5ad", dest = "h5seurat", overwrite = FALSE, layers=FALSE)
# 
# # read the h5seurat file
# pbmc <- LoadH5Seurat("/home/jmitchel/data/covid_data/COVID19_subset.h5seurat")


## trying to read in the subsetted mtx file from scanpy
mat <- Matrix::readMM(file = "/home/jmitchel/data/covid_data/COVID19_subset.mtx")

mat[1:5,1:5]

print(dim(mat))

lib_sizes <- rowSums(mat)

max(lib_sizes)
boxplot(lib_sizes)

pbmc_meta <- read.csv(file = "/home/jmitchel/data/covid_data/COVID19_meta.csv")
rownames(pbmc_meta) <- pbmc_meta[,1]
pbmc_meta[,1] <- NULL

mat <- t(mat)
dim(mat)
colnames(mat) <- rownames(pbmc_meta)

# adding gene names
features <- read.csv(file = "/home/jmitchel/data/covid_data/COVID19_genes.csv")
features <- features[,2]
rownames(mat) <- features

# checking counts for RPS4Y1 of M vs F donor to make sure gene ordering correct
pbmc_meta$Sex[100000]
mat['RPS4Y1',100000]
mat['RPS4Y1',100]

# create seurat object
pbmc <- CreateSeuratObject(counts = mat, meta.data=pbmc_meta)

rm(mat)
gc()

# saveRDS(pbmc,file="/home/jmitchel/data/covid_data/COVID19_subsetted_seurat.rds",compress = "xz")


# library(Seurat)
# pbmc <- readRDS(file="/home/jmitchel/data/covid_data/COVID19_subsetted_seurat.rds")










