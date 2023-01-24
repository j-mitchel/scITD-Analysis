
library(Seurat)
library(readxl)

all_direc <- list.dirs('/home/jmitchel/data/pediatric_lupus/scRNA_raw')
all_direc <- all_direc[2:length(all_direc)]
names(all_direc) <- sapply(all_direc,function(x){
  strsplit(x,split='/')[[1]][[7]]
})

dat <- Read10X(
  all_direc,
  gene.column = 2,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE)

pbmc <- CreateSeuratObject(counts=dat)

colnames(pbmc@meta.data)[1] <- 'donors'
Idents(pbmc)<-rep(1,ncol(dat))

# add meta data
clin_vars <- read_excel('/home/jmitchel/data/pediatric_lupus/pediatric_lupus_clinical.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$Names
clin_vars$Names <- NULL
tmp_meta <- clin_vars[as.character(pbmc@meta.data$donors),]
rownames(tmp_meta) <- rownames(pbmc@meta.data)
pbmc@meta.data <- cbind.data.frame(pbmc@meta.data,tmp_meta)


## subsetting to just pediatric cases/controls for now
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$Groups %in% c('cSLE','cHD')]
pbmc_pediatric <- subset(pbmc,cells=cells_keep)
pbmc_pediatric@meta.data$donors <- as.character(pbmc_pediatric@meta.data$donors)

# save as h5ad
library(SeuratDisk)
# SaveH5Seurat(pbmc_pediatric, filename = "/home/jmitchel/data/pediatric_lupus/processed/pbmc_pediatric.h5Seurat")
# Convert("/home/jmitchel/data/pediatric_lupus/processed/pbmc_pediatric.h5Seurat", dest = "h5ad")


### at this point, I just followed their preprocessing pipeline on github: https://github.com/dnehar/SingleCells_SLE_paper/blob/master/cSLE_scanpy_Pipeline.ipynb
Convert("/home/jmitchel/data/pediatric_lupus/processed/pbmc_pediatric_clean_doublet_removed_batch_corrected.h5ad", dest = "h5seurat", overwrite = FALSE)
pbmc_annotated <- LoadH5Seurat("/home/jmitchel/data/pediatric_lupus/processed/pbmc_pediatric_clean_doublet_removed_batch_corrected.h5seurat", assays='counts')

# since the scanpy version doesn't have raw counts anymore, I'll load the original UMI matrix
# and subset to proper cells and add the metadata, and projections to it
pbmc <- LoadH5Seurat("/home/jmitchel/data/pediatric_lupus/processed/pbmc_pediatric.h5Seurat")

# making sure its counts again
pbmc@assays$RNA[1:50,1:10]

# limit to correct cells
pbmc <- subset(pbmc, cells=colnames(pbmc_annotated))

# adding reductions from annotated data
pbmc@reductions$pca <- pbmc_annotated@reductions$pca
pbmc@reductions$umap <- pbmc_annotated@reductions$umap

# adding metadata from annotated data
pbmc@meta.data <- pbmc_annotated@meta.data

# show dim reduction
DimPlot(pbmc, reduction = "umap", group.by = 'clusters_fine')

# # save seurat object
# saveRDS(pbmc,
#         file='/home/jmitchel/data/pediatric_lupus/processed/pbmc_pediatric_clean_annotated_seurat.rds',
#         compress = "xz")

pbmc <- readRDS(file='/home/jmitchel/data/pediatric_lupus/processed/pbmc_pediatric_clean_annotated_seurat.rds')









##### old code below, not used any more

# pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, 
#         group.by = NULL)
# 
# pbmc <- subset(pbmc, subset = nFeature_RNA > 400 & nFeature_RNA < 2500 & percent.mt < 20)
# 
# print(dim(pbmc))
# 
# # standard Seurat pipeline
# pbmc = NormalizeData(pbmc)
# pbmc = FindVariableFeatures(pbmc, verbose = F)
# pbmc = ScaleData(pbmc, vars.to.regress = c("nFeature_RNA", "percent.mt"),
#                       verbose = F)
# pbmc = RunPCA(pbmc, verbose = F, npcs = 20)
# pbmc = RunUMAP(pbmc, dims = 1:10, verbose = F)
# 
# DimPlot(pbmc, reduction = "umap")
# 
# # running doubletFinder
# nExp <- round(ncol(data.filt) * 0.04) 
# data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)


