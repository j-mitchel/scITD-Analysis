
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
# need to do this in a few steps because the matrix is too large
norm_mat_all <- log1p(sweep(pbmc[['RNA']]@counts[,1:50000],
                         MARGIN=2,
                         lib_sizes[1:50000]*all_nf[1:50000],FUN='/') * 10000)
from_cell <- 50001
to_cell <- 100000
while (to_cell <= ncol(pbmc[['RNA']]@counts)) {
  print(from_cell)
  norm_mat <- log1p(sweep(pbmc[['RNA']]@counts[,from_cell:to_cell],
                           MARGIN=2,
                           lib_sizes[from_cell:to_cell]*all_nf[from_cell:to_cell],FUN='/') * 10000)
  norm_mat_all <- cbind(norm_mat_all,norm_mat)
  
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

# Get regular normalization of the data
pbmc <- NormalizeData(pbmc)

pbmc@assays$RNA <- pbmc@assays$RNA[,samp_ndx] # because didnt change order yet

# # save the cleaned, normalized data
# saveRDS(norm_mat_all,file='/home/jmitchel/data/lupus_data/lupus_counts_trim_mean.rds',compress = "xz")
# saveRDS(pbmc@meta.data,file='/home/jmitchel/data/lupus_data/lupus_meta_trim_mean.rds',compress = "xz")
# saveRDS(pbmc@assays$RNA,file='/home/jmitchel/data/lupus_data/lupus_counts_regular_norm.rds',compress = "xz")


# create meta data matrix with only some of the variables
pbmc_meta <- readRDS('/home/jmitchel/data/lupus_data/lupus_meta_all_vars.rds')

# select and rename columns of metadata
metadata_cols=c("Genotype.ID",
                "SLE_status",
                "Status",
                "cg_cov",
                "sex",
                "age",
                "batch_cov",
                "Processing_Cohort")
metadata_col_nm=c('donors',
                  'SLE_status',
                  'Status',
                  'ctypes',
                  'sex',
                  'age',
                  'pool',
                  'processing')
pbmc_meta <- pbmc_meta[,metadata_cols]
colnames(pbmc_meta) <- metadata_col_nm

# saveRDS(pbmc_meta,file='/home/jmitchel/data/lupus_data/lupus_meta_select_vars.rds',compress = "xz")








# for creating smaller dataset with balanced batches only

# read the h5seurat file
pbmc <- LoadH5Seurat("/home/jmitchel/data/lupus_data/Lupus_study_adjusted.h5seurat")

# qc analyses
pbmc[['nCount_RNA']] <- colSums(pbmc[['RNA']]@counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# remove cells with high mitochondrial expression
pbmc <- subset(pbmc, subset = percent.mt < 10)

# get donor status counts for each batch
batches <- unique(as.character(pbmc@meta.data$batch_cov))
for (b in batches) {
  tmp <- pbmc@meta.data[pbmc@meta.data$batch_cov==b,]
  print(b)
  print(table(tmp$Status))
}

# initial pools to remove (ones with only one donor status represented)
# "dmx_YE110", "dmx_count_BH7YT2DMXX_YE_0907", "dmx_count_AH7TNHDMXX_YE_8-30", "dmx_count_AHCM2CDMXX_YE_0831"

batch_rem <- c("dmx_YE110", "dmx_count_BH7YT2DMXX_YE_0907", 
               "dmx_count_AH7TNHDMXX_YE_8-30", "dmx_count_AHCM2CDMXX_YE_0831")
cells_rem <- rownames(pbmc@meta.data[pbmc@meta.data$batch_cov %in% batch_rem,])
pbmc <- subset(pbmc, cells = cells_rem, invert=TRUE)

# next I'll remove the flare processing batches because my sims indicate that if this batch has any batch effect
# it will be very challenging to distinguish real from noise factors
proc_rem <- c("L2")
cells_rem <- rownames(pbmc@meta.data[pbmc@meta.data$Processing_Cohort %in% proc_rem,])
pbmc <- subset(pbmc, cells = cells_rem, invert=TRUE)

# should now only contain cells from L1_UCSF and L3
table(pbmc@meta.data$Processing_Cohort)

# starting with batches that have fewest healthy donor samples represented
# if those donors have samples in other batches, then remove those
# using greedy process of removing donors to minimize imbalance in sets
donors <- unique(as.character(pbmc@meta.data$Genotype.ID))
print(length(donors))
for (d in donors) {
  print(d)
  # subset to the donor
  tmp <- pbmc@meta.data[pbmc@meta.data$Genotype.ID==d,]
  # get number of batches donor is in
  if (sum(table(tmp$batch_cov)>0)>1) {
    # get the donor's status
    d_stat <- as.character(tmp$Status[1])
    
    # get fractions of statuses in each batch
    all_fracs <- list()
    batch_pres <- names(table(tmp$batch_cov))[table(tmp$batch_cov)>0]
    for (b in batch_pres) {
      tmp2 <- pbmc@meta.data[pbmc@meta.data$batch_cov==b,]
      fracs <- table(tmp2$Status)/sum(table(tmp2$Status))
      all_fracs[[b]] <- fracs[d_stat]
    }
    
    # keep sample where there is lowest fraction
    all_fracs_vec <- unlist(all_fracs)
    names(all_fracs_vec) <- names(all_fracs)
    ndx_rem <- which(all_fracs_vec != min(all_fracs_vec))
    batches_rem <- names(all_fracs_vec)[ndx_rem]
    samps_rem <- sapply(batches_rem,function(x) {
      paste0(d,x)
    })
    
    cells_rem <- rownames(pbmc@meta.data[pbmc@meta.data$ind_cov_batch_cov %in% samps_rem,])
    pbmc <- subset(pbmc, cells = cells_rem, invert=TRUE)
  }
}

# checking that still have same number of donors
donors <- unique(as.character(pbmc@meta.data$Genotype.ID))
print(length(donors))

# saving subsetted seurat object so can revisit this step later and try different normalization schemes
saveRDS(pbmc,file='/home/jmitchel/data/lupus_data/lupus_subsetted_seurat.rds',compress = "xz")





pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat.rds')

# check that donors only in one batch
donors <- unique(pbmc@meta.data$Genotype.ID)
for (d in donors) {
  # subset to the donor
  tmp <- pbmc@meta.data[pbmc@meta.data$Genotype.ID==d,]
  # get number of batches donor is in
  if (sum(table(tmp$batch_cov)>0)>1) {
    print(sum(table(tmp$batch_cov)>0)>1)
  }
}

# get status breakdown for each batch
batches <- unique(pbmc@meta.data$batch_cov)
for (b in batches) {
  tmp <- pbmc@meta.data[pbmc@meta.data$batch_cov==b,]
  tmp <- tmp[,c('Genotype.ID','Status')]
  tmp <- unique(tmp)
  rownames(tmp) <- tmp$Genotype.ID
  print(table(tmp$Status))
}


# need to normalize the data and correct batch class
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat.rds')

# Get regular normalization of the data if using h5
pbmc <- NormalizeData(pbmc)

# convert batch to factor for meta association analysis
pbmc@meta.data$batch_cov <- factor(pbmc@meta.data$batch_cov,levels=unique(pbmc@meta.data$batch_cov))

# save new file
saveRDS(pbmc,file='/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v2.rds',compress = "xz")






