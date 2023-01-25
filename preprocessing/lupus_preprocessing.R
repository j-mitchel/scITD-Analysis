library(Seurat)
library(SeuratDisk)

# Download the H5AD file here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174188
# resave it with gzip compression and free of NaN values as done in preprocessing/resave_sle_data.ipynb

# convert H5AD to h5seurat
Convert("/home/jmitchel/data/temp_test_data/GSE174188_CLUES1_adjusted2.h5ad", dest = "h5seurat", overwrite = FALSE)

# read the h5seurat file
pbmc <- LoadH5Seurat("/home/jmitchel/data/temp_test_data/GSE174188_CLUES1_adjusted2.h5seurat")

# convert age to numeric
pbmc@meta.data$Age <- as.numeric(as.character(pbmc@meta.data$Age))

# remove NA levels from metadata where not present
vars_rem_all <- c('batch_cov','ind_cov','Processing_Cohort','louvain','cg_cov',
                  'ind_cov_batch_cov','pop_cov','Status','SLE_status')
for (var_rem in vars_rem_all) {
  meta_var_counts <- table(pbmc@meta.data[,var_rem])
  pbmc@meta.data[,var_rem] <- factor(pbmc@meta.data[,var_rem],levels=names(meta_var_counts)[meta_var_counts>0])
}

# re-level processing cohort variable
levels(pbmc@meta.data$Processing_Cohort) <- c('L1_Broad', 'L1_UCSF', 'L2', 'L3')


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
donors <- unique(as.character(pbmc@meta.data$ind_cov))
print(length(donors))
for (d in donors) {
  print(d)
  # subset to the donor
  tmp <- pbmc@meta.data[pbmc@meta.data$ind_cov==d,]
  # get number of batches donor is in
  if (sum(table(tmp$batch_cov)>0)>1) {
    # get the donor's status
    d_stat <- as.character(tmp$SLE_status[1])
    
    # get fractions of statuses in each batch
    all_fracs <- list()
    batch_pres <- names(table(tmp$batch_cov))[table(tmp$batch_cov)>0]
    for (b in batch_pres) {
      tmp2 <- pbmc@meta.data[pbmc@meta.data$batch_cov==b,]
      fracs <- table(tmp2$SLE_status)/sum(table(tmp2$SLE_status))
      all_fracs[[b]] <- fracs[d_stat]
    }
    
    # keep sample where there is lowest fraction
    all_fracs_vec <- unlist(all_fracs)
    names(all_fracs_vec) <- names(all_fracs)
    ndx_rem <- which(all_fracs_vec != min(all_fracs_vec))
    batches_rem <- names(all_fracs_vec)[ndx_rem]
    samps_rem <- sapply(batches_rem,function(x) {
      paste0(d,':',x)
    })
    
    cells_rem <- rownames(pbmc@meta.data[pbmc@meta.data$ind_cov_batch_cov %in% samps_rem,])
    pbmc <- subset(pbmc, cells = cells_rem, invert=TRUE)
  }
}

# checking that still have same number of donors
donors <- unique(as.character(pbmc@meta.data$ind_cov))
print(length(donors))

# # saving subsetted seurat object so can revisit this step later and try different normalization schemes
# saveRDS(pbmc,file='/home/jmitchel/data/lupus_data/lupus_subsetted_seurat.rds',compress = "xz")





# pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat.rds')

# check that donors only in one batch
donors <- unique(pbmc@meta.data$ind_cov)
for (d in donors) {
  # subset to the donor
  tmp <- pbmc@meta.data[pbmc@meta.data$ind_cov==d,]
  # get number of batches donor is in
  if (sum(table(tmp$batch_cov)>0)>1) {
    print(sum(table(tmp$batch_cov)>0)>1)
  }
}

# get status breakdown for each batch
batches <- unique(pbmc@meta.data$batch_cov)
for (b in batches) {
  tmp <- pbmc@meta.data[pbmc@meta.data$batch_cov==b,]
  tmp <- tmp[,c('ind_cov','SLE_status')]
  tmp <- unique(tmp)
  rownames(tmp) <- tmp$ind_cov
  print(table(tmp$SLE_status))
}


# # need to normalize the data and correct batch class
# pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat.rds')

# convert batch to factor for meta association analysis
pbmc@meta.data$batch_cov <- factor(pbmc@meta.data$batch_cov,levels=unique(pbmc@meta.data$batch_cov))

# # save new file
# saveRDS(pbmc,file='/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v2.rds',compress = "xz")





# ## need to remove cells with large fraction of UMIs from 1 gene as these are unreliable

# pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v2.rds')
count_data <- pbmc@assays$RNA@counts

thresh <- .3
ncells <- ncol(count_data)
csub <- c(seq(1,ncells,50000),ncells)
all_cells_keep <- c()
for (i in 1:(length(csub)-1)) {
  print(i)
  tmp <- count_data[,(csub[i]+1):csub[i+1]]
  
  maxvals <- apply(tmp,MARGIN=2,FUN=max)
  lib_sizes <- colSums(tmp)
  fracs <- maxvals/lib_sizes
  
  cells_keep <- names(fracs)[fracs < thresh]
  
  all_cells_keep <- c(all_cells_keep,cells_keep)
  
  print(paste0(length(cells_keep),' cells kept'))
  print('')
  
}

pbmc <- subset(pbmc,cells = all_cells_keep)
dim(pbmc@assays$RNA@counts)
dim(pbmc@meta.data)

# rename some variables just so it's consistent with what I did before downstream
ndx_change <- which(colnames(pbmc@meta.data)=='Sex')
colnames(pbmc@meta.data)[ndx_change] <- 'sex'

ndx_change <- which(colnames(pbmc@meta.data)=='pop_cov')
colnames(pbmc@meta.data)[ndx_change] <- 'Ethnicity'

# add the original subclusters
subc_orig <- readRDS(url('http://pklab.med.harvard.edu/jonathan/data/original_subcluster_labels.rds'))
pbmc@meta.data$ct_cov <- subc_orig

## save the final object used in downstream analyses
# saveRDS(pbmc,file='/home/jmitchel/data/temp_test_data/lupus_subsetted_seurat_v3.rds',compress = "xz")






