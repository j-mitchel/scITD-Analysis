
library(scITD)
library(Seurat)

# the data were downloaded from here:
# https://cellxgene.cziscience.com/collections/7d7cabfd-1d1f-40af-96b7-26a0825a306d
pbmc <- readRDS('/home/jmitchel/data/covid_data_ye/local.rds')

# it seems in their paper they say the Case patients negative for COVID19 often had
# other bacteria present so they're likely not really C19 patients

## need to remove the covid negative cases and limit each sample to the first time point D0 where severity was measured
## if I need more samples, can use the the time series samples each as separate samples,
## though I should use NIH ordinal measure for severity
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$COVID_status!='Negative']
pbmc_sub <- subset(pbmc,cells=cells_keep)
pbmc_sub@meta.data[1:100,c('COVID_status','timepoint')]

cells_keep <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$timepoint %in% c(0,'None')]
pbmc_sub <- subset(pbmc_sub,cells=cells_keep)


## need to see that donor samples only in 1 batch...

meta <- pbmc_sub@meta.data[,c('donor','batch','pool','well')]
meta <- unique(meta)
table(meta)

# I'm unable to just pick one donor per batch at well-level because there would be very few
# cells per donor per cell type that way. For now I'll use pool-level batch correction.
# Still need to check it's 1 donor per pool. If this doesn't give me good results, I could
# try doing batch correction at well-level on individual cells.

meta <- pbmc_sub@meta.data[,c('donor','pool')]
meta <- unique(meta)
table(meta)

meta <- pbmc_sub@meta.data[,c('donor','pool','COVID_status')]
meta <- unique(meta)
meta <- meta[meta$COVID_status=='Healthy',c('donor','pool')]
table(meta$pool)

meta <- pbmc_sub@meta.data[,c('donor','pool','COVID_status')]
meta <- unique(meta)
meta <- meta[meta$COVID_status=='Positive',c('donor','pool')]
table(meta$pool)

meta <- pbmc_sub@meta.data[,c('donor','pool','COVID_status')]
meta <- unique(meta)
meta <- meta[meta$COVID_status=='Healthy',c('donor','pool')]
meta$donor <- factor(meta$donor,levels=unique(meta$donor))
meta$pool <- factor(meta$pool,levels=unique(meta$pool))
table(meta)

# remove pool 10 altogether
# keep ICC-C-0177 pool 9 sample
# keep ICC_C_0003 pool 7
# keep ICC_C_0001 pool 3

cells_keep <- c()
for (i in 1:nrow(pbmc_sub@meta.data)){
  cn <- rownames(pbmc_sub@meta.data)[i]
  if (pbmc_sub@meta.data$pool[i]==10) {
    next
  } else if (pbmc_sub@meta.data$donor[i]=='ICC-C-0177' && pbmc_sub@meta.data$pool[i]!=9) {
    next
  } else if (pbmc_sub@meta.data$donor[i]=='ICC_C_0003' && pbmc_sub@meta.data$pool[i]!=7) {
    next
  } else if (pbmc_sub@meta.data$donor[i]=='ICC_C_0001' && pbmc_sub@meta.data$pool[i]!=3) {
    next
  } else {
    cells_keep <- c(cells_keep,cn)
  }
}

pbmc_sub <- subset(pbmc_sub,cells=cells_keep)
test_cells <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$donor=='ICC_C_0003']
test <- subset(pbmc_sub,cells=test_cells)

# need to convert gene names to gene symbols...
feature.names <- readRDS('/home/jmitchel/data/van_der_wijst/genes.rds')
colnames(feature.names) <- c('ensg','hugo')
rownames(feature.names) <- feature.names$ensg

dim(pbmc_sub)
new_names <- feature.names[rownames(pbmc_sub),2]
sum(is.na(new_names))
ndx_keep <- which(!is.na(new_names))
pbmc_counts <- pbmc_sub@assays$RNA@counts[ndx_keep,]
new_names <- feature.names[rownames(pbmc_counts),2]
# get rid of duplicates
dup1=duplicated(new_names)
dup2=duplicated(new_names,fromLast=TRUE)
ndx1 <- which(dup1)
ndx2 <- which(dup2)
ndx_keep <- which(!(new_names %in% unique(new_names[ndx1])))
pbmc_counts <- pbmc_counts[ndx_keep,]
rownames(pbmc_counts) <- feature.names[rownames(pbmc_counts),2]

pbmc_meta <- pbmc_sub@meta.data
pbmc_sub2 <- CreateSeuratObject(pbmc_counts, meta.data=pbmc_meta)



saveRDS(pbmc_sub2,file='/home/jmitchel/data/covid_data_ye/covid_subset.rds')






