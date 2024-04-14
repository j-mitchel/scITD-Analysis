
library(scTensor)
library(Seurat)
library("AnnotationHub")
library("LRBaseDbi")
library(SingleCellExperiment)
library(Homo.sapiens)
library(scTGIF)
library(devtools)
load_all('/home/jmitchel/scITD/')

# load up the lupus dataset: see preprocessing/lupus_preprocessing.R 
# for code used to generate this object
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')

# converting shorthand cell type names to full names
new_names <- sapply(as.character(pbmc@meta.data$cg_cov), function(x){
  if (x=='cM') {
    return('cMono')
  } else if (x=='ncM') {
    return('ncMono')
  } else if (x=='T4') {
    return('Th')
  } else if (x=='T8') {
    return('Tc')
  } else {
    return(x)
  }
})
names(new_names) <- NULL
pbmc@meta.data$cg_cov <- factor(new_names,levels=unique(new_names))


# subset data to SLE patients only
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$Status=='Managed']
pbmc <- subset(pbmc,cells = cells_keep)

# subset data to just cell types analyzed with scITD
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$cg_cov %in% c("B","NK","Th","Tc","cDC",
                                                                    "cMono","ncMono")]
pbmc <- subset(pbmc,cells = cells_keep)

## converting gene names to ncbi gene id
rowID <- rownames(pbmc)
# saveRDS(rowID,file='/home/jmitchel/sle_genes_tmp.rds')

# Convert gene symbols to NCBI Gene IDs
gene_conv <- readRDS('/home/jmitchel/data/lupus_data/sle_genes_ncbi.rds')

match_ndx <- match(rowID,gene_conv[,'HGNC symbol'])
gene_nm_new <- gene_conv[match_ndx,2]

# remove genes with na gene names from expression data
g_keep_ndx <- which(!is.na(gene_nm_new))
pbmc@assays$RNA@counts <- pbmc@assays$RNA@counts[g_keep_ndx,]
pbmc@assays$RNA@data <- pbmc@assays$RNA@data[g_keep_ndx,]

# removing na elements from ncbi ids
gene_nm_new <- gene_nm_new[!is.na(gene_nm_new)]

# now removing duplicated elements
g_ndx_keep <- which(!duplicated(gene_nm_new))
pbmc@assays$RNA@counts <- pbmc@assays$RNA@counts[g_ndx_keep,]
pbmc@assays$RNA@data <- pbmc@assays$RNA@data[g_ndx_keep,]
gene_nm_new <- gene_nm_new[g_ndx_keep]

# change names of the genes to the ncbi ids
rownames(pbmc@assays$RNA@counts) <- gene_nm_new
rownames(pbmc@assays$RNA@data) <- gene_nm_new

# convert to single-cell experiment object needed for scTensor
pbmc.sce <- as.SingleCellExperiment(pbmc)



## set up lr database to use
ah <- AnnotationHub() # needed to install BiocFileCache from github
# ah <- readRDS('/home/jmitchel/data/lupus_data/sle_ah_file.rds') # saved from locally in case it errors
mcols(ah)

dbfile <- query(ah, c("LRBaseDb", "Homo sapiens", "v002"))[[1]]

LRBase.Hsa.eg.db <- LRBaseDbi::LRBaseDb(dbfile)
key_HSA <- keys(LRBase.Hsa.eg.db, keytype="GENEID_L")
head(select(LRBase.Hsa.eg.db, keys=key_HSA,
            columns=c("GENEID_L", "GENEID_R"), keytype="GENEID_L"))

# saving lr pairs for later scITD run
lr_pairs <- select(LRBase.Hsa.eg.db, keys=key_HSA,
                   columns=c("GENEID_L", "GENEID_R"), keytype="GENEID_L")


## set parameters 
cellCellSetting(pbmc.sce, LRBase.Hsa.eg.db, as.character(pbmc@meta.data$cg_cov)) # third parameter is cell barcodes

# run rank determination
(rks <- cellCellRanks(pbmc.sce))
print(rks$selected)

# run method
pbmc.sce <- cellCellDecomp(pbmc.sce, ranks=rks$selected) # ranks are 3,3

# extract results
l_ct_f <- pbmc.sce@metadata[["sctensor"]][["ligand"]] # ligand factors, by cell types
r_ct_f <- pbmc.sce@metadata[["sctensor"]][["receptor"]] # receptor factors, by cell types
core_tensor <- pbmc.sce@metadata[["sctensor"]][["lrpair"]] # ligand factors, by receptor factors, by lr pairs

# # saving decomposition results
# saveRDS(list(l_ct_f,r_ct_f,core_tensor),file='/home/jmitchel/data/lupus_data/scTensor_res.rds')

lr_pairs_in_tensor <- dimnames(core_tensor@data)[[3]]







### now to run scITD with same lr pairs

# loading up scITD data
pbmc_container <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_container_w_decomp.rds')

# converting LR table to have gene names
lr_pairs2 <- lapply(lr_pairs_in_tensor,function(x){
  mypair <- unlist(strsplit(x,split='_')[[1]])
  
  ndx_gene <- which(gene_conv[,'NCBI gene (formerly Entrezgene) ID']==mypair[1])
  lig <- gene_conv[ndx_gene,'HGNC symbol']
  
  ndx_gene <- which(gene_conv[,'NCBI gene (formerly Entrezgene) ID']==mypair[2])
  rec <- gene_conv[ndx_gene,'HGNC symbol']
  
  mypair <- as.data.frame(t(c(lig,rec)))
  return(mypair)
})

# just removing the two cases with multiple receptors
tester=sapply(lr_pairs2,function(x){
  if (length(x)==2) {
    return(FALSE)
  } else {
    return(TRUE)
  }
})
ndx_keep <- which(!tester)
lr_pairs3 <- lr_pairs2[ndx_keep]
lr_pairs3 <- do.call(rbind.data.frame,lr_pairs3)

# infer active LR interactions
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs3, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='pool')
scITD_scale_pb_scT <- pbmc_container$scale_pb_extra
# saveRDS(scITD_scale_pb_scT,file='/home/jmitchel/data/lupus_data/scITD_scTensor_pb_scaled.rds')

sft_thresh <- c(12,14,12,10,12,9,12)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

lr_hmap <- compute_LR_interact(pbmc_container, lr_pairs3, sig_thresh=.00000000005,
                               percentile_exp_rec=0.85, add_ld_fact_sig=TRUE)

lr_hmap

# saveRDS(pbmc_container$lr_res,file='/home/jmitchel/data/lupus_data/scTensor_scITD_lr_res.rds')



# first, binarize scTensor factors
for (i in 1:nrow(l_ct_f)) {
  med_val <- median(l_ct_f[i,])
  sd_val <- sd(l_ct_f[i,])
  ct_rem_ndx <- which(l_ct_f[i,]<(med_val+sd_val))
  l_ct_f[i,ct_rem_ndx] <- 0
}

for (i in 1:nrow(r_ct_f)) {
  med_val <- median(r_ct_f[i,])
  sd_val <- sd(r_ct_f[i,])
  ct_rem_ndx <- which(r_ct_f[i,]<(med_val+sd_val))
  r_ct_f[i,ct_rem_ndx] <- 0
}

# next select the top 200 ligands with largest max values in core tensor
lr_max_core <- c()
for (i in 1:length(lr_pairs_in_tensor)) {
  max_val <- max(core_tensor@data[,,i])
  lr_max_core <- c(lr_max_core,max_val)
}
ndx_top <- order(lr_max_core,decreasing=TRUE)[1:200] #indices of top 200 lr pairs


# for each lr pair, select biggest factor interaction from core tensor
sig_channels <- c()
for (i in ndx_top) {
  pair_current <- lr_pairs_in_tensor[i]
  
  lig <- strsplit(pair_current,split='_')[[1]][[1]]
  rec <- strsplit(pair_current,split='_')[[1]][[2]]
  
  # convert lig and rec back to gene symbols
  ndx_gene <- which(gene_conv[,'NCBI gene (formerly Entrezgene) ID']==lig)
  lig <- gene_conv[ndx_gene,'HGNC symbol']
  ndx_gene <- which(gene_conv[,'NCBI gene (formerly Entrezgene) ID']==rec)
  rec <- gene_conv[ndx_gene,'HGNC symbol']
  
  ndx_best <- which(core_tensor@data[,,i]==max(core_tensor@data[,,i]),arr.ind = TRUE) # indicates lig_factor, rec_factor
  
  ct_ndx <- which(l_ct_f[ndx_best[1],] != 0)
  lig_ct <- colnames(l_ct_f)[ct_ndx]
  
  ct_ndx <- which(r_ct_f[ndx_best[2],] != 0)
  rec_ct <- colnames(r_ct_f)[ct_ndx]
  
  # get combinations of ligands and receptors
  ct_combos <- expand.grid(lig_ct, rec_ct)
  ct_combos[,1] <- as.character(ct_combos[,1])
  ct_combos[,2] <- as.character(ct_combos[,2])
  ct_combos <- ct_combos[ct_combos[,1] != ct_combos[,2],]
  
  for (ct_ndx in 1:nrow(ct_combos)) {
    ct_combo <- ct_combos[ct_ndx,]
    source_ct <- as.character(ct_combo[1])
    rec_ct <- as.character(ct_combo[2])
    cur_channel <- paste(lig,source_ct,rec,rec_ct,sep = '_')
    sig_channels <- c(sig_channels,cur_channel)
  }
}

sig_channels <- unique(sig_channels)
length(sig_channels)



### now choosing the same number of top scITD channels
N_select <- length(sig_channels)

mat <- pbmc_container$lr_res
mat_df <- as.data.frame(as.table(mat))
colnames(mat_df) <- c("Row", "Column", "Value")

# Sort the data frame by 'Value' column to get smallest elements
sorted_df <- mat_df[order(mat_df$Value), ]

# Select top N smallest elements
top_N_smallest <- sorted_df[1:N_select, ]
top_N_smallest$Column <- as.character(top_N_smallest$Column)
top_N_smallest$Row <- as.character(top_N_smallest$Row)
top_N_smallest$rec_ct <- sapply(top_N_smallest$Column,function(x){
  strsplit(x,split='_')[[1]][[1]]
})
scITD_top_channels <- paste0(top_N_smallest$Row,'_',top_N_smallest$rec_ct)

# calculate jaccard overlap
jacc_scTensor <- length(intersect(scITD_top_channels,sig_channels)) / length(union(scITD_top_channels,sig_channels))

# get jaccard from random channels
all_ct <- c("B","NK","Th","Tc","cDC",
            "cMono","ncMono")
ct_combos <- expand.grid(all_ct, all_ct)
ct_combos[,1] <- as.character(ct_combos[,1])
ct_combos[,2] <- as.character(ct_combos[,2])
ct_combos <- ct_combos[ct_combos[,1] != ct_combos[,2],]
total_dat_sub_channels <- c()
for (i in 1:nrow(lr_pairs3)) {
  for (j in 1:nrow(ct_combos)) {
    cur_channel <- paste(lr_pairs3$V1[i],ct_combos[j,1],lr_pairs3$V2[i],ct_combos[j,2],sep = '_')
    total_dat_sub_channels <- c(total_dat_sub_channels,cur_channel)
  }
}

total_dat_sub_channels <- unique(total_dat_sub_channels)
print(length(total_dat_sub_channels))

rj_all_scTensor <- sapply(1:10000,function(i) {
  ch_samp <- sample(total_dat_sub_channels,N_select)
  in_both <- intersect(scITD_top_channels,ch_samp)
  in_either <- union(scITD_top_channels,ch_samp)
  rand_jaccard <- length(in_both)/length(in_either)
})

p_enr_scTensor <- sum(rj_all_scTensor>jacc_scTensor)/length(rj_all_scTensor) # pval=1e-04

res_lst <- c(mean(rj_all_scTensor),sd(rj_all_scTensor),jacc_scTensor)
# saveRDS(res_lst,file='/home/jmitchel/data/lupus_data/scTensor_stats.rds')

















