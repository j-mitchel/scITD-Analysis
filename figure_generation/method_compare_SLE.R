
library(MASS)
library(Seurat)
library(DIALOGUE)
library(MOFAcellulaR)
library(dplyr)
library(tidyr)
library(pROC)
library(PRROC)
library(ComplexHeatmap)
library(readxl)
library(devtools)
load_all('/home/jmitchel/scITD')
library(reticulate)
reticulate::use_condaenv("sandbox", required=TRUE)


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

param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc","cDC",
                                               "cMono","ncMono"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(seurat_obj=pbmc,
                                     params=param_list,
                                     metadata_cols=c('ind_cov_batch_cov',
                                                     "SLE_status",
                                                     "Status",
                                                     "cg_cov",
                                                     "sex",
                                                     "Age",
                                                     "batch_cov",
                                                     "Processing_Cohort",
                                                     "Ethnicity",
                                                     "ind_cov"),
                                     metadata_col_nm=c('donors',
                                                       'SLE_status',
                                                       'Status',
                                                       'ctypes',
                                                       'sex',
                                                       'Age',
                                                       'pool',
                                                       'processing',
                                                       'Ethnicity',
                                                       'ind_cov'))


pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20),
                                 tucker_type = 'regular', rotation_type = 'hybrid')

# flip sign of F1 so high ISG expression is positive instead of negative (signs are arbitrary)
pbmc_container$tucker_results[[1]][,1] <- pbmc_container$tucker_results[[1]][,1] * -1
pbmc_container$tucker_results[[2]][1,] <- pbmc_container$tucker_results[[2]][1,] * -1
pbmc_container$projection_data[[1]][1,] <- pbmc_container$projection_data[[1]][1,] * -1














##### now running dialogue
ctypes_use = c("B","NK","Th","Tc","cDC",
               "cMono","ncMono")

ct_list <- list()
for (ct in ctypes_use) {
  cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$cg_cov==ct]
  pbmc_sub <- subset(pbmc,cells=cells_keep)
  samples <- pbmc_sub@meta.data[,'ind_cov_batch_cov',drop=TRUE]
  pbmc_sub <- NormalizeData(pbmc_sub)
  tpm <- pbmc_sub@assays$RNA@data
  
  
  # running combat before computing pcs, since it is used before the other approaches as well
  pbmc_sub <- FindVariableFeatures(pbmc_sub, selection.method = "vst", nfeatures = 1000)
  tpm <- tpm[VariableFeatures(object = pbmc_sub),]
  tpm <- as.matrix(tpm)
  
  tpm_scale <- t(scale(t(tpm)))
  combat_res <- sva::ComBat(tpm_scale,pbmc_sub@meta.data$batch_cov)
  pbmc_sub@assays$RNA@scale.data <- combat_res
  pbmc_sub <- RunPCA(pbmc_sub, features = VariableFeatures(object = pbmc_sub))
  
  X <- pbmc_sub@reductions[["pca"]]@cell.embeddings[,c(1:10)]
  metadata <- as.data.frame(pbmc_sub@meta.data[,c("SLE_status",
                                                  "Status",
                                                  "cg_cov",
                                                  "sex",
                                                  "Age",
                                                  "batch_cov",
                                                  "Processing_Cohort",
                                                  "Ethnicity","nCount_RNA")])
  colnames(metadata)[9] <- 'cellQ'
  ct_dat <- make.cell.type(name = ct,tpm,samples,X,metadata,cellQ = metadata$cellQ)
  ct_list[[ct]] <- ct_dat
}

param <- DLG.get.param(k = 7,
                       results.dir = "/home/jmitchel/data/dialogue_results/",
                       conf = c("batch_cov","cellQ"),
                       find.genes = FALSE)

dial_res <- DIALOGUE.run(rA = ct_list, # list of cell.type objects
                         main = "sle.data.v2",
                         param = param,
                         plot.flag = FALSE)


## read in saved dialogue results
# dial_res <- readRDS('/home/jmitchel/data/dialogue_results/DLG.full.output_sle.data.v3.rds')
dial_res <- readRDS('/home/jmitchel/data/dialogue_results/DIALOGUE1_sle.data.v3.rds')


# computing average donor scores across included cell types per MCP
dial_dscores <- list()
for (i in 1:length(dial_res$MCPs)) {
  ct_include <- names(dial_res$MCPs[[i]])
  ct_include <- sapply(ct_include,function(x){
    strsplit(x,split='.',fixed = TRUE)[[1]][[1]]
  })
  ct_include <- unique(ct_include)
  av_scores <- list()
  for (ct in ct_include) {
    ct_scores <- dial_res[["scores"]][[ct]][,c(paste0('MCP',i),'samples')]
    colnames(ct_scores)[1] <- 'MCP'
    if (sum(ct_scores$MCP==0)==nrow(ct_scores)) {
      next
    }
    dscore_means <- ct_scores %>%
      group_by(samples) %>%
      summarize(mean = mean(MCP,na.rm=TRUE))
    dscore_means$mean <- scale(dscore_means$mean)
    av_scores[[ct]] <- dscore_means
  }
  d_intersect <- av_scores[[1]]$samples
  for (ct in names(av_scores)) {
    d_intersect <- intersect(d_intersect,av_scores[[ct]]$samples)
  }
  for (ct in names(av_scores)) {
    dscore_means <- as.data.frame(av_scores[[ct]])
    rownames(dscore_means) <- dscore_means$samples
    dscore_means <- dscore_means[d_intersect,'mean']
    av_scores[[ct]] <- dscore_means
  }
  
  av_av_scores <- as.data.frame(do.call(cbind, av_scores))
  av_av_scores <- rowMeans(av_av_scores)
  names(av_av_scores) <- d_intersect
  dial_dscores[[i]] <- av_av_scores
}
d_intersect <- names(dial_dscores[[1]])
for (i in 1:length(dial_dscores)) {
  d_intersect <- intersect(d_intersect,names(dial_dscores[[i]]))
}
for (i in 1:length(dial_dscores)) {
  dial_dscores[[i]] <- dial_dscores[[i]][d_intersect]
}
dial_dscores <- as.data.frame(do.call(cbind, dial_dscores))


### alternate way of calculating dscores for dialogue
dscores <- lapply(names(dial_res[["sample.PCs"]]), function(i) dial_res[["sample.PCs"]][[i]] %*% dial_res$cca$ws[[i]])
dscores_cross_ct_av <- matrix(0,ncol=ncol(dscores[[1]]),nrow=nrow(dscores[[1]]))
for (i in 1:length(dscores)) {
  dscores_cross_ct_av <- dscores_cross_ct_av + dscores[[i]]
}
dscores_cross_ct_av <- dscores_cross_ct_av / length(dscores)
dial_dscores <- dscores_cross_ct_av



scITD_dsc <- pbmc_container[["tucker_results"]][[1]]
cormat <- abs(cor(scITD_dsc,dial_dscores[rownames(scITD_dsc),])) # scITD dsc are the rows here and mofa's are columns
rownames(cormat) <- paste0('scITD_',1:nrow(cormat))
colnames(cormat) <- paste0('DIALOGUE_',1:ncol(cormat))
cormat <- t(cormat)
col_fun = colorRamp2(c(0, 1), c("white", "red"))
myhmap <- Heatmap(cormat,name = "Pearson r",
                  cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  column_names_gp = gpar(fontsize = 10),
                  row_names_gp = gpar(fontsize = 10),
                  col = col_fun,border=TRUE, show_column_names=TRUE,
                  show_row_names=TRUE,show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  column_names_side = "top",
                  row_names_side = "left",  
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid::grid.text(sprintf("%.2f", cormat[i, j]), x, y, gp = gpar(fontsize = 10))
                  })

myhmap

### Figure 3b middle
pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/scITD_dialogue_dscores.pdf", useDingbats = FALSE,
    width = 8, height = 3.5)
myhmap
dev.off()









##### running mofa
ctypes_use = c("B","NK","Th","Tc","cDC",
               "cMono","ncMono")

# convert seurat object to single cell experiment type
pbmc.sce <- as.SingleCellExperiment(pbmc)

dat_pb <- scuttle::summarizeAssayByGroup(pbmc.sce,pbmc.sce@colData[,c('ind_cov','cg_cov')],statistics='sum')
pb_meta <- as.data.frame(dat_pb@colData@listData)
pb_counts <- dat_pb@assays@data@listData[["sum"]]

colnames(pb_meta) <- c('donor_id','cell_type','cell_counts')
concat_nms <- paste0(pb_meta$cell_type,'_',pb_meta$donor_id)
rownames(pb_meta) <- concat_nms
colnames(pb_counts) <- concat_nms

pb_obj <- MOFAcellulaR::create_init_exp(counts = pb_counts,  coldata = pb_meta)

ct_list <- MOFAcellulaR::filt_profiles(pb_dat = pb_obj,
                                       cts = ctypes_use,
                                       ncells = 20, # don't remove any samples for not having enough cells
                                       counts_col = "cell_counts", # This refers to the column name in testcoldata where the number of cells per profile was stored
                                       ct_col = "cell_type") # This refers to the column name in testcoldata where the cell-type label was stored

### filtering - might need to adjust accordingly
# remove cell types that have too few samples remaining
ct_list <- MOFAcellulaR::filt_views_bysamples(pb_dat_list = ct_list,
                                              nsamples = 2)


ct_list <- MOFAcellulaR::filt_gex_byexpr(pb_dat_list = ct_list,
                                         min.count = 5, # Modify!!
                                         min.prop = 0.25) # Modify!!

ct_list <- filt_views_bygenes(pb_dat_list = ct_list,
                              ngenes = 15)

ct_list <- filt_samples_bycov(pb_dat_list = ct_list,
                              prop_coverage = 0.9)
###

ct_list <- MOFAcellulaR::tmm_trns(pb_dat_list = ct_list,
                                  scale_factor = 1000000)

ct_list <- MOFAcellulaR::filt_gex_byhvg(pb_dat_list = ct_list,
                                        prior_hvg = NULL,
                                        var.threshold = 0)


##### replacing normalized expression with batch corrected expression so can compare factors to scITD
pbmc_container <- parse_data_by_ctypes(pbmc_container)
for (ct in pbmc_container$experiment_params$ctypes_use) {
  ctype_sub <- pbmc_container$scMinimal_ctype[[ct]]
  
  # identify donors with few cells
  donor_counts <- table(ctype_sub$metadata$donors)
  donors_keep <- names(donor_counts)[donor_counts > 1]
  
  # subset on donors
  ctype_sub <- subset_scMinimal(ctype_sub, donors_use = donors_keep)
}
pbmc_container <- get_pseudobulk(pbmc_container)
pbmc_container <- normalize_pseudobulk(pbmc_container, method='trim', scale_factor=1000000)
pbmc_container <- get_normalized_variance(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method='norm_var_pvals', thresh=1)
pbmc_container <- apply_combat(pbmc_container,batch_var='pool')

pbmc_container <- get_donor_meta(pbmc_container,c('ind_cov','pool'),only_analyzed = FALSE)
d_meta <- pbmc_container$donor_metadata

for (ct in ctypes_use) {
  dat_to_replace <- ct_list[[ct]]@assays@data@listData[["logcounts"]]
  dat_replace_with <- pbmc_container$scMinimal_ctype[[ct]]$pseudobulk
  dat_replace_with <- t(dat_replace_with)
  colnames(dat_replace_with) <- d_meta[colnames(dat_replace_with),'ind_cov']
  colnames(dat_replace_with) <- paste0(ct,'_',colnames(dat_replace_with))
  dat_replace_with <- dat_replace_with[rownames(dat_to_replace),colnames(dat_to_replace)]
  ct_list[[ct]]@assays@data@listData[["logcounts"]] <- dat_replace_with
}

## convert the list to a MOFA object
multiview_dat <- pb_dat2MOFA(pb_dat_list = ct_list, 
                             sample_column = "donor_id")


### running the MOFA model
MOFAobject <- MOFA2::create_mofa(multiview_dat)

data_opts <- MOFA2::get_default_data_options(MOFAobject)
train_opts <- MOFA2::get_default_training_options(MOFAobject)
model_opts <- MOFA2::get_default_model_options(MOFAobject)

# This avoids the regularization of multicellular programs per cell type.
# This avoids less sparse gene weights
model_opts$spikeslab_weights <- FALSE 

# Define the number of factors needed
model_opts$num_factors <- 7

# Prepare MOFA model:
MOFAobject <- MOFA2::prepare_mofa(object = MOFAobject,
                                  data_options = data_opts,
                                  model_options = model_opts,
                                  training_options = train_opts)

mofa_model <- MOFA2::run_mofa(MOFAobject, outfile=NULL)




## plotting the donor scores correlation from MOFA to scITD
d_meta <- unique(pb_meta[,'donor_id',drop=FALSE])
rownames(d_meta) <- d_meta$donor_id
colnames(d_meta)[1] <- c('sample')
rownames(d_meta) <- NULL
mofa_dsc <- MOFAcellulaR::get_tidy_factors(model = mofa_model,
                                           metadata = d_meta,
                                           factor = 'all',
                                           sample_id_column = "sample")

# need to get donor id conversions from long to short names
pbmc_container <- get_donor_meta(pbmc_container,'ind_cov',only_analyzed = FALSE)
nm_conv <- as.character(pbmc_container$donor_metadata[,'ind_cov'])
names(nm_conv) <- rownames(pbmc_container$donor_metadata)
nm_conv <- setNames(names(nm_conv), nm_conv)

mofa_dsc_melt <- mofa_dsc %>%
  pivot_wider(names_from = Factor, values_from = value)

mofa_dsc_melt <- as.data.frame(mofa_dsc_melt)

rownames(mofa_dsc_melt) <- mofa_dsc_melt$sample
mofa_dsc_melt$sample <- NULL

scITD_dsc <- pbmc_container[["tucker_results"]][[1]]
rownames(mofa_dsc_melt) <- nm_conv[rownames(mofa_dsc_melt)]
cormat <- abs(cor(scITD_dsc,mofa_dsc_melt[rownames(scITD_dsc),])) # scITD dsc are the rows here and mofa's are columns
rownames(cormat) <- paste0('scITD_',1:nrow(cormat))
colnames(cormat) <- paste0('MOFA_',1:ncol(cormat))
cormat <- t(cormat)

col_fun = colorRamp2(c(0, 1), c("white", "red"))
myhmap <- Heatmap(cormat,name = "Pearson r",
                  cluster_columns = FALSE,
                  cluster_rows = FALSE,
                  column_names_gp = gpar(fontsize = 10),
                  row_names_gp = gpar(fontsize = 10),
                  col = col_fun,border=TRUE, show_column_names=TRUE,
                  show_row_names=TRUE,show_row_dend = FALSE,
                  show_column_dend = FALSE,
                  column_names_side = "top",
                  row_names_side = "left",  
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid::grid.text(sprintf("%.2f", cormat[i, j]), x, y, gp = gpar(fontsize = 10))
                  })

### Figure 3b top
pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/scITD_mofa_dscores.pdf", useDingbats = FALSE,
    width = 8, height = 3.5)
myhmap
dev.off()




### now comparing scITD with PCA
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')
scITD_dsc <- pbmc_container[["tucker_results"]][[1]]

pbmc_container <- pca_unfolded(pbmc_container, 7)
pbmc_container <- get_lm_pvals(pbmc_container)
PCA_dsc <- pbmc_container[["tucker_results"]][[1]]

cormat <- abs(cor(scITD_dsc,PCA_dsc[rownames(scITD_dsc),])) # scITD dsc are the rows here and mofa's are columns
rownames(cormat) <- paste0('scITD_',1:nrow(cormat))
colnames(cormat) <- paste0('PCA_',1:ncol(cormat))

col_fun = colorRamp2(c(0, 1), c("white", "red"))
Heatmap(cormat,name = "Pearson r",
        column_title = 'scITD-PCA donor scores correlations',
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        col = col_fun,border=TRUE, show_column_names=TRUE,
        show_row_names=TRUE,show_row_dend = FALSE,
        show_column_dend = FALSE,
        column_title_side = "top",
        row_names_side = "left",  
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid::grid.text(sprintf("%.2f", cormat[i, j]), x, y, gp = gpar(fontsize = 10))
        })




### testing applying ica rotation to pca to see how close that makes it to scITD
PCA_dsc_ica <- ica::icafast(PCA_dsc,ncol(PCA_dsc),center=FALSE,alg='def')$S #ica
PCA_dsc_ica <- stats::varimax(PCA_dsc,eps = 1e-15)$loadings #varimax

# trying to get rotation matrix using loadings and counter rotating dscores
PCA_lds <- pbmc_container[["tucker_results"]][[2]]
vari_res <- stats::varimax(PCA_lds,eps = 1e-15)$loadings #varimax


cormat <- abs(cor(scITD_dsc,PCA_dsc_ica[rownames(scITD_dsc),])) # scITD dsc are the rows here and mofa's are columns
rownames(cormat) <- paste0('scITD_',1:nrow(cormat))
colnames(cormat) <- paste0('PCA_',1:ncol(cormat))
cormat <- t(cormat)
col_fun = colorRamp2(c(0, 1), c("white", "red"))
Heatmap(cormat,name = "Pearson r",
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        col = col_fun,border=TRUE, show_column_names=TRUE,
        show_row_names=TRUE,show_row_dend = FALSE,
        show_column_dend = FALSE,
        column_names_side = "top",
        row_names_side = "left",  
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid::grid.text(sprintf("%.2f", cormat[i, j]), x, y, gp = gpar(fontsize = 10))
        })

pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/scITD_pca_dscores.pdf", useDingbats = FALSE,
    width = 8, height = 3.5)
myhmap
dev.off()




# ## saving dscore matrices
# saveRDS(list(scITD_dsc,dial_dscores,mofa_dsc_melt,PCA_dsc),'/home/jmitchel/data/scITD_sim_res/SLE_dsc_compare.rds')







##### clinical association tests
dscore_dat <- readRDS('/home/jmitchel/data/scITD_sim_res/SLE_dsc_compare.rds')

# reload the main sle dataset to get additional metadata
pbmc_container <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_container_w_decomp.rds')

scITD_dsc <- dscore_dat[[1]]
scITD_dsc[,1] <- scITD_dsc[,1] * (-1)
dial_dscores <- dscore_dat[[2]]
mofa_dsc_melt <- dscore_dat[[3]]
PCA_dsc <- dscore_dat[[4]]

clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_categorical.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL

meds <- read_excel('/home/jmitchel/data/lupus_data/SLE_meds_cleaned.xlsx')
meds <- as.data.frame(meds)
rownames(meds) <- meds[,'Sample ID']
meds[,'Sample ID'] <- NULL

# make all NA into zeros, since no 0 are put in the table
meds[is.na(meds)] <- 0

# getting age, ethnicity, and sex covariates
pbmc_container <- get_donor_meta(pbmc_container,additional_meta = c('Age','sex','Ethnicity'),only_analyzed = FALSE)
d_covars <- pbmc_container$donor_metadata
d_covars$donors <- NULL

# trim donor IDs
trim_names <- sapply(rownames(scITD_dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
rownames(scITD_dsc) <- trim_names

trim_names <- sapply(rownames(PCA_dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
rownames(PCA_dsc) <- trim_names

trim_names <- sapply(rownames(dial_dscores), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
rownames(dial_dscores) <- trim_names

trim_names <- sapply(rownames(mofa_dsc_melt),function(x){
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
rownames(mofa_dsc_melt) <- trim_names

trim_names <- sapply(rownames(d_covars), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
rownames(d_covars) <- trim_names

dscore_dat <- list(scITD_dsc,dial_dscores,mofa_dsc_melt,PCA_dsc)
names(dscore_dat) <- c('scITD','dialogue','mofa','PCA')


get_assoc_auc <- function(var_dat_use,dat,factor_use,var_use) {
  d_both <- intersect(rownames(dat),rownames(var_dat_use))
  tmp <- cbind.data.frame(dat[d_both,factor_use],var_dat_use[d_both,var_use])
  colnames(tmp) <- c('dsc','myvar')
  if (sum(unique(tmp$myvar) %in% c(0,1))!=2) {
    tmp$myvar <- factor(tmp$myvar,levels=unique(tmp$myvar))
    levels(tmp$myvar) <- c(1,0)
  }
  pROC <- roc(tmp$myvar,tmp$dsc,
              smoothed = TRUE,
              plot=FALSE, AUC=TRUE)
  auc <- pROC[["auc"]]
  
  ## getting bootstrap se for the auc here
  boot_auc_all <- c()
  for (i in 1:1000) {
    tmp2 <- tmp[sample(1:nrow(tmp),replace = TRUE),]
    if (all(tmp2$myvar==tmp2$myvar[1])) {
      next
    }
    pROC <- roc(tmp2$myvar,tmp2$dsc,
                smoothed = TRUE,
                plot=FALSE, AUC=TRUE)
    boot_auc <- pROC[["auc"]]
    boot_auc_all <- c(boot_auc_all,boot_auc)
  }
  auc_se <- sd(boot_auc_all)
  return(list(auc,auc_se))
}

get_assoc_neph_abs_auc <- function(var_dat_use,dat,factor_use) {
  d_both <- intersect(rownames(dat),rownames(var_dat_use))
  d_both <- intersect(d_both,rownames(scITD_dsc))
  d_both <- intersect(d_both,rownames(meds))
  tmp <- cbind.data.frame(dat[d_both,factor_use],scITD_dsc[d_both,1],meds[d_both,'prednisone'])
  colnames(tmp) <- c('dsc','dsc_f1','prednisone')
  tmp <- cbind.data.frame(tmp,var_dat_use[d_both,'crflupusneph'])
  colnames(tmp)[4] <- 'neph'
  
  # subset to only donors with high IFN factor
  tmp <- tmp[tmp$dsc_f1>0,]
  
  tmp$neph <- as.factor(tmp$neph)
  
  # now compute AUC separately for donors on/off prednisone and average them
  pROC <- roc(tmp$neph,tmp$dsc,
              smoothed = TRUE,
              plot=FALSE, AUC=TRUE)
  auc <- pROC[["auc"]]
  
  ## getting bootstrap se for the auc here
  boot_auc_all <- c()
  for (i in 1:1000) {
    tmp_resamp <- tmp[sample(1:nrow(tmp),replace = TRUE),]
    
    pROC <- roc(tmp_resamp$neph,tmp_resamp$dsc,
                smoothed = TRUE,
                plot=FALSE, AUC=TRUE)
    boot_auc <- pROC[["auc"]]
    
    boot_auc_all <- c(boot_auc_all,boot_auc)
  }
  auc_se <- sd(boot_auc_all)
  
  return(list(auc,auc_se))
}

get_assoc_rsq <- function(var_dat_use,dat,factor_use,var_use) {
  ### comparing strength of the f1-autoantibody association
  d_both <- intersect(rownames(dat),rownames(var_dat_use))
  tmp <- cbind.data.frame(dat[d_both,factor_use],var_dat_use[d_both,var_use])
  colnames(tmp) <- c('dsc','myvar')
  lmres <- summary(lm(myvar~dsc,data=tmp))
  rsq <- lmres$r.squared
  
  ## getting bootstrap se for the auc here
  boot_rsq_all <- c()
  for (i in 1:1000) {
    tmp2 <- tmp[sample(1:nrow(tmp),replace = TRUE),]
    lmres <- summary(lm(myvar~dsc,data=tmp2))
    boot_rsq <- lmres$r.squared
    boot_rsq_all <- c(boot_rsq_all,boot_rsq)
  }
  rsq_se <- sd(boot_rsq_all)
  
  return(list(rsq,rsq_se))
}


### making it into general fn
get_res_to_plot <- function(var_dat_use,var_to_test,scITD_factor_use,expected_best) {
  ## first get results for anti dsdna associations
  res_plot <- data.frame(matrix(ncol=5,nrow=length(dscore_dat)))
  rownames(res_plot) <- names(dscore_dat)
  colnames(res_plot) <- c('AUC','AUC_se','method','best_factor','var_tested')
  ## computing the results for scITD
  if (var_to_test=='crflupusneph') {
    auc_scITD <- get_assoc_neph_abs_auc(var_dat_use,scITD_dsc,scITD_factor_use)
  } else if (var_to_test=='Age') {
    auc_scITD <- get_assoc_rsq(var_dat_use,scITD_dsc,scITD_factor_use,var_to_test)
  } else {
    auc_scITD <- get_assoc_auc(var_dat_use,scITD_dsc,scITD_factor_use,var_to_test)
  }
  res_plot['scITD',] <- c(auc_scITD[[1]],auc_scITD[[2]],'scITD',scITD_factor_use,var_to_test)
  for (method_ndx in 2:length(dscore_dat)) {
    method_dsc <- dscore_dat[[method_ndx]]
    method_nm <- names(dscore_dat)[method_ndx]
    method_all_auc <- c()
    method_all_auc_se <- c()
    for (i in 1:7) {
      if (var_to_test=='crflupusneph') {
        auc_res <- get_assoc_neph_abs_auc(var_dat_use,method_dsc,i)
      } else if (var_to_test=='Age') {
        auc_res <- get_assoc_rsq(var_dat_use,method_dsc,i,var_to_test)
      } else {
        auc_res <- get_assoc_auc(var_dat_use,method_dsc,i,var_to_test)
      }
      
      method_all_auc <- c(method_all_auc,auc_res[[1]])
      method_all_auc_se <- c(method_all_auc_se,auc_res[[2]])
    }
    ndx_best <- which(method_all_auc==max(method_all_auc))
    method_max_auc <- method_all_auc[ndx_best]
    method_max_auc_se <- method_all_auc_se[ndx_best]
    res_plot[method_nm,'AUC'] <- method_max_auc
    res_plot[method_nm,'AUC_se'] <- method_max_auc_se
    res_plot[method_nm,'method'] <- method_nm
    res_plot[method_nm,'best_factor'] <- ndx_best
    res_plot[method_nm,'var_tested'] <- var_to_test
  }
  
  # expected matching factors for scITD, dial, mofa, pca in that order
  match_exp <- c()
  for (i in 1:nrow(res_plot)) {
    method_matches_exp <- res_plot$best_factor[i]==expected_best[i]
    match_exp <- c(match_exp,method_matches_exp)
  }
  res_plot$match_exp <- match_exp
  return(res_plot)
}

res_plot_dsdna <- get_res_to_plot(clin_vars,'acrantidsdna',1,c(1,1,1,1))
res_plot_neph <- get_res_to_plot(clin_vars,'crflupusneph',2,c(2,3,2,2))
res_plot_pred <- get_res_to_plot(meds,'prednisone',3,c(3,3,4,4))
res_plot_eth <- get_res_to_plot(d_covars,'Ethnicity',5,c(5,4,6,5))
res_plot_s <- get_res_to_plot(d_covars,'sex',6,c(6,7,7,6))
res_plot_age <- get_res_to_plot(d_covars,'Age',4,c(4,5,3,3))

colnames(res_plot_age)[1:2] <- c('rsq','rsq_se')


rownames(res_plot_dsdna) <- NULL
rownames(res_plot_neph) <- NULL
rownames(res_plot_pred) <- NULL
rownames(res_plot_age) <- NULL
rownames(res_plot_eth) <- NULL
rownames(res_plot_s) <- NULL

res_plot_dsdna$AUC <- as.numeric(res_plot_dsdna$AUC)
res_plot_neph$AUC <- as.numeric(res_plot_neph$AUC)
res_plot_pred$AUC <- as.numeric(res_plot_pred$AUC)
res_plot_age$rsq <- as.numeric(res_plot_age$rsq)
res_plot_eth$AUC <- as.numeric(res_plot_eth$AUC)
res_plot_s$AUC <- as.numeric(res_plot_s$AUC)

res_plot_dsdna$AUC_se <- as.numeric(res_plot_dsdna$AUC_se)
res_plot_neph$AUC_se <- as.numeric(res_plot_neph$AUC_se)
res_plot_pred$AUC_se <- as.numeric(res_plot_pred$AUC_se)
res_plot_age$rsq_se <- as.numeric(res_plot_age$rsq_se)
res_plot_eth$AUC_se <- as.numeric(res_plot_eth$AUC_se)
res_plot_s$AUC_se <- as.numeric(res_plot_s$AUC_se)

library(ggpattern)

annotate_match <- function(res_plot) {
  res_plot$match_exp <- sapply(res_plot$match_exp,function(x){
    if (x) {
      return('match')
    } else {
      return('not_match')
    }
  })
  return(res_plot)
}

res_plot_dsdna <- annotate_match(res_plot_dsdna)
res_plot_neph <- annotate_match(res_plot_neph)
res_plot_pred <- annotate_match(res_plot_pred)
res_plot_age <- annotate_match(res_plot_age)
res_plot_eth <- annotate_match(res_plot_eth)
res_plot_s <- annotate_match(res_plot_s)


## removing the pca bars because it's getting too complicated to show everything
res_plot_dsdna <- res_plot_dsdna[1:3,]
res_plot_neph <- res_plot_neph[1:3,]
res_plot_pred <- res_plot_pred[1:3,]
res_plot_age <- res_plot_age[1:3,]
res_plot_eth <- res_plot_eth[1:3,]
res_plot_s <- res_plot_s[1:3,]

res_plot_dsdna$method <- factor(res_plot_dsdna$method,levels=c('scITD','dialogue','mofa'))
res_plot_neph$method <- factor(res_plot_neph$method,levels=c('scITD','dialogue','mofa'))
res_plot_pred$method <- factor(res_plot_pred$method,levels=c('scITD','dialogue','mofa'))
res_plot_age$method <- factor(res_plot_age$method,levels=c('scITD','dialogue','mofa'))
res_plot_eth$method <- factor(res_plot_eth$method,levels=c('scITD','dialogue','mofa'))
res_plot_s$method <- factor(res_plot_s$method,levels=c('scITD','dialogue','mofa'))

res_plot_dsdna$var_tested <- 'antidsdna'
res_plot_neph$var_tested <- 'nephritis'
res_plot_pred$var_tested <- 'prednisone'
res_plot_age$var_tested <- 'age'
res_plot_eth$var_tested <- 'ethnicity'
res_plot_s$var_tested <- 'sex'


p1 <- ggplot(res_plot_dsdna,aes(x=var_tested,y=AUC,fill=method,pattern=match_exp)) +
  geom_bar(stat="identity", position=position_dodge(),width=0.8) +
  geom_bar_pattern(stat="identity",
                   position = position_dodge(),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.001,
                   pattern_spacing = 0.1,
                   pattern_key_scale_factor = 0.6,
                   width=0.8) +
  geom_errorbar(
    aes(ymin = AUC-AUC_se, ymax = AUC+AUC_se), 
    width = 0.2, position = position_dodge(0.8)
  ) +
  scale_pattern_manual(values = c(not_match = "stripe", match = "none")) +
  ylim(0,1) +
  xlab('') +
  ylab('AUC') +
  theme_classic() + 
  theme(legend.position="none")
p1


p2 <- ggplot(res_plot_neph,aes(x=var_tested,y=AUC,fill=method,pattern=match_exp)) +
  geom_bar(stat="identity", position=position_dodge(),width=0.8) +
  geom_bar_pattern(stat="identity",
                   position = position_dodge(),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.001,
                   pattern_spacing = 0.1,
                   pattern_key_scale_factor = 0.6,
                   width=0.8) +
  geom_errorbar(
    aes(ymin = AUC-AUC_se, ymax = AUC+AUC_se), 
    width = 0.2, position = position_dodge(0.8)
  ) +
  scale_pattern_manual(values = c(not_match = "stripe", match = "none")) +
  ylim(0,1) +
  xlab('') +
  ylab('AUC') +
  theme_classic() + 
  theme(legend.position="none")
p2



p3 <- ggplot(res_plot_pred,aes(x=var_tested,y=AUC,fill=method,pattern=match_exp)) +
  geom_bar(stat="identity", position=position_dodge(),width=0.8) +
  geom_bar_pattern(stat="identity",
                   position = position_dodge(),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.001,
                   pattern_spacing = 0.1,
                   pattern_key_scale_factor = 0.6,
                   width=0.8) +
  geom_errorbar(
    aes(ymin = AUC-AUC_se, ymax = AUC+AUC_se), 
    width = 0.2, position = position_dodge(0.8)
  ) +
  scale_pattern_manual(values = c(not_match = "stripe", match = "none")) +
  ylim(0,1) +
  xlab('') +
  ylab('AUC') +
  theme_classic() + 
  theme(legend.position="none")
p3


p4 <- ggplot(res_plot_age,aes(x=var_tested,y=rsq,fill=method,pattern=match_exp)) +
  geom_bar(stat="identity", position=position_dodge(),width=0.8) +
  geom_bar_pattern(stat="identity",
                   position = position_dodge(),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.001,
                   pattern_spacing = 0.1,
                   pattern_key_scale_factor = 0.6,
                   width=0.8) +
  geom_errorbar(
    aes(ymin = rsq-rsq_se, ymax = rsq+rsq_se), 
    width = 0.2, position = position_dodge(0.8)
  ) +
  scale_pattern_manual(values = c(not_match = "stripe", match = "none")) +
  xlab('') +
  ylab('R-squared') +
  theme_classic() + 
  theme(legend.position="none")
p4



p5 <- ggplot(res_plot_eth,aes(x=var_tested,y=AUC,fill=method,pattern=match_exp)) +
  geom_bar(stat="identity", position=position_dodge(),width=0.8) +
  geom_bar_pattern(stat="identity",
                   position = position_dodge(),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.001,
                   pattern_spacing = 0.1,
                   pattern_key_scale_factor = 0.6,
                   width=0.8) +
  geom_errorbar(
    aes(ymin = AUC-AUC_se, ymax = AUC+AUC_se), 
    width = 0.2, position = position_dodge(0.8)
  ) +
  scale_pattern_manual(values = c(not_match = "stripe", match = "none")) +
  ylim(0,1) +
  xlab('') +
  ylab('AUC') +
  theme_classic() + 
  theme(legend.position="none")
p5

p6 <- ggplot(res_plot_s,aes(x=var_tested,y=AUC,fill=method,pattern=match_exp)) +
  geom_bar(stat="identity", position=position_dodge(),width=0.8) +
  geom_bar_pattern(stat="identity",
                   position = position_dodge(),
                   color = "black", 
                   pattern_fill = "black",
                   pattern_angle = 45,
                   pattern_density = 0.001,
                   pattern_spacing = 0.1,
                   pattern_key_scale_factor = 0.6,
                   width=0.8) +
  geom_errorbar(
    aes(ymin = AUC-AUC_se, ymax = AUC+AUC_se), 
    width = 0.2, position = position_dodge(0.8)
  ) +
  scale_pattern_manual(values = c(not_match = "stripe", match = "none")) +
  ylim(0,1) +
  xlab('') +
  ylab('AUC') +
  theme_classic()
p6




library(cowplot)

### Figure 3b bottom
fig <- plot_grid(plotlist=list(p1,p2,p3,p4,p5,p6),nrow=1,rel_widths = c(1,1,1,1,1,1.5))
pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/method_compare_meta6.pdf", useDingbats = FALSE,
    width = 15, height = 5)
fig
dev.off()





