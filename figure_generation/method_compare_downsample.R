.libPaths(c("/home/jmitchel/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))
.libPaths(c("/home/jmitchel/R/x86_64-pc-linux-gnu-library/4.1", .libPaths()))
.libPaths(c("/home/jmitchel/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))

library(MASS)
library(tictoc)
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

### functions to run each tool
# should take seurat object as input and return a matrix of donor scores
run_scITD_single <- function(seurat_obj,use_PCA=FALSE) {
  param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc","cDC",
                                                 "cMono","ncMono"),
                                  ncores = 30, rand_seed = 10)
  
  pbmc_container <- make_new_container(seurat_obj=seurat_obj,
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
  
  tic()
  if (use_PCA) {
    pbmc_container <- pca_unfolded(pbmc_container, 7)
  } else {
    pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20),
                                     tucker_type = 'regular', rotation_type = 'hybrid')
  }
  time_lapse <- toc()
  runtime <- time_lapse[["toc"]][["elapsed"]]-time_lapse[["tic"]][["elapsed"]]
  
  dscores <- pbmc_container$tucker_results[[1]]
  
  return(list(dscores,runtime))
}




run_dialogue_single <- function(seurat_obj) {
  ##### now running dialogue
  ctypes_use = c("B","NK","Th","Tc","cDC",
                 "cMono","ncMono")
  
  ct_list <- list()
  for (ct in ctypes_use) {
    print(ct)
    cells_keep <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$cg_cov==ct]
    pbmc_sub <- subset(seurat_obj,cells=cells_keep)
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
  
  tic()
  dial_res <- DIALOGUE.run(rA = ct_list, # list of cell.type objects
                           main = "sle.data.downsamp",
                           param = param,
                           plot.flag = FALSE)
  
  time_lapse <- toc()
  runtime <- time_lapse[["toc"]][["elapsed"]]-time_lapse[["tic"]][["elapsed"]]
  
  dscores <- lapply(names(dial_res[["sample.PCs"]]), function(i) dial_res[["sample.PCs"]][[i]] %*% dial_res$cca$ws[[i]])
  dscores_cross_ct_av <- matrix(0,ncol=ncol(dscores[[1]]),nrow=nrow(dscores[[1]]))
  for (i in 1:length(dscores)) {
    dscores_cross_ct_av <- dscores_cross_ct_av + dscores[[i]]
  }
  dscores_cross_ct_av <- dscores_cross_ct_av / length(dscores)
  dial_dscores <- dscores_cross_ct_av
  dscores <- dial_dscores
  
  file.remove('/home/jmitchel/data/dialogue_results/sle.data.downsamp.rds')
  unlink("/home/jmitchel/data/dialogue_results/sle.data.downsamp/",recursive=TRUE)
  
  return(list(dscores,runtime))
}


run_MOFA_single <- function(seurat_obj) {
  
  ctypes_use = c("B","NK","Th","Tc","cDC",
                 "cMono","ncMono")
  
  # convert seurat object to single cell experiment type
  pbmc.sce <- as.SingleCellExperiment(seurat_obj)
  
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
  ## using scITD just to get the batch corrected data
  param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc","cDC",
                                                 "cMono","ncMono"),
                                  ncores = 30, rand_seed = 10)
  
  pbmc_container <- make_new_container(seurat_obj=seurat_obj,
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
  
  tic()
  mofa_model <- MOFA2::run_mofa(MOFAobject, outfile=NULL)
  time_lapse <- toc()
  runtime <- time_lapse[["toc"]][["elapsed"]]-time_lapse[["tic"]][["elapsed"]]
  
  
  
  ## plotting the donor scores correlatino from MOFA to scITD
  d_meta <- unique(pb_meta[,'donor_id',drop=FALSE])
  rownames(d_meta) <- d_meta$donor_id
  colnames(d_meta)[1] <- c('sample')
  rownames(d_meta) <- NULL
  mofa_dsc <- MOFAcellulaR::get_tidy_factors(model = mofa_model,
                                             metadata = d_meta,
                                             factor = 'all',
                                             sample_id_column = "sample")
  
  pbmc_container <- get_donor_meta(pbmc_container,'ind_cov',only_analyzed = FALSE)
  nm_conv <- as.character(pbmc_container$donor_metadata[,'ind_cov'])
  names(nm_conv) <- rownames(pbmc_container$donor_metadata)
  nm_conv <- setNames(names(nm_conv), nm_conv)
  
  mofa_dsc_melt <- mofa_dsc %>%
    pivot_wider(names_from = Factor, values_from = value)
  
  mofa_dsc_melt <- as.data.frame(mofa_dsc_melt)
  
  rownames(mofa_dsc_melt) <- mofa_dsc_melt$sample
  mofa_dsc_melt$sample <- NULL
  
  rownames(mofa_dsc_melt) <- nm_conv[rownames(mofa_dsc_melt)]
  
  dscores <- mofa_dsc_melt
  return(list(dscores,runtime))
}


get_max_cors <- function(orig_dsc,new_dsc) {
  all_max_cors <- c()
  d_both <- intersect(rownames(orig_dsc),rownames(new_dsc))
  for (i in 1:ncol(orig_dsc)) {
    all_cors <- c()
    for (j in 1:ncol(new_dsc)) {
      corval <- abs(cor(orig_dsc[d_both,i],new_dsc[d_both,j]))
      all_cors <- c(all_cors,corval)
    }
    max_cor <- max(all_cors)
    all_max_cors <- c(all_max_cors,max_cor)
  }
  return(all_max_cors)
}





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





# load up dscores generated from full data
dsc_all <- readRDS('/home/jmitchel/data/scITD_sim_res/SLE_dsc_compare.rds')
scITD_dsc_full <- dsc_all[[1]]
dial_dsc_full <- dsc_all[[2]]
mofa_dsc_full <- dsc_all[[3]]
pca_dsc_full <- dsc_all[[4]]

# now iteratively downsample and rerun the tools
n_iter <- 10
rand_seeds <- sample(1:100000000,n_iter, replace=FALSE)
donors_all <- as.character(unique(pbmc@meta.data$ind_cov_batch_cov))
res_all <- data.frame(matrix(ncol=3,nrow=0))
runtime_res <- data.frame(matrix(ncol=2,nrow=0))
for (i in 1:n_iter) {
  print(i)
  # downsample data by donors
  dsamp_seed <- rand_seeds[i]
  set.seed(dsamp_seed)
  d_keep <- sample(donors_all, size=round(length(donors_all)*.85),replace = FALSE)
  cells_keep <- rownames(pbmc@meta.data)[as.character(pbmc@meta.data$ind_cov_batch_cov) %in% d_keep]
  pbmc_down <- subset(pbmc,cells=cells_keep)
  
  scITD_dsc <- run_scITD_single(pbmc_down,use_PCA=FALSE)
  dial_dsc <- run_dialogue_single(pbmc_down)
  mofa_dsc <- run_MOFA_single(pbmc_down)
  pca_dsc <- run_scITD_single(pbmc_down,use_PCA=TRUE)
  
  scITD_runtime <- scITD_dsc[[2]]
  dial_runtime <- dial_dsc[[2]]
  mofa_runtime <- mofa_dsc[[2]]
  pca_runtime <- pca_dsc[[2]]
  
  scITD_dsc <- scITD_dsc[[1]]
  dial_dsc <- dial_dsc[[1]]
  mofa_dsc <- mofa_dsc[[1]]
  pca_dsc <- pca_dsc[[1]]
  
  # get max correlations for each original factor to all new factors
  scITD_max_cors <- get_max_cors(scITD_dsc_full,scITD_dsc)
  print(scITD_max_cors)
  dial_max_cors <- get_max_cors(dial_dsc_full,dial_dsc)
  mofa_max_cors <- get_max_cors(mofa_dsc_full,mofa_dsc)
  pca_max_cors <- get_max_cors(pca_dsc_full,pca_dsc)

  # store result
  max_cor_all_methods <- c(scITD_max_cors,dial_max_cors,mofa_max_cors,pca_max_cors)
  method_indicator <- c(rep('scITD',7),rep('dialogue',7),rep('mofa',7),rep('PCA',7))
  factor_indicator <- rep(c(1:7),4)
  res_all <- rbind.data.frame(res_all,cbind.data.frame(max_cor_all_methods,method_indicator,factor_indicator))

  method_runtimes <- c(scITD_runtime,dial_runtime,mofa_runtime,pca_runtime)
  method_indicator <- c('scITD','dialogue','mofa','PCA')
  runtime_res <- rbind.data.frame(runtime_res,cbind.data.frame(method_runtimes,method_indicator))
}

# saveRDS(list(res_all,runtime_res),file='/home/jmitchel/data/scITD_sim_res/SLE_downsamp_compare.rds') # accidentally downsampled 95% of donors - also mistakenly resampled same donors
# saveRDS(list(res_all,runtime_res),file='/home/jmitchel/data/scITD_sim_res/SLE_downsamp_compare2.rds') # downsampling 85% of donors - actually still had resampling problem
# saveRDS(list(res_all,runtime_res),file='/home/jmitchel/data/scITD_sim_res/SLE_downsamp_compare3.rds') # downsampling 85% of donors - fixed resampling problem




### plotting results
library(ggplot2)
downsamp_res <- readRDS(file='/home/jmitchel/data/scITD_sim_res/SLE_downsamp_compare3.rds')
res_all <- downsamp_res[[1]]
runtime_res <- downsamp_res[[2]]

res_all$factor_indicator <- as.factor(res_all$factor_indicator)
colnames(res_all)[2] <- 'method'
p1 <- ggplot(res_all,aes(x=factor_indicator,y=max_cor_all_methods,color=method)) +
  geom_boxplot() +
  xlab('Factor from full decomposition') +
  ylab('Max cor (to downsampled)') +
  theme_bw()
p1

### Figure 3c
pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/method_compare_downsample.pdf", useDingbats = FALSE,
    width = 9, height = 3)
p1
dev.off()




colnames(runtime_res)[2] <- 'method'
p2 <- ggplot(runtime_res,aes(x=method,y=method_runtimes)) +
  geom_boxplot() +
  xlab('') +
  ylab('Runtime (seconds)') +
  theme_bw()
p2

runtime_res_no_dial <- runtime_res[runtime_res$method!='dialogue',]
p3 <- ggplot(runtime_res_no_dial,aes(x=method,y=method_runtimes)) +
  geom_boxplot() +
  xlab('') +
  ylab('Runtime (seconds)') +
  theme_bw()
p3

### Figure 3d left
pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/method_compare_runtime1.pdf", useDingbats = FALSE,
    width = 4, height = 3.5)
p2
dev.off()

### Figure 3d right
pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/method_compare_runtime2.pdf", useDingbats = FALSE,
    width = 2.5, height = 2)
p3
dev.off()












