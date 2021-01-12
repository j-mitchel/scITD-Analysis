
library(scITD)
library(edgeR)
library(Seurat)
library(SeuratDisk)

### for loading data from h5ad file...
# or read the h5seurat file
pbmc <- LoadH5Seurat("/home/jmitchel/data/lupus_data/Lupus_study_adjusted.h5seurat")

# or load up the subsetted dataset
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v2.rds')

# Get regular normalization of the data if using h5
pbmc <- NormalizeData(pbmc)

# convert batch to factor for meta association analysis
pbmc@meta.data$batch_cov <- factor(pbmc@meta.data$batch_cov,levels=unique(pbmc@meta.data$batch_cov))

# prep the container
pbmc_scMinimal <- seurat_to_scMinimal(pbmc,normalize_counts=FALSE,
                                      metadata_cols=c('ind_cov_batch_cov',
                                                      "Genotype.ID",
                                                      "SLE_status",
                                                      "Status",
                                                      "cg_cov",
                                                      "sex",
                                                      "age",
                                                      "batch_cov",
                                                      "Processing_Cohort"),
                                      metadata_col_nm=c('donors',
                                                        'donor_genotype',
                                                        'SLE_status',
                                                        'Status',
                                                        'ctypes',
                                                        'sex',
                                                        'age',
                                                        'pool',
                                                        'processing'))


# get median number cells per donor
ctypes <- levels(pbmc_scMinimal$metadata$ctypes)
for (ctype in ctypes) {
  tmp <- pbmc_scMinimal$metadata[pbmc_scMinimal$metadata$ctypes==ctype,]
  print(ctype)
  print(median(table(tmp$donors)))
}

# remove original pbmc object because it takes up a lot of RAM
rm(pbmc)
gc()


# load the cleaned data
pbmc_norm_counts <- readRDS('/home/jmitchel/data/lupus_data/lupus_counts_trim_mean.rds')
pbmc_norm_counts <- readRDS('/home/jmitchel/data/lupus_data/lupus_counts_regular_norm.rds')
pbmc_meta <- readRDS('/home/jmitchel/data/lupus_data/lupus_meta_select_vars.rds')

# trying to remove the one batch with only one donor because I think thats causing batch correct issues
ndx_keep <- which(pbmc_meta$pool!='dmx_flare2')
pbmc_norm_counts <- pbmc_norm_counts[,ndx_keep]
pbmc_meta <- pbmc_meta[ndx_keep,]
pbmc_meta$pool <- factor(pbmc_meta$pool,levels=unique(pbmc_meta$pool))

pbmc_scMinimal <- instantiate_scMinimal(data_sparse=pbmc_norm_counts, meta_data=pbmc_meta)

# pbmc_container <- make_new_container(pbmc_scMinimal,
#                                      ctypes_use = c("B","NK","T4","T8","cDC",
#                                                     "cM","ncM"),
#                                      scale_var = TRUE,
#                                      var_scale_power = 1.25,
#                                      tucker_type = 'sparse', rotation_type = 'ica',
#                                      ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = c("NK","T4","T8","cM"),
                                     scale_var = TRUE,
                                     var_scale_power = 1.5,
                                     tucker_type = 'regular', rotation_type = 'ica',
                                     ncores = 30, rand_seed = 10)

rm(pbmc_scMinimal)
gc()


# finish prepping the data
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=500)

# trying with my most recent batch correct method
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=1000)
# pbmc_container <- apply_cellwise_batch_correct(pbmc_container,'pool')

# trying mixed model for batch correct
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="anova", thresh=.05)

# trying deseq for variable gene selection
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="deseq", thresh=.001)
saved_vargenes<- pbmc_container$all_vargenes
pbmc_container$all_vargenes <- saved_vargenes #done after recomputing everything
pbmc_container <- reduce_to_vargenes(pbmc_container)

# apply batch correction
for (ct in pbmc_container$experiment_params$ctypes_use) {
  scMinimal <- pbmc_container$scMinimal_ctype[[ct]]
  
  metadata <- scMinimal$metadata
  modcombat <- model.matrix(~1, data=metadata)
  tmp <- sva::ComBat(dat=scMinimal$data_sparse,
                     batch=metadata[,'pool'],
                     mod=modcombat, par.prior=TRUE,
                     prior.plots=FALSE)
  

  tmp <- Matrix(tmp, sparse = TRUE)
  
  # need to replace both the metadata and the data_sparse since this field is used to get vargenes
  scMinimal$data_sparse <- tmp
}
# now need to rerun everything with normalized data and manually limit to these genes


# testing out batch correction on collapsed data before var gene selection
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- collapse_by_donors(pbmc_container, shuffle=FALSE)
pbmc_container <- apply_pseudobulk_batch_correct(pbmc_container,'pool')
pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=1500) #make necessary changes
## probably need more changes to make it work with run tucker ica because this recalls collapse by donors
## maybe need to have the apply pseudobulk fn put the data in data sparse. I think I do something similar with
## the gtex data so take a look at that though that may have been relating to vargenes fn only

# how many in saved vargenes not in new vargenes
sum(saved_vargenes %in% pbmc_container$all_vargenes)


# trying the above method but without batch correction first...
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- collapse_by_donors(pbmc_container, shuffle=FALSE)
# need to set collapsed matrices as data_sparse matrices
for (ct in pbmc_container$experiment_params$ctypes_use) {
  scMinimal <- pbmc_container$scMinimal_ctype[[ct]]
  
  # need metadata at donor level
  metadata <- unique(scMinimal$metadata)
  rownames(metadata) <- metadata$donors
  metadata <- metadata[rownames(scMinimal$data_means),]
  
  # need to replace both the metadata and the data_sparse since this field is used to get vargenes
  scMinimal$data_sparse <- Matrix(t(scMinimal$data_means), sparse = TRUE)
  scMinimal$metadata <- metadata
}
pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=1500) #make necessary changes
# now applying batch correct after selecting vargenes
for (ct in pbmc_container$experiment_params$ctypes_use) {
  scMinimal <- pbmc_container$scMinimal_ctype[[ct]]
  metadata <- scMinimal$metadata
  modcombat <- model.matrix(~1, data=metadata)
  tmp <- sva::ComBat(dat=scMinimal$data_sparse,
                     batch=metadata[,'pool'],
                     mod=modcombat, par.prior=TRUE,
                     prior.plots=FALSE)
  

  tmp <- Matrix(tmp, sparse = TRUE)
  
  # need to replace both the metadata and the data_sparse since this field is used to get vargenes
  scMinimal$data_sparse <- tmp
}
saved_vargenes <- pbmc_container$all_vargenes


# determine appropriate variance scaling parameter
pbmc_container <- optimize_var_scale_power(pbmc_container, min_ranks_test=c(5,8,5),
                                           max_ranks_test=c(10,15,5),
                                           min_power_test=1.5,
                                           max_power_test=2.25)
pbmc_container$plots$var_scale_plot

# determine appropriate ranks to use for decomposition
pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(10,18,7),
                                         method='svd', num_iter=5, shuffle_level='cells')
pbmc_container$plots$rank_determination_plot

# run tucker
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(15,25,4), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,10,4), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(9,20,4), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,10,4), batch_var='pool')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(12,16,4), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(12,24,4), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(14,20,4), shuffle=FALSE)

# get metadata associations for donor scores plot
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','pool'))

# plot donor scores first by clustering by sex
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status','pool','processing'),
                                    cluster_by_meta='pool', show_donor_ids = FALSE,add_meta_associations=T)
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status','pool','processing'),
                                    show_donor_ids = FALSE, add_meta_associations=T)
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('Status'),
                                    cluster_by_meta='Status', show_donor_ids = FALSE,
                                    add_meta_associations=T)
pbmc_container <- plot_donor_matrix(pbmc_container,show_donor_ids = FALSE, add_meta_associations=T)
pbmc_container$plots$donor_matrix

# compare do decomp with batch effect
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(15,25,4), shuffle=FALSE)
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','pool'))
tucker_res1 <- pbmc_container$tucker_results
meta_anno1 <- pbmc_container$meta_associations
# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(8,16,4), shuffle=FALSE, batch_var='pool')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(12,24,4), shuffle=FALSE, batch_var='pool')
# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(10,20,4), shuffle=FALSE, batch_var='pool')
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','pool'))
tucker_res2 <- pbmc_container$tucker_results
meta_anno2 <- pbmc_container$meta_associations
decomp_names <- c('no combat','with combat')
myhmap <- compare_decompositions(tucker_res1,tucker_res2,decomp_names,meta_anno1,meta_anno2,comparison='dscores')
myhmap

# plot meta associations
pbmc_container <- plot_meta_associations(pbmc_container,vars_test=c('sex','pool'),anno_side='top')
pbmc_container <- plot_meta_associations(pbmc_container,vars_test=c('pool'))
draw(pbmc_container$plots$meta_associations)

# plot donor scores with associations together
pbmc_container <- plot_donor_matrix(pbmc_container,show_donor_ids = FALSE,add_meta_associations=TRUE)

# plot donor scores with row clustering
pbmc_container <- plot_donor_matrix(pbmc_container, show_donor_ids = FALSE)
pbmc_container$plots$donor_matrix

# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, n_fibers=250, n_iter=500)

# show all loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=FALSE, 
                                           nonsig_to_zero=FALSE,
                                           display_genes=FALSE,
                                           gene_callouts=FALSE)
render_all_lds_plots(pbmc_container, n_rows=2)

# show all loadings plots sig genes only
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE, 
                                           nonsig_to_zero=TRUE, sig_thresh=0.01,
                                           display_genes=FALSE,
                                           gene_callouts=TRUE)
render_all_lds_plots(pbmc_container, n_rows=2)



# plot single loadings plot
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=5, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                    pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.001, display_genes=FALSE, 
                    gene_callouts=TRUE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, show_xlab=TRUE,
                    show_var_eplained=TRUE, show_all_legends=TRUE)
pbmc_container$plots$single_lds_plot



# get donor sig genes plots
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=1, 
                                       top_n_per_ctype=15)
pbmc_container$plots$donor_sig_genes$Factor1

# get donor sig genes plots
container <- plot_donor_sig_genes(container, factor_select=1, 
                                       top_n_per_ctype=15)
container$plots$donor_sig_genes$Factor1

