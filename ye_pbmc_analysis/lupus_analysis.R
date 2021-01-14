
library(scITD)
library(Seurat)

# load up the subsetted dataset
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v2.rds')

# prep the container
pbmc_scMinimal <- seurat_to_scMinimal(pbmc,normalize_counts=FALSE,
                                      metadata_cols=c('ind_cov_batch_cov',
                                                      "SLE_status",
                                                      "Status",
                                                      "cg_cov",
                                                      "sex",
                                                      "age",
                                                      "batch_cov",
                                                      "Processing_Cohort"),
                                      metadata_col_nm=c('donors',
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


pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = c("B","NK","T4","T8","cDC",
                                                    "cM","ncM"),
                                     scale_var = TRUE,
                                     var_scale_power = 1.5,
                                     tucker_type = 'regular', rotation_type = 'ica',
                                     ncores = 30, rand_seed = 10)

# pbmc_container <- make_new_container(pbmc_scMinimal,
#                                      ctypes_use = c("NK","T4","T8","cM"),
#                                      scale_var = TRUE,
#                                      var_scale_power = 1.5,
#                                      tucker_type = 'regular', rotation_type = 'ica',
#                                      ncores = 30, rand_seed = 10)

rm(pbmc_scMinimal)
gc()


# finish prepping the data
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=500)


# determine appropriate variance scaling parameter
pbmc_container <- optimize_var_scale_power(pbmc_container, min_ranks_test=c(5,8,5),
                                           max_ranks_test=c(10,15,5),
                                           min_power_test=1.5,
                                           max_power_test=2.25)
pbmc_container$plots$var_scale_plot

# determine appropriate ranks to use for decomposition
pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(30,210,7),
                                         method='svd', num_iter=5, shuffle_level='cells',
                                         batch_var='pool')
pbmc_container$plots$rank_determination_plot

# run tucker
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(13,20,7), shuffle=FALSE, batch_var='pool')

# get metadata associations for donor scores plot
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','pool','Status'))

# plot donor scores first by clustering by sex
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status','pool','processing'),
                                    show_donor_ids = FALSE,add_meta_associations=T,cluster_by_meta='Status')
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status','pool','processing'),
                                    show_donor_ids = FALSE, add_meta_associations=T)
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('Status'),
                                    cluster_by_meta='Status', show_donor_ids = FALSE,
                                    add_meta_associations=T)
pbmc_container <- plot_donor_matrix(pbmc_container,show_donor_ids = FALSE, add_meta_associations=T)
pbmc_container$plots$donor_matrix

# get significant genes via jackstraw
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
                    show_var_explained=TRUE, show_all_legends=TRUE)
pbmc_container$plots$single_lds_plot



# get sig genes plots in donor-centric manner
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=1, 
                                       top_n_per_ctype=15)
pbmc_container$plots$donor_sig_genes$Factor1







