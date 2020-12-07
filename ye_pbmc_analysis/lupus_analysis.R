
library(scITD)
library(edgeR)
library(Seurat)
library(SeuratDisk)

# load the cleaned data
readRDS()

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

pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = c("B","NK","T4","T8","cDC",
                                                    "cM","ncM"),
                                     scale_var = TRUE,
                                     var_scale_power = 1.25,
                                     tucker_type = 'sparse', rotation_type = 'ica',
                                     ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = c("NK","T4","T8","cM"),
                                     scale_var = TRUE,
                                     var_scale_power = 1.5,
                                     tucker_type = 'sparse', rotation_type = 'ica',
                                     ncores = 30, rand_seed = 10)

rm(pbmc_scMinimal)
gc()


# finish prepping the data
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=1000)

# trying mixed model for batch correct
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="anova", thresh=.05)


# testing out batch correction on collapsed data before var gene selection
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- collapse_by_donors(pbmc_container, shuffle=FALSE)
pbmc_container <- apply_pseudobulk_batch_correct(pbmc_container,'pool')
pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=1000) #make necessary changes
## probably need more changes to make it work with run tucker ica because this recalls collapse by donors
## maybe need to have the apply pseudobulk fn put the data in data sparse. I think I do something similar with
## the gtex data so take a look at that though that may have been relating to vargenes fn only

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
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,10,4), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,4), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,10,7), shuffle=FALSE)


# plot donor scores first by clustering by sex
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status','pool'),
                                    cluster_by_meta='sex', show_donor_ids = FALSE)
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status','pool'),
                                    show_donor_ids = FALSE)
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('Status','pool'),
                                    cluster_by_meta='pool', show_donor_ids = FALSE)
pbmc_container$plots$donor_matrix

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

