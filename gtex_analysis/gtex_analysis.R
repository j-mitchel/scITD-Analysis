
library(scITD)

# read in the processed data
gtex_tpm_sub_transform <- readRDS(file='/home/jmitchel/data/gtex/7_tissues_counts.rds')
gtex_meta <- readRDS(file='/home/jmitchel/data/gtex/7_tissues_meta.rds')
feature.names.final <- readRDS(file='/home/jmitchel/data/gtex/genes.rds')


## basic scITD pipeline with some slight adjustments for the bulk RNA data
ttypes_use <- unique(gtex_meta$ctypes)
gtex_scMinimal <- instantiate_scMinimal(data_sparse=gtex_tpm_sub_transform,
                                        meta_data=gtex_meta)
gtex_container <- make_new_container(gtex_scMinimal,
                                     ctypes_use = ttypes_use,
                                     gn_convert = feature.names.final,
                                     scale_var = TRUE,
                                     var_scale_power = 1,
                                     tucker_type = 'sparse', 
                                     rotation_type = 'ica',
                                     ncores = 30, rand_seed = 10)

# setting donor_min_cells to 0 because in this case there is only one sample per
# donor per tissue type, whereas usually, there are several cells per donor per cell type
gtex_container <- get_ctype_data(gtex_container,donor_min_cells=0)

# using norm_var method because the anova and empir methods only work with single cell data
gtex_container <- get_ctype_vargenes(gtex_container, method="norm_var", thresh=500)

# get variance scaling power parameter
gtex_container <- optimize_var_scale_power(gtex_container,min_ranks_test=c(2,10,7),
                                           max_ranks_test=c(7,15,7),
                                           min_power_test=0.25,max_power_test=1)
gtex_container$plots$var_scale_plot

# adjust var scale param as needed
gtex_container <- set_experiment_params(gtex_container, var_scale_power = 1)

# determine appropriate ranks
gtex_container <- determine_ranks_tucker(gtex_container, max_ranks_test=c(12,20,7),
                                         method='svd', shuffle_level='tensor', num_iter=5)
gtex_container$plots$rank_determination_plot

# run tucker
gtex_container <- run_tucker_ica(gtex_container, ranks=c(10,20,7), shuffle=FALSE)

# plot donor scores
gtex_container <- plot_donor_matrix(gtex_container, meta_vars=c('sex'),
                                    cluster_by_meta='sex')
gtex_container$plots$donor_matrix

# get significant genes
gtex_container <- run_jackstraw(gtex_container, n_fibers=100, n_iter=500)


# check that significant genes correspond with large magnitude loadings
gtex_container <- get_all_lds_factor_plots(gtex_container, use_sig_only=FALSE, 
                                           nonsig_to_zero=FALSE, 
                                           annot='sig_genes',
                                           sig_thresh=0.05, 
                                           display_genes=FALSE,
                                           gene_callouts=FALSE)

render_all_lds_plots(gtex_container, n_rows=3)


# plot loadings of significant genes only and add gene callouts
gtex_container <- get_all_lds_factor_plots(gtex_container, use_sig_only=TRUE, 
                                           nonsig_to_zero=TRUE, 
                                           annot='sig_genes',
                                           sig_thresh=0.05, 
                                           display_genes=FALSE,
                                           gene_callouts=TRUE)

render_all_lds_plots(gtex_container, n_rows=3)

# run gsea for an interesting factor
gtex_container <- run_gsea_one_factor(gtex_container, factor_select=9, method="fgsea",
                                      thresh=0.005, db_use="GO", num_iter=10000)
gtex_container$plots$gsea[['Factor2']][['up']]

# plot associations with meta data
gtex_container <- plot_meta_associations(gtex_container)
gtex_container$plots$meta_associations



