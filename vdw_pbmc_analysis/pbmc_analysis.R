
library(scITD)

# counts matrix
pbmc_counts <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_counts.rds')

# meta data matrix
pbmc_meta <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_meta.rds')

# ensembl to gene name conversions
feature.names <- readRDS('/home/jmitchel/data/van_der_wijst/genes.rds')

# put data in project container
pbmc_scMinimal <- instantiate_scMinimal(count_data=pbmc_counts, meta_data=pbmc_meta)
pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = c("CD4+ T", "CD8+ T", "cMonocyte", "CD56(dim) NK", "B"),
                                     gn_convert = feature.names, scale_var = TRUE,
                                     var_scale_power = 1.25,
                                     tucker_type = 'sparse', rotation_type = 'ica',
                                     ncores = 30, rand_seed = 10)

# get sex meta data
pbmc_container <- identify_sex_metadata(pbmc_container)

# finish prepping the data
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="empir", thresh=0.01)

# determine appropriate variance scaling parameter
pbmc_container <- optimize_var_scale_power(pbmc_container, min_ranks_test=c(5,8,5),
                                           max_ranks_test=c(10,15,5),
                                           min_power_test=1.5,
                                           max_power_test=2.25)
pbmc_container$plots$var_scale_plot

# determine appropriate ranks to use for decomposition
pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(6,10,5),
                                         method='svd', num_iter=5, shuffle_level='cells')
pbmc_container$plots$rank_determination_plot

# run tucker
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5), shuffle=FALSE)

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes'),
                                    cluster_by_meta='sex', show_donor_ids = FALSE)
pbmc_container$plots$donor_matrix

# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, n_fibers=100, n_iter=500)

# show that significant genes correspond to large magnitude loadings
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=FALSE, 
                                           nonsig_to_zero=FALSE, 
                                           annot='sig_genes',
                                           sig_thresh=0.05, 
                                           display_genes=FALSE,
                                           gene_callouts=FALSE)
render_all_lds_plots(pbmc_container, n_rows=2)

# plot loadings with significant genes only and gene callouts
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE, 
                                           nonsig_to_zero=TRUE, annot='none',
                                           sig_thresh=0.05, display_genes=FALSE,
                                           gene_callouts=TRUE,
                                           callout_n_gene_per_ctype=5, 
                                           callout_ctypes=list(c(NULL),
                                                               c('CD4+ T',
                                                                 'cMonocyte'),
                                                               c(NULL),
                                                               c('CD8+ T'),
                                                               c(NULL)))
render_all_lds_plots(pbmc_container, n_rows=2)

### plot_donor_sig_genes and run_gsea_one_factor should be made automatic for all
### all factors and resulting plots should be auto formatted...
# generate some plots of scaled expression for top loadings genes of a factor
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=2, 
                                       top_n_per_ctype=c(30,30), 
                                       ctypes_use=c('CD4+ T','cMonocyte'))
pbmc_container$plots$donor_sig_genes$Factor2

pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=4, 
                                       top_n_per_ctype=c(5,20), 
                                       ctypes_use=c('CD4+ T','CD8+ T'))
pbmc_container$plots$donor_sig_genes$Factor4


# run gsea for a few factors
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.01,
                    db_use="GO", collapse_paths=TRUE)
pbmc_container$plots$gsea$Factor2

pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.01,
                                      db_use="GO", collapse_paths=TRUE)
pbmc_container$plots$gsea$Factor4


# run cell subtype proportion analysis
pbmc_container <- get_subtype_prop_associations(pbmc_container,max_res=1.5,
                                                stat_type='adj_pval',
                                                integration_var='lanes')
container$plots$subtype_prop_factor_associations

# generate subcluster plots and association significance for individual subtypes
pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('CD4+ T','CD56(dim) NK',
                                                  'cMonocyte','CD8+ T'),
                                         res=c(.7,.6,.5,.5),
                                         factors=c(2,2,2,4))

render_subtype_plots(pbmc_container)


# get associations between proportions of major cell types and factors
pbmc_container <- get_ctype_prop_associations(pbmc_container,
                                              stat_type='adj_pval')
pbmc_container$plots$ctype_prop_factor_associations


# run stability analysis
pbmc_container <- run_stability_analysis(pbmc_container, downsample_ratio=0.9,
                                         n_iter=500)

# look at number of significant genes for additional factors
pbmc_container <- get_min_sig_genes(pbmc_container,donor_rank_range=c(4:9),thresh=0.05)
pbmc_container$plots$min_sig_genes


