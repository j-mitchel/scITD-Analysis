library(Seurat)

# load up the subsetted dataset
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')

# set up project parameters
param_list <- initialize_params(ctypes_use = c("B","NK","T4","T8","cDC",
                                               "cM","ncM"),
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
                                                     "Processing_Cohort"),
                                     metadata_col_nm=c('donors',
                                                       'SLE_status',
                                                       'Status',
                                                       'ctypes',
                                                       'sex',
                                                       'Age',
                                                       'pool',
                                                       'processing'))


pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                              scale_var = TRUE, var_scale_power = 1.5)

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(20,30,7),
                                 tucker_type = 'regular', rotation_type = 'ica')

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('pool','processing'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('pool','processing'),
                                    cluster_by_meta='pool',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='rsq')

pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_dscores.pdf", useDingbats = FALSE,
    width = 7, height = 6.5)
pbmc_container$plots$donor_matrix
dev.off()

pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=FALSE,
                                           nonsig_to_zero=FALSE,
                                           display_genes=FALSE,
                                           gene_callouts=FALSE,
                                           show_var_explained = FALSE)

# 2,4,6,15
pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_f2_lds.pdf", useDingbats = FALSE,
    width = 4, height = 5)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["2"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["2"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_f4_lds.pdf", useDingbats = FALSE,
    width = 4, height = 5)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["4"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["4"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_f6_lds.pdf", useDingbats = FALSE,
    width = 4, height = 5)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["6"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["6"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_f15_lds.pdf", useDingbats = FALSE,
    width = 4, height = 5)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["15"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["15"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()



# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(20,30,7), n_fibers=100, n_iter=500,
                                tucker_type='regular', rotation_type='ica')

# saveRDS(pbmc_container[["gene_score_associations"]],file='/home/jmitchel/data/lupus_data/lupus_no_combat_jackstraw.rds')
pbmc_container[["gene_score_associations"]] <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_no_combat_jackstraw.rds')

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=15, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.15, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.15, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)

# run gsea
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=15, method="hypergeometric", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=15,direc='up',thresh=.05)


pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="hypergeometric", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=15,direc='up',thresh=.05)


