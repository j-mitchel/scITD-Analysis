
library(scITD)
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


# ## initial test with no batch correction
# pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
#                               norm_method='trim', scale_factor=10000,
#                               vargenes_method='norm_var_pvals', vargenes_thresh=.05,
#                               scale_var = TRUE, var_scale_power = 1.5)
# 
# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(20,30,7),
#                                  tucker_type = 'regular', rotation_type = 'ica')
# 
# # get factor-meta data associations
# pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','pool','processing'))
# 
# # plot donor scores
# pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','pool','processing'),
#                                     cluster_by_meta='pool',
#                                     show_donor_ids = FALSE,
#                                     add_meta_associations=TRUE)
# 
# pdf(file = "/home/jmitchel/figures/for_paper/lupus_no_combat_dscores.pdf", useDingbats = FALSE,
#     width = 9, height = 8)
# pbmc_container$plots$donor_matrix
# dev.off()



## now with batch correction applied
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                              scale_var = TRUE, var_scale_power = 1.5,
                              batch_var='pool')

# get assistance with rank determination
pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(20,40,7),
                                         shuffle_level='tensor', shuffle_within=NULL,
                                         num_iter=10,
                                         norm_method='trim',
                                         scale_factor=10000,
                                         scale_var=TRUE,
                                         var_scale_power=1.5,
                                         batch_var='pool')
pdf(file = "/home/jmitchel/figures/for_paper/lupus_combat_rank_det_tensor.pdf", useDingbats = FALSE,
    width = 9, height = 9)
pbmc_container$plots$rank_determination_plot
dev.off()


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(10,20,7),
                                 tucker_type = 'regular', rotation_type = 'ica')


# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','pool','processing','Status'))

# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status'),
                                    cluster_by_meta = 'Status',
                                    show_donor_ids = FALSE,
                                    add_meta_associations=TRUE)

# pdf(file = "/home/jmitchel/figures/for_paper/lupus_combat_dscores_by_status.pdf", useDingbats = FALSE,
#     width = 7, height = 8)
pbmc_container$plots$donor_matrix
dev.off()

# plot donor scores as naturally clustered
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status'),
                                    show_donor_ids = FALSE,
                                    add_meta_associations=TRUE)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_combat_dscores.pdf", useDingbats = FALSE,
    width = 7, height = 8)
pbmc_container$plots$donor_matrix
dev.off()


pbmc_container <- plot_scores_by_meta(pbmc_container,'Status')

pdf(file = "/home/jmitchel/figures/for_paper/lupus_status_associations.pdf", useDingbats = FALSE,
    width = 11, height = 7)
pbmc_container$plots$indv_meta_scores_associations
dev.off()


# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(10,20,7), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')

saveRDS(pbmc_container[["gene_score_associations"]],file='/home/jmitchel/data/lupus_data/lupus_jackstraw.rds')

# get loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=0.02,
                                           display_genes=FALSE,
                                           gene_callouts=TRUE,
                                           callout_n_gene_per_ctype=5)

# makes tall rendering
myfig <- render_multi_plots(pbmc_container,data_type='loadings')
myfig

pdf(file = "/home/jmitchel/figures/for_paper/lupus_combat_loadings.pdf", useDingbats = FALSE,
    width = 18, height = 18)
myfig
dev.off()

# makes wide rendering
myfig <- render_multi_plots(pbmc_container,data_type='loadings',max_cols=5)
myfig

pdf(file = "/home/jmitchel/figures/for_paper/lupus_combat_loadings.pdf", useDingbats = FALSE,
    width = 27, height = 11)
myfig
dev.off()



# Factor 5 deep dive
# run gsea 
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=5, method="hypergeometric", thresh=0.001,
                                      db_use=c("GO"), collapse_paths=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f5_enrich.pdf", useDingbats = FALSE,
    width = 9, height = 16)
pbmc_container[["plots"]][["gsea"]][["5"]]
dev.off()

oldplt <- pbmc_container[["plots"]][["gsea"]][["5"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f5_dsig_genes.pdf", useDingbats = FALSE,
    width = 9, height = 13)
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=5,
                                       top_n_per_ctype=8, show_donor_labels=FALSE,
                                       additional_meta='Status')
dev.off()


# Factor 2 deep dive
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=2,
                                       top_n_per_ctype=8, show_donor_labels=FALSE,
                                       additional_meta='Status')
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f2_dsig_genes.pdf", useDingbats = FALSE,
    width = 9, height = 13)
pbmc_container$plots$donor_sig_genes[['2']]
dev.off()

pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.01,
                                      db_use=c("GO"), collapse_paths=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f2_gsea.pdf", useDingbats = FALSE,
    width = 7, height = 9)
pbmc_container[["plots"]][["gsea"]][["2"]]
dev.off()


# Factor 4 deep dive
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=4,
                                       top_n_per_ctype=8, show_donor_labels=FALSE,
                                       additional_meta='Status')
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f4_dsig_genes.pdf", useDingbats = FALSE,
    width = 9, height = 13)
pbmc_container$plots$donor_sig_genes[['4']]
dev.off()

pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=4, method="fgsea", thresh=0.01,
                                      db_use=c("GO"), collapse_paths=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f4_gsea.pdf", useDingbats = FALSE,
    width = 7, height = 9)
pbmc_container[["plots"]][["gsea"]][["4"]]
dev.off()

# Factor 1 deep dive
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=1,
                                       top_n_per_ctype=8, show_donor_labels=FALSE,
                                       additional_meta='Status')
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f1_dsig_genes.pdf", useDingbats = FALSE,
    width = 9, height = 13)
pbmc_container$plots$donor_sig_genes[['1']]
dev.off()

pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=1, method="fgsea", thresh=0.01,
                                      db_use=c("GO"), collapse_paths=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f1_gsea.pdf", useDingbats = FALSE,
    width = 7, height = 9)
pbmc_container[["plots"]][["gsea"]][["1"]]
dev.off()


pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=3,
                                       top_n_per_ctype=8, show_donor_labels=FALSE,
                                       additional_meta='Status')
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f3_dsig_genes.pdf", useDingbats = FALSE,
    width = 9, height = 13)
pbmc_container$plots$donor_sig_genes[['3']]
dev.off()

pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=3, method="fgsea", thresh=0.01,
                                      db_use=c("GO"), collapse_paths=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f3_gsea.pdf", useDingbats = FALSE,
    width = 7, height = 18)
pbmc_container[["plots"]][["gsea"]][["3"]]
dev.off()

myplot <- plot_gsea_hmap(pbmc_container,factor_select=1,thresh=.05)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f1_gsea.pdf", useDingbats = FALSE,
    width = 9, height = 16)
myplot
dev.off()


pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=6,
                                       top_n_per_ctype=8, show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['6']]




pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=7, ctypes_use=c('cM','T4','T8'),
                                       top_n_per_ctype=c(5,5,30), show_donor_labels=FALSE,
                                       additional_meta='Status')

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f7_dsig_genes.pdf", useDingbats = FALSE,
    width = 9, height = 13)
pbmc_container$plots$donor_sig_genes[['7']]
dev.off()

pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=7, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f7_gsea.pdf", useDingbats = FALSE,
    width = 7, height = 10)
pbmc_container[["plots"]][["gsea"]][["7"]]
dev.off()












# add conos object for cell proportion analysis
con <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_conos.rds')
pbmc_container$embedding <- con
rm(con)
gc()


# saveRDS(pbmc_container[["gsea_results"]],file='/home/jmitchel/data/lupus_data/lupus_gsea_results.rds')
# saveRDS(pbmc_container[["gsea_results"]],file='/home/jmitchel/data/lupus_data/lupus_plots.rds')

# myplots <- readRDS('/home/jmitchel/data/lupus_data/lupus_plots.rds')
# mygsea <- readRDS('/home/jmitchel/data/lupus_data/lupus_gsea_results.rds')



pbmc_container <- get_subtype_prop_associations(pbmc_container, max_res=.8, stat_type='adj_pval')

pdf(file = "/home/jmitchel/figures/for_paper/lupus_subtype_associations.pdf", useDingbats = FALSE,
    width = 9, height = 9)
pbmc_container$plots$subtype_prop_factor_associations
dev.off()


# in case below fn errors in the middle, need to save these objects
# save original embedding
orig_embed <- pbmc_container$embedding[["embedding"]]
pbmc_container$embedding[["embedding"]] <- orig_embed

# save original cluster labels
orig_clusts <- pbmc_container$embedding$clusters$leiden$groups
pbmc_container$embedding$clusters$leiden$groups <- orig_clusts

pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('cM','cDC','T4','B'),
                                         res=c(.7,.5,.8,.6),
                                         factors=c(5,5,5,5))


pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('cM'),
                                         res=c(.7),
                                         factors=c(5))

subc_fig <- render_subtype_plots_v2(pbmc_container)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f5_cM_subtypes.pdf", useDingbats = FALSE,
    width = 4, height = 9)
subc_fig
dev.off()

pbmc_container[["plots"]][["subc_plots"]][["cM_res:0.7"]] <- NULL


pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('T4'),
                                         res=c(.6),
                                         factors=c(4))

subc_fig <- render_subtype_plots_v2(pbmc_container)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f4_T4_subtypes.pdf", useDingbats = FALSE,
    width = 4, height = 9)
subc_fig
dev.off()

pbmc_container[["plots"]][["subc_plots"]][["T4_res:0.6"]] <- NULL


pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('T8'),
                                         res=c(.5),
                                         factors=c(3))

subc_fig <- render_subtype_plots_v2(pbmc_container)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f3_T8_subtypes.pdf", useDingbats = FALSE,
    width = 4, height = 9)
subc_fig
dev.off()

pbmc_container[["plots"]][["subc_plots"]][["T8_res:0.5"]] <- NULL


pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('T8'),
                                         res=c(.5),
                                         factors=c(1))

subc_fig <- render_subtype_plots_v2(pbmc_container)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f1_T8_subtypes.pdf", useDingbats = FALSE,
    width = 4, height = 9)
subc_fig
dev.off()

pbmc_container[["plots"]][["subc_plots"]][["T8_res:0.5"]] <- NULL


pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('T8'),
                                         res=c(.5),
                                         factors=c(2))

subc_fig <- render_subtype_plots_v2(pbmc_container)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f2_T8_subtypes.pdf", useDingbats = FALSE,
    width = 4, height = 9)
subc_fig
dev.off()

pbmc_container[["plots"]][["subc_plots"]][["T8_res:0.5"]] <- NULL


pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('T8'),
                                         res=c(.5),
                                         factors=c(7))

subc_fig <- render_subtype_plots_v2(pbmc_container)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f7_T8_subtypes.pdf", useDingbats = FALSE,
    width = 4, height = 9)
subc_fig
dev.off()


pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('T4'),
                                         res=c(.7),
                                         factors=c(2))

subc_fig <- render_subtype_plots_v2(pbmc_container)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f7_T8_subtypes.pdf", useDingbats = FALSE,
    width = 4, height = 9)
subc_fig
dev.off()


pbmc_container[["plots"]][["subc_plots"]][["T4_res:0.7"]] <- NULL

pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('cDC'),
                                         res=c(.5),
                                         factors=c(5))

subc_fig <- render_subtype_plots_v2(pbmc_container)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f5_cDC_subtypes.pdf", useDingbats = FALSE,
    width = 4, height = 9)
subc_fig
dev.off()


# now do whole ctype proportion association tests
pbmc_container <- get_ctype_prop_associations(pbmc_container,'adj_pval',n_col=3)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_major_ctype_associations.pdf", useDingbats = FALSE,
    width = 7, height = 7)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()


# create main embedding for lupus data
pbmc_container$embedding$plotGraph(alpha=0.1)



# create UMAP using donor scores
pbmc_container <- plot_donor_umap(pbmc_container, color_by_factor=c(1:10), 
                                  color_by_meta=NULL, n_col=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_dscores_umap.pdf", useDingbats = FALSE,
    width = 15, height = 5)
pbmc_container$plots$dscores_umap
dev.off()





### new workflow for subtype analysis

# add conos object for cell proportion analysis
con <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_conos.rds')
pbmc_container$embedding <- con
rm(con)
gc()

# large number of cores seems to hamper some stuff below
pbmc_container$embedding$n.cores <- 5

# save original embedding
orig_embed1 <- pbmc_container$embedding[["embedding"]]
pbmc_container$embedding[["embedding"]] <- orig_embed

# save original cluster labels
orig_clusts1 <- pbmc_container$embedding$clusters$leiden$groups
pbmc_container$embedding$clusters$leiden$groups <- orig_clusts


pbmc_container <- get_subtype_prop_associations(pbmc_container, max_res=.7, stat_type='adj_pval',
                                                min_cells_group=200)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_subtype_prop_associations.pdf", useDingbats = FALSE,
    width = 9, height = 9)
pbmc_container$plots$subtype_prop_factor_associations
dev.off()

saveRDS(pbmc_container$plots$subtype_prop_factor_associations,file='/home/jmitchel/data/lupus_data/lupus_sub_associations_plot.rds')
saveRDS(pbmc_container$subclusters,file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')



### generate figure with all ctype information for all ctypes/factors
# first determine what resolution of each ctype to choose
all_ctypes=c('cM','cM',
             'ncM','ncM',
             'cDC','cDC',
             'B','B',
             'T4','T4',
             'T8','T8',
             'NK','NK')
all_res=c(.5,.7,
          .5,.6,
          .6,.7,
          .6,.7,
          .6,.7,
          .5,.6,
          .6,.7)


# get all subcluster umaps
pbmc_container <- get_subclust_umap(pbmc_container,all_ctypes=all_ctypes,
                                    all_res=all_res,n_col=4)

pdf(file = "/home/jmitchel/test.pdf", useDingbats = FALSE,
    width = 15, height = 15)
pbmc_container$subc_umap_fig
dev.off()

all_ctypes=c('cM',
             'ncM',
             'cDC',
             'B',
             'T4',
             'T8',
             'NK')
all_res=c(.7,
          .5,
          .6,
          .6,
          .6,
          .5,
          .7)

pbmc_container <- get_subclust_enr_fig(pbmc_container,all_ctypes,all_res)


pdf(file = "/home/jmitchel/figures/for_paper/lupus_all_subc_fig3.pdf", useDingbats = FALSE,
    width = 30, height = 18)
pbmc_container$plots$subc_fig
dev.off()

# saving the DE heatmaps specifically
saveRDS(tmp,file='/home/jmitchel/data/lupus_data/lupus_sub_de_plots.rds')
saveRDS(pbmc_container[["plots"]][["subtype_de"]],file='/home/jmitchel/data/lupus_data/lupus_sub_de_plots_no_auc_lim.rds')








# old pipeline

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
                                    show_donor_ids = FALSE,add_meta_associations=TRUE)
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status','pool','processing'),
                                    show_donor_ids = FALSE,add_meta_associations=TRUE,cluster_by_meta='Status')
pbmc_container$plots$donor_matrix

# get significant genes via jackstraw
pbmc_container <- run_jackstraw(pbmc_container, n_fibers=250, n_iter=500)

# show all loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=FALSE,
                                           nonsig_to_zero=FALSE,
                                           display_genes=FALSE,
                                           gene_callouts=FALSE)
render_multi_plots(pbmc_container, n_rows=4, 'loadings')

# show all loadings plots sig genes only
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE, sig_thresh=0.01,
                                           display_genes=FALSE,
                                           gene_callouts=TRUE,
                                           callout_n_gene_per_ctype=3)
render_multi_plots(pbmc_container, n_rows=3, 'loadings')

# plot Factor 2 loadings plot
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                    pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.01, display_genes=FALSE,
                    gene_callouts=TRUE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, show_xlab=TRUE,
                    show_var_explained=TRUE)


# get sig genes plots in donor-centric manner for factor 2
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=2,
                                       top_n_per_ctype=10)

# run gsea for factor 2
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.01,
                                      db_use="GO", collapse_paths=TRUE)
pbmc_container$plots$gsea$`2`

# run cell subtype proportion analysis
pbmc_container <- get_subtype_prop_associations(pbmc_container,max_res=1.5,
                                                stat_type='adj_pval',
                                                integration_var='pool')
container$plots$subtype_prop_factor_associations

# generate subcluster plots and association significance for individual subtypes
pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('CD4+ T','CD56(dim) NK',
                                                  'cMonocyte','CD8+ T'),
                                         res=c(.7,.6,.5,.5),
                                         factors=c(2,2,2,4))

render_subtype_plots(pbmc_container)




