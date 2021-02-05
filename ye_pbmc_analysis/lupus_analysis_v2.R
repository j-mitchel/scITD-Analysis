
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


# testing out logit version of plot_scores_by_meta
pbmc_container <- plot_scores_by_meta(pbmc_container,'Status')

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_status_associations.pdf", useDingbats = FALSE,
    width = 11, height = 7)
pbmc_container$plots$indv_meta_scores_associations
dev.off()



# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(10,20,7), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')

# saveRDS(pbmc_container[["gene_score_associations"]],file='/home/jmitchel/data/lupus_data/lupus_jackstraw.rds')
pbmc_container[["gene_score_associations"]] <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_jackstraw.rds')


# get loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=0.02,
                                           display_genes=FALSE,
                                           gene_callouts=TRUE,
                                           callout_n_gene_per_ctype=8)

# render loadings figure
myfig <- render_multi_plots(pbmc_container,data_type='loadings',max_cols=5)
myfig

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_loadings.pdf", useDingbats = FALSE,
    width = 30, height = 14)
myfig
dev.off()




# generating all dsig genes plots
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=1,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['1']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=2,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['2']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=3,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['3']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=4,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['4']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=5,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['5']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=6,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['6']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=7,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=c(10,50,5,10,5,5,5),
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['7']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=8,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['8']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=9,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=7,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['9']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=10,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=c(5,5,5,50,5,5,5),
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['10']]

dsig_fig <- render_multi_plots(pbmc_container,data_type='dgenes',max_cols=5)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_dsig_genes.pdf", useDingbats = FALSE,
    width = 35, height = 35)
dsig_fig
dev.off()


# generate dsig genes plot with donor scores heatmap appended so we can see that
# stripes correspond with donors from the other interferon pathways


pbmc_container <- plot_donor_sig_genes_v2(pbmc_container, factor_select=5,
                                          ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                          top_n_per_ctype=15, show_donor_labels=FALSE,
                                          additional_meta='Status')
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_dsig_f5_w_scores.pdf", useDingbats = FALSE,
    width = 20, height = 17)
draw(pbmc_container$plots$donor_sig_genes[['5']], column_title='Donors', column_title_side='bottom')
dev.off()


# generate gsea plots
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=1, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=TRUE)
pbmc_container[["plots"]][["gsea"]][["1"]]
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=TRUE)
pbmc_container[["plots"]][["gsea"]][["2"]]
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=3, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=TRUE)
pbmc_container[["plots"]][["gsea"]][["3"]]
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=4, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=TRUE)
pbmc_container[["plots"]][["gsea"]][["4"]]
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=5, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=TRUE)
pbmc_container[["plots"]][["gsea"]][["5"]]
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=6, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=TRUE)
pbmc_container[["plots"]][["gsea"]][["6"]]
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=7, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=TRUE)
pbmc_container[["plots"]][["gsea"]][["7"]]
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=8, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=TRUE)
pbmc_container[["plots"]][["gsea"]][["8"]]
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=9, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=TRUE)
pbmc_container[["plots"]][["gsea"]][["9"]]
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=10, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=TRUE)
pbmc_container[["plots"]][["gsea"]][["10"]]


gsea_fig <- render_multi_plots(pbmc_container,data_type='gsea',max_cols=3)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_gsea.pdf", useDingbats = FALSE,
    width = 30, height = 39)
gsea_fig
dev.off()

# saveRDS(pbmc_container[["gsea_results"]],file='/home/jmitchel/data/lupus_data/lupus_gsea_results.rds')



## now for cell subtype proportion analysis
# add conos object for cell proportion analysis
con <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_conos.rds')
pbmc_container$embedding <- con
rm(con)
gc()


# in case below fn errors in the middle, need to save these objects
orig_embed <- pbmc_container$embedding[["embedding"]]
orig_clusts <- pbmc_container$embedding$clusters$leiden$groups

# to recover original embedding/cell assignments
pbmc_container$embedding[["embedding"]] <- orig_embed
pbmc_container$embedding$clusters$leiden$groups <- orig_clusts


# large number of cores seems to hamper some stuff below
pbmc_container$embedding$n.cores <- 5

pbmc_container <- get_subtype_prop_associations(pbmc_container, max_res=.9, stat_type='adj_pval',
                                                min_cells_group=200)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_subtype_prop_associations.pdf", useDingbats = FALSE,
    width = 9, height = 9)
pbmc_container$plots$subtype_prop_factor_associations
dev.off()

# saveRDS(pbmc_container$plots$subtype_prop_factor_associations,file='/home/jmitchel/data/lupus_data/lupus_sub_associations_plot.rds')
# saveRDS(pbmc_container$subclusters,file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')



### generate figure with all ctype information for all ctypes/factors
# first determine what resolution of each ctype to choose
all_ctypes=c('cM','cM','cM','cM',
             'ncM','ncM','ncM','ncM',
             'cDC','cDC','cDC','cDC',
             'B','B','B','B',
             'T4','T4','T4','T4',
             'T8','T8','T8','T8',
             'NK','NK','NK','NK')
all_res=c(.5,.7,.8,.9,
          .5,.6,.7,.8,
          .5,.7,.8,.9,
          .6,.7,.8,.9,
          .6,.7,.8,.9,
          .5,.6,.8,.9,
          .6,.7,.8,.9)


# get all subcluster umaps
pbmc_container <- get_subclust_umap(pbmc_container,all_ctypes=all_ctypes,
                                    all_res=all_res,n_col=4)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/all_subc_umaps.pdf", useDingbats = FALSE,
    width = 20, height = 25)
pbmc_container$subc_umap_fig
dev.off()

all_ctypes=c('cM',
             'ncM',
             'cDC',
             'B',
             'T4',
             'T8',
             'NK')
all_res=c(.5,
          .6,
          .5,
          .8,
          .6,
          .6,
          .6)

pbmc_container <- get_subclust_enr_fig(pbmc_container,all_ctypes,all_res)


pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_all_subc_fig2.pdf", useDingbats = FALSE,
    width = 30, height = 16)
pbmc_container$plots$subc_fig
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
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_dscores_umap.pdf", useDingbats = FALSE,
    width = 16, height = 5)
pbmc_container$plots$dscores_umap
dev.off()

# make UMAP using donor meta data
pbmc_container <- plot_donor_umap(pbmc_container, 
                                  color_by_meta=c('Status','sex','pool'), n_col=3)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_meta_umap.pdf", useDingbats = FALSE,
    width = 16, height = 5)
pbmc_container$plots$dscores_umap
dev.off()




# compare similar factors to find out what's similar/different
pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,5), direction=c('up','down'),
                                  compare_type='same', sig_thresh=0.02)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_f2_f5_same_up.pdf", useDingbats = FALSE,
    width = 10, height = 14)
pbmc_container$plots$comparisons[['2_5']]
dev.off()

pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,5), direction=c('up','down'),
                                  compare_type='different', sig_thresh=0.02)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_f2_f5_different_up.pdf", useDingbats = FALSE,
    width = 10, height = 14)
pbmc_container$plots$comparisons[['2_5']]
dev.off()

pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,5), direction=c('down','up'),
                                  compare_type='same', sig_thresh=0.02)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_f2_f5_same_down.pdf", useDingbats = FALSE,
    width = 10, height = 14)
pbmc_container$plots$comparisons[['2_5']]
dev.off()

pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,5), direction=c('down','up'),
                                  compare_type='different', sig_thresh=0.02)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_f2_f5_different_down.pdf", useDingbats = FALSE,
    width = 10, height = 35)
pbmc_container$plots$comparisons[['2_5']]
dev.off()


