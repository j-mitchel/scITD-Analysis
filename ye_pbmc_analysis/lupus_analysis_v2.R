
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

# # get assistance with rank determination
# pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(20,40,7),
#                                          shuffle_level='tensor', shuffle_within=NULL,
#                                          num_iter=10,
#                                          norm_method='trim',
#                                          scale_factor=10000,
#                                          scale_var=TRUE,
#                                          var_scale_power=1.5,
#                                          batch_var='pool')
# pdf(file = "/home/jmitchel/figures/for_paper/lupus_combat_rank_det_tensor.pdf", useDingbats = FALSE,
#     width = 9, height = 9)
# pbmc_container$plots$rank_determination_plot
# dev.off()


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(10,20,7),
                                 tucker_type = 'regular', rotation_type = 'ica')


# get factor-meta data associations
# pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Status'),
#                                         stat_use='rsq')
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Status'),
                                        stat_use='pval')

# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status'),
                                    cluster_by_meta = 'Status',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# pdf(file = "/home/jmitchel/figures/for_paper/lupus_combat_dscores_by_status.pdf", useDingbats = FALSE,
#     width = 7, height = 8)
pbmc_container$plots$donor_matrix
dev.off()

# plot donor scores as naturally clustered
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status'),
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

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



# to create GO similarity map for a cell type
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=5, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_f5_gsea_double2.pdf", useDingbats = FALSE,
    width = 20, height = 15)
cool <- plot_gsea_hmap_w_similarity(pbmc_container,factor_select=5,direc='down',thresh=.05)
dev.off()

# make plot of just one cluster
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=1)



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

pbmc_container <- get_subtype_prop_associations(pbmc_container, max_res=.5, stat_type='adj_pval',
                                                min_cells_group=200)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_subtype_prop_associations.pdf", useDingbats = FALSE,
    width = 9, height = 9)
pbmc_container$plots$subtype_prop_factor_associations
dev.off()

# saveRDS(pbmc_container$plots$subtype_prop_factor_associations,file='/home/jmitchel/data/lupus_data/lupus_sub_associations_plot.rds')
# saveRDS(pbmc_container$subclusters,file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')
pbmc_container$subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')


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
          .5,.6,.8,.9,
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


pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_all_subc_total_cells.pdf", useDingbats = FALSE,
    width = 30, height = 16)
pbmc_container$plots$subc_fig
dev.off()


# saveRDS(pbmc_container[["plots"]][["subtype_de"]],file='/home/jmitchel/data/lupus_data/lupus_subc_de_plots.rds')
pbmc_container[["plots"]][["subtype_de"]] <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subc_de_plots.rds')



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



# create "vignettes" for individual findings with the subtype proportion analysis

# CD4 
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_de.pdf", useDingbats = FALSE,
    width = 5, height = 7)
grid::grid.draw(pbmc_container[["plots"]][["subtype_de"]][['T4:0.6']])
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_umap.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container$plots$subc_umaps[['T4:0.6']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='T4',factor_use=4)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_f4_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['T4']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='T4',factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_f5_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['T4']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='T4',factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_f2_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['T4']]
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'T4',0.6,subtype=1,factor_use=4)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_f4_sub1_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'T4',0.6,subtype=1,factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_f2_sub1_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='T4',factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_f5_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['T4']]
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'T4',0.6,subtype=3,factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_f5_sub3_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'T4',0.6,subtype=2,factor_use=4)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_f4_sub2_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'T4',0.6,subtype=1,factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_f5_sub1_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()


# CD8 
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD8_de.pdf", useDingbats = FALSE,
    width = 5, height = 7)
grid::grid.draw(pbmc_container[["plots"]][["subtype_de"]][['T8:0.6']])
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD8_umap.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container$plots$subc_umaps[['T8:0.6']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='T8',factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD8_f2_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['T8']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='T8',factor_use=1)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD8_f1_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['T8']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='T8',factor_use=3)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD8_f3_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['T8']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='T8',factor_use=7)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD8_f7_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['T8']]
dev.off()


myplot <- get_subclust_enr_dotplot(pbmc_container,'T8',0.6,subtype=2,factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD8_f2_sub2_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'T8',0.6,subtype=1,factor_use=1)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD8_f1_sub1_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'T8',0.6,subtype=1,factor_use=7)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD8_f7_sub1_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'T8',0.6,subtype=3,factor_use=3)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD8_f3_sub3_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()



# NK 
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/NK_de.pdf", useDingbats = FALSE,
    width = 5, height = 7)
grid::grid.draw(pbmc_container[["plots"]][["subtype_de"]][['NK:0.6']])
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/NK_umap.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container$plots$subc_umaps[['NK:0.6']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='NK',factor_use=4)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/NK_f4_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['NK']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='NK',factor_use=9)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/NK_f9_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['NK']]
dev.off()


myplot <- get_subclust_enr_dotplot(pbmc_container,'NK',0.6,subtype=1,factor_use=4)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/NK_f4_sub1_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()



# B 
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/B_de.pdf", useDingbats = FALSE,
    width = 5, height = 7)
grid::grid.draw(pbmc_container[["plots"]][["subtype_de"]][['B:0.8']])
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/B_umap.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container$plots$subc_umaps[['B:0.8']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='B',factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/B_f2_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['B']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='B',factor_use=4)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/B_f4_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['B']]
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'B',0.8,subtype=4,factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/B_f2_sub4_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'B',0.8,subtype=1,factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/B_f2_sub1_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'B',0.8,subtype=3,factor_use=4)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/B_f4_sub3_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()


# cDC
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cDC_de.pdf", useDingbats = FALSE,
    width = 5, height = 7)
grid::grid.draw(pbmc_container[["plots"]][["subtype_de"]][['cDC:0.5']])
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cDC_umap.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container$plots$subc_umaps[['cDC:0.5']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='cDC',factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cDC_f5_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['cDC']]
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'cDC',0.5,subtype=2,factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cDC_f5_sub2_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'cDC',0.5,subtype=1,factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cDC_f5_sub1_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()


# cM  
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cM_de.pdf", useDingbats = FALSE,
    width = 5, height = 7)
grid::grid.draw(pbmc_container[["plots"]][["subtype_de"]][['cM:0.5']])
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cM_umap.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container$plots$subc_umaps[['cM:0.5']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='cM',factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cM_f2_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['cM']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='cM',factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cM_f5_collapsed_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['cM']] + scale_x_discrete(labels= c('cM_1','cM_2, cM_4','cM_3','cM_5'))
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'cM',0.5,subtype=4,factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cM_f2_sub4_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'cM',0.5,subtype=2,factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cM_f2_sub2_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'cM',0.5,subtype=4,factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cM_f2_sub4_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'cM',0.5,subtype=2,factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cM_f5_sub2_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'cM',0.5,subtype=4,factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cM_f5_sub4_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()




# ncM 
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/ncM_de.pdf", useDingbats = FALSE,
    width = 5, height = 7)
grid::grid.draw(pbmc_container[["plots"]][["subtype_de"]][['ncM:0.6']])
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/ncM_umap.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container$plots$subc_umaps[['ncM:0.6']]
dev.off()

pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='ncM',factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/ncM_f5_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['ncM']]
dev.off()


myplot <- get_subclust_enr_dotplot(pbmc_container,'ncM',0.6,subtype=3,factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/ncM_f5_sub3_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()




# combine cM subclusters 2 and 4
old <- pbmc_container[["subclusters"]][["cM"]][["res:0.5"]]
tmp <- pbmc_container[["subclusters"]][["cM"]][["res:0.5"]]
tmp <- sapply(tmp,function(x) {
    if (x==4) {
        return(2)
    } else {
        return(x)
    }
})

tmp <- sapply(tmp,function(x) {
    if (x==5) {
        return(4)
    } else {
        return(x)
    }
})

pbmc_container[["subclusters"]][["cM"]][["res:0.5"]] <- tmp
# recreated bplot to show F5 has an overall decrease in HLA expressing cM

# also recreated dotplot see below
myplot <- get_subclust_enr_dotplot(pbmc_container,'cM',0.5,subtype=2,factor_use=5)
myplot <- myplot + ggtitle('cM_2 + cM_4 Proportions')
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cM_f5_sub24_dot.pdf", useDingbats = FALSE,
    width = 5, height = 4.5)
myplot
dev.off()


## running LR analysis with new functions
# prep for new LR analysis
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
    rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
    return(rname)
})
lr_pairs$interaction_name <- NULL

pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                              var_scale_power=1.5, batch_var='pool')
sft_thresh <- c(3,3,2,2,2,2,2)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=5, sig_thresh=0.05, percentile_exp_rec=.9)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_f5.pdf", useDingbats = FALSE,
    width = 15, height = 13)
pbmc_container$plots$lr_analysis[['Factor5']]
dev.off()


# getting GO enrichment HMAPs for modules
ctypes <- c('NK','cM','T4','cDC','ncM')
modules <- c(8,4,7,3,2)

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use='TF')
mod_enr

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('GO'))
mod_enr

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('Hallmark'))
mod_enr


ctypes <- c('NK','cDC','ncM','cM','T4','T8')
modules <- c(2,1,1,1,2,3)

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.1, db_use=c('TF'))
mod_enr

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('Hallmark'))
mod_enr



# running gsea to see what functionally distinguishes these cd4 subtype 2 genes
mod_genes <- c('ITGB1', 'S100A4', 'S100A10', 'ANXA1', 'LGALS1', 'KLRB1', 'CRIP1', 'SH3BGRL3')
db_use='immuno'
pvals[order(pvals,decreasing=FALSE)][1:10]
## top 3 hits are all naive vs memory cd4 cell down, so these are likely markers of memory cd4 cells

# now will try for the T4_3 groups of genes
mod_genes <- c('JUN', 'JUNB', 'CD69', 'FOS', 'DUSP1', 'CXCR4', 'TSC22D3', 'LEPROTL1', 'ZFP36L2', 'PNRC1')
db_use='immuno'
db_use='ctype'
pvals[order(pvals,decreasing=FALSE)][1:30]
## top hits are memory cd4 sets as well as rosiglitazone (PPARg) set and others






## creating loadings heatmaps with gene callouts for select gene sets.
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=5, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_EXOCYTOSIS")
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_INTERLEUKIN_1_MEDIATED_SIGNALING_PATHWAY") #2nd in cM, ncM, cDC
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY") #2nd in cM, ncM, cDC
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_POSITIVE_REGULATION_OF_HYDROLASE_ACTIVITY") #2nd in t4, t8
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_EXOCYTOSIS") #2nd in cM
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_POSITIVE_REGULATION_OF_MEMBRANE_PERMEABILITY") #2nd in t4
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_CYTOSKELETON_ORGANIZATION") # down in cM
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II") # down in cM (SAME THING WE SEE WITH SUBTYPE DEPLETION)
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_RESPONSE_TO_HORMONE")
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_RESPONSE_TO_CORTICOSTEROID")
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_RESPONSE_TO_LIPID")
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",'GO_OSSIFICATION')
gset_cmap <- c('blue',
               'orange')

gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_INTERLEUKIN_1_MEDIATED_SIGNALING_PATHWAY",
           "GO_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
           "GO_POSITIVE_REGULATION_OF_MEMBRANE_PERMEABILITY",
           "GO_CYTOSKELETON_ORGANIZATION",
           "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II")

gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
           "GO_POSITIVE_REGULATION_OF_MEMBRANE_PERMEABILITY",
           "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II",
           "GO_RESPONSE_TO_HORMONE")

gset_cmap <- c('blue',
               'orange',
               'red',
               'dark green',
               'black')
names(gset_cmap) <- gsets

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/f5_lds2.pdf", useDingbats = FALSE,
    width = 9, height = 12)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=5, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gsets, le_set_colormap=gset_cmap, le_set_num_per=5, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)
dev.off()



# will try to identify good summary gene sets to plot on the loadings heatmap
cool <- plot_gsea_hmap_w_similarity(pbmc_container,factor_select=5,direc='up',thresh=.05)
cool <- plot_gsea_hmap_w_similarity(pbmc_container,factor_select=5,direc='down',thresh=.05)
dev.off()
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=1)
dev.off()


gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",'GO_INTERLEUKIN_2_PRODUCTION')
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",'GO_EXOCYTOSIS')
gsets <- c('GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I')
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",'GO_RESPONSE_TO_INTERFERON_GAMMA')
gset_cmap <- c('blue')
gset_cmap <- c('blue','orange')
names(gset_cmap) <- gsets

gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           'GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I',
           'GO_MYELOID_CELL_DIFFERENTIATION')
gset_cmap <- c('blue',
               'orange',
               'dark green')
names(gset_cmap) <- gsets

# pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/f2_lds.pdf", useDingbats = FALSE,
#     width = 9, height = 10)
dev.off()
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gsets, le_set_colormap=gset_cmap, le_set_num_per=6, show_le_legend=TRUE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)
dev.off()


# compare similar factors to find out what's similar/different
pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,5), direction=c('up','down'),
                                  compare_type='different', sig_thresh=0.05)
pbmc_container$plots$comparisons[['2_5']]
diff_enr <- get_compare_go_enrich(pbmc_container,'B',-1)
print(diff_enr[order(diff_enr,decreasing=F)][1:10])


# now picking some gene sets for f4
cool <- plot_gsea_hmap_w_similarity(pbmc_container,factor_select=4,direc='up',thresh=.05)
plot_gsea_sub(pbmc_container,factor_select=4,direc='up',thresh=.05,clust_select=4)

gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           'GO_CELLULAR_EXTRAVASATION',
           'GO_REGULATION_OF_WNT_SIGNALING_PATHWAY')
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           'GO_CELLULAR_EXTRAVASATION',
           'GO_IMPORT_INTO_CELL')
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           'GO_CELLULAR_EXTRAVASATION',
           'GO_OSSIFICATION')

gset_cmap <- c('blue',
               'orange',
               'dark green')
names(gset_cmap) <- gsets

dev.off()
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=4, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gsets, le_set_colormap=gset_cmap, le_set_num_per=6, show_le_legend=TRUE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)
dev.off()

color_sets <- c('blue',
                'dark orange',
                'red',
                'dark gray',
                'black',
                'forest green',
                'magenta',
                'brown',
                'purple')
sets_plot <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
               'GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I',
               'GO_MYELOID_CELL_DIFFERENTIATION',
               'GO_CELLULAR_EXTRAVASATION',
               "GO_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
               "GO_POSITIVE_REGULATION_OF_MEMBRANE_PERMEABILITY",
               "GO_RESPONSE_TO_HORMONE",
               "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II",
               'GO_OSSIFICATION')
color_map <- color_sets
names(color_map) <- sets_plot


##### creating loadings + gene sets for JUST FACTOR 5
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=5, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
           "GO_POSITIVE_REGULATION_OF_MEMBRANE_PERMEABILITY",
           "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II",
           "GO_RESPONSE_TO_HORMONE")

gset_cmap <- color_map[match(gsets,names(color_map))]

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=5, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gsets, le_set_colormap=gset_cmap, le_set_num_per=4, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)


hm_list <- plot_select_sets(pbmc_container, factors_all=c(5), 
                            sets_plot=gsets,
                            thresh=.05,
                            color_sets=gset_cmap)

p1 <- pbmc_container$plots$all_lds_plots[['5']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["5"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/f5_lds_gsets.pdf", useDingbats = FALSE,
    width = 12, height = 12)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()



##### creating loadings + gene sets for factors 2,4,5. Still sort of need to do it one at a time, but can
# generate the select gene sets all together for each...

hm_list <- plot_select_sets(pbmc_container, factors_all=c(5), 
                            sets_plot=sets_plot,
                            thresh=.05,
                            color_sets=color_sets)

p1 <- pbmc_container$plots$all_lds_plots[['5']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["5"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/f5_lds_gsets_combo.pdf", useDingbats = FALSE,
    width = 12, height = 12)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()

# now doing factor 4
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=4, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           'GO_CELLULAR_EXTRAVASATION',
           'GO_OSSIFICATION')

gset_cmap <- color_map[match(gsets,names(color_map))]

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=4, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gsets, le_set_colormap=gset_cmap, le_set_num_per=6, show_le_legend=TRUE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)


hm_list <- plot_select_sets(pbmc_container, factors_all=c(4), 
                            sets_plot=sets_plot,
                            thresh=.05,
                            color_sets=color_sets)

p1 <- pbmc_container$plots$all_lds_plots[['4']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["4"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/f4_lds_gsets_combo.pdf", useDingbats = FALSE,
    width = 12, height = 12)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()




# now doing factor 2
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           'GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I',
           'GO_MYELOID_CELL_DIFFERENTIATION')

gset_cmap <- color_map[match(gsets,names(color_map))]

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gsets, le_set_colormap=gset_cmap, le_set_num_per=6, show_le_legend=TRUE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)


hm_list <- plot_select_sets(pbmc_container, factors_all=c(2), 
                            sets_plot=sets_plot,
                            thresh=.05,
                            color_sets=color_sets)

p1 <- pbmc_container$plots$all_lds_plots[['2']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["2"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/f2_lds_gsets_combo.pdf", useDingbats = FALSE,
    width = 12, height = 12)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()




# seeing what gene sets look good for factor 1 if I decide to include it in presentation
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=1, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
           "GO_CELL_KILLING",
           "GO_INTERLEUKIN_8_PRODUCTION")

gset_cmap <- c('blue',
                'dark orange',
                'red',
                'dark gray')

gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_CELL_KILLING",
           "GO_INTERLEUKIN_8_PRODUCTION")

gset_cmap <- c('blue',
               'dark orange',
               'red')

names(gset_cmap) <- gsets

dev.off()
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gsets, le_set_colormap=gset_cmap, le_set_num_per=7, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)

# will try to identify good summary gene sets to plot on the loadings heatmap
cool <- plot_gsea_hmap_w_similarity(pbmc_container,factor_select=1,direc='up',thresh=.05)
cool <- plot_gsea_hmap_w_similarity(pbmc_container,factor_select=1,direc='down',thresh=.05)
dev.off()
plot_gsea_sub(pbmc_container,factor_select=1,direc='down',thresh=.05,clust_select=1)
dev.off()





####### regenerating above stuff with new all gene sets to include
color_sets <- c('blue',
                'salmon',
                'turquoise',
                'dark orange',
                'red',
                'dark gray',
                'black',
                'forest green',
                'magenta',
                'brown',
                'purple')
sets_plot <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
               "GO_CELL_KILLING",
               "GO_INTERLEUKIN_8_PRODUCTION",
               'GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I',
               'GO_MYELOID_CELL_DIFFERENTIATION',
               'GO_CELLULAR_EXTRAVASATION',
               "GO_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
               "GO_POSITIVE_REGULATION_OF_MEMBRANE_PERMEABILITY",
               "GO_RESPONSE_TO_HORMONE",
               "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II",
               'GO_OSSIFICATION')
color_map <- color_sets
names(color_map) <- sets_plot

## F5
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
           "GO_POSITIVE_REGULATION_OF_MEMBRANE_PERMEABILITY",
           "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II",
           "GO_RESPONSE_TO_HORMONE")

gset_cmap <- color_map[match(gsets,names(color_map))]

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=5, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gsets, le_set_colormap=gset_cmap, le_set_num_per=4, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)


hm_list <- plot_select_sets(pbmc_container, factors_all=c(5), 
                            sets_plot=sets_plot,
                            thresh=.05,
                            color_sets=color_sets)

p1 <- pbmc_container$plots$all_lds_plots[['5']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["5"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/f5_lds_gsets_combo2.pdf", useDingbats = FALSE,
    width = 12, height = 12)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()


## F1
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_CELL_KILLING",
           "GO_INTERLEUKIN_8_PRODUCTION")
gsets <- c("GO_CELL_KILLING")

gset_cmap <- color_map[match(gsets,names(color_map))]

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gsets, le_set_colormap=gset_cmap, le_set_num_per=10, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)


hm_list <- plot_select_sets(pbmc_container, factors_all=c(1), 
                            sets_plot=sets_plot,
                            thresh=.05,
                            color_sets=color_sets)

p1 <- pbmc_container$plots$all_lds_plots[['1']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["1"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/f1_lds_gsets_combo2.pdf", useDingbats = FALSE,
    width = 12, height = 12)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()

## F2
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           'GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I',
           'GO_MYELOID_CELL_DIFFERENTIATION')

gset_cmap <- color_map[match(gsets,names(color_map))]

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gsets, le_set_colormap=gset_cmap, le_set_num_per=6, show_le_legend=TRUE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)


hm_list <- plot_select_sets(pbmc_container, factors_all=c(2), 
                            sets_plot=sets_plot,
                            thresh=.05,
                            color_sets=color_sets)

p1 <- pbmc_container$plots$all_lds_plots[['2']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["2"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/f2_lds_gsets_combo2.pdf", useDingbats = FALSE,
    width = 12, height = 12)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()


## F4
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           'GO_CELLULAR_EXTRAVASATION',
           'GO_OSSIFICATION')

gset_cmap <- color_map[match(gsets,names(color_map))]

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=4, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gsets, le_set_colormap=gset_cmap, le_set_num_per=6, show_le_legend=TRUE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)


hm_list <- plot_select_sets(pbmc_container, factors_all=c(4), 
                            sets_plot=sets_plot,
                            thresh=.05,
                            color_sets=color_sets)

p1 <- pbmc_container$plots$all_lds_plots[['4']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["4"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/f4_lds_gsets_combo2.pdf", useDingbats = FALSE,
    width = 12, height = 12)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()







# create a zoom in plot for my slides
cool <- plot_gsea_hmap_w_similarity(pbmc_container,factor_select=1,direc='up',thresh=.05)
cool <- plot_gsea_hmap_w_similarity(pbmc_container,factor_select=5,direc='down',thresh=.05)
dev.off()
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=1)
dev.off()


# trying LR analysis for factor 4 to see if anything interesting pops out that could explain the T4 depletion
# prep for new LR analysis
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
    rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
    return(rname)
})
lr_pairs$interaction_name <- NULL

pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=1.5, batch_var='pool')
sft_thresh <- c(3,3,2,2,2,2,2)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=4, sig_thresh=0.05, percentile_exp_rec=.9)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_f4_2.pdf", useDingbats = FALSE,
    width = 15, height = 17)
pbmc_container$plots$lr_analysis[['Factor4']]
dev.off()


ctypes <- c('T4','T4')
modules <- c(3,4)

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('BioCarta'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('GO'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('TF')) 
# KLF7 regulates cell proliferation, differentiation, survival, HES2 can be activated by PI3K pathway which is activated by ICOS
# STAT5A signaling is linked to ICOS, SAFB2 involved in cell cycle regulation, TEAD2 involved in proliferation, DIDO1 does apoptosis

mod_enr

mod_genes <- pbmc_container[["module_genes"]][["T4"]]
mygenes <- names(mod_genes)[mod_genes==3]
write.csv(mygenes,file='/home/jmitchel/test_genes.csv',row.names=F,quote=F)




## now trying with the signed network
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=1.5, batch_var='pool')
sft_thresh <- c(9,9,7,4,7,4,6)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=4, sig_thresh=0.05, percentile_exp_rec=.9)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_test.pdf", useDingbats = FALSE,
    width = 15, height = 17)
pbmc_container$plots$lr_analysis[['Factor4']]
dev.off()

ctypes <- rep('T4',2)
modules <- c(5,6)

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.1, db_use=c('KEGG'))
mod_enr
# DYRK1A could be interesting as it can potentially regulate cell proliferation. Specifically, it inhibits apoptosis
# and promotes survival. The fact that we have this upregulated may suggest a reaction to the low T4 numbers...
# WNT could be a potential upstream signal says this article:
# https://jcs.biologists.org/content/125/3/561
# We see enrichment of WNT signaling which could be an interesting connection. We also see enrichment
# of apoptosis gene sets in the downregulated modules for the lupus donors

ctypes <- rep('T4',5)
modules <- c(1,3,4,5,6)
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.1, db_use=c('Hallmark'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.1, db_use=c('GO'))
mod_enr

# another interesting one is the ICOSL interaction with T4_3
# we have enrichemnt of targets for KAT2A and CREB3L4 TFs
# KAT2A is interesting because it is linked to NOTCH signagling and we also have
# a ligand correlated that binds to NOTCH receptors. see here
# https://www.genecards.org/cgi-bin/carddisp.pl?gene=KAT2A



#### looking at directionality of signed and unsigned sets
## now trying with the signed network
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=1.5, batch_var='pool')
sft_thresh <- c(9,9,6,4,7,4,6)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

# get a module from cM
mod_genes <- pbmc_container[["module_genes"]][["T4"]]
mygenes <- names(mod_genes)[mod_genes==3]

# now compute correlation hmap of expression of these genes
pb <- pbmc_container$scMinimal_ctype[['T4']]$pseudobulk
# check if they are all (or mostly) positively correlated with one another

mygenes <- mygenes[mygenes %in% colnames(pb)]

Heatmap(cor(pb[,mygenes]))

# get the module eigengene
MEs <- pbmc_container[["module_eigengenes"]][['T4']]
ME <- MEs[,3]
names(ME) <- rownames(MEs)

# make heatmap of genes by donors with annotation for donor eigengene expression
# myannot <- ComplexHeatmap::rowAnnotation(me_val=anno_simple(ME),
#                                                      show_annotation_name=TRUE,
#                                                      col = list(cell_types = mycol))
col_fun = colorRamp2(c(min(ME), 0, max(ME)), c("blue", "white", "red"))
ra <- rowAnnotation(me_val=ME,col = list(me_val=col_fun),
                        border=TRUE,show_legend=TRUE)
# myannot <- ComplexHeatmap::rowAnnotation(me_val=anno_simple(ME),
#                                          show_annotation_name=TRUE, col=list(me_val=col_fun))
Heatmap(pb[,mygenes],right_annotation=ra)

## now add a row annotation indicating the genes in certain gene sets
m_df <- data.frame()
m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                    category = "C5", subcategory = "BP"))
my_pathways = split(m_df$gene_symbol, f = m_df$gs_name)
gset <- my_pathways[['GO_CELL_CYCLE']]

m_df <- data.frame()
m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                    category = "C2", subcategory = "CP:BIOCARTA"))
my_pathways = split(m_df$gene_symbol, f = m_df$gs_name)
gset <- my_pathways[['BIOCARTA_TCR_PATHWAY']]

in_set <- sapply(mygenes,function(x) {
    if (x %in% gset) {
        return(1)
    } else {
        return(0)
    }
})
names(in_set) <- mygenes
ba <- HeatmapAnnotation(gs=in_set,
                    border=TRUE,show_legend=TRUE)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/f4_T4_mod3_tcr.pdf", useDingbats = FALSE,
    width = 40, height = 17)
Heatmap(pb[,mygenes],
        right_annotation=ra,
        bottom_annotation=ba)
dev.off()

# plotting levels of 1 CDK cell cycle inhibitory gene to check its direction
# CDKN1A
dsc <- pbmc_container$tucker_results[[1]][,4]
pb <- pbmc_container$scMinimal_ctype[['T4']]$pseudobulk[,'CD3D']
pb <- pbmc_container$scMinimal_ctype[['T4']]$pseudobulk[,'KLF3']
pb <- pbmc_container$scMinimal_ctype[['T4']]$pseudobulk[,'CNOT1']
pb <- pbmc_container$scMinimal_ctype[['T4']]$pseudobulk[,'WHAMM']
pb <- pbmc_container$scMinimal_ctype[['T4']]$pseudobulk[,'SDCBP']
pb <- pbmc_container$scMinimal_ctype[['T4']]$pseudobulk[,'DCTN6']

cool <- as.data.frame(cbind(dsc,pb[names(dsc)]))
colnames(cool) <- c('dsc','pb')
class(cool$dsc)
ggplot(cool,aes(x=dsc,y=pb)) +
    geom_point()

# # get donors with largest positive scores
# top_d <- ME[order(ME,decreasing=TRUE)][1:10]
# 
# # get donors with lowest scores
# bot_d <- ME[order(ME,decreasing=FALSE)][1:10]
# 
# # do these donors have higher expression of the module genes consistently?
# pb[names(top_d),mygenes[1]]
# 
# pb[names(bot_d),mygenes[1]]

# loop through the module genes to see which ones have highest cor with the factor
all_cors <- c()
for (g in mygenes) {
    mycor <- cor(pbmc_container$scMinimal_ctype[['T4']]$pseudobulk[names(dsc),g],dsc)
    all_cors <- c(all_cors,mycor)
}
names(all_cors) <- mygenes
all_cors <- all_cors[order(all_cors,decreasing=T)]
all_cors <- all_cors[order(all_cors,decreasing=F)]
f_cors <- all_cors[1:25]

## actually I think it makes more sense to look for gnes most correlated with the ligand, ICOSLG or THBS1

all_cors <- c()
for (g in mygenes) {
    mycor <- cor(pbmc_container$scMinimal_ctype[['T4']]$pseudobulk[names(dsc),g],pbmc_container$scMinimal_ctype[['cM']]$pseudobulk[names(dsc),'ICOSLG'])
    all_cors <- c(all_cors,mycor)
}
names(all_cors) <- mygenes
all_cors <- all_cors[order(all_cors,decreasing=T)]
all_cors <- all_cors[order(all_cors,decreasing=F)]
cool <- all_cors[1:30]
cool[names(cool)%in%names(f_cors)] #get ones highly associated with the factor

## it definitely looks like more of the genes highly correlated go in the expected direction (pro-proliferative ones are higher when ICOLG is high and
# antiproliferative ones are high when ICOS is low)

# now I'll look for a similar thing with THBS1
# firs to plot overall module correlation though
plot_mod_and_lig(pbmc_container,factor_select=4,mod_ct='T4',mod=3,lig_ct='cM',lig='THBS1')

all_cors <- c()
for (g in mygenes) {
    mycor <- cor(pbmc_container$scMinimal_ctype[['T4']]$pseudobulk[names(dsc),g],pbmc_container$scMinimal_ctype[['cM']]$pseudobulk[names(dsc),'THBS1'])
    all_cors <- c(all_cors,mycor)
}
names(all_cors) <- mygenes
all_cors <- all_cors[order(all_cors,decreasing=T)]
all_cors <- all_cors[order(all_cors,decreasing=F)]
all_cors[1:10]



## looking at whether some of the strongest linked ICOSLG targets are cell subtype dependent
# genes included: CNOT, DYNLT1, SBDS, CETN3
subclusts <- pbmc_container[["subclusters"]][["T4"]][["res:0.6"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
subclusts[1]
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
    return(paste0('T4','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
sum(subclusts==6)

# set up project parameters
param_list <- initialize_params(ctypes_use = c("T4_1","T4_2","T4_3"),
                                ncores = 30, rand_seed = 10)

pbmc_sub_container <- make_new_container(seurat_obj=pbmc_sub,
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

pbmc_sub_container <- parse_data_by_ctypes(pbmc_sub_container)

pbmc_sub_container <- clean_data(pbmc_sub_container, donor_min_cells=3, gene_min_cells=10)

pbmc_sub_container <- get_pseudobulk(pbmc_sub_container)

pbmc_sub_container <- normalize_pseudobulk(pbmc_sub_container, method='trim', scale_factor=10000)

# now look at expresion of these genes in T4 subtypes
# genes included: CNOT1, DYNLT1, SBDS, CETN3
d_included <- rownames(pbmc_sub_container$scMinimal_ctype[[1]]$pseudobulk)
dsc2 <- dsc[d_included]

pb <- pbmc_sub_container$scMinimal_ctype[['T4_1']]$pseudobulk[names(dsc2),'CNOT1']
mycor1 <- cor(dsc2,pb)
tmp <- as.data.frame(cbind(dsc2,pb))
colnames(tmp) <- c('dsc','myexp')
p1 <- ggplot(tmp,aes(x=dsc,y=myexp)) +
    geom_point() +
    xlab(paste0('factor4',' donor score')) +
    ylab('CNOT1 expression in T4_1')

gn_test <- c('CNOT1', 'DYNLT1', 'SBDS', 'CETN3')
sub_test <- c("T4_1","T4_2","T4_3")
myres <- as.data.frame(matrix(ncol=4,nrow=3))
colnames(myres) <- gn_test
rownames(myres) <- sub_test
for (sub in sub_test) {
    for (g in gn_test) {
        pb <- pbmc_sub_container$scMinimal_ctype[[sub]]$pseudobulk[names(dsc2),g]
        mycor1 <- cor(dsc2,pb)
        myres[sub,g] <- mycor1
    }
}

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
Heatmap(as.matrix(myres), name = 'expres-dsc cor',
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col=col_fun,
        row_names_side='left',
        column_names_side='top',
        column_names_rot = 40,
        row_title = 'T4 subtype',
        border=TRUE,
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
            grid.text(round(myres[i,j],2), x, y)
        })

# see if ICOS is expressed on naive T4 cells T4_1
pbmc_sub_container$scMinimal_ctype[['T4_1']]$pseudobulk[names(dsc2),'ICOS']
## yeah it looks like the vast majority of donors express ICOS in naive T4 cells and the others

### testing for interaction between age and Status in determining T4_1 proportions
# using age cutoff of <=48 yrs and >48 years as in that abstract
container <- pbmc_container
ctype <- 'T4'
res <- .6
resolution_name <- paste0('res:',as.character(res))
subclusts <- container$subclusters[[ctype]][[resolution_name]]

# append large cell type name to subclusters
subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})

# limit cells in subclusts to those that we actually have scores for
donor_scores <- container$tucker_results[[1]]
donor_vec <- container$scMinimal_full$metadata[names(subclusts),'donors']
subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]

# make subtype association plot
subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
scMinimal <- container$scMinimal_ctype[[ctype]]
sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

# get donor proportions of subclusters
donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)

head(donor_props)

# get meta data for donor status
meta <- container$scMinimal_full$metadata[,c('donors','Status','Age')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL
head(meta)

dsc <- container$tucker_results[[1]][,4]
tmp <- as.data.frame(cbind(donor_props[names(dsc),1],dsc,as.character(meta[names(dsc),'Status']),meta[names(dsc),'Age']))
colnames(tmp) <- c('prop','dsc','Status','Age')
head(tmp)
tmp$prop <- as.numeric(tmp$prop)
tmp$Age <- as.numeric(tmp$Age)

tmp$age_category <- sapply(tmp$Age,function(x) {
    if (x > 48) {
        return('age > 48')
    } else {
        return('age <= 48')
    }
})

ggplot(tmp,aes(x=as.factor(age_category),y=prop,color=as.factor(Status))) +
    geom_violin() +
    xlab('Age Category') +
    ylab('Naive CD4 Proportion') +
    labs(color = "Status")

tmp$age_category <- as.factor(tmp$age_category)
old <- as.data.frame(tmp[tmp$Age>48,])
young <- as.data.frame(tmp[tmp$Age<=48,])

t.test(prop~Status,data=old)
t.test(prop~Status,data=young)

t.test(Age~Status,data=tmp)


lmr <- lm(tmp$prop~tmp$Age)
summary(lmr)

plot(tmp$prop~tmp$Age)

tmp$dsc <- as.numeric(tmp$dsc)
lm1 <- lm(prop~Age,data=tmp)
lm2 <- lm(prop~Age+dsc,data=tmp)
anova(lm1,lm2)

ggplot(tmp,aes(x=Age,y=prop,color=Status)) +
    geom_point() +
    xlab('Age') +
    ylab('Naive CD4 Proportion') +
    labs(color = 'Status')

lm1 <- lm(prop~dsc,data=tmp)
lm2 <- lm(prop~Age+dsc,data=tmp)
anova(lm1,lm2)

test <- tmp[tmp$Status=='Healthy',]
lm3 <- lm(Age~dsc,data=test)
summary(lm3)


test <- tmp[tmp$Status=='Healthy',]
ggplot(test,aes(x=Age,y=prop)) +
    geom_point(color='#F8766D') +
    xlab('Age') +
    ylab('Naive CD4 Proportion') +
    ggtitle('Healthy Donors') +
    theme(plot.title = element_text(hjust = 0.5))

test <- tmp[tmp$Status=='Managed',]
ggplot(test,aes(x=Age,y=prop)) +
    geom_point(color='#00BFC4') +
    xlab('Age') +
    ylab('Naive CD4 Proportion') +
    ggtitle('Managed Donors') +
    theme(plot.title = element_text(hjust = 0.5))


test <- tmp[tmp$Status=='Healthy',]
ggplot(test,aes(x=dsc,y=Age)) +
    geom_point(color='#F8766D') +
    xlab('Factor 4 Donor Score') +
    ylab('Age') +
    ggtitle('Healthy Donors') +
    theme(plot.title = element_text(hjust = 0.5))

test <- tmp[tmp$Status=='Managed',]
ggplot(test,aes(x=dsc,y=Age)) +
    geom_point(color='#00BFC4') +
    xlab('Factor 4 Donor Score') +
    ylab('Age') +
    ggtitle('Managed Donors') +
    theme(plot.title = element_text(hjust = 0.5))


test <- tmp[tmp$Status=='Managed',]
lm3 <- lm(Age~dsc,data=test)
summary(lm3)


test <- tmp[tmp$Status=='Managed',]
lm3 <- lm(Age~prop,data=test)
summary(lm3)


## plotting IL16 levels vs donor score for factor 5
meta <- container$scMinimal_full$metadata[,c('donors','Status','Age')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL
head(meta)


plot_mod_and_lig(pbmc_container,factor_select=5,mod_ct='ncM',mod=4,lig_ct='ncM',lig='TNF')
factor_select=4
lig_ct='T4'
lig='IL16'
# lig_ct='ncM'
# lig='TNF'
mod=4
mod_ct='cM'
dsc <- container$tucker_results[[1]][,factor_select]
lig_exp <- container$scMinimal_ctype[[lig_ct]]$pseudobulk[,lig]
MEs <- container[["module_eigengenes"]][[mod_ct]]
ME <- MEs[,mod]
names(ME) <- rownames(MEs)

tmp <- as.data.frame(cbind(dsc[names(ME)],ME,lig_exp[names(ME)]))
colnames(tmp) <- c('dsc','ME','lig_exp')
tmp$Status <- as.character(meta[names(ME),'Status'])
head(tmp)
class(tmp$dsc)

mycor1 <- cor(tmp$dsc,tmp$lig_exp)
p1 <- ggplot(tmp,aes(x=dsc,y=lig_exp,color=Status)) +
    geom_point() +
    xlab(paste0('factor',factor_select,' donor score')) +
    ylab(paste0(lig,' expression in ',lig_ct)) +
    annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
             label=paste0('pearson r = ',round(mycor1,digits=3)))



## regenerating the LR analysis vignette trend plots for presentation
plot_mod_and_lig(pbmc_container,factor_select=4,mod_ct='T4',mod=3,lig_ct='cM',lig='ICOSLG')










### testing out plotting all proportions for all donors
container <- pbmc_container
ctype <- 'T4'
res <- .6
resolution_name <- paste0('res:',as.character(res))
subclusts <- container$subclusters[[ctype]][[resolution_name]]

# append large cell type name to subclusters
subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})

# limit cells in subclusts to those that we actually have scores for
donor_scores <- container$tucker_results[[1]]
donor_vec <- container$scMinimal_full$metadata[names(subclusts),'donors']
subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]

# make subtype association plot
subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
scMinimal <- container$scMinimal_ctype[[ctype]]
sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

# get donor proportions of subclusters
donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)

## new stuff

# order donor props by dscore
dsc <- donor_scores[,1]
dsc <- dsc[order(dsc,decreasing=T)]
donor_props <- donor_props[names(dsc),]

col_fun = colorRamp2(c(0, 1), c("white", "red"))
Heatmap(donor_props, name='subtype proportions',
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col = col_fun,
        show_row_names = FALSE,
        border=TRUE)

# now try selecting just the top and bottom donors to show
nrow(donor_props)
top_bot_d_props <- donor_props[c(1:10,162:171),]
# top_bot_d_props <- donor_props[c(1:20,152:171),]

# need to add dsc annotation
dsc_sub <- dsc[c(1:10,162:171)]
# dsc_sub <- dsc[c(1:20,152:171)]
col_fun2 = colorRamp2(c(min(dsc_sub),0, max(dsc_sub)), c("blue","white", "red"))
la = rowAnnotation(F4_dscore=dsc_sub,col=list(F4_dscore=col_fun2))

# make new colnames
new_col_names <- sapply(1:ncol(top_bot_d_props),function(x) {
    paste0(ctype,'_',x)
})

col_fun = colorRamp2(c(0, 1), c("white", "red"))
Heatmap(top_bot_d_props, name='subtype proportions',
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col = col_fun,
        column_labels=new_col_names,
        show_row_names = FALSE,
        row_title = "Donors",
        column_title = "CD4 T Cell Subpopulations",
        column_title_side = "bottom",
        border=TRUE,
        left_annotation=la)


# make dotplot showing T4_4 IFN subtype proportions for factor 1 donors
myplot <- get_subclust_enr_dotplot(pbmc_container,'T4',0.6,subtype=4,factor_use=1)









