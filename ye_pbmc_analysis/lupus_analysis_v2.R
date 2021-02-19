
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
                                       top_n_per_ctype=c(5,10,5,5,5,5,5),
                                       show_donor_labels=FALSE,
                                       additional_meta='Status',
                                       add_genes=c('IFI6','ISG15','MX1','XAF1'))
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



# to create GO similarity map for a cell type
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=5, method="hypergeometric", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=5, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=7, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=4, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
pbmc_container[["plots"]][["gsea"]][["7"]]

tmp <- pbmc_container[["gsea_results"]][["5"]][["up"]][["cM"]]
tmp <- pbmc_container[["gsea_results"]][["5"]][["down"]][["cDC"]]
tmp <- tmp[tmp<.05]
tmp <- names(tmp)
mat <- simplifyEnrichment::GO_similarity(tmp,ont='BP')
simplifyEnrichment::simplifyGO(mat)

# or get all sets in up group
all_go <- c()
for (ct in pbmc_container$experiment_params$ctypes_use) {
    tmp <- pbmc_container[["gsea_results"]][["5"]][["down"]][[ct]]
    tmp <- tmp[tmp<.05]
    tmp <- names(tmp)
    all_go <- c(all_go,tmp)
}
all_go <- unique(all_go)
mat <- GO_similarity(all_go,ont='BP')
simplifyGO(mat)
simplifyGO(mat,method="dynamicTreeCut")
simplifyGO(mat,method="apcluster")
simplifyGO(mat,method="hdbscan")
simplifyGO(mat,method="kmeans")

cool <- grid.grabExpr(simplifyGO(mat,method="kmeans"))

cl = simplifyEnrichment::binary_cut(mat)
ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 80))

# testing ht_clusters fn
go_id = simplifyEnrichment::random_GO(500)
mat <- simplifyEnrichment::GO_similarity(go_id,ont='BP')
cl = simplifyEnrichment::binary_cut(mat)
cool <- ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 80))

# general test
library(simplifyEnrichment)
set.seed(888)
go_id = random_GO(500)
mat = GO_similarity(go_id,ont='BP')
df = simplifyGO(mat)

# testing new function to create side by side heatmaps

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_f5_gsea_double2.pdf", useDingbats = FALSE,
    width = 20, height = 15)
cool <- plot_gsea_hmap_w_similarity(pbmc_container,factor_select=2,direc='up',thresh=.05)
dev.off()

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


pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,5), direction=c('down','up'),
                                  compare_type='different', sig_thresh=0.05)
pbmc_container$plots$comparisons[['2_5']]
diff_enr <- get_compare_go_enrich(pbmc_container,'cM',-1)
print(diff_enr[order(diff_enr,decreasing=F)][1:10])


pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,4), direction=c('up','down'),
                                  compare_type='different', sig_thresh=0.05)
pbmc_container$plots$comparisons[['2_4']]
diff_enr <- get_compare_go_enrich(pbmc_container,'T8',-1)
print(diff_enr[order(diff_enr,decreasing=F)][1:15])


# create "vignettes" for individual findings with the subtype proportion analysis

# CD4 naive cell depletion
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





# CD8 granzyme cell expansion
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



# NK granzyme cell expansion
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



# B cell expansion/depletion
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


# cDC cell expansion/depletion
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


# cM cell expansion/depletion 
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




# ncM cell expansion/depletion 
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

### for testing purposes only
# need to look at subtype proportions vs dscores
# Want to:
# double check associations
# see what the deal is with an overlapping top donors for contradictory factors
container <- pbmc_container
ctype <- 'cM'
res <- 0.5
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


subtype_associations <- get_indv_subtype_associations(pbmc_container,donor_props,5)




subtype <- donor_props[,5,drop=FALSE]

# add a second column with value of 1 - first column
subtype <- cbind(subtype,1-subtype)

# get balances
donor_balances <- coda.base::coordinates(subtype)



# for trying coordinates on counts
clusts <- subclusts_num
metadata <- sub_meta_tmp

# got donor_props count version from compute_donor_props
donor_balances <- coda.base::coordinates(donor_props)




# append dscores for factor 4
donor_props2 <- cbind(donor_props,donor_scores[rownames(donor_props),4])
colnames(donor_props2)[ncol(donor_props2)] <- 'dsc'

ggplot(as.data.frame(donor_props2),aes(x=as.numeric(dsc),y=as.numeric(K2))) +
    geom_point()

# append disease status
meta <- unique(container$scMinimal_full$metadata[,c('donors','Status')])
rownames(meta) <- meta$donors
donor_props2 <- cbind(donor_props2,as.character(meta[rownames(donor_props2),'Status']))
colnames(donor_props2)[ncol(donor_props2)] <- 'Status'

ggplot(as.data.frame(donor_props2),aes(x=as.numeric(dsc),y=as.numeric(K3),color=as.factor(Status))) +
    geom_point()



# see if any high loading donors from f1 are low loading in f3
top_f1 <- donor_scores[,1]
top_f1 <- top_f1[order(top_f1,decreasing=TRUE)][1:10]

top_f3 <- donor_scores[,3]
top_f3 <- top_f3[order(top_f3,decreasing=FALSE)][1:10]

sum(names(top_f3) %in% names(top_f1))

top_f3[names(top_f3) %in% names(top_f1)]







## trying to compute donor props using total cell numbers
clusts <- subclusts_num
metadata <- sub_meta_tmp

names(clusts) <- metadata[names(clusts),"donors"]
all_donors <- unique(as.character(metadata$donors))

# store results in df
donor_props <- data.frame(matrix(0,ncol=length(unique(clusts)),nrow = length(all_donors)))
colnames(donor_props) <- sapply(1:ncol(donor_props),function(x) {
    paste0('K',as.character(x))
})
rownames(donor_props) <- all_donors
for (d in all_donors) {
    tmp_clusts <- clusts[names(clusts)==d]
    counts <- table(tmp_clusts)
    names(counts) <- sapply(names(counts),function(x) {
        paste0('K',as.character(x))
    })
    for (j in 1:length(counts)) {
        donor_props[d,names(counts)[j]] <- counts[j]
    }
}
donor_props <- donor_props + 1 #adding pseudocount to avoid infinities when make balances
# donor_props <- t(apply(donor_props, 1, function(i) i/sum(i))) # counts -> props

new_totals <- table(container$scMinimal_full$metadata$donors)
donor_props <- t(sweep(t(donor_props),MARGIN=2,new_totals[rownames(donor_props)],FUN='/'))
subtype_associations <- get_indv_subtype_associations(container,donor_props,1)
subtype_associations





# I want to explore some of the differences a bit to see why some new things are significant and old are not...
container <- pbmc_container
ctype <- 'T8'
res <- 0.6
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

# append dscores for factor 4
donor_props2 <- cbind(donor_props,donor_scores[rownames(donor_props),1])
colnames(donor_props2)[ncol(donor_props2)] <- 'dsc'

ggplot(as.data.frame(donor_props2),aes(x=as.numeric(dsc),y=as.numeric(K3))) +
    geom_point()

# append disease status
meta <- unique(container$scMinimal_full$metadata[,c('donors','Status')])
rownames(meta) <- meta$donors
donor_props2 <- cbind(donor_props2,as.character(meta[rownames(donor_props2),'Status']))
colnames(donor_props2)[ncol(donor_props2)] <- 'Status'

ggplot(as.data.frame(donor_props2),aes(x=as.numeric(dsc),y=as.numeric(K1),color=as.factor(Status))) +
    geom_point()


# trying different way to get balances
subtype_associations <- get_indv_subtype_associations(container,donor_props,5)

lmres <- lm(as.numeric(dsc)~as.numeric(K3),data=as.data.frame(donor_props2))
summary(lmres)

j <- 3
tmp <- donor_props[,j,drop=FALSE]
donor_props <- donor_props[,-j]
donor_props <- cbind(donor_props,tmp)
# donor_balances <- coda.base::coordinates(donor_props)
donor_balances <- compositions::ilr(donor_props)
rownames(donor_balances) <- rownames(donor_props)
donor_balances <- donor_balances[,ncol(donor_balances),drop=FALSE]
donor_props2 <- cbind(donor_balances,donor_scores[rownames(donor_balances),4])
colnames(donor_props2)[ncol(donor_props2)] <- 'dsc'
meta <- unique(container$scMinimal_full$metadata[,c('donors','Status')])
rownames(meta) <- meta$donors
donor_props2 <- cbind(donor_props2,as.character(meta[rownames(donor_props2),'Status']))
colnames(donor_props2)[ncol(donor_props2)] <- 'Status'
colnames(donor_props2)[1] <- 'ilr4'

ggplot(as.data.frame(donor_props2),aes(x=as.numeric(dsc),y=as.numeric(ilr4),color=as.factor(Status))) +
    geom_point()

head(donor_props2)

donor_props2 <- as.data.frame(donor_props2)
donor_props2$dsc <- as.numeric(donor_props2$dsc)
donor_props2$ilr4 <- as.numeric(donor_props2$ilr4)
rownames(donor_props2)[order(donor_props2[,'ilr4'],decreasing=F)][1:10]


container <- get_subclust_enr_hmap(container,all_ctypes,all_res,1:10)
container$plots$subc_enr_hmap


# getting donors with top scores for factor 7
tmp <- container$tucker_results[[1]][,7]
cool <- names(tmp)[tmp>.1]

test <- plot_dscore_enr(pbmc_container,factor_use=6,meta_var='Status')
test




# for testing enrichment of anergic pathway among factorss associated with naive cd4 depletion
my_pathways2 <- my_pathways[c('GSE46242_TH1_VS_ANERGIC_TH1_CD4_TCELL_DN','GSE46242_TH1_VS_ANERGIC_TH1_CD4_TCELL_UP','GSE5960_TH1_VS_ANERGIC_TH1_DN','GSE5960_TH1_VS_ANERGIC_TH1_UP')]

plt <- plotEnrichment(mypaths[[meta_vals[i]]],
                      myranks) + labs(title=paste0('',' - Factor ',as.character(6)))
plt <- plt +
    annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
             label=paste0('adj pval: ',
                          round(fgseaRes[fgseaRes$pathway==meta_vals[i],'padj'],digits=4)))


tmp <- fgsea_res[3,8][[1]][[1]]


pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=6,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=c(5,5,5,5,5,5,5),
                                       show_donor_labels=FALSE,
                                       additional_meta='Status',
                                       add_genes=tmp)

pbmc_container$plots$donor_sig_genes[['6']]

# a few caveats:
# -the relatively borderline pvalue
# -the gene set is from mouse
# -




# LR interaction analysis
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/NicheNet-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('from','to')]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_f1.pdf", useDingbats = FALSE,
    width = 9, height = 9)
tmp <- get_LR_interact(container,lr_pairs,1)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_f2.pdf", useDingbats = FALSE,
    width = 9, height = 9)
tmp <- get_LR_interact(container,lr_pairs,2)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_f4.pdf", useDingbats = FALSE,
    width = 9, height = 9)
tmp <- get_LR_interact(container,lr_pairs,4)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_f5.pdf", useDingbats = FALSE,
    width = 9, height = 9)
tmp <- get_LR_interact(container,lr_pairs,5)
dev.off()











#### exploring new idea for LR interactions
container <- pbmc_container

factor_select <- 5

### prep stuff from lr fn
ctypes_use <- container$experiment_params$ctypes_use

# extract significance of all genes in all ctypes and put in list
sig_vectors <- get_significance_vectors(container,
                                        factor_select, ctypes_use)
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

# set 0 pvals to the min nonzero pval and take -log10
min_nonz <- min(sig_df[sig_df!=0])
sig_df[sig_df==0] <- min_nonz
sig_df <- -log10(sig_df)

# sign sig_df by loading
ldngs <- container$tucker_results[[2]]
genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})
sr_col <- ldngs[factor_select,]
tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)
tmp_casted_num <- tmp_casted_num[rownames(sig_df),colnames(sig_df)]
neg_mask <- tmp_casted_num < 0
sig_df[neg_mask] <- sig_df[neg_mask] * -1
###

lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')

# first identify ligands significantly associated with dscore for a factor
ligs <- lr_pairs[,'ligand']
ligs <- unique(ligs)

sig_ligs <- list()
for (lig in ligs) {
    if (lig %in% rownames(sig_df)) {
        for (ct in ctypes_use) {
            sig_val <- sig_df[lig,ct]
            if (abs(sig_val) > -log10(.01)) {
                sig_ligs[[paste0(lig,'_',ct)]] <- sig_val
            }
        } 
    }
}

print(sig_ligs)

# pick a ligand, see if receptor(s) present in any cell type
mylig <- 'RETN'
mylig <- 'LAIR1'
mylig <- 'THBS1'
mylig <- 'ADM'

myrec <- lr_pairs[lr_pairs$ligand==mylig,]
recs <- lapply(myrec$interaction_name,function(x) {
    tmp <- strsplit(x,split='_')[[1]]
    myrecs <- tmp[2:length(tmp)]
    return(myrecs)
})

print(recs)

# container <- get_pseudobulk(container)
# container <- normalize_pseudobulk(container, method='trim', scale_factor=10000)

lig_pos <- FALSE
r_thresh <- .01
for (j in 1:length(recs)) {
    rs <- recs[[j]]
    num_in_df <- sum(rs %in% rownames(sig_df))
    if (num_in_df == length(rs)) {
        for (ct in ctypes_use) {
            # need to use pseudobulked/normalized expression (but not scaled!)
            pb <- container$scMinimal_ctype[[ct]]$pseudobulk
            checks <- list()
            for (r in rs) {
                # determine direction of ligand expressing donors
                if (lig_pos) {
                    dsc <- container$tucker_results[[1]]
                    tmp <- dsc[,factor_select]
                    tmp <- tmp[order(tmp,decreasing=TRUE)]
                    top_n <- names(tmp)[1:10]
                    d_exp <- sum(pb[top_n,r] > r_thresh)
                    if (d_exp == 10) {
                        checks[[r]] <- TRUE
                    } else {
                        checks[[r]] <- FALSE
                    }
                } else {
                    dsc <- container$tucker_results[[1]]
                    tmp <- dsc[,factor_select]
                    tmp <- tmp[order(tmp,decreasing=FALSE)]
                    top_n <- names(tmp)[1:10]
                    d_exp <- sum(pb[top_n,r] > r_thresh)
                    if (d_exp == 10) {
                        checks[[r]] <- TRUE
                    } else {
                        checks[[r]] <- FALSE
                    }
                }
            }
            if (sum(unlist(checks))==length(checks)) {
                print(rs)
                print(ct)
            }
        }
    }
}


genes_test <- rownames(sig_df)[sig_df[,'ncM']<log10(.01)]
print(genes_test[1:30])

genes_test <- rownames(sig_df)[sig_df[,'cM']<log10(.01)]
print(genes_test[1:30])

genes_test <- rownames(sig_df)[sig_df[,'T4']<log10(.01)]
print(genes_test[1:30])

sig_df_tmp <- sig_df[,colnames(sig_df)!='cM']
sig_genes_other <- rownames(sig_df_tmp)[rowSums(sig_df_tmp<log10(.01)) > 0]
genes_test <- genes_test[!(genes_test %in% sig_genes_other)]

sig_df_tmp <- sig_df[,colnames(sig_df)!='T4']
sig_genes_other <- rownames(sig_df_tmp)[rowSums(sig_df_tmp<log10(.01)) > 0]
genes_test <- genes_test[!(genes_test %in% sig_genes_other)]

dsc <- container$tucker_results[[1]]
tmp <- dsc[,factor_select]

pb <- container$scMinimal_ctype[['ncM']]$pseudobulk
pb <- container$scMinimal_ctype[['cM']]$pseudobulk
pb <- container$scMinimal_ctype[['T4']]$pseudobulk

tmp2 <- cbind(tmp,pb[names(tmp),'HMGB2'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'IL1R2'])
plot(tmp2[,1],tmp2[,2])

# can look for genes with high spearman correlation to HMGB2
res <- list()
for (g in genes_test) {
    tmp2 <- cbind(pb[names(tmp),'UGCG'],pb[names(tmp),g])
    mycor <- cor(tmp2,method='spearman')[1,2]
    res[[g]] <- mycor
}
res <- unlist(res)
res <- res[order(res,decreasing=TRUE)]
print(res[1:10])

tmp2 <- cbind(tmp,pb[names(tmp),'FKBP5'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'ZFAND5'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'ETS2'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'CXCR4'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'ISG15'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'IFI6'])
plot(tmp2[,1],tmp2[,2])

# trying to normalize by receptor levels
rlevs <- pb[names(tmp),'CD36']
sum(rlevs==0)
tmp2 <- cbind(tmp,pb[names(tmp),'SESN1']/rlevs)
plot(tmp2[,1],tmp2[,2])



# negative controls
tmp2 <- cbind(tmp,pb[names(tmp),'ISG15'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'IFI6'])
plot(tmp2[,1],tmp2[,2])



# looking at ligand expression
pb <- container$scMinimal_ctype[['cDC']]$pseudobulk

tmp2 <- cbind(tmp,pb[names(tmp),'RETN'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'THBS1'])
plot(tmp2[,1],tmp2[,2])

pb <- container$scMinimal_ctype[['cM']]$pseudobulk
tmp2 <- cbind(tmp,pb[names(tmp),'ADM'])
plot(tmp2[,1],tmp2[,2])

pb <- container$scMinimal_ctype[['cM']]$pseudobulk
tmp2 <- cbind(tmp,pb[names(tmp),'ZBTB16'])
plot(tmp2[,1],tmp2[,2])


# look at combined ligand levels for all expressing ctypes
pb <- container$scMinimal_ctype[['cM']]$pseudobulk

pb1 <- container$scMinimal_ctype[['cDC']]$pseudobulk
pb2 <- container$scMinimal_ctype[['cM']]$pseudobulk
pb3 <- container$scMinimal_ctype[['ncM']]$pseudobulk

tmp2 <- cbind(tmp,pb1[names(tmp),'ADM']+pb2[names(tmp),'ADM']+pb3[names(tmp),'ADM'])
plot(tmp2[,1],tmp2[,2])


res <- list()
for (i in 1:ncol(pb)) {
    tmp3 <- cbind(pb[names(tmp),i],tmp2[,2])
    mycor <- cor(tmp3,method='pearson')[1,2]
    res[[colnames(pb)[i]]] <- mycor
}
res <- unlist(res)
res <- res[order(res,decreasing=TRUE)]
print(res[1:10])




pb1 <- container$scMinimal_ctype[['T8']]$pseudobulk
pb2 <- container$scMinimal_ctype[['T4']]$pseudobulk

tmp2 <- cbind(tmp,pb1[names(tmp),'CD70']+pb2[names(tmp),'CD70'])
plot(tmp2[,1],tmp2[,2])


# look at receptor levels
pb <- container$scMinimal_ctype[['T4']]$pseudobulk

tmp2 <- cbind(tmp,pb[names(tmp),'CD36'])
plot(tmp2[,1],tmp2[,2])



# incorporate ratio to receptor levels
# use total ligand expression levels from the different cell types (maybe)
# show plots for "negative controls" things like ISG15, which are associated with 
# score but shouldn't show the expected trend


# looking to see if CALCRL explains ZBTB16 levels

pb <- container$scMinimal_ctype[['cM']]$pseudobulk
tmp2 <- cbind(pb[names(tmp)[1:12],'CALCRL'],pb[names(tmp)[1:12],'ZBTB16'])
plot(tmp2[,1],tmp2[,2])

# doesnt really appear to explain it... 
# now I'll do unbiased search of receptors that might explain it
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/NicheNet-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('from','to')]
myres <- list()
rs_test <- unique(lr_pairs$to)
for (r in rs_test) {
    if (r %in% colnames(pb)) {
        tmp2 <- as.data.frame(cbind(pb[names(tmp)[1:12],r],pb[names(tmp)[1:12],'ZBTB16']))
        colnames(tmp2) <- c('rec','sig')
        lmres <- lm(sig~rec,data=tmp2)
        lmres <- summary(lmres)
        # myres[[r]] <-  lmres$fstatistic[[1]]
        myres[[r]] <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
    }
}

myres <- unlist(myres)
myres <- p.adjust(myres,method='fdr')
myres[order(myres,decreasing=F)][1:10]
# myres[order(myres,decreasing=TRUE)][1:10]

tmp2 <- cbind(pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'ZBTB16'])
plot(tmp2[,1],tmp2[,2])


# I also showed increased levels of THBS1 in these patients and this happens to be the
# ligand for the CD47 gene!
# I wonder if I can use this to estimate Kd values for signal transduction since I have
# relative ligand and receptor levels and saturation


# to gain further evidence, I should see if CD47 correlates with residuals of other genes with the same patterns
# will try FKBP5

tmp2 <- cbind(pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'FKBP5'])
plot(tmp2[,1],tmp2[,2])

# it looks pretty okay actually
# whats the pvalue though?

tmp2 <- as.data.frame(cbind(pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'FKBP5']))
colnames(tmp2) <- c('rec','sig')
lmres <- lm(sig~rec,data=tmp2)
lmres <- summary(lmres)


# trying also for CCND3
tmp2 <- cbind(pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'CCND3'])
plot(tmp2[,1],tmp2[,2])






