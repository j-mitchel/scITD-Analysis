
library(scITD)
library(Seurat)

# load up the subsetted dataset
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_flare_v3.rds')

# set up project parameters
param_list <- initialize_params(ctypes_use = c("B","NK","T4","T8",
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

# get median number cells per donor
ctypes <- levels(pbmc_container$scMinimal_full$metadata$ctypes)
for (ctype in ctypes) {
  tmp <- pbmc_container$scMinimal_full$metadata[pbmc_container$scMinimal_full$metadata$ctypes==ctype,]
  print(ctype)
  print(table(tmp$donors))
}
# we would lose about 9 donors if continue to use cDC cells because there are ~9 with < 10 cDC cells

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=10, gene_min_cells=10,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                              scale_var = TRUE, var_scale_power = 1.5,
                              batch_var='pool')

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,10,6),
                                 tucker_type = 'regular', rotation_type = 'ica')

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Status'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status'),
                                    cluster_by_meta='Status',
                                    show_donor_ids = FALSE,
                                    add_meta_associations=TRUE)

pdf(file = "/home/jmitchel/figures/for_paper/flare/flare_dscores.pdf", useDingbats = FALSE,
    width = 6, height = 6)
pbmc_container$plots$donor_matrix
dev.off()



pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(12,20,6),
                                         shuffle_level='tensor', shuffle_within=NULL,
                                         num_iter=10,
                                         norm_method='trim',
                                         scale_factor=10000,
                                         scale_var=TRUE,
                                         var_scale_power=1.5,
                                         batch_var='pool')
pdf(file = "/home/jmitchel/figures/for_paper/flare/flare_rank_determination.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$rank_determination_plot
dev.off()

# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(5,10,6), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')

# get loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=0.02,
                                           display_genes=FALSE,
                                           gene_callouts=TRUE,
                                           callout_n_gene_per_ctype=8)

# render loadings figure
myfig <- render_multi_plots(pbmc_container,data_type='loadings',max_cols=3)
myfig

pdf(file = "/home/jmitchel/figures/for_paper/flare/flare_loadings.pdf", useDingbats = FALSE,
    width = 17, height = 14)
myfig
dev.off()

# saveRDS(pbmc_container[["gene_score_associations"]],file='/home/jmitchel/data/lupus_data/flare_jackstraw.rds')



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

gsea_fig <- render_multi_plots(pbmc_container,data_type='gsea',max_cols=3)

pdf(file = "/home/jmitchel/figures/for_paper/flare/flare_gsea.pdf", useDingbats = FALSE,
    width = 16, height = 10)
gsea_fig
dev.off()


# cell proportion analysis
pbmc_container <- get_subtype_prop_associations(pbmc_container, max_res=.9, stat_type='adj_pval',
                                                min_cells_group=200, integration_var='pool')

pdf(file = "/home/jmitchel/figures/for_paper/flare/flare_subtype_prop_associations.pdf", useDingbats = FALSE,
    width = 8, height = 9)
pbmc_container$plots$subtype_prop_factor_associations
dev.off()


# saveRDS(pbmc_container$embedding,file='/home/jmitchel/data/lupus_data/flare_conos.rds')
# saveRDS(pbmc_container$subclusters,file='/home/jmitchel/data/lupus_data/flare_subcluster_data.rds')
# saveRDS(pbmc_container$plots$subtype_prop_factor_associations,file='/home/jmitchel/data/lupus_data/flare_sub_associations_plot.rds')


# in case below fn errors in the middle, need to save these objects
orig_embed <- pbmc_container$embedding[["embedding"]]
orig_clusts <- pbmc_container$embedding$clusters$leiden$groups

# to recover original embedding/cell assignments
pbmc_container$embedding[["embedding"]] <- orig_embed
pbmc_container$embedding$clusters$leiden$groups <- orig_clusts



### generate figure with all ctype information for all ctypes/factors
# first determine what resolution of each ctype to choose
all_ctypes=c('cM','cM','cM',
             'ncM','ncM',
             'B','B',
             'T4','T4',
             'T8','T8',
             'NK','NK','NK')
all_res=c(.5,.7,.9,
          .5,.6,
          .7,.8,
          .6,.9,
          .5,.8,
          .6,.7,.8)

pbmc_container <- get_subclust_umap(pbmc_container,all_ctypes=all_ctypes,
                                    all_res=all_res,n_col=4)

pdf(file = "/home/jmitchel/figures/for_paper/flare/flare_all_subc_umaps.pdf", useDingbats = FALSE,
    width = 20, height = 20)
pbmc_container$subc_umap_fig
dev.off()


all_ctypes=c('cM',
             'ncM',
             'B',
             'T4',
             'T8',
             'NK')
all_res=c(.5,
          .5,
          .7,
          .6,
          .5,
          .6)

pbmc_container[["embedding"]][["n.cores"]] <- 5

pbmc_container <- get_subclust_enr_fig(pbmc_container,all_ctypes,all_res)


pdf(file = "/home/jmitchel/figures/for_paper/flare/flare_all_subc_fig5.pdf", useDingbats = FALSE,
    width = 36, height = 16)
pbmc_container$plots$subc_fig
dev.off()










