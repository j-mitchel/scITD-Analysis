library(Seurat)

# load up the subsetted dataset
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$Status=='Managed']
pbmc <- subset(pbmc,cells = cells_keep)

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
# # pdf(file = "/home/jmitchel/figures/for_paper/lupus_combat_rank_det_tensor.pdf", useDingbats = FALSE,
# #     width = 9, height = 9)
# pbmc_container$plots$rank_determination_plot
# # dev.off()


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,20,7),
                                 tucker_type = 'regular', rotation_type = 'ica')

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,25,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,26,7),
                                 tucker_type = 'regular', rotation_type = 'ica')

# get factor-meta data associations
# pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Status'),
#                                         stat_use='rsq')
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing'),
                                        stat_use='pval')

# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# pdf(file = "/home/jmitchel/figures/for_paper/sle_only_dscores.pdf", useDingbats = FALSE,
#     width = 5, height = 6)
pbmc_container$plots$donor_matrix
dev.off()

# plot donor scores as naturally clustered
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status'),
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# pdf(file = "/home/jmitchel/figures/for_paper/lupus_combat_dscores.pdf", useDingbats = FALSE,
#     width = 7, height = 8)
pbmc_container$plots$donor_matrix
dev.off()


# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(5,20,7), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')

# saveRDS(pbmc_container[["gene_score_associations"]],file='/home/jmitchel/data/lupus_data/sle_only_jackstraw.rds')
pbmc_container[["gene_score_associations"]] <- readRDS(file='/home/jmitchel/data/lupus_data/sle_only_jackstraw.rds')


# get loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=0.01,
                                           display_genes=FALSE,
                                           gene_callouts=TRUE,
                                           callout_n_gene_per_ctype=5)

# render loadings figure
myfig <- render_multi_plots(pbmc_container,data_type='loadings',max_cols=3)
myfig

pdf(file = "/home/jmitchel/figures/for_paper/sle_only_loadings.pdf", useDingbats = FALSE,
    width = 13, height = 12)
myfig
dev.off()

# plotting IFN response gene expression against f1
d_exp <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,'ISG20']
dsc <- pbmc_container$tucker_results[[1]][,1]
tmp <- cbind.data.frame(dsc[names(d_exp)],d_exp)
colnames(tmp) <- c('dscore','expres')

pdf(file = "/home/jmitchel/figures/for_paper/sle_only_f1_IFN.pdf", useDingbats = FALSE,
    width = 4, height = 3.25)
ggplot(tmp,aes(x=dscore,y=expres)) +
  geom_point() +
  ylab('ISG20 expression (CD4+ T)') +
  xlab('Factor 1 donor scores') +
  theme_bw()
dev.off()


pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=1, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=1,direc='down',thresh=.05)
plot_gsea_sub(pbmc_container,factor_select=1,direc='down',thresh=.05,clust_select=8)

gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_LYMPHOCYTE_ACTIVATION",
           "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
           "GO_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
           "GO_INTERLEUKIN_1_MEDIATED_SIGNALING_PATHWAY",
           "GO_NIK_NF_KAPPAB_SIGNALING",
           "GO_INTERLEUKIN_8_PRODUCTION",
           "GO_CELL_CYCLE",
           "GO_INTERLEUKIN_1_BETA_PRODUCTION",
           "GO_INTERLEUKIN_10_PRODUCTION",
           "GO_RNA_PHOSPHODIESTER_BOND_HYDROLYSIS",
           "GO_SECRETION",
           "GO_APOPTOTIC_PROCESS",
           "GO_METAL_ION_HOMEOSTASIS")

gset_cmap <- rep('black',length(gsets))

names(gset_cmap) <- gsets

pdf(file = "/home/jmitchel/figures/for_paper/sle_only_f1_go.pdf", useDingbats = FALSE,
    width = 10, height = 8)
hm_list <- plot_select_sets(pbmc_container, 1, gsets, thresh=.05, color_sets=gset_cmap, 
                            cl_rows=TRUE)
dev.off()




pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=4, method="hypergeometric", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=4,direc='up',thresh=.05)




## running LR analysis with new functions
# prep for new LR analysis
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
  rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
  return(rname)
})
lr_pairs$interaction_name <- NULL


# factor 5 LR analysis
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=1.5, batch_var='pool')
sft_thresh <- c(3,3,3,3,3,3,3)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=1, 
                                      sig_thresh=0.01, percentile_exp_rec=.8,
                                      show_rec_sig=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper/sle_only_lr_f1.pdf", useDingbats = FALSE,
    width = 7.75, height = 7.75)
pbmc_container$plots$lr_analysis[['Factor1']]
dev.off()


ctypes <- c('B','B')
modules <- c(1,2)
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('TF'))
pdf(file = "/home/jmitchel/figures/for_paper/sle_only_lr_tf_mods.pdf", useDingbats = FALSE,
    width = 4, height = 2.75)
mod_enr
dev.off()


lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=1,mod_ct='B',mod=1,lig_ct='cM',lig='TNFSF13B')
pdf(file = "/home/jmitchel/figures/for_paper/sle_only_lig_mod.pdf", useDingbats = FALSE,
    width = 7.5, height = 4.75)
lig_mod_fact
dev.off()

pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=4, 
                                      sig_thresh=0.05, percentile_exp_rec=.99,
                                      show_rec_sig=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper/sle_only_lr_f4.pdf", useDingbats = FALSE,
    width = 7.75, height = 7.75)
pbmc_container$plots$lr_analysis[['Factor4']]
dev.off()







## testing whether any factors are associated with IFN after limiting them to just sle patients
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Status','Age')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL
head(meta)

d_keep <- rownames(meta)[meta$Status=='Managed']
# d_keep <- rownames(meta)

d_exp <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,'MX1']
d_exp <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,'ISG20']
dsc <- pbmc_container$tucker_results[[1]][,1]
# tmp <- cbind.data.frame(dsc[d_keep],d_exp[d_keep])
tmp <- cbind.data.frame(dsc[names(d_exp)],d_exp)
colnames(tmp) <- c('dscore','expres')
plot(tmp$dscore,tmp$expres)
summary(lm(expres~dscore,data=tmp))


trim_names <- sapply(names(d_exp), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(d_exp) <- trim_names

duplicated(trim_names)


## now look at association with sledaiscore...

trim_names <- sapply(names(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(dsc) <- trim_names

colnames(clin_vars)

d_both <- names(dsc)

tmp <- cbind.data.frame(clin_vars[d_both,'sledaiscore'],dsc[d_both])
colnames(tmp) <- c('sledai','dscore')
plot(tmp$dscore,tmp$sledai)
summary(lm(sledai~dscore,data=tmp))



tmp$sledai <- as.factor(tmp$sledai)
m <- polr(sledai ~ ., data = tmp, Hess=TRUE, method='probit')

## view a summary of the model
ctable <- coef(summary(m))
## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
pval <- p[1:10]









# compare similar factors to find out what's similar/different
pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,4), direction=c('up','down'),
                                  compare_type='different', sig_thresh=0.05)
pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,4), direction=c('up','down'),
                                  compare_type='same', sig_thresh=0.05)
pbmc_container$plots$comparisons[['2_4']]
diff_enr <- get_compare_go_enrich(pbmc_container,'cDC',1)
diff_enr <- get_compare_go_enrich(pbmc_container,'ncM',1)
diff_enr <- get_compare_go_enrich(pbmc_container,'NK',1)
diff_enr <- get_compare_go_enrich(pbmc_container,'T8',1)
diff_enr <- get_compare_go_enrich(pbmc_container,'T4',1)
diff_enr <- get_compare_go_enrich(pbmc_container,'B',-1)
print(diff_enr[order(diff_enr,decreasing=F)][1:20])






ctypes <- c('B','NK','cDC','T4','T8','cM','ncM')
modules <- c(6,7,3,4,6,3,2)
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('GO'))
mod_enr



