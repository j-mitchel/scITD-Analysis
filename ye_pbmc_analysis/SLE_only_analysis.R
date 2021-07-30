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
                                                     "Processing_Cohort",
                                                     "Ethnicity"),
                                     metadata_col_nm=c('donors',
                                                       'SLE_status',
                                                       'Status',
                                                       'ctypes',
                                                       'sex',
                                                       'Age',
                                                       'pool',
                                                       'processing',
                                                       'Ethnicity'))


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
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=0,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
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
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(9,20,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(8,23,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(10,40,7),
                                 tucker_type = 'regular', rotation_type = 'ica')

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(12,25,7),
                                 tucker_type = 'sparse', rotation_type = 'ica',
                                 sparsity=sqrt(2))

# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,30,7),
#                                  tucker_type = 'regular', rotation_type = 'varimax')

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,25,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(17,40,7),
                                 tucker_type = 'regular', rotation_type = 'ica')

# get factor-meta data associations
# pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Status'),
#                                         stat_use='rsq')
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Ethnicity'),
                                        stat_use='pval')

# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','processing'),
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
d_exp <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,'IFI6']
dsc <- pbmc_container$tucker_results[[1]][,1]
tmp <- cbind.data.frame(dsc[names(d_exp)],d_exp)
colnames(tmp) <- c('dscore','expres')

# add regression line
lmres <- lm(expres~dscore,data=tmp)
line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

pdf(file = "/home/jmitchel/figures/for_paper/sle_only_f1_IFN.pdf", useDingbats = FALSE,
    width = 4, height = 3.25)
ggplot(tmp,aes(x=dscore,y=expres)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  ylab('IFI6 expression (CD4+ T)') +
  xlab('Factor 1 donor scores') +
  theme_bw()
dev.off()


pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=3, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=3,direc='up',thresh=.05)
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=3)

pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=1, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
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
                                      db_use=c("GO"))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=4,direc='up',thresh=.05)




## factor 3 gsea
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=3, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=3,direc='down',thresh=.05)
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=3)

gsets <- c("GO_B_CELL_RECEPTOR_SIGNALING_PATHWAY",
           "GO_RESPONSE_TO_INTERFERON_GAMMA",
           "GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_REGULATION_OF_RESPONSE_TO_BIOTIC_STIMULUS",
           "GO_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY",
           "GO_INTERLEUKIN_1_PRODUCTION",
           "GO_PROTEIN_POLYUBIQUITINATION",
           "GO_EXOCYTOSIS",
           "GO_CELL_MOTILITY",
           "GO_CELL_CYCLE_G2_M_PHASE_TRANSITION",
           "GO_CELL_ACTIVATION",
           "GO_RESPONSE_TO_STEROID_HORMONE",
           "GO_CELL_GROWTH",
           "GO_CIRCADIAN_RHYTHM",
           "GO_REGULATION_OF_CELL_DEATH",
           "GO_ENZYME_LINKED_RECEPTOR_PROTEIN_SIGNALING_PATHWAY",
           "GO_POSITIVE_REGULATION_OF_DEVELOPMENTAL_PROCESS")

gset_cmap <- rep('black',length(gsets))

names(gset_cmap) <- gsets

pdf(file = "/home/jmitchel/figures/for_paper/sle_only_f3_go.pdf", useDingbats = FALSE,
    width = 9, height = 4)
hm_list <- plot_select_sets(pbmc_container, 3, gsets, color_sets=gset_cmap, 
                            cl_rows=TRUE, h_w=c(8.25,9), myfontsize=7.15)
dev.off()






## running LR analysis with new functions
pbmc_container <- get_lm_pvals(pbmc_container)

# prep for new LR analysis
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
  rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
  return(rname)
})
lr_pairs$interaction_name <- NULL


# pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
#                                    var_scale_power=1.5, batch_var='pool')
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='pool')
# pbmc_container <- prep_LR_interact_v2(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
#                                    var_scale_power=.5, batch_var='pool')
# sft_thresh <- c(3,3,3,2,3,3,3) # for unsigned version
# sft_thresh <- c(10,10,7,5,9,4,7) # for signed version
# sft_thresh <- c(3,3,3,2,3,3,3) # for unsigned version
# sft_thresh <- c(12,14,12,10,12,9,12) # for signed version
sft_thresh <- c(12,14,12,10,12,9,12) # for signed version
# sft_thresh <- c(12,12,12,10,12,9,10) # for signed version
# sft_thresh <- c(12,12,16,8,20,7,10) # for signed version with prep v2
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

# pbmc_container <- compute_LR_interact_v2(pbmc_container, lr_pairs, factor_select=2, 
#                                       sig_thresh=0.05, percentile_exp_rec=.9,
#                                       show_rec_sig=FALSE)
pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=1, 
                                      sig_thresh=0.01, percentile_exp_rec=.85,
                                      show_rec_sig=F)
pbmc_container$plots$lr_analysis[['Factor1']]

ctypes <- c('T4','T4','T4','T4','T8','T8','T8')
modules <- c(3,6,7,10,3,5,6)
ctypes <- c('T4','T4','T4','T8','T8')
modules <- c(3,5,6,4,5)
ctypes <- c('T4','T4','T4')
modules <- c(4,5,6)
ctypes <- c('T8','T8')
modules <- c(4,5)
ctypes <- c('T4')
modules <- c(6)
ctypes <- c('B','B')
modules <- c(1,2)
ctypes <- c('B')
modules <- c(1)
ctypes <- c('T8')
modules <- c(5)

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('TF'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('GO'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('Hallmark'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('BioCarta'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('Reactome'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('KEGG'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('GO','Reactome','KEGG'))
pdf(file = "/home/jmitchel/figures/for_paper/test.pdf", useDingbats = FALSE,
    width = 7, height = 10)
mod_enr
dev.off()

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,mod_ct='T8',mod=4,lig_ct='cM',lig='ICOSLG')
lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,mod_ct='ncM',mod=7,lig_ct='T8',lig='CCL13')
lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,mod_ct='cM',mod=5,lig_ct='T4',lig='TNFSF10')
lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=1,mod_ct='B',mod=1,lig_ct='cM',lig='TNFSF13B')

lig_mod_fact

## looking at expression of cell cycle arrest genes
# first looking into T4_3 module - it is supposedly upregulated in F2
mod_genes <- pbmc_container[["module_genes"]][["T4"]]
mod_genes <- pbmc_container[["module_genes"]][["T8"]]
mod_genes <- mod_genes[mod_genes==2]
mod_genes <- mod_genes[mod_genes==5]
mod_genes <- names(mod_genes)
gset <- my_pathways[['GO_CELL_CYCLE_ARREST']]
gset <- my_pathways[['GO_CELL_CYCLE']]
gset <- my_pathways[['GO_REGULATION_OF_CELL_CYCLE']]
gset <- my_pathways[['GO_NEGATIVE_REGULATION_OF_CELL_CYCLE']]
gset <- my_pathways[['GO_MITOTIC_CELL_CYCLE']]
mod_genes_gset <- mod_genes[mod_genes %in% gset]
mod_genes_gset
mod_genes_gset1=mod_genes_gset
mod_genes_gset2=mod_genes_gset
mod_genes_gset2[!(mod_genes_gset2 %in% mod_genes_gset1)]



# plotting IFN response gene expression against f1
mygene <- 'TNFSF8'
myfactor <- 3
d_exp <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,mygene]
d_exp <- pbmc_container[["scMinimal_ctype"]][['cM']][["pseudobulk"]][,mygene]
# d_exp <- pbmc_container[["scMinimal_ctype"]][['cM']][["pseudobulk"]][,mygene]
dsc <- pbmc_container$tucker_results[[1]][,myfactor]
tmp <- cbind.data.frame(dsc[names(d_exp)],d_exp)
colnames(tmp) <- c('dscore','expres')

# add regression line
lmres <- lm(expres~dscore,data=tmp)
line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

ggplot(tmp,aes(x=dscore,y=expres)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  ylab('IFI6 expression (CD4+ T)') +
  xlab('Factor 1 donor scores') +
  theme_bw()



print(9)











pdf(file = "/home/jmitchel/figures/for_paper/sle_only_lr_f1.pdf", useDingbats = FALSE,
    width = 7.75, height = 7.75)
pbmc_container$plots$lr_analysis[['Factor1']]
dev.off()


ctypes <- c('B','B')
modules <- c(1,7)
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('TF'))
# pdf(file = "/home/jmitchel/figures/for_paper/sle_only_lr_tf_mods.pdf", useDingbats = FALSE,
#     width = 4, height = 2.75)
mod_enr
dev.off()

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.1, db_use=c('GO'))
pdf(file = "/home/jmitchel/figures/for_paper/sle_only_lr_go_mods.pdf", useDingbats = FALSE,
    width = 6, height = 8)
mod_enr
dev.off()


lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=1,mod_ct='B',mod=1,lig_ct='cM',lig='TNFSF13B')
pdf(file = "/home/jmitchel/figures/for_paper/sle_only_lig_mod.pdf", useDingbats = FALSE,
    width = 7.5, height = 4.75)
lig_mod_fact
dev.off()


lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=1,mod_ct='B',mod=7,lig_ct='cM',lig='TNFSF13B')
pdf(file = "/home/jmitchel/figures/for_paper/sle_only_lig_mod_B7.pdf", useDingbats = FALSE,
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




# LR analysis for f3
pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=3, 
                                      sig_thresh=0.01, percentile_exp_rec=.90,
                                      show_rec_sig=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper/sle_only_lr_f3.pdf", useDingbats = FALSE,
    width = 7.75, height = 7.75)
pbmc_container$plots$lr_analysis[['Factor3']]
dev.off()


ctypes <- c('T4','T8')
modules <- c(3,2)
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('GO'))
pdf(file = "/home/jmitchel/figures/for_paper/sle_only_lr_f3_go_mods.pdf", useDingbats = FALSE,
    width = 5.5, height = 5)
mod_enr
dev.off()

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=3,mod_ct='T8',mod=2,lig_ct='cM',lig='ICOSLG')
pdf(file = "/home/jmitchel/figures/for_paper/sle_only_lr_f3_lig_mod.pdf", useDingbats = FALSE,
    width = 7.5, height = 4.75)
lig_mod_fact
dev.off()


# seeing how well L + R predicts module and if interaction does even better
dsc <- container$tucker_results[[1]][,3]
lig_exp <- container$scale_pb_extra[["T8"]][,"PTPRC"]
rec_exp <- container$scale_pb_extra[["B"]][,"CD22"]
MEs <- container[["module_eigengenes"]][["B"]]
ME <- MEs[,5]
names(ME) <- rownames(MEs)
tmp <- as.data.frame(cbind(dsc[names(ME)],ME,lig_exp[names(ME)],rec_exp[names(ME)]))
colnames(tmp) <- c('dsc','ME','lig_expr','rec_expr')

lm1 = lm(ME ~ lig_expr, data=tmp)
lm2 = lm(ME ~ lig_expr + rec_expr, data=tmp)
anova_res <- anova(lm1,lm2)
anova_pval <- anova_res$`Pr(>F)`[2]
anova_pval







# showing factor 1 association with IFN
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



## testing if those on prednisone prioritized by f2 are different from those prioritized by f4
tmp <- cbind.data.frame(dsc[,2],clin_vars[,'prednisone'])
colnames(tmp) <- c('dscore','cvar')
tmp <- tmp[tmp$cvar==1,]
tmp <- tmp[tmp$dscore>.05,]
f2_pred <- rownames(tmp)

tmp <- cbind.data.frame(dsc[,4],clin_vars[,'prednisone'])
colnames(tmp) <- c('dscore','cvar')
tmp <- tmp[tmp$cvar==1,]
tmp <- tmp[tmp$dscore<(-.05),]
f4_pred <- rownames(tmp)

intersect(f2_pred,f4_pred)
union(f2_pred,f4_pred)

print(f2_pred)
print(f4_pred)


f4_only <- f4_pred[!(f4_pred %in% f2_pred)]

# does f4_only have lower disease activity or one of the other symptoms compared to f2 pred donors?
all_pvals <- c()
for (cvar in colnames(clin_vars)) {
  print(cvar)

  # do hypergeometric test
  group2 <- sum(clin_vars[f4_only,cvar])
  group1 <- sum(clin_vars[f2_pred,cvar])
  
  total <- length(c(f4_only,f2_pred))
  
  pval <- phyper(group1-1, group1+group2, total-(group1+group2), length(f2_pred), lower.tail= FALSE)
  all_pvals <- c(all_pvals,pval)
}
all_pvals <- p.adjust(all_pvals,method='fdr')

all_pvals <- all_pvals[!is.na(all_pvals)]





## DE testing to see if there are differences between donors with different symptoms
colnames(clin_vars)
myvar <- 'crflupusneph' # ones tested: crfmucuulcers, crflupusneph, sledairash
meta <- clin_vars[,myvar,drop=FALSE]
pb <- pbmc_container[["scMinimal_ctype"]][["T4"]][["pseudobulk"]]
trim_names <- sapply(rownames(pb), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
rownames(pb) <- trim_names


pvals <- c()
for (j in 1:ncol(pb)) {
  tmp <- cbind.data.frame(meta,pb[rownames(meta),j])
  colnames(tmp) <- c('cvar','expr')
  
  tres <- t.test(tmp[tmp$cvar==1,2],tmp[tmp$cvar==0,2])
  pv <- tres$p.value
  pvals <- c(pvals,pv)
}
pvals <- p.adjust(pvals,method='fdr')
pvals[order(pvals)][1:10]














## seeing if factor 1 score is more highly associated with anti-dsDNA compared to IFN genes alone
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"acrantidsdna"])]
tmp <- as.data.frame(cbind(dsc[d_keep,1], clin_vars[d_keep,"acrantidsdna"]))
tmp <- as.data.frame(cbind(dsc[d_keep,1], clin_vars[d_keep,"sledaiscore"]))
colnames(tmp) <- c('dscore','cvar')
# tmp$cvar <- as.factor(tmp$cvar)
lmres <- lm(dscore~cvar,tmp)
summary(lmres)

# now trying with expression of IFI6 in T4
pb <- pbmc_container[["scMinimal_ctype"]][["T4"]][["pseudobulk"]][,'MX1']
pb <- pbmc_container[["scMinimal_ctype"]][["T4"]][["pseudobulk"]][,'IFI6']
pb <- pbmc_container[["scMinimal_ctype"]][["T4"]][["pseudobulk"]][,'MX2']
trim_names <- sapply(names(pb), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
names(pb) <- trim_names

tmp <- as.data.frame(cbind(pb[d_keep], clin_vars[d_keep,"acrantidsdna"]))
tmp <- as.data.frame(cbind(pb[d_keep], clin_vars[d_keep,"sledaiscore"]))
colnames(tmp) <- c('dscore','cvar')
# tmp$cvar <- as.factor(tmp$cvar)
lmres <- lm(dscore~cvar,tmp)
summary(lmres)

fmod <- glm(cvar~dscore, data=tmp, family = "binomial") ##"full" mod
nmod <- glm(cvar~1, data=tmp, family = 'binomial') ##"null" mod
a_res <- anova(nmod, fmod, test = 'Chisq')
pval <- a_res$`Pr(>Chi)`[2]
pval










