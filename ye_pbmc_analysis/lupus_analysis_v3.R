

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


# form tensor with batch correction applied
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                              scale_var = TRUE, var_scale_power = 1.5,
                              batch_var='pool')


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(10,20,7),
                                 tucker_type = 'regular', rotation_type = 'ica')


# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Status'),
                                        stat_use='pval')

# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# pdf(file = "/home/jmitchel/figures/for_paper/lupus_dscores_v2.pdf", useDingbats = FALSE,
#     width = 6, height = 7)
pbmc_container$plots$donor_matrix
# dev.off()


# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(10,20,7), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')

# saveRDS(pbmc_container[["gene_score_associations"]],file='/home/jmitchel/data/lupus_data/lupus_jackstraw.rds')
pbmc_container[["gene_score_associations"]] <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_jackstraw.rds')


# run gsea for a f4
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=4, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=4,direc='down',thresh=.05)

# investigate several clusters of enriched sets
plot_gsea_sub(pbmc_container,factor_select=4,direc='down',thresh=.05,clust_select=1)
plot_gsea_sub(pbmc_container,factor_select=4,direc='down',thresh=.05,clust_select=2)
plot_gsea_sub(pbmc_container,factor_select=4,direc='down',thresh=.05,clust_select=3)
plot_gsea_sub(pbmc_container,factor_select=4,direc='down',thresh=.05,clust_select=4)
plot_gsea_sub(pbmc_container,factor_select=4,direc='down',thresh=.05,clust_select=5)
plot_gsea_sub(pbmc_container,factor_select=4,direc='down',thresh=.05,clust_select=6)
plot_gsea_sub(pbmc_container,factor_select=4,direc='down',thresh=.05,clust_select=7)
plot_gsea_sub(pbmc_container,factor_select=4,direc='down',thresh=.05,clust_select=8)
plot_gsea_sub(pbmc_container,factor_select=4,direc='down',thresh=.05,clust_select=9)

plot_gsea_hmap_w_similarity(pbmc_container,factor_select=4,direc='up',thresh=.05)
plot_gsea_sub(pbmc_container,factor_select=4,direc='up',thresh=.05,clust_select=1)
plot_gsea_sub(pbmc_container,factor_select=4,direc='up',thresh=.05,clust_select=2)
plot_gsea_sub(pbmc_container,factor_select=4,direc='up',thresh=.05,clust_select=3)
plot_gsea_sub(pbmc_container,factor_select=4,direc='up',thresh=.05,clust_select=4)


## f4 sets to show on loading hmap
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_OSSIFICATION",
           "GO_CELLULAR_EXTRAVASATION",
           'GO_MYELOID_LEUKOCYTE_ACTIVATION',
           "GO_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY",
           'GO_LYMPHOCYTE_MIGRATION',
           'GO_RESPONSE_TO_TUMOR_NECROSIS_FACTOR',
           'GO_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION',
           'GO_NEGATIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY',
           'GO_SECRETION',
           'GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I',
           'GO_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY',
           "GO_INTERLEUKIN_4_PRODUCTION")

gset_cmap <- c('blue',
               'orange',
               'forest green',
               'purple',
               'black',
               'black',
               'black',
               'black',
               'black',
               'black',
               'black',
               'black',
               'black')

names(gset_cmap) <- gsets

gset_cmap_sub <- gset_cmap[gset_cmap!='black']
gset_sub <- names(gset_cmap_sub)

dev.off()
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=4, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gset_sub, le_set_colormap=gset_cmap_sub, le_set_num_per=6, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)
dev.off()
hm_list <- plot_select_sets(pbmc_container, 4, gsets, thresh=.05, color_sets=gset_cmap, 
                            cl_rows=TRUE)

p1 <- pbmc_container$plots$all_lds_plots[['4']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["4"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f4_lds_go.pdf", useDingbats = FALSE,
    width = 7.5, height = 8.5)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()


# plot dscore status violin plot for factor 4
dscores <- pbmc_container[["tucker_results"]][[1]]
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Status')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

tmp <- cbind.data.frame(dscores[,4],meta[rownames(dscores),1])
colnames(tmp) <- c('dsc','Status')
f4_status_plot <- ggplot(tmp,aes(x=Status,y=dsc)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.65, binwidth = .01) +
  ylab('Factor 4 Donor Score') +
  xlab('SLE Status') +
  geom_hline(yintercept = 0, linetype="dashed", 
               color = "gray", size=1.5) +
  coord_flip() +
  theme_bw()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f4_status.pdf", useDingbats = FALSE,
    width = 4.5, height = 2)
f4_status_plot
dev.off()

# get the pvalue for status difference
pbmc_container[["meta_associations"]]['Status','Factor4']


# run gsea for a f7
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=7, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=7,direc='up',thresh=.05)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=7,direc='down',thresh=.05)

# investigate several clusters of enriched sets
plot_gsea_sub(pbmc_container,factor_select=7,direc='up',thresh=.05,clust_select=1)
plot_gsea_sub(pbmc_container,factor_select=7,direc='up',thresh=.05,clust_select=2)
plot_gsea_sub(pbmc_container,factor_select=7,direc='up',thresh=.05,clust_select=3)
plot_gsea_sub(pbmc_container,factor_select=7,direc='up',thresh=.05,clust_select=4)



## f7 sets to show on loading hmap
gsets <- c("GO_IMMUNE_EFFECTOR_PROCESS",
           "GO_INTERLEUKIN_8_PRODUCTION",
           "GO_CELL_CELL_ADHESION",
           "GO_LYMPHOCYTE_MIGRATION",
           "GO_REGULATION_OF_CELL_DEVELOPMENT",
           "GO_EXOCYTOSIS",
           "GO_CELL_KILLING")

gset_cmap <- c('blue',
               'orange',
               'forest green',
               'purple',
               'black',
               'black',
               'black')

names(gset_cmap) <- gsets

gset_cmap_sub <- gset_cmap[gset_cmap!='black']
gset_sub <- names(gset_cmap_sub)

dev.off()
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=7, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gset_sub, le_set_colormap=gset_cmap_sub, le_set_num_per=6, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)
dev.off()
hm_list <- plot_select_sets(pbmc_container, 7, gsets, thresh=.05, color_sets=gset_cmap, 
                            cl_rows=TRUE)

p1 <- pbmc_container$plots$all_lds_plots[['7']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["7"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f7_lds_go.pdf", useDingbats = FALSE,
    width = 7.5, height = 8.5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f7_lds_go.pdf", useDingbats = FALSE,
    width = 9.5, height = 12.25)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()


# save gsea results
# saveRDS(pbmc_container[["gsea_res_full"]],file='/home/jmitchel/data/lupus_data/lupus_v2_gsea_f4_f7_res_full.rds')
# saveRDS(pbmc_container[["gsea_results"]],file='/home/jmitchel/data/lupus_data/lupus_v2_gsea_f4_f7_results.rds')

# load gsea results
pbmc_container[["gsea_res_full"]] <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_v2_gsea_f4_f7_res_full.rds')
pbmc_container[["gsea_results"]] <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_v2_gsea_f4_f7_results.rds')




# run gsea for a f2
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=2,direc='up',thresh=.05)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=2,direc='down',thresh=.05)

# investigate several clusters of enriched sets
plot_gsea_sub(pbmc_container,factor_select=2,direc='up',thresh=.05,clust_select=1)
plot_gsea_sub(pbmc_container,factor_select=2,direc='up',thresh=.05,clust_select=2)
plot_gsea_sub(pbmc_container,factor_select=2,direc='up',thresh=.05,clust_select=3)
plot_gsea_sub(pbmc_container,factor_select=2,direc='up',thresh=.05,clust_select=4)
plot_gsea_sub(pbmc_container,factor_select=2,direc='up',thresh=.05,clust_select=5)
plot_gsea_sub(pbmc_container,factor_select=2,direc='up',thresh=.05,clust_select=6)
plot_gsea_sub(pbmc_container,factor_select=2,direc='up',thresh=.05,clust_select=7)
plot_gsea_sub(pbmc_container,factor_select=2,direc='up',thresh=.05,clust_select=8)
plot_gsea_sub(pbmc_container,factor_select=2,direc='up',thresh=.05,clust_select=9)
plot_gsea_sub(pbmc_container,factor_select=2,direc='up',thresh=.05,clust_select=10)
plot_gsea_sub(pbmc_container,factor_select=2,direc='up',thresh=.05,clust_select=11)
plot_gsea_sub(pbmc_container,factor_select=2,direc='up',thresh=.05,clust_select=17)

plot_gsea_sub(pbmc_container,factor_select=2,direc='down',thresh=.05,clust_select=1)
plot_gsea_sub(pbmc_container,factor_select=2,direc='down',thresh=.05,clust_select=2)


## f2 sets to show on loading hmap
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_ENDOTHELIAL_CELL_APOPTOTIC_PROCESS",
           "GO_CELLULAR_KETONE_METABOLIC_PROCESS",
           "GO_ERBB_SIGNALING_PATHWAY",
           "GO_APOPTOTIC_PROCESS",
           "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
           "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II",
           "GO_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE",
           "GO_NEGATIVE_REGULATION_OF_CELL_DEATH",
           "GO_INTERLEUKIN_1_MEDIATED_SIGNALING_PATHWAY",
           "GO_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
           "GO_POSITIVE_REGULATION_OF_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING",
           "GO_NEGATIVE_REGULATION_OF_CELL_DIFFERENTIATION",
           "GO_NEGATIVE_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION",
           "GO_NUCLEIC_ACID_PHOSPHODIESTER_BOND_HYDROLYSIS",
           "GO_INTERLEUKIN_2_PRODUCTION",
           "GO_INTERLEUKIN_4_PRODUCTION",
           "GO_INTERLEUKIN_1_BETA_PRODUCTION",
           "GO_TYPE_I_INTERFERON_PRODUCTION",
           "GO_SECRETION")

gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_ENDOTHELIAL_CELL_APOPTOTIC_PROCESS",
           "GO_CELLULAR_KETONE_METABOLIC_PROCESS",
           "GO_ERBB_SIGNALING_PATHWAY",
           "GO_APOPTOTIC_PROCESS",
           "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
           "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II",
           "GO_ACTIVATION_OF_INNATE_IMMUNE_RESPONSE",
           "GO_NEGATIVE_REGULATION_OF_CELL_DEATH",
           "GO_INTERLEUKIN_1_MEDIATED_SIGNALING_PATHWAY",
           "GO_RESPONSE_TO_TUMOR_NECROSIS_FACTOR",
           "GO_POSITIVE_REGULATION_OF_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING",
           "GO_NEGATIVE_REGULATION_OF_CELL_DIFFERENTIATION",
           "GO_NEGATIVE_REGULATION_OF_CELL_CYCLE_G2_M_PHASE_TRANSITION",
           "GO_NUCLEIC_ACID_PHOSPHODIESTER_BOND_HYDROLYSIS",
           "GO_INTERLEUKIN_2_PRODUCTION",
           "GO_INTERLEUKIN_4_PRODUCTION",
           "GO_INTERLEUKIN_1_BETA_PRODUCTION",
           "GO_RESPONSE_TO_INTERFERON_GAMMA",
           "GO_SECRETION")

gset_cmap <- c('blue',
               'orange',
               'forest green',
               'purple',
               rep('black',16))

names(gset_cmap) <- gsets

gset_cmap_sub <- gset_cmap[gset_cmap!='black']
gset_sub <- names(gset_cmap_sub)

dev.off()
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gset_sub, le_set_colormap=gset_cmap_sub, le_set_num_per=6, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)
dev.off()
hm_list <- plot_select_sets(pbmc_container, 2, gsets, thresh=.05, color_sets=gset_cmap, 
                            cl_rows=TRUE)

p1 <- pbmc_container$plots$all_lds_plots[['2']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["2"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f2_lds_go_V2.pdf", useDingbats = FALSE,
    width = 9.5, height = 12.25)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()

# plot dscore status violin plot for factor 2
dscores <- pbmc_container[["tucker_results"]][[1]]
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Status')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

tmp <- cbind.data.frame(dscores[,2],meta[rownames(dscores),1])
colnames(tmp) <- c('dsc','Status')
f2_status_plot <- ggplot(tmp,aes(x=Status,y=dsc)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, binwidth = .01) +
  ylab('Factor 2 Donor Score') +
  xlab('SLE Status') +
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "gray", size=1.5) +
  coord_flip() +
  theme_bw()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f2_status.pdf", useDingbats = FALSE,
    width = 5, height = 3)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_f2_status.pdf", useDingbats = FALSE,
    width = 4.5, height = 3.5)
f2_status_plot
dev.off()






## get loadings and status plot for f5
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=5, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=5,direc='up',thresh=.05)
plot_gsea_sub(pbmc_container,factor_select=5,direc='up',thresh=.05,clust_select=1)
plot_gsea_sub(pbmc_container,factor_select=5,direc='up',thresh=.05,clust_select=2)
plot_gsea_sub(pbmc_container,factor_select=5,direc='up',thresh=.05,clust_select=3)
plot_gsea_sub(pbmc_container,factor_select=5,direc='up',thresh=.05,clust_select=4)

plot_gsea_hmap_w_similarity(pbmc_container,factor_select=5,direc='down',thresh=.05)
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=1)
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=2)
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=3)
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=4)
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=5)
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=6)
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=7)
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=8)
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=9)
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=10)
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=11)
plot_gsea_sub(pbmc_container,factor_select=5,direc='down',thresh=.05,clust_select=20)

## f5 sets to show on loading hmap
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
           "GO_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
           "GO_POSITIVE_REGULATION_OF_MEMBRANE_PERMEABILITY",
           "GO_RESPONSE_TO_CORTICOSTEROID",
           "GO_RESPONSE_TO_LIPID",
           "GO_RESPONSE_TO_HORMONE",
           "GO_NEGATIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS",
           "GO_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING",
           "GO_APOPTOTIC_PROCESS",
           "GO_TYPE_I_INTERFERON_PRODUCTION",
           "GO_NEGATIVE_REGULATION_OF_TYPE_I_INTERFERON_PRODUCTION",
           "GO_POSITIVE_REGULATION_OF_PROTEIN_KINASE_ACTIVITY",
           "GO_SECRETION",
           "GO_REGULATION_OF_CELL_POPULATION_PROLIFERATION",
           "GO_REGULATION_OF_CELL_ADHESION")

gset_cmap <- c('blue',
               'brown',
               'forest green',
               'purple',
               rep('black',length(gsets)-4))

names(gset_cmap) <- gsets

gset_cmap_sub <- gset_cmap[gset_cmap!='black']
gset_sub <- names(gset_cmap_sub)

dev.off()
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=5, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gset_sub, le_set_colormap=gset_cmap_sub, le_set_num_per=6, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)
dev.off()
hm_list <- plot_select_sets(pbmc_container, 5, gsets, thresh=.05, color_sets=gset_cmap, 
                            cl_rows=TRUE)

p1 <- pbmc_container$plots$all_lds_plots[['5']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["5"]]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f5_lds_go.pdf", useDingbats = FALSE,
    width = 10, height = 11.5)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()

# plot dscore status violin plot for factor 5
dscores <- pbmc_container[["tucker_results"]][[1]]
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Status')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

tmp <- cbind.data.frame(dscores[,5],meta[rownames(dscores),1])
colnames(tmp) <- c('dsc','Status')
f5_status_plot <- ggplot(tmp,aes(x=Status,y=dsc)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, binwidth = .01) +
  ylab('Factor 5 Donor Score') +
  xlab('SLE Status') +
  geom_hline(yintercept = 0, linetype="dashed", 
             color = "gray", size=1.5) +
  coord_flip() +
  theme_bw()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f5_status.pdf", useDingbats = FALSE,
    width = 4.5, height = 3.5)
f5_status_plot
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


# factor 5 LR analysis
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=1.5, batch_var='pool')
sft_thresh <- c(3,3,2,2,2,2,2)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=5, 
                                      sig_thresh=0.05, percentile_exp_rec=.9,
                                      show_rec_sig=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_LR_f5.pdf", useDingbats = FALSE,
    width = 7.75, height = 7.75)
pbmc_container$plots$lr_analysis[['Factor5']]
dev.off()


# getting GO enrichment HMAPs for modules
ctypes <- c('cM','T4','cDC','ncM','NK','B')
modules <- c(3,5,3,2,8,6)

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use='TF')
pdf(file = "/home/jmitchel/figures/for_paper/lupus_LR_f5_TF.pdf", useDingbats = FALSE,
    width = 6, height = 3)
mod_enr
dev.off()

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('GO'))
pdf(file = "/home/jmitchel/figures/for_paper/lupus_LR_f5_GO.pdf", useDingbats = FALSE,
    width = 6, height = 9)
mod_enr
dev.off()


lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=5,mod_ct='cM',mod=3,lig_ct='T4',lig='TNFSF8')
pdf(file = "/home/jmitchel/figures/for_paper/lupus_LR_lig_mod_fact.pdf", useDingbats = FALSE,
    width = 5.5, height = 4.5)
lig_mod_fact
dev.off()

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('Hallmark'))
pdf(file = "/home/jmitchel/figures/for_paper/lupus_LR_f5_Hallmark.pdf", useDingbats = FALSE,
    width = 4.5, height = 3)
mod_enr
dev.off()


# factor 4 LR analysis
pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=4, 
                                      sig_thresh=0.05, percentile_exp_rec=.9,
                                      show_rec_sig=T)

pbmc_container$plots$lr_analysis[['Factor4']]

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=4,mod_ct='T8',mod=2,lig_ct='cM',lig='ICOSLG')
pdf(file = "/home/jmitchel/figures/for_paper/lupus_LR__f4_lig_mod_fact.pdf", useDingbats = FALSE,
    width = 5.5, height = 4.5)
lig_mod_fact
dev.off()


# get enriched gene sets as additional evidence of the ICOSLG find
ctypes <- c('T8','T4')
modules <- c(2,3)

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use='TF')
pdf(file = "/home/jmitchel/figures/for_paper/lupus_LR_f4_TF.pdf", useDingbats = FALSE,
    width = 4, height = 6)
mod_enr
dev.off()

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('GO'))
pdf(file = "/home/jmitchel/figures/for_paper/lupus_LR_f4_GO.pdf", useDingbats = FALSE,
    width = 5.25, height = 3.75)
mod_enr
dev.off()

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('BioCarta'))
pdf(file = "/home/jmitchel/figures/for_paper/lupus_LR_f4_BioCarta.pdf", useDingbats = FALSE,
    width = 4, height = 2)
mod_enr
dev.off()



# investigate whether ICOS module is localized to a subtype
pbmc_container$subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')
subclusts <- pbmc_container[["subclusters"]][["T8"]][["res:0.6"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('T8','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg

# set up project parameters
param_list <- initialize_params(ctypes_use = c("T8_1","T8_2","T8_3"),
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
# dsc <- pbmc_container$tucker_results[[1]][,4]
d_included <- rownames(pbmc_sub_container$scMinimal_ctype[[1]]$pseudobulk)
# dsc2 <- dsc[d_included]

# trying to use ligands expression instead
lig_exp <- pbmc_container$scMinimal_ctype[['cM']]$pseudobulk[,'ICOSLG']
lig_exp2 <- lig_exp[d_included]

gn_test <- c('CNOT1', 'DYNLT1', 'SBDS', 'CETN3')
sub_test <- c("T8_1","T8_2","T8_3")
myres <- as.data.frame(matrix(ncol=4,nrow=3))
colnames(myres) <- gn_test
rownames(myres) <- sub_test
for (sub in sub_test) {
  for (g in gn_test) {
    pb <- pbmc_sub_container$scMinimal_ctype[[sub]]$pseudobulk[names(lig_exp2),g]
    mycor1 <- cor(lig_exp2,pb)
    myres[sub,g] <- mycor1
  }
}

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
subtype_hmap <- Heatmap(as.matrix(myres), name = 'lig-target cor',
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        col=col_fun,
        row_names_side='left',
        column_names_side='top',
        column_names_rot = 40,
        row_title = 'T8 subtype',
        border=TRUE,
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(round(myres[i,j],2), x, y)
        })

pdf(file = "/home/jmitchel/figures/for_paper/lupus_LR_f4_subtype_hmap.pdf", useDingbats = FALSE,
    width = 5, height = 3.75)
subtype_hmap
dev.off()




# # look at MIF association with f5 specific modules
# plot_mod_and_lig(pbmc_container,factor_select=5,mod_ct='cM',mod=3,lig_ct='T8',lig='MIF')
# 
# # my guess is that it correlates even better with those lupus donors that have had highest dose of prednisone
# library("readxl")
# clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_meds_cleaned.xlsx')
# clin_vars <- as.data.frame(clin_vars)
# rownames(clin_vars) <- clin_vars[,'Sample ID']
# clin_vars[,'Sample ID'] <- NULL
# 
# # make all NA into zeros
# clin_vars[is.na(clin_vars)] <- 0
# 
# # separate out pred dose as it's the only continuous variable here
# pred_dose <- clin_vars[,'pred_dose',drop=FALSE]
# 
# mif_exp <- pbmc_container$scMinimal_ctype[['T8']]$pseudobulk[,'MIF']
# trim_names <- sapply(names(mif_exp), function(x) {
#   strsplit(x,split='_')[[1]][[1]]
# })
# names(mif_exp) <- trim_names
# 
# names_both <- intersect(rownames(pred_dose),names(mif_exp))
# tmp <- cbind.data.frame(mif_exp[names_both],pred_dose[names_both,])
# summary(lm(tmp[,1]~tmp[,2]))
# plot(tmp[,1],tmp[,2])
# 
# 
# # trying to limit donors to those on pred. See if get inverse trend with mif vs dscore
# clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_meds_cleaned.xlsx')
# clin_vars <- as.data.frame(clin_vars)
# rownames(clin_vars) <- clin_vars[,'Sample ID']
# clin_vars[,'Sample ID'] <- NULL
# clin_vars[is.na(clin_vars)] <- 0
# clin_vars <- clin_vars[names_both,]
# 
# dsc <- pbmc_container$tucker_results[[1]][,5]
# trim_names <- sapply(names(dsc), function(x) {
#   strsplit(x,split='_')[[1]][[1]]
# })
# names(dsc) <- trim_names
# tmp <- cbind.data.frame(mif_exp[names_both],dsc[names_both],clin_vars[,'prednisone'])
# 
# tmp <- tmp[tmp[,3]==1,]
# 
# summary(lm(tmp[,1]~tmp[,2]))
# plot(tmp[,2],tmp[,1])
# 
# MEs <- pbmc_container[["module_eigengenes"]][['cM']]
# ME <- MEs[,3]
# names(ME) <- rownames(MEs)
# names(ME) <- sapply(names(ME), function(x) {
#   strsplit(x,split='_')[[1]][[1]]
# })
# tmp <- cbind.data.frame(mif_exp[names_both],ME[names_both],clin_vars[,'prednisone'])
# tmp <- tmp[tmp[,3]==1,]
# plot(tmp[,2],tmp[,1])





















## now for cell subtype proportion analysis
# add conos object for cell proportion analysis
con <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_conos.rds')
pbmc_container$embedding <- con
rm(con)
gc()

# in case below fn errors in the middle, need to save these objects
orig_embed <- pbmc_container$embedding[["embedding"]]
orig_clusts <- pbmc_container$embedding$clusters$leiden$groups

# # to recover original embedding/cell assignments
# pbmc_container$embedding[["embedding"]] <- orig_embed
# pbmc_container$embedding$clusters$leiden$groups <- orig_clusts


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

# previously saved the de hmaps at this point
pbmc_container[["plots"]][["subtype_de"]] <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subc_de_plots.rds')










## render cell subtype DE heatmaps
## using seurat for to make them as dotplots

# T8
subclusts <- pbmc_container[["subclusters"]][["T8"]][["res:0.6"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('T8','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('GZMH','FGFBP2','GNLY','GZMB','NKG7',
            'GZMA','CCL5','CST7','LGALS1','KLRD1',
            'DUSP2','GZMK','LYAR','CMC1','PIK3R1',
            'CD69','ZFP36L2','CD8B','NOSIP','CCR7',
            'EIF3E','C6orf48','PRKCQ-AS1','LTB','LDHB',
            'SELL')
t8_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()

# T4
subclusts <- pbmc_container[["subclusters"]][["T4"]][["res:0.6"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('T4','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('CCR7','FHIT','AIF1','NUCB2','SELL','PRKCQ-AS1','LINC00861','CD7',
            'RGS10','LEF1','ITGB1','S100A4','ANXA1','S100A10','S100A11','LGALS1',
            'KLRB1','EMP3','SH3BGRL3','CRIP1','JUN','JUNB','CD69','FOS','DUSP1',
            'CXCR4','TSC22D3','LEPROTL1','ZFP36L2','PNRC1','IFI6','ISG15','IFI44L',
            'LY6E','MT2A','MX1','IFIT3','IFITM1','EIF2AK2','ISG20','NEAT1',
            'JUND','MT-ND4L','SYNE2','ANKRD11','ETS1','DDX17','NKTR','MT-ND5','C1QA',
            'C1QB','COX5B','PLAC8','PSAP')
t4_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()


# NK
subclusts <- pbmc_container[["subclusters"]][["NK"]][["res:0.6"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('NK','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('FGFBP2','SPON2','GZMH','FCGR3A','LGALS1','PRF1','LAIR2','CCL4',
            'GZMB','IGFBP7','CMC1','CXCR4','DUSP2','CD160','PIK3R1','XCL2',
            'MAP3K8','CCL3','ZFP36','JUNB','GZMK','XCL1','SELL','KLRC1','CD44',
            'NFKBIA','COTL1','C1orf56','CDC42SE1','HNRNPH1','APOBEC3C','MDM4',
            'B4GALT1','CDC42','CTNNB1','SAR1A','SET')
nk_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()



# cM
subclusts <- pbmc_container[["subclusters"]][["cM"]][["res:0.5"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('cM','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('IFI6','ISG15','LY6E','IFI44L','IFITM3','IFI44','MX1','MT2A',
            'S100A12','EPSTI1','HLA-DPB1','HLA-DPA1','HLA-DMA','HLA-DQB1',
            'HLA-DQA1','LGALS2','CPVL','HLA-DRB1','SLC25A5','EEF1B2','ALOX5AP',
            'MGST1','RBP7','VCAN','PLBD1','CDA','METTL9','RETN','VNN2',
            'APOBEC3A','MARCKS','PSME2','IL8','IL1B','G0S2','CCL3','EREG',
            'TNFAIP3','NFKBIZ','SOD2','IER3','NAMPT')
cM_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()


# ncM
subclusts <- pbmc_container[["subclusters"]][["ncM"]][["res:0.6"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('ncM','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('TESC','VMO1','LYPD2','PPM1N','CKB','CDKN1C','ICAM4','SOD1','MEG3',
            'IFI6','ISG15','IFI44L','APOBEC3A','EPSTI1','IFIT3','TNFSF10',
            'PLSCR1','PLAC8','MX1','C1QA','C1QB','C1QC','HLA-DQA1','HLA-DMA',
            'VAMP8','HLA-DQB1','VAMP5','HLA-DMB','MS4A6A','JUND','LGALS2','G0S2',
            'VCAN','GPX1','MT-ND4L')
ncM_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()



# cDC
subclusts <- pbmc_container[["subclusters"]][["cDC"]][["res:0.5"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('cDC','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('VCAN','FCN1','S100A12','CD14','RNASE2','MS4A6A','CFD','CD36',
            'S100A8','SERPINA1','CD1C','FCER1A','ENHO','IL2RG','NDRG2','BASP1',
            'FCGR2B','CD1E','ARL4C','CLEC10A','CLEC9A','C1orf54','IRF8',
            'DNASE1L3','CPNE3','WDFY4','HLA-DOB','IDO1','GYPC')
cDC_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()


# B
subclusts <- pbmc_container[["subclusters"]][["B"]][["res:0.8"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('B','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('GPR183','CRIP1','COTL1','EMP3','CLECL1','LGALS1','GAPDH','S100A6',
            'S100A10','IL4R','PPAPDC1B','FCER2','MEF2C','ADAM28','TCL1A','BIRC3',
            'STAG3','HVCN1','VPREB3','CD79B','IGLL5','SNX29P2','FCRLA','PPP1R14A',
            'MZB1','CD72','IER2','CD69','FOS','FOSB','JUN','NFKBIA','ZFP36','JUNB',
            'YPEL5','CD83','C1orf56','HNRNPH1','CDC42SE1','MDM4','B4GALT1',
            'PPP3CA','APOBEC3C','CDC42','CTNNB1','FGD2')
b_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()


all_dot1 <- list(b_dot,cDC_dot,cM_dot,ncM_dot)
all_dot2 <- list(NULL,nk_dot,t4_dot,t8_dot,NULL)
pg1 <- plot_grid(plotlist = all_dot1,nrow=1)
pg2 <- plot_grid(plotlist = all_dot2,nrow=1,rel_widths = c(.2,.5,.5,.5,.2))
pg3 <- plot_grid(pg1,pg2,nrow=2)


pdf(file = "/home/jmitchel/figures/for_paper/lupus_all_subc_dot.pdf", useDingbats = FALSE,
    width = 25, height = 15)
pg3
dev.off()


# needed to make t4 plot taller so resaving it here
pdf(file = "/home/jmitchel/figures/for_paper/lupus_t4_dot.pdf", useDingbats = FALSE,
    width = 7, height = 11)
t4_dot
dev.off()





# create main embedding for lupus data
library(RColorBrewer)
mycolors <- brewer.pal(n = 9, name = "Set1")
mycolors <- mycolors[c(1:5,7,9)]
tmp <- pbmc_container$embedding$plotGraph(alpha=0.1)
tmp <- tmp + theme(panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black", size=1, fill=NA)) +
  scale_color_manual(values=mycolors)
  # scale_color_brewer(palette="Set1")
tmp$layers[[2]] <- NULL


# pdf(file = "/home/jmitchel/figures/for_paper/lupus_embedding.pdf", useDingbats = FALSE,
#     width = 7, height = 7)
# saved a jpeg 550 x 400 dimensions
tmp
dev.off()




# make barplots and dotplots showing associations of subtype proportions with factors

# CD4 T cell factor 4 associations
pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='T4',factor_use=4)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_f4_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['T4']]
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'T4',0.6,subtype=1,factor_use=4)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD4_f4_sub1_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()


# B cell factor 4 associations
pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='B',factor_use=4)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/B_f4_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['B']]
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'B',0.8,subtype=3,factor_use=4)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/B_f4_sub3_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()


# B cell factor 2 associations
pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='B',factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/B_f2_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['B']]
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'B',0.8,subtype=4,factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/B_f2_sub4_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()


# cM factor 2 associations
pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='cM',factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/cM_f2_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['cM']]
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


# ncM factor 2 associations
pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='ncM',factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/ncM_f2_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['ncM']]
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'ncM',0.6,subtype=2,factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/ncM_f2_sub2_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()


# CD8 factor 7 associations
pbmc_container <- get_subclust_enr_bplot(pbmc_container,ctype='T8',factor_use=7)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD8_f7_bplot.pdf", useDingbats = FALSE,
    width = 4, height = 3)
pbmc_container$plots$subc_bplots[['T8']]
dev.off()

myplot <- get_subclust_enr_dotplot(pbmc_container,'T8',0.6,subtype=1,factor_use=7)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/vignettes/CD8_f7_sub1_dot.pdf", useDingbats = FALSE,
    width = 4.5, height = 4)
myplot
dev.off()




## getting age association and naive T4 association with healthy vs SLE donors
# first compute donor proportions
ctype <- 'T4'
res <- .6
resolution_name <- paste0('res:',as.character(res))
subclusts <- pbmc_container$subclusters[[ctype]][[resolution_name]]

# append large cell type name to subclusters
subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})

# limit cells in subclusts to those that we actually have scores for
donor_scores <- pbmc_container$tucker_results[[1]]
donor_vec <- pbmc_container$scMinimal_full$metadata[names(subclusts),'donors']
subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]

# make subtype association plot
subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
scMinimal <- pbmc_container$scMinimal_ctype[[ctype]]
sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

# get donor proportions of subclusters
donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)

head(donor_props)

# now get donor status meta data
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Status','Age')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL
head(meta)

dsc <- pbmc_container$tucker_results[[1]][,4]
tmp <- as.data.frame(cbind(donor_props[names(dsc),1],dsc,as.character(meta[names(dsc),'Status']),meta[names(dsc),'Age']))
colnames(tmp) <- c('prop','dsc','Status','Age')
head(tmp)
tmp$prop <- as.numeric(tmp$prop)
tmp$Age <- as.numeric(tmp$Age)
tmp$dsc <- as.numeric(tmp$dsc)


# test associations for healthy donors only
library(cowplot)
test <- tmp[tmp$Status=='Healthy',]
h_prop <- ggplot(test,aes(x=Age,y=prop)) +
  geom_point(color='#F8766D') +
  xlab('Age') +
  ylab('T4_1 Proportion') +
  ggtitle('Healthy Donors') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

h_fact <- ggplot(test,aes(x=dsc,y=Age)) +
  geom_point(color='#F8766D') +
  xlab('Factor 4 Donor Score') +
  ylab('Age') +
  ggtitle('Healthy Donors') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

lm_h_prop <- lm(Age~prop,data=test)
summary(lm_h_prop)

lm_h_fact <- lm(Age~dsc,data=test)
summary(lm_h_fact)


# test associations for SLE donors only
test <- tmp[tmp$Status=='Managed',]
s_prop <- ggplot(test,aes(x=Age,y=prop)) +
  geom_point(color='#00BFC4') +
  xlab('Age') +
  ylab('T4_1 Proportion') +
  ggtitle('Managed Donors') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

s_fact <- ggplot(test,aes(x=dsc,y=Age)) +
  geom_point(color='#00BFC4') +
  xlab('Factor 4 Donor Score') +
  ylab('Age') +
  ggtitle('Managed Donors') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

lm_s_prop <- lm(Age~prop,data=test)
summary(lm_s_prop)

lm_s_fact <- lm(Age~dsc,data=test)
summary(lm_s_fact)

age_prop_plots <- plot_grid(h_fact,h_prop,s_fact,s_prop,nrow=2,scale=.95)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_f4_age_prop_plots.pdf", useDingbats = FALSE,
    width = 7.5, height = 4.75)
age_prop_plots
dev.off()



# now do whole ctype proportion association tests
pbmc_container <- get_ctype_prop_associations(pbmc_container,'adj_pval',n_col=3)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_major_ctype_associations.pdf", useDingbats = FALSE,
    width = 17.5, height = 11.5)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()







## testing whether TLR7 is up in f4 low score donors
# parse data by cell type
pbmc_container <- parse_data_by_ctypes(pbmc_container)

# clean counts
pbmc_container <- clean_data(pbmc_container, donor_min_cells=5, gene_min_cells=1)

# collapse data to donor-level
pbmc_container <- get_pseudobulk(pbmc_container)

# normalize data
pbmc_container <- normalize_pseudobulk(pbmc_container, method='trim', scale_factor=10000)

# apply batch correction
pbmc_container <- apply_combat(pbmc_container,batch_var='pool')

# look at expresion of TLR7 
'TLR7' %in% colnames(pbmc_container[["scMinimal_ctype"]][["cM"]][["pseudobulk"]])

gene_test <- 'TLR7'
gene_test <- 'TLR8'
ctype <- 'cM'
ctype <- 'ncM'
d_exp <- pbmc_container[["scMinimal_ctype"]][[ctype]][["pseudobulk"]][,gene_test]
dsc <- pbmc_container$tucker_results[[1]][,4]
tmp <- cbind.data.frame(dsc,d_exp[names(dsc)])
colnames(tmp) <- c('dscore','expres')
plot(tmp$dscore,tmp$expres)
summary(lm(expres~dscore,data=tmp))








## testing whether any factors are associated with IFN after limiting them to just sle patients
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Status','Age')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL
head(meta)

d_keep <- rownames(meta)[meta$Status=='Managed']
d_keep <- d_keep[d_keep %in% rownames(pbmc_container$tucker_results[[1]])]
# d_keep <- rownames(meta)

d_exp <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,'MX1']
d_exp <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,'ISG20']
dsc <- pbmc_container$tucker_results[[1]][,2]
tmp <- cbind.data.frame(dsc[d_keep],d_exp[d_keep])
colnames(tmp) <- c('dscore','expres')
plot(tmp$dscore,tmp$expres)
summary(lm(expres~dscore,data=tmp))


pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=3,
                                       top_n_per_ctype=6, show_donor_labels=F)


pbmc_container[["plots"]][["donor_sig_genes"]][[2]]



# try plotting IFI44L levels vs sledai score
trim_names <- sapply(names(d_exp), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(d_exp) <- trim_names

dsc <- dsc[d_keep]
trim_names <- sapply(names(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(dsc) <- trim_names

colnames(clin_vars)

# d_both <- rownames(clin_vars)[rownames(clin_vars) %in% names(d_exp)]
d_both <- trim_names

tmp <- cbind.data.frame(clin_vars[d_both,'sledaiscore'],d_exp[d_both])
tmp <- cbind.data.frame(clin_vars[d_both,'sledaiscore'],dsc[d_both])
colnames(tmp) <- c('dscore','expres')
plot(tmp$expres,tmp$dscore)
summary(lm(expres~dscore,data=tmp))





# compare similar factors to find out what's similar/different
pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,4), direction=c('up','down'),
                                  compare_type='different', sig_thresh=0.05)
pbmc_container$plots$comparisons[['2_4']]
diff_enr <- get_compare_go_enrich(pbmc_container,'cDC',1)
diff_enr <- get_compare_go_enrich(pbmc_container,'ncM',1)
diff_enr <- get_compare_go_enrich(pbmc_container,'NK',1)
diff_enr <- get_compare_go_enrich(pbmc_container,'T8',-1)
diff_enr <- get_compare_go_enrich(pbmc_container,'T4',-1)
diff_enr <- get_compare_go_enrich(pbmc_container,'B',-1)
print(diff_enr[order(diff_enr,decreasing=F)][1:20])





# computing average number of cells per donor per cell type.
for (ct in pbmc_container$experiment_params$ctypes_use) {
  print(ct)
  print(median(table(pbmc_container$scMinimal_ctype[[ct]]$metadata$donors)))
}


























