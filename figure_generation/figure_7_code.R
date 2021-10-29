library(scITD)
library(Seurat)
library(RColorBrewer)
library(ggplot2)

pbmc <- readRDS(file="/home/jmitchel/data/covid_data_uk/haniffa21_subset.rds")

# remove the patient with malignant cells
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$sample_id!='BGCV10_CV0198']
pbmc <- subset(pbmc,cells=cells_keep)

# just keep in mind there is a patient with covid but negative test...

## trying to remove the lps patients
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$Status!='LPS']
pbmc <- subset(pbmc,cells=cells_keep)
##


## saving the dataset without lps or cancer patients as this is used in the analysis below.
# saveRDS(pbmc,file="/home/jmitchel/data/covid_data_uk/haniffa21_subset_no_lps.rds")
pbmc <- readRDS(file="/home/jmitchel/data/covid_data_uk/haniffa21_subset_no_lps.rds")



# collapsing cell subtypes and converting cell type names 
new_names <- sapply(as.character(pbmc@meta.data$full_clustering), function(x){
  if (x=='B_exhausted' || x=='B_immature' || x=='B_naive' || x=='B_non-switched_memory' || x=='B_switched_memory') {
    return('B')
  } else if (x=='CD14_mono' || x=='CD83_CD14_mono') {
    return('cMono')
  } else if (x=='C1_CD16_mono' || x=='CD16_mono') {
    return('ncMono')
  } else if (x=='CD4.CM' || x=='CD4.EM' || x=='CD4.IL22' || x=='CD4.Naive' || x=='CD4.Prolif' || x=='CD4.Tfh' || x=='CD4.Th1') {
    return('Th')
  } else if (x=='CD8.EM' || x=='CD8.Naive' || x=='CD8.Prolif' || x=='CD8.TE') {
    return('Tc')
  } else if (x=='NK_16hi' || x=='NK_56hi' || x=='NK_prolif') {
    return('NK')
  } else {
    return(x)
  }
})

names(new_names) <- NULL
pbmc@meta.data$initial_clustering <- factor(new_names,levels=unique(new_names))

param_list <- initialize_params(ctypes_use = c("B","Tc","Th","NK","cMono","Platelets"), ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(count_data = pbmc@assays$raw@counts,
                                     meta_data = pbmc@meta.data,
                                     params=param_list,
                                     metadata_cols=c('patient_id',
                                                     "Sex",
                                                     "Age_interval",
                                                     "initial_clustering",
                                                     "Status",
                                                     "Status_on_day_collection",
                                                     "Status_on_day_collection_summary",
                                                     'Days_from_onset',
                                                     "Site",
                                                     "Smoker",
                                                     'Worst_Clinical_Status',
                                                     'Outcome',
                                                     'Swab_result',
                                                     'time_after_LPS'),
                                     metadata_col_nm=c('donors',
                                                       'sex',
                                                       'age',
                                                       'ctypes',
                                                       'status',
                                                       'status_on_day_collection',
                                                       'status_on_day_collection_summary',
                                                       'days_from_onset',
                                                       'site',
                                                       "smoker",
                                                       'worst_Clinical_Status',
                                                       'outcome',
                                                       'swab_result',
                                                       'time_after_LPS'))

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=2, 
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.00001,
                              batch_var = 'site',
                              scale_var = TRUE, var_scale_power = .5) 

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(9,26), # current used
                                 tucker_type = 'regular', rotation_type = 'hybrid') # good with >100 num 5

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','age','status',
                                                                   'status_on_day_collection_summary',
                                                                   'site',"smoker"),stat_use='pval')

# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval',h_w=c(10,8))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_dscores.pdf", useDingbats = FALSE,
#     width = 8, height = 8)
pbmc_container$plots$donor_matrix
# dev.off()












##### plotting some metadata associations
## factor 1 associations first
f_test <- get_one_factor(pbmc_container,1)
dsc <- f_test[[1]]
pbmc_container <- get_donor_meta(pbmc_container,additional_meta = c('status_on_day_collection',
                                                                    'status_on_day_collection_summary',
                                                                    'worst_Clinical_Status',
                                                                    'swab_result',
                                                                    'status',
                                                                    'outcome',
                                                                    'site','age'),only_analyzed = FALSE)

tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'age'],
                        pbmc_container$donor_metadata[rownames(dsc),'status'])
colnames(tmp) <- c('dscore','age','status')

mycol = brewer.pal(n = 8, name = "Dark2")

tmp2 <- tmp
tmp2$status <- rep('violin',nrow(tmp2))
p <- ggplot(tmp,aes(x=age,y=dscore,fill=status)) +
  geom_violin(data=tmp2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.85, binwidth = .01) +
  ylab('Factor 1 Donor Score') +
  xlab('Age range') +
  coord_flip() +
  scale_fill_manual(values=c(mycol[5], mycol[2], 'light gray')) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=26))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f1_dsc_age.pdf", useDingbats = FALSE,
#     width = 12.5, height = 6.5)
p
# dev.off()


tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'status_on_day_collection_summary'],
                        pbmc_container$donor_metadata[rownames(dsc),'status'])
colnames(tmp) <- c('dscore','status_on_day_collection_summary','status')
# order severity levels appropriately
tmp$status_on_day_collection_summary <- factor(tmp$status_on_day_collection_summary,levels=c('Healthy','Asymptomatic','Mild','Moderate','Severe','Critical'))


tmp2 <- tmp
tmp2$status <- rep('violin',nrow(tmp2))
p <- ggplot(tmp,aes(x=status_on_day_collection_summary,y=dscore,fill=status)) +
  geom_violin(data=tmp2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.85, binwidth = .01) +
  ylab('Factor 3 Donor Score') +
  xlab('Severity on collection day') +
  coord_flip() +
  scale_fill_manual(values=c(mycol[5], mycol[2], 'light gray')) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=26))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f1_dsc_severity.pdf", useDingbats = FALSE,
#     width = 12.5, height = 6.5)
p
# dev.off()



## now plotting factor 3 associations
f_test <- get_one_factor(pbmc_container,3)
dsc <- f_test[[1]]
pbmc_container <- get_donor_meta(pbmc_container,additional_meta = c('status_on_day_collection',
                                                                    'status_on_day_collection_summary',
                                                                    'worst_Clinical_Status',
                                                                    'swab_result',
                                                                    'status',
                                                                    'outcome',
                                                                    'site','age'),only_analyzed = F)

tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'status_on_day_collection_summary'],
                        pbmc_container$donor_metadata[rownames(dsc),'status'])
colnames(tmp) <- c('dscore','status_on_day_collection_summary','status')
# order severity levels appropriately
tmp$status_on_day_collection_summary <- factor(tmp$status_on_day_collection_summary,levels=c('Healthy','Asymptomatic','Mild','Moderate','Severe','Critical'))


tmp2 <- tmp
tmp2$status <- rep('violin',nrow(tmp2))
p <- ggplot(tmp,aes(x=status_on_day_collection_summary,y=dscore,fill=status)) +
  geom_violin(data=tmp2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.85, binwidth = .01) +
  ylab('Factor 3 Donor Score') +
  xlab('Severity on collection day') +
  coord_flip() +
  scale_fill_manual(values=c(mycol[5], mycol[2], 'light gray')) +
  theme_bw() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=26))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f3_dsc_severity.pdf", useDingbats = FALSE,
#     width = 12.5, height = 8.5)
p
# dev.off()
















##### plotting IFI6 to double check that critical donors really don't have high IFN expression
# not shown in paper
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','status_on_day_collection_summary')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

d_exp <- pbmc_container[["scMinimal_ctype"]][['CD4']][["pseudobulk"]][,'IFI6']
dsc <- pbmc_container$tucker_results[[1]][,1]
tmp <- cbind.data.frame(dsc[names(d_exp)],d_exp,meta[names(d_exp),1])
colnames(tmp) <- c('dscore','expres','status')

# add regression line
lmres <- lm(expres~dscore,data=tmp)
line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

# pdf(file = "/home/jmitchel/figures/for_paper_v2/f1_IFN_status.pdf", useDingbats = FALSE,
#     width = 5, height = 2.75)
# mycol <- RColorBrewer::brewer.pal(n = 3, name = "Accent")
mycol <- RColorBrewer::brewer.pal(n = 7, name = "Dark2")
ggplot(tmp,aes(x=dscore,y=expres,color=status)) +
  geom_point(alpha = 0.75,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy,color='line')) +
  scale_color_manual(values=mycol) +
  ylab('IFI6 expression (CD4+ T)') +
  xlab('Factor 1 donor scores') +
  theme_bw()
# dev.off()





##### running gsea
# run gsea for f1
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=1, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))

plot_gsea_hmap_w_similarity(pbmc_container,factor_select=1,direc='up',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=1,direc='down',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))

plot_gsea_sub(pbmc_container,thresh=.001,clust_select=1)

## f1 sets to show on loading hmap
gsets <- c("GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
           "GOBP_RESPONSE_TO_INTERFERON_GAMMA",
           "GOBP_MYELOID_LEUKOCYTE_ACTIVATION",
           "GOBP_PLATELET_DEGRANULATION")

gset_cmap <- c('forest green',
               'black',
               'orange',
               'brown')

names(gset_cmap) <- gsets

gset_cmap_sub <- gset_cmap[gset_cmap!='black']
gset_sub <- names(gset_cmap_sub)

dev.off()
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gset_sub, le_set_colormap=gset_cmap_sub, le_set_num_per=7, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                      clust_method='mcquitty', h_w=c(9,6.5))
dev.off()
hm_list <- plot_select_sets(pbmc_container, 1, gsets, color_sets=gset_cmap, 
                            cl_rows=F, myfontsize=6.5, h_w=c(3,6.5))

p1 <- pbmc_container$plots$all_lds_plots[['1']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["1"]]

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f1_lds_go.pdf", useDingbats = FALSE,
    width = 12, height = 10)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()






# run gsea for f3
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=3, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))

plot_gsea_hmap_w_similarity(pbmc_container,factor_select=3,direc='up',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=3,direc='down',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))

plot_gsea_sub(pbmc_container,thresh=.05,clust_select=25)

## f3 sets to show on loading hmap
gsets <- c("GOBP_CELL_CYCLE",
           "GOBP_PLATELET_ACTIVATION",
           "GOBP_COAGULATION",
           "GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
           "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION",
           "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
           "GOBP_APOPTOTIC_PROCESS",
           "GOBP_REGULATION_OF_GROWTH",
           "GOBP_REGULATION_OF_STEROID_METABOLIC_PROCESS",
           "GOBP_LEUKOCYTE_DIFFERENTIATION",
           "GOBP_EPITHELIAL_CELL_APOPTOTIC_PROCESS",
           "GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_GROWTH_FACTOR_STIMULUS",
           "GOBP_CELL_ACTIVATION")

gsets <- c("GOBP_CELL_CYCLE",
           "GOBP_PLATELET_ACTIVATION",
           "GOBP_COAGULATION",
           "GOBP_LEUKOCYTE_DIFFERENTIATION",
           "GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_GROWTH_FACTOR_STIMULUS",
           "GOBP_REGULATION_OF_STEROID_METABOLIC_PROCESS",
           "GOBP_EPITHELIAL_CELL_APOPTOTIC_PROCESS",
           "GOBP_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
           "GOBP_CELL_ACTIVATION",
           "GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION",
           "GOBP_REGULATION_OF_GROWTH",
           "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
           "GOBP_APOPTOTIC_PROCESS")

gset_cmap <- c('turquoise',
               'forest green',
               'black',
               'brown',
               'black',
               'black',
               'black',
               'black',
               'black',
               'black',
               'black',
               'orange',
               'black')

names(gset_cmap) <- gsets

gset_cmap_sub <- gset_cmap[gset_cmap!='black']
gset_sub <- names(gset_cmap_sub)

dev.off()
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=3, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gset_sub, le_set_colormap=gset_cmap_sub, le_set_num_per=7, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                      clust_method='mcquitty', h_w=c(9,6.5))
dev.off()
hm_list <- plot_select_sets(pbmc_container, 3, gsets, color_sets=gset_cmap, 
                            cl_rows=FALSE, myfontsize=6.25, h_w=c(6,6.5))

p1 <- pbmc_container$plots$all_lds_plots[['3']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["3"]]

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f3_lds_go.pdf", useDingbats = FALSE,
    width = 12, height = 10)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()










##### get umap plots for figure
# keep just donors and cell types analyzed
d_kept <- rownames(pbmc_container$scMinimal_ctype[[1]]$pseudobulk)
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$patient_id %in% d_kept]
pbmc <- subset(pbmc,cells=cells_keep)

ct_kept <- pbmc_container$experiment_params$ctypes_use
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$initial_clustering %in% ct_kept]
pbmc <- subset(pbmc,cells=cells_keep)

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_umap_major_clustering.pdf", useDingbats = FALSE,
#     width = 6, height = 4)
DimPlot(pbmc, reduction = "umap", group.by = 'initial_clustering', label = FALSE) + ggtitle('Cell types analyzed')
# dev.off()

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_umap_subclusters.pdf", useDingbats = FALSE,
#     width = 8, height = 4)
DimPlot(pbmc, reduction = "umap", group.by = 'full_clustering', label = FALSE) + ggtitle('Cell subtypes')
# dev.off()

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_umap_donors.pdf", useDingbats = FALSE,
#     width = 5, height = 4)
DimPlot(pbmc, reduction = "umap", group.by = 'patient_id') + NoLegend() + ggtitle('Colored by donor')
# dev.off()









##### get median number cells per donor and total number of cells used
d_keep <- rownames(pbmc_container$scMinimal_ctype[[1]]$pseudobulk)
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$patient_id %in% d_keep]
pbmc_sub <- subset(pbmc,cells=cells_keep)
ctypes <- c("B","Tc","Th","NK","cMono","Platelets")
cells_keep <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$initial_clustering %in% ctypes]
pbmc_sub2 <- subset(pbmc_sub,cells=cells_keep)
pbmc_sub2@meta.data$patient_id <- factor(as.character(pbmc_sub2@meta.data$patient_id),levels=unique(as.character(pbmc_sub2@meta.data$patient_id)))
for (ctype in ctypes) {
  tmp <- pbmc_sub2@meta.data[pbmc@meta.data$initial_clustering==ctype,]
  print(ctype)
  print(median(table(tmp$patient_id)))
}











