library(Seurat)

pbmc <- readRDS(file="/home/jmitchel/data/covid_data_uk/haniffa21_subset.rds")
pbmc <- readRDS(file="/home/jmitchel/data/covid_data_uk/haniffa21_subset_no_lps.rds")

# remove the patient with malignant cells
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$sample_id!='BGCV10_CV0198']
pbmc <- subset(pbmc,cells=cells_keep)

# just keep in mind there is a patient with covid but negative test...

## trying to remove the lps patients
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$Status!='LPS']
pbmc <- subset(pbmc,cells=cells_keep)
##

# saveRDS(pbmc,file="/home/jmitchel/data/covid_data_uk/haniffa21_subset_no_lps.rds")


## trying to limit dataset to just sle patients
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$Status=='Covid']
pbmc <- subset(pbmc,cells=cells_keep)
##



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


# colnames(pbmc@meta.data)
# 
# table(pbmc@meta.data$full_clustering)
# table(pbmc@meta.data$initial_clustering)
# 
# # look at number of unique donors
# length(unique(pbmc@meta.data$patient_id))
# 
# param_list <- initialize_params(ctypes_use = c("B_cell","CD4","CD8","CD14","CD16",
#                                                "DCs","MAIT","NK_16hi","NK_56hi",
#                                                "Plasmablast","Platelets","Treg","gdT",
#                                                "pDC"), ncores = 30, rand_seed = 10)
# 
# param_list <- initialize_params(ctypes_use = c("B_cell","CD4","CD8","CD14","CD16",
#                                                "DCs","NK_16hi","NK_56hi",
#                                                "Platelets","Treg"), ncores = 30, rand_seed = 10) #cur
# 
# param_list <- initialize_params(ctypes_use = c("B_cell","CD4","CD8","CD14","CD16",
#                                                "NK_16hi","NK_56hi",
#                                                "Platelets"), ncores = 30, rand_seed = 10)
# 
# param_list <- initialize_params(ctypes_use = c("B_cell","CD4","CD8","CD14",
#                                                "NK_16hi"), ncores = 30, rand_seed = 10)
# 
# param_list <- initialize_params(ctypes_use = c("B_cell","CD4","CD8","CD14",
#                                                "NK_56hi"), ncores = 30, rand_seed = 10) #4alt
# 
# param_list <- initialize_params(ctypes_use = c("B_cell","CD4","CD8","CD14",
#                                                "NK_56hi","NK_16hi"), ncores = 30, rand_seed = 10) #4alt2
# 
# param_list <- initialize_params(ctypes_use = c("B_cell","CD4","CD8"), ncores = 30, rand_seed = 10)
# 
# param_list <- initialize_params(ctypes_use = c("CD4","CD14"), ncores = 30, rand_seed = 10)
# 
# param_list <- initialize_params(ctypes_use = c("B_cell","CD4","NK_56hi"), ncores = 30, rand_seed = 10) #####
# 
# param_list <- initialize_params(ctypes_use = c("CD4","CD8","B_cell","NK_56hi","NK_16hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10)
# 
# param_list <- initialize_params(ctypes_use = c("MAIT","gdT","CD14","Treg","DCs","CD4","CD8","B_cell","NK_56hi","NK_16hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10) #starting point for set discovery
# 
# ### ones with over 100 donors
# param_list <- initialize_params(ctypes_use = c("MAIT","CD4","CD8","B_cell",
#                                                "NK_56hi","NK_16hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10)
# 
# param_list <- initialize_params(ctypes_use = c("MAIT",'Treg',"CD4","CD8","B_cell",
#                                                "NK_56hi","NK_16hi"), 
#                                 ncores = 30, rand_seed = 10)
# 
# param_list <- initialize_params(ctypes_use = c('Treg',"CD4","CD8","B_cell",
#                                                "NK_56hi","NK_16hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10) 
# 
# param_list <- initialize_params(ctypes_use = c("MAIT",'Treg',"CD8","B_cell",
#                                                "NK_56hi","NK_16hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10) 
# 
# param_list <- initialize_params(ctypes_use = c("MAIT",'Treg',"CD4","CD8","B_cell",
#                                                "NK_16hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10) #5
# 
# param_list <- initialize_params(ctypes_use = c("MAIT","Treg","CD4","CD8","B_cell",
#                                                "NK_56hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10) 
# 
# param_list <- initialize_params(ctypes_use = c("MAIT","CD14","Treg","CD4","CD8",
#                                                "NK_16hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10) 
# 
# param_list <- initialize_params(ctypes_use = c("MAIT","Treg","CD4","B_cell",
#                                                "NK_56hi", "NK_16hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10) 
# 
# param_list <- initialize_params(ctypes_use = c("CD14","Treg","CD4","CD8",
#                                                "NK_56hi","NK_16hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10) #9
# 
# param_list <- initialize_params(ctypes_use = c("CD14","Treg","CD4","CD8","B_cell",
#                                                "NK_16hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10)
# 
# param_list <- initialize_params(ctypes_use = c("MAIT","Treg","CD4","CD8",
#                                                "NK_56hi", "NK_16hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10) #11
# 
# param_list <- initialize_params(ctypes_use = c("MAIT","CD14","NK_56hi"), 
#                                 ncores = 30, rand_seed = 10) #11
# 
# ## current
# param_list <- initialize_params(ctypes_use = c("CD4","CD8","B_cell","CD14","NK_56hi","Platelets"), 
#                                 ncores = 30, rand_seed = 10) #11

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

# determine_ctype_set(pbmc_container,min_d_intersect=110,n_ctypes=7,donor_min_cells=5,max_iter=400)

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=2, 
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.00001,
                              batch_var = 'site',
                              scale_var = TRUE, var_scale_power = .5) # current
# pbmc_container <- form_tensor(pbmc_container, donor_min_cells=2, 
#                               norm_method='trim', scale_factor=10000,
#                               vargenes_method='norm_var_pvals', vargenes_thresh=.05,
#                               scale_var = TRUE, var_scale_power = .5)
# pbmc_container <- form_tensor(pbmc_container, donor_min_cells=2,
#                               norm_method='trim', scale_factor=10000,
#                               vargenes_method='norm_var', vargenes_thresh=750,
#                               batch_var = 'site',
#                               scale_var = TRUE, var_scale_power = .5)

# pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, 
#                               norm_method='trim', scale_factor=10000,
#                               vargenes_method='norm_var_pvals', vargenes_thresh=.001,
#                               batch_var = 'site',
#                               scale_var = TRUE, var_scale_power = .5)

# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,20),
#                                  tucker_type = 'regular', rotation_type = 'ica_dsc')
# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(14,30),
#                                  tucker_type = 'regular', rotation_type = 'hybrid') # works well w ct group 4 or 5
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(9,26), # current used
                                 tucker_type = 'regular', rotation_type = 'hybrid') # good with >100 num 5
# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(30,40),
#                                  tucker_type = 'regular', rotation_type = 'hybrid') # good with >100 num 5
# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(14,40),
#                                  tucker_type = 'regular', rotation_type = 'hybrid') # good with >100 num 5
# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(10,30),
#                                  tucker_type = 'regular', rotation_type = 'hybrid')  # works well with new ctype naming

# pbmc_container <- run_stability_analysis(pbmc_container, ranks=c(9,26), tucker_type='regular',
#                                          rotation_type='hybrid',  sparsity=sqrt(2),
#                                          subset_type='subset', sub_prop=.75,
#                                          n_iterations=10)
# pbmc_container$plots$stability_plot_dsc
# pbmc_container$plots$stability_plot_lds

# get factor-meta data associations
# pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','age',
#                                                                    'status','status_on_day_collection',
#                                                                    'status_on_day_collection_summary',
#                                                                    'site',"smoker",'worst_Clinical_Status',
#                                                                    'outcome','swab_result'),stat_use='pval')
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','age','status',
                                                                   'status_on_day_collection_summary',
                                                                   'site',"smoker"),stat_use='pval')
# pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','age',
#                                                                    'status_on_day_collection',
#                                                                    'status_on_day_collection_summary',
#                                                                    'site',"smoker",'worst_Clinical_Status',
#                                                                    'outcome','swab_result'),stat_use='pval')

# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval',h_w=c(10,8))
# pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars = 'sex',
#                                     cluster_by_meta = 'sex',
#                                     show_donor_ids = FALSE,
#                                     add_meta_associations='pval',h_w=c(10,10))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_dscores.pdf", useDingbats = FALSE,
#     width = 8, height = 8)
pbmc_container$plots$donor_matrix
dev.off()



f_test <- get_one_factor(pbmc_container,3)
dsc <- f_test[[1]]
pbmc_container <- get_donor_meta(pbmc_container,additional_meta = c('status_on_day_collection',
                                                                    'status_on_day_collection_summary',
                                                                    'worst_Clinical_Status',
                                                                    'swab_result',
                                                                    'status',
                                                                    'outcome',
                                                                    'site','age'),only_analyzed = F)
tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'status_on_day_collection'])
tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'status_on_day_collection_summary'])
tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'worst_Clinical_Status'])
tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'swab_result'])
tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'status'])
tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'age'])
tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'outcome'])
tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'age'],
                        pbmc_container$donor_metadata[rownames(dsc),'status_on_day_collection_summary'])
tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'age'],
                        pbmc_container$donor_metadata[rownames(dsc),'status'])
tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'status_on_day_collection_summary'],
                        pbmc_container$donor_metadata[rownames(dsc),'status'])
colnames(tmp) <- c('dscore','status_on_day_collection_summary','status')
colnames(tmp) <- c('dscore','age','status')

tmp$status_on_day_collection_summary <- factor(tmp$status_on_day_collection_summary,levels=c('Healthy','Asymptomatic','Mild','Moderate','Severe','Critical'))

library(RColorBrewer)
mycol = brewer.pal(n = 8, name = "Dark2")


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
p

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
p

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f1_dsc_age.pdf", useDingbats = FALSE,
    width = 12.5, height = 6.5)
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f1_dsc_severity.pdf", useDingbats = FALSE,
    width = 12.5, height = 6.5)
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f3_dsc_severity.pdf", useDingbats = FALSE,
    width = 12.5, height = 8.5)
p
dev.off()

ggplot(tmp,aes(x=age,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 9 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()


colnames(tmp) <- c('dscore','cvar')
tmp$cvar <- factor(tmp$cvar,levels=c('Healthy','Asymptomatic','Mild','Moderate','Severe','Critical'))

ggplot(tmp,aes(x=cvar,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 9 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()




# # checking for RPS4Y1 differences between M/F
# meta <- pbmc_container$scMinimal_full$metadata[,c('donors','sex')]
# meta <- unique(meta)
# rownames(meta) <- meta$donors
# meta$donors <- NULL
# 
# pb <- pbmc_container$scMinimal_ctype[['B']]$pseudobulk[,'XIST']
# tmp <- cbind.data.frame(meta[names(pb),1],pb)
# colnames(tmp) <- c('sex','expr')
# ggplot(tmp,aes(x=as.factor(sex),y=expr)) +
#   geom_boxplot()



## look at how the LPS profile differs from the one that distinguishes covid by its severity or healthy vs covid


## it's possible that the factor distinguishing the healthy from negative samples is actually a batch factor since there
# was different notation used between the two sites. This would be a batch effect that only took place in healthy donors for some reason...

# try a decomp with just the covid patients to see if I can more easily extract the outcome associated factor




# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_lds_f1.pdf", useDingbats = FALSE,
    width = 7, height = 7)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(12,10))
dev.off()
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_lds_f3.pdf", useDingbats = FALSE,
    width = 7, height = 7)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=3, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=3, h_w=c(12,10))
dev.off()

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=5, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=3, h_w=c(12,10))

# plotting IFI6 to see if critical donors really don't have high IFN expression (maybe due to autoantibodies?)
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
dev.off()




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
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=3, method="fgsea", thresh=0.05,
                                      db_use=c("KEGG",'BioCarta', 'Reactome'))

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

gset_cmap <- c('turquoise',
               'forest green',
               'black',
               'black',
               'black',
               'orange',
               'black',
               'black',
               'black',
               'brown',
               'black',
               'black',
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
                            cl_rows=T, myfontsize=6.25, h_w=c(6,6.5))

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



## get umap plots for figure
# keep just donors and cell types analyzed
d_kept <- rownames(pbmc_container$scMinimal_ctype[[1]]$pseudobulk)
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$patient_id %in% d_kept]
pbmc <- subset(pbmc,cells=cells_keep)

ct_kept <- pbmc_container$experiment_params$ctypes_use
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$initial_clustering %in% ct_kept]
pbmc <- subset(pbmc,cells=cells_keep)

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_umap_major_clustering.pdf", useDingbats = FALSE,
    width = 6, height = 4)
DimPlot(pbmc, reduction = "umap", group.by = 'initial_clustering', label = FALSE) + ggtitle('Cell types analyzed')
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_umap_subclusters.pdf", useDingbats = FALSE,
    width = 8, height = 4)
DimPlot(pbmc, reduction = "umap", group.by = 'full_clustering', label = FALSE) + ggtitle('Cell subtypes')
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_umap_donors.pdf", useDingbats = FALSE,
    width = 5, height = 4)
DimPlot(pbmc, reduction = "umap", group.by = 'patient_id') + NoLegend() + ggtitle('Colored by donor')
dev.off()






### get median number cells per donor and total number of cells used
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











