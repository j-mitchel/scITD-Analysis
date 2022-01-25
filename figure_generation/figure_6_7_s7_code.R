# using current dev version of scITD instead of renv version...

library(Seurat)
library(RColorBrewer)
library(ggplot2)

# # This chunk run originally but faster to use saved data below...
# # load up covid PBMC dataset: see preprocessing/covid_uk_preprocessing.ipynb for code used to generate this object
# pbmc <- readRDS(file="/home/jmitchel/data/covid_data_uk/haniffa21_subset.rds")
# 
# # remove the patient with malignant cells
# cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$sample_id!='BGCV10_CV0198']
# pbmc <- subset(pbmc,cells=cells_keep)
# 
# # just keep in mind there is a patient with covid but negative test...
# 
# ## trying to remove the lps patients
# cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$Status!='LPS']
# pbmc <- subset(pbmc,cells=cells_keep)
# ##


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


# ## testing markers of plasmablast like B cell population
# b_cells <- rownames(pbmc@meta.data)[pbmc@meta.data$initial_clustering=='B']
# pbmc <- subset(pbmc,cells=b_cells)
# cells_keep1 <- rownames(pbmc@reductions[["umap"]]@cell.embeddings)[pbmc@reductions[["umap"]]@cell.embeddings[,1]>-1]
# cells_keep2 <- rownames(pbmc@reductions[["umap"]]@cell.embeddings)[pbmc@reductions[["umap"]]@cell.embeddings[,2]>.5]
# cells_keep <- unique(c(cells_keep1,cells_keep2))
# cells_not_keep <- rownames(pbmc@meta.data)[!(rownames(pbmc@meta.data) %in% cells_keep)]
# pbmc@meta.data$new_meta <- NULL
# pbmc@meta.data[cells_keep,'new_meta'] <- 'good'
# pbmc@meta.data[cells_not_keep,'new_meta'] <- 'bad'
# Idents(pbmc) <- pbmc@meta.data$new_meta
# monocyte.de.markers <- FindMarkers(pbmc, ident.1 = "good", ident.2 = "bad")
# head(monocyte.de.markers)
# # it the plasmablast like cells underexpresses CD19 and AB-CD19 and overexpress MKI67

# remove plasmablast-like B cells
cells_keep1 <- rownames(pbmc@reductions[["umap"]]@cell.embeddings)[pbmc@reductions[["umap"]]@cell.embeddings[,1]>-1]
cells_keep2 <- rownames(pbmc@reductions[["umap"]]@cell.embeddings)[pbmc@reductions[["umap"]]@cell.embeddings[,2]>.5]
cells_keep <- unique(c(cells_keep1,cells_keep2))

pbmc <- subset(pbmc,cells=cells_keep)




# remove antibody "genes"
g_ndx_rem <- which(pbmc@assays[["RNA"]]@meta.features=='Antibody Capture')
g_ndx_keep <- c(1:nrow(pbmc@assays[["raw"]]@counts))[!(c(1:nrow(pbmc@assays[["raw"]]@counts)) %in% g_ndx_rem)]
dim(pbmc)
pbmc <- CreateSeuratObject(pbmc@assays[["raw"]]@counts[g_ndx_keep,], project = "SeuratProject", assay = "raw",
                               meta.data = pbmc@meta.data)
dim(pbmc)


# starting analysis
param_list <- initialize_params(ctypes_use = c("B","Tc","Th","NK","cMono"), ncores = 30, rand_seed = 10)

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
                              vargenes_method='norm_var_pvals', vargenes_thresh=.0001,
                              batch_var = 'site',
                              scale_var = TRUE, var_scale_power = .5) 

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(9,38), 
                                 tucker_type = 'regular', rotation_type = 'hybrid') # best with vargenes_thresh=.0001 also .01

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','age',
                                                                   'status_on_day_collection_summary',
                                                                   'site',"smoker"),stat_use='pval')
# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval',h_w=c(10,8))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_dscores2.pdf", useDingbats = FALSE,
#     width = 8, height = 8)
pbmc_container$plots$donor_matrix
# dev.off()

# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)











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

tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'status_on_day_collection_summary'],
                        pbmc_container$donor_metadata[rownames(dsc),'status'])
colnames(tmp) <- c('dscore','status_on_day_collection_summary','status')
# order severity levels appropriately
tmp$status_on_day_collection_summary <- factor(tmp$status_on_day_collection_summary,levels=c('Healthy','Asymptomatic','Mild','Moderate','Severe','Critical'))

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

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f1_dsc_severity.pdf", useDingbats = FALSE,
#     width = 12.5, height = 6.5)
p
# dev.off()


# plotting F1 against status and coloring the severe critical patients differently
tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'status_on_day_collection_summary'],
                        pbmc_container$donor_metadata[rownames(dsc),'status'])
colnames(tmp) <- c('dscore','status_on_day_collection_summary','status')
# order severity levels appropriately
tmp$status_on_day_collection_summary <- factor(tmp$status_on_day_collection_summary,levels=c('Healthy','Asymptomatic','Mild','Moderate','Severe','Critical'))
tmp$is_critical <- sapply(as.character(tmp$status_on_day_collection_summary),function(x){
  if (x=='Critical') {
    return('Critical condition')
  } else {
    return('Other')
  }
})
tmp$is_critical <- factor(tmp$is_critical,levels=c('Other','Critical condition'))

tmp2 <- tmp
tmp2$is_critical <- rep('violin',nrow(tmp2))
p <- ggplot(tmp,aes(x=status,y=dscore,fill=is_critical)) +
  geom_violin(data=tmp2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.25, binwidth = .01) +
  ylab('Factor 1 Donor Score') +
  xlab('COVID-19 Status') +
  coord_flip() +
  scale_fill_manual(values=c('red', 'black', 'light gray')) +
  theme_bw() +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f1_status.pdf", useDingbats = FALSE,
#     width = 12.5, height = 8.5)
p
# dev.off()




# plotting F1 against status of three groups, healthy, covid, covid-critical
tmp <- cbind.data.frame(dsc,pbmc_container$donor_metadata[rownames(dsc),'status_on_day_collection_summary'],
                        pbmc_container$donor_metadata[rownames(dsc),'status'])
colnames(tmp) <- c('dscore','status_on_day_collection_summary','status')
# order severity levels appropriately
tmp$status_on_day_collection_summary <- factor(tmp$status_on_day_collection_summary,levels=c('Healthy','Asymptomatic','Mild','Moderate','Severe','Critical'))
tmp$is_critical <- sapply(as.character(tmp$status_on_day_collection_summary),function(x){
  if (x=='Critical') {
    return('COVID-19-critical')
  } else if (x=='Healthy') {
    return(x)
  } else {
    return('COVID-19')
  }
})
tmp$is_critical <- factor(tmp$is_critical,levels=c('Healthy','COVID-19','COVID-19-critical'))

tmp2 <- tmp
tmp2$status <- rep('violin',nrow(tmp2))
p <- ggplot(tmp,aes(x=is_critical,y=dscore,fill=status)) +
  geom_violin(data=tmp2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.25, binwidth = .01) +
  ylab('Factor 1 Donor Score') +
  xlab('COVID-19 Status') +
  coord_flip() +
  scale_fill_manual(values=c(mycol[5], mycol[2], 'light gray')) +
  theme_bw() +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26))

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f1_status2.pdf", useDingbats = FALSE,
    width = 12.5, height = 8.5)
p
dev.off()





## getting enrichment p-values for f1 with status
enr_fig <- plot_dscore_enr(pbmc_container, factor_use=1, meta_var='status')
enr_fig


## trying enrichment without healthy donors
container=pbmc_container
meta_var='status_on_day_collection_summary'
factor_use=1
meta <- unique(container$scMinimal_full$metadata[,c('donors',meta_var)])
rownames(meta) <- meta$donors
meta$is_critical <- sapply(as.character(meta$status_on_day_collection_summary),function(x){
  if (x=='Critical') {
    return('COVID-19-critical')
  } else if (x=='Healthy') {
    return('Healthy')
  } else {return('COVID-19')}
})
meta_var <- 'is_critical'

meta_vals <- unlist(unique(as.character(meta[,meta_var])))
d_use <- intersect(rownames(container$tucker_results[[1]]),rownames(meta))
meta <- meta[d_use,]
mypaths <- list()
for (i in 1:length(meta_vals)) {
  mypaths[[meta_vals[i]]] <- rownames(meta)[meta[,meta_var]==meta_vals[i]]
}

myranks <- container$tucker_results[[1]][d_use,factor_use]

fgseaRes <- fgsea::fgsea(pathways = mypaths,
                         stats    = myranks,
                         minSize  = 0,
                         maxSize  = 5000)

plt_lst <- list()
for (i in 1:length(meta_vals)) {
  plt <- fgsea::plotEnrichment(mypaths[[meta_vals[i]]],
                               myranks) + labs(title=paste0(meta_vals[i],' - Factor ',as.character(factor_use)))
  plt <- plt +
    annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
             label=paste0('adj pval: ',
                          round(fgseaRes[fgseaRes$pathway==meta_vals[i],'padj'],digits=4)))
  
  plt_lst[[i]] <- plt
}

fig <- cowplot::plot_grid(plotlist=plt_lst,nrow=1)
fig











## now plotting factor 2 associations
f_test <- get_one_factor(pbmc_container,2)
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

mycol = brewer.pal(n = 8, name = "Dark2")

tmp2 <- tmp
tmp2$status <- rep('violin',nrow(tmp2))
p <- ggplot(tmp,aes(x=status_on_day_collection_summary,y=dscore,fill=status)) +
  geom_violin(data=tmp2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.25, binwidth = .01) +
  ylab('Factor 2 Donor Score') +
  xlab('Severity on collection day') +
  coord_flip() +
  scale_fill_manual(values=c(mycol[5], mycol[2], 'light gray')) +
  theme_bw() +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f2_severity2.pdf", useDingbats = FALSE,
#     width = 12.5, height = 8.5)
p
# dev.off()















##### plotting IFI6 to double check that critical donors really don't have high IFN expression
# not shown in paper
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','status_on_day_collection_summary')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

d_exp <- pbmc_container[["scMinimal_ctype"]][['Th']][["pseudobulk"]][,'IFI6']
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

plot_gsea_sub(pbmc_container,thresh=.05,clust_select=1)

## f1 sets to show on loading hmap
gsets <- c("GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
           "GOBP_RESPONSE_TO_INTERFERON_GAMMA")

gset_cmap <- c('forest green',
               'black')

#### testing out adding some more gene sets
gsets <- c("GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
           "GOBP_RESPONSE_TO_INTERFERON_GAMMA",
           "GOBP_PROTEOLYSIS",
           "GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
           "GOBP_SECRETION",
           "GOBP_REGULATION_OF_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION",
           "GOBP_CELL_CYCLE")

gset_cmap <- c('forest green',
               'black',
               'black',
               'black',
               'black',
               'black',
               'black')

####

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
                            cl_rows=F, myfontsize=6.5, h_w=c(4.5,6.5))

p1 <- pbmc_container$plots$all_lds_plots[['1']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["1"]]

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f1_lds_go2.pdf", useDingbats = FALSE,
    width = 12, height = 10)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()

# to save just the gene sets component
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f1_lds_go3.pdf", useDingbats = FALSE,
    width = 12, height = 10)
hm_list
dev.off()





# run gsea for f2
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f2_all_go_up.pdf", useDingbats = FALSE,
    width = 12, height = 10)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=2,direc='up',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))
dev.off()
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f2_all_go_down.pdf", useDingbats = FALSE,
    width = 12, height = 10)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=2,direc='down',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))
dev.off()

plot_gsea_sub(pbmc_container,thresh=.05,clust_select=20)

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
               'black',
               'black',
               'black',
               'black')

gsets <- c("GOBP_CELL_CYCLE",
           "GOBP_POSITIVE_REGULATION_OF_MEMBRANE_PERMEABILITY",
           "GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
           "GOBP_LYMPHOCYTE_MEDIATED_IMMUNITY",
           "GOBP_CELL_ACTIVATION",
           "GOBP_APOPTOTIC_PROCESS",
           "GOBP_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR",
           "GOBP_CELL_KILLING")


gset_cmap <- c('blue',
               'green',
               'orange',
               'purple',
               'black',
               'black')
gset_cmap <- c('black',
               'black',
               'orange',
               'black',
               'black',
               'black')

gsets <- c("GOBP_CELL_CYCLE",
           "GOBP_CELL_KILLING",
           "GOBP_MYELOID_LEUKOCYTE_MEDIATED_IMMUNITY",
           "GOBP_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR",
           'GOBP_CELLULAR_CARBOHYDRATE_METABOLIC_PROCESS',
           'GOBP_REGULATION_OF_CHEMOTAXIS',
           'GOBP_JNK_CASCADE',
           'GOBP_PROTEIN_KINASE_B_SIGNALING')
gset_cmap <- c('black',
               'black',
               'red',
               'orange',
               'green',
               'black',
               'black',
               'black')
names(gset_cmap) <- gsets

gset_cmap_sub <- gset_cmap[gset_cmap!='black']
gset_sub <- names(gset_cmap_sub)

dev.off()
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gset_sub, le_set_colormap=gset_cmap_sub, le_set_num_per=5, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                      clust_method='mcquitty', h_w=c(9,6.5))
dev.off()
hm_list <- plot_select_sets(pbmc_container, 2, gsets, color_sets=gset_cmap, 
                            cl_rows=FALSE, myfontsize=6.25, h_w=c(6,6.5))

p1 <- pbmc_container$plots$all_lds_plots[['2']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["2"]]

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f2_lds_go2.pdf", useDingbats = FALSE,
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

# since I created a new seurat object at the top, I lost their original umap...
# will reload pbmc to get their umap back
pbmc_tmp <- readRDS(file="/home/jmitchel/data/covid_data_uk/haniffa21_subset_no_lps.rds")
cells_keep <- rownames(pbmc@meta.data)
pbmc_tmp <- subset(pbmc_tmp,cells=cells_keep)

new_names <- sapply(as.character(pbmc_tmp@meta.data$full_clustering), function(x){
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
pbmc_tmp@meta.data$initial_clustering <- factor(new_names,levels=unique(new_names))

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_umap_major_clustering2.pdf", useDingbats = FALSE,
    width = 6, height = 4)
DimPlot(pbmc_tmp, reduction = "umap", group.by = 'initial_clustering', label = FALSE) + ggtitle('Cell types analyzed')
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_umap_subclusters2.pdf", useDingbats = FALSE,
    width = 8, height = 4)
DimPlot(pbmc_tmp, reduction = "umap", group.by = 'full_clustering', label = FALSE) + ggtitle('Cell subtypes')
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_umap_donors2.pdf", useDingbats = FALSE,
    width = 5, height = 4)
DimPlot(pbmc_tmp, reduction = "umap", group.by = 'patient_id') + NoLegend() + ggtitle('Colored by donor')
dev.off()

rm(pbmc_tmp)
gc()







##### get median number cells per donor and total number of cells used
d_keep <- rownames(pbmc_container$scMinimal_ctype[[1]]$pseudobulk)
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$patient_id %in% d_keep]
pbmc_sub <- subset(pbmc,cells=cells_keep)
ctypes <- c("B","Tc","Th","NK","cMono")
cells_keep <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$initial_clustering %in% ctypes]
pbmc_sub2 <- subset(pbmc_sub,cells=cells_keep)
pbmc_sub2@meta.data$patient_id <- factor(as.character(pbmc_sub2@meta.data$patient_id),levels=unique(as.character(pbmc_sub2@meta.data$patient_id)))
for (ctype in ctypes) {
  tmp <- pbmc_sub2@meta.data[pbmc@meta.data$initial_clustering==ctype,]
  print(ctype)
  print(median(table(tmp$patient_id)))
}




















##### prepping the SLE data for comparing IFN factors
# will compare to full SLE decomp as in fig 2
# load up the subsetted dataset: see preprocessing/lupus_preprocessing.R for code used to generate this object
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')

# converting shorthand cell type
new_names <- sapply(as.character(pbmc@meta.data$cg_cov), function(x){
  if (x=='cM') {
    return('cMono')
  } else if (x=='ncM') {
    return('ncMono')
  } else if (x=='T4') {
    return('Th')
  } else if (x=='T8') {
    return('Tc')
  } else {
    return(x)
  }
})
names(new_names) <- NULL
pbmc@meta.data$cg_cov <- factor(new_names,levels=unique(new_names))

# set up project parameters
param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc","cDC",
                                               "cMono","ncMono"),
                                ncores = 30, rand_seed = 10)

pbmc_container_SLE <- make_new_container(seurat_obj=pbmc,
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


pbmc_container_SLE <- form_tensor(pbmc_container_SLE, donor_min_cells=20, 
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')


pbmc_container_SLE <- run_tucker_ica(pbmc_container_SLE, ranks=c(7,20),
                                 tucker_type = 'regular', rotation_type = 'hybrid')

# get factor-meta data associations
pbmc_container_SLE <- get_meta_associations(pbmc_container_SLE,vars_test=c('sex','Age','pool','processing','Status','Ethnicity'),
                                        stat_use='pval')

# plot donor scores by status
pbmc_container_SLE <- plot_donor_matrix(pbmc_container_SLE, meta_vars=c('sex'),
                                    cluster_by_meta = 'sex',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# pdf(file = "/home/jmitchel/figures/for_paper_v2/lupus_dscores_full2.pdf", useDingbats = FALSE,
#     width = 6, height = 7)
pbmc_container_SLE$plots$donor_matrix
# dev.off()












############################### setting up the analysis in one pipeline
get_tstats <- function(container, n.cores = container$experiment_params$ncores) {
  tensor_data <- container$tensor_data
  tucker_results <- container$tucker_results
  
  if (is.null(tucker_results)) {
    stop('Need to run run_tucker_ica() first')
  }
  
  n_genes <- length(tensor_data[[2]])
  n_ctypes <- length(tensor_data[[3]])
  all_pvals <- data.frame(matrix(ncol=3,nrow=0))
  if (is.null(container$experiment_params$ncores)){
    n.cores = 1
  }
  all_pvals <- mclapply(1:n_genes, function(i) {
    gene <- tensor_data[[2]][i]
    gene_res <- list()
    for (j in 1:n_ctypes) {
      ctype <- tensor_data[[3]][j]
      gene_res[[ctype]] <- list()
      for (k in 1:ncol(tucker_results[[1]])) {
        tmp_fiber <- tensor_data[[4]][,i,j]
        
        # if expression is 0 for all donors just skip
        if (sum(tmp_fiber==0)==length(tmp_fiber)) {
          gene_res[[ctype]][[as.character(k)]] <- NA
        } else {
          df_test <- as.data.frame(cbind(tmp_fiber, tucker_results[[1]][,k]))
          colnames(df_test) <- c('fiber','factor')
          lmres <- lm(factor~fiber,df_test)
          
          x <- summary(lmres)
          tstat <- x$coefficients['fiber','t value']
          
          gene_res[[ctype]][[as.character(k)]] <- tstat
        }
      }
    }
    return(gene_res)
  }, mc.cores = n.cores)
  
  names(all_pvals) <- tensor_data[[2]]
  
  # unpack the list
  all_pvals <- unlist(all_pvals)
  
  # set NA values to 0
  all_pvals[is.na(all_pvals)] <- 0
  
  new_names <- sapply(names(all_pvals),function(x) {
    tmp <- strsplit(x,split = '.', fixed = TRUE)[[1]]
    if (length(tmp)==3) {
      return(paste0(tmp[[1]],'.',tmp[[2]],'.',tmp[[3]]))
    } else if (length(tmp)==4) {
      return(paste0(tmp[[1]],'.',tmp[[2]],'.',tmp[[3]],'.',tmp[[4]]))
    } else if (length(tmp)==5) {
      return(paste0(tmp[[1]],'.',tmp[[2]],'.',tmp[[3]],'.',tmp[[4]],'.',tmp[[5]]))
    } else if (length(tmp)==6) {
      return(paste0(tmp[[1]],'.',tmp[[2]],'.',tmp[[3]],'.',tmp[[4]],'.',tmp[[5]],'.',tmp[[6]]))
    }
  })
  names(new_names) <- NULL
  names(all_pvals) <- new_names
  container[["gene_score_associations"]] <- all_pvals
  return(container)
}

## reprocess both to include all genes
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=2, 
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=1,
                              batch_var = 'site',
                              scale_var = TRUE, var_scale_power = .5) 

pbmc_container_SLE <- form_tensor(pbmc_container_SLE, donor_min_cells=20, 
                                  norm_method='trim', scale_factor=10000,
                                  vargenes_method='norm_var_pvals', vargenes_thresh=1,
                                  scale_var = TRUE, var_scale_power = .5,
                                  batch_var='pool')

## getting gene significance pvals
pbmc_container <- get_lm_pvals(pbmc_container)
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=1, pbmc_container$experiment_params$ctypes_use)
sig_df_covid2 <- t(as.data.frame(do.call(rbind, sig_vectors)))

# for sle data
pbmc_container_SLE <- get_lm_pvals(pbmc_container_SLE)
sig_vectors <- get_significance_vectors(pbmc_container_SLE,
                                        factor_select=1, pbmc_container_SLE$experiment_params$ctypes_use)
sig_df_SLE2 <- t(as.data.frame(do.call(rbind, sig_vectors)))

## getting tstats
# for covid data
pbmc_container <- get_tstats(pbmc_container)
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=1, pbmc_container$experiment_params$ctypes_use)
sig_df_covid <- t(as.data.frame(do.call(rbind, sig_vectors)))

# for sle data
pbmc_container_SLE <- get_tstats(pbmc_container_SLE)
sig_vectors <- get_significance_vectors(pbmc_container_SLE,
                                        factor_select=1, pbmc_container_SLE$experiment_params$ctypes_use)
sig_df_SLE <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df_SLE <- sig_df_SLE * -1 # reverse it so increased ISG is positive
genes_both <- intersect(rownames(sig_df_covid),rownames(sig_df_SLE))

tmp1 <- cbind.data.frame(sig_df_covid[,1],sig_df_covid2[rownames(sig_df_covid),1])
colnames(tmp1) <- c('tstat','pval')
tmp1 <- tmp1[order(abs(tmp1$tstat),decreasing=TRUE),]
ndx_thr <- which(tmp1$pval<.01)
tmp1[508:515,]
covid_thresh <- abs(tmp1[max(ndx_thr),'tstat'])

tmp1 <- cbind.data.frame(sig_df_SLE[,1],sig_df_SLE2[rownames(sig_df_SLE),1])
colnames(tmp1) <- c('tstat','pval')
tmp1 <- tmp1[order(abs(tmp1$tstat),decreasing=TRUE),]
ndx_thr <- which(tmp1$pval<.01)
tmp1[1:100,]
sle_thresh <- abs(tmp1[max(ndx_thr),'tstat'])





## trying to use clusterProfiler for gene overrepresentation analysis
# will use union of IFN upregulated genes as background
library(clusterProfiler)
library(org.Hs.eg.db)

## NK cells
covid_sub <- sig_df_covid[genes_both,'NK']
covid_sig_up <- names(covid_sub)[covid_sub>covid_thresh] 
covid_sig_down <- names(covid_sub)[covid_sub<(-covid_thresh)] 
covid_sig_all <- names(covid_sub)[abs(covid_sub)>covid_thresh] 

sle_sub <- sig_df_SLE[genes_both,'NK']
sle_sig_up <- names(sle_sub)[sle_sub>sle_thresh]
sle_sig_down <- names(sle_sub)[sle_sub<(-sle_thresh)]
sle_sig_all <- names(sle_sub)[abs(sle_sub)>sle_thresh] 

bg <- unique(c(covid_sig_up,sle_sig_up))
covid_unique_up <- covid_sig_up[!(covid_sig_up %in% sle_sig_up)]
sle_unique_up <- sle_sig_up[!(sle_sig_up %in% covid_sig_up)]

## writing the genes to csvs
sig_both <- intersect(covid_sig_up,sle_sig_up)
write.csv(sig_both,file='/home/jmitchel/data/gene_lists/NK_sig_both.csv')
write.csv(covid_unique_up,file='/home/jmitchel/data/gene_lists/NK_sig_covid.csv')
write.csv(sle_unique_up,file='/home/jmitchel/data/gene_lists/NK_sig_sle.csv')

# bg <- unique(c(covid_sig_down,sle_sig_down))
# covid_unique_up <- covid_sig_down[!(covid_sig_down %in% sle_sig_down)]
# sle_unique_up <- sle_sig_down[!(sle_sig_down %in% covid_sig_down)]

bg <- unique(c(covid_sig_all,sle_sig_all))
covid_unique_up <- covid_sig_all[!(covid_sig_all %in% sle_sig_all)]
sle_unique_up <- sle_sig_all[!(sle_sig_all %in% covid_sig_all)]

enr_res_c1 <- enrichGO(covid_unique_up,
                    keyType = 'SYMBOL',
                    universe=bg,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    minGSSize = 10)

enr_res_c1@result[1:10,c('Description','p.adjust')]

enr_res_s1 <- enrichGO(sle_unique_up, #increased maxGSSize for proteolysis enrichment calculation...
                    keyType = 'SYMBOL',
                    universe=bg,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    minGSSize = 10)

enr_res_s1@result[1:10,c('Description','p.adjust')]
# sets to show later:
# leukocyte cell-cell adhesion, tissue development, apoptotic process, cell activation, leukocyte migration,
# carboxylic acid metabolic process, Fc receptor signaling pathway
pth <- c('leukocyte cell-cell adhesion', 'tissue development', 'apoptotic process', 
         'cell activation', 'leukocyte migration',
         'carboxylic acid metabolic process', 'Fc receptor signaling pathway')
enr_res_s1@result[enr_res_s1@result$Description %in% pth,]

# proteolysis p-value: .016

# save relevant genes
tmp <- as.data.frame(enr_res_s1@result[enr_res_s1@result$Description %in% pth,'geneID'])
write.csv(tmp,file='/home/jmitchel/data/gene_lists/NK_gsets.csv')


## other genes to potentially show
# TNFAIP8L2, SELL, CX3CR1, CD74, HAVCR2, HLA-DRB1, HLA-DPA1, HLA-DPB1,
# CRTAM, IL2RA, ITGB7, LGALS3, PYCARD, LGALS9B, ADA
# ENO1, PSMB2, PSMB2, PSMA5, ELOVL6, HPGD, CD74, MTHFD1L, FABP5, AZIN1, PFKP, 
# ACSL5, GAPDH, PSME1, PSMA6, PSMA4, PSMB10

## Th cells
covid_sub <- sig_df_covid[genes_both,'Th']
covid_sig_up <- names(covid_sub)[covid_sub>covid_thresh] 
covid_sig_down <- names(covid_sub)[covid_sub<(-covid_thresh)] 

sle_sub <- sig_df_SLE[genes_both,'Th']
sle_sig_up <- names(sle_sub)[sle_sub>sle_thresh]
sle_sig_down <- names(sle_sub)[sle_sub<(-sle_thresh)]

bg <- unique(c(covid_sig_up,sle_sig_up))
covid_unique_up <- covid_sig_up[!(covid_sig_up %in% sle_sig_up)]
sle_unique_up <- sle_sig_up[!(sle_sig_up %in% covid_sig_up)]
# bg <- unique(c(covid_sig_down,sle_sig_down))
# covid_unique_up <- covid_sig_down[!(covid_sig_down %in% sle_sig_down)]
# sle_unique_up <- sle_sig_down[!(sle_sig_down %in% covid_sig_down)]

## writing the genes to csvs
sig_both <- intersect(covid_sig_up,sle_sig_up)
write.csv(sig_both,file='/home/jmitchel/data/gene_lists/Th_sig_both.csv')
write.csv(covid_unique_up,file='/home/jmitchel/data/gene_lists/Th_sig_covid.csv')
write.csv(sle_unique_up,file='/home/jmitchel/data/gene_lists/Th_sig_sle.csv')

enr_res_c2 <- enrichGO(covid_unique_up,
                    keyType = 'SYMBOL',
                    universe=bg,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    minGSSize = 10)

enr_res_c2@result[1:10,c('Description','p.adjust')]

enr_res_s2 <- enrichGO(sle_unique_up,
                    keyType = 'SYMBOL',
                    universe=bg,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    minGSSize = 10)

enr_res_s2@result[1:40,c('Description','p.adjust')]

#proteolysis p-value: .033

## Tc cells
covid_sub <- sig_df_covid[genes_both,'Tc']
covid_sig_up <- names(covid_sub)[covid_sub>covid_thresh] 
covid_sig_down <- names(covid_sub)[covid_sub<(-covid_thresh)] 

sle_sub <- sig_df_SLE[genes_both,'Tc']
sle_sig_up <- names(sle_sub)[sle_sub>sle_thresh]
sle_sig_down <- names(sle_sub)[sle_sub<(-sle_thresh)]

# length(intersect(covid_sig_up,sle_sig_up))

bg <- unique(c(covid_sig_up,sle_sig_up))
covid_unique_up <- covid_sig_up[!(covid_sig_up %in% sle_sig_up)]
sle_unique_up <- sle_sig_up[!(sle_sig_up %in% covid_sig_up)]
# bg <- unique(c(covid_sig_down,sle_sig_down))
# covid_unique_up <- covid_sig_down[!(covid_sig_down %in% sle_sig_down)]
# sle_unique_up <- sle_sig_down[!(sle_sig_down %in% covid_sig_down)]

## writing the genes to csvs
sig_both <- intersect(covid_sig_up,sle_sig_up)
write.csv(sig_both,file='/home/jmitchel/data/gene_lists/Tc_sig_both.csv')
write.csv(covid_unique_up,file='/home/jmitchel/data/gene_lists/Tc_sig_covid.csv')
write.csv(sle_unique_up,file='/home/jmitchel/data/gene_lists/Tc_sig_sle.csv')

enr_res_c3 <- enrichGO(covid_unique_up,
                    keyType = 'SYMBOL',
                    universe=bg,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    minGSSize = 10)

enr_res_c3@result[1:10,c('Description','p.adjust')]

enr_res_s3 <- enrichGO(sle_unique_up,
                    keyType = 'SYMBOL',
                    universe=bg,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    minGSSize = 10)

enr_res_s3@result[1:10,c('Description','p.adjust')]
# sets to show later:
pth <- c('ATP metabolic process','carboxylic acid metabolic process','mitotic cell cycle process',
         'Fc receptor signaling pathway')
enr_res_s3@result[enr_res_s3@result$Description %in% pth,]

# proteolysis pval: .51

# save relevant genes
tmp <- as.data.frame(enr_res_s3@result[enr_res_s3@result$Description %in% pth,'geneID'])
write.csv(tmp,file='/home/jmitchel/data/gene_lists/Tc_gsets.csv')

# genes to show
# ENO1, SDHC, COA6, NDUFB3, UQCRC1, NDUFB4, NDUFS
# SDHB, PSMB2, CTPS1, SCP2, CYP2J2, PSMA5, MDH1, UGP2, ATIC, NIT2, QDPR, QKI, PSMB1, AIMP


## B cells
covid_sub <- sig_df_covid[genes_both,'B']
covid_sig_up <- names(covid_sub)[covid_sub>covid_thresh] 
covid_sig_down <- names(covid_sub)[covid_sub<(-covid_thresh)] 

sle_sub <- sig_df_SLE[genes_both,'B']
sle_sig_up <- names(sle_sub)[sle_sub>sle_thresh]
sle_sig_down <- names(sle_sub)[sle_sub<(-sle_thresh)]

bg <- unique(c(covid_sig_up,sle_sig_up))
covid_unique_up <- covid_sig_up[!(covid_sig_up %in% sle_sig_up)]
sle_unique_up <- sle_sig_up[!(sle_sig_up %in% covid_sig_up)]
# bg <- unique(c(covid_sig_down,sle_sig_down))
# covid_unique_up <- covid_sig_down[!(covid_sig_down %in% sle_sig_down)]
# sle_unique_up <- sle_sig_down[!(sle_sig_down %in% covid_sig_down)]

## writing the genes to csvs
sig_both <- intersect(covid_sig_up,sle_sig_up)
write.csv(sig_both,file='/home/jmitchel/data/gene_lists/B_sig_both.csv')
write.csv(covid_unique_up,file='/home/jmitchel/data/gene_lists/B_sig_covid.csv')
write.csv(sle_unique_up,file='/home/jmitchel/data/gene_lists/B_sig_sle.csv')

enr_res_c4 <- enrichGO(covid_unique_up,
                    keyType = 'SYMBOL',
                    universe=bg,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    minGSSize = 10)

enr_res_c4@result[1:10,c('Description','p.adjust')]

enr_res_s4 <- enrichGO(sle_unique_up,
                    keyType = 'SYMBOL',
                    universe=bg,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    minGSSize = 10)

enr_res_s4@result[1:10,c('Description','p.adjust')]
# sets to show later:
pth <- c('regulation of MAP kinase activity')
enr_res_s4@result[enr_res_s4@result$Description %in% pth,]

#proteolysis pval: .97

# save relevant genes
tmp <- as.data.frame(enr_res_s4@result[enr_res_s4@result$Description %in% pth,'geneID'])
write.csv(tmp,file='/home/jmitchel/data/gene_lists/B_gsets.csv')

# genes to show
# RGS2, AIDA, CXCR4, DUSP1, STK38, DNAJA1, GNG3, DUSP5, DUSP6, RASGRP1, UBB, GADD45B


## cMono cells
covid_sub <- sig_df_covid[genes_both,'cMono']
covid_sig_up <- names(covid_sub)[covid_sub>covid_thresh] 
covid_sig_down <- names(covid_sub)[covid_sub<(-covid_thresh)] 
covid_sig_all <- names(covid_sub)[abs(covid_sub)>covid_thresh] 

sle_sub <- sig_df_SLE[genes_both,'cMono']
sle_sig_up <- names(sle_sub)[sle_sub>sle_thresh]
sle_sig_down <- names(sle_sub)[sle_sub<(-sle_thresh)]
sle_sig_all <- names(sle_sub)[abs(sle_sub)>sle_thresh] 

bg <- unique(c(covid_sig_up,sle_sig_up))
covid_unique_up <- covid_sig_up[!(covid_sig_up %in% sle_sig_up)]
sle_unique_up <- sle_sig_up[!(sle_sig_up %in% covid_sig_up)]

# bg <- unique(c(covid_sig_down,sle_sig_down))
# covid_unique_up <- covid_sig_down[!(covid_sig_down %in% sle_sig_down)]
# sle_unique_up <- sle_sig_down[!(sle_sig_down %in% covid_sig_down)]

bg <- unique(c(covid_sig_all,sle_sig_all))
covid_unique_up <- covid_sig_all[!(covid_sig_all %in% sle_sig_all)]
sle_unique_up <- sle_sig_all[!(sle_sig_all %in% covid_sig_all)]

## writing the genes to csvs
sig_both <- intersect(covid_sig_up,sle_sig_up)
write.csv(sig_both,file='/home/jmitchel/data/gene_lists/cMono_sig_both.csv')
write.csv(covid_unique_up,file='/home/jmitchel/data/gene_lists/cMono_sig_covid.csv')
write.csv(sle_unique_up,file='/home/jmitchel/data/gene_lists/cMono_sig_sle.csv')

enr_res_c5 <- enrichGO(covid_unique_up,
                    keyType = 'SYMBOL',
                    universe=bg,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    minGSSize = 10)

enr_res_c5@result[1:10,c('Description','p.adjust')]

enr_res_s5 <- enrichGO(sle_unique_up,
                    keyType = 'SYMBOL',
                    universe=bg,
                    OrgDb = org.Hs.eg.db,
                    ont = "BP",
                    minGSSize = 10)

enr_res_s5@result[1:10,c('Description','p.adjust')]
# sets to show later:
pth <- c('actin filament organization')
enr_res_s5@result[enr_res_s5@result$Description %in% pth,]

# save relevant genes
tmp <- as.data.frame(enr_res_s5@result[enr_res_s5@result$Description %in% pth,'geneID'])
write.csv(tmp,file='/home/jmitchel/data/gene_lists/cMono_gsets.csv')

#proteolysis pval: .089

# genes to show
# CAPZB, CDC42, RHOC, S100A10, TPM3, HAX1, ARPC5, ARF1, PLEK, TMSB10, CAPG, ACTR3,
# WIPF1, RHOA, CTNNA1, AIF1, EZR, RAC1, DBNL, CAPZA2


## making per cell type heatmap of enrichment results
all_sets <- unique(c('actin filament organization','regulation of MAP kinase activity',
              'ATP metabolic process','carboxylic acid metabolic process','mitotic cell cycle process',
              'Fc receptor signaling pathway','leukocyte cell-cell adhesion', 'tissue development', 'apoptotic process', 
              'cell activation', 'leukocyte migration',
              'carboxylic acid metabolic process', 'Fc receptor signaling pathway'))

enr_res <- matrix(1,nrow=length(all_sets),ncol=5)
rownames(enr_res) <- all_sets
colnames(enr_res) <- c('NK','Th','Tc','B','cMono')
res_lst <- list(enr_res_s1,enr_res_s2,enr_res_s3,enr_res_s4,enr_res_s5)
# populate results matrix
for (i in 1:ncol(enr_res)) {
  print(i)
  print('')
  res <- res_lst[[i]]@result
  for (j in 1:nrow(enr_res)) {
    print(j)
    gs <- rownames(enr_res)[j]
    if (gs %in% res$Description) {
      padj <- res[res$Description %in% gs,'p.adjust']
      enr_res[j,i] <- padj
    }
  }
}

col_fun <- colorRamp2(c(-log10(.15), 4), c("white", "red"))
h_w <- c(10,10)
enr_res <- -log10(enr_res)
myhmap <- Heatmap(as.matrix(enr_res), name = '-log10(padj)',
                  cluster_rows = TRUE,
                  cluster_columns = TRUE,
                  show_row_dend = FALSE, show_column_dend = FALSE,
                  column_names_gp = gpar(fontsize = 12),
                  col = col_fun,
                  row_title_gp = gpar(fontsize = 12),
                  column_title = 'Cell Types',
                  column_title_side = "bottom",
                  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                  row_title = 'Gene Sets',
                  row_title_side = "left",
                  border=TRUE,
                  row_names_side = "right",
                  row_names_gp = gpar(fontsize = 12),
                  show_heatmap_legend = TRUE,
                  width = unit(h_w[2], "cm"),
                  height = unit(h_w[1], "cm"),
                  cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                    if (abs(enr_res[i,j]) > -log10(.001)) {
                      grid.text('***', x, y, gp = gpar(fontface='bold'))
                    } else if (abs(enr_res[i,j]) > -log10(.01)) {
                      grid.text('**', x, y, gp = gpar(fontface='bold'))
                    } else if (abs(enr_res[i,j]) > -log10(.05)) {
                      grid.text('*', x, y, gp = gpar(fontface='bold'))
                    } else {
                      grid.text('', x, y)
                    }
                  })

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_sle_comp.pdf", useDingbats = FALSE,
#     width = 9, height = 7)
myhmap
dev.off()









require("ggrepel")

tmp <- cbind.data.frame(sig_df_covid[genes_both,'cMono'],sig_df_SLE[genes_both,'cMono'])
# adding gene labels for callouts
colnames(tmp) <- c('covid','sle')
tmp$callouts <- rep("",nrow(tmp))
tmp['MX1','callouts'] <- 'MX1'
tmp['IFI6','callouts'] <- 'IFI6'
tmp['XAF1','callouts'] <- 'XAF1'
tmp['IFIT3','callouts'] <- 'IFIT3'
g_include <- c('S100A10','PLEK','RHOA','EZR')
tmp[g_include,'callouts'] <- g_include
tmp$is_callout <- sapply(tmp$callouts,function(x){
  if (x=='') {
    return(1)
  } else {
    return(2)
  }
})

p1 <- ggplot(tmp,aes(x=covid,y=sle,color=as.factor(is_callout))) +
  geom_point(alpha=.2,pch=19) +
  xlab('Factor 1 T-statistics (Covid dataset)') +
  ylab('Factor 1 T-statistics (SLE dataset)') +
  geom_vline(xintercept=0, linetype="solid", 
             color = "black", size=1) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  geom_hline(yintercept=sle_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_hline(yintercept=-sle_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_vline(xintercept=covid_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_vline(xintercept=-covid_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  ggtitle('cMonocytes') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_manual(values = c("black","red")) +
  geom_text_repel(
    data = tmp,
    aes(label = callouts),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps=2000
  )

p1

tmp <- cbind.data.frame(sig_df_covid[genes_both,'Th'],sig_df_SLE[genes_both,'Th'])
# adding gene labels for callouts
colnames(tmp) <- c('covid','sle')
tmp$callouts <- rep("",nrow(tmp))
tmp['MX1','callouts'] <- 'MX1'
tmp['IFI6','callouts'] <- 'IFI6'
tmp['XAF1','callouts'] <- 'XAF1'
tmp['IFIT3','callouts'] <- 'IFIT3'
tmp$is_callout <- sapply(tmp$callouts,function(x){
  if (x=='') {
    return(1)
  } else {
    return(2)
  }
})
p2 <- ggplot(tmp,aes(x=covid,y=sle,color=as.factor(is_callout))) +
  geom_point(alpha=.2,pch=19) +
  xlab('Factor 1 T-statistics (Covid dataset)') +
  ylab('Factor 1 T-statistics (SLE dataset)') +
  geom_vline(xintercept=0, linetype="solid", 
             color = "black", size=1) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  geom_hline(yintercept=sle_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_hline(yintercept=-sle_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_vline(xintercept=covid_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_vline(xintercept=-covid_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  ggtitle('Th cells') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_manual(values = c("black","red")) +
  geom_text_repel(
    data = tmp,
    aes(label = callouts),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps=2000
  )

p2


tmp <- cbind.data.frame(sig_df_covid[genes_both,'Tc'],sig_df_SLE[genes_both,'Tc'])
# adding gene labels for callouts
colnames(tmp) <- c('covid','sle')
tmp$callouts <- rep("",nrow(tmp))
tmp['MX1','callouts'] <- 'MX1'
tmp['IFI6','callouts'] <- 'IFI6'
tmp['XAF1','callouts'] <- 'XAF1'
tmp['IFIT3','callouts'] <- 'IFIT3'
g_include <- c('PSMB2','PSMA5','MDH1','UGP2')
tmp[g_include,'callouts'] <- g_include
tmp$is_callout <- sapply(tmp$callouts,function(x){
  if (x=='') {
    return(1)
  } else {
    return(2)
  }
})
p3 <- ggplot(tmp,aes(x=covid,y=sle,color=as.factor(is_callout))) +
  geom_point(alpha=.2,pch=19) +
  xlab('Factor 1 T-statistics (Covid dataset)') +
  ylab('Factor 1 T-statistics (SLE dataset)') +
  geom_vline(xintercept=0, linetype="solid", 
             color = "black", size=1) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  geom_hline(yintercept=sle_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_hline(yintercept=-sle_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_vline(xintercept=covid_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_vline(xintercept=-covid_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  ggtitle('Tc cells') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_manual(values = c("black","red")) +
  geom_text_repel(
    data = tmp,
    aes(label = callouts),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps=2000
  )

p3


tmp <- cbind.data.frame(sig_df_covid[genes_both,'B'],sig_df_SLE[genes_both,'B'])
# adding gene labels for callouts
colnames(tmp) <- c('covid','sle')
tmp$callouts <- rep("",nrow(tmp))
tmp['MX1','callouts'] <- 'MX1'
tmp['IFI6','callouts'] <- 'IFI6'
tmp['XAF1','callouts'] <- 'XAF1'
tmp['IFIT3','callouts'] <- 'IFIT3'
g_include <- c('DUSP5','DUSP1','UBB')
tmp[g_include,'callouts'] <- g_include
tmp$is_callout <- sapply(tmp$callouts,function(x){
  if (x=='') {
    return(1)
  } else {
    return(2)
  }
})
p4 <- ggplot(tmp,aes(x=covid,y=sle,color=as.factor(is_callout))) +
  geom_point(alpha=.2,pch=19) +
  xlab('Factor 1 T-statistics (Covid dataset)') +
  ylab('Factor 1 T-statistics (SLE dataset)') +
  geom_vline(xintercept=0, linetype="solid", 
             color = "black", size=1) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  geom_hline(yintercept=sle_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_hline(yintercept=-sle_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_vline(xintercept=covid_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_vline(xintercept=-covid_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  ggtitle('B cells') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_manual(values = c("black","red")) +
  geom_text_repel(
    data = tmp,
    aes(label = callouts),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps=2000
  )

p4

tmp <- cbind.data.frame(sig_df_covid[genes_both,'NK'],sig_df_SLE[genes_both,'NK'])
# adding gene labels for callouts
colnames(tmp) <- c('covid','sle')
tmp$callouts <- rep("",nrow(tmp))
tmp['MX1','callouts'] <- 'MX1'
tmp['IFI6','callouts'] <- 'IFI6'
tmp['XAF1','callouts'] <- 'XAF1'
tmp['IFIT3','callouts'] <- 'IFIT3'
g_include <- c('HAVCR2','ITGB7','LGALS3','PSME1')
tmp[g_include,'callouts'] <- g_include
tmp$is_callout <- sapply(tmp$callouts,function(x){
  if (x=='') {
    return(1)
  } else {
    return(2)
  }
})
p5 <- ggplot(tmp,aes(x=covid,y=sle,color=as.factor(is_callout))) +
  geom_point(alpha=.2,pch=19) +
  xlab('Factor 1 T-statistics (Covid dataset)') +
  ylab('Factor 1 T-statistics (SLE dataset)') +
  geom_vline(xintercept=0, linetype="solid", 
             color = "black", size=1) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  geom_hline(yintercept=sle_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_hline(yintercept=-sle_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_vline(xintercept=covid_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  geom_vline(xintercept=-covid_thresh, linetype="dashed", 
             color = "blue", size=.5, alpha=.5) +
  ggtitle('NK cells') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  scale_color_manual(values = c("black","red")) +
  geom_text_repel(
    data = tmp,
    aes(label = callouts),
    size = 5,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps=2000
  )

p5

library(cowplot)
# row1 <- plot_grid(p5,p3,p2,nrow=1)
# row2 <- plot_grid(plotlist=list(NULL,NULL,p4,p1,NULL,NULL),
#                   nrow=1,rel_widths = c(.25,.25,1,1.25,.5))
# 
# fig <- plot_grid(row1,row2,ncol=1)
# 
# fig <- plot_grid(p5,p3,p2,p4,p1,ncol=2,align = 'hv')

fig <- plot_grid(p3,p1,nrow=1)

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_sle_dotplots2.pdf", useDingbats = FALSE,
    width = 9.5, height = 6)
fig
dev.off()


fig <- plot_grid(p2,p4,nrow=1)

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_sle_dotplots3.pdf", useDingbats = FALSE,
    width = 9.5, height = 6)
fig
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_sle_dotplots4.pdf", useDingbats = FALSE,
    width = 4.75, height = 6)
p5
dev.off()





## cor vals for comparing the two covid datasets
# cMono: .77, Th: .68, Tc: .56, B:, .65, NK: .61

cool1 <- get_one_factor(pbmc_container,1)[[2]]
cool2 <- get_one_factor(pbmc_container_SLE,1)[[2]]
gb <- intersect(rownames(cool1),rownames(cool2))
cor(cool1[gb,1],cool2[gb,1])
cor(cool1[gb,2],cool2[gb,2])
cor(cool1[gb,3],cool2[gb,3])
cor(cool1[gb,4],cool2[gb,4])
cor(cool1[gb,5],cool2[gb,5])





## getting number of same significant genes for table
get_num_sig <- function(ctype) {
  covid_sub <- sig_df_covid2[genes_both,ctype]
  covid_sig <- names(covid_sub)[covid_sub<.01] 
  
  sle_sub <- sig_df_SLE2[genes_both,ctype]
  sle_sig <- names(sle_sub)[sle_sub<.01]
  
  print(length(intersect(covid_sig,sle_sig)))
  print(length(covid_sig[!(covid_sig %in% sle_sig)]))
  print(length(sle_sig[!(sle_sig %in% covid_sig)]))
}

get_num_sig('Th')
get_num_sig('Tc')
get_num_sig('NK')
get_num_sig('B')
get_num_sig('cMono')













## showing enrichment of proteolysis process in the thing

# from inside run_fgsea because need these values to make the enrichment plot
get_paths_ranks <- function(container, factor_select, ctype, db_use="GO", signed=TRUE, 
                      max_gs_size=500, ncores=container$experiment_params$ncores) {
  donor_scores <- container$tucker_results[[1]]
  
  # select mean exp data for one cell type
  tnsr_slice <- container$scMinimal_ctype[[ctype]]$pseudobulk
  tnsr_slice <- scale(tnsr_slice, center=TRUE) # rescaling to unit variance
  
  # get transformed expression for each gene by summing d_score * scaled exp
  exp_vals <- sapply(1:ncol(tnsr_slice), function(j) {
    if (signed) {
      exp_transform <- tnsr_slice[,j] * donor_scores[rownames(tnsr_slice),factor_select]
      de_val <- sum(exp_transform)
    } else {
      # testing out using undirected statistics
      exp_transform <- tnsr_slice[,j] * donor_scores[rownames(tnsr_slice),factor_select]
      de_val <- abs(sum(exp_transform))
    }
    
    return(de_val)
  })
  
  names(exp_vals) <- convert_gn(container,colnames(tnsr_slice))
  
  # remove duplicate genes
  ndx_remove <- duplicated(names(exp_vals)) | duplicated(names(exp_vals), fromLast = TRUE)
  exp_vals <- exp_vals[!ndx_remove]
  
  m_df <- data.frame()
  for (db in db_use) {
    if (db == "GO") {
      # select the GO Biological Processes group of gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C5", subcategory = "BP"))
    } else if (db == "GOCC") {
      # select the GOCC gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C5", subcategory = "CC"))
    } else if (db == "Reactome") {
      # select the Reactome gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C2", subcategory = "CP:REACTOME"))
    } else if (db == "KEGG") {
      # select the KEGG gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C2", subcategory = "CP:KEGG"))
    } else if (db == "BioCarta") {
      # select the BioCarts gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C2", subcategory = "CP:BIOCARTA"))
    } else if (db == "Hallmark") {
      # select the BioCarts gene sets
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "H"))
    } else if (db == "TF") {
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C3", subcategory = "TFT:GTRD"))
      m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                          category = "C3", subcategory = "TFT:TFT_Legacy"))
    }
  }
  
  my_pathways <- split(m_df$gene_symbol, f = m_df$gs_name)
  
  return(list(my_pathways,exp_vals))
}

# for covid data first
go_covid <- run_fgsea(pbmc_container, factor_select=1, ctype='Tc', 
                      db_use="GOCC", signed=TRUE, min_gs_size=5, ncores=10)
covid_paths_ranks <- get_paths_ranks(pbmc_container, factor_select=1, ctype='Tc', 
                            db_use="GOCC", signed=TRUE, ncores=10)
mypaths <- covid_paths_ranks[[1]]
myranks <- covid_paths_ranks[[2]]
myranks <- myranks[!is.na(myranks)]

plt_covid <- fgsea::plotEnrichment(mypaths[['GOCC_ENDOPEPTIDASE_COMPLEX']],
                             myranks) + labs(title='ENDOPEPTIDASE_COMPLEX - COVID-19 Factor 1')

plt_covid <- plt_covid +
  annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
           label=paste0('p-value: ',signif(go_covid[go_covid$pathway=='GOCC_ENDOPEPTIDASE_COMPLEX','pval'],digits=2)))

# flip sign of F1 so high ISG expression is positive instead of negative (signs are arbitrary)
pbmc_container_SLE$tucker_results[[1]][,1] <- pbmc_container_SLE$tucker_results[[1]][,1] * -1
pbmc_container_SLE$tucker_results[[2]][1,] <- pbmc_container_SLE$tucker_results[[2]][1,] * -1
go_SLE <- run_fgsea(pbmc_container_SLE, factor_select=1, ctype='Tc', 
                    db_use="GOCC", signed=TRUE, min_gs_size=5, ncores=10)
SLE_paths_ranks <- get_paths_ranks(pbmc_container_SLE, factor_select=1, ctype='Tc', 
                                     db_use="GOCC", signed=TRUE, ncores=10)
mypaths <- SLE_paths_ranks[[1]]
myranks <- SLE_paths_ranks[[2]]
myranks <- myranks[!is.na(myranks)]

plt_sle <- fgsea::plotEnrichment(mypaths[['GOCC_ENDOPEPTIDASE_COMPLEX']],
                             myranks) + labs(title='ENDOPEPTIDASE_COMPLEX - SLE Factor 1')

plt_sle <- plt_sle +
  annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
           label=paste0('p-value: ',signif(go_SLE[go_SLE$pathway=='GOCC_ENDOPEPTIDASE_COMPLEX','pval'],digits=2)))

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_proteolysis.pdf", useDingbats = FALSE,
    width = 3.5, height = 3)
plt_covid
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_proteolysis.pdf", useDingbats = FALSE,
    width = 3.5, height = 3)
plt_sle
dev.off()






## plotting fraction of SLE/covid-only genes in proteolysis set
pathways <- msigdbr::msigdbr(species = "Homo sapiens",category = "C5", subcategory = "CC")
my_pathways <- split(pathways$gene_symbol, f = pathways$gs_name)
prot_set <- my_pathways[['GOCC_ENDOPEPTIDASE_COMPLEX']]
# prot_set <- my_pathways[['GOBP_PROTEOLYSIS']]

get_fracs_pval <- function(ctype) {
  covid_sub <- sig_df_covid[genes_both,ctype]
  covid_sig_up <- names(covid_sub)[covid_sub>covid_thresh] 
  
  sle_sub <- sig_df_SLE[genes_both,ctype]
  sle_sig_up <- names(sle_sub)[sle_sub>sle_thresh]
  
  covid_unique_up <- covid_sig_up[!(covid_sig_up %in% sle_sig_up)]
  sle_unique_up <- sle_sig_up[!(sle_sig_up %in% covid_sig_up)]
  bg <- unique(c(covid_unique_up,sle_unique_up))
  
  covid_frac <- sum(covid_unique_up %in% prot_set) / length(covid_unique_up)
  SLE_frac <- sum(sle_unique_up %in% prot_set) / length(sle_unique_up)
  
  # get enrichment p-value
  enr_res <- enrichGO(sle_unique_up,
                      keyType = 'SYMBOL',
                      universe=bg,
                      OrgDb = org.Hs.eg.db,ont = "CC",
                      minGSSize = 2)
  pval <- enr_res@result[enr_res@result$Description=='endopeptidase complex','pvalue']
  return(list(covid_frac,SLE_frac,pval))
}

res <- matrix(nrow=10,ncol=3)
colnames(res) <- c('ctype','dataset','frac')

res1 <- get_fracs_pval('Tc')
res[1,] <- c('Tc','COVID-19',res1[[1]])
res[2,] <- c('Tc','SLE',res1[[2]])
res1[[3]]

res2 <- get_fracs_pval('Th')
res[3,] <- c('Th','COVID-19',res2[[1]])
res[4,] <- c('Th','SLE',res2[[2]])
res2[[3]]

res3 <- get_fracs_pval('NK')
res[5,] <- c('NK','COVID-19',res3[[1]])
res[6,] <- c('NK','SLE',res3[[2]])
res3[[3]]

res4 <- get_fracs_pval('B')
res[7,] <- c('B','COVID-19',res4[[1]])
res[8,] <- c('B','SLE',res4[[2]])
res4[[3]]

res5 <- get_fracs_pval('cMono')
res[9,] <- c('cMono','COVID-19',res5[[1]])
res[10,] <- c('cMono','SLE',res5[[2]])
res5[[3]]

res <- as.data.frame(res)

res$ctype <- as.factor(res$ctype)
res$dataset <- as.factor(res$dataset)
res$frac <- as.numeric(res$frac)

# adding very small value where 0, so it shows up on the plot
res$frac[res$frac==0] <- .0001

# plotting results
p <- ggplot(res,aes(x=ctype,y=frac,fill=dataset)) +
  geom_bar(stat="identity", position=position_dodge(preserve = "single")) +
  ylab('Fraction of dataset-specific F1 genes') +
  xlab ('') +
  theme_bw()

pdf(file = "/home/jmitchel/figures/for_paper_v2/proteasome_comp.pdf", useDingbats = FALSE,
    width = 6, height = 2.8)
p
dev.off()







# ## computing and comparing IFN TF activity and targets between covid and SLE
# # first, rerun lm to get pvalues in correct slot
# pbmc_container <- get_lm_pvals(pbmc_container)
# pbmc_container_SLE <- get_lm_pvals(pbmc_container_SLE)
# 
# # run unsigned TF gsea for each dataset
# go_covid_Tc <- run_hypergeometric_gsea(pbmc_container, factor_select=1, ctype='cMono', up_down='unsigned',
#                         thresh=0.01, db_use="TF")
# go_SLE_Tc <- run_hypergeometric_gsea(pbmc_container_SLE, factor_select=1, ctype='cMono', up_down='unsigned',
#                                        thresh=0.01, db_use="TF")
# go_covid_Tc_2 <- run_hypergeometric_gsea(pbmc_ye, factor_select=1, ctype='cMono', up_down='unsigned',
#                                      thresh=0.25, db_use="TF")
# 
# go_covid_Tc[order(go_covid_Tc)][1:10]
# go_SLE_Tc[order(go_SLE_Tc)][1:10]
# go_covid_Tc_2[order(go_covid_Tc_2)][1:10]
# 
# 
# # extracting leading edge genes
# regulons <- dorothea_hs %>%
#   filter(confidence %in% c("A", "B", "C"))
# colnames(regulons)[1] <- c('gs_name')
# colnames(regulons)[3] <- c('gene_symbol')
# 
# sig_both <- intersect(names(go_covid_Tc)[go_covid_Tc<.05],names(go_SLE_Tc)[go_SLE_Tc<.05])
# sig_both
# all_targets <- regulons$gene_symbol[regulons$gs_name=='STAT2']
# le_1 <- all_targets[all_targets %in% covid_sig_all]
# le_1 <- le_1[le_1 %in% genes_both]
# le_2 <- all_targets[all_targets %in% sle_sig_all]
# le_2 <- le_2[le_2 %in% genes_both]
# length(intersect(le_1,le_2))/length(union(le_1,le_2))
# 
# 
# go_covid <- run_fgsea(pbmc_container, factor_select=1, ctype='cMono', 
#                     db_use="TF", signed=FALSE, ncores=10)
# go_SLE <- run_fgsea(pbmc_container_SLE, factor_select=1, ctype='cMono', 
#                        db_use="TF", signed=FALSE, ncores=10)
# go_ye <- run_fgsea(pbmc_ye, factor_select=1, ctype='cMono', 
#                     db_use="TF", signed=FALSE, ncores=10)
# go_covid <- go_covid[go_covid$padj<.05,]
# go_SLE <- go_SLE[go_SLE$padj<.05,]
# go_ye <- go_ye[go_ye$padj<.05,]
# 
# sig_both <- intersect(go_covid$pathway,go_SLE$pathway)
# sig_both <- intersect(go_covid$pathway,go_ye$pathway)
# sig_both
# # calculate jaccard coefficient of overlap
# le_1 <- go_covid[go_covid$pathway=='PSMB5_TARGET_GENES','leadingEdge'][[1]][[1]]
# le_1 <- le_1[le_1 %in% genes_both]
# 
# le_2 <- go_SLE[go_SLE$pathway=='IRF1','leadingEdge'][[1]][[1]]
# le_2 <- le_2[le_2 %in% genes_both]
# 
# le_2 <- go_ye[go_ye$pathway=='PSMB5_TARGET_GENES','leadingEdge'][[1]][[1]]
# le_2 <- le_2[le_2 %in% genes_both]
# 
# length(intersect(le_1,le_2))/length(union(le_1,le_2))
# 
# # probably should have removed genes w 0 expression for significance estimation...
# # for validation, should try this with the Jimmie Ye dataset and should get higher matched leadingEdge sets
# # also critical that I make sure the same sets of genes are being compared...
# 
# # determine consistently activated TFs
# 
# # get their leading edge genes
# 
# # compute jaccard coefficients
# 
# # get top 5 TFs with highest target gene overlap
# 
# # get top 5 TFs with least target gene overlap



















##### before comparing the two covid datasets, I need to reform tensor and run
# decomp using only top od genes
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=2, 
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.0001,
                              batch_var = 'site',
                              scale_var = TRUE, var_scale_power = .5) 

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(9,38), 
                                 tucker_type = 'regular', rotation_type = 'hybrid') # best with vargenes_thresh=.0001 also .01

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','age',
                                                                   'status_on_day_collection_summary',
                                                                   'site',"smoker"),stat_use='pval')
# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval',h_w=c(10,8))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_dscores2.pdf", useDingbats = FALSE,
#     width = 8, height = 8)
pbmc_container$plots$donor_matrix
# dev.off()

# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)


##### prepping jimmie's covid data for comparison
pbmc <- readRDS(file='/home/jmitchel/data/covid_data_ye/covid_subset.rds')
# making pool variable not have batches without representation
pbmc@meta.data$pool <- factor(pbmc@meta.data$pool,levels=unique(pbmc@meta.data$pool))

# changing cell type names to match others used
new_names <- sapply(as.character(pbmc@meta.data$ct1), function(x){
  if (x=='cM') {
    return('cMono')
  } else if (x=='ncM') {
    return('ncMono')
  } else if (x=='T4') {
    return('Th')
  } else if (x=='T8') {
    return('Tc')
  } else {
    return(x)
  }
})
names(new_names) <- NULL
pbmc@meta.data$ct1 <- factor(new_names,levels=unique(new_names))

param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc",
                                               "cMono"),
                                ncores = 30, rand_seed = 10)
pbmc_ye <- make_new_container(seurat_obj=pbmc,
                                     params=param_list,
                                     metadata_cols=c('donor',
                                                     "COVID_status",
                                                     'COVID_severity',
                                                     "COVID_severity_merged",
                                                     "ct1",
                                                     "sex",
                                                     "pool",
                                                     "ethnicity",
                                                     "respiratory_support_D0",
                                                     'onset_to_D0_days',
                                                     'intubated_days',
                                                     'admission_to_discharge',
                                                     'death',
                                                     'pulmonary_infection',
                                                     'non_pulmonary_infection',
                                                     'WBC_count1',
                                                     'WBC_count2',
                                                     'WBC_count3',
                                                     'Monocyte_count',
                                                     'respiratory_support',
                                                     'NIH_clinical',
                                                     'NIH_ordinal',
                                                     'admission_level'),
                                     metadata_col_nm=c('donors',
                                                       'covid_status',
                                                       'covid_severity',
                                                       'covid_severity_merged',
                                                       'ctypes',
                                                       'sex',
                                                       'pool',
                                                       'Ethnicity',
                                                       'respiratory_support_D0',
                                                       'onset_to_D0_days',
                                                       'intubated_days',
                                                       'admission_to_discharge',
                                                       'death',
                                                       'pulmonary_infection',
                                                       'non_pulmonary_infection',
                                                       'WBC_count1',
                                                       'WBC_count2',
                                                       'WBC_count3',
                                                       'Monocyte_count',
                                                       'respiratory_support',
                                                       'NIH_clinical',
                                                       'NIH_ordinal',
                                                       'admission_level'))

pbmc_ye <- form_tensor(pbmc_ye, donor_min_cells=2, 
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')

pbmc_ye <- run_tucker_ica(pbmc_ye, ranks=c(10,30), # tested well with comparison below
                                 tucker_type = 'regular', rotation_type = 'hybrid')

# get factor-meta data associations
pbmc_ye <- get_meta_associations(pbmc_ye,vars_test=c('sex','pool',
                                                                   'covid_status',
                                                                   'covid_severity',
                                                                   'covid_severity_merged',
                                                                   'death'),
                                        stat_use='pval')

# plot donor scores by status
pbmc_ye <- plot_donor_matrix(pbmc_ye, meta_vars=c('sex'),
                                    cluster_by_meta = 'sex',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# pdf(file = "/home/jmitchel/figures/for_paper_v2/lupus_dscores_full2.pdf", useDingbats = FALSE,
#     width = 6, height = 7)
pbmc_ye$plots$donor_matrix
# dev.off()

pbmc_ye <- get_lm_pvals(pbmc_ye)

# getting count of covid and healthy patients
pbmc_ye <- get_donor_meta(pbmc_ye,'covid_status')
table(pbmc_ye$donor_metadata$covid_status)

# getting # cells used from the major clusters
d_kept <- rownames(pbmc_ye$scMinimal_ctype[[1]]$pseudobulk)
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$donor %in% d_kept]
pbmc <- subset(pbmc,cells=cells_keep)

ct_kept <- pbmc_ye$experiment_params$ctypes_use
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$ct1 %in% ct_kept]
pbmc <- subset(pbmc,cells=cells_keep)
dim(pbmc)


##### get median number cells per donor and total number of cells used
d_keep <- rownames(pbmc_ye$scMinimal_ctype[[1]]$pseudobulk)
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$donor %in% d_keep]
pbmc_sub <- subset(pbmc,cells=cells_keep)
ctypes <- c("B","Tc","Th","NK","cMono")
cells_keep <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$ct1 %in% ctypes]
pbmc_sub2 <- subset(pbmc_sub,cells=cells_keep)
pbmc_sub2@meta.data$donor <- factor(as.character(pbmc_sub2@meta.data$donor),levels=unique(as.character(pbmc_sub2@meta.data$donor)))
for (ctype in ctypes) {
  tmp <- pbmc_sub2@meta.data[pbmc@meta.data$ct1==ctype,]
  print(ctype)
  print(median(table(tmp$donor)))
}





## projecting UK factors onto Jimmie covid data
pbmc_ye <- project_new_data(pbmc_ye,pbmc_container)
cor(pbmc_ye$projected_scores,pbmc_ye$tucker_results[[1]])

# now plotting projected scores versus metadata
dsc <- pbmc_ye$projected_scores[,2,drop=FALSE]
pbmc_ye <- get_donor_meta(pbmc_ye,additional_meta = c('covid_status',
                                                                    'covid_severity_merged'),only_analyzed = TRUE)
tmp <- cbind.data.frame(dsc,pbmc_ye$donor_metadata[rownames(dsc),'covid_severity_merged'],
                        pbmc_ye$donor_metadata[rownames(dsc),'covid_status'])
colnames(tmp) <- c('dscore','severity','status')
levels(tmp$severity)
levels(tmp$severity) <- c('Critical','Healthy','Moderate','Severe','Critical_anti-IFN','Negative')
tmp$severity <- factor(tmp$severity,levels=c('Healthy','Moderate','Severe','Critical','Critical_anti-IFN'))

mycol <- RColorBrewer::brewer.pal(n = 7, name = "Dark2")

tmp2 <- tmp
tmp2$status <- rep('violin',nrow(tmp2))
p <- ggplot(tmp,aes(x=severity,y=dscore,fill=status)) +
  geom_violin(data=tmp2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.85, binwidth = .01) +
  ylab('Projected Factor 2 Scores') +
  xlab('Severity on collection day') +
  coord_flip() +
  scale_fill_manual(values=c(mycol[2], mycol[5], 'light gray')) +
  theme_bw() +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f2_proj.pdf", useDingbats = FALSE,
#     width = 12.5, height = 8.5)
p
dev.off()

## calculate significance of the association
dsc <- pbmc_ye$projected_scores[,2,drop=FALSE]
tmp <- cbind.data.frame(dsc,
                        pbmc_ye$donor_metadata[rownames(dsc),'covid_severity_merged'])
colnames(tmp) <- c('proj_dsc','status')
lmres <- summary(lm(proj_dsc~status,data=tmp))
lmres
pv1 <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)

pv2 <- pbmc_container[["meta_associations"]]['status_on_day_collection_summary',2]


## combining severity association p-values for meta analysis
library(metap)
sumlog(c(pv1,pv2))
# p =  4.342494e-17 






##### now comparing loadings values for the corresponding severity associated factors
decomp1 <- get_one_factor(pbmc_container,2)
decomp2 <- get_one_factor(pbmc_ye,5)

lds1 <- decomp1[[2]]
lds2 <- decomp2[[2]]

genes_both <- intersect(rownames(lds1),rownames(lds2))
cor(lds1[genes_both,'cMono'],lds2[genes_both,'cMono'])
cor(lds1[genes_both,'B'],lds2[genes_both,'B'])
cor(lds1[genes_both,'Th'],lds2[genes_both,'Th'])
cor(lds1[genes_both,'Tc'],lds2[genes_both,'Tc'])
cor(lds1[genes_both,'NK'],lds2[genes_both,'NK'])

uk_lds <- get_one_factor(pbmc_container,2)[[2]]
# sig_vectors <- get_significance_vectors(pbmc_container,
#                                         factor_select=2, pbmc_container$experiment_params$ctypes_use)
# uk_lds <- t(as.data.frame(do.call(rbind, sig_vectors)))
ctypes <- pbmc_ye$experiment_params$ctypes_use
dsc <- pbmc_ye$tucker_results[[1]]
# matrix to store results
myres <- matrix(ncol=length(ctypes),nrow=ncol(dsc))
colnames(myres) <- ctypes
# loop through ye factors
for (f in 1:ncol(dsc)) {
  ye_lds <- get_one_factor(pbmc_ye,f)[[2]]
  # sig_vectors <- get_significance_vectors(pbmc_ye,
  #                                         factor_select=f, pbmc_ye$experiment_params$ctypes_use)
  # ye_lds <- t(as.data.frame(do.call(rbind, sig_vectors)))
  genes_both <- intersect(rownames(uk_lds),rownames(ye_lds))
  # loop through cell types
  for (ct in ctypes) {
    lds_cor <- cor(uk_lds[genes_both,ct],ye_lds[genes_both,ct],method='spearman')
    # uk_sub <- uk_lds[genes_both,ct,drop=F]
    # sig_uk <- rownames(uk_sub)[uk_sub[,1]<.01]
    # ye_sub <- ye_lds[genes_both,ct,drop=F]
    # sig_ye <- rownames(ye_sub)[ye_sub[,1]<.01]
    # sig_intersect <- intersect(sig_uk,sig_ye)
    # sig_union <- union(sig_uk,sig_ye)
    # lds_cor <- length(sig_intersect)/length(sig_union)
    myres[f,ct] <- lds_cor
  }
}

# COULD ALSO TRY USING JACCARD COEFFICIENT OF SIGNIFICANT GENES PER CTYPE INSTEAD OF COR
# multiply by -1 because ye ones are arbitrarily in opposite direction
# title: F2 associations w vdw factors
myres <- myres*-1
rownames(myres) <- sapply(c(1:nrow(myres)),function(x){paste0('Factor ',as.character(x))})

# add right annotation
col_fun_annot = colorRamp2(c(0, -log10(.05), 5), c("white", "white", "forest green"))
logpv <- -log10(pbmc_ye$meta_associations)
logpv <- logpv['covid_severity_merged',,drop=FALSE]
rownames(logpv)[1] <- 'severity'
ra <- rowAnnotation('-log_10_pval'=t(logpv),col = list('-log_10_pval'=col_fun_annot),
                        border=TRUE)

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
# col_fun <- colorRamp2(c(0, max(myres)), c("white", "red"))
h_w <- c(12,8)
myhmap <- Heatmap(as.matrix(myres), name = 'spearman cor',
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  show_row_dend = FALSE, show_column_dend = FALSE,
                  column_names_gp = gpar(fontsize = 12),
                  col = col_fun,
                  row_title_gp = gpar(fontsize = 12),
                  column_title = 'Cell Types',
                  column_title_side = "top",
                  column_names_side = "top",
                  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                  row_title = 'Factors (van der Wijst dataset)',
                  row_title_side = "left",
                  border=TRUE,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = 12),
                  show_heatmap_legend = TRUE,
                  width = unit(h_w[2], "cm"),
                  height = unit(h_w[1], "cm"),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid::grid.text(sprintf("%.2f", myres[i, j]), x, y, gp = gpar(fontsize = 10))
                  },
                  right_annotation = ra)

pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_severity_comp.pdf", useDingbats = FALSE,
    width = 9, height = 7)
myhmap
dev.off()








##### running LR analysis
# iTalk db
library(iTALK)
lr_pairs <- database[,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')]
colnames(lr_pairs) <- c('ligand','receptor')
lr_pairs <- unique(lr_pairs)

# infer active LR interactions
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='site')
sft_thresh <- c(12,12,12,12,10) #cchat
sft_thresh <- c(12,12,12,12,10) #italk
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

lr_hmap <- compute_LR_interact(pbmc_container, lr_pairs, sig_thresh=.0001,
                               percentile_exp_rec=0.85, add_ld_fact_sig=TRUE)

pdf(file = "/home/jmitchel/figures/test2.pdf", useDingbats = FALSE,
    width = 20, height = 85)
lr_hmap
dev.off()

# trying a smaller sig_thresh
sig_thresh=.0000001
myres_mat <- pbmc_container$lr_res # at 285
container=pbmc_container
myhmap1
pdf(file = "/home/jmitchel/figures/test2.pdf", useDingbats = FALSE,
    width = 14, height = 55)
myhmap1
dev.off()

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,
                                 mod_ct='NK',mod=12,lig_ct='Th',lig='PF4')
lig_mod_fact

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,
                                 mod_ct='Th',mod=2,lig_ct='Tc',lig='IL16')
lig_mod_fact

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,
                                 mod_ct='cMono',mod=13,lig_ct='Th',lig='PECAM1')
lig_mod_fact

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,
                                 mod_ct='cMono',mod=3,lig_ct='Th',lig='IL16')
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_IL16_trio2.pdf", useDingbats = FALSE,
    width = 5, height = 5.5)
lig_mod_fact 
dev.off()
# probably my most interesting/consistent hit so far...
# can check if it's replicated in Jimmie's dataset

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,
                                 mod_ct='cMono',mod=15,lig_ct='Th',lig='IL16')
lig_mod_fact

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,
                                 mod_ct='cMono',mod=6,lig_ct='Th',lig='IL16')
lig_mod_fact

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,
                                 mod_ct='cMono',mod=14,lig_ct='Th',lig='IL16')
lig_mod_fact

ctypes <- c('cMono')
# modules <- c(3)
modules <- c(14)
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.009, 
                                 db_use=c('GO'),
                                 max_plt_pval=.009,h_w=c(12,6))
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_IL16_enr.pdf", useDingbats = FALSE,
    width = 6, height = 7)
mod_enr
dev.off()

new_mods <- pbmc_container[["module_genes"]][["cMono"]]
new_mods[new_mods %in% c(2,3,4,6,9,15,14)] <- 3
new_mods[new_mods %in% c(3,6,15)] <- 3
old_mods <- pbmc_container[["module_genes"]][["cMono"]]
pbmc_container[["module_genes"]][["cMono"]] <- new_mods
pbmc_container[["module_genes"]][["cMono"]] <- old_mods

# getting p-values 
pbmc_container$lr_res['IL16_Th_CD4','cMono_m14']
which(rownames(myres_mat)=='IL16_Th_CD4')
10**(-fact_res2[299,2])

#### computing enrichment of modules with the nichenet scores
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

## ICOSLG target module first
test = ligand_target_matrix[,'IL16']
t_mod = pbmc_container[["module_genes"]][["cMono"]]
mymod <- c(2,3,4,6,9,15,14)
mymod <- c(3)
g_in_mod <- names(t_mod)[t_mod%in%mymod]
g_not_mod <- names(t_mod)[!(t_mod%in%mymod)]
tmp <- cbind.data.frame(c(g_in_mod,g_not_mod),
                        c(rep('cMono_m3',length(g_in_mod)),rep('other',length(g_not_mod))))
colnames(tmp) <- c('gn','in_mod')
tmp$in_mod <- factor(tmp$in_mod,levels=c('cMono_m3','other'))
target_scores <- test[tmp$gn]
tmp$target_scores <- target_scores
tmp <- tmp[which(!is.na(target_scores)),]
p <- ggplot(tmp,aes(x=as.factor(in_mod),y=target_scores)) +
  geom_boxplot(notch=TRUE) +
  ylab('PPBP NicheNet regulatory potential') +
  xlab('') +
  theme_bw()
test_res <- wilcox.test(target_scores~in_mod,data=tmp)
print(test_res$p.value)

## p-val for IL16 w cMono 3 is 0.000893



pdf(file = "/home/jmitchel/figures/for_paper_v2/ICOSLG_NicheNet_enr2.pdf", useDingbats = FALSE,
    width = 4, height = 3.5)
p
dev.off()



## plotting association between IL16 and HLA-DRA1
old_mods[old_mods==14]
old_mods[old_mods==3]
d_exp1 <- pbmc_container[["scale_pb_extra"]][['Th']][,'IL16']
d_exp2 <- pbmc_container[["no_scale_pb_extra"]][['cMono']][,'HLA-DPA1']
d_exp2 <- pbmc_container[["no_scale_pb_extra"]][['cMono']][,'SCARF1']
plot(d_exp1,d_exp2)
test <- colnames(pbmc_container[["no_scale_pb_extra"]][['cMono']])
test2 <- sapply(test,function(x){strsplit(x,split='HLA-DRA')[[1]]})
test2=unlist(test2)
which(test2=="")
for (i in test2) {
  if (length(test2[[i]]>1)) {
    print(test2[[i]])
  }
}








##### testing for subtype association shifts with F2
my_ctype <- 'Th'
subc <- pbmc@meta.data[pbmc@meta.data$initial_clustering==my_ctype,'full_clustering']
subc_lev <- unique(subc)
subc <- factor(subc,levels=subc_lev)
levels(subc) <- 1:length(subc_lev)
names(subc) <- rownames(pbmc@meta.data)[pbmc@meta.data$initial_clustering==my_ctype]
# convert subtype names to numbers
pbmc_container$subclusters[[my_ctype]][['res:0.5']] <- subc
pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype=my_ctype,res=.5,n_col=2)
pbmc_container$plots$ctype_prop_factor_associations

dotplot <- get_subclust_enr_dotplot(pbmc_container,ctype='Th',res=.5,subtype=1,factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_Th_prop.pdf", useDingbats = FALSE,
    width = 3.5, height = 4)
dotplot + ylim(0,.5)
dev.off()

## all Tc except naive change
dotplot <- get_subclust_enr_dotplot(pbmc_container,ctype='Tc',res=.5,subtype=1,factor_use=2)
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_Tc_prop.pdf", useDingbats = FALSE,
    width = 3.5, height = 4)
dotplot
dev.off()

## recording p-values of the prolif subtypes
dotplot <- get_subclust_enr_dotplot(pbmc_container,ctype='Th',res=.5,subtype=6,factor_use=2)
dotplot <- dotplot + scale_y_continuous(c(0,.10))
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_Th_prolif.pdf", useDingbats = FALSE,
    width = 3.5, height = 2.5)
dotplot # 3.73781501836407e-10
dev.off()

dotplot <- get_subclust_enr_dotplot(pbmc_container,ctype='Tc',res=.5,subtype=4,factor_use=2)
dotplot <- dotplot + scale_y_continuous(c(0,.15))
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_Tc_prolif.pdf", useDingbats = FALSE,
    width = 3.5, height = 2.5)
dotplot # 1.48574233120848e-13
dev.off()

dotplot <- get_subclust_enr_dotplot(pbmc_container,ctype='NK',res=.5,subtype=2,factor_use=2)
dotplot <- dotplot + scale_y_continuous(c(0,.3))
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_NK_prolif.pdf", useDingbats = FALSE,
    width = 3.5, height = 2.5)
dotplot # 1.44605467115252e-16
dev.off()

subc_lev










## running the full lr analysis for jimmie's covid data
pbmc_ye <- prep_LR_interact(pbmc_ye, lr_pairs, norm_method='trim', scale_factor=10000,
                            var_scale_power=.5, batch_var='pool')
sft_thresh <- c(12,14,10,14,12) #italk
pbmc_ye <- get_gene_modules(pbmc_ye,sft_thresh)

lr_hmap <- compute_LR_interact(pbmc_ye, lr_pairs, sig_thresh=.0001,
                               percentile_exp_rec=0.85, add_ld_fact_sig=TRUE)

pdf(file = "/home/jmitchel/figures/test3.pdf", useDingbats = FALSE,
    width = 20, height = 85)
lr_hmap
dev.off()

# make sure that dscores flipped to match those of Stephenson et al. dataset
pbmc_ye$tucker_results[[1]][,5] <- pbmc_ye$tucker_results[[1]][,5]*-1
lig_mod_fact <- plot_mod_and_lig(pbmc_ye,factor_select=5,
                                 mod_ct='cMono',mod=4,lig_ct='Th',lig='IL16')
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_jimmie_IL16_trio.pdf", useDingbats = FALSE,
    width = 5, height = 5.5)
lig_mod_fact
dev.off()

ctypes <- c('cMono')
modules <- c(4)
mod_enr <- plot_multi_module_enr(pbmc_ye, ctypes, modules, sig_thresh=.005, 
                                 db_use=c('GO'),
                                 max_plt_pval=.005,h_w=c(12,6))
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_jimmie_IL16_GO.pdf", useDingbats = FALSE,
    width = 6, height = 7)
mod_enr
dev.off()

tmp <- pbmc_ye[["module_genes"]][["cMono"]]
tmp[tmp==4]

sig_thresh=.0001
myres_mat <- pbmc_ye$lr_res # at 285
container=pbmc_ye
pbmc_ye$lr_res['IL16_Th_CD4','cMono_m4']
which(rownames(myres_mat)=='IL16_Th_CD4')
10**(-fact_res2[331,4])




# checking subcluster shifts in Jimmie's data
my_ctype <- 'NK'
subc <- pbmc@meta.data[pbmc@meta.data$ct1==my_ctype,'ct3']
subc_lev <- unique(subc)
subc <- factor(subc,levels=subc_lev)
levels(subc) <- 1:length(subc_lev)
names(subc) <- rownames(pbmc@meta.data)[pbmc@meta.data$ct1==my_ctype]
# convert subtype names to numbers
pbmc_ye$subclusters[[my_ctype]][['res:0.5']] <- subc
pbmc_ye <- get_ctype_subc_prop_associations(pbmc_ye,ctype=my_ctype,res=.5,n_col=2)
pbmc_ye$plots$ctype_prop_factor_associations

dotplot <- get_subclust_enr_dotplot(pbmc_ye,ctype='NK',res=.5,subtype=4,factor_use=5)
pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_Th_prop.pdf", useDingbats = FALSE,
    width = 3.5, height = 4)
dotplot + ylim(0,.5)
dev.off()

### seems like not the same trends but we don't really have the exact same populations
# Prolif subpops also sig: Tc prolif early 0.019, Th prolif early 0.045, NK prolif 0.036






# 
# f_test <- get_one_factor(pbmc_ye,5)
# dsc <- f_test[[1]]
# pbmc_ye <- get_donor_meta(pbmc_ye,additional_meta = c('covid_status','covid_severity_merged'),only_analyzed = F)
# 
# tmp <- cbind.data.frame(dsc,pbmc_ye$donor_metadata[rownames(dsc),'covid_severity_merged'],
#                         pbmc_ye$donor_metadata[rownames(dsc),'covid_status'])
# colnames(tmp) <- c('dscore','status_on_day_collection_summary','status')
# # order severity levels appropriately
# tmp$status_on_day_collection_summary <- factor(tmp$status_on_day_collection_summary,levels=c('Healthy','Moderate','Severe','Critical'))
# 
# mycol = brewer.pal(n = 8, name = "Dark2")
# 
# tmp2 <- tmp
# tmp2$status <- rep('violin',nrow(tmp2))
# p <- ggplot(tmp,aes(x=status_on_day_collection_summary,y=dscore,fill=status)) +
#   geom_violin(data=tmp2) +
#   geom_dotplot(binaxis='y', stackdir='center', dotsize=1.25, binwidth = .01) +
#   ylab('Factor 2 Donor Score') +
#   xlab('Severity on collection day') +
#   coord_flip() +
#   scale_fill_manual(values=c(mycol[5], mycol[2], 'light gray')) +
#   theme_bw() +
#   theme(axis.text=element_text(size=24),
#         axis.title=element_text(size=26))
# 
# # pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f2_severity2.pdf", useDingbats = FALSE,
# #     width = 12.5, height = 8.5)
# p
















### paused analysis just before line 858
# looking at genes that are IFN upregulated in factor 1 SLE













