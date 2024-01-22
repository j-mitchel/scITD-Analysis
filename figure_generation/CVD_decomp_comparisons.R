library(Seurat)
library(cowplot)
library(scITD)
library(cacoa)

# load up the lupus dataset: see preprocessing/lupus_preprocessing.R 
# for code used to generate this object
pbmc_sle <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')

# converting shorthand cell type
new_names <- sapply(as.character(pbmc_sle@meta.data$cg_cov), function(x){
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
pbmc_sle@meta.data$cg_cov <- factor(new_names,levels=unique(new_names))

# subset data to SLE patients only
cells_keep <- rownames(pbmc_sle@meta.data)[pbmc_sle@meta.data$Status=='Managed']
pbmc_sle <- subset(pbmc_sle,cells = cells_keep)

param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc","cDC",
                                               "cMono","ncMono"),
                                ncores = 30, rand_seed = 10)

pbmc_container_SLE <- make_new_container(seurat_obj=pbmc_sle,
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

# flip sign of F1 so high ISG expression is positive instead of negative (signs are arbitrary)
pbmc_container_SLE$tucker_results[[1]][,1] <- pbmc_container_SLE$tucker_results[[1]][,1] * -1
pbmc_container_SLE$tucker_results[[2]][1,] <- pbmc_container_SLE$tucker_results[[2]][1,] * -1
pbmc_container_SLE$projection_data[[1]][1,] <- pbmc_container_SLE$projection_data[[1]][1,] * -1

pbmc_container_SLE <- get_meta_associations(pbmc_container_SLE,vars_test=c('sex','Age','pool','processing','Ethnicity'),
                                            stat_use='pval')

## plot donor score
pbmc_container_SLE <- plot_donor_matrix(pbmc_container_SLE,
                                        show_donor_ids = FALSE,
                                        add_meta_associations='pval')

pbmc_container_SLE$plots$donor_matrix








# load up the covid dataset: see preprocessing/covid_uk_preprocessing.ipynb
# for code used to generate this object
pbmc_covid <- readRDS(file="/home/jmitchel/data/covid_data_uk/haniffa21_subset_no_lps.rds")

# collapsing cell subtypes and converting cell type names 
new_names <- sapply(as.character(pbmc_covid@meta.data$full_clustering), function(x){
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
pbmc_covid@meta.data$initial_clustering <- factor(new_names,levels=unique(new_names))

# remove plasmablast-like B cells
cells_keep1 <- rownames(pbmc_covid@reductions[["umap"]]@cell.embeddings)[pbmc_covid@reductions[["umap"]]@cell.embeddings[,1]>-1]
cells_keep2 <- rownames(pbmc_covid@reductions[["umap"]]@cell.embeddings)[pbmc_covid@reductions[["umap"]]@cell.embeddings[,2]>.5]
cells_keep <- unique(c(cells_keep1,cells_keep2))

pbmc_covid <- subset(pbmc_covid,cells=cells_keep)



# remove antibody "genes"
g_ndx_rem <- which(pbmc_covid@assays[["RNA"]]@meta.features=='Antibody Capture')
g_ndx_keep <- c(1:nrow(pbmc_covid@assays[["raw"]]@counts))[!(c(1:nrow(pbmc_covid@assays[["raw"]]@counts)) %in% g_ndx_rem)]
dim(pbmc_covid)
pbmc_covid <- CreateSeuratObject(pbmc_covid@assays[["raw"]]@counts[g_ndx_keep,], project = "SeuratProject", assay = "raw",
                                 meta.data = pbmc_covid@meta.data)
dim(pbmc_covid)


# starting analysis
param_list <- initialize_params(ctypes_use = c("B","Tc","Th","NK","cMono"), ncores = 30, rand_seed = 10)

pbmc_container_covid <- make_new_container(count_data = pbmc_covid@assays$raw@counts,
                                           meta_data = pbmc_covid@meta.data,
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

pbmc_container_covid <- form_tensor(pbmc_container_covid, donor_min_cells=2, 
                                    norm_method='trim', scale_factor=10000,
                                    vargenes_method='norm_var_pvals', vargenes_thresh=.0001,
                                    batch_var = 'site',
                                    scale_var = TRUE, var_scale_power = .5) 

pbmc_container_covid <- run_tucker_ica(pbmc_container_covid, ranks=c(9,38), 
                                       tucker_type = 'regular', rotation_type = 'hybrid') # best with vargenes_thresh=.0001 also .01

# get factor-meta data associations
pbmc_container_covid <- get_meta_associations(pbmc_container_covid,vars_test=c('sex','age',
                                                                               'status_on_day_collection_summary',
                                                                               'site'),stat_use='pval')
# plot donor scores by status
pbmc_container_covid <- plot_donor_matrix(pbmc_container_covid,
                                          show_donor_ids = FALSE,
                                          add_meta_associations='pval',h_w=c(10,8))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_dscores2.pdf", useDingbats = FALSE,
#     width = 8, height = 8)
pbmc_container_covid$plots$donor_matrix
# dev.off()


## flipping signs of f3 and f5
pbmc_container_covid$tucker_results[[1]][,3] <- pbmc_container_covid$tucker_results[[1]][,3] * -1
pbmc_container_covid$tucker_results[[2]][3,] <- pbmc_container_covid$tucker_results[[2]][3,] * -1
pbmc_container_covid$tucker_results[[1]][,5] <- pbmc_container_covid$tucker_results[[1]][,5] * -1
pbmc_container_covid$tucker_results[[2]][5,] <- pbmc_container_covid$tucker_results[[2]][5,] * -1









#### plotting the dimplot for only donors and cell types analyzed
d_kept <- rownames(pbmc_container_covid$scMinimal_ctype[[1]]$pseudobulk)
cells_keep <- rownames(pbmc_covid@meta.data)[pbmc_covid@meta.data$patient_id %in% d_kept]
pbmc_covid <- subset(pbmc_covid,cells=cells_keep)

ct_kept <- pbmc_container_covid$experiment_params$ctypes_use
cells_keep <- rownames(pbmc_covid@meta.data)[pbmc_covid@meta.data$initial_clustering %in% ct_kept]
pbmc_covid <- subset(pbmc_covid,cells=cells_keep)

# need to reload data because lost the umap when did CreateSeuratObject
# see preprocessing/covid_uk_preprocessing.ipynb
# for code used to generate this object
pbmc_tmp <- readRDS(file="/home/jmitchel/data/covid_data_uk/haniffa21_subset_no_lps.rds")
cells_keep <- rownames(pbmc_covid@meta.data)
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

### Figure 5a
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/covid_umap_major_clustering3.pdf", useDingbats = FALSE,
#     width = 15, height = 13)
DimPlot(pbmc_tmp, reduction = "umap", group.by = 'initial_clustering', label = TRUE) + ggtitle('Cell types analyzed') + NoLegend()
# dev.off()

###
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/covid_umap_subclusters3.pdf", useDingbats = FALSE,
#     width = 15, height = 13)
DimPlot(pbmc_tmp, reduction = "umap", group.by = 'full_clustering', label = TRUE) + ggtitle('Cell subtypes') + NoLegend()
# dev.off()



### now to project all SLE factors onto covid data
pbmc_container_covid <- project_new_data(pbmc_container_covid,pbmc_container_SLE)

res_orig <- cor(pbmc_container_covid[["projected_scores"]],pbmc_container_covid$tucker_results[[1]])
colnames(res_orig) <- sapply(c(1:ncol(res_orig)),function(x){
  paste0('Factor ',x)
})
rownames(res_orig) <- sapply(c(1:nrow(res_orig)),function(x){
  paste0('Projected ',x)
})

# use metadata as annotation for the hmap
col_fun_annot = colorRamp2(c(0, -log10(.1), 5), c("white", "white", "forest green"))
# removing site variable since nothing significant
pbmc_container_covid$meta_associations <- pbmc_container_covid$meta_associations[c('sex','age','status_on_day_collection_summary'),]
logpv <- -log10(pbmc_container_covid$meta_associations)
ba <- HeatmapAnnotation('-log_10_pval'=t(logpv),col = list('-log_10_pval'=col_fun_annot),
                        border=TRUE,annotation_name_side="right")

# make heatmap of this
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
cor_hmap <- Heatmap(res_orig, name = "Pearson r",
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    column_names_gp = gpar(fontsize = 10),
                    row_names_gp = gpar(fontsize = 10),
                    col = col_fun,border=TRUE, show_column_names=TRUE,
                    show_row_names=TRUE,show_row_dend = FALSE,
                    show_column_dend = FALSE,
                    column_title_side = "bottom",
                    row_names_side = "left",  
                    bottom_annotation = ba,
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid::grid.text(sprintf("%.2f", res_orig[i, j]), x, y, gp = gpar(fontsize = 10))
                    })

### Figure S9a
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/sle_covid_projection4.pdf", useDingbats = FALSE,
#     width = 8, height = 6.5)
cor_hmap
# dev.off()








# plotting F1 against status of three groups, healthy, covid, covid-critical
f_test <- get_one_factor(pbmc_container_covid,1)
dsc <- f_test[[1]]
pbmc_container_covid <- get_donor_meta(pbmc_container_covid,additional_meta = c('status_on_day_collection',
                                                                    'status_on_day_collection_summary',
                                                                    'worst_Clinical_Status',
                                                                    'swab_result',
                                                                    'status',
                                                                    'outcome',
                                                                    'site','age'),only_analyzed = FALSE)

tmp <- cbind.data.frame(dsc,pbmc_container_covid$donor_metadata[rownames(dsc),'status_on_day_collection_summary'],
                        pbmc_container_covid$donor_metadata[rownames(dsc),'status'])
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
mycol = brewer.pal(n = 8, name = "Dark2")

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

### Figure S9b
# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f1_status2.pdf", useDingbats = FALSE,
#     width = 12.5, height = 8.5)
p
# dev.off()



## compute enrichment p-values for covid-critical/non-critical patients in high/low dscores
meta_var='status_on_day_collection_summary'
factor_use=1
meta <- unique(pbmc_container_covid$scMinimal_full$metadata[,c('donors',meta_var)])
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
d_use <- intersect(rownames(pbmc_container_covid$tucker_results[[1]]),rownames(meta))
meta <- meta[d_use,]
mypaths <- list()
for (i in 1:length(meta_vals)) {
  mypaths[[meta_vals[i]]] <- rownames(meta)[meta[,meta_var]==meta_vals[i]]
}

myranks <- pbmc_container_covid$tucker_results[[1]][d_use,factor_use]

fgseaRes <- fgsea::fgsea(pathways = mypaths,
                         stats    = myranks,
                         minSize  = 0,
                         maxSize  = 5000)

print(fgseaRes)

# now regressing out the age association first
pbmc_container_covid <- get_donor_meta(pbmc_container_covid,additional_meta = c('status_on_day_collection_summary','status','age'))
head(pbmc_container_covid$donor_metadata)
meta <- pbmc_container_covid$donor_metadata
f <- 1
f_test <- get_one_factor(pbmc_container_covid,f)
dsc <- f_test[[1]]
tmp <- cbind.data.frame(dsc[rownames(meta),1],
                        meta$age)
colnames(tmp) <- c('dscore','age')
lmres <- lm(dscore~age,data=tmp)
dscore_resid <- lmres$residuals

fgseaRes2 <- fgsea::fgsea(pathways = mypaths,
                          stats    = dscore_resid,
                          minSize  = 0,
                          maxSize  = 5000)

print(fgseaRes2)



## now plotting factor 2 association with severity
f_test <- get_one_factor(pbmc_container_covid,2)
dsc <- f_test[[1]]
pbmc_container_covid <- get_donor_meta(pbmc_container_covid,additional_meta = c('status_on_day_collection',
                                                                    'status_on_day_collection_summary',
                                                                    'worst_Clinical_Status',
                                                                    'swab_result',
                                                                    'status',
                                                                    'outcome',
                                                                    'site','age'),only_analyzed = F)

tmp <- cbind.data.frame(dsc,pbmc_container_covid$donor_metadata[rownames(dsc),'status_on_day_collection_summary'],
                        pbmc_container_covid$donor_metadata[rownames(dsc),'status'])
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

### Figure 5b
# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f2_severity2.pdf", useDingbats = FALSE,
#     width = 12.5, height = 8.5)
p
# dev.off()





#### getting loadings plots
get_ct_specific_genes <- function(container,ct,factor_select,thresh=.05,signed=NULL) {
  ldngs <- get_one_factor(container,factor_select)[[2]]
  
  sig_vectors <- get_significance_vectors(container,
                                          factor_select, colnames(ldngs))
  # convert list to df
  sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
  
  # order df same way as in ldngs
  sig_df <- sig_df[rownames(ldngs),colnames(ldngs)]
  
  # loop through cell type columns and update mask vector for all genes
  gene_mask <- rep(TRUE,nrow(sig_df))
  for (i in 1:ncol(sig_df)) {
    if (colnames(sig_df)[i] %in% ct) {
      gene_mask <- gene_mask & sig_df[,i] < thresh
    } else {
      gene_mask <- gene_mask & !(sig_df[,i] < thresh)
    }
  }
  
  specific_genes <- names(gene_mask)[gene_mask]
  
  # now shorten this list to get strongest
  top_genes <- ldngs[specific_genes,ct,drop=FALSE]
  top_genes <- rowMeans(top_genes)
  
  if (is.null(signed)) {
    top_genes <- top_genes[order(abs(top_genes),decreasing = TRUE)]
  } else if (signed=='neg') {
    top_genes <- top_genes[order(top_genes,decreasing = FALSE)]
    top_genes <- top_genes[top_genes<0]
  } else if (signed=='pos') {
    top_genes <- top_genes[order(top_genes,decreasing = TRUE)]
    top_genes <- top_genes[top_genes>0]
  }
  
  return(top_genes)
}

# get significant genes
pbmc_container_covid <- get_lm_pvals(pbmc_container_covid)

t_all_spec_gns <- names(get_ct_specific_genes(pbmc_container_covid,c('Th'),1,thresh = .05,signed = 'neg'))
cm_spec_gns <- names(get_ct_specific_genes(pbmc_container_covid,c('cMono'),1,thresh = .05,signed = 'pos'))
cm_all_spec_gns <- names(get_ct_specific_genes(pbmc_container_covid,c('cMono','ncMono','cDC'),1,thresh = .01,signed = 'neg'))
th_spec_gns <- names(get_ct_specific_genes(pbmc_container_covid,c('Th'),1,thresh = .01,signed = 'neg'))
b_spec_gns <- names(get_ct_specific_genes(pbmc_container_covid,c('B'),1,thresh = .05,signed = 'pos'))

g_show <- c('IFI6','ISG15','MX1','USP18',
            'JUP','LILRB2','FOXP3','CD70')

### Figure S9c
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/covid_f1_lds.pdf", useDingbats = FALSE,
#     width = 5, height = 5)
pbmc_container_covid <- plot_loadings_annot(pbmc_container_covid, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                            pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                            gene_callouts=TRUE, specific_callouts=g_show,
                                            show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                            clust_method='complete', h_w=c(5,6.5))
# dev.off()





## making the factor 2 loadings plot
pbmc_container_covid <- run_gsea_one_factor(pbmc_container_covid, factor_select=2, method="fgsea", thresh=0.05,
                                            db_use=c("GO"))

## making the factor 2 loadings plot
gsets <- c("GOBP_CELL_CYCLE",
           "GOBP_CELL_KILLING")
gset_cmap <- c('red',
               'black')
names(gset_cmap) <- gsets

gset_cmap_sub <- gset_cmap[gset_cmap!='black']
gset_sub <- names(gset_cmap_sub)

pbmc_container_covid <- plot_loadings_annot(pbmc_container_covid, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                            pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.01, display_genes=FALSE,
                                            gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                            le_set_callouts=gset_sub, le_set_colormap=gset_cmap_sub, le_set_num_per=10, show_le_legend=FALSE,
                                            show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                            clust_method='mcquitty', h_w=c(9,6.5))

# extracting some proliferation markers that showed up in the loadings
Idents(pbmc_covid) <- pbmc_covid@meta.data$full_clustering
prolif.de.markers1 <- FindMarkers(pbmc_covid, ident.1 = c("CD8.Prolif"), ident.2 = c('CD8.Naive','CD8.EM','CD8.TE'), min.pct = 0.5, logfc.threshold = log(2))
prolif.de.markers2 <- FindMarkers(pbmc_covid, ident.1 = c("CD4.Prolif"), ident.2 = c('CD4.Naive','CD4.IL22','CD4.CM','CD4.EM','CD4.Th1','CD4.Tfh'),min.pct = 0.5,logfc.threshold = log(2))
prolif.de.markers3 <- FindMarkers(pbmc_covid, ident.1 = c("NK_prolif"), ident.2 = c('NK_16hi','NK_56hi'), min.pct = 0.5,logfc.threshold = log(2))

p_g1 <- intersect(rownames(prolif.de.markers1)[1:100],rownames(prolif.de.markers2)[1:100])
p_g1 <- intersect(p_g1,rownames(prolif.de.markers3)[1:100])
print(p_g1)

t_all_spec_gns <- names(get_ct_specific_genes(pbmc_container_covid,c('Th','Tc','NK'),2,thresh = .05,signed = 'pos'))
t_all_spec_gns
t_all_spec_gns[t_all_spec_gns%in%rownames(prolif.de.markers1)]
t_all_spec_gns[t_all_spec_gns%in%rownames(prolif.de.markers2)]
t_all_spec_gns[t_all_spec_gns%in%rownames(prolif.de.markers3)]

g_show <- c('RAD51','UBE2C','CENPF','CDK1','HLA-DRB1','HLA-DPB1','HLA-DRA')

### Figure 5c
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/covid_f2_lds.pdf", useDingbats = FALSE,
#     width = 5, height = 5)
pbmc_container_covid <- plot_loadings_annot(pbmc_container_covid, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                            pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                            gene_callouts=TRUE, specific_callouts=g_show,
                                            show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                            clust_method='complete', h_w=c(5,6.5))
# dev.off()






##### testing for subtype association shifts with F2
my_ctype <- 'Th'
subc <- pbmc_covid@meta.data[pbmc_covid@meta.data$initial_clustering==my_ctype,'full_clustering']
subc_lev <- unique(subc)
subc <- factor(subc,levels=subc_lev)
levels(subc) <- 1:length(subc_lev)
names(subc) <- rownames(pbmc_covid@meta.data)[pbmc_covid@meta.data$initial_clustering==my_ctype]
# convert subtype names to numbers
pbmc_container_covid$subclusters[[my_ctype]][['res:0.5']] <- subc
pbmc_container_covid <- get_ctype_subc_prop_associations(pbmc_container_covid,ctype=my_ctype,res=.5,n_col=2)
pbmc_container_covid$plots$ctype_prop_factor_associations
dotplot <- get_subclust_enr_dotplot(pbmc_container_covid,ctype='Th',res=.5,subtype=1,factor_use=2)

### Figure S8a left
# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_Th_prop.pdf", useDingbats = FALSE,
#     width = 3.5, height = 4)
dotplot + ylim(0,.5)
# dev.off()

my_ctype <- 'Tc'
subc <- pbmc_covid@meta.data[pbmc_covid@meta.data$initial_clustering==my_ctype,'full_clustering']
subc_lev <- unique(subc)
subc <- factor(subc,levels=subc_lev)
levels(subc) <- 1:length(subc_lev)
names(subc) <- rownames(pbmc_covid@meta.data)[pbmc_covid@meta.data$initial_clustering==my_ctype]
# convert subtype names to numbers
pbmc_container_covid$subclusters[[my_ctype]][['res:0.5']] <- subc
pbmc_container_covid <- get_ctype_subc_prop_associations(pbmc_container_covid,ctype=my_ctype,res=.5,n_col=2)
pbmc_container_covid$plots$ctype_prop_factor_associations
dotplot <- get_subclust_enr_dotplot(pbmc_container_covid,ctype='Tc',res=.5,subtype=1,factor_use=2)

### Figure S8a right
# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_Tc_prop.pdf", useDingbats = FALSE,
#     width = 3.5, height = 4)
dotplot
# dev.off()

## recording p-values of the prolif subtypes
dotplot <- get_subclust_enr_dotplot(pbmc_container_covid,ctype='Th',res=.5,subtype=6,factor_use=2)
dotplot <- dotplot + scale_y_continuous(c(0,.10))

### Figure 5e top right
# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_Th_prolif.pdf", useDingbats = FALSE,
#     width = 3.5, height = 2.5)
dotplot # 3.73781501836407e-10
# dev.off()

dotplot <- get_subclust_enr_dotplot(pbmc_container_covid,ctype='Tc',res=.5,subtype=4,factor_use=2)
dotplot <- dotplot + scale_y_continuous(c(0,.15))

### Figure 5e bottom
# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_Tc_prolif.pdf", useDingbats = FALSE,
#     width = 3.5, height = 2.5)
dotplot # 1.48574233120848e-13
# dev.off()


# need to add the subcluster data for NK
my_ctype <- 'NK'
subc <- pbmc_covid@meta.data[pbmc_covid@meta.data$initial_clustering==my_ctype,'full_clustering']
subc_lev <- unique(subc)
subc <- factor(subc,levels=subc_lev)
levels(subc) <- 1:length(subc_lev)
names(subc) <- rownames(pbmc_covid@meta.data)[pbmc_covid@meta.data$initial_clustering==my_ctype]
# convert subtype names to numbers
pbmc_container_covid$subclusters[[my_ctype]][['res:0.5']] <- subc
dotplot <- get_subclust_enr_dotplot(pbmc_container_covid,ctype='NK',res=.5,subtype=2,factor_use=2)
dotplot <- dotplot + scale_y_continuous(c(0,.3))

### Figure 5e top left
# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_NK_prolif.pdf", useDingbats = FALSE,
#     width = 3.5, height = 2.5)
dotplot # 1.44605467115252e-16
# dev.off()












####### new factor comparison approach to compare covid vs sle cognate factors
get_betas <- function(container, all_f_test, donor_min_cells=2,
                      norm_method='trim', scale_factor=10000,
                      batch_var = 'Batch', get_pvals=FALSE) {
  # Need to run a modified script to get pseudobulk expression without scaling, so can remove lowly expressed genes
  container <- parse_data_by_ctypes(container)
  container <- clean_data(container, donor_min_cells=donor_min_cells)
  container <- get_pseudobulk(container)
  container <- normalize_pseudobulk(container, method=norm_method, scale_factor=scale_factor)
  
  # remove lowly expressed genes here
  for (ct in container$experiment_params$ctypes_use) {
    # select only genes expressed to some amount in at least 5% of donors
    donor_thresh <- round(nrow(container$scMinimal_ctype[[ct]]$pseudobulk) * .05)
    g_keep <- colSums(container$scMinimal_ctype[[ct]]$pseudobulk>0) > donor_thresh
    container$scMinimal_ctype[[ct]]$pseudobulk <- container$scMinimal_ctype[[ct]]$pseudobulk[,g_keep]
  }
  
  container <- apply_combat(container,batch_var=batch_var)
  
  # unit scale each gene in each cell type
  for (ct in container$experiment_params$ctypes_use) {
    container$scMinimal_ctype[[ct]]$pseudobulk <- scale(container$scMinimal_ctype[[ct]]$pseudobulk)
  }
  
  # need to unit scale donor scores for each factor, so can compare betas
  container$tucker_results[[1]] <- scale(container$tucker_results[[1]])
  
  ## now compute expression-dsc associations for all factors and all cell types
  total_res <- list() # stores the results for all factors
  total_res_pv <- list()
  for (f_test in all_f_test) {
    print(f_test)
    dsc <- container$tucker_results[[1]][,f_test,drop=FALSE]
    
    f_res <- list() # stores the result from a single factor
    f_res_pv <- list() 
    
    # loop through cell types
    for (ct in container$experiment_params$ctypes_use) {
      print(ct)
      ct_res <- c()
      ct_res_pv <- c()
      
      pb <- container$scMinimal_ctype[[ct]]$pseudobulk
      # loop through genes
      for (g_ndx in 1:ncol(pb)) {
        tmp <- cbind.data.frame(dsc,pb[rownames(dsc),g_ndx])
        colnames(tmp) <- c('dscore','expr')
        lmres <- summary(lm(expr~dscore,data=tmp))
        beta <- lmres$coefficients['dscore','Estimate']
        ct_res <- c(ct_res,beta)
        
        pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
        ct_res_pv <- c(ct_res_pv,pval)
      }
      names(ct_res) <- colnames(pb)
      f_res[[ct]] <- ct_res
      f_res_pv[[ct]] <- p.adjust(ct_res_pv,method='fdr')
    }
    total_res[[f_test]] <- f_res
    total_res_pv[[f_test]] <- f_res_pv
  }
  
  if (get_pvals) {
    return(list(total_res,total_res_pv))
  } else {
    return(total_res)
  }
}

betas_covid <- get_betas(pbmc_container_covid, all_f_test=c(1), donor_min_cells=2,
                       norm_method='trim', scale_factor=10000,
                       batch_var = 'site')

betas_sle <- get_betas(pbmc_container_SLE, all_f_test=c(1), donor_min_cells=20,
                         norm_method='trim', scale_factor=10000,
                         batch_var='pool')


# function to get differences in betas, run gsea on these, and plot the result
get_diff_gsea <- function(betas1,betas2,cao,f_use1,f_use2,plot=TRUE,plot_n=25) {
  # limit to intersection of genes we have betas for in both datasets
  for (ct in names(betas1[[f_use1]])) {
    genes1 <- names(betas1[[f_use1]][[ct]])
    genes2 <- names(betas2[[f_use2]][[ct]])
    gene_intersect <- intersect(genes1,genes2)
    betas1[[f_use1]][[ct]] <- betas1[[f_use1]][[ct]][gene_intersect]
    betas2[[f_use2]][[ct]] <- betas2[[f_use2]][[ct]][gene_intersect]
  }
  
  # for each cell type take the difference across datasets
  betas_diff <- list()
  for (ct in names(betas1[[f_use1]])) {
    betas_diff[[ct]] <- betas1[[f_use1]][[ct]] - betas2[[f_use2]][[ct]]
  }
  
  all_res <- lapply(1:length(betas_diff),function(x){
    ct <- names(betas_diff)[x]
    diff_vals <- betas_diff[[ct]]
    res <- as.data.frame(matrix(ncol=3,nrow=length(diff_vals)))
    colnames(res) <- c('padj','Z','Gene')
    res$Gene <- names(diff_vals)
    
    # res$Z <- diff_vals
    # res$Z <- abs(diff_vals) # TESTING ONLY
    res$Z <- abs(diff_vals**2) # TESTING ONLY
    # res$Z <- abs(diff_vals**2.5) # TESTING ONLY
    
    res$padj <- rep(1,length(diff_vals))
    # sort by absolute value score
    res <- res[order(abs(res$Z),decreasing = TRUE),]
    res_outer <- list(res)
    names(res_outer) <- 'res'
    return(res_outer)
  })
  names(all_res) <- names(betas_diff)
  
  cao[["test.results"]][["de"]] <- all_res
  
  # cao$estimateOntology(type="GSEA", org.db=org.Hs.eg.db::org.Hs.eg.db, verbose=TRUE,
  #                      n.cores=5, ignore.cache=TRUE)
  cao$estimateOntology(type="GSEA", org.db=org.Hs.eg.db::org.Hs.eg.db, verbose=TRUE,
                       n.cores=5, ignore.cache=TRUE,
                       scoreType = "pos")
  
  # cao$estimateOntology(type="GO", org.db=org.Hs.eg.db::org.Hs.eg.db, verbose=TRUE,
  #                      n.cores=5, ignore.cache=TRUE)
  
  if (!plot) {
    return(cao)
  } else {
    ex_words <- c('regulation', 'process', 'cell', 'positive', 'negative')
    p <- cao$plotOntologyHeatmapCollapsed(
      name="GSEA", genes="up", n=plot_n, clust.method="ward.D", size.range=c(1, 4),
      exclude.words=ex_words, p.adj=.05, q.value=.2
    )
    
    return(p)
  }
  
}


# need some empty cacoa object, so this is arbitrary data that's not use except for its data structure
cao <- cacoa::Cacoa$new(
  NULL, sample.groups=as.factor(c('ctrl','test')), 
  cell.groups=as.factor(c('ct1','ct2')),
  sample.per.cell=as.factor(c('1','2')),
  target.level='test', ref.level='ctrl', n.cores=5, verbose=FALSE
)


# this analysis needs to be done per factor
p1 <- get_diff_gsea(betas_covid,betas_sle,cao,f_use1=1,f_use2=1,plot_n=10)

### Figure S9e
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/gsea_compare_covid_sle3.pdf", useDingbats = FALSE,
#     width = 6, height = 5)
p1
# dev.off()


p1 <- get_diff_gsea(betas_covid,betas_sle,cao,f_use1=1,f_use2=1,plot=F)
ct <- 'Tc'
tmp=p1[["test.results"]][["GSEA"]][["res"]][[ct]][["BP"]]@result
rownames(tmp) <- tmp$Description
print(tmp[1:60,c('Description','p.adjust')])


### making comparative enrichment plots of genes in "generation of precursor metabolites and energy" gene set
m_df <- msigdbr::msigdbr(species = "Homo sapiens",
                                    category = "C5", subcategory = "BP")
my_pathways <- split(m_df$gene_symbol, f = m_df$gs_name)

mypaths <- my_pathways['GOBP_ATP_METABOLIC_PROCESS']
myranks <- betas_covid[[1]][['Tc']]
fgseaRes_cvd <- fgsea::fgsea(pathways = mypaths,
                         stats    = myranks,
                         eps = 0,
                         minSize  = 0,
                         maxSize  = 5000)
print(fgseaRes_cvd)

plt_covid <- fgsea::plotEnrichment(my_pathways[['GOBP_ATP_METABOLIC_PROCESS']],
                                   myranks) + labs(title='ATP metabolic process genes - CVD_F1')
plt_covid <- plt_covid + theme_bw()

myranks <- betas_sle[[1]][['Tc']]
fgseaRes_sle <- fgsea::fgsea(pathways = mypaths,
                         stats    = myranks,
                         eps = 0,
                         minSize  = 0,
                         maxSize  = 5000)
print(fgseaRes_sle)

plt_sle <- fgsea::plotEnrichment(my_pathways[['GOBP_ATP_METABOLIC_PROCESS']],
                                 myranks) + labs(title='ATP metabolic process genes - aSLE_F1')
plt_sle <- plt_sle + theme_bw()

fig <- cowplot::plot_grid(plt_sle,plt_covid,nrow=2)

### Figure S9f
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/ATP_process_fgsea.pdf", useDingbats = FALSE,
#     width = 5, height = 8.5)
fig
# dev.off()















## plotting correlations across genes for gene-factor associations
betas_covid <- get_betas(pbmc_container_covid, all_f_test=c(1,3,5,6,7), donor_min_cells=2,
                         norm_method='trim', scale_factor=10000,
                         batch_var = 'site', get_pvals = TRUE)

betas_sle <- get_betas(pbmc_container_SLE, all_f_test=c(1,2,3,4,5), donor_min_cells=20,
                       norm_method='trim', scale_factor=10000,
                       batch_var='pool', get_pvals = TRUE)


beta_vals_cvd <- betas_covid[[1]]
beta_vals_sle <- betas_sle[[1]]
p_vals_cvd <- betas_covid[[2]]
p_vals_sle <- betas_sle[[2]]

ctypes_all <- c('B','NK','Th','Tc','cMono')
cor_res <- matrix(nrow=length(beta_vals_sle),ncol=length(ctypes_all))
rownames(cor_res) <- c('aSLE_F1:CVD_F1','aSLE_F2:CVD_F3','aSLE_F5:CVD_F5','aSLE_F3:CVD_F6','aSLE_F4:CVD_F7')
colnames(cor_res) <- ctypes_all
covid_f_map <- c(1,3,5,6,7)
sle_f_map <- c(1,2,5,3,4)
for (i in 1:length(beta_vals_sle)) {
  betas_sle_f <- beta_vals_sle[[sle_f_map[i]]]
  betas_cvd_f <- beta_vals_cvd[[covid_f_map[i]]]
  
  p_sle_f <- p_vals_sle[[sle_f_map[i]]]
  p_cvd_f <- p_vals_cvd[[covid_f_map[i]]]
  
  # limit all to the cell types found in both
  betas_sle_f <- betas_sle_f[ctypes_all]
  betas_cvd_f <- betas_cvd_f[ctypes_all]
  p_sle_f <- p_sle_f[ctypes_all]
  p_cvd_f <- p_cvd_f[ctypes_all]
  
  
  # loop through cell types
  for (ct in names(betas_sle_f)) {
    betas_sle_f_ct <- betas_sle_f[[ct]]
    betas_cvd_f_ct <- betas_cvd_f[[ct]]
    p_sle_f_ct <- p_sle_f[[ct]]
    p_cvd_f_ct <- p_cvd_f[[ct]]
    
    names(p_sle_f_ct) <- names(betas_sle_f_ct)
    names(p_cvd_f_ct) <- names(betas_cvd_f_ct)
    
    # reduce genes to the intersection of genes tested in both
    g_both <- intersect(names(betas_sle_f_ct),names(betas_cvd_f_ct))
    betas_sle_f_ct <- betas_sle_f_ct[g_both]
    betas_cvd_f_ct <- betas_cvd_f_ct[g_both]
    p_sle_f_ct <- p_sle_f_ct[g_both]
    p_cvd_f_ct <- p_cvd_f_ct[g_both]
    
    # get union of significant genes
    g_keep1 <- names(p_sle_f_ct)[p_sle_f_ct<.01]
    g_keep2 <- names(p_cvd_f_ct)[p_cvd_f_ct<.01]
    g_keep <- unique(c(g_keep1,g_keep2))
    # print(i)
    # print(ct)
    # print(length(g_keep))
    
    if (length(g_keep)<10) {
      cor_res[i,ct] <- NA
    } else {
      betas_sle_f_ct <- betas_sle_f_ct[g_keep]
      betas_cvd_f_ct <- betas_cvd_f_ct[g_keep]
      
      mycor <- cor(betas_sle_f_ct,betas_cvd_f_ct)
      cor_res[i,ct] <- mycor
    }
  }
}

col_fun <- colorRamp2(c(0, 1), c("white", "red"))
h_w <- c(5,7)
myhmap <- Heatmap(as.matrix(cor_res), name = 'Pearson cor',
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
                  border=TRUE,
                  row_names_side = "left",
                  na_col = "light gray",
                  row_names_gp = gpar(fontsize = 12),
                  show_heatmap_legend = TRUE,
                  width = unit(h_w[2], "cm"),
                  height = unit(h_w[1], "cm"),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid::grid.text(sprintf("%.2f", as.matrix(cor_res)[i, j]), x, y, gp = gpar(fontsize = 10))
                  })

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/factor_assoc_cors_sle_cvd_v2.pdf", useDingbats = FALSE,
#     width = 6, height = 5)
myhmap
# dev.off()













####### trying subsampling to match age/sex characteristics of the populations

pbmc_container_covid <- get_donor_meta(pbmc_container_covid,c('age','sex'))
pbmc_container_SLE <- get_donor_meta(pbmc_container_SLE,c('Age','sex'))


age_levels <- levels(pbmc_container_covid$donor_metadata$age)
age_ranges <- list(c(20,29),c(30,39),c(40,49),c(50,59),c(60,69),c(70,79),c(80,89),c(90,99))
sle_age_fact <- sapply(pbmc_container_SLE$donor_metadata$Age,function(x){
  for (i in 1:length(age_ranges)) {
    if (x >= age_ranges[[i]][1] && x <= age_ranges[[i]][2]) {
      return(age_levels[i])
    }
  }
})
pbmc_container_SLE$donor_metadata$age_levels <- sle_age_fact
head(pbmc_container_SLE$donor_metadata)

age_bins_sle <- table(pbmc_container_SLE$donor_metadata$age_levels) / nrow(pbmc_container_SLE$donor_metadata)
age_bins_cvd <- table(pbmc_container_covid$donor_metadata$age) / nrow(pbmc_container_covid$donor_metadata)

# adding a 90-99 bin range for sle since no one in that bin
age_bins_sle <- c(age_bins_sle,0)

c1 <- c(age_bins_sle,age_bins_cvd)
c2 <- c(rep('sle',length(age_bins_sle)),rep('covid',length(age_bins_cvd)))
c3 <- rep(names(age_bins_cvd),2)
tmp <- cbind.data.frame(c1,c2,c3)
colnames(tmp) <- c('age_frac','disease','age_bin')

tmp1 <- tmp[tmp$disease=='sle',]
tmp2 <- tmp[tmp$disease=='covid',]

p1 <- ggplot(tmp1,aes(x=as.factor(age_bin),y=age_frac,fill=disease)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#E69F00")) +
  theme_bw()

p2 <- ggplot(tmp2,aes(x=as.factor(age_bin),y=age_frac,fill=disease)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#56B4E9")) +
  theme_bw()

f1 <- plot_grid(p1,p2,ncol=1)


### also showing sample sex distribution for the two datasets
sx_bins_sle <- table(pbmc_container_SLE$donor_metadata$sex) / nrow(pbmc_container_SLE$donor_metadata)
sx_bins_cvd <- table(pbmc_container_covid$donor_metadata$sex) / nrow(pbmc_container_covid$donor_metadata)

c1 <- c(sx_bins_sle,sx_bins_cvd)
c2 <- c(rep('sle',length(sx_bins_sle)),rep('covid',length(sx_bins_cvd)))
c3 <- rep(names(sx_bins_cvd),2)
tmp <- cbind.data.frame(c1,c2,c3)
colnames(tmp) <- c('sx_frac','disease','sx_bin')

tmp1 <- tmp[tmp$disease=='sle',]
tmp2 <- tmp[tmp$disease=='covid',]

p1 <- ggplot(tmp1,aes(x=as.factor(sx_bin),y=sx_frac,fill=disease)) +
  geom_bar(position="dodge", stat="identity", width=.75) +
  scale_fill_manual(values=c("#E69F00")) +
  ylim(0,1) +
  theme_bw()

p2 <- ggplot(tmp2,aes(x=as.factor(sx_bin),y=sx_frac,fill=disease)) +
  geom_bar(position="dodge", stat="identity", width=.75) +
  scale_fill_manual(values=c("#56B4E9")) +
  ylim(0,1) +
  theme_bw()

f2 <- plot_grid(p1,p2,ncol=1)

f3 <- plot_grid(f1,f2,ncol=2,rel_widths = c(.5,.25))
f3

### Figure S10a
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/covid_sle_sex_age_unbalanced.pdf", useDingbats = FALSE,
#     width = 9, height = 6.5)
f3
# dev.off()













### trying to downsample sex first, then age after
# now downsampling males from covid dataset
pbmc_container_covid$donor_metadata_sub <- pbmc_container_covid$donor_metadata
group1 <- pbmc_container_covid$donor_metadata_sub$donors[pbmc_container_covid$donor_metadata_sub$sex=='Male']
group_other <- pbmc_container_covid$donor_metadata_sub$donors[!(pbmc_container_covid$donor_metadata_sub$sex %in% c('Male'))]

group1 <- sample(group1,round(length(group1)*.1))

d_final <- c(group1,group_other)

pbmc_container_covid$donor_metadata_sub <- pbmc_container_covid$donor_metadata_sub[pbmc_container_covid$donor_metadata_sub$donors%in%d_final,]

age_bins_sle <- table(pbmc_container_SLE$donor_metadata$age_levels) / nrow(pbmc_container_SLE$donor_metadata)
age_bins_cvd <- table(pbmc_container_covid$donor_metadata_sub$age) / nrow(pbmc_container_covid$donor_metadata_sub)

# adding a 90-99 bin range for sle since no one in that bin
age_bins_sle <- c(age_bins_sle,0)

c1 <- c(age_bins_sle,age_bins_cvd)
c2 <- c(rep('sle',length(age_bins_sle)),rep('covid',length(age_bins_cvd)))
c3 <- rep(names(age_bins_cvd),2)
tmp <- cbind.data.frame(c1,c2,c3)
colnames(tmp) <- c('age_frac','disease','age_bin')

tmp1 <- tmp[tmp$disease=='sle',]
tmp2 <- tmp[tmp$disease=='covid',]

p1 <- ggplot(tmp1,aes(x=as.factor(age_bin),y=age_frac,fill=disease)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#E69F00")) +
  theme_bw()

p2 <- ggplot(tmp2,aes(x=as.factor(age_bin),y=age_frac,fill=disease)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#56B4E9")) +
  theme_bw()

f1 <- plot_grid(p1,p2,ncol=1)


### also showing sample sex distribution for the two datasets
sx_bins_sle <- table(pbmc_container_SLE$donor_metadata$sex) / nrow(pbmc_container_SLE$donor_metadata)
sx_bins_cvd <- table(pbmc_container_covid$donor_metadata_sub$sex) / nrow(pbmc_container_covid$donor_metadata_sub)

c1 <- c(sx_bins_sle,sx_bins_cvd)
c2 <- c(rep('sle',length(sx_bins_sle)),rep('covid',length(sx_bins_cvd)))
c3 <- rep(names(sx_bins_cvd),2)
tmp <- cbind.data.frame(c1,c2,c3)
colnames(tmp) <- c('sx_frac','disease','sx_bin')

tmp1 <- tmp[tmp$disease=='sle',]
tmp2 <- tmp[tmp$disease=='covid',]

p1 <- ggplot(tmp1,aes(x=as.factor(sx_bin),y=sx_frac,fill=disease)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#E69F00")) +
  ylim(0,1) +
  theme_bw()

p2 <- ggplot(tmp2,aes(x=as.factor(sx_bin),y=sx_frac,fill=disease)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#56B4E9")) +
  ylim(0,1) +
  theme_bw()

f2 <- plot_grid(p1,p2,ncol=1)

f4 <- plot_grid(f1,f2,ncol=2)
f4


## now downsampling by age
# sle dataset first
pbmc_container_SLE$donor_metadata_sub <- pbmc_container_SLE$donor_metadata
group1 <- pbmc_container_SLE$donor_metadata$donors[pbmc_container_SLE$donor_metadata$age_levels=='(20, 29]']
group2 <- pbmc_container_SLE$donor_metadata$donors[pbmc_container_SLE$donor_metadata$age_levels=='(30, 39]']
group3 <- pbmc_container_SLE$donor_metadata$donors[pbmc_container_SLE$donor_metadata$age_levels=='(40, 49]']
group_other <- pbmc_container_SLE$donor_metadata$donors[!(pbmc_container_SLE$donor_metadata$age_levels %in% c('(20, 29]','(30, 39]','(40, 49]'))]

group1 <- sample(group1,round(length(group1)/2))
group2 <- sample(group2,round(length(group2)/3.5))
group3 <- sample(group3,round(length(group3)*3/5))

d_final <- c(group1,group2,group3,group_other)

pbmc_container_SLE$donor_metadata_sub <- pbmc_container_SLE$donor_metadata_sub[pbmc_container_SLE$donor_metadata_sub$donors%in%d_final,]

# now covid data
pbmc_container_covid$donor_metadata_sub_sub <- pbmc_container_covid$donor_metadata_sub
group2 <- pbmc_container_covid$donor_metadata_sub$donors[pbmc_container_covid$donor_metadata_sub$age=='(70, 79]']
group3 <- pbmc_container_covid$donor_metadata_sub$donors[pbmc_container_covid$donor_metadata_sub$age=='(80, 89]']
group_other <- pbmc_container_covid$donor_metadata_sub$donors[!(pbmc_container_covid$donor_metadata_sub$age %in% c('(70, 79]','(80, 89]'))]

group2 <- sample(group2,round(length(group2)/4))
group3 <- sample(group3,round(length(group3)/5))

d_final <- c(group2,group3,group_other)

pbmc_container_covid$donor_metadata_sub_sub <- pbmc_container_covid$donor_metadata_sub_sub[pbmc_container_covid$donor_metadata_sub_sub$donors%in%d_final,]


age_bins_sle <- table(pbmc_container_SLE$donor_metadata_sub$age_levels) / nrow(pbmc_container_SLE$donor_metadata_sub)
age_bins_cvd <- table(pbmc_container_covid$donor_metadata_sub_sub$age) / nrow(pbmc_container_covid$donor_metadata_sub_sub)

# adding a 90-99 bin range for sle since no one in that bin
age_bins_sle <- c(age_bins_sle,0)

c1 <- c(age_bins_sle,age_bins_cvd)
c2 <- c(rep('sle',length(age_bins_sle)),rep('covid',length(age_bins_cvd)))
c3 <- rep(names(age_bins_cvd),2)
tmp <- cbind.data.frame(c1,c2,c3)
colnames(tmp) <- c('age_frac','disease','age_bin')

tmp1 <- tmp[tmp$disease=='sle',]
tmp2 <- tmp[tmp$disease=='covid',]

p1 <- ggplot(tmp1,aes(x=as.factor(age_bin),y=age_frac,fill=disease)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#E69F00")) +
  theme_bw()

p2 <- ggplot(tmp2,aes(x=as.factor(age_bin),y=age_frac,fill=disease)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("#56B4E9")) +
  theme_bw()

f1 <- plot_grid(p1,p2,ncol=1)


### also showing sample sex distribution for the two datasets
sx_bins_sle <- table(pbmc_container_SLE$donor_metadata_sub$sex) / nrow(pbmc_container_SLE$donor_metadata_sub)
sx_bins_cvd <- table(pbmc_container_covid$donor_metadata_sub_sub$sex) / nrow(pbmc_container_covid$donor_metadata_sub_sub)

c1 <- c(sx_bins_sle,sx_bins_cvd)
c2 <- c(rep('sle',length(sx_bins_sle)),rep('covid',length(sx_bins_cvd)))
c3 <- rep(names(sx_bins_cvd),2)
tmp <- cbind.data.frame(c1,c2,c3)
colnames(tmp) <- c('sx_frac','disease','sx_bin')

tmp1 <- tmp[tmp$disease=='sle',]
tmp2 <- tmp[tmp$disease=='covid',]

p1 <- ggplot(tmp1,aes(x=as.factor(sx_bin),y=sx_frac,fill=disease)) +
  geom_bar(position="dodge", stat="identity", width=.75) +
  scale_fill_manual(values=c("#E69F00")) +
  ylim(0,1) +
  theme_bw()

p2 <- ggplot(tmp2,aes(x=as.factor(sx_bin),y=sx_frac,fill=disease)) +
  geom_bar(position="dodge", stat="identity", width=.75) +
  scale_fill_manual(values=c("#56B4E9")) +
  ylim(0,1) +
  theme_bw()

f2 <- plot_grid(p1,p2,ncol=1)

f4 <- plot_grid(f1,f2,ncol=2,rel_widths = c(.5,.25))
f4

### Figure S10b
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/covid_sle_sex_age_balanced.pdf", useDingbats = FALSE,
#     width = 9, height = 6.5)
f4
# dev.off()

### saving donor subsample lists so I don't need to recompute these
cvd_dnr_keep <- pbmc_container_covid$donor_metadata_sub_sub$donors
sle_dnr_keep <- pbmc_container_SLE$donor_metadata_sub$donors

all_dnr_keep <- list(cvd_dnr_keep,sle_dnr_keep)
# saveRDS(all_dnr_keep,file='/home/jmitchel/data/lupus_data/donors_age_sex_match.rds')

all_dnr_keep <- readRDS(file='/home/jmitchel/data/lupus_data/donors_age_sex_match.rds')
print(length(unique(all_dnr_keep[[1]])))
print(length(unique(all_dnr_keep[[2]])))










# now get betas again from the original datasets
betas_covid_full <- get_betas(pbmc_container_covid, all_f_test=c(1), donor_min_cells=2,
                         norm_method='trim', scale_factor=10000,
                         batch_var = 'site')

betas_sle_full <- get_betas(pbmc_container_SLE, all_f_test=c(1), donor_min_cells=20,
                       norm_method='trim', scale_factor=10000,
                       batch_var='pool')

betas1 <- betas_covid_full
betas2 <- betas_sle_full
f_use1 <- 1
f_use2 <- 1
# compute delta betas
for (ct in names(betas1[[f_use1]])) {
  genes1 <- names(betas1[[f_use1]][[ct]])
  genes2 <- names(betas2[[f_use2]][[ct]])
  gene_intersect <- intersect(genes1,genes2)
  betas1[[f_use1]][[ct]] <- betas1[[f_use1]][[ct]][gene_intersect]
  betas2[[f_use2]][[ct]] <- betas2[[f_use2]][[ct]][gene_intersect]
}

# for each cell type take the difference across datasets
betas_diff <- list()
for (ct in names(betas1[[f_use1]])) {
  betas_diff[[ct]] <- betas1[[f_use1]][[ct]] - betas2[[f_use2]][[ct]]
}

betas_diff_full <- betas_diff



## now subset the data objects to just the sex/age balanced set of donors
# subsetting cells first
# covid data
meta <- pbmc_container_covid[["scMinimal_full"]][["metadata"]]
counts <- pbmc_container_covid[["scMinimal_full"]][["count_data"]]
cells_keep <- rownames(meta)[meta$donors%in%cvd_dnr_keep]
meta <- meta[cells_keep,]
counts <- counts[,cells_keep]
pbmc_container_covid[["scMinimal_full"]][["metadata"]] <- meta
pbmc_container_covid[["scMinimal_full"]][["count_data"]] <- counts

# sle data
meta <- pbmc_container_SLE[["scMinimal_full"]][["metadata"]]
counts <- pbmc_container_SLE[["scMinimal_full"]][["count_data"]]
cells_keep <- rownames(meta)[meta$donors%in%sle_dnr_keep]
meta <- meta[cells_keep,]
counts <- counts[,cells_keep]
pbmc_container_SLE[["scMinimal_full"]][["metadata"]] <- meta
pbmc_container_SLE[["scMinimal_full"]][["count_data"]] <- counts

# now subsetting donor scores 
# covid
pbmc_container_covid$tucker_results[[1]] <- pbmc_container_covid$tucker_results[[1]][as.character(cvd_dnr_keep),]
# sle
pbmc_container_SLE$tucker_results[[1]] <- pbmc_container_SLE$tucker_results[[1]][as.character(sle_dnr_keep),]



## now rerunning get betas
betas_covid_sub <- get_betas(pbmc_container_covid, all_f_test=c(1), donor_min_cells=2,
                              norm_method='trim', scale_factor=10000,
                              batch_var = 'site')

betas_sle_sub <- get_betas(pbmc_container_SLE, all_f_test=c(1), donor_min_cells=20,
                            norm_method='trim', scale_factor=10000,
                            batch_var='pool')

betas1 <- betas_covid_sub
betas2 <- betas_sle_sub
f_use1 <- 1
f_use2 <- 1
# compute delta betas
for (ct in names(betas1[[f_use1]])) {
  genes1 <- names(betas1[[f_use1]][[ct]])
  genes2 <- names(betas2[[f_use2]][[ct]])
  gene_intersect <- intersect(genes1,genes2)
  betas1[[f_use1]][[ct]] <- betas1[[f_use1]][[ct]][gene_intersect]
  betas2[[f_use2]][[ct]] <- betas2[[f_use2]][[ct]][gene_intersect]
}

# for each cell type take the difference across datasets
betas_diff <- list()
for (ct in names(betas1[[f_use1]])) {
  betas_diff[[ct]] <- betas1[[f_use1]][[ct]] - betas2[[f_use2]][[ct]]
}

betas_diff_sub <- betas_diff

betas_all <- list(betas_diff_full,betas_diff_sub)
# saveRDS(betas_all,file='/home/jmitchel/data/lupus_data/sle_covid_matched_diff_betas.rds')
betas_all <- readRDS(file='/home/jmitchel/data/lupus_data/sle_covid_matched_diff_betas.rds')

betas_diff_full <- betas_all[[1]]
betas_diff_sub <- betas_all[[2]]





### to run gsea using full data
all_res <- lapply(1:length(betas_diff_full),function(x){
  ct <- names(betas_diff_full)[x]
  diff_vals <- betas_diff_full[[ct]]
  res <- as.data.frame(matrix(ncol=3,nrow=length(diff_vals)))
  colnames(res) <- c('padj','Z','Gene')
  res$Gene <- names(diff_vals)

  res$Z <- abs(diff_vals**2) # TESTING ONLY
  # res$Z <- abs(diff_vals**2.5) # TESTING ONLY
  
  res$padj <- rep(1,length(diff_vals))
  
  # sort by absolute value score
  res <- res[order(abs(res$Z),decreasing = TRUE),]
  res_outer <- list(res)
  names(res_outer) <- 'res'
  return(res_outer)
})
names(all_res) <- names(betas_diff_full)

cao <- cacoa::Cacoa$new(
  NULL, sample.groups=as.factor(c('ctrl','test')), 
  cell.groups=as.factor(c('ct1','ct2')),
  sample.per.cell=as.factor(c('1','2')),
  target.level='test', ref.level='ctrl', n.cores=5, verbose=FALSE
)
cao[["test.results"]][["de"]] <- all_res

cao$estimateOntology(type="GSEA", org.db=org.Hs.eg.db::org.Hs.eg.db, verbose=TRUE,
                     n.cores=5, ignore.cache=TRUE,
                     scoreType = "pos")

gsea_full <- cao[["test.results"]][["GSEA"]][["res"]]

### to run gsea using subset data
all_res <- lapply(1:length(betas_diff_sub),function(x){
  ct <- names(betas_diff_sub)[x]
  diff_vals <- betas_diff_sub[[ct]]
  res <- as.data.frame(matrix(ncol=3,nrow=length(diff_vals)))
  colnames(res) <- c('padj','Z','Gene')
  res$Gene <- names(diff_vals)
  
  res$Z <- abs(diff_vals**2) # TESTING ONLY
  # res$Z <- abs(diff_vals**2.5) # TESTING ONLY
  
  res$padj <- rep(1,length(diff_vals))
  
  # sort by absolute value score
  res <- res[order(abs(res$Z),decreasing = TRUE),]
  res_outer <- list(res)
  names(res_outer) <- 'res'
  return(res_outer)
})
names(all_res) <- names(betas_diff_sub)

cao <- cacoa::Cacoa$new(
  NULL, sample.groups=as.factor(c('ctrl','test')), 
  cell.groups=as.factor(c('ct1','ct2')),
  sample.per.cell=as.factor(c('1','2')),
  target.level='test', ref.level='ctrl', n.cores=5, verbose=FALSE
)
cao[["test.results"]][["de"]] <- all_res

cao$estimateOntology(type="GSEA", org.db=org.Hs.eg.db::org.Hs.eg.db, verbose=TRUE,
                     n.cores=5, ignore.cache=TRUE,
                     scoreType = "pos")

gsea_sub <- cao[["test.results"]][["GSEA"]][["res"]]


# compute fraction of original gene sets still significant after subsetting

for (ct in names(gsea_full)) {
  gsea_full_ct <- gsea_full[[ct]][['BP']]
  gsea_sub_ct <- gsea_sub[[ct]][['BP']]
  
  gsea_full_ct <- gsea_full_ct[gsea_full_ct$pvalue<.05,]
  gsea_sub_ct <- gsea_sub_ct[gsea_sub_ct$pvalue<.05,]
  orig_sets <- gsea_full_ct$Description
  sub_sets <- gsea_sub_ct$Description
  sum(orig_sets %in% sub_sets) / length(orig_sets)
}


# to get out the leading edge genes and look at the delta betas for them from full vs subset data

total_res <- NULL
for (ct in names(gsea_full)) {
  gsea_full_ct <- gsea_full[[ct]][['BP']]
  if (is.null(total_res)) {
    total_res <- rbind.data.frame(cbind.data.frame(gsea_full_ct@result,rep(ct,nrow(gsea_full_ct@result))))
  } else {
    total_res <- rbind.data.frame(total_res,cbind.data.frame(gsea_full_ct@result,rep(ct,nrow(gsea_full_ct@result))))
  }
}
colnames(total_res)[ncol(total_res)] <- 'cell_type'
total_res$p.adjust2 <- p.adjust(total_res$pvalue,method='fdr')

total_res_sub <- total_res[total_res$cell_type=='NK',]
sum(total_res_sub$p.adjust<.05)
sum(total_res_sub$p.adjust2<.05)

total_res_full <- total_res
total_res_full_ct <- total_res_sub
total_res_full_ct <- total_res_full_ct[total_res_full_ct$p.adjust2<.05,]

total_res <- NULL
for (ct in names(gsea_sub)) {
  gsea_full_ct <- gsea_sub[[ct]][['BP']]
  if (is.null(total_res)) {
    total_res <- rbind.data.frame(cbind.data.frame(gsea_full_ct@result,rep(ct,nrow(gsea_full_ct@result))))
  } else {
    total_res <- rbind.data.frame(total_res,cbind.data.frame(gsea_full_ct@result,rep(ct,nrow(gsea_full_ct@result))))
  }
}
colnames(total_res)[ncol(total_res)] <- 'cell_type'
total_res$p.adjust2 <- p.adjust(total_res$pvalue,method='fdr')

total_res_sub <- total_res[total_res$cell_type=='NK',]
sum(total_res_sub$p.adjust<.05)
sum(total_res_sub$p.adjust2<.05)

total_res_sub <- total_res_sub[total_res_sub$p.adjust2<.5,]

sum(total_res_sub$Description %in% total_res_full_ct$Description)




# to get out the leading edge genes and look at the delta betas for them from full vs subset data
# make df with columns as full delta beta, subset delta beta, gene_type (leading edge/other)
df1_total <- NULL
df2_total <- NULL
for (ct in names(cao[["test.results"]][["GSEA"]][["res"]])) {
  gsea_res <- cao[["test.results"]][["GSEA"]][["res"]][[ct]][['BP']]
  gsea_res_sig <- gsea_res[gsea_res$p.adjust<.05,]
  genes <- gsea_res_sig$core_enrichment
  genes <- sapply(genes,function(x){
    unlist(strsplit(x,split='/'))
  })
  genes <- unlist(genes)
  genes <- unique(genes)

  genes1 <- genes[genes %in% names(betas_diff_full[[ct]])]
  genes2 <- genes[genes %in% names(betas_diff_sub[[ct]])]
  genes <- intersect(genes1,genes2)
  
  other_genes1 <- names(betas_diff_full[[ct]])[!(betas_diff_full[[ct]] %in% genes)]
  other_genes2 <- names(betas_diff_sub[[ct]])[!(betas_diff_sub[[ct]] %in% genes)]
  other_genes <- intersect(other_genes1,other_genes2)
  
  # trying without squaring...
  df1 <- cbind.data.frame(betas_diff_full[[ct]][genes],betas_diff_sub[[ct]][genes],rep('le_genes',length(genes)))
  colnames(df1) <- c('full','sub','type')
  df2 <- cbind.data.frame(betas_diff_full[[ct]][other_genes],
                          betas_diff_sub[[ct]][other_genes],rep('other_genes',length(other_genes)))
  colnames(df2) <- c('full','sub','type')
  # df3 <- rbind.data.frame(df1,df2)
  
  if (is.null(df1_total)) {
    df1_total <- df1
  } else {
    df1_total <- rbind.data.frame(df1_total,df1)
  }
  
  if (is.null(df2_total)) {
    df2_total <- df2
  } else {
    df2_total <- rbind.data.frame(df2_total,df2)
  }
  
}


df2_total$type <- as.factor(df2_total$type)
levels(df2_total$type)
levels(df2_total$type) <- c('other')
df1_total$type <- as.factor(df1_total$type)
levels(df1_total$type)
levels(df1_total$type) <- c('leading edge')

library(ggrastr)
p <- ggplot(df2_total,aes(x=full,y=sub,color=type)) +
  ggrastr::rasterise(geom_point(alpha=.1)) +
  ggrastr::rasterise(geom_point(data=df1_total,aes(x=full,y=sub,color=type),alpha=.1)) +
  geom_abline(intercept = 0, slope = 1, colour = 'blue') +
  geom_hline(yintercept = 0, colour = 'black') +
  geom_vline(xintercept = 0, colour = 'black') +
  scale_color_manual(values=c( "red","black")) +
  xlab('Full dataset - delta beta') +
  ylab('Balanced dataset - delta beta') +
  theme_bw(base_size = 17) +
  guides(color=guide_legend(title="Gene category"))

### Figure S10d
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/full_vs_balanced_delta_betas.pdf", useDingbats = FALSE,
#     width = 8, height = 5)
p
# dev.off()
