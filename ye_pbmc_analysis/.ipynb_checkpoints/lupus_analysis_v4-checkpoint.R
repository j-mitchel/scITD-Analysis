
library(Seurat)

# load up the subsetted dataset
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')


# converting shorthand cell type names to full names
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
# param_list <- initialize_params(ctypes_use = c("B","NK","T4","T8","cDC",
#                                                "cM","ncM"),
#                                 ncores = 30, rand_seed = 10)
param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc","cDC",
                                               "cMono","ncMono"),
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


pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=0,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')



pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(9,20,7),
                                 tucker_type = 'regular', rotation_type = 'ica') #for dmat rot

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Status','Ethnicity'),
                                        stat_use='pval')

# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    cluster_by_meta = 'sex',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','processing','pool'),
                                    cluster_by_meta = 'pool',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

pdf(file = "/home/jmitchel/figures/for_paper_v2/lupus_dscores_full.pdf", useDingbats = FALSE,
    width = 6, height = 7)
pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 6, height = 7)
pbmc_container$plots$donor_matrix
dev.off()



# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)



# get loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=0.02,
                                           display_genes=FALSE,
                                           gene_callouts=TRUE,
                                           callout_n_gene_per_ctype=10)

myfig <- render_multi_plots(pbmc_container,data_type='loadings',max_cols=3)

pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 18, height = 23)
myfig
dev.off()


# run gsea for a f1
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=1, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 14, height = 14)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=1,direc='down',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))
dev.off()

pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 14, height = 14)
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=7)
dev.off()

## f4 sets to show on loading hmap
gsets <- c("GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
           "GOBP_RESPONSE_TO_INTERFERON_GAMMA",
           "GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY",
           "GOBP_PATTERN_RECOGNITION_RECEPTOR_SIGNALING_PATHWAY",
           "GOBP_RECEPTOR_SIGNALING_PATHWAY_VIA_STAT",
           "GOBP_INTERLEUKIN_1_PRODUCTION",
           "GOBP_MYELOID_LEUKOCYTE_ACTIVATION",
           "GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_PROLIFERATION",
           "GOBP_REGULATORY_T_CELL_DIFFERENTIATION")

gset_cmap <- c('blue',
               'black',
               'black',
               'black',
               'black',
               'black',
               'orange',
               'black',
               'forest green')

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
                            cl_rows=F, myfontsize=6.5, h_w=c(6,6.5))

p1 <- pbmc_container$plots$all_lds_plots[['1']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["1"]]

pdf(file = "/home/jmitchel/figures/for_paper_v2/lupus_f1_lds_go.pdf", useDingbats = FALSE,
    width = 12, height = 10)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()














# plotting IFN response gene expression against f1
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Status')]
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

pdf(file = "/home/jmitchel/figures/for_paper_v2/f1_IFN_status.pdf", useDingbats = FALSE,
    width = 5, height = 2.75)
# mycol <- RColorBrewer::brewer.pal(n = 3, name = "Accent")
mycol <- RColorBrewer::brewer.pal(n = 6, name = "Dark2")
ggplot(tmp,aes(x=dscore,y=expres,color=status)) +
  geom_point(alpha = 0.75,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy,color='line')) +
  scale_color_manual(values=c(mycol[1],"#000000",mycol[6])) +
  ylab('IFI6 expression (CD4+ T)') +
  xlab('Factor 1 donor scores') +
  theme_bw()
dev.off()





## checking TLR7 expression in cMono and ncMono

test <- get_one_factor(pbmc_container,1)
dsc <- test[[1]]
lds <- test[[2]]

sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=1, colnames(lds))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

# limit to just the genes in tmp_casted_num
sig_df <- sig_df[rownames(lds),colnames(lds)]

sig_df['TLR7',]

dsc <- dsc[,1]
dsc <- dsc[order(dsc)]

meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Status','sex')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

tmp <- cbind.data.frame(meta[names(dsc),],dsc)
tmp[1:50,]



## looking at distributions for some factors using the different rotations
test=get_one_factor(pbmc_container,3)

pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 5, height = 3.75)
hist(test[[1]][,1])
dev.off()

pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 5, height = 3.75)
qqnorm(test[[1]][,1], pch = 1, frame = FALSE)
qqline(test[[1]][,1], col = "steelblue", lwd = 2)
dev.off()




