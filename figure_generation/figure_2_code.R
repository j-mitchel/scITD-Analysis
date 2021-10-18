library(scITD)
library(Seurat)
library(ggplot2)
library(coda.base)
library(RColorBrewer)

# load up the subsetted dataset
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


pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, 
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20),
                                 tucker_type = 'regular', rotation_type = 'hybrid')

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Status','Ethnicity'),
                                        stat_use='pval')

# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    cluster_by_meta = 'sex',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# pdf(file = "/home/jmitchel/figures/for_paper_v2/lupus_dscores_full2.pdf", useDingbats = FALSE,
#     width = 6, height = 7)
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


# run gsea for a f1
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=1, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=1,direc='down',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=7)

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

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, 
                                      le_set_callouts=gset_sub, le_set_colormap=gset_cmap_sub, le_set_num_per=7, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                      clust_method='mcquitty', h_w=c(9,6.5))

hm_list <- plot_select_sets(pbmc_container, 1, gsets, color_sets=gset_cmap, 
                            cl_rows=F, myfontsize=6.5, h_w=c(6,6.5))

p1 <- pbmc_container$plots$all_lds_plots[['1']]
p2 <- p1 %v% hm_list[[1]]

pd <- pbmc_container[["plots"]][["all_legends"]][["1"]]

pdf(file = "/home/jmitchel/figures/for_paper_v2/lupus_f1_lds_go2.pdf", useDingbats = FALSE,
    width = 12, height = 10)
draw(p2,annotation_legend_list = pd,
     legend_grouping = "original", annotation_legend_side = "left",
     heatmap_legend_list = hm_list[[2]], heatmap_legend_side = "left",
     newpage=TRUE, auto_adjust = FALSE)
dev.off()




##### plotting IFN response gene expression against f1
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

mycol <- RColorBrewer::brewer.pal(n = 6, name = "Dark2")
ifi6_dot <- ggplot(tmp,aes(x=dscore,y=expres,color=status)) +
  geom_point(alpha = 0.75,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy,color='line')) +
  scale_color_manual(values=c(mycol[1],"#000000",mycol[6])) +
  ylab('IFI6 expression (CD4+ T)') +
  xlab('Factor 1 donor scores') +
  theme_bw()

pdf(file = "/home/jmitchel/figures/for_paper_v2/f1_IFN_status2.pdf", useDingbats = FALSE,
    width = 5, height = 2.75)
ifi6_dot
dev.off()


# get association p-value for IFI6
print(pbmc_container[["gene_score_associations"]]['IFI6.Th.1'])








##### plotting Treg proportions with healthy/SLE labels
pbmc_container$subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')
t4_sub <- colnames(pbmc_container$scMinimal_ctype[['Th']]$count_data)

tmp <- pbmc@meta.data[t4_sub,'ct_cov']

ndx_mark <- which(tmp=='T4reg')
ndx_other <- which(tmp!='T4reg')

tmp <- as.character(tmp)
tmp[ndx_mark] <- 1
tmp[ndx_other] <- 2

pbmc_container[["subclusters"]][["T4"]][["res:0.6"]] <- as.numeric(tmp)
names(pbmc_container[["subclusters"]][["T4"]][["res:0.6"]]) <- t4_sub

ctype <- 'T4'
res <- 0.6
subtype=1
factor_use=1
ctype_cur='Th'

resolution_name <- paste0('res:',as.character(res))
subclusts <- pbmc_container$subclusters[[ctype]][[resolution_name]]

# append large cell type name to subclusters
subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})

# limit cells in subclusts to those that we actually have scores for
donor_scores <- pbmc_container$tucker_results[[1]]
cell_intersect <- intersect(names(subclusts),rownames(pbmc_container$scMinimal_full$metadata))
donor_vec <- pbmc_container$scMinimal_full$metadata[cell_intersect,'donors']
subclusts <- subclusts[cell_intersect]
subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]

# make subtype association plot
subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
scMinimal <- pbmc_container$scMinimal_ctype[[ctype_cur]]
sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

# get donor proportions of subclusters
donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)
donor_props <- donor_props[,subtype,drop=FALSE]
colnames(donor_props) <- 'prop'

# append dscores for factor 1
donor_props2 <- cbind(donor_props,donor_scores[rownames(donor_props),factor_use])
colnames(donor_props2)[ncol(donor_props2)] <- 'dsc'

# append disease status
meta <- unique(pbmc_container$scMinimal_full$metadata[,c('donors','Status')])
rownames(meta) <- meta$donors
donor_props2 <- cbind(donor_props2,as.character(meta[rownames(donor_props2),'Status']))
colnames(donor_props2)[ncol(donor_props2)] <- 'Status'

donor_props2 <- as.data.frame(donor_props2)
donor_props2$dsc <- as.numeric(donor_props2$dsc)
donor_props2$prop <- as.numeric(donor_props2$prop)
donor_props2$Status <- as.factor(donor_props2$Status)

lmres <- lm(prop~dsc,data=donor_props2)
line_range <- seq(min(donor_props2$dsc),max(donor_props2$dsc),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
line_df <- cbind.data.frame(line_df,rep('1',nrow(line_df)))
colnames(line_df) <- c('myx','myy','Status')

p <- ggplot(donor_props2,aes(x=dsc,y=prop,color=Status)) +
  geom_point(alpha = 0.75,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  scale_color_manual(values=c("#000000",mycol[1],mycol[6])) +
  xlab(paste0('Factor ',as.character(factor_use),' Donor Score')) +
  ylab(paste0('Proportion Treg/Th')) +
  ylim(0,1) +
  labs(color = "Status") +
  ggtitle(paste0(ctype,'_',as.character(subtype),' Proportions')) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  ylim(0,.25) + ggtitle('')

pdf(file = "/home/jmitchel/figures/for_paper_v2/Treg_props2.pdf", useDingbats = FALSE,
    width = 5, height = 3.25)
p
dev.off()

donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)
donor_balances <- coda.base::coordinates(donor_props)
rownames(donor_balances) <- rownames(donor_props)

f1 <- get_one_factor(pbmc_container,1)
f1_dsc <- f1[[1]]
tmp <- cbind.data.frame(f1_dsc[rownames(donor_balances),1,drop=FALSE],donor_balances)
colnames(tmp) <- c('dsc','my_balance')
lmres <- summary(lm(my_balance~dsc,data=tmp))
pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
print(pval)










##### plotting embedding of all cells from SLE dataset
# conos object from embedding prep file
con <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_conos2.rds')
tmp <- con$plotGraph(alpha=0.1)
mycolors <- brewer.pal(n = 9, name = "Set1")
mycolors <- mycolors[c(1:5,7,9)]
tmp <- tmp + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black", size=1, fill=NA)) +
  scale_color_manual(values=mycolors) +
  scale_y_reverse() +
  scale_x_reverse()
# scale_color_brewer(palette="Set1")
tmp$layers[[2]] <- NULL


# pdf(file = "/home/jmitchel/figures/for_paper/lupus_embedding.pdf", useDingbats = FALSE,
#     width = 7, height = 7)
# saved a jpeg 550 x 400 dimensions
tmp
# dev.off()
