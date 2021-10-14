
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


pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, 
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')



pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20),
                                 tucker_type = 'regular', rotation_type = 'hybrid')
pca_unfolded(pbmc_container,9) # testing using PCA on unfolded tensor instead of Tucker

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(9,20),
                                 tucker_type = 'regular', rotation_type = 'ica_dsc') #for dmat rot

# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20),
#                                  tucker_type = 'regular', rotation_type = 'ica_dsc') 
# 
# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(25,37),
#                                  tucker_type = 'regular', rotation_type = 'ica_dsc') 

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




# testing for increase in Treg cell subset out of Th cells
pbmc_container$subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')

test <- pbmc@meta.data[names(t4_sub),'ct_cov']
ndx_keep <- which(test!='T4naive')
test <- test[ndx_keep]

tmp <- test
ndx_mark <- which(test=='T4reg')
ndx_other <- which(test!='T4reg')

tmp <- as.character(tmp)
tmp[ndx_mark] <- 1
tmp[ndx_other] <- 2
names(tmp) <- names(t4_sub)[ndx_keep]

pbmc_container[["subclusters"]][["T4"]][["res:0.6"]] <- as.numeric(tmp)
names(pbmc_container[["subclusters"]][["T4"]][["res:0.6"]]) <- names(t4_sub)[ndx_keep]

myplot <- get_subclust_enr_dotplot(pbmc_container,'T4',0.6,subtype=1,factor_use=1,'Th')
pdf(file = "/home/jmitchel/figures/for_paper_v2/Treg_props.pdf", useDingbats = FALSE,
    width = 4.5, height = 3.75)
myplot + ylim(0,.35) + ggtitle('Treg proportions')
dev.off()

container <- pbmc_container
ctype <- 'T4'
res <- .6
resolution_name <- paste0('res:',as.character(res))
subclusts <- container$subclusters[[ctype]][[resolution_name]]

# append large cell type name to subclusters
subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})

# limit cells in subclusts to those that we actually have scores for
donor_scores <- container$tucker_results[[1]]
donor_vec <- container$scMinimal_full$metadata[names(subclusts),'donors']
subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]

# make subtype association plot
subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
scMinimal <- container$scMinimal_ctype[['Th']]
sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

# get donor proportions of subclusters
donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)

f1 <- get_one_factor(container,1)
f1_dsc <- f1[[1]]
tmp <- cbind.data.frame(f1_dsc[rownames(donor_props),1],donor_props)
colnames(tmp) <- c('dsc','Treg','Tother')
lmres <- summary(lm(Treg~dsc,data=tmp))
pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
print(pval)



## testing whether there are any residual gc content associations after batch correction
all_pv <- c()
for (i in 1:25) {
  print(i)
  pv <- test_gc_association(pbmc_container,my_factor=i,b_direc='up',comp_type='any_up')
  all_pv <- c(all_pv,pv[[1]])
  names(all_pv)[length(all_pv)] <- i
  pv <- test_gc_association(pbmc_container,my_factor=i,b_direc='down',comp_type='any_up')
  all_pv <- c(all_pv,pv[[1]])
  names(all_pv)[length(all_pv)] <- i
}
all_pv2 <- p.adjust(all_pv,method='fdr')
all_pv2
min(all_pv2)

for (i in 1:50) {
  if (i %% 2 == 1) {
    print(ceil(i/2))
  }
  print(all_pv2[i])
}


# testing whether there are enriched gene sets among the high gc content genes
# at some arbitrary cutoff
head(pbmc_container[["gen_dat"]])
test <- pbmc_container[["gen_dat"]]
dim(test)
up_gc <- test$hgnc_symbol[test$gc>.42]
length(up_gc)
down_gc <- test$hgnc_symbol[test$gc<=.42]
length(down_gc)

db_use <- 'GO'
m_df <- data.frame()
for (db in db_use) {
  if (db == "GO") {
    # select the GO Biological Processes group of gene sets
    m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                        category = "C5", subcategory = "BP"))
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
  }
}

my_pathways = split(m_df$gene_symbol, f = m_df$gs_name)
# my_pathways <- split(m_df$gene_symbol, f = m_df$gs_exact_source)

all_genes <- unique(c(up_gc,down_gc))
total_num_genes <- length(all_genes)
# sig_genes <- up_gc
sig_genes <- down_gc


pvals <- c()
for (i in 1:length(my_pathways)) {
  pth <- my_pathways[[i]]
  pth_name <- names(my_pathways)[i]
  
  # A: total num genes in pathway in tmp_casted_num
  pth_in_df <- pth[which(pth %in% all_genes)]
  num_pth_in_df <- length(pth_in_df)
  
  # if set is too small continue
  if (num_pth_in_df < 15) {
    next
  }
  
  # B: number of genes from A in sig_genes
  num_in_sig <- sum(pth_in_df %in% sig_genes)
  
  # compute pvalue
  pval <- stats::phyper(num_in_sig-1, num_pth_in_df, total_num_genes - num_pth_in_df,
                        length(sig_genes), lower.tail = FALSE) # I double checked this is right
  pvals[pth_name] <- pval
}
padj <- p.adjust(pvals,method='fdr')

padj[order(padj,decreasing=FALSE)][1:5]















# plotting Treg proportions with healthy/SLE labels
pbmc_container$subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')
t4_sub <- colnames(pbmc_container$scMinimal_ctype[['Th']]$count_data)

test <- pbmc@meta.data[t4_sub,'ct_cov']
# ndx_keep <- which(test!='T4naive')
# test <- test[ndx_keep]

tmp <- test
ndx_mark <- which(test=='T4reg')
ndx_other <- which(test!='T4reg')

tmp <- as.character(tmp)
tmp[ndx_mark] <- 1
tmp[ndx_other] <- 2
# names(tmp) <- t4_sub[ndx_keep]

pbmc_container[["subclusters"]][["T4"]][["res:0.6"]] <- as.numeric(tmp)
names(pbmc_container[["subclusters"]][["T4"]][["res:0.6"]]) <- t4_sub

container <- pbmc_container
ctype <- 'T4'
res <- 0.6
subtype=1
factor_use=1
ctype_cur='Th'

resolution_name <- paste0('res:',as.character(res))
subclusts <- container$subclusters[[ctype]][[resolution_name]]

# append large cell type name to subclusters
subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})

# limit cells in subclusts to those that we actually have scores for
donor_scores <- container$tucker_results[[1]]
cell_intersect <- intersect(names(subclusts),rownames(container$scMinimal_full$metadata))
donor_vec <- container$scMinimal_full$metadata[cell_intersect,'donors']
subclusts <- subclusts[cell_intersect]
subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]

# make subtype association plot
subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
scMinimal <- container$scMinimal_ctype[[ctype_cur]]
sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

# get donor proportions of subclusters
donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)
donor_props <- donor_props[,subtype,drop=FALSE]
colnames(donor_props) <- 'prop'

# append dscores for factor 4
donor_props2 <- cbind(donor_props,donor_scores[rownames(donor_props),factor_use])
colnames(donor_props2)[ncol(donor_props2)] <- 'dsc'

# append disease status
meta <- unique(container$scMinimal_full$metadata[,c('donors','Status')])
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
  xlab(paste0('Factor ',as.character(factor_use),' Donor Score')) +
  ylab(paste0('Proportion Treg/Th')) +
  ylim(0,1) +
  labs(color = "Status") +
  ggtitle(paste0(ctype,'_',as.character(subtype),' Proportions')) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  scale_color_manual(values = c("Healthy" = '#F8766D',
                                "Managed" = '#00BFC4',
                                "1" = "black")) + 
  ylim(0,.25) + ggtitle('')

pdf(file = "/home/jmitchel/figures/for_paper_v2/Treg_props.pdf", useDingbats = FALSE,
    width = 5, height = 3.25)
p
dev.off()



### get median number cells per donor and total number of cells used
d_keep <- rownames(pbmc_container$scMinimal_ctype[[1]]$pseudobulk)
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$ind_cov_batch_cov %in% d_keep]
pbmc_sub <- subset(pbmc,cells=cells_keep)
ctypes <- c("B","NK","Th","Tc","cDC","cMono","ncMono")
cells_keep <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$cg_cov %in% ctypes]
pbmc_sub2 <- subset(pbmc_sub,cells=cells_keep)
pbmc_sub2@meta.data$ind_cov_batch_cov <- factor(as.character(pbmc_sub2@meta.data$ind_cov_batch_cov),levels=unique(as.character(pbmc_sub2@meta.data$ind_cov_batch_cov)))
for (ctype in ctypes) {
  tmp <- pbmc_sub2@meta.data[pbmc@meta.data$cg_cov==ctype,]
  print(ctype)
  print(median(table(tmp$ind_cov_batch_cov)))
}










