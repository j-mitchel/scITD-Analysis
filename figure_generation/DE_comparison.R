
library(Seurat)
library(ggplot2)
library(ggrastr)
library(readxl)
library(ComplexHeatmap)
library(MASS)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(sccore)
library(reshape2)
library(ggbeeswarm)
library(ggrepel)
library(devtools)
load_all('/home/jmitchel/scITD/')

# load up the lupus dataset: see preprocessing/lupus_preprocessing.R 
# for code used to generate this object
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


# subset data to SLE patients only
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$Status=='Managed']
pbmc <- subset(pbmc,cells = cells_keep)

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

# flip sign of F1 so high ISG expression is positive instead of negative (signs are arbitrary)
pbmc_container$tucker_results[[1]][,1] <- pbmc_container$tucker_results[[1]][,1] * -1
pbmc_container$tucker_results[[2]][1,] <- pbmc_container$tucker_results[[2]][1,] * -1
pbmc_container$projection_data[[1]][1,] <- pbmc_container$projection_data[[1]][1,] * -1



#### recomputing tensor with all expressed genes, not just ones used in the tensor
pbmc_container_full <- make_new_container(seurat_obj=pbmc,
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


pbmc_container_full <- form_tensor(pbmc_container_full, donor_min_cells=20,
                                   norm_method='trim', scale_factor=10000,
                                   vargenes_method='norm_var_pvals', vargenes_thresh=.5,
                                   scale_var = TRUE, var_scale_power = .5,
                                   batch_var='pool')

pbmc_container_full$tucker_results <- pbmc_container$tucker_results

# get significant genes
pbmc_container_full <- get_lm_pvals(pbmc_container_full)
pbmc_container <- pbmc_container_full

## using ordinal variables
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_ordinal.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL

# get tucker donor scores to test
dsc <- pbmc_container$tucker_results[[1]]

## get donors IDs matching to the clinical variable IDs
trim_names <- sapply(rownames(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})

# get donors in both dataframes
d_both <- rownames(clin_vars)[rownames(clin_vars) %in% trim_names]
trim_names <- trim_names[trim_names %in% d_both]

de_genes <- c()
meta_var <- 'sledaiscore'
for (ct in pbmc_container$experiment_params$ctypes_use) {
  print(ct)
  # getting out pseudobulk for a cell type
  pb <- pbmc_container$scMinimal_ctype[[ct]]$pseudobulk
  rownames(pb) <- trim_names[rownames(pb)]
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  
  # run DE t-test
  ct_res <- plapply(1:ncol(pb),function(j) {
    tmp <- cbind.data.frame(pb[,j],clin_vars[rownames(pb),meta_var])
    colnames(tmp) <- c('expr','meta')
    # trying linear test instead
    lmres <- summary(lm(meta~expr,data=tmp))
    if (all(tmp$expr==0)) {
      tres_pval <- NaN
    } else {
      tres_pval <- lmres$coefficients['expr','Pr(>|t|)']
    }
    names(tres_pval) <- colnames(pb)[j]
    return(tres_pval)
  },mc.preschedule=TRUE,n.cores=20,progress=TRUE)
  
  de_genes <- c(de_genes,ct_res)
}
de_genes <- unlist(de_genes)

# correct pvals
pv_thresh <- .1
de_genes_adj <- p.adjust(de_genes,method = 'fdr')
de_genes_adj <- de_genes_adj[!is.na(de_genes_adj)]
de_genes_id <- names(de_genes_adj)[de_genes_adj<pv_thresh]
print(length(de_genes_id))

## generate correlation matrix for DE significant genes
# first need to append expression for all cell types
pb_all <- matrix(nrow=length(trim_names),ncol=0)
for (ct in pbmc_container$experiment_params$ctypes_use) {
  # getting out pseudobulk for a cell type
  pb <- pbmc_container$scMinimal_ctype[[ct]]$pseudobulk
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  g_both <- intersect(colnames(pb),de_genes_id)
  if (length(g_both)==0) {
    next
  }
  pb <- pb[names(trim_names),g_both,drop=FALSE]
  pb_all <- cbind(pb_all,pb)
}
de_cormat <- abs(cor(pb_all,method = 'pearson'))


# cluster genes and run gsea
hres <- hclust(as.dist(1-de_cormat),method = 'ward.D2')
hclusts <- cutree(hres, k = 4)
table(hclusts)

## manually reordering rows and cols by clustering
de_cormat2 <- de_cormat[hres[["order"]],hres[["order"]]]
hclusts <- hclusts[hres[["order"]]]

g_ct_cev_de <- sapply(colnames(de_cormat2),function(x){
  return(strsplit(x,split='_')[[1]][[2]])
})

# set up color schemes
group_col <- c('darkblue','darkgreen','darkorange','red4')
names(group_col) <- c(1,2,3,4)
annot_col <- brewer.pal(7, 'Set3')
names(annot_col) <- unique(g_ct_cev_de)
col_fun = colorRamp2(c(0, 1), c("white", "red"))

# create annotations
ra_de <- ComplexHeatmap::rowAnnotation(ctype = g_ct_cev_de, show_annotation_name=FALSE, col = list(ctype = annot_col))
ca_grp <- ComplexHeatmap::HeatmapAnnotation(hcluster = hclusts,
                                            col = list(hcluster = group_col),
                                            show_legend = TRUE)

hmap_de <- Heatmap(de_cormat2,name = "pearson r",
                   cluster_columns = FALSE,
                   col = col_fun,
                   cluster_rows = FALSE,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   show_row_dend = FALSE,
                   show_column_dend = FALSE,
                   top_annotation=ca_grp,
                   right_annotation=ra_de,
                   border = TRUE,
                   column_title = 'SLEDAI DE genes',
                   column_title_gp = gpar(fontsize = 20))
hmap_de

bg <- names(de_genes)[!is.na(de_genes)]
bg <- sapply(bg,function(x){
  strsplit(x,split='_')[[1]][[1]]
})
bg <- unique(bg)

c_g <- names(hclusts)[hclusts %in% c(3,4)]
c_g <- sapply(c_g,function(x){
  strsplit(x,split='_')[[1]][[1]]
})

go_res_c <- enrichGO(gene = c_g,
                      universe = bg,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'SYMBOL',
                      ont           = 'BP',
                      pAdjustMethod = "BH",
                      pvalueCutoff  = .01,
                      qvalueCutoff  = .05,
                      minGSSize=10,
                      maxGSSize=500)
go_res_c <- go_res_c@result
### some significant gene sets for cluster 3 include:
# "catecholamine secretion"
# "regulation of chromosome segregation"

## now trying for just cluster 4
c_g <- names(hclusts)[hclusts==4]
c_g <- sapply(c_g,function(x){
  strsplit(x,split='_')[[1]][[1]]
})

go_res_c <- enrichGO(gene = c_g,
                     universe = bg,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = 'BP',
                     pAdjustMethod = "BH",
                     pvalueCutoff  = .01,
                     qvalueCutoff  = .05,
                     minGSSize=10,
                     maxGSSize=500)
go_res_c <- go_res_c@result


## now getting correlation matrix for scITD F1 significant genes
f1_pvals <- get_one_factor_gene_pvals(pbmc_container,1)

# trying to just get top n gene_ct pairs to match number in DE
pv_ord <- order(f1_pvals,decreasing=FALSE)[1:ncol(de_cormat)]
top_indexes <- as.data.frame(matrix(0, length(pv_ord), 2))
colnames(top_indexes) <- c("Row", "Column")
top_indexes$Row <- as.integer((pv_ord-1) %% nrow(f1_pvals)) + 1
top_indexes$Column <- as.integer((pv_ord-1) %/% nrow(f1_pvals)) + 1
f1_pvals2 <- matrix(1,nrow=nrow(f1_pvals),ncol=ncol(f1_pvals))
for (i in 1:nrow(top_indexes)) {
  f1_pvals2[top_indexes$Row[i],top_indexes$Column[i]] <-   f1_pvals[top_indexes$Row[i],top_indexes$Column[i]]
}
colnames(f1_pvals2) <- colnames(f1_pvals)
rownames(f1_pvals2) <- rownames(f1_pvals)

pb_all <- matrix(nrow=length(trim_names),ncol=0)
for (ct in pbmc_container$experiment_params$ctypes_use) {
  # getting out pseudobulk for a cell type
  pb <- pbmc_container$scMinimal_ctype[[ct]]$pseudobulk
  pb <- pb[names(trim_names),rownames(f1_pvals2)[f1_pvals2[,ct]<pv_thresh]]
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  # append pb to existing pb
  pb_all <- cbind(pb_all,pb)
}

scITD_cormat <- abs(cor(pb_all,method='pearson'))



### running go for significant genes by either method
de_genes_ct <- sapply(de_genes_id,function(x) {
  strsplit(x,split='_')[[1]][[2]]
})
de_genes_only <- sapply(de_genes_id,function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
scITD_genes_id <- colnames(scITD_cormat)
scITD_genes_ct <- sapply(scITD_genes_id,function(x) {
  strsplit(x,split='_')[[1]][[2]]
})
scITD_genes_only <- sapply(scITD_genes_id,function(x) {
  strsplit(x,split='_')[[1]][[1]]
})


de_genes_unq <- unique(de_genes_only)
scITD_genes_unq <- unique(scITD_genes_only)

go_bp_res_de <- enrichGO(gene = de_genes_unq,
                      universe = bg,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'SYMBOL',
                      ont           = 'BP',
                      pAdjustMethod = "BH",
                      pvalueCutoff  = .01,
                      qvalueCutoff  = .05,
                      minGSSize=10,
                      maxGSSize=1500)


go_bp_res_de <- go_bp_res_de@result

go_bp_res_f1 <- enrichGO(gene = scITD_genes_unq,
                      universe = bg,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'SYMBOL',
                      ont           = 'BP',
                      pAdjustMethod = "BH",
                      pvalueCutoff  = .01,
                      qvalueCutoff  = .05,
                      minGSSize=10,
                      maxGSSize=1500)


go_bp_res_f1 <- go_bp_res_f1@result

## seeing if there are significant sets in scITD not in de
go_res_f1_sub <- go_bp_res_f1$Description[go_bp_res_f1$p.adjust<.1]
go_res_de_ns <- go_bp_res_de[go_bp_res_de$pvalue>.1,]
go_res_de_ns <- go_res_de_ns$Description[go_res_de_ns$Count<=4]
go_res_f1_sub[go_res_f1_sub%in%go_res_de_ns]


## response to type I interferon enrichment 
# for de: padj=1.39829e-05
# for de cluster 1: padj=3.934150e-09
# for scITD: padj=1.660958e-16


go_mf_res_de <- enrichGO(gene = de_genes_unq,
                         universe = bg,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = 'MF',
                         pAdjustMethod = "BH",
                         pvalueCutoff  = .01,
                         qvalueCutoff  = .05,
                         minGSSize=10,
                         maxGSSize=1500)


go_mf_res_de <- go_mf_res_de@result

go_mf_res_f1 <- enrichGO(gene = scITD_genes_unq,
                         universe = bg,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = 'MF',
                         pAdjustMethod = "BH",
                         pvalueCutoff  = .01,
                         qvalueCutoff  = .05,
                         minGSSize=10,
                         maxGSSize=1500)


go_mf_res_f1 <- go_mf_res_f1@result

go_res_f1_sub <- go_mf_res_f1$Description[go_mf_res_f1$p.adjust<.05]
go_res_de_ns <- go_mf_res_de[go_mf_res_de$pvalue>.1,]
go_res_de_ns <- go_res_de_ns$Description[go_res_de_ns$Count<=2]
go_res_f1_sub[go_res_f1_sub%in%go_res_de_ns]


### computing overlap for either de set or our comparable top set with the in vitro ifn stimulation data
load("/home/jmitchel/data/lupus_data/cell.type.diffexp.RData")

B.expressed.res.sig <- cd19.expressed.res$featureData.symbol[cd19.expressed.res$padj<.05]
nk.expressed.res.sig <- nk.expressed.res$featureData.symbol[nk.expressed.res$padj<.05]
cM.expressed.res.sig <- cd14.expressed.res$featureData.symbol[cd14.expressed.res$padj<.05]
cd4.expressed.res.sig <- cd4.expressed.res$featureData.symbol[cd4.expressed.res$padj<.05]
cd8.expressed.res.sig <- cd8.expressed.res$featureData.symbol[cd8.expressed.res$padj<.05]
ncM.expressed.res.sig <- ncmono.expressed.res$featureData.symbol[ncmono.expressed.res$padj<.05]

ct_test <- list(B.expressed.res.sig,nk.expressed.res.sig,cM.expressed.res.sig,cd4.expressed.res.sig,cd8.expressed.res.sig,ncM.expressed.res.sig)
ct_test_full <- list(cd19.expressed.res,nk.expressed.res,cd14.expressed.res,cd4.expressed.res,cd8.expressed.res,ncmono.expressed.res)
names(ct_test) <- c('B','NK','cMono','Th','Tc','ncMono')
names(ct_test_full) <- c('B','NK','cMono','Th','Tc','ncMono')

total_intersect_de <- 0
total_intersect_scITD <- 0
num_sig_in_dat <- 0
num_tested_in_dat <- 0
intersect_genes_de <- c()
intersect_genes_scITD <- c()
for (ct in names(ct_test)) {
  iv_res <- ct_test[[ct]]
  de_res <- de_genes_only[de_genes_ct==ct]
  de_intersect <- intersect(de_res,iv_res)
  ilen_de <- length(de_intersect)
  de_intersect <- paste0(de_intersect,'_',ct)
  intersect_genes_de <- c(intersect_genes_de,de_intersect)
  total_intersect_de <- total_intersect_de + ilen_de
  scITD_res <- scITD_genes_only[scITD_genes_ct==ct]
  scITD_intersect <- intersect(scITD_res,iv_res)
  ilen_scITD <- length(scITD_intersect)
  scITD_intersect <- paste0(scITD_intersect,'_',ct)
  intersect_genes_scITD <- c(intersect_genes_scITD,scITD_intersect)
  total_intersect_scITD <- total_intersect_scITD + ilen_scITD
  num_sig_in_dat <- num_sig_in_dat + length(iv_res)
  num_tested_in_dat <- num_tested_in_dat + length(intersect(ct_test_full[[ct]]$featureData.symbol,bg))
}
print(total_intersect_de)
print(total_intersect_scITD)
tmp <- cbind.data.frame(c(total_intersect_de,total_intersect_scITD),c('DE','scITD_F1'))
colnames(tmp) <- c('count','method')
tmp$count <- as.numeric(tmp$count)
p <- ggplot(tmp,aes(x=method,y=count)) +
  geom_bar(stat = 'identity') +
  ggtitle('Overlap with genes in IFN stimulation DE') +
  ylab(paste0('# genes (out of ',as.character(ncol(de_cormat)),')')) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/de_f1_IFN_overlap.pdf", useDingbats = FALSE,
#     width = 4, height = 2.25)
p
dev.off()


### hypergeometric test enrichment pvals for in vitro ifn data gene set
num_bg <- length(de_genes[!is.na(de_genes)]) # population size
gset_sig <- num_sig_in_dat # number of success in population
sample_size <- ncol(de_cormat) # sample size
total_intersect_scITD # number of successes
total_intersect_de # number of successes
pval_scITD <- phyper(total_intersect_scITD-1, sample_size, num_bg-sample_size, gset_sig, lower.tail=FALSE)
pval_de <- phyper(total_intersect_de-1, sample_size, num_bg-sample_size, gset_sig, lower.tail=FALSE)
# pval scITD = 1.9e-118
# pval de = .0016


### plotting hmaps for de and scITD with annotations for
# gene sets, clusters, cell types, and significance in IFN DE
genes_lab_de1 <- go_bp_res_de[go_bp_res_de$Description=='cellular response to type I interferon','geneID']
genes_lab_de2 <- go_bp_res_de[go_bp_res_de$Description=="catecholamine secretion",'geneID']
genes_lab_de3 <- go_bp_res_de[go_bp_res_de$Description=="chromosome segregation",'geneID']
# genes_lab_de3 <- go_bp_res_de[go_bp_res_de$Description=="regulation of chromosome segregation",'geneID']
genes_lab_de4 <- go_bp_res_de[go_bp_res_de$Description=="organic hydroxy compound metabolic process",'geneID']
genes_lab_de5 <- go_mf_res_de[go_mf_res_de$Description=="exonuclease activity, active with either ribo- or deoxyribonucleic acids and producing 5'-phosphomonoesters" ,'geneID']

genes_lab_de1 <- unlist(strsplit(genes_lab_de1,split='/'))
genes_lab_de2 <- unlist(strsplit(genes_lab_de2,split='/'))
genes_lab_de3 <- unlist(strsplit(genes_lab_de3,split='/'))
genes_lab_de4 <- unlist(strsplit(genes_lab_de4,split='/'))
genes_lab_de5 <- unlist(strsplit(genes_lab_de5,split='/'))

genes_lab_de <- cbind.data.frame(c(de_genes_id %in% intersect_genes_de),
                                 c(de_genes_only %in% genes_lab_de1),
                                 c(de_genes_only %in% genes_lab_de2),
                                 c(de_genes_only %in% genes_lab_de3),
                                 c(de_genes_only %in% genes_lab_de4),
                                 c(de_genes_only %in% genes_lab_de5))
rownames(genes_lab_de) <- de_genes_id
colnames(genes_lab_de) <- c('IFN_stim_DE','GO_IFN','GO_catecholamine','GO_chrom_seg','GO_hydroxy_metab','GO_exonuclease')
genes_lab_de$IFN_stim_DE <- factor(genes_lab_de$IFN_stim_DE,levels=c(TRUE,FALSE))
genes_lab_de$GO_IFN <- factor(genes_lab_de$GO_IFN,levels=c(TRUE,FALSE))
genes_lab_de$GO_catecholamine <- factor(genes_lab_de$GO_catecholamine,levels=c(TRUE,FALSE))
genes_lab_de$GO_chrom_seg <- factor(genes_lab_de$GO_chrom_seg,levels=c(TRUE,FALSE))
genes_lab_de$GO_hydroxy_metab <- factor(genes_lab_de$GO_hydroxy_metab,levels=c(TRUE,FALSE))
genes_lab_de$GO_exonuclease <- factor(genes_lab_de$GO_exonuclease,levels=c(TRUE,FALSE))
levels(genes_lab_de$IFN_stim_DE) <- c('in_set','not_in_set')
levels(genes_lab_de$GO_IFN) <- c('in_set','not_in_set')
levels(genes_lab_de$GO_catecholamine) <- c('in_set','not_in_set')
levels(genes_lab_de$GO_chrom_seg) <- c('in_set','not_in_set')
levels(genes_lab_de$GO_hydroxy_metab) <- c('in_set','not_in_set')
levels(genes_lab_de$GO_exonuclease) <- c('in_set','not_in_set')

# ordering the genes the same as in the preclustered cor matrix
genes_lab_de <- genes_lab_de[rownames(de_cormat2),]

de_col <- c('red','gray95')
names(de_col) <- c('in_set','not_in_set')
ca_de <- ComplexHeatmap::HeatmapAnnotation(IFN_stim_DE = genes_lab_de$IFN_stim_DE,
                                           GO_IFN = genes_lab_de$GO_IFN,
                                           GO_catecholamine = genes_lab_de$GO_catecholamine,
                                           GO_chrom_seg = genes_lab_de$GO_chrom_seg,
                                           GO_hydroxy_metab = genes_lab_de$GO_hydroxy_metab,
                                           GO_exonuclease = genes_lab_de$GO_exonuclease,
                                           col = list(IFN_stim_DE = de_col,GO_IFN = de_col,GO_catecholamine = de_col,GO_chrom_seg = de_col,GO_hydroxy_metab = de_col,GO_exonuclease = de_col),
                                           show_legend = FALSE,show_annotation_name = TRUE,
                                           annotation_name_side = 'left')

col_fun = colorRamp2(c(0, 1), c("white", "red"))
hmap_de <- Heatmap(de_cormat2,name = "pearson r",
                   cluster_columns = FALSE,
                   col = col_fun,
                   cluster_rows = FALSE,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   show_row_dend = FALSE,
                   show_column_dend = FALSE,
                   bottom_annotation = ca_de,
                   right_annotation=ra_de,
                   top_annotation=ca_grp,
                   border = TRUE,
                   column_title = 'SLEDAI DE genes',
                   column_title_gp = gpar(fontsize = 14),
                   use_raster = TRUE)

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/DE_hmap_gset3.pdf", useDingbats = FALSE,
#     width = 7, height = 4)
hmap_de
dev.off()


# now for scITD result
genes_lab_sc1 <- go_bp_res_f1[go_bp_res_f1$Description=='cellular response to type I interferon','geneID']
genes_lab_sc2 <- go_bp_res_f1[go_bp_res_f1$Description=="catecholamine secretion",'geneID']
genes_lab_sc3 <- go_bp_res_f1[go_bp_res_f1$Description=="chromosome segregation",'geneID']
genes_lab_sc4 <- go_bp_res_f1[go_bp_res_f1$Description=="organic hydroxy compound metabolic process",'geneID']
genes_lab_sc5 <- go_mf_res_f1[go_mf_res_f1$Description=="exonuclease activity, active with either ribo- or deoxyribonucleic acids and producing 5'-phosphomonoesters" ,'geneID']

genes_lab_sc1 <- unlist(strsplit(genes_lab_sc1,split='/'))
genes_lab_sc2 <- unlist(strsplit(genes_lab_sc2,split='/'))
genes_lab_sc3 <- unlist(strsplit(genes_lab_sc3,split='/'))
genes_lab_sc4 <- unlist(strsplit(genes_lab_sc4,split='/'))
genes_lab_sc5 <- unlist(strsplit(genes_lab_sc5,split='/'))

genes_lab_sc <- cbind.data.frame(c(scITD_genes_id %in% intersect_genes_scITD),
                                 c(scITD_genes_only %in% genes_lab_sc1),
                                 c(scITD_genes_only %in% genes_lab_sc2),
                                 c(scITD_genes_only %in% genes_lab_sc3),
                                 c(scITD_genes_only %in% genes_lab_sc4),
                                 c(scITD_genes_only %in% genes_lab_sc5))
rownames(genes_lab_sc) <- scITD_genes_id
colnames(genes_lab_sc) <- c('IFN_stim_DE','GO_IFN','GO_catecholamine','GO_chrom_seg','GO_hydroxy_metab','GO_exonuclease')
genes_lab_sc$IFN_stim_DE <- factor(genes_lab_sc$IFN_stim_DE,levels=c(TRUE,FALSE))
genes_lab_sc$GO_IFN <- factor(genes_lab_sc$GO_IFN,levels=c(TRUE,FALSE))
genes_lab_sc$GO_catecholamine <- factor(genes_lab_sc$GO_catecholamine,levels=c(TRUE,FALSE))
genes_lab_sc$GO_chrom_seg <- factor(genes_lab_sc$GO_chrom_seg,levels=c(TRUE,FALSE))
genes_lab_sc$GO_hydroxy_metab <- factor(genes_lab_sc$GO_hydroxy_metab,levels=c(TRUE,FALSE))
genes_lab_sc$GO_exonuclease <- factor(genes_lab_sc$GO_exonuclease,levels=c(TRUE,FALSE))
levels(genes_lab_sc$IFN_stim_DE) <- c('in_set','not_in_set')
levels(genes_lab_sc$GO_IFN) <- c('in_set','not_in_set')
levels(genes_lab_sc$GO_catecholamine) <- c('in_set','not_in_set')
levels(genes_lab_sc$GO_chrom_seg) <- c('in_set','not_in_set')
levels(genes_lab_sc$GO_hydroxy_metab) <- c('in_set','not_in_set')
levels(genes_lab_sc$GO_exonuclease) <- c('in_set','not_in_set')

# ordering the genes the same as in the preclustered cor matrix
genes_lab_sc <- genes_lab_sc[rownames(scITD_cormat),]

de_col <- c('red','gray95')
names(de_col) <- c('in_set','not_in_set')
ca_scITD <- ComplexHeatmap::HeatmapAnnotation(IFN_stim_DE = genes_lab_sc$IFN_stim_DE,
                                              GO_IFN = genes_lab_sc$GO_IFN,
                                              GO_catecholamine = genes_lab_sc$GO_catecholamine,
                                              GO_chrom_seg = genes_lab_sc$GO_chrom_seg,
                                              GO_hydroxy_metab = genes_lab_sc$GO_hydroxy_metab,
                                              GO_exonuclease = genes_lab_sc$GO_exonuclease,
                                           col = list(IFN_stim_DE = de_col,GO_IFN = de_col,GO_catecholamine = de_col,GO_chrom_seg = de_col,GO_hydroxy_metab = de_col,GO_exonuclease = de_col),
                                           show_legend = FALSE,show_annotation_name = TRUE,
                                           annotation_name_side = 'left')

g_ct_cev_f1 <- sapply(colnames(scITD_cormat),function(x){
  return(strsplit(x,split='_')[[1]][[2]])
})

ra_f1 <- ComplexHeatmap::rowAnnotation(ctype = g_ct_cev_f1, show_annotation_name=FALSE, 
                                       col = list(ctype = annot_col))
clust_method <- 'ward.D2'
scITD_hmap <- Heatmap(scITD_cormat,name = "pearson r",
                      col = col_fun,
                      clustering_method_rows = clust_method,
                      clustering_method_columns = clust_method,
                      cluster_columns = TRUE,
                      cluster_rows = TRUE,
                      show_row_names = FALSE,
                      show_column_names = FALSE,
                      show_row_dend = FALSE,
                      show_column_dend = FALSE,
                      right_annotation=ra_f1,
                      bottom_annotation = ca_scITD,
                      column_title = 'scITD F1 significant genes',
                      column_title_gp = gpar(fontsize = 14),
                      border=TRUE,
                      use_raster = TRUE)




# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/F1_cor_hmap_gset3.pdf", useDingbats = FALSE,
#     width = 6, height = 4)
scITD_hmap
dev.off()








### comparing gene pair correlations between de and f1 for matched number of hits
diag(scITD_cormat) <- NA
tmp1 <- as.data.frame(c(scITD_cormat))
tmp1$method <- 'scITD'
colnames(tmp1)[1] <- 'pearson_r'
diag(de_cormat) <- NA
tmp2 <- as.data.frame(c(de_cormat))
tmp2$method <- 'DE'
colnames(tmp2)[1] <- 'pearson_r'
tmp <- rbind.data.frame(tmp1,tmp2)

p <- ggplot(tmp, aes(x=pearson_r, color=method)) +
  geom_histogram(fill="white", position="dodge",bins = 100) +
  # geom_histogram(data=tmp[tmp$method=='DE',],fill="white", position="dodge",bins = 100,aes(y = after_stat(count / sum(count)))) +
  # geom_histogram(data=tmp[tmp$method=='scITD',],fill="white", position="dodge",bins = 100,aes(y = after_stat(count / sum(count)))) +
  theme_classic()

pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/pairwise_cor_de_f1.pdf", useDingbats = FALSE,
    width = 5, height = 2.5)
p
dev.off()








## run DE also for prednisone
## using meds variables
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_meds_cleaned.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars[,'Sample ID']
clin_vars[,'Sample ID'] <- NULL
clin_vars[is.na(clin_vars)] <- 0
##

# get tucker donor scores to test
dsc <- pbmc_container$tucker_results[[1]]

## get donors IDs matching to the clinical variable IDs
trim_names <- sapply(rownames(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})

# get donors in both dataframes
d_both <- rownames(clin_vars)[rownames(clin_vars) %in% trim_names]
trim_names <- trim_names[trim_names %in% d_both]

##### running DE analysis with the same data that's in the tensor used in scITD
de_genes <- c()
meta_var <- 'prednisone'
for (ct in pbmc_container$experiment_params$ctypes_use) {
  print(ct)
  # getting out pseudobulk for a cell type
  pb <- pbmc_container$scMinimal_ctype[[ct]]$pseudobulk
  rownames(pb) <- trim_names[rownames(pb)]
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  
  # run DE t-test
  ct_res <- plapply(1:ncol(pb),function(j) {
    tmp <- cbind.data.frame(pb[,j],clin_vars[rownames(pb),meta_var])
    colnames(tmp) <- c('expr','meta')
    g1 <- tmp$expr[tmp$meta==0]
    g2 <- tmp$expr[tmp$meta==1]
    tres <- wilcox.test(g1,g2)
    tres <- t.test(g1,g2)
    tres_pval <- tres$p.value
    names(tres_pval) <- colnames(pb)[j]
    return(tres_pval)
  },mc.preschedule=TRUE,n.cores=20,progress=TRUE)
  
  de_genes <- c(de_genes,ct_res)
}
de_genes <- unlist(de_genes)

# correct pvals
pv_thresh <- .1
de_genes_adj <- p.adjust(de_genes,method = 'fdr')
de_genes_adj <- de_genes_adj[!is.na(de_genes_adj)]
de_genes_id <- names(de_genes_adj)[de_genes_adj<pv_thresh]
print(length(de_genes_id))


f3_pvals <- get_one_factor_gene_pvals(pbmc_container,3)
f3_pvals_vec <- c()
for (j in 1:ncol(f3_pvals)) {
  ct_vec <- f3_pvals[,j]
  names(ct_vec) <- paste0(rownames(f3_pvals),'_',colnames(f3_pvals)[j])
  f3_pvals_vec <- c(f3_pvals_vec,ct_vec)
}

g_both <- intersect(names(de_genes_adj),names(f3_pvals_vec))
tmp <- cbind.data.frame(de_genes_adj[g_both],f3_pvals_vec[g_both])
colnames(tmp) <- c('de_padj','f3_padj')
tmp$de_padj <- -log10(tmp$de_padj)
tmp$f3_padj <- -log10(tmp$f3_padj)
tmp$ctype <- sapply(rownames(tmp),function(x){
  strsplit(x,split='_')[[1]][[2]]
})
tmp$gene <- sapply(rownames(tmp),function(x){
  strsplit(x,split='_')[[1]][[1]]
})

m_df <- data.frame()
m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                    category = "C5", subcategory = "BP"))
my_pathways <- split(m_df$gene_symbol, f = m_df$gs_name)

tmp_sub <- tmp[tmp$ctype=='cMono',]
mystats <- tmp_sub$de_padj
names(mystats) <- tmp_sub$gene
fgsea_res_de <- fgsea::fgsea(pathways = my_pathways,
                               stats = mystats,
                               minSize=10,
                               maxSize=500,
                               eps=0,
                               gseaParam=0,
                               scoreType = "pos")

mystats <- tmp_sub$f3_padj
names(mystats) <- tmp_sub$gene
fgsea_res_f3 <- fgsea::fgsea(pathways = my_pathways,
                               stats = mystats,
                               minSize=10,
                               maxSize=500,
                               eps=0,
                               gseaParam=0,
                               scoreType = "pos")

## GOBP_RESPONSE_TO_HORMONE pvals
# de pval=.026
# f3 pval=.00036


tmp$gr_gset <- tmp$gene %in% my_pathways[['GOBP_RESPONSE_TO_HORMONE']]
g_to_label <- fgsea_res_f3[fgsea_res_f3$pathway=='GOBP_RESPONSE_TO_HORMONE','leadingEdge'][[1]][[1]][1:15]
tmp$glab <- NA
tmp[tmp$gene %in% g_to_label & tmp$ctype=='cMono','glab'] <- tmp$gene[tmp$gene %in% g_to_label & tmp$ctype=='cMono']
myColors <- c('grey80','red4')
p <- ggplot(tmp,aes(x=f3_padj,y=de_padj,color=gr_gset,label=glab)) +
  geom_point_rast(alpha=.6) +
  geom_text_repel(size=3) +
  geom_hline(yintercept=-log10(.01),color='red',linetype='dashed') +
  geom_vline(xintercept=-log10(.01),color='red',linetype='dashed') +
  scale_colour_manual(name = "grp",values = myColors) +
  xlab('scITD F3 -log10(padj)') +
  ylab('Prednisone DE -log10(padj)') +
  theme_bw()
  
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/pred_de_vs_f3_pvals.pdf", useDingbats = FALSE,
#     width = 6, height = 4)
p
dev.off()



### making a heatmap of expression to show why many of these genes not de
dsc <- dsc[names(trim_names),]
rownames(dsc) <- trim_names

# g_to_label2 <- c('FKBP5', 'ZBTB16', 'TNFAIP3', 'PDK4', 'TSC22D3')
# g_to_label <- fgsea_res_f3[fgsea_res_f3$pathway=='GOBP_RESPONSE_TO_HORMONE','leadingEdge'][[1]][[1]][1:15]
# g_to_label <- fgsea_res_f3[fgsea_res_f3$pathway=='GOBP_CELLULAR_RESPONSE_TO_CORTICOSTEROID_STIMULUS','leadingEdge'][[1]][[1]]
pb <- pbmc_container$scMinimal_ctype[['cMono']]$pseudobulk
pb <- pb[names(trim_names),g_to_label]
rownames(pb) <- trim_names
clin_vars_use <- clin_vars[trim_names,'prednisone']
pb <- pb[clin_vars_use==1,]
f3_dsc <- dsc[rownames(pb),3] 
pb <- t(pb)
ndx_order <- order(f3_dsc,decreasing = FALSE)
pb <- pb[,ndx_order]
f3_dsc <- f3_dsc[ndx_order]

col_fun_dsc = colorRamp2(c(min(f3_dsc),0, max(f3_dsc)), c("darkgreen","white", "darkorange"))
ca_dsc <- ComplexHeatmap::HeatmapAnnotation(f3_dsc = f3_dsc,
                                            col = list(f3_dsc = col_fun_dsc),
                                            show_legend = TRUE)

col_fun_expr = colorRamp2(c(min(pb),0, max(pb)), c("darkblue","white", "darkred"))
hmap_expr <- Heatmap(pb,name = "expression",
                   cluster_columns = FALSE,
                   col = col_fun_expr,
                   cluster_rows = TRUE,
                   show_row_names = TRUE,
                   show_column_names = FALSE,
                   show_row_dend = FALSE,
                   show_column_dend = FALSE,
                   top_annotation=ca_dsc,
                   border = TRUE,
                   column_title = 'Hormone response gene expression\nin prednisone patients',
                   column_title_gp = gpar(fontsize = 14))

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/pred_hormone_expr2.pdf", useDingbats = FALSE,
#     width = 5.5, height = 4)
hmap_expr
dev.off()



### plotting expression of these genes for prednisone positive patients vs not
g_to_label <- c('KLF9','IRS2','HMGB2','ZFP36L2')
pb <- pbmc_container$scMinimal_ctype[['cMono']]$pseudobulk
pb <- pb[names(trim_names),g_to_label]
rownames(pb) <- trim_names
melted_pb <- melt(pb, id.vars = "Gene", variable.name = "Donor", value.name = "Expression")
colnames(melted_pb) <- c('donor','gene','expression')
melted_pb$donor <- as.character(melted_pb$donor)
melted_pb$pred <- clin_vars[melted_pb$donor,'prednisone']
melted_pb$pred <- factor(melted_pb$pred,levels=c(1,0))
levels(melted_pb$pred) <- c('prednisone','no prednisone')
p <- ggplot(melted_pb,aes(x=gene,y=expression,color=pred)) +
  geom_beeswarm(dodge.width=.75,cex=2,corral = 'wrap',corral.width = .35) +
  theme_bw()

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/pred_group_expr.pdf", useDingbats = FALSE,
#     width = 7, height = 4)
p
dev.off()











