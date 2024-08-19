
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

# # or load up the already generated object
# pbmc_container <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_container_w_decomp.rds')


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
pbmc_container_full$projection_data <- pbmc_container$projection_data


# get significant genes
pbmc_container_full <- get_lm_pvals(pbmc_container_full)

## using ordinal variables
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_ordinal.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL

# get tucker donor scores to test
dsc <- pbmc_container_full$tucker_results[[1]]

## get donors IDs matching to the clinical variable IDs
trim_names <- sapply(rownames(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})

# get donors in both dataframes
d_both <- rownames(clin_vars)[rownames(clin_vars) %in% trim_names]
trim_names <- trim_names[trim_names %in% d_both]

de_genes_sledai <- c()
de_genes_sledai_z <- c()
meta_var <- 'sledaiscore'
for (ct in pbmc_container_full$experiment_params$ctypes_use) {
  print(ct)
  # getting out pseudobulk for a cell type
  pb <- pbmc_container_full$scMinimal_ctype[[ct]]$pseudobulk
  rownames(pb) <- trim_names[rownames(pb)]
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  
  # run DE test
  ct_res <- plapply(1:ncol(pb),function(j) {
    tmp <- cbind.data.frame(pb[,j],clin_vars[rownames(pb),meta_var])
    colnames(tmp) <- c('expr','meta')
    # with linear model
    lmres <- summary(lm(meta~expr,data=tmp))
    if (all(tmp$expr==0)) {
      tres_pval <- NaN
      tres_z <- NaN
    } else {
      tres_pval <- lmres$coefficients['expr','Pr(>|t|)']
      tres_z <- lmres$coefficients['expr','t value']
    }
    names(tres_pval) <- colnames(pb)[j]
    names(tres_z) <- colnames(pb)[j]
    return(list(tres_pval,tres_z))
  },mc.preschedule=TRUE,n.cores=20,progress=TRUE)
  
  ct_res_pv <- lapply(1:length(ct_res),function(x){
    return(ct_res[[x]][[1]])
  })
  ct_res_z <- lapply(1:length(ct_res),function(x){
    return(ct_res[[x]][[2]])
  })
  
  de_genes_sledai <- c(de_genes_sledai,ct_res_pv)
  de_genes_sledai_z <- c(de_genes_sledai_z,ct_res_z)
  
}
de_genes_sledai <- unlist(de_genes_sledai)
de_genes_sledai_z <- unlist(de_genes_sledai_z)

# correct pvals
pv_thresh <- .1
de_genes_sledai_adj <- p.adjust(de_genes_sledai,method = 'fdr')
de_genes_sledai_adj <- de_genes_sledai_adj[!is.na(de_genes_sledai_adj)]
de_genes_sledai_id <- names(de_genes_sledai_adj)[de_genes_sledai_adj<pv_thresh]
print(length(de_genes_sledai_id))

## generate correlation matrix for DE significant genes
# first need to append expression for all cell types
pb_all <- matrix(nrow=length(trim_names),ncol=0)
for (ct in pbmc_container_full$experiment_params$ctypes_use) {
  # getting out pseudobulk for a cell type
  pb <- pbmc_container_full$scMinimal_ctype[[ct]]$pseudobulk
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  g_both <- intersect(colnames(pb),de_genes_sledai_id)
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

bg <- names(de_genes_sledai)[!is.na(de_genes_sledai)]
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
f1_pvals <- get_one_factor_gene_pvals(pbmc_container_full,1)

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
for (ct in pbmc_container_full$experiment_params$ctypes_use) {
  # getting out pseudobulk for a cell type
  pb <- pbmc_container_full$scMinimal_ctype[[ct]]$pseudobulk
  pb <- pb[names(trim_names),rownames(f1_pvals2)[f1_pvals2[,ct]<pv_thresh]]
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  # append pb to existing pb
  pb_all <- cbind(pb_all,pb)
}

scITD_cormat <- abs(cor(pb_all,method='pearson'))



### running go for significant genes by either method
de_genes_ct <- sapply(de_genes_sledai_id,function(x) {
  strsplit(x,split='_')[[1]][[2]]
})
de_genes_only <- sapply(de_genes_sledai_id,function(x) {
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

### Figure S3j
pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/de_f1_IFN_overlap.pdf", useDingbats = FALSE,
    width = 4, height = 3)
p
dev.off()


### hypergeometric test enrichment pvals for in vitro ifn data gene set
num_bg <- length(de_genes_sledai[!is.na(de_genes_sledai)]) # population size
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

genes_lab_de <- cbind.data.frame(c(de_genes_sledai_id %in% intersect_genes_de),
                                 c(de_genes_only %in% genes_lab_de1),
                                 c(de_genes_only %in% genes_lab_de2),
                                 c(de_genes_only %in% genes_lab_de3),
                                 c(de_genes_only %in% genes_lab_de4),
                                 c(de_genes_only %in% genes_lab_de5))
rownames(genes_lab_de) <- de_genes_sledai_id
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

## Figure S3g
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



## Figure S3h
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
  theme_classic()

## Figure S3i
pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/pairwise_cor_de_f1.pdf", useDingbats = FALSE,
    width = 5, height = 4)
p
dev.off()



## making a version of the histogram but only with pairs of different genes from different cell types
scITD_cormat2 <- scITD_cormat
de_cormat2 <- de_cormat

for (i in 1:nrow(scITD_cormat2)) {
  g_ct1 <- strsplit(rownames(scITD_cormat2)[i],split='_')[[1]]
  g1 <- g_ct1[[1]]
  ct1 <- g_ct1[[2]]
  for (j in 1:ncol(scITD_cormat2)) {
    g_ct2 <- strsplit(colnames(scITD_cormat2)[j],split='_')[[1]]
    g2 <- g_ct2[[1]]
    ct2 <- g_ct2[[2]]
    if (ct1==ct2 | g1==g2) {
      scITD_cormat2[i,j] <- NA
    }
  }
}

for (i in 1:nrow(de_cormat2)) {
  g_ct1 <- strsplit(rownames(de_cormat2)[i],split='_')[[1]]
  g1 <- g_ct1[[1]]
  ct1 <- g_ct1[[2]]
  for (j in 1:ncol(de_cormat2)) {
    g_ct2 <- strsplit(colnames(de_cormat2)[j],split='_')[[1]]
    g2 <- g_ct2[[1]]
    ct2 <- g_ct2[[2]]
    if (ct1==ct2 | g1==g2) {
      de_cormat2[i,j] <- NA
    }
  }
}


diag(scITD_cormat2) <- NA
tmp1 <- as.data.frame(c(scITD_cormat2))
tmp1$method <- 'scITD'
colnames(tmp1)[1] <- 'pearson_r'
diag(de_cormat2) <- NA
tmp2 <- as.data.frame(c(de_cormat2))
tmp2$method <- 'DE'
colnames(tmp2)[1] <- 'pearson_r'
tmp <- rbind.data.frame(tmp1,tmp2)

p <- ggplot(tmp, aes(x=pearson_r, color=method)) +
  geom_histogram(fill="white", position="dodge",bins = 100) +
  theme_classic()

pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/pairwise_cor_de_f1_diff_only.pdf", useDingbats = FALSE,
    width = 5, height = 2.5)
p
dev.off()











## loading categorical variables to test nephritis
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_categorical.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL

# get tucker donor scores to test
dsc <- pbmc_container_full$tucker_results[[1]]

## get donors IDs matching to the clinical variable IDs
trim_names <- sapply(rownames(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})

# get donors in both dataframes
d_both <- rownames(clin_vars)[rownames(clin_vars) %in% trim_names]
trim_names <- trim_names[trim_names %in% d_both]

##### running DE analysis with the same data that's in the tensor used in scITD
de_genes_neph <- c()
de_genes_neph_z <- c()
meta_var <- 'crflupusneph'
for (ct in pbmc_container_full$experiment_params$ctypes_use) {
  print(ct)
  # getting out pseudobulk for a cell type
  pb <- pbmc_container_full$scMinimal_ctype[[ct]]$pseudobulk
  rownames(pb) <- trim_names[rownames(pb)]
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  
  # run DE t-test
  ct_res <- plapply(1:ncol(pb),function(j) {
    tmp <- cbind.data.frame(pb[,j],clin_vars[rownames(pb),meta_var])
    colnames(tmp) <- c('expr','meta')
    g1 <- tmp$expr[tmp$meta==0]
    g2 <- tmp$expr[tmp$meta==1]
    tres <- t.test(g1,g2)
    tres_pval <- tres$p.value
    tres_z <- tres[["statistic"]][["t"]]

    names(tres_pval) <- colnames(pb)[j]
    names(tres_z) <- colnames(pb)[j]
    
    return(list(tres_pval,tres_z))
  },mc.preschedule=TRUE,n.cores=20,progress=TRUE)
  
  ct_res_pv <- lapply(1:length(ct_res),function(x){
    return(ct_res[[x]][[1]])
  })
  ct_res_z <- lapply(1:length(ct_res),function(x){
    return(ct_res[[x]][[2]])
  })
  
  de_genes_neph <- c(de_genes_neph,ct_res_pv)
  de_genes_neph_z <- c(de_genes_neph_z,ct_res_z)
  
}
de_genes_neph <- unlist(de_genes_neph)
de_genes_neph_z <- unlist(de_genes_neph_z)


# correct pvals
pv_thresh <- .1
de_genes_neph_adj <- p.adjust(de_genes_neph,method = 'fdr')
de_genes_neph_adj <- de_genes_neph_adj[!is.na(de_genes_neph_adj)]
de_genes_neph_id <- names(de_genes_neph_adj)[de_genes_neph_adj<pv_thresh]
print(length(de_genes_neph_id))






## run DE also for prednisone
## using meds variables
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_meds_cleaned.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars[,'Sample ID']
clin_vars[,'Sample ID'] <- NULL
clin_vars[is.na(clin_vars)] <- 0
##

# get tucker donor scores to test
dsc <- pbmc_container_full$tucker_results[[1]]

## get donors IDs matching to the clinical variable IDs
trim_names <- sapply(rownames(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})

# get donors in both dataframes
d_both <- rownames(clin_vars)[rownames(clin_vars) %in% trim_names]
trim_names <- trim_names[trim_names %in% d_both]

##### running DE analysis with the same data that's in the tensor used in scITD
de_genes_pred <- c()
de_genes_pred_z <- c()
meta_var <- 'prednisone'
for (ct in pbmc_container_full$experiment_params$ctypes_use) {
  print(ct)
  # getting out pseudobulk for a cell type
  pb <- pbmc_container_full$scMinimal_ctype[[ct]]$pseudobulk
  rownames(pb) <- trim_names[rownames(pb)]
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  
  # run DE t-test
  ct_res <- plapply(1:ncol(pb),function(j) {
    tmp <- cbind.data.frame(pb[,j],clin_vars[rownames(pb),meta_var])
    colnames(tmp) <- c('expr','meta')
    g1 <- tmp$expr[tmp$meta==0]
    g2 <- tmp$expr[tmp$meta==1]
    tres <- t.test(g1,g2)
    tres_pval <- tres$p.value
    tres_z <- tres[["statistic"]][["t"]]

    names(tres_pval) <- colnames(pb)[j]
    names(tres_z) <- colnames(pb)[j]
    
    return(list(tres_pval,tres_z))
  },mc.preschedule=TRUE,n.cores=20,progress=TRUE)
  
  ct_res_pv <- lapply(1:length(ct_res),function(x){
    return(ct_res[[x]][[1]])
  })
  ct_res_z <- lapply(1:length(ct_res),function(x){
    return(ct_res[[x]][[2]])
  })
  
  de_genes_pred <- c(de_genes_pred,ct_res_pv)
  de_genes_pred_z <- c(de_genes_pred_z,ct_res_z)
  
}
de_genes_pred <- unlist(de_genes_pred)
de_genes_pred_z <- unlist(de_genes_pred_z)


# correct pvals
pv_thresh <- .1
de_genes_pred_adj <- p.adjust(de_genes_pred,method = 'fdr')
de_genes_pred_adj <- de_genes_pred_adj[!is.na(de_genes_pred_adj)]
de_genes_pred_id <- names(de_genes_pred_adj)[de_genes_pred_adj<pv_thresh]
print(length(de_genes_pred_id))


f3_pvals <- get_one_factor_gene_pvals(pbmc_container_full,3)
f3_pvals_vec <- c()
for (j in 1:ncol(f3_pvals)) {
  ct_vec <- f3_pvals[,j]
  names(ct_vec) <- paste0(rownames(f3_pvals),'_',colnames(f3_pvals)[j])
  f3_pvals_vec <- c(f3_pvals_vec,ct_vec)
}

g_both <- intersect(names(de_genes_pred_adj),names(f3_pvals_vec))
tmp <- cbind.data.frame(de_genes_pred_adj[g_both],f3_pvals_vec[g_both])
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
  
## Figure S5f
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
g_to_label <- g_to_label[!(g_to_label %in% 'ZFP36L1')] # just removing it because it's the only one going the opposite direction, not an issue though
pb <- pbmc_container_full$scMinimal_ctype[['cMono']]$pseudobulk
pb <- pb[names(trim_names),g_to_label]
rownames(pb) <- trim_names
clin_vars_use <- clin_vars[trim_names,'prednisone']
names(clin_vars_use) <- trim_names
# pb <- pb[clin_vars_use==1,]
f3_dsc <- dsc[rownames(pb),3] 
pb <- t(pb)
ndx_order <- order(f3_dsc,decreasing = FALSE)
pb <- pb[,ndx_order]
f3_dsc <- f3_dsc[ndx_order]
clin_vars_use <- clin_vars_use[ndx_order]

col_fun_dsc = colorRamp2(c(min(f3_dsc),0, max(f3_dsc)), c("darkgreen","white", "darkorange"))
ca_dsc <- ComplexHeatmap::HeatmapAnnotation(f3_dsc = f3_dsc,
                                            col = list(f3_dsc = col_fun_dsc),
                                            show_legend = TRUE)

clin_vars_use <- factor(clin_vars_use,levels=c(1,0))
levels(clin_vars_use) <- c('Prednisone','No prednisone')
col_fun_pred <- c('yellow','green')
names(col_fun_pred) <- c('Prednisone','No prednisone')
pred_use <- ComplexHeatmap::HeatmapAnnotation(pred_use = clin_vars_use,
                                            col = list(pred_use = col_fun_pred),
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
                   bottom_annotation = pred_use,
                   border = TRUE,
                   column_title = 'Hormone response gene expression',
                   column_title_gp = gpar(fontsize = 14))

## Figure S5h
pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/pred_hormone_expr3.pdf", useDingbats = FALSE,
    width = 5.5, height = 3.5)
hmap_expr
dev.off()



### plotting expression of these genes for prednisone positive patients vs not
g_to_label <- c('KLF9','IRS2','HMGB2','ZFP36L2')
pb <- pbmc_container_full$scMinimal_ctype[['cMono']]$pseudobulk
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

## Figure S5g
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/pred_group_expr.pdf", useDingbats = FALSE,
#     width = 7, height = 4)
p
dev.off()























######### evaluating how SLEDAI DE vs aSLE_F1 effect sizes compare with those of the pediatric dataset
# load up the pediatric lupus dataset: see preprocessing/pediatric_lupus_preprocessing.R 
# for code used to generate this object
pbmc <- readRDS(file='/home/jmitchel/data/pediatric_lupus/processed/pbmc_pediatric_clean_annotated_seurat.rds')

### change types of some metadata
pbmc@meta.data$SLEDAI <- as.numeric(as.character(pbmc@meta.data$SLEDAI))
pbmc@meta.data$ESR <- as.numeric(as.character(pbmc@meta.data$ESR))
pbmc@meta.data$WBC <- as.numeric(as.character(pbmc@meta.data$WBC))
pbmc@meta.data$RBC <- as.numeric(as.character(pbmc@meta.data$RBC))
pbmc@meta.data$MONOCYTE_per <- as.numeric(as.character(pbmc@meta.data$MONOCYTE_per))
pbmc@meta.data$NEUTROPHIL_per <- as.numeric(as.character(pbmc@meta.data$NEUTROPHIL_per))
pbmc@meta.data$LYMPHOCYTE_per <- as.numeric(as.character(pbmc@meta.data$LYMPHOCYTE_per))
pbmc@meta.data$HGB <- as.numeric(as.character(pbmc@meta.data$HGB))
pbmc@meta.data$HCT <- as.numeric(as.character(pbmc@meta.data$HCT))
pbmc@meta.data$PLATELETS <- as.numeric(as.character(pbmc@meta.data$PLATELETS))
pbmc@meta.data$NEU_ABS <- as.numeric(as.character(pbmc@meta.data$NEU_ABS))
pbmc@meta.data$LYM_ABS <- as.numeric(as.character(pbmc@meta.data$LYM_ABS))
pbmc@meta.data$CREATININE <- as.numeric(as.character(pbmc@meta.data$CREATININE))
pbmc@meta.data$ALBUMIN <- as.numeric(as.character(pbmc@meta.data$ALBUMIN))
pbmc@meta.data$C3 <- as.numeric(as.character(pbmc@meta.data$C3))
pbmc@meta.data$C4 <- as.numeric(as.character(pbmc@meta.data$C4))
pbmc@meta.data$ALT <- as.numeric(as.character(pbmc@meta.data$ALT))
pbmc@meta.data$AST <- as.numeric(as.character(pbmc@meta.data$AST))
pbmc@meta.data$ALD <- as.numeric(as.character(pbmc@meta.data$ALD))
pbmc@meta.data$LDH <- as.numeric(as.character(pbmc@meta.data$LDH))
pbmc@meta.data$DSDNA_ratio <- NULL
# pbmc@meta.data$DSDNA <- NULL
levels(pbmc@meta.data$DSDNA)
levels(pbmc@meta.data$DSDNA) <- c('detected','detected','detected','ND','detected','none detected')

# removing clusters PC_pDC, Eryth, and Mgk
cells_keep <- rownames(pbmc@meta.data)[!(pbmc@meta.data$clusters_fine %in% c('PC_pDC','Eryth','Mgk'))]
pbmc <- subset(pbmc,cells=cells_keep)

meta_cols <- colnames(pbmc@meta.data)[10:55]

# need to reset factor levels for clustering and donors to be the unique values
pbmc@meta.data$donors <- factor(pbmc@meta.data$donors,levels=unique(pbmc@meta.data$donors))
pbmc@meta.data$clusters_coarse <- factor(pbmc@meta.data$clusters_coarse,levels=unique(pbmc@meta.data$clusters_coarse))

# making all metadata ND values to be NA values...
pbmc@meta.data[pbmc@meta.data=='ND'] <- NA

# set up project parameters
param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc",
                                               "cMono","ncMono"),
                                ncores = 30, rand_seed = 10)

pbmc_container_pediatric <- make_new_container(seurat_obj=pbmc,
                                     params=param_list,
                                     metadata_cols=c('donors',
                                                     'clusters_coarse',
                                                     "Batch",
                                                     "Age",
                                                     "Gender",
                                                     "Ethnicity",
                                                     "Groups",
                                                     meta_cols),
                                     metadata_col_nm=c('donors',
                                                       'ctypes',
                                                       'Batch',
                                                       'Age',
                                                       'sex',
                                                       'Ethnicity',
                                                       'SLE_status',
                                                       meta_cols))

pbmc_container_pediatric <- form_tensor(pbmc_container_pediatric, donor_min_cells=2,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='Batch')


pbmc_container_pediatric <- run_tucker_ica(pbmc_container_pediatric, ranks=c(3,10),
                                 tucker_type = 'regular', rotation_type = 'hybrid')

pbmc_container_pediatric <- get_meta_associations(pbmc_container_pediatric,vars_test=c('sex','SLE_status','Age','Batch','Ethnicity'),
                                        stat_use='pval')

# plot donor scores by status
pbmc_container_pediatric <- plot_donor_matrix(pbmc_container_pediatric,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

pbmc_container_pediatric$plots$donor_matrix



# project all original factors onto the pediatric dataset
pbmc_container_pediatric <- project_new_data(pbmc_container_pediatric,pbmc_container)


### comparing projected scores to the pediatric lupus decomposition dscores

# reversing the signs of factors so that they match up with those from the main sle analysis
pbmc_container_pediatric$tucker_results[[1]] <- pbmc_container_pediatric$tucker_results[[1]] * -1
pbmc_container_pediatric$tucker_results[[2]] <- pbmc_container_pediatric$tucker_results[[2]] * -1

cor_mat <- cor(pbmc_container_pediatric[["projected_scores"]],pbmc_container_pediatric$tucker_results[[1]])
colnames(cor_mat) <- sapply(c(1:ncol(cor_mat)),function(x){
  paste0('Factor ',x)
})
rownames(cor_mat) <- sapply(c(1:nrow(cor_mat)),function(x){
  paste0('Projected ',x)
})

# make heatmap of this
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
cor_hmap <- Heatmap(cor_mat, name = "Pearson r",
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    column_names_gp = gpar(fontsize = 10),
                    row_names_gp = gpar(fontsize = 10),
                    col = col_fun,border=TRUE, show_column_names=TRUE,
                    show_row_names=TRUE,show_row_dend = FALSE,
                    show_column_dend = FALSE,
                    column_title_side = "bottom",
                    row_names_side = "left",
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid::grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 10))
                    })

cor_hmap


# recompute tensor and gene associations using all genes
pbmc_container_full_pediatric <- make_new_container(seurat_obj=pbmc,
                                     params=param_list,
                                     metadata_cols=c('donors',
                                                     'clusters_coarse',
                                                     "Batch",
                                                     "Age",
                                                     "Gender",
                                                     "Ethnicity",
                                                     "Groups",
                                                     meta_cols),
                                     metadata_col_nm=c('donors',
                                                       'ctypes',
                                                       'Batch',
                                                       'Age',
                                                       'sex',
                                                       'Ethnicity',
                                                       'SLE_status',
                                                       meta_cols))

pbmc_container_full_pediatric <- form_tensor(pbmc_container_full_pediatric, donor_min_cells=2,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.5,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='Batch')

pbmc_container_full_pediatric$tucker_results <- pbmc_container_pediatric$tucker_results
pbmc_container_full_pediatric$projection_data <- pbmc_container_pediatric$projection_data

# get significant genes
pbmc_container_full_pediatric <- get_lm_pvals(pbmc_container_full_pediatric)


# now compute SLEDAI DE for pSLE dataset
pbmc_container_full_pediatric <- get_donor_meta(pbmc_container_full_pediatric, additional_meta = meta_cols)
pbmc_container_full_pediatric$donor_metadata[1:10,1:10]
table(pbmc_container_full_pediatric$donor_metadata$SLEDAI)

de_genes_sledai_pediatric <- c()
de_genes_z_sledai_pediatric <- c()
meta_var <- 'SLEDAI'
for (ct in pbmc_container_full_pediatric$experiment_params$ctypes_use) {
  print(ct)
  # getting out pseudobulk for a cell type
  pb <- pbmc_container_full_pediatric$scMinimal_ctype[[ct]]$pseudobulk
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  
  # run DE test
  ct_res <- plapply(1:ncol(pb),function(j) tryCatch({
    tmp <- cbind.data.frame(pb[,j],pbmc_container_full_pediatric$donor_metadata[rownames(pb),meta_var])
    colnames(tmp) <- c('expr','meta')
    tmp <- tmp[!is.na(tmp$meta),]
    # with linear model
    lmres <- summary(lm(meta~expr,data=tmp))
    if (all(tmp$expr==0) | all(tmp$expr==tmp$expr[1])) {
      tres_pval <- NaN
      tres_z <- NaN
    } else {
      tres_pval <- lmres$coefficients['expr','Pr(>|t|)']
      tres_z <- lmres$coefficients['expr','t value']
    }
    names(tres_pval) <- colnames(pb)[j]
    names(tres_z) <- colnames(pb)[j]
    return(list(tres_pval,tres_z))
  },error=function(e) paste0('error_index_',j)),progress = TRUE,n.cores = 20,mc.preschedule = TRUE)
  
  ct_res_pv <- lapply(1:length(ct_res),function(x){
    return(ct_res[[x]][[1]])
  })
  ct_res_z <- lapply(1:length(ct_res),function(x){
    return(ct_res[[x]][[2]])
  })
  
  de_genes_sledai_pediatric <- c(de_genes_sledai_pediatric,ct_res_pv)
  de_genes_z_sledai_pediatric <- c(de_genes_z_sledai_pediatric,ct_res_z)
  
}
de_genes_sledai_pediatric <- unlist(de_genes_sledai_pediatric)
de_genes_z_sledai_pediatric <- unlist(de_genes_z_sledai_pediatric)

# correct pvals
pv_thresh <- .1
de_genes_sledai_adj_pediatric <- p.adjust(de_genes_sledai_pediatric,method = 'fdr')
de_genes_sledai_adj_pediatric <- de_genes_sledai_adj_pediatric[!is.na(de_genes_sledai_adj_pediatric)]
de_genes_sledai_id_pediatric <- names(de_genes_sledai_adj_pediatric)[de_genes_sledai_adj_pediatric<pv_thresh]
print(length(de_genes_sledai_id_pediatric))


de_genes_neph_pediatric <- c()
de_genes_z_neph_pediatric <- c()
meta_var <- 'Neph_all'
for (ct in pbmc_container_full_pediatric$experiment_params$ctypes_use) {
  print(ct)
  # getting out pseudobulk for a cell type
  pb <- pbmc_container_full_pediatric$scMinimal_ctype[[ct]]$pseudobulk
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  
  # run DE t-test
  ct_res <- plapply(1:ncol(pb),function(j) {
    tmp <- cbind.data.frame(pb[,j],pbmc_container_full_pediatric$donor_metadata[rownames(pb),meta_var])
    colnames(tmp) <- c('expr','meta')
    tmp <- tmp[!is.na(tmp$meta),]
    if (all(tmp$expr==0) | all(tmp$expr==tmp$expr[1])) {
      tres_pval <- NaN
      tres_z <- NaN
    } else {
      g1 <- tmp$expr[tmp$meta==0]
      g2 <- tmp$expr[tmp$meta==1]
      tres <- t.test(g1,g2)
      tres_pval <- tres$p.value
      tres_z <- tres[["statistic"]][["t"]]
    }
    
    names(tres_pval) <- colnames(pb)[j]
    names(tres_z) <- colnames(pb)[j]
    
    return(list(tres_pval,tres_z))
  },mc.preschedule=TRUE,n.cores=20,progress=TRUE)
  
  ct_res_pv <- lapply(1:length(ct_res),function(x){
    return(ct_res[[x]][[1]])
  })
  ct_res_z <- lapply(1:length(ct_res),function(x){
    return(ct_res[[x]][[2]])
  })
  
  de_genes_neph_pediatric <- c(de_genes_neph_pediatric,ct_res_pv)
  de_genes_z_neph_pediatric <- c(de_genes_z_neph_pediatric,ct_res_z)
  
}
de_genes_neph_pediatric <- unlist(de_genes_neph_pediatric)
de_genes_z_neph_pediatric <- unlist(de_genes_z_neph_pediatric)

# correct pvals
pv_thresh <- .1
de_genes_neph_adj_pediatric <- p.adjust(de_genes_neph_pediatric,method = 'fdr')
de_genes_neph_adj_pediatric <- de_genes_neph_adj_pediatric[!is.na(de_genes_neph_adj_pediatric)]
de_genes_neph_id_pediatric <- names(de_genes_neph_adj_pediatric)[de_genes_neph_adj_pediatric<pv_thresh]
print(length(de_genes_neph_id_pediatric))





de_genes_pred_pediatric <- c()
de_genes_z_pred_pediatric <- c()
meta_var <- 'OS'
for (ct in pbmc_container_full_pediatric$experiment_params$ctypes_use) {
  print(ct)
  # getting out pseudobulk for a cell type
  pb <- pbmc_container_full_pediatric$scMinimal_ctype[[ct]]$pseudobulk
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  
  # run DE test
  # run DE t-test
  ct_res <- plapply(1:ncol(pb),function(j) {
    tmp <- cbind.data.frame(pb[,j],pbmc_container_full_pediatric$donor_metadata[rownames(pb),meta_var])
    colnames(tmp) <- c('expr','meta')
    tmp <- tmp[!is.na(tmp$meta),]
    if (all(tmp$expr==0) | all(tmp$expr==tmp$expr[1])) {
      tres_pval <- NaN
      tres_z <- NaN
    } else {
      g1 <- tmp$expr[tmp$meta==0]
      g2 <- tmp$expr[tmp$meta==1]
      tres <- t.test(g1,g2)
      tres_pval <- tres$p.value
      tres_z <- tres[["statistic"]][["t"]]
    }
    
    names(tres_pval) <- colnames(pb)[j]
    names(tres_z) <- colnames(pb)[j]
    
    return(list(tres_pval,tres_z))
  },mc.preschedule=TRUE,n.cores=20,progress=TRUE)
  
  ct_res_pv <- lapply(1:length(ct_res),function(x){
    return(ct_res[[x]][[1]])
  })
  ct_res_z <- lapply(1:length(ct_res),function(x){
    return(ct_res[[x]][[2]])
  })
  
  de_genes_pred_pediatric <- c(de_genes_pred_pediatric,ct_res_pv)
  de_genes_z_pred_pediatric <- c(de_genes_z_pred_pediatric,ct_res_z)
  
}
de_genes_pred_pediatric <- unlist(de_genes_pred_pediatric)
de_genes_z_pred_pediatric <- unlist(de_genes_z_pred_pediatric)

# correct pvals
pv_thresh <- .1
de_genes_pred_adj_pediatric <- p.adjust(de_genes_pred_pediatric,method = 'fdr')
de_genes_pred_adj_pediatric <- de_genes_pred_adj_pediatric[!is.na(de_genes_pred_adj_pediatric)]
de_genes_pred_id_pediatric <- names(de_genes_pred_adj_pediatric)[de_genes_pred_adj_pediatric<pv_thresh]
print(length(de_genes_pred_id_pediatric))





### plot effect sizes across datasets
print(de_genes_adj[1:10])
print(de_genes_adj_pediatric[1:10])
print(de_genes_sledai_z[1:10])
print(de_genes_z_sledai_pediatric[1:10])

genes_both <- intersect(names(de_genes_sledai_z),names(de_genes_z_sledai_pediatric))
tmp <- cbind.data.frame(de_genes_z_sledai_pediatric[genes_both],de_genes_sledai_z[genes_both])
colnames(tmp) <- c('pSLE','aSLE')
lmres1 <- summary(lm(pSLE~aSLE,data=tmp))
p1 <- ggplot(tmp,aes(x=aSLE,y=pSLE)) +
  geom_point_rast(alpha=.1,size=.5,pch=21, colour='black') +
  ggtitle('SLEDAI DE') +
  xlab('aSLE Z-score') +
  ylab('pSLE Z-score') +
  geom_abline(slope=lmres1$coefficients['aSLE','Estimate'],intercept=lmres1$coefficients[1,'Estimate'],color='black') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
p1

## Figure S3k
pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/DE_cross_dataset_SLEDAI.pdf", useDingbats = FALSE,
    width = 4, height = 2.75)
p1
dev.off()


genes_both <- intersect(names(de_genes_neph_z),names(de_genes_z_neph_pediatric))
tmp <- cbind.data.frame(de_genes_z_neph_pediatric[genes_both],de_genes_neph_z[genes_both])
colnames(tmp) <- c('pSLE','aSLE')
lmres2 <- summary(lm(pSLE~aSLE,data=tmp))
p2 <- ggplot(tmp,aes(x=aSLE,y=pSLE)) +
  geom_point_rast(alpha=.1,size=.5,pch=21, colour='black') +
  ggtitle('Lupus nephritis DE') +
  xlab('aSLE Z-score') +
  ylab('pSLE Z-score') +
  geom_abline(slope=lmres2$coefficients['aSLE','Estimate'],intercept=lmres2$coefficients[1,'Estimate'],color='black') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
p2

## Figure S3l
pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/DE_cross_dataset_neph.pdf", useDingbats = FALSE,
    width = 4, height = 2.75)
p2
dev.off()


genes_both <- intersect(names(de_genes_pred_z),names(de_genes_z_pred_pediatric))
tmp <- cbind.data.frame(de_genes_z_pred_pediatric[genes_both],de_genes_pred_z[genes_both])
colnames(tmp) <- c('pSLE','aSLE')
lmres3 <- summary(lm(pSLE~aSLE,data=tmp))
p3 <- ggplot(tmp,aes(x=aSLE,y=pSLE)) +
  geom_point_rast(alpha=.1,size=.5,pch=21, colour='black') +
  ggtitle('Prednisone DE') +
  xlab('aSLE Z-score') +
  ylab('pSLE Z-score') +
  geom_abline(slope=lmres3$coefficients['aSLE','Estimate'],intercept=lmres3$coefficients[1,'Estimate'],color='black') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
p3

## Figure S3m
pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/DE_cross_dataset_pred.pdf", useDingbats = FALSE,
    width = 4, height = 2.75)
p3
dev.off()

# get gene-factor associations for each dataset

get_z <- function(container, f_test) {
  dsc <- container$tucker_results[[1]][,f_test,drop=FALSE]
  f_res_z <- list() # stores the result from a single factor
  f_res_p <- list() # stores the result from a single factor
  
  # loop through cell types
  for (ct in container$experiment_params$ctypes_use) {
    print(ct)
    pb <- container$scMinimal_ctype[[ct]]$pseudobulk
    # loop through genes
    ct_res <- plapply(1:ncol(pb), function(g_ndx) {
      tmp <- cbind.data.frame(dsc,pb[rownames(dsc),g_ndx])
      colnames(tmp) <- c('dscore','expr')
      lmres <- summary(lm(expr~dscore,data=tmp))
      z_val <- lmres$coefficients['dscore','t value']
      p_val <- lmres$coefficients['dscore','Pr(>|t|)']
      names(z_val) <- paste0(colnames(pb)[g_ndx],'_',ct)
      names(p_val) <- paste0(colnames(pb)[g_ndx],'_',ct)
      return(list(z_val,p_val))
    },mc.preschedule=TRUE,n.cores=20,progress=TRUE)
    ct_res_z <- lapply(ct_res,function(x) {
      return(x[[1]])
    })
    ct_res_p <- lapply(ct_res,function(x) {
      return(x[[2]])
    })
    f_res_z[[ct]] <- unlist(ct_res_z)
    f_res_p[[ct]] <- unlist(ct_res_p)
  }
  return(list(f_res_z,f_res_p))
}


plot_z <- function(aSLE_res,pSLE_res,f_tested) {
  aSLE_z <- unlist(aSLE_res[[1]])
  aSLE_p <- unlist(aSLE_res[[2]])
  pSLE_z <- unlist(pSLE_res[[1]])
  pSLE_p <- unlist(pSLE_res[[2]])
  
  aSLE_z <- aSLE_z[!is.na(aSLE_z)]
  aSLE_p <- aSLE_p[!is.na(aSLE_p)]
  pSLE_z <- pSLE_z[!is.na(pSLE_z)]
  pSLE_p <- pSLE_p[!is.na(pSLE_p)]
  
  aSLE_p <- p.adjust(aSLE_p,method='fdr')
  pSLE_p <- p.adjust(pSLE_p,method='fdr')
  g_sig_pSLE <- names(pSLE_p)[pSLE_p < .001]
  g_sig_aSLE <- names(aSLE_p)[aSLE_p < .001]
  
  genes_both <- intersect(names(aSLE_z),names(pSLE_z))
  tmp <- cbind.data.frame(pSLE_z[genes_both],aSLE_z[genes_both])
  colnames(tmp) <- c('pSLE','aSLE')
  lmres <- summary(lm(pSLE~aSLE,data=tmp))
  p2 <- ggplot(tmp,aes(x=aSLE,y=pSLE)) +
    geom_point_rast(alpha=.1,size=.5,pch=21, colour='black') +
    ggtitle(paste0('scITD F',f_tested,'-gene associations')) +
    xlab(paste0('aSLE_F',f_tested ,' Z-scores')) +
    ylab(paste0('pSLE_F',f_tested, ' Z-scores')) +
    geom_abline(slope=lmres$coefficients['aSLE','Estimate'],intercept=lmres$coefficients[1,'Estimate'],color='black') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  return(list(p2,lmres))
}

aSLE_res <- get_z(pbmc_container_full,1)
pSLE_res <- get_z(pbmc_container_full_pediatric,1)
f1_plot <- plot_z(aSLE_res,pSLE_res,1)

lmres4 <- f1_plot[[2]]
f1_plot <- f1_plot[[1]]

aSLE_res2 <- get_z(pbmc_container_full,2)
pSLE_res2 <- get_z(pbmc_container_full_pediatric,2)
f2_plot <- plot_z(aSLE_res2,pSLE_res2,2)

lmres5 <- f2_plot[[2]]
f2_plot <- f2_plot[[1]]

aSLE_res3 <- get_z(pbmc_container_full,3)
pSLE_res3 <- get_z(pbmc_container_full_pediatric,3)
f3_plot <- plot_z(aSLE_res3,pSLE_res3,3)

lmres6 <- f3_plot[[2]]
f3_plot <- f3_plot[[1]]

## Figure S3k
pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/scITD_cross_dataset_f1.pdf", useDingbats = FALSE,
    width = 4, height = 2.75)
f1_plot
dev.off()

## Figure S3l
pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/scITD_cross_dataset_f2.pdf", useDingbats = FALSE,
    width = 4, height = 2.75)
f2_plot
dev.off()

## Figure S3m
pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/scITD_cross_dataset_f3.pdf", useDingbats = FALSE,
    width = 4, height = 2.75)
f3_plot
dev.off()

plotlist <- list(p1,p2,p3,f1_plot,f2_plot,f3_plot)
lmlist <- list(lmres1,lmres2,lmres3,lmres4,lmres5,lmres6)
# saveRDS(list(plotlist,lmlist),file='/home/jmitchel/data/lupus_data/factor_de_comparison_plots.rds')

lmres1 # rsq 0.1238
lmres2 # rsq 0.01399
lmres3 # rsq 0.03007
lmres4 # rsq 0.5228
lmres5 # rsq 0.2973
lmres6 # rsq 0.2306


# f3_plot +
#   geom_point_rast(alpha=.8,size=.5,pch=21, colour='black')




##### plotting heatmap with DE genes only
pbmc_container <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_container_w_decomp.rds')
pbmc_container <- get_lm_pvals(pbmc_container)


# appending the sledai variable
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_ordinal.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid

dsc <- pbmc_container$tucker_results[[1]]

#
trim_names <- sapply(rownames(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})

orig_names <- names(trim_names)
names(orig_names) <- trim_names
clin_vars <- clin_vars[rownames(clin_vars) %in% names(orig_names),]
rownames(clin_vars) <- orig_names[rownames(clin_vars)]

for (ct in pbmc_container$experiment_params$ctypes_use) {
  ct_meta <- pbmc_container$scMinimal_ctype[[ct]]$metadata
  ct_meta$sledai <- sapply(as.character(ct_meta$donors),function(x){
    return(clin_vars[x,'sledaiscore'])
  })
  pbmc_container$scMinimal_ctype[[ct]]$metadata <- ct_meta
}

full_meta <- pbmc_container$scMinimal_full$metadata
full_meta$sledai <- sapply(as.character(full_meta$donors),function(x){
  return(clin_vars[x,'sledaiscore'])
})
pbmc_container$scMinimal_full$metadata <- full_meta

pbmc_container <- get_donor_meta(pbmc_container,additional_meta='sledai')
head(pbmc_container$donor_metadata)

pbmc_container <- get_de_res(pbmc_container,'sledai')
  
pbmc_container <- plot_loadings_vs_de(pbmc_container, 1, 'SLEDAI')
hm_all <- pbmc_container$hmap_comparison

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/f1_vs_de_hmap.pdf", useDingbats = FALSE,
#     width = 8, height = 5)
hm_all
dev.off()




