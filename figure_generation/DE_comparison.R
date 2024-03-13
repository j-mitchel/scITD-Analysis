
library(Seurat)
library(ggplot2)
library(readxl)
library(ComplexHeatmap)
library(MASS)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(sccore)
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

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Ethnicity'),
                                        stat_use='pval')




## plot donor scores to make sure it looks the same as before
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

pbmc_container$plots$donor_matrix


# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)


# ## load up the clinical variables to use for DE
# clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_categorical.xlsx')
# clin_vars <- as.data.frame(clin_vars)
# rownames(clin_vars) <- clin_vars$subjectid
# clin_vars$subjectid <- NULL
# clin_vars$abs_either <- sapply(1:nrow(clin_vars),function(x){
#   if (clin_vars$acrantidsdna[x]==1 | clin_vars$acrantismith[x]==1) {
#     return(1)
#   } else {
#     return(0)
#   }
# })

## using ordinal variables
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_ordinal.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL
# clin_vars$sledaiscore_bin <- sapply(clin_vars$sledaiscore,function(x){
#   if (x>4) {
#     return(1)
#   } else {
#     return(0)
#   }
# })

# ## to use a combined version of antiantibodies and sledscore
# clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_categorical.xlsx')
# clin_vars <- as.data.frame(clin_vars)
# rownames(clin_vars) <- clin_vars$subjectid
# clin_vars$subjectid <- NULL
# clin_vars_ord <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_ordinal.xlsx')
# clin_vars_ord <- as.data.frame(clin_vars_ord)
# rownames(clin_vars_ord) <- clin_vars_ord$subjectid
# clin_vars_ord$subjectid <- NULL
# clin_vars_ord$sledaiscore_bin <- sapply(clin_vars_ord$sledaiscore,function(x){
#   if (x>6) {
#     return(1)
#   } else {
#     return(0)
#   }
# })
# clin_vars$sledaiscore_bin <- clin_vars_ord[rownames(clin_vars),'sledaiscore_bin']
# clin_vars$abs_either <- sapply(1:nrow(clin_vars),function(x){
#   if (clin_vars$acrantidsdna[x]==1 | clin_vars$acrantismith[x]==1 | clin_vars$sledaiscore_bin[x]==1) {
#     return(1)
#   } else {
#     return(0)
#   }
# })


# ## using meds variables
# clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_meds_cleaned.xlsx')
# clin_vars <- as.data.frame(clin_vars)
# rownames(clin_vars) <- clin_vars[,'Sample ID']
# clin_vars[,'Sample ID'] <- NULL
# clin_vars[is.na(clin_vars)] <- 0
# ##

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
# meta_var <- 'acrantidsdna'
# meta_var <- 'abs_either'
# meta_var <- 'prednisone'
# meta_var <- 'sledaiscore_bin'
meta_var <- 'sledaiscore'
for (ct in pbmc_container$experiment_params$ctypes_use) {
  print(ct)
  # getting out pseudobulk for a cell type
  pb <- pbmc_container$scMinimal_ctype[[ct]]$pseudobulk
  rownames(pb) <- trim_names[rownames(pb)]
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  
  # run DE t-test
  for (j in 1:ncol(pb)) {
    # append metadata variable
    tmp <- cbind.data.frame(pb[,j],clin_vars[rownames(pb),meta_var])
    colnames(tmp) <- c('expr','meta')
    # # run test
    # g1 <- tmp$expr[tmp$meta==0]
    # g2 <- tmp$expr[tmp$meta==1]
    # tres <- wilcox.test(g1,g2)
    # tres <- t.test(g1,g2)
    # tres_pval <- tres$p.value
    
    # trying linear test instead
    lmres <- summary(lm(meta~expr,data=tmp))
    if (all(tmp$expr==0)) {
      tres_pval <- NaN
    } else {
      tres_pval <- lmres$coefficients['expr','Pr(>|t|)']
    }
    
    # ## using the ordinal regression
    # tmp$meta <- as.factor(tmp$meta)
    # if (all(tmp$expr==0)) {
    #   tres_pval <- NaN
    # } else {
    #   tres_pval <- tryCatch({
    #     m <- polr(meta ~ ., data = tmp, Hess=TRUE, method='probit')
    #     ## view a summary of the model
    #     ctable <- coef(summary(m))
    #     ## calculate and store p values
    #     p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
    #     tres_pval <- p[1]
    #   },error = function(cond) {
    #     return(NA)
    #   })
    # }
    
    names(tres_pval) <- colnames(pb)[j]
    de_genes <- c(de_genes,tres_pval)
  }
}
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

g_ct_cev_de <- sapply(colnames(pb_all),function(x){
  return(strsplit(x,split='_')[[1]][[2]])
})

annot_col <- brewer.pal(7, 'Set3')
names(annot_col) <- unique(g_ct_cev_de)
ra_de <- ComplexHeatmap::rowAnnotation(df = g_ct_cev_de, show_annotation_name=FALSE, col = list(df = annot_col))
col_fun = colorRamp2(c(0, 1), c("white", "red"))
# clust_method <- 'single'
clust_method <- 'complete'
hmap_de <- Heatmap(de_cormat,name = "pearson r",
        cluster_columns = TRUE,
        col = col_fun,
        clustering_method_rows = clust_method,
        clustering_method_columns = clust_method,
        cluster_rows = TRUE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        right_annotation=ra_de,
        border = TRUE,
        column_title = 'anti-dsDNA DE genes',
        column_title_gp = gpar(fontsize = 20))

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/SLEDAI_cor_hmap.pdf", useDingbats = FALSE,
#     width = 4, height = 3.5)
hmap_de
dev.off()




### now generate the heatmap for the f1 significant genes
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
  # pb <- pb[names(trim_names),rownames(f1_pvals)[f1_pvals[,ct]<pv_thresh]]
  pb <- pb[names(trim_names),rownames(f1_pvals2)[f1_pvals2[,ct]<pv_thresh]]
  colnames(pb) <- paste0(colnames(pb),'_',ct)
  # append pb to existing pb
  pb_all <- cbind(pb_all,pb)
}

## to downsample genes randomly
# colsamp <- sample(1:ncol(pb_all),ncol(de_cormat))
# colsamp <- sample(1:ncol(pb_all),ncol(pb_all)/2)
# pb_all <- pb_all[,colsamp]

scITD_cormat <- abs(cor(pb_all,method='pearson'))

g_ct_cev_f1 <- sapply(colnames(pb_all),function(x){
  return(strsplit(x,split='_')[[1]][[2]])
})

ra_f1 <- ComplexHeatmap::rowAnnotation(df = g_ct_cev_f1, show_annotation_name=FALSE, col = list(df = annot_col))

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
        column_title = 'F1 significant genes',
        column_title_gp = gpar(fontsize = 20))




# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/F1_cor_hmap.pdf", useDingbats = FALSE,
#     width = 4, height = 3.5)
scITD_hmap
dev.off()

print(dim(scITD_cormat))



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

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/pairwise_cor_de_f1.pdf", useDingbats = FALSE,
#     width = 6, height = 3.5)
p
dev.off()


### computing jaccard overlap for either de set or our comparable top set with the in vitro ifn stimulation data
load("/home/jmitchel/data/lupus_data/cell.type.diffexp.RData")

B.expressed.res.sig <- cd19.expressed.res$featureData.symbol[cd19.expressed.res$padj<.05]
nk.expressed.res.sig <- nk.expressed.res$featureData.symbol[nk.expressed.res$padj<.05]
cM.expressed.res.sig <- cd14.expressed.res$featureData.symbol[cd14.expressed.res$padj<.05]
cd4.expressed.res.sig <- cd4.expressed.res$featureData.symbol[cd4.expressed.res$padj<.05]
cd8.expressed.res.sig <- cd8.expressed.res$featureData.symbol[cd8.expressed.res$padj<.05]
ncM.expressed.res.sig <- ncmono.expressed.res$featureData.symbol[ncmono.expressed.res$padj<.05]

ct_test <- list(B.expressed.res.sig,nk.expressed.res.sig,cM.expressed.res.sig,cd4.expressed.res.sig,cd8.expressed.res.sig,ncM.expressed.res.sig)
names(ct_test) <- c('B','NK','cMono','Th','Tc','ncMono')
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
total_intersect_de <- 0
total_intersect_scITD <- 0
for (ct in names(ct_test)) {
  iv_res <- ct_test[[ct]]
  de_res <- de_genes_only[de_genes_ct==ct]
  ilen_de <- length(intersect(de_res,iv_res))
  total_intersect_de <- total_intersect_de + ilen_de
  scITD_res <- scITD_genes_only[scITD_genes_ct==ct]
  ilen_scITD <- length(intersect(scITD_res,iv_res))
  total_intersect_scITD <- total_intersect_scITD + ilen_scITD
}
total_intersect_de
total_intersect_scITD
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
#     width = 4, height = 3.5)
p
dev.off()



## plotting enriched gene sets in each set of genes against background
# should show two example gene sets with the genes colored and highlighted on the correlation plot to indicate it's a separate process
# and ideally it should be a process that typically isn't thought to be caused by IFN or to cause IFN itself.


de_genes_unq <- unique(de_genes_only)
scITD_genes_unq <- unique(scITD_genes_only)
bg <- names(de_genes)[!is.na(de_genes)]
bg <- sapply(bg,function(x){
  strsplit(x,split='_')[[1]][[1]]
})

go_res_de <- enrichGO(gene = de_genes_unq,
                   universe = bg,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = 'CC',
                 pAdjustMethod = "BH",
                 pvalueCutoff  = .01,
                 qvalueCutoff  = .05,
                 minGSSize=10,
                 maxGSSize=1500)


go_res_de <- go_res_de@result

go_res_f1 <- enrichGO(gene = scITD_genes_unq,
                      universe = bg,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'SYMBOL',
                      ont           = 'CC',
                      pAdjustMethod = "BH",
                      pvalueCutoff  = .01,
                      qvalueCutoff  = .05,
                      minGSSize=10,
                      maxGSSize=1500)


go_res_f1 <- go_res_f1@result

### looking at genes sets sig in one but not the other
go_res_de_sub <- go_res_de$Description[go_res_de$p.adjust<.05]
go_res_f1_sub <- go_res_f1$Description[go_res_f1$p.adjust<.05]
go_res_de_sub[!(go_res_de_sub%in%go_res_f1_sub)]
go_res_f1_sub[!(go_res_f1_sub%in%go_res_de_sub)]

go_res_de_ns <- go_res_de[go_res_de$pvalue>.1,]
go_res_de_ns <- go_res_de_ns$Description[go_res_de_ns$Count<=2]
go_res_f1_sub[go_res_f1_sub%in%go_res_de_ns]



genes_lab_de1 <- go_res_de[go_res_de$Description=='cellular response to type I interferon','geneID']
genes_lab_de1 <- go_res_de[go_res_de$Description=="exonuclease activity, active with either ribo- or deoxyribonucleic acids and producing 5'-phosphomonoesters" ,'geneID']
# genes_lab_de2 <- go_res_de[go_res_de$Description=='response to cadmium ion','geneID']
# genes_lab_de2 <- go_res_de[go_res_de$Description=="response to zinc ion",'geneID']
genes_lab_de2 <- go_res_de[go_res_de$Description=="catecholamine secretion",'geneID']
genes_lab_de2 <- go_res_de[go_res_de$Description=="regulation of chromosome segregation",'geneID']
# genes_lab_de2 <- go_res_de[go_res_de$Description=='chromosome segregation','geneID']
# genes_lab_de2 <- go_res_de[go_res_de$Description=='regulation of chromosome segregation','geneID']
# genes_lab_de2 <- go_res_de[go_res_de$Description=='hematopoietic progenitor cell differentiation','geneID']
genes_lab_scITD1 <- go_res_f1[go_res_f1$Description=='cellular response to type I interferon','geneID']


genes_lab_de1 <- unlist(strsplit(genes_lab_de1,split='/'))
genes_lab_de2 <- unlist(strsplit(genes_lab_de2,split='/'))
genes_lab_de <- cbind.data.frame(c(de_genes_only %in% genes_lab_de1),c(de_genes_only %in% genes_lab_de2))
rownames(genes_lab_de) <- de_genes_id
colnames(genes_lab_de) <- c('GO_IFN','GO_cadmium')
genes_lab_de$GO_IFN <- factor(genes_lab_de$GO_IFN,levels=c(TRUE,FALSE))
genes_lab_de$GO_cadmium <- factor(genes_lab_de$GO_cadmium,levels=c(TRUE,FALSE))
levels(genes_lab_de$GO_IFN) <- c('in_set','not_in_set')
levels(genes_lab_de$GO_cadmium) <- c('in_set','not_in_set')

hres <- hclust(as.dist(1-de_cormat),method = 'ward.D2')

## manually reordering rows and cols by clustering
de_cormat2 <- de_cormat[hres[["order"]],hres[["order"]]]
hclusts <- hclusts[hres[["order"]]]

genes_lab_de <- genes_lab_de[hres[["order"]],]


de_col <- c('darkblue','gray95')
names(de_col) <- c('in_set','not_in_set')
ca_de <- ComplexHeatmap::HeatmapAnnotation(df_ifn = genes_lab_de$GO_IFN, df_cadmium = genes_lab_de$GO_cadmium,
                                           col = list(df_ifn = de_col, df_cadmium = de_col),
                                           show_legend = FALSE,show_annotation_name = c(TRUE, TRUE),
                                           annotation_name_side = c("left", "left"))

col_fun = colorRamp2(c(0, 1), c("white", "red"))
hmap_de <- Heatmap(de_cormat2,name = "pearson r",
                   cluster_columns = FALSE,
                   col = col_fun,
                   cluster_rows = FALSE,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   show_row_dend = FALSE,
                   show_column_dend = FALSE,
                   # right_annotation=ra_de,
                   bottom_annotation=ca_de,
                   border = TRUE,
                   column_title = 'anti-dsDNA DE genes',
                   column_title_gp = gpar(fontsize = 20))

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/SLEDAI_cor_hmap_gset.pdf", useDingbats = FALSE,
#     width = 4, height = 3.5)
hmap_de
dev.off()

### looking to see what genes are specifically in the separate clusters to try go with
# hres <- hclust(as.dist(1/de_cormat),method = 'complete')
# hres <- hclust(as.dist(1/de_cormat),method = 'ward.D2')
hres <- hclust(as.dist(1-de_cormat),method = 'ward.D2')
# hres <- hclust(as.dist(1-de_cormat),method = 'ward.D2')
hclusts <- cutree(hres, k = 5)
hclusts <- cutree(hres, k = 3)
table(hclusts)

group_col <- c('darkblue','darkgreen','darkorange','darkred','lightblue','yellow')
names(group_col) <- c(1,2,3,4,5,6)

group_col <- c('darkblue','darkgreen','darkorange','darkred','lightblue')
names(group_col) <- c(1,2,3,4,5)

group_col <- c('darkblue','darkgreen','darkorange','darkred')
names(group_col) <- c(1,2,3,4)

group_col <- c('darkblue','darkgreen','darkorange')
names(group_col) <- c(1,2,3)

group_col <- c('darkblue','darkorange')
names(group_col) <- c(1,2)

## manually reordering rows and cols by clustering
de_cormat2 <- de_cormat[hres[["order"]],hres[["order"]]]
hclusts <- hclusts[hres[["order"]]]

ca_grp <- ComplexHeatmap::HeatmapAnnotation(df = hclusts,
                                           col = list(df = group_col),
                                           show_legend = TRUE)

col_fun = colorRamp2(c(0, 1), c("white", "red"))
# clust_method <- 'single'
# clust_method <- 'complete'
hmap_de <- Heatmap(de_cormat2,name = "pearson r",
                   cluster_columns = FALSE,
                   col = col_fun,
                   cluster_rows = FALSE,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   show_row_dend = FALSE,
                   show_column_dend = FALSE,
                   bottom_annotation=ca_grp,
                   border = TRUE,
                   column_title = 'anti-dsDNA DE genes',
                   column_title_gp = gpar(fontsize = 20))
hmap_de

c2_g <- names(hclusts)[hclusts==2]
c3_g <- names(hclusts)[hclusts==3]
c2_g <- sapply(c2_g,function(x){
  strsplit(x,split='_')[[1]][[1]]
})
c3_g <- sapply(c3_g,function(x){
  strsplit(x,split='_')[[1]][[1]]
})

go_res_c2 <- enrichGO(gene = c2_g,
                      universe = bg,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'SYMBOL',
                      ont           = 'BP',
                      pAdjustMethod = "BH",
                      pvalueCutoff  = .01,
                      qvalueCutoff  = .05,
                      minGSSize=10,
                      maxGSSize=500)
go_res_c2 <- go_res_c2@result

go_res_c3 <- enrichGO(gene = c3_g,
                      universe = bg,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'SYMBOL',
                      ont           = 'BP',
                      pAdjustMethod = "BH",
                      pvalueCutoff  = .01,
                      qvalueCutoff  = .05,
                      minGSSize=10,
                      maxGSSize=300)
go_res_c3 <- go_res_c3@result


genes_lab_scITD1 <- unlist(strsplit(genes_lab_scITD1,split='/'))


















#### recomputing associated genes using all genes

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


