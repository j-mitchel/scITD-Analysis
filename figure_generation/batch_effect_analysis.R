
library(Seurat)
library(ComplexHeatmap)
library(readxl)
library(MASS)
library(limma)
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

pbmc_container <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_container_w_decomp.rds')





#### now running it without batch correction
pbmc_container_uncorrected <- make_new_container(seurat_obj=pbmc,
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


pbmc_container_uncorrected <- form_tensor(pbmc_container_uncorrected, donor_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5)

# pbmc_container_uncorrected <- run_tucker_ica(pbmc_container_uncorrected, ranks=c(7,20),
#                                  tucker_type = 'regular', rotation_type = 'hybrid')

# pbmc_container_uncorrected <- run_tucker_ica(pbmc_container_uncorrected, ranks=c(12,36),
#                                              tucker_type = 'regular', rotation_type = 'hybrid')

pbmc_container_uncorrected <- run_tucker_ica(pbmc_container_uncorrected, ranks=c(13,40),
                                             tucker_type = 'regular', rotation_type = 'hybrid')


### plot donor scores
pbmc_container_uncorrected <- get_meta_associations(pbmc_container_uncorrected,vars_test=c('sex','Age','pool','processing','Ethnicity'),
                                        stat_use='pval')

## plot donor scores to make sure it looks the same as before
pbmc_container_uncorrected <- plot_donor_matrix(pbmc_container_uncorrected,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval',
                                    cluster_by_meta='pool',
                                    meta_vars='pool')

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/batch_dscores.pdf", useDingbats = FALSE,
#     width = 7, height = 6)
pbmc_container_uncorrected$plots$donor_matrix
# dev.off()



### adding the same metadata variables to the top of the plot
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_categorical.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL

# get tucker donor scores to test
dsc_uncorrected <- pbmc_container_uncorrected$tucker_results[[1]]

## get donors in both dsc and in clin_vars
# trim donor IDs in dsc
trim_names <- sapply(rownames(dsc_uncorrected), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
old_names <- rownames(dsc_uncorrected)
names(old_names) <- trim_names
rownames(dsc_uncorrected) <- trim_names

# get donors in both dataframes
d_both <- rownames(clin_vars)[rownames(clin_vars) %in% rownames(dsc_uncorrected)]

# limit both dataframes to just the intersection of donors and in the same order
dsc_uncorrected <- dsc_uncorrected[d_both,]
clin_vars <- clin_vars[d_both,]

## calculate associations
all_pvals <- c()
f_tested <- c()
c_tested <- c()
# loop through the variables to test
for (j in 1:ncol(clin_vars)) {
  print(j)
  # loop through factors
  for (f in 1:ncol(dsc_uncorrected)) {
    d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]
    
    tmp <- as.data.frame(cbind(dsc_uncorrected[d_keep,f], clin_vars[d_keep,j]))
    colnames(tmp) <- c('dscore','cvar')
    
    # force cvar to be factor
    tmp$cvar <- as.factor(tmp$cvar)
    
    # if smallest level has less thatn n donors skip this one
    if (min(table(tmp$cvar)) < 20) {
      next
    }
    
    # logistic regression model
    fmod <- glm(cvar~dscore, data=tmp, family = "binomial") ##"full" mod
    nmod <- glm(cvar~1, data=tmp, family = 'binomial') ##"null" mod
    a_res <- anova(nmod, fmod, test = 'Chisq')
    pval <- a_res$`Pr(>Chi)`[2]
    
    all_pvals <- c(all_pvals,pval)
    f_tested <- c(f_tested,f)
    c_tested <- c(c_tested,colnames(clin_vars)[j])
    
  }
}
all_pvals <- p.adjust(all_pvals,method='fdr')
all_pvals[order(all_pvals,decreasing=FALSE)][1:10]
f_tested[order(all_pvals,decreasing=FALSE)][1:10]
c_tested[order(all_pvals,decreasing=FALSE)][1:10]

res_add <- matrix(NA,nrow=4,ncol=ncol(dsc_uncorrected))

# get a row of padj for autoantibody associations and others
ndx_extract <- which(c_tested=='acrantidsdna')
f_ordering <- f_tested[ndx_extract]
c_pv <- all_pvals[ndx_extract]
ndx_reorder <- order(f_ordering)
print(f_ordering[ndx_reorder])
print(c_pv[ndx_reorder])
res_add[1,] <- c_pv[ndx_reorder]

ndx_extract <- which(c_tested=='acrantismith')
f_ordering <- f_tested[ndx_extract]
c_pv <- all_pvals[ndx_extract]
ndx_reorder <- order(f_ordering)
print(f_ordering[ndx_reorder])
print(c_pv[ndx_reorder])
res_add[2,] <- c_pv[ndx_reorder]




##### testing factors against the ordinal variables
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_ordinal.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL

clin_vars$sliccmalignancy <- NULL
clin_vars$lupusseverityindex <- NULL
clin_vars$smokestat <- NULL
clin_vars$acrcsum <- NULL
clin_vars$sliccavasnec <- NULL
clin_vars$slicccva <- NULL

clin_vars <- clin_vars[d_both,]

all_pvals <- c()
f_tested <- c()
c_tested <- c()
for (j in 1:ncol(clin_vars)) {
  print(j)
  
  # get donors in clin var that don't have an NA value
  d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]
  
  tmp <- as.data.frame(cbind(dsc_uncorrected[d_keep,], clin_vars[d_keep,j]))
  colnames(tmp)[1:ncol(dsc_uncorrected)] <- sapply(1:ncol(dsc_uncorrected),function(x){paste0('Factor',x)})
  colnames(tmp)[ncol(dsc_uncorrected)+1] <- 'cvar'
  
  # using ordinal logistic regression
  tmp$cvar <- as.factor(tmp$cvar)
  m <- polr(cvar ~ ., data = tmp, Hess=TRUE, method='probit')
  
  ## view a summary of the model
  ctable <- coef(summary(m))
  ## calculate and store p values
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  pval <- p[1:ncol(dsc_uncorrected)]
  
  all_pvals <- c(all_pvals,pval)
  f_tested <- c(f_tested,1:ncol(dsc_uncorrected))
  c_tested <- c(c_tested,rep(colnames(clin_vars)[j],ncol(dsc_uncorrected)))
}
all_pvals <- p.adjust(all_pvals,method='fdr')
all_pvals[order(all_pvals,decreasing=FALSE)][1:10]
f_tested[order(all_pvals,decreasing=FALSE)][1:10]
c_tested[order(all_pvals,decreasing=FALSE)][1:10]


# get a row of padj for SLEDAI associations
ndx_extract <- which(c_tested=='sledaiscore')
f_ordering <- f_tested[ndx_extract]
c_pv <- all_pvals[ndx_extract]
ndx_reorder <- order(f_ordering)
print(f_ordering[ndx_reorder])
print(c_pv[ndx_reorder])
res_add[3,] <- c_pv[ndx_reorder]


clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_meds_cleaned.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars[,'Sample ID']
clin_vars[,'Sample ID'] <- NULL

# make all NA into zeros, since no 0 are put in the table
clin_vars[is.na(clin_vars)] <- 0

# separate out pred dose as it's the only continuous variable here
pred_dose <- clin_vars[,'pred_dose',drop=FALSE]
clin_vars[,'pred_dose'] <- NULL

# make sure there are no columns of all zeros
colSums(clin_vars)

# need to remove a few columns that have only 1 or 0 donors on the med
clin_vars[,c('solumedrol','rx_abatacept','rx_cyclophosphamide','rx_etanercept',
             'rx_IGG','rx_leflunomide','rx_rituximab','rx_sulfasalazine')] <- NULL

# get tucker donor scores to test
dsc_uncorrected <- pbmc_container_uncorrected$tucker_results[[1]]

## get donors in both dsc and in clin_vars
# trim donor IDs in dsc
trim_names <- sapply(rownames(dsc_uncorrected), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
old_names <- rownames(dsc_uncorrected)
names(old_names) <- trim_names
rownames(dsc_uncorrected) <- trim_names

# get donors in both dataframes
d_both <- rownames(clin_vars)[rownames(clin_vars) %in% rownames(dsc_uncorrected)]

# limit both dataframes to just the intersection of donors and in the same order
dsc_uncorrected <- dsc_uncorrected[d_both,]
clin_vars <- clin_vars[d_both,]

all_pvals <- c()
f_tested <- c()
c_tested <- c()
# loop through the variables to test
for (j in 1:ncol(clin_vars)) {
  print(j)
  # loop through factors
  for (f in 1:ncol(dsc_uncorrected)) {
    # get donors in clin var that don't have an NA value
    d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]
    
    tmp <- as.data.frame(cbind(dsc_uncorrected[d_keep,f], clin_vars[d_keep,j]))
    colnames(tmp) <- c('dscore','cvar')
    
    # force cvar to be factor
    tmp$cvar <- as.factor(tmp$cvar)
    
    # trying with logistic regression model
    fmod <- glm(cvar~dscore, data=tmp, family = "binomial") ##"full" mod
    nmod <- glm(cvar~1, data=tmp, family = 'binomial') ##"null" mod
    a_res <- anova(nmod, fmod, test = 'Chisq')
    pval <- a_res$`Pr(>Chi)`[2]
    
    all_pvals <- c(all_pvals,pval)
    f_tested <- c(f_tested,f)
    c_tested <- c(c_tested,colnames(clin_vars)[j])
  }
}
all_pvals <- p.adjust(all_pvals,method='fdr')
all_pvals[order(all_pvals,decreasing=FALSE)]
f_tested[order(all_pvals,decreasing=FALSE)]
c_tested[order(all_pvals,decreasing=FALSE)]

ndx_extract <- which(c_tested=='prednisone')
f_ordering <- f_tested[ndx_extract]
c_pv <- all_pvals[ndx_extract]
ndx_reorder <- order(f_ordering)
print(f_ordering[ndx_reorder])
print(c_pv[ndx_reorder])
res_add[4,] <- c_pv[ndx_reorder]

rownames(res_add) <- c('anti_dsDNA','anti_Smith','SLEDAI','Prednisone')

pbmc_container_uncorrected$meta_associations <- rbind(pbmc_container_uncorrected$meta_associations,res_add)

# reorder rows
pbmc_container_uncorrected$meta_associations <- pbmc_container_uncorrected$meta_associations[c('sex','Age','pool',
                                                                       'processing','Ethnicity',
                                                                       'anti_dsDNA','anti_Smith',
                                                                       'SLEDAI','Prednisone'),]

## plot donor score
pbmc_container_uncorrected <- plot_donor_matrix(pbmc_container_uncorrected,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval',
                                    meta_vars=c('pool'),
                                    cluster_by_meta='pool')

pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/uncorrected_dscores2.pdf", useDingbats = FALSE,
    width = 7.5, height = 6)
pbmc_container_uncorrected$plots$donor_matrix
dev.off()





### plot correlations between factor scores from corrected to uncorrected
dsc <- pbmc_container$tucker_results[[1]]
dsc_uncorrected <- pbmc_container_uncorrected$tucker_results[[1]]
cormat <- abs(cor(dsc,dsc_uncorrected))
colnames(cormat) <- paste0('uncorrected_F',1:ncol(cormat))
rownames(cormat) <- paste0('corrected_F',1:nrow(cormat))

col_fun = colorRamp2(c(0, 1), c("white", "red"))
hmap <- Heatmap(cormat,name = "pearson r",
                cluster_columns = FALSE,
                col = col_fun,
                cluster_rows = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                column_names_side = 'top',
                row_names_side = 'left',
                column_names_rot = 30,
                column_names_gp = grid::gpar(fontsize = 8),
                row_names_gp = grid::gpar(fontsize = 8),
                border = TRUE,
                column_title = 'Donor scores comparison',
                column_title_gp = gpar(fontsize = 14),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid::grid.text(sprintf("%.2f", cormat[i, j]), x, y, gp = gpar(fontsize = 10))
                })


pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/corrected_vs_uncorrected_dscores2.pdf", useDingbats = FALSE,
    width = 7, height = 4)
hmap
dev.off()








### trying with limma
pbmc_container_cor2 <- make_new_container(seurat_obj=pbmc,
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

pbmc_container_cor2 <- parse_data_by_ctypes(pbmc_container_cor2)
pbmc_container_cor2 <- clean_data(pbmc_container_cor2, donor_min_cells=20)
pbmc_container_cor2 <- get_pseudobulk(pbmc_container_cor2)
pbmc_container_cor2 <- normalize_pseudobulk(pbmc_container_cor2, method='trim', scale_factor=10000)
pbmc_container_cor2 <- get_normalized_variance(pbmc_container_cor2)
pbmc_container_cor2 <- get_ctype_vargenes(pbmc_container_cor2, method='norm_var_pvals', thresh=.15)
pbmc_container_cor2 <- scale_variance(pbmc_container_cor2,var_scale_power=.5)

### run different batch correction
for (ct in pbmc_container_cor2$experiment_params$ctypes_use) {
  pb <- pbmc_container_cor2$scMinimal_ctype[[ct]]$pseudobulk
  metadata <- unique(pbmc_container_cor2$scMinimal_ctype[[ct]]$metadata[,c('donors', 'pool')])
  rownames(metadata) <- metadata$donors
  metadata <- metadata[rownames(pb),]
  pb_cln <- removeBatchEffect(t(pb), batch=metadata$pool,
                    design=matrix(1,ncol(t(pb)),1))
  pbmc_container_cor2$scMinimal_ctype[[ct]]$pseudobulk <- t(pb_cln)
}

pbmc_container_cor2 <- stack_tensor(pbmc_container_cor2)
pbmc_container_cor2 <- run_tucker_ica(pbmc_container_cor2, ranks=c(7,20),
                                      tucker_type = 'regular', rotation_type = 'hybrid')


### plot donor scores
pbmc_container_cor2 <- get_meta_associations(pbmc_container_cor2,vars_test=c('sex','Age','pool','processing','Ethnicity'),
                                             stat_use='pval')

## plot donor scores to make sure it looks the same as before
pbmc_container_cor2 <- plot_donor_matrix(pbmc_container_cor2,
                                         show_donor_ids = FALSE,
                                         add_meta_associations='pval',
                                         cluster_by_meta='pool',
                                         meta_vars='pool')

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/limma_dscores.pdf", useDingbats = FALSE,
#     width = 7, height = 6)
pbmc_container_cor2$plots$donor_matrix
dev.off()



### plot correlations between factor scores from combat to limma
dsc <- pbmc_container$tucker_results[[1]]
dsc_limma <- pbmc_container_cor2$tucker_results[[1]]
cormat <- abs(cor(dsc,dsc_limma))
colnames(cormat) <- paste0('limma_F',1:ncol(cormat))
rownames(cormat) <- paste0('combat_F',1:nrow(cormat))

col_fun = colorRamp2(c(0, 1), c("white", "red"))
hmap <- Heatmap(cormat,name = "pearson r",
                cluster_columns = FALSE,
                col = col_fun,
                cluster_rows = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                column_names_side = 'top',
                row_names_side = 'left',
                column_names_rot = 30,
                column_names_gp = grid::gpar(fontsize = 8),
                row_names_gp = grid::gpar(fontsize = 8),
                border = TRUE,
                column_title = 'Donor scores comparison',
                column_title_gp = gpar(fontsize = 14),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid::grid.text(sprintf("%.2f", cormat[i, j]), x, y, gp = gpar(fontsize = 10))
                })


# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/combat_vs_limma_dscores.pdf", useDingbats = FALSE,
#     width = 4.5, height = 3.5)
hmap
dev.off()






##### doing a comparison with different normalization
pbmc_container_reg_norm <- make_new_container(seurat_obj=pbmc,
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


pbmc_container_reg_norm <- form_tensor(pbmc_container_reg_norm, donor_min_cells=20,
                              norm_method='regular', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')

pbmc_container_reg_norm <- run_tucker_ica(pbmc_container_reg_norm, ranks=c(7,20),
                                 tucker_type = 'regular', rotation_type = 'hybrid')


### plot correlations between factor scores from trim to regular normalization
dsc <- pbmc_container$tucker_results[[1]]
dsc_reg <- pbmc_container_reg_norm$tucker_results[[1]]
cormat <- abs(cor(dsc,dsc_reg))
colnames(cormat) <- paste0('regnorm_F',1:ncol(cormat))
rownames(cormat) <- paste0('trimm_F',1:nrow(cormat))

col_fun = colorRamp2(c(0, 1), c("white", "red"))
hmap <- Heatmap(cormat,name = "pearson r",
                cluster_columns = FALSE,
                col = col_fun,
                cluster_rows = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                column_names_side = 'top',
                row_names_side = 'left',
                column_names_rot = 30,
                column_names_gp = grid::gpar(fontsize = 8),
                row_names_gp = grid::gpar(fontsize = 8),
                border = TRUE,
                column_title = 'Donor scores comparison',
                column_title_gp = gpar(fontsize = 14),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid::grid.text(sprintf("%.2f", cormat[i, j]), x, y, gp = gpar(fontsize = 10))
                })


# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/norm_compare_dscores.pdf", useDingbats = FALSE,
#     width = 4.5, height = 3.5)
hmap
dev.off()



### now comparing loadings
lds <- pbmc_container$tucker_results[[2]]
lds_reg <- pbmc_container_reg_norm$tucker_results[[2]]
g_both <- intersect(colnames(lds_reg),colnames(lds))
lds <- lds[,g_both]
lds_reg <- lds_reg[,g_both]
cormat <- abs(cor(t(lds),t(lds_reg)))
colnames(cormat) <- paste0('regnorm_F',1:ncol(cormat))
rownames(cormat) <- paste0('trimm_F',1:nrow(cormat))

col_fun = colorRamp2(c(0, 1), c("white", "red"))
hmap <- Heatmap(cormat,name = "pearson r",
                cluster_columns = FALSE,
                col = col_fun,
                cluster_rows = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                column_names_side = 'top',
                row_names_side = 'left',
                column_names_rot = 30,
                column_names_gp = grid::gpar(fontsize = 8),
                row_names_gp = grid::gpar(fontsize = 8),
                border = TRUE,
                column_title = 'Loadings comparison',
                column_title_gp = gpar(fontsize = 14),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid::grid.text(sprintf("%.2f", cormat[i, j]), x, y, gp = gpar(fontsize = 10))
                })

hmap


