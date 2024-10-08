library(Seurat)
library(cowplot)
library(scITD)
library(devtools)
library(cacoa)
library(readxl)
library(ggplot2)

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


# putting the two umaps side by side
p1 <- DimPlot(pbmc, reduction = "umap", group.by = 'clusters_fine', label = TRUE, repel = TRUE) +
  NoLegend() +
  ggtitle('Cell subtypes')
p2 <- DimPlot(pbmc, reduction = "umap", group.by = 'clusters_coarse', label = TRUE, repel = TRUE) +
  NoLegend() +
  ggtitle('Cell types analyzed')

fig <- plot_grid(p2,p1)

### Figure S2e
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/pediatric_umap.pdf", useDingbats = FALSE,
#     width = 8, height = 4)
fig
# dev.off()


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

pbmc_container <- make_new_container(seurat_obj=pbmc,
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

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=2,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='Batch')


pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(3,10),n_iterations=50,subset_type='subset', sub_prop=.85)
pbmc_container$plots$stability_plot_dsc

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(3,10),
                                 tucker_type = 'regular', rotation_type = 'hybrid')

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','SLE_status','Age','Batch','Ethnicity'),
                                        stat_use='pval')

# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

pbmc_container$plots$donor_matrix









### load up the lupus dataset to do projections: see preprocessing/lupus_preprocessing.R 
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




# project all original factors onto the pediatric dataset
pbmc_container <- project_new_data(pbmc_container,pbmc_container_SLE)

pbmc_container <- get_donor_meta(pbmc_container, additional_meta = c(meta_cols,'Ethnicity'))


colnames(pbmc_container$donor_metadata)
table(pbmc_container$donor_metadata$Neph_all)
table(pbmc_container$donor_metadata$KIDNEY)
table(pbmc_container$donor_metadata$SERUM)
table(pbmc_container$donor_metadata$Memb_Neph)
table(pbmc_container$donor_metadata$ALBUMIN)
table(pbmc_container$donor_metadata$CREATININE)
table(pbmc_container$donor_metadata$low_com)

# listing out all factor-metadata combos to test associations for
factors_test <- c(1,1,1,1,1,2,2,2,2,3)
meta_vars <- c('SLEDAI','dsDNA','MONOCYTE_per','LYMPHOCYTE_per','LYM_ABS',
               'Neph_all','KIDNEY','CREATININE','ALBUMIN',
               'OS')

all_pvals <- c()
for (i in 1:length(factors_test)) {
  f <- factors_test[i]
  mv <- meta_vars[i]
  tmp <- cbind.data.frame(pbmc_container$donor_metadata[,mv],pbmc_container$projected_scores[rownames(pbmc_container$donor_metadata),f])
  colnames(tmp) <- c('cvar','dsc')

  if (is.factor(tmp$cvar)) {
    fmod <- glm(cvar~dsc, data=tmp, family = "binomial") ##"full" mod
    nmod <- glm(cvar~1, data=tmp, family = 'binomial') ##"null" mod
    a_res <- anova(nmod, fmod, test = 'Chisq')
    pval <- a_res$`Pr(>Chi)`[2]
  } else {
    x <- summary(lm(cvar~dsc,data = tmp))
    pval <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)
  }
  all_pvals <- c(all_pvals,pval)

}
p_adj <- p.adjust(all_pvals,method='fdr')
print(all_pvals)
print(p_adj)

ndx_order <- order(p_adj,decreasing=FALSE)

print(factors_test[ndx_order])
print(meta_vars[ndx_order])
print(all_pvals[ndx_order])
print(p_adj[ndx_order])


# computing the meta analysis pvalues with their corresponding trait p-values in the main analysis
library(metap)
# not shown - rerunning stuff from aSLE clinical association tests to get the original raw p-values
# prednisone and OS
p1 <- sumlog(c(1.738293e-08,3.147647e-02))$p # = 1.22159e-08
# sledaiscore and SLEDAI
p2 <- sumlog(c(0.0003896732,3.022751e-02))$p # = 0.0001454594
# nephritis incidence and creatinine
p3 <- sumlog(c(.01,3.515598e-02))$p # =  0.003147561
# nephritis incidence and neph_all
p4 <- sumlog(c(.01,4.627291e-02))$p # =  0.004015734
# lymphocyte % and LYM_ABS
p5 <- sumlog(c(1.541e-05,8.139915e-05))$p # = 2.696454e-08
# acrantidsdna and dsDNA
p6 <- sumlog(c(0.000342,.3383254))$p # = .001

# manually doing FDR correction, taking into account all tests but applied just to these 5
p4_adj <- p4 * length(factors_test) / 6 # 6th smallest nominal pval
print(p4_adj)
p3_adj <- p3 * length(factors_test) / 5 # 5th smallest nominal pval
print(p3_adj)
p6_adj <- p6 * length(factors_test) / 4 # 4th smallest nominal pval
print(p6_adj)
p2_adj <- p2 * length(factors_test) / 3 # 3rd smallest nominal pval
print(p2_adj)
p5_adj <- p5 * length(factors_test) / 2 # 2nd smallest nominal pval
print(p5_adj)
p1_adj <- p1 * length(factors_test) / 1 # 1st smallest nominal pval
print(p1_adj)
# all are monotonically decreasing




############### for the above, need to recompute association statistic just for naive Th cells ########
subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds') # from aSLE_decomp_associations.R file
subclusts <- subclusters[["T4"]][["res:0.6"]]
c_keep <- colnames(pbmc_container_SLE[["scMinimal_ctype"]][["Th"]][["count_data"]])
subclusts <- subclusts[c_keep]
subclusts <- factor(subclusts,levels=unique(subclusts))
sub_meta_tmp <- pbmc_container_SLE$scMinimal_ctype[['Th']]$metadata
donor_props <- compute_donor_props(subclusts,sub_meta_tmp)
dsc <- pbmc_container_SLE$tucker_results[[1]][,1,drop=F]
tmp <- cbind.data.frame(donor_props[rownames(dsc),1],dsc)
colnames(tmp) <- c('props','dscore')
plot(tmp$dscore,tmp$props)
summary(lm(tmp$dscore~tmp$props))







## plotting some of these results

## LYM_ABS
tmp <- pbmc_container$projected_scores[,1]
meta_var <- 'LYM_ABS'
tmp2 <- cbind.data.frame(pbmc_container$donor_metadata[,meta_var],
                         tmp[rownames(pbmc_container$donor_metadata)])
colnames(tmp2) <- c('cvar','dsc')

# adding best fit line
lmres <- lm(cvar~dsc,data=tmp2)
line_range <- seq(min(tmp2$dsc),max(tmp2$dsc),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

p1 <- ggplot(tmp2,aes(x=dsc,y=cvar)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Projected Factor 1 Scores') +
  ylab('Lymphocyte Count') +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/proj_F1_lymph_abs.pdf", useDingbats = FALSE,
#     width = 4, height = 3.25)
p1
# dev.off()


## SLEDAI
tmp <- pbmc_container$projected_scores[,1]
meta_var <- 'SLEDAI'
tmp2 <- cbind.data.frame(pbmc_container$donor_metadata[,meta_var],
                         tmp[rownames(pbmc_container$donor_metadata)])
colnames(tmp2) <- c('cvar','dsc')

# adding best fit line
lmres <- lm(cvar~dsc,data=tmp2)
line_range <- seq(min(tmp2$dsc),max(tmp2$dsc),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

p2 <- ggplot(tmp2,aes(x=dsc,y=cvar)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Projected Factor 1 Scores') +
  ylab('SLEDAI Score') +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/proj_F1_SLEDAI.pdf", useDingbats = FALSE,
#     width = 4, height = 3.25)
p2
# dev.off()



## CREATININE
tmp <- pbmc_container$projected_scores[,2]
meta_var <- 'CREATININE'
tmp2 <- cbind.data.frame(pbmc_container$donor_metadata[,meta_var],
                         tmp[rownames(pbmc_container$donor_metadata)])
colnames(tmp2) <- c('cvar','dsc')

# adding best fit line
lmres <- lm(cvar~dsc,data=tmp2)
line_range <- seq(min(tmp2$dsc),max(tmp2$dsc),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

p3 <- ggplot(tmp2,aes(x=dsc,y=cvar)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Projected Factor 2 Scores') +
  ylab('Creatinine') +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))

### Figure S4b
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/proj_F2_CREATININE.pdf", useDingbats = FALSE,
#     width = 4, height = 3.25)
p3
# dev.off()




## Neph_all
tmp <- pbmc_container$projected_scores[,2]
meta_var <- 'Neph_all'
tmp2 <- cbind.data.frame(pbmc_container$donor_metadata[,meta_var],
                         tmp[rownames(pbmc_container$donor_metadata)])
colnames(tmp2) <- c('cvar','dsc')

# remove na donors
tmp2 <- tmp2[!is.na(tmp2$cvar),]

levels(tmp2$cvar) <- c('No nephritis','Nephritis','ND')
tmp2$cvar <- factor(tmp2$cvar,levels=c('Nephritis','No nephritis')) # reverse order of levels

p4 <- ggplot(tmp2,aes(x=cvar,y=dsc)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5, binwidth = .015) +
  ylab('Projected Factor 2 Scores') +
  xlab('') +
  coord_flip() +
  theme_bw()

### Figure S4c
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/proj_F2_Neph_all.pdf", useDingbats = FALSE,
#     width = 4, height = 3.25)
p4
# dev.off()





## OS
tmp <- pbmc_container$projected_scores[,3]
meta_var <- 'OS'
tmp2 <- cbind.data.frame(pbmc_container$donor_metadata[,meta_var],
                         tmp[rownames(pbmc_container$donor_metadata)])
colnames(tmp2) <- c('cvar','dsc')

# remove na donors
tmp2 <- tmp2[!is.na(tmp2$cvar),]

levels(tmp2$cvar) <- c('No oral steroids','Oral steroids','ND')

p5 <- ggplot(tmp2,aes(x=cvar,y=dsc)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5, binwidth = .015) +
  ylab('Projected Factor 3 Scores') +
  xlab("") +
  coord_flip() +
  theme_bw()

### Figure S5b
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/proj_F3_OS.pdf", useDingbats = FALSE,
#     width = 4, height = 3.25)
p5
# dev.off()








####### storing the effect sizes, errors and fisher meta p-values for bar plots
res_plot_lym <- matrix(nrow=3,ncol=4)
colnames(res_plot_lym) <- c('study','beta','beta_se','adj_pval')

res_plot_sledai <- matrix(nrow=2,ncol=4)
colnames(res_plot_sledai) <- c('study','beta','beta_se','adj_pval')

res_plot_dsdna <- matrix(nrow=2,ncol=4)
colnames(res_plot_dsdna) <- c('study','beta','beta_se','adj_pval')


## dsdna
tmp <- pbmc_container$projected_scores[,1]
meta_var <- 'dsDNA'
tmp2 <- cbind.data.frame(pbmc_container$donor_metadata[,meta_var],
                         tmp[rownames(pbmc_container$donor_metadata)])
colnames(tmp2) <- c('cvar','dsc')
tmp2 <- tmp2[!is.na(tmp2$cvar),]
tmp2[,'dsc'] <- scale(tmp2[,'dsc'])
tmp2 <- as.data.frame(tmp2)
fmod <- glm(cvar~dsc, data=tmp2, family = "binomial") ##"full" mod
lmres <- summary(fmod)
ndx_retrieve <- which(meta_vars==meta_var)

res_plot_dsdna[1,] <- c('pediatric',lmres$coefficients['dsc','Estimate'],
                        lmres$coefficients['dsc','Std. Error'],
                        p_adj[ndx_retrieve])

# now to recompute betas from main data
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_categorical.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL

# get tucker donor scores to test
dsc <- pbmc_container_SLE$tucker_results[[1]]

trim_names <- sapply(rownames(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
rownames(dsc) <- trim_names
# get donors in both dataframes
d_both <- rownames(clin_vars)[rownames(clin_vars) %in% rownames(dsc)]
dsc <- dsc[d_both,]
clin_vars <- clin_vars[d_both,]
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,'acrantidsdna'])]
tmp <- as.data.frame(cbind(dsc[d_keep,1], clin_vars[d_keep,'acrantidsdna']))
colnames(tmp) <- c('dscore','cvar')
tmp[,'dscore'] <- scale(tmp[,'dscore'])
tmp <- as.data.frame(tmp)
# force cvar to be factor
tmp$cvar <- as.factor(tmp$cvar)
fmod <- glm(cvar~dscore, data=tmp, family = "binomial") ##"full" mod
lmres <- summary(fmod)
res_plot_dsdna[2,] <- c('adult',lmres$coefficients['dscore','Estimate'],
                        lmres$coefficients['dscore','Std. Error'],
                        .012) # padj from figure_2_code_v2



## LYM_ABS
tmp <- pbmc_container$projected_scores[,1]
meta_var <- 'LYM_ABS'
tmp2 <- cbind.data.frame(pbmc_container$donor_metadata[,meta_var],
                         tmp[rownames(pbmc_container$donor_metadata)])
colnames(tmp2) <- c('cvar','dsc')
tmp2 <- tmp2[!is.na(tmp2$cvar),]
tmp2 <- as.data.frame(scale(tmp2))
lmres <- lm(cvar~dsc,data=tmp2)
lmres <- summary(lmres)
ndx_retrieve <- which(meta_vars==meta_var)

res_plot_lym[1,] <- c('pediatric',lmres$coefficients['dsc','Estimate'],
                      lmres$coefficients['dsc','Std. Error'],
                      p_adj[ndx_retrieve])

### get original adult SLE effect size and significance for this association
# need to make sure the full data is limited to the cells used in analysis
all_cells <- c()
for (ct in pbmc_container_SLE$experiment_params$ctypes_use) {
  cells_in_ctype <- rownames(pbmc_container_SLE$scMinimal_ctype[[ct]]$metadata)
  all_cells <- c(all_cells,cells_in_ctype)
}

pbmc_container_SLE$scMinimal_full$metadata <- pbmc_container_SLE$scMinimal_full$metadata[all_cells,]
pbmc_container_SLE$scMinimal_full$count_data <- pbmc_container_SLE$scMinimal_full$count_data[,all_cells]

scMinimal <- pbmc_container_SLE$scMinimal_full
donor_scores <- pbmc_container_SLE$tucker_results[[1]]
metadata <- scMinimal$metadata

# map cell types to numbers temporarily
all_ctypes <- unique(as.character(metadata$ctypes)) # index of this is the mapping
cell_clusters <- sapply(as.character(metadata$ctypes),function(x){
  return(which(all_ctypes %in% x))
})
names(cell_clusters) <- rownames(metadata)

# get matrix of donor proportions of different cell types
donor_props <- compute_donor_props(cell_clusters,metadata)
colnames(donor_props) <- pbmc_container_SLE$experiment_params$ctypes_use
# donor_props <- cbind(donor_props,rowSums(donor_props[,c('Th','Tc','B','NK')]))
# colnames(donor_props)[ncol(donor_props)] <- c('lym')
dsc <- get_one_factor(pbmc_container_SLE,1)[[1]]
tmp <- cbind.data.frame(dsc,donor_props[rownames(dsc),'Th'])
colnames(tmp)[2] <- c('Th')

### plotting this result
# adding best fit line
lmres <- lm(Th~dsc,data=tmp)
line_range <- seq(min(tmp$dsc),max(tmp$dsc),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

p <- ggplot(tmp,aes(x=dsc,y=Th)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Factor 1 Scores') +
  ylab('Overall Th Fractions') +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))

# now adding the beta and significance
tmp2 <- as.data.frame(scale(tmp))
lmres <- summary(lm(Th~dsc,data=tmp2))

res_plot_lym[2,] <- c('adult',lmres$coefficients['dsc','Estimate'],
                      lmres$coefficients['dsc','Std. Error'],
                      lmres$coefficients['dsc','Pr(>|t|)'])

###

all_cells <- c()
for (ct in pbmc_container$experiment_params$ctypes_use) {
  cells_in_ctype <- rownames(pbmc_container$scMinimal_ctype[[ct]]$metadata)
  all_cells <- c(all_cells,cells_in_ctype)
}

pbmc_container$scMinimal_full$metadata <- pbmc_container$scMinimal_full$metadata[all_cells,]
pbmc_container$scMinimal_full$count_data <- pbmc_container$scMinimal_full$count_data[,all_cells]

scMinimal <- pbmc_container$scMinimal_full
donor_scores <- pbmc_container$tucker_results[[1]]
metadata <- scMinimal$metadata

# map cell types to numbers temporarily
all_ctypes <- unique(as.character(metadata$ctypes)) # index of this is the mapping
cell_clusters <- sapply(as.character(metadata$ctypes),function(x){
  return(which(all_ctypes %in% x))
})
names(cell_clusters) <- rownames(metadata)

# get matrix of donor proportions of different cell types
donor_props <- compute_donor_props(cell_clusters,metadata)
colnames(donor_props) <- pbmc_container$experiment_params$ctypes_use
# donor_props <- cbind(donor_props,rowSums(donor_props[,c('Th','Tc','B','NK')]))
# colnames(donor_props)[ncol(donor_props)] <- c('lym')
dsc <- pbmc_container$projected_scores[,1,drop=F]
tmp <- cbind.data.frame(dsc,donor_props[rownames(dsc),'Th'])
colnames(tmp)[2] <- c('Th')

# now adding the beta and significance
tmp2 <- as.data.frame(scale(tmp))
lmres <- summary(lm(Th~dsc,data=tmp2))

res_plot_lym[3,] <- c('pediatric',lmres$coefficients['dsc','Estimate'],
                      lmres$coefficients['dsc','Std. Error'],
                      lmres$coefficients['dsc','Pr(>|t|)'])


## SLEDAI
tmp <- pbmc_container$projected_scores[,1]
meta_var <- 'SLEDAI'
tmp2 <- cbind.data.frame(pbmc_container$donor_metadata[,meta_var],
                         tmp[rownames(pbmc_container$donor_metadata)])
colnames(tmp2) <- c('cvar','dsc')
tmp2 <- tmp2[!is.na(tmp2$cvar),]
tmp2 <- as.data.frame(scale(tmp2))
lmres <- lm(cvar~dsc,data=tmp2)
lmres <- summary(lmres)
ndx_retrieve <- which(meta_vars==meta_var)

res_plot_sledai[1,] <- c('pediatric',lmres$coefficients['dsc','Estimate'],
                         lmres$coefficients['dsc','Std. Error'],
                         p_adj[ndx_retrieve])

# recompute adult SLE SLEDIA effect size and such
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_ordinal.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid

# get tucker donor scores to test
dsc <- pbmc_container_SLE$tucker_results[[1]]

## get donors in both dsc and in clin_vars
# trim donor IDs in dsc
trim_names <- sapply(rownames(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()

clin_var_sub <- clin_vars[,'sledaiscore',drop=F]
d_both <- intersect(trim_names,rownames(clin_var_sub))
tmp <- cbind.data.frame(dsc[,1],clin_var_sub[trim_names,])
colnames(tmp) <- c('dsc','SLEDAI')
tmp <- as.data.frame(scale(tmp))
lmres <- summary(lm(SLEDAI~dsc,data=tmp))

res_plot_sledai[2,] <- c('adult',lmres$coefficients['dsc','Estimate'],
                         lmres$coefficients['dsc','Std. Error'],
                         .008) # padj from figure_2_covde_v2

## creating the full bar plots
res_plot_sledai <- as.data.frame(res_plot_sledai)
res_plot_sledai$beta <- as.numeric(res_plot_sledai$beta)
res_plot_sledai$beta_se <- as.numeric(res_plot_sledai$beta_se)
res_plot_sledai$adj_pval <- as.numeric(res_plot_sledai$adj_pval)
res_plot_dsdna <- as.data.frame(res_plot_dsdna)
res_plot_dsdna$beta <- as.numeric(res_plot_dsdna$beta)
res_plot_dsdna$beta_se <- as.numeric(res_plot_dsdna$beta_se)
res_plot_dsdna$adj_pval <- as.numeric(res_plot_dsdna$adj_pval)
res_plot_lym <- as.data.frame(res_plot_lym)
res_plot_lym$beta <- as.numeric(res_plot_lym$beta)
res_plot_lym$beta_se <- as.numeric(res_plot_lym$beta_se)
res_plot_lym$adj_pval <- as.numeric(res_plot_lym$adj_pval)

res_plot <- rbind.data.frame(res_plot_sledai,res_plot_dsdna,res_plot_lym)
res_plot$var_lab <- c('sledai','sledai','dsDNA','dsDNA','lym_abs','lym','lym')

res_plot$var_lab <- factor(res_plot$var_lab,levels=c('sledai','dsDNA','lym','lym_abs'))

res_plot2 <- res_plot
res_plot2 <- rbind.data.frame(res_plot2,c('adult',0,0,1,'lym_abs'))
res_plot2$beta <- as.numeric(res_plot2$beta)
res_plot2$beta_se <- as.numeric(res_plot2$beta_se)
res_plot2$adj_pval <- as.numeric(res_plot2$adj_pval)

fig <- ggplot(res_plot2, aes(x=var_lab, y=beta, fill=study)) +
  geom_bar(stat="identity", position=position_dodge(),width=.75) +
  geom_errorbar(aes(ymin=beta-beta_se, ymax=beta+beta_se), width=.25,
                position=position_dodge(.75)) +
  xlab('') +
  # ylim(c(-1,.8)) +
  geom_hline(yintercept=0) +
  theme_classic()

### Figure S3f
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/F1_meta_barplot2.pdf", useDingbats = FALSE,
#     width = 7.5, height = 4)
fig
# dev.off()























### comparing projected scores to the pediatric lupus decomposition dscores

# reversing the signs of factors so that they match up with those from the main sle analysis
pbmc_container$tucker_results[[1]] <- pbmc_container$tucker_results[[1]] * -1
pbmc_container$tucker_results[[2]] <- pbmc_container$tucker_results[[2]] * -1

cor_mat <- cor(pbmc_container[["projected_scores"]],pbmc_container$tucker_results[[1]])
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

### Figure S3e
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/pediatric_comp_proj_decomp.pdf", useDingbats = FALSE,
#     width = 4.5, height = 5.5)
cor_hmap
# dev.off()











##### checking that the LR interactions are also replicated...
## first rerunning LR analysis for my main SLE dataset
pbmc_container_SLE <- form_tensor(pbmc_container_SLE, donor_min_cells=20,
                                  norm_method='trim', scale_factor=10000,
                                  vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                                  scale_var = TRUE, var_scale_power = .5,
                                  batch_var='pool')

# using cellchat database
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
  rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
  return(rname)
})
lr_pairs$interaction_name <- NULL

# infer active LR interactions
pbmc_container_SLE <- prep_LR_interact(pbmc_container_SLE, lr_pairs, norm_method='trim', scale_factor=10000,
                                       var_scale_power=.5, batch_var='pool')
sft_thresh <- c(12,14,12,10,12,9,12)
pbmc_container_SLE <- get_gene_modules(pbmc_container_SLE,sft_thresh)

# reshowing original highlighted LR hits
lig_mod_fact <- plot_mod_and_lig(pbmc_container_SLE,factor_select=2,mod_ct='Th',mod=5,lig_ct='cMono',lig='ICOSLG')
lig_mod_fact

lig_mod_fact <- plot_mod_and_lig(pbmc_container_SLE,factor_select=3,mod_ct='Th',mod=9,lig_ct='cDC',lig='THBS1')
lig_mod_fact


# reform pediatric tensor with all genes
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=2, 
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=1,
                              batch_var = 'Batch',
                              scale_var = TRUE, var_scale_power = .5)

# Look at ICOSLG interaction
# getting genes in the module previously used
th_modules <- pbmc_container_SLE[["module_genes"]][["Th"]]
ndx_keep <- which(th_modules==5)
mymod <- names(th_modules)[ndx_keep]

pb <- pbmc_container[["scMinimal_ctype"]][["Th"]][["pseudobulk"]]
mod_pres <- mymod[mymod %in% colnames(pb)]
pb_mod <- pb[,mod_pres]
eg <- prcomp(pb_mod,scale. = FALSE)
eg <- eg$x[,1,drop=FALSE]

lig_exp <- pbmc_container[["scMinimal_ctype"]][["cMono"]][["pseudobulk"]][,'ICOSLG']

dsc <- pbmc_container$tucker_results[[1]][,2]
tmp <- cbind.data.frame(dsc,eg[names(dsc),1])
colnames(tmp) <- c('l_exp','mod_vals')
cor(tmp)
# eigengene pc is anti-correlated with the factor, so we'll just flip it's arbitrary direction
eg <- eg * -1

# now to plot module vs ligand expression
tmp <- cbind.data.frame(lig_exp,eg[names(lig_exp),1])
colnames(tmp) <- c('lig_exp','ME')
cor(tmp)

lmres <- lm(ME~lig_exp,data=tmp)
line_range <- seq(min(tmp$lig_exp),max(tmp$lig_exp),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')
x <- summary(lmres)
pval <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)

p3 <- ggplot(tmp,aes(x=lig_exp,y=ME)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab(paste0('ICOSLG expression in cMono')) +
  ylab(paste0('Th_m5 module expression')) +
  annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
           label=paste0('p-value = ',
                        as.character(signif(as.numeric(format(pval,scientific=TRUE)),digits=2)))) +
  theme_bw()
p3

# now to plot ligand expression versus factor 2 expression
tmp <- cbind.data.frame(lig_exp,dsc[names(lig_exp)])
colnames(tmp) <- c('l_exp','dscores')
cor(tmp)

lmres <- lm(l_exp~dscores,data=tmp)
line_range <- seq(min(tmp$dscores),max(tmp$dscores),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')
x <- summary(lmres)
pval <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)

p1 <- ggplot(tmp,aes(x=dscores,y=l_exp)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab(paste0('Factor 2 donor score')) +
  ylab(paste0('ICOSLG expression in cMono')) +
  annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
           label=paste0('p-value = ',
                        as.character(signif(as.numeric(format(pval,scientific=TRUE)),digits=2)))) +
  theme_bw()
p1

fig <- plot_grid(p3,p1)

### Figure S4e
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/ICOSLG_ped_rep.pdf", useDingbats = FALSE,
#     width = 5.75, height = 2.5)
fig
# dev.off()



# THBS1 interaction

# need to reform tensor to include dendritic cells
pbmc_container[["experiment_params"]][["ctypes_use"]] <- c("B","NK","Th","Tc", 
                                                           "cMono","ncMono","cDC")
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=2, 
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=1,
                              batch_var = 'Batch',
                              scale_var = TRUE, var_scale_power = .5)

th_modules <- pbmc_container_SLE[["module_genes"]][["Th"]]
ndx_keep <- which(th_modules==9)
mymod <- names(th_modules)[ndx_keep]

pb <- pbmc_container[["scMinimal_ctype"]][["Th"]][["pseudobulk"]]
mod_pres <- mymod[mymod %in% colnames(pb)]
pb_mod <- pb[,mod_pres]
eg <- prcomp(pb_mod,scale. = FALSE)
eg <- eg$x[,1,drop=FALSE]

lig_exp <- pbmc_container[["scMinimal_ctype"]][["cDC"]][["pseudobulk"]][,'THBS1']

dsc <- pbmc_container$tucker_results[[1]][,3]
tmp <- cbind.data.frame(dsc[rownames(eg)],eg[,1])
colnames(tmp) <- c('l_exp','mod_vals')
cor(tmp)

# now to plot module vs ligand expression
tmp <- cbind.data.frame(lig_exp,eg[names(lig_exp),1])
colnames(tmp) <- c('lig_exp','ME')
cor(tmp)

lmres <- lm(ME~lig_exp,data=tmp)
line_range <- seq(min(tmp$lig_exp),max(tmp$lig_exp),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')
x <- summary(lmres)
pval <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)

p3 <- ggplot(tmp,aes(x=lig_exp,y=ME)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab(paste0('THBS1 expression in cDC')) +
  ylab(paste0('Th_m9 module expression')) +
  annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
           label=paste0('p-value = ',
                        as.character(signif(as.numeric(format(pval,scientific=TRUE)),digits=2)))) +
  theme_bw()
p3

# now to plot ligand expression versus factor 2 expression
tmp <- cbind.data.frame(lig_exp,dsc[names(lig_exp)])
colnames(tmp) <- c('l_exp','dscores')
cor(tmp)

lmres <- lm(l_exp~dscores,data=tmp)
line_range <- seq(min(tmp$dscores),max(tmp$dscores),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')
x <- summary(lmres)
pval <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)

p1 <- ggplot(tmp,aes(x=dscores,y=l_exp)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab(paste0('Factor 3 donor score')) +
  ylab(paste0('THBS1 expression in cDC')) +
  annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
           label=paste0('p-value = ',
                        as.character(signif(as.numeric(format(pval,scientific=TRUE)),digits=2)))) +
  theme_bw()
p1

fig <- plot_grid(p3,p1)

### Figure S5d
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/THBS1_ped_rep.pdf", useDingbats = FALSE,
#     width = 5.75, height = 2.5)
fig
# dev.off()

