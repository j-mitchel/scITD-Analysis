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

### Figure S2G
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

pbmc_container <- get_donor_meta(pbmc_container, additional_meta = meta_cols)


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

### Figure 3C
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

### Figure S3A
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

### Figure 4C
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

### Figure 2F
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

### Figure S2H
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/pediatric_comp_proj_decomp.pdf", useDingbats = FALSE,
#     width = 4.5, height = 5.5)
cor_hmap
# dev.off()









####### now doing gsea of genes that are similar/different in each decomposition
# need to decide whether to take intersection of genes from both decompositions
# if a gene is DE in one dataset but not another, it would likely be a HVG in one dataset but not another
# maybe the best way to do this is just to take the slope of the lm coefficient as the "statistic"
# this way, I am able to include all genes (or at least all genes expressed above some threshold)

get_betas <- function(container, all_f_test, donor_min_cells=2,
                      norm_method='trim', scale_factor=10000,
                      batch_var = 'Batch', get_pvals=FALSE, add_error=FALSE) {
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
  total_res_pv <- list() # stores the results for all factors
  total_res_error <- list()
  for (f_test in all_f_test) {
    print(f_test)
    dsc <- container$tucker_results[[1]][,f_test,drop=FALSE]

    f_res <- list() # stores the result from a single factor
    f_res_pv <- list() # stores the result from a single factor
    f_res_error <- list()
    
    # loop through cell types
    for (ct in container$experiment_params$ctypes_use) {
      print(ct)
      ct_res <- c()
      ct_res_pv <- c()
      ct_res_error <- c()
      pb <- container$scMinimal_ctype[[ct]]$pseudobulk
      # loop through genes
      for (g_ndx in 1:ncol(pb)) {
        tmp <- cbind.data.frame(dsc,pb[rownames(dsc),g_ndx])
        colnames(tmp) <- c('dscore','expr')
        lmres <- summary(lm(expr~dscore,data=tmp))
        beta <- lmres$coefficients['dscore','Estimate']
        ct_res <- c(ct_res,beta)
        
        bse <- lmres$coefficients['dscore','Std. Error']
        ct_res_error <- c(ct_res_error,bse)
        
        pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
        ct_res_pv <- c(ct_res_pv,pval)
        
      }
      names(ct_res) <- colnames(pb)
      names(ct_res_error) <- colnames(pb)
      f_res[[ct]] <- ct_res
      f_res_error[[ct]] <- ct_res_error
      f_res_pv[[ct]] <- p.adjust(ct_res_pv,method='fdr')
    }
    total_res[[f_test]] <- f_res
    total_res_error[[f_test]] <- f_res_error
    total_res_pv[[f_test]] <- f_res_pv
  }

  if (get_pvals) {
    return(list(total_res,total_res_pv))
  } else if (add_error) {
    return(list(total_res,total_res_error))
  } else {
    return(total_res)
  }
}

betas_ped <- get_betas(pbmc_container, all_f_test=c(1,2,3), donor_min_cells=2,
          norm_method='trim', scale_factor=10000,
          batch_var = 'Batch')

betas_adult <- get_betas(pbmc_container_SLE, all_f_test=c(1,2,3), donor_min_cells=20,
          norm_method='trim', scale_factor=10000,
          batch_var='pool')

betas_all <- list(betas_ped,betas_adult)
# saveRDS(betas_all,file='/home/jmitchel/data/lupus_data/ped_adult_betas.rds')




# function to get differences in betas, run gsea on these, and plot the result
get_diff_gsea <- function(betas_ped,betas_adult,cao,f_use1,f_use2,plot=TRUE,plot_n=25) {
  # limit to intersection of genes we have betas for in both datasets
  for (ct in names(betas_ped[[f_use1]])) {
    genes1 <- names(betas_ped[[f_use1]][[ct]])
    genes2 <- names(betas_adult[[f_use2]][[ct]])
    gene_intersect <- intersect(genes1,genes2)
    betas_ped[[f_use1]][[ct]] <- betas_ped[[f_use1]][[ct]][gene_intersect]
    betas_adult[[f_use2]][[ct]] <- betas_adult[[f_use2]][[ct]][gene_intersect]
  }

  # for each cell type take the difference across datasets
  betas_diff <- list()
  for (ct in names(betas_ped[[f_use1]])) {
    betas_diff[[ct]] <- betas_ped[[f_use1]][[ct]] - betas_adult[[f_use2]][[ct]]
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

# reload up the betas generated above
betas_all <- readRDS(file='/home/jmitchel/data/lupus_data/ped_adult_betas.rds')

betas_ped <- betas_all[[1]]
betas_adult <- betas_all[[2]]

# this analysis needs to be done per factor
p1 <- get_diff_gsea(betas_ped,betas_adult,cao,f_use1=1,f_use2=1,plot_n=15)
p2 <- get_diff_gsea(betas_ped,betas_adult,cao,f_use=2,f_use2=2,plot_n=5)
p3 <- get_diff_gsea(betas_ped,betas_adult,cao,f_use=3,f_use2=3,plot_n=5)

### Figure 5C
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/covid_sle_f1_compare_gsea.pdf", useDingbats = FALSE,
#     width = 6, height = 4)
p1
# dev.off()

### Figure 5E top
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/covid_sle_f2_compare_gsea.pdf", useDingbats = FALSE,
#     width = 6, height = 4)
p2
# dev.off()

### Figure 5E bottom
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/covid_sle_f3_compare_gsea.pdf", useDingbats = FALSE,
#     width = 6, height = 4)
p3
# dev.off()

library(cowplot)
fig <- cowplot::plot_grid(p1,p2,p3,ncol=2)
fig




##### plotting betas of some leading edge genes with bars

# first need to rerun get_betas but with error vals
betas_ped <- get_betas(pbmc_container, all_f_test=c(1,2,3), donor_min_cells=2,
                       norm_method='trim', scale_factor=10000,
                       batch_var = 'Batch', add_error=TRUE)

betas_adult <- get_betas(pbmc_container_SLE, all_f_test=c(1,2,3), donor_min_cells=20,
                         norm_method='trim', scale_factor=10000,
                         batch_var='pool', add_error=TRUE)

betas_all <- list(betas_ped,betas_adult)
# saveRDS(betas_all,file='/home/jmitchel/data/lupus_data/ped_adult_betas_errors.rds')

genes_plot <- c('BATF','CD22','CXCR5','PRF1','IL4R','SOX4','CCNA2','KIFC1')
cell_types <- c('B','B','B','Tc','Tc','Tc','cMono','cMono')
res <- data.frame(matrix(nrow=length(genes_plot),ncol=3))
rownames(res) <- genes_plot
colnames(res) <- c('gene','beta_ped','beta_adult')
res_se <- res
colnames(res_se) <- c('gene','se_ped','se_adult')
for (i in 1:length(genes_plot)) {
  gene <- genes_plot[i]
  ct <- cell_types[i]
  b_ped <- betas_ped[[1]][[1]][[ct]][[gene]]
  b_adult <- betas_adult[[1]][[1]][[ct]][[gene]]
  res[gene,] <- c(gene,b_ped,b_adult)
  
  ###### NEED TO ADD [[1]] WHERE RELEVANT
  bse_ped <- betas_ped[[2]][[1]][[ct]][[gene]]
  bse_adult <- betas_adult[[2]][[1]][[ct]][[gene]]
  res_se[gene,] <- c(gene,bse_ped,bse_adult)
}

res2 <- melt(res, id = c('gene')) 
colnames(res2)[2:3] <- c('dataset','beta')
res2$beta <- as.numeric(res2$beta)
res2$gene <- factor(genes_plot,levels = genes_plot)

res2_se <- melt(res_se, id = c('gene')) 
colnames(res2_se)[2:3] <- c('dataset','beta_se')
res2_se$beta_se <- as.numeric(res2_se$beta_se)
res2_se$gene <- factor(genes_plot,levels = genes_plot)

res2 <- cbind.data.frame(res2,res2_se)
res2[,4:5] <- NULL

p1 <- ggplot(res2,aes(x=gene,y=beta,fill=dataset)) +
  geom_bar(position="dodge", stat="identity", width = 0.75) +
  geom_errorbar(aes(ymin=beta-beta_se, ymax=beta+beta_se), width=.2,position=position_dodge(.75)) +
  geom_hline(yintercept = 0) +
  coord_flip() +
  theme_minimal()

### Figure 5D
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/covid_sle_f1_select_genes2.pdf", useDingbats = FALSE,
#     width = 6, height = 5)
p1
# dev.off()





### negative control analysis
n_iter <- 30
myres <- lapply(1:n_iter,function(x) {
  d_all <- unique(as.character(pbmc_sle@meta.data$ind_cov))
  down_samp_n <- round(.9 * length(d_all))
  d_keep <- sample(d_all,down_samp_n)
  cells_keep <- rownames(pbmc_sle@meta.data)[as.character(pbmc_sle@meta.data$ind_cov) %in% d_keep]
  pbmc_sle_sub <- subset(pbmc_sle,cells = cells_keep)
  
  param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc","cDC",
                                                 "cMono","ncMono"),
                                  ncores = 30, rand_seed = 10)
  
  pbmc_container_SLE2 <- make_new_container(seurat_obj=pbmc_sle_sub,
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
  
  
  pbmc_container_SLE2 <- form_tensor(pbmc_container_SLE2, donor_min_cells=20,
                                     norm_method='trim', scale_factor=10000,
                                     vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                                     scale_var = TRUE, var_scale_power = .5,
                                     batch_var='pool')
  
  
  pbmc_container_SLE2 <- run_tucker_ica(pbmc_container_SLE2, ranks=c(7,20),
                                        tucker_type = 'regular', rotation_type = 'hybrid')
  
  # test if need to flip sign
  t1 <- get_one_factor(pbmc_container_SLE,1)[[1]]
  t2 <- get_one_factor(pbmc_container_SLE2,1)[[1]]
  mycor <- cor(t2,t1[rownames(t2),])
  
  if (abs(mycor)<.9) {
    return(NA)
  } else if (mycor < -.9) {
    # flip F1 signs
    pbmc_container_SLE2$tucker_results[[1]][,1] <- pbmc_container_SLE2$tucker_results[[1]][,1] * -1
    pbmc_container_SLE2$tucker_results[[2]][1,] <- pbmc_container_SLE2$tucker_results[[2]][1,] * -1
  }
  
  
  betas_adult2 <- get_betas(pbmc_container_SLE2, all_f_test=c(1), donor_min_cells=20,
                            norm_method='trim', scale_factor=10000,
                            batch_var='pool')
  
  cao_tmp <- get_diff_gsea(betas_adult,betas_adult2,cao,f_use1=1,f_use2=1,plot=FALSE)
  ctypes <- names(cao_tmp$test.results$GSEA$res)
  ct_res <- list()
  for (ct in ctypes) {
    n_sig_padj_01 <- sum(cao_tmp$test.results$GSEA$res[[ct]][['BP']]$p.adjust<.01)
    n_sig_padj_05 <- sum(cao_tmp$test.results$GSEA$res[[ct]][['BP']]$p.adjust<.05)
    n_sig_nom_01 <- sum(cao_tmp$test.results$GSEA$res[[ct]][['BP']]$pvalue<.01)
    n_sig_nom_05 <- sum(cao_tmp$test.results$GSEA$res[[ct]][['BP']]$pvalue<.05)
    n_total <- nrow(cao_tmp$test.results$GSEA$res[[ct]][['BP']])
    ct_res[[ct]] <- c(n_sig_padj_01,n_sig_padj_05,n_sig_nom_01,n_sig_nom_05,n_total)
  }
  return(list(cao_tmp,ct_res))
})

# saveRDS(myres,file='/home/jmitchel/data/lupus_data/adult_v_ped_f_compare_gsea_null.rds')

myres <- readRDS('/home/jmitchel/data/lupus_data/adult_v_ped_f_compare_gsea_null.rds') # generated above


# making two barplots
# - fraction of nominal p<.05, p<.01 and Padj<.01, Padj<.05 for each cell type

ctypes_all <- names(myres[[1]][[2]])

prep_dat1 <- as.data.frame(matrix(0,nrow=length(ctypes_all),ncol=4))
rownames(prep_dat1) <- ctypes_all
colnames(prep_dat1) <- c('sig_padj_01','sig_padj_05','sig_nom_01','sig_nom_05')

prep_dat2 <- as.data.frame(matrix(0,nrow=length(ctypes_all),ncol=4))
rownames(prep_dat2) <- ctypes_all
colnames(prep_dat2) <- c('padj_01_sd','padj_05_sd','p_01_sd','p_05_sd')

dat_list <- list()

for (i in 1:length(myres)) {
  iter_res <- myres[[i]]
  iter_res <- iter_res[[2]]
  
  # loop through cell types and append fractions of null positive results
  for (ct in names(iter_res)) {
    sig_padj_01 <- iter_res[[ct]][1] / iter_res[[ct]][5]
    sig_padj_05 <- iter_res[[ct]][2] / iter_res[[ct]][5]
    sig_nom_01 <- iter_res[[ct]][3] / iter_res[[ct]][5]
    sig_nom_05 <- iter_res[[ct]][4] / iter_res[[ct]][5]
    
    dat_list[[ct]][['sig_padj_01']] <- c(dat_list[[ct]][['sig_padj_01']],sig_padj_01)
    dat_list[[ct]][['sig_padj_05']] <- c(dat_list[[ct]][['sig_padj_05']],sig_padj_05)
    dat_list[[ct]][['sig_nom_01']] <- c(dat_list[[ct]][['sig_nom_01']],sig_nom_01)
    dat_list[[ct]][['sig_nom_05']] <- c(dat_list[[ct]][['sig_nom_05']],sig_nom_05)
    
    # prep_dat[ct,c(2:5)] <- prep_dat[ct,c(2:5)] + c(sig_padj_01,sig_padj_05,sig_nom_01,sig_nom_05)
  }
}

for (ct in ctypes_all) {
  ct_res <- dat_list[[ct]]
  for (i in 1:length(ct_res)) {
    my_sd <- sd(ct_res[[i]])
    my_mean <- mean(ct_res[[i]])
    prep_dat1[ct,i] <- my_mean 
    prep_dat2[ct,i] <- my_sd
  }
}

# append the two dataframes

prep_dat3 <- cbind.data.frame(prep_dat1,prep_dat2)
prep_dat3 <- cbind.data.frame(prep_dat3,ctypes_all)

# need to melt df by stacking the different thresholds
c1 <- c(prep_dat3$sig_padj_01,prep_dat3$sig_padj_05,prep_dat3$sig_nom_01,prep_dat3$sig_nom_05)
c2 <- c(prep_dat3$padj_01_sd,prep_dat3$padj_05_sd,prep_dat3$p_01_sd,prep_dat3$p_05_sd)
c3 <- rep(ctypes_all,4)
c4 <- c(rep('sig_padj_01',7),rep('sig_padj_05',7),rep('sig_nom_01',7),rep('sig_nom_05',7))
prep_dat4 <- cbind.data.frame(c1,c2,c3,c4)
colnames(prep_dat4) <- c('av','stdev','ctype','type')

# results columns: cell type, fraction p<.05, sd of fraction<.05, ...
p <- ggplot(prep_dat4,aes(x=as.factor(type),y=av,fill=ctype)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=av-stdev, ymax=av+stdev), width=.2,position=position_dodge(.9)) +
  theme_bw()

### Figure 5F
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/gsea_compare_beta_bplot.pdf", useDingbats = FALSE,
#     width = 8, height = 4)
p
# dev.off()






#### computing gene pattern similarity across factors from the two datasets
# need to recompute it but with get_pvals=TRUE
betas_ped <- get_betas(pbmc_container, all_f_test=c(1,2,3), donor_min_cells=2,
                       norm_method='trim', scale_factor=10000,
                       batch_var = 'Batch', get_pvals = TRUE)

betas_adult <- get_betas(pbmc_container_SLE, all_f_test=c(1,2,3), donor_min_cells=20,
                         norm_method='trim', scale_factor=10000,
                         batch_var='pool', get_pvals = TRUE)

beta_vals_ped <- betas_ped[[1]]
beta_vals_adult <- betas_adult[[1]]
p_vals_ped <- betas_ped[[2]]
p_vals_adult <- betas_adult[[2]]

ctypes_all <- c('B','NK','Th','Tc','cMono','ncMono')
cor_res <- matrix(nrow=length(beta_vals_ped),ncol=length(ctypes_all))
rownames(cor_res) <- c('aSLE_F1:pSLE_F1','aSLE_F2:pSLE_F2','aSLE_F3:pSLE_F3')
colnames(cor_res) <- ctypes_all
for (i in 1:length(beta_vals_adult)) {
  betas_adult_f <- beta_vals_adult[[i]]
  betas_ped_f <- beta_vals_ped[[i]]
  
  p_adult_f <- p_vals_adult[[i]]
  p_ped_f <- p_vals_ped[[i]]
  
  # limit all to the cell types found in both
  betas_adult_f <- betas_adult_f[ctypes_all]
  betas_ped_f <- betas_ped_f[ctypes_all]
  p_adult_f <- p_adult_f[ctypes_all]
  p_ped_f <- p_ped_f[ctypes_all]
  
  
  # loop through cell types
  for (ct in names(betas_adult_f)) {
    betas_adult_f_ct <- betas_adult_f[[ct]]
    betas_ped_f_ct <- betas_ped_f[[ct]]
    p_adult_f_ct <- p_adult_f[[ct]]
    p_ped_f_ct <- p_ped_f[[ct]]
    
    names(p_adult_f_ct) <- names(betas_adult_f_ct)
    names(p_ped_f_ct) <- names(betas_ped_f_ct)
    
    # reduce genes to the intersection of genes tested in both
    g_both <- intersect(names(betas_adult_f_ct),names(betas_ped_f_ct))
    betas_adult_f_ct <- betas_adult_f_ct[g_both]
    betas_ped_f_ct <- betas_ped_f_ct[g_both]
    p_adult_f_ct <- p_adult_f_ct[g_both]
    p_ped_f_ct <- p_ped_f_ct[g_both]
    
    # get union of significant genes
    g_keep1 <- names(p_adult_f_ct)[p_adult_f_ct<.01]
    g_keep2 <- names(p_ped_f_ct)[p_ped_f_ct<.01]
    g_keep <- unique(c(g_keep1,g_keep2))
    
    betas_adult_f_ct <- betas_adult_f_ct[g_keep]
    betas_ped_f_ct <- betas_ped_f_ct[g_keep]
    
    mycor <- cor(betas_adult_f_ct,betas_ped_f_ct)
    cor_res[i,ct] <- mycor
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
                  row_names_gp = gpar(fontsize = 12),
                  show_heatmap_legend = TRUE,
                  width = unit(h_w[2], "cm"),
                  height = unit(h_w[1], "cm"),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid::grid.text(sprintf("%.2f", as.matrix(cor_res)[i, j]), x, y, gp = gpar(fontsize = 10))
                  })

### Figure 5A
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/factor_assoc_cors_v2.pdf", useDingbats = FALSE,
#     width = 6, height = 4)
myhmap
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

### Figure S3D
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

### Figure 4F
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/THBS1_ped_rep.pdf", useDingbats = FALSE,
#     width = 5.75, height = 2.5)
fig
# dev.off()

