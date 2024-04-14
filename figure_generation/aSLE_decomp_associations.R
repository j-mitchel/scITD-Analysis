library(Seurat)
library(conos)
library(ggplot2)
library(coda.base)
library(RColorBrewer)
library(readxl)
library(MASS)
library(ggrastr)
library(devtools)
load_all('/home/jmitchel/scITD/')

# conos object from preprocessing/embedding_prep.R file
con <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_conos2.rds')
con[["embeddings"]][["UMAP"]] <- con[["embeddings"]][["UMAP"]][cells_keep,]
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

### Figure 2a
# pdf(file = "/home/jmitchel/figures/for_paper/lupus_embedding.pdf", useDingbats = FALSE,
#     width = 7, height = 7)
# saved a jpeg 550 x 400 dimensions
tmp
# dev.off()


## make umap showing subtype annotations
c_use <- names(con$clusters$leiden$groups)
con$clusters$leiden$groups <- pbmc@meta.data[c_use,'ct_cov']
names(con$clusters$leiden$groups) <- c_use
tmp <- con$plotGraph(alpha=0.1)
# mycolors <- brewer.pal(n = 9, name = "Set1")
mycolors <- c(brewer.pal(n = 5, name = "Set1"),brewer.pal(n = 8, name = "Set2"))
# mycolors <- mycolors[c(1:5,7,9)]
tmp <- tmp + theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black", size=1, fill=NA)) +
  scale_color_manual(values=mycolors) +
  scale_y_reverse() +
  scale_x_reverse()
# scale_color_brewer(palette="Set1")
tmp$layers[[2]] <- NULL

### Figure S2a
# saved a jpeg 550 x 400 dimensions
tmp
# dev.off()














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

pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(14,20),
                                         shuffle_level='cells', shuffle_within=NULL,
                                         num_iter=50, batch_var='pool',
                                         norm_method='trim',
                                         scale_factor=10000,
                                         scale_var=TRUE,
                                         var_scale_power=.5)

### Figure S2c
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/sle_only_rank_determination2.pdf", useDingbats = FALSE,
#     width = 9, height = 7)
pbmc_container$plots$rank_determination_plot
# dev.off()

pbmc_container <- get_min_sig_genes(pbmc_container, donor_rank_range=c(2:14), gene_ranks=20,
                                    use_lm=TRUE, tucker_type='regular',
                                    rotation_type='hybrid',
                                    n.cores = 5, thresh=0.05)

p <- pbmc_container[["plots"]][["min_sig_genes"]]
p <- p + 
  xlab('Total number of factors') +
  scale_x_continuous(breaks = seq(0, 14, by = 1)) +
  theme_bw()
p

### 
# pdf("/home/jmitchel/figures/scITD_revision_figs2/sle_only_min_sig_genes.pdf", width = 4, height = 4)
p
# dev.off()


##### running stability analysis
set.seed(1234)
pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(7,20),n_iterations=500,subset_type='subset', sub_prop=.85)

### Figure S2d left
# pdf(file = "/home/jmitchel/figures/for_paper_v2/stability_dsc2.pdf", useDingbats = FALSE,
#     width = 4, height = 3.5)
pbmc_container$plots$stability_plot_dsc
# dev.off()

### Figure S2d right
# pdf(file = "/home/jmitchel/figures/for_paper_v2/stability_lds2.pdf", useDingbats = FALSE,
#     width = 4, height = 3.5)
pbmc_container$plots$stability_plot_lds
# dev.off()


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

# saveRDS(pbmc_container,file='/home/jmitchel/data/lupus_data/lupus_container_w_decomp.rds')
pbmc_container <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_container_w_decomp.rds')


#### computing and adding relevant clinical associations to the metadata associations
# load data of categorical variables
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_categorical.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL

# get tucker donor scores to test
dsc <- pbmc_container$tucker_results[[1]]

## get donors in both dsc and in clin_vars
# trim donor IDs in dsc
trim_names <- sapply(rownames(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
old_names <- rownames(dsc)
names(old_names) <- trim_names
rownames(dsc) <- trim_names

# get donors in both dataframes
d_both <- rownames(clin_vars)[rownames(clin_vars) %in% rownames(dsc)]

# limit both dataframes to just the intersection of donors and in the same order
dsc <- dsc[d_both,]
clin_vars <- clin_vars[d_both,]

#### apply other checks to clin vars before testing!!
ndx_rem <- c()
for (j in 1:ncol(clin_vars)) {
  d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]
  tmp <- clin_vars[d_keep,j]
  num_levs <- length(unique(tmp))
  if (num_levs==1) {
    ndx_rem <- c(ndx_rem,j)
  }
}
clin_vars <- clin_vars[,-ndx_rem]
####

all_pvals <- c()
f_tested <- c()
c_tested <- c()
# loop through the variables to test
for (j in 1:ncol(clin_vars)) {
  print(j)
  # loop through factors
  for (f in 1:ncol(dsc)) {
    d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]
    
    tmp <- as.data.frame(cbind(dsc[d_keep,f], clin_vars[d_keep,j]))
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
print(all_pvals[order(all_pvals,decreasing=FALSE)][1:40])
print(f_tested[order(all_pvals,decreasing=FALSE)][1:40])
print(c_tested[order(all_pvals,decreasing=FALSE)][1:40])


##### now plotting some clinical associations for the categorical variables

### plot anti-dsDNA antibody hits for f1 sle only 
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"acrantidsdna"])]
tmp <- as.data.frame(cbind(dsc[d_keep,1], clin_vars[d_keep,"acrantidsdna"]))
colnames(tmp) <- c('dscore','cvar')

# force cvar to be factor
tmp$cvar <- as.factor(tmp$cvar)

tmp$cvar_word <- sapply(tmp$cvar,function(x){
  if (x==1) {
    return('anti-dsDNA')
  } else {
    return('No anti-dsDNA')
  }
})
tmp$cvar_word <- as.factor(tmp$cvar_word)

# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_antidsdna3.pdf", useDingbats = FALSE,
#     width = 3.75, height = 2.75)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.25, binwidth = .01) +
  ylab('Factor 1 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
# dev.off()



### plot anti-smith antibody hits for f1 sle only
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"acrantismith"])]
tmp <- as.data.frame(cbind(dsc[d_keep,1], clin_vars[d_keep,"acrantismith"]))
colnames(tmp) <- c('dscore','cvar')

# force cvar to be factor
tmp$cvar <- as.factor(tmp$cvar)

tmp$cvar_word <- sapply(tmp$cvar,function(x){
  if (x==1) {
    return('anti-smith')
  } else {
    return('No anti-smith')
  }
})
tmp$cvar_word <- as.factor(tmp$cvar_word)

# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_antismith3.pdf", useDingbats = FALSE,
#     width = 3.75, height = 2.75)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.25, binwidth = .01) +
  ylab('Factor 1 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
dev.off()




#### make empty matrix to add these associations to dscores metadata heatmap
res_add <- matrix(NA,nrow=5,ncol=ncol(dsc))

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





### calculating nephritis association using an interaction test with F1
# including prednisone variable as covariate
clin_vars2 <- read_excel('/home/jmitchel/data/lupus_data/SLE_meds_cleaned.xlsx')
clin_vars2 <- as.data.frame(clin_vars2)
rownames(clin_vars2) <- clin_vars2[,'Sample ID']
clin_vars2[,'Sample ID'] <- NULL
# make all NA into zeros, since no 0 are put in the table
clin_vars2[is.na(clin_vars2)] <- 0

tmp <- as.data.frame(cbind(dsc[,1],dsc[,2],clin_vars2[rownames(dsc),'pred_dose'],clin_vars[,'crflupusneph'],clin_vars[,'acrantidsdna'],dmeta[rownames(dsc),c('sex','Ethnicity','Age')]))
tmp <- tmp[order(tmp[,1],decreasing=TRUE),]
colnames(tmp)[1:5] <-  c('dscore1','dscore2','pred','ln','dsdna')

tmp$ln <- as.factor(tmp$ln)
tmp$dsdna <- as.factor(tmp$dsdna)

# logistic regression model
nmod <- glm(ln~dscore1+dscore2+pred+Ethnicity, data=tmp, family = 'binomial') ##"null" mod
fmod <- glm(ln~dscore1+dscore2+pred+Ethnicity+dscore2*dscore1, data=tmp, family = "binomial") ##"full" mod

a_res <- anova(nmod, fmod, test = 'Chisq')
pval <- a_res$`Pr(>Chi)`[2]
pval
res_add[3,2] <- pval


# plotting the interaction effect
tmp$ln <- factor(tmp$ln,levels=c(1,0))
levels(tmp$ln) <- c('Nephritis','No nephritis')

f1_thresh <- 0
tmp1 <- tmp[tmp$dscore1>f1_thresh,]
tmp2 <- tmp[tmp$dscore1<f1_thresh,]
tmp1 <- tmp1[tmp1$pred>0,]
tmp2 <- tmp2[tmp2$pred>0,]

p1 <- ggplot(tmp1,aes(x=ln,y=dscore2)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.25, binwidth = .01) +
  xlab('') +
  ylab('aSLE_F2 Scores') +
  coord_flip() +
  theme_bw() +
  ggtitle('Donors with high aSLE_F1 scores') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(-.28,.28)
p1

p2 <- ggplot(tmp2,aes(x=ln,y=dscore2)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.25, binwidth = .01) +
  xlab('') +
  ylab('aSLE_F2 Scores') +
  coord_flip() +
  theme_bw() +
  ggtitle('Donors with low aSLE_F1 scores') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(-.28,.28)
p2

tmp1 <- tmp[tmp$dscore1>f1_thresh,]
tmp2 <- tmp[tmp$dscore1<f1_thresh,]
tmp1 <- tmp1[tmp1$pred==0,]
tmp2 <- tmp2[tmp2$pred==0,]

p3 <- ggplot(tmp1,aes(x=ln,y=dscore2)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.25, binwidth = .01) +
  xlab('') +
  ylab('aSLE_F2 Scores') +
  coord_flip() +
  theme_bw() +
  ggtitle('Donors with high aSLE_F1 scores') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(-.28,.28)
p3

p4 <- ggplot(tmp2,aes(x=ln,y=dscore2)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.25, binwidth = .01) +
  xlab('') +
  ylab('aSLE_F2 Scores') +
  coord_flip() +
  theme_bw() +
  ggtitle('Donors with low aSLE_F1 scores') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(-.28,.28)
p4

fig <- cowplot::plot_grid(p1,p2,p3,p4,nrow=2,ncol=2,align = 'hv')
fig

### Figure S4a
pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/sle_f2_nephritis_interaction.pdf", useDingbats = FALSE,
    width = 8, height = 4)
fig
dev.off()



# co-occurrance of LN and dsdna with factor 2 with window test
old_dsc <- dsc
tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantidsdna']))
tmp <- tmp[order(tmp[,1],decreasing=TRUE),]
colnames(tmp) <-  c('dscore','ln','dsdna')

# remove rows that don't have dsdna==1
tmp <- tmp[tmp$dsdna==1,]
window_size <- 19
stored_counts <- c()
for (i in 1:(nrow(tmp)-window_size+1)) {
  tests <- tmp[i:(i+window_size-1),'ln']
  stored_counts <- c(stored_counts,sum(tests==1))
}
dscores <- tmp$dscore[(floor(window_size/2)+1):(nrow(tmp)-floor(window_size/2))]

plot_df <- cbind.data.frame(stored_counts,dscores)
plot_df$stored_counts <- (plot_df$stored_counts / window_size) * 100

### 
# pdf(file = "/home/jmitchel/figures/for_paper_v2/LN_antidsDNA_link2.pdf", useDingbats = FALSE,
#     width = 4, height = 3.25)
ggplot(plot_df,aes(x=dscores,y=stored_counts)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  xlab('Factor 2 Donor Score (window center)') +
  ylab('Percent in window with LN') +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))
# dev.off()


# calculating p-value for the association
myfstat <- cor(stored_counts,dscores,method='spearman')

all_scores <- c()
for (testndx in 1:10000) {
  
  # shuffle donor scores
  for (j in 1:ncol(dsc)) {
    dsc[,j] <- sample(dsc[,j])
  }
  
  # testing co-occurrance of LN and dsdna with factor 2
  tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantidsdna']))
  tmp <- tmp[order(tmp[,1],decreasing=TRUE),]
  colnames(tmp) <-  c('dscore','ln','dsdna')
  
  tmp <- tmp[tmp$dsdna==1,]
  
  stored_counts <- c()
  for (i in 1:(nrow(tmp)-window_size+1)) {
    tests <- tmp[i:(i+window_size-1),'ln']
    stored_counts <- c(stored_counts,sum(tests==1))
  }
  dscores <- tmp$dscore[(floor(window_size/2)+1):(nrow(tmp)-floor(window_size/2))]
  lmres <- lm(stored_counts~dscores)
  lmres <- summary(lmres)
  fstat <- cor(stored_counts,dscores,method='spearman')
  all_scores <- c(all_scores,fstat)
}

pval <- sum(all_scores>myfstat)/10000
print(pval)

dsc <- old_dsc


# res_add[3,2] <- pval











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
  
  tmp <- as.data.frame(cbind(dsc[d_keep,], clin_vars[d_keep,j]))
  colnames(tmp)[1:ncol(dsc)] <- sapply(1:ncol(dsc),function(x){paste0('Factor',x)})
  colnames(tmp)[ncol(dsc)+1] <- 'cvar'
  
  # using ordinal logistic regression
  tmp$cvar <- as.factor(tmp$cvar)
  m <- polr(cvar ~ ., data = tmp, Hess=TRUE, method='probit')
  
  ## view a summary of the model
  ctable <- coef(summary(m))
  ## calculate and store p values
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  pval <- p[1:ncol(dsc)]
  
  all_pvals <- c(all_pvals,pval)
  f_tested <- c(f_tested,1:ncol(dsc))
  c_tested <- c(c_tested,rep(colnames(clin_vars)[j],ncol(dsc)))
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
res_add[4,] <- c_pv[ndx_reorder]




##### testing factors against meds
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
dsc <- pbmc_container$tucker_results[[1]]

## get donors in both dsc and in clin_vars
# trim donor IDs in dsc
trim_names <- sapply(rownames(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
old_names <- rownames(dsc)
names(old_names) <- trim_names
rownames(dsc) <- trim_names

# get donors in both dataframes
d_both <- rownames(clin_vars)[rownames(clin_vars) %in% rownames(dsc)]

# limit both dataframes to just the intersection of donors and in the same order
dsc <- dsc[d_both,]
clin_vars <- clin_vars[d_both,]

all_pvals <- c()
f_tested <- c()
c_tested <- c()
# loop through the variables to test
for (j in 1:ncol(clin_vars)) {
  print(j)
  # loop through factors
  for (f in 1:ncol(dsc)) {
    # get donors in clin var that don't have an NA value
    d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]
    
    tmp <- as.data.frame(cbind(dsc[d_keep,f], clin_vars[d_keep,j]))
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




## now plotting the prednisone associations
tmp <- cbind.data.frame(dsc[,3],clin_vars[,'prednisone'])
colnames(tmp) <- c('dscore','cvar')
tmp$cvar_word <- sapply(tmp$cvar,function(x){
  if (x==1) {
    return('On Prednisone')
  } else {
    return('Not on Prednisone')
  }
})
tmp$cvar_word <- as.factor(tmp$cvar_word)

### Figure 2c
# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_prednisone_binary2.pdf", useDingbats = FALSE,
#     width = 4.5, height = 3.5)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 3 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
# dev.off()



## remove outlier for lm
tmp <- cbind.data.frame(dsc[,3],pred_dose[d_both,1])
colnames(tmp) <- c('dscore','cvar')
row_max <- which(tmp$cvar==max(tmp$cvar))
tmp2 <- tmp[-row_max,]
lmres <- lm(cvar~dscore,data=tmp2)
summary(lmres)

line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

# Figure S5a
# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_prednisone_dose2.pdf", useDingbats = FALSE,
#     width = 6, height = 4)
ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point(alpha = 0.3,pch=19,size=3) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Factor 3 Donor Score') +
  ylab('Prednisone Dose (mg)') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))
# dev.off()





# get a row of padj for prednisone associations
ndx_extract <- which(c_tested=='prednisone')
f_ordering <- f_tested[ndx_extract]
c_pv <- all_pvals[ndx_extract]
ndx_reorder <- order(f_ordering)
print(f_ordering[ndx_reorder])
print(c_pv[ndx_reorder])
res_add[5,] <- c_pv[ndx_reorder]

rownames(res_add) <- c('anti_dsDNA','anti_Smith','Nephritis','SLEDAI','Prednisone')

pbmc_container$meta_associations <- rbind(pbmc_container$meta_associations,res_add)

# reorder rows
pbmc_container$meta_associations <- pbmc_container$meta_associations[c('sex','Age',
                                                                       'processing','Ethnicity',
                                                                       'anti_dsDNA','anti_Smith',
                                                                       'SLEDAI','Prednisone','Nephritis'),]

## plot donor score
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval',
                                    meta_vars=c('sex'))

### Figure 2b
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/sle_only_dscores2.pdf", useDingbats = FALSE,
#     width = 4, height = 4)
pbmc_container$plots$donor_matrix
# dev.off()












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

### for hmap below I temporarily changed these parameters in run_gsea scITD code to alter the size
#   sim_hmap_res <- ht_clusters(mat, cl, word_cloud_grob_param = list(max_width = 70),
#   exclude_words=exclude_words, fontsize_range = c(10, 15))

###
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/F1_gsea_summary.pdf", useDingbats = FALSE,
#     width = 17, height = 10)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=1,direc='up',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))
# dev.off()

## looking at some gene sets that were significant
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=3)

## f1 sets to show on loading hmap
# gsets <- c("GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
#            "GOBP_RESPONSE_TO_INTERFERON_GAMMA",
#            "GOBP_PROTEOLYSIS",
#            "GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY",
#            "GOBP_INTERLEUKIN_1_PRODUCTION",
#            "GOBP_TUMOR_NECROSIS_FACTOR_SUPERFAMILY_CYTOKINE_PRODUCTION",
#            "GOBP_MYELOID_LEUKOCYTE_ACTIVATION",
#            "GOBP_INTERLEUKIN_6_PRODUCTION",
#            "GOBP_APOPTOTIC_CELL_CLEARANCE",
#            "GOBP_REGULATION_OF_LEUKOCYTE_PROLIFERATION",
#            "GOBP_APOPTOTIC_SIGNALING_PATHWAY",
#            "GOBP_REGULATION_OF_CELL_CYCLE",
#            "GOBP_NEGATIVE_REGULATION_OF_GROWTH",
#            "GOBP_VESICLE_BUDDING_FROM_MEMBRANE",
#            "GOBP_REGULATION_OF_NIK_NF_KAPPAB_SIGNALING")

gsets <- c("GOBP_RESPONSE_TO_TYPE_I_INTERFERON",
           "GOBP_MYELOID_LEUKOCYTE_ACTIVATION",
           "GOBP_REGULATION_OF_NIK_NF_KAPPAB_SIGNALING")

gset_cmap <- c('blue',
               'maroon',
               'forest green')
#######


names(gset_cmap) <- gsets

gset_cmap_sub <- gset_cmap[gset_cmap!='black']
gset_sub <- names(gset_cmap_sub)

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/sle_only_f1_lds.pdf", useDingbats = FALSE,
#     width = 12, height = 7)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=8, callout_ctypes=NULL, 
                                      le_set_callouts=gset_sub, le_set_colormap=gset_cmap_sub, le_set_num_per=7, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                      clust_method='complete', h_w=c(9,6.5))
# dev.off()

hm_list <- plot_select_sets(pbmc_container, 1, gsets, color_sets=gset_cmap, 
                            cl_rows=F, myfontsize=6.5, h_w=c(2,6.5))

hm_list

pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/sle_f1_sel_gsets.pdf", useDingbats = FALSE,
    width = 10, height = 3)
hm_list
dev.off()

### Figure 2d
## plotting just specific genes h_w=c(9,6.5) original
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/sle_only_f1_lds_v2.pdf", useDingbats = FALSE,
#     width = 10, height = 18)
g_show <- c('IFI6','ISG15','MX1','XAF1','USP18',
            'JUP','BATF2','STOM','LILRB2','PTPN22','IL18R1','RTKN2','CALR','FOXP3')
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, specific_callouts=g_show,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                      clust_method='complete', h_w=c(18,13))
# dev.off()




pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))

gsets <- c('GOBP_CELL_CYCLE',
           'GOBP_CELL_MIGRATION',
           'GOBP_P38MAPK_CASCADE')

gset_cmap <- c('blue','maroon','forest green')

names(gset_cmap) <- gsets

gset_cmap_sub <- gset_cmap[gset_cmap!='black']
gset_sub <- names(gset_cmap_sub)

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=8, callout_ctypes=NULL, 
                                      le_set_callouts=gset_sub, le_set_colormap=gset_cmap_sub, le_set_num_per=7, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                      clust_method='complete', h_w=c(9,6.5))

hm_list <- plot_select_sets(pbmc_container, 2, gsets, color_sets=gset_cmap, 
                            cl_rows=F, myfontsize=6.5, h_w=c(2,6.5))

hm_list

pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/sle_f2_sel_gsets.pdf", useDingbats = FALSE,
    width = 10, height = 3)
hm_list
dev.off()




pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=3, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))

gsets <- c('GOBP_RESPONSE_TO_HORMONE',
           'GOBP_POSITIVE_REGULATION_OF_CELL_DIFFERENTIATION',
           'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II')

gset_cmap <- c('blue','maroon','forest green')

names(gset_cmap) <- gsets

gset_cmap_sub <- gset_cmap[gset_cmap!='black']
gset_sub <- names(gset_cmap_sub)

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=3, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=8, callout_ctypes=NULL, 
                                      le_set_callouts=gset_sub, le_set_colormap=gset_cmap_sub, le_set_num_per=7, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                      clust_method='complete', h_w=c(9,6.5))

hm_list <- plot_select_sets(pbmc_container, 3, gsets, color_sets=gset_cmap, 
                            cl_rows=F, myfontsize=6.5, h_w=c(2,6.5))

hm_list

pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/sle_f3_sel_gsets.pdf", useDingbats = FALSE,
    width = 10, height = 3)
hm_list
dev.off()







##### project F1 onto whole dataset including healthy donors
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

param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc","cDC",
                                               "cMono","ncMono"),
                                ncores = 30, rand_seed = 10)

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
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')

# do the projection
pbmc_container_full <- project_new_data(pbmc_container_full,pbmc_container)

pbmc_container_full <- get_donor_meta(pbmc_container_full, additional_meta = 'Status')

tmp <- cbind.data.frame(pbmc_container_full$donor_metadata[,2],
                        pbmc_container_full$projected_scores[rownames(pbmc_container_full$donor_metadata),1])
colnames(tmp) <- c('status','proj_sc')
tmp$status <- factor(tmp$status)
levels(tmp$status) <- c('Healthy','SLE')
p <- ggplot(tmp,aes(x=status,y=proj_sc)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  xlab('') +
  ylab('Projected Factor 1 Scores') +
  coord_flip() +
  theme_bw()

### Figure S3a
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/proj_f1_status.pdf", useDingbats = FALSE,
#     width = 4, height = 3)
p
# dev.off()


# compute the status association p-value
lmres <- lm(proj_sc~status,data=tmp)
summary(lmres)














###### ploting the Th proportions and Treg proportion associations
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
dsc <- get_one_factor(pbmc_container,1)[[1]]
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
  geom_point(alpha = 0.5,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Factor 1 Donor Score') +
  ylab('Proportion Th/PBMCs') +
  theme_bw() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))

# now adding the beta and significance
tmp2 <- as.data.frame(scale(tmp))
lmres <- summary(lm(Th~dsc,data=tmp2))
print(lmres)

### Figure S3b
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/F1_Th_props.pdf", useDingbats = FALSE,
#     width = 4, height = 3)
p
# dev.off()


# now plotting the Treg association
my_ctype <- 'Th'
subc <- pbmc@meta.data[pbmc@meta.data$cg_cov==my_ctype,'ct_cov']
subc_lev <- unique(subc)
print(subc_lev)
subc <- factor(subc,levels=subc_lev)
levels(subc) <- 1:length(subc_lev)
names(subc) <- rownames(pbmc@meta.data)[pbmc@meta.data$cg_cov==my_ctype]
# convert subtype names to numbers
pbmc_container$subclusters[[my_ctype]][['res:0.5']] <- subc
pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype=my_ctype,res=.5,n_col=2)
dotplot <- get_subclust_enr_dotplot(pbmc_container,ctype='Th',res=.5,subtype=1,factor_use=1)

dotplot <- dotplot + ylim(0,.3) + ylab('Proportion Treg/Th') + 
  theme(plot.title = element_blank())

### Figure S3c
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/F1_T4reg_props.pdf", useDingbats = FALSE,
#     width = 4, height = 3)
dotplot
# dev.off()






### generating a plot to show actual expression values for top/bot donors
plot_donor_sig_genes_override <- function(container, factor_select, top_n_per_ctype,
                                          genes_plot, ct_in_hmap, ctypes_use=NULL, show_donor_labels=FALSE,
                                          additional_meta=NULL, add_genes=NULL, top_bot_n=NULL) {
  
  # extract tensor information
  tensor_data <- container$tensor_data
  donor_nm <- tensor_data[[1]]
  gene_nm  <- tensor_data[[2]]
  ctype_nm  <- tensor_data[[3]]
  tnsr <- tensor_data[[4]]
  
  # get the loadings matrix
  ldngs <- container$tucker_results[[2]]
  
  # break down a factor from the loadings matrix
  genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
  ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})
  
  sr_col <- ldngs[factor_select,]
  
  tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)
  
  # extract the genes to show
  if (is.null(ctypes_use)) {
    ctypes <- container$experiment_params$ctypes_use
  } else {
    ctypes <- ctypes_use
  }
  sig_vecs <- get_significance_vectors(container,factor_select,ctypes)
  ct_in_hmap <- factor(ct_in_hmap)
  
  # unfold tensor along donor mode
  donor_unfold <- rTensor::k_unfold(rTensor::as.tensor(tnsr),1)@data
  
  gn_ctype_cnames <- c()
  for (ct in ctype_nm) {
    for (gn in gene_nm) {
      gn_ctype_cnames <- c(gn_ctype_cnames,paste0(gn,"_",ct))
    }
  }
  
  colnames(donor_unfold) <- gn_ctype_cnames
  rownames(donor_unfold) <- donor_nm
  
  # testing out scaling the data to unit variance
  donor_unfold <- scale(donor_unfold)
  
  # subset data to just genes to plot
  donor_unfold_sub <- donor_unfold[,genes_plot]
  donor_unfold_sub <- t(donor_unfold_sub)
  
  # reorder donors by their score for the factor
  donor_scores <- container$tucker_results[[1]]
  donor_scores <- donor_scores[,factor_select]
  donor_unfold_sub <- donor_unfold_sub[,order(donor_scores)]
  donor_scores <- donor_scores[order(donor_scores)]
  
  donor_scores <- unlist(donor_scores)
  col_fun2 = circlize::colorRamp2(c(min(donor_scores), 0, max(donor_scores)), c("purple", "white", "green"))
  if (!is.null(top_bot_n)) {
    donor_scores <- donor_scores[c(1:top_bot_n,(length(donor_scores)-top_bot_n+1):length(donor_scores))]
  }
  ha <- ComplexHeatmap::HeatmapAnnotation(score = donor_scores,col=list(score=col_fun2),
                                          show_annotation_name=FALSE)
  
  
  if (!is.null(additional_meta)) {
    meta <- container$scMinimal_full$metadata[,c('donors',additional_meta)]
    meta <- unique(meta)
    rownames(meta) <- meta$donors
    meta$donors <- NULL
    meta <- meta[colnames(donor_unfold_sub),,drop=FALSE]
    
    # make all columns of meta to be factors
    for (i in 1:ncol(meta)) {
      meta[,i] <- factor(unlist(meta[,i]),levels=unique(unlist(meta[,i]))[order(unique(unlist(meta[,i])))])
    }
    
    set.seed(30)
    if (length(levels(meta)) < 3) {
      mycol <- RColorBrewer::brewer.pal(n = 3, name = "Paired")
    } else {
      mycol <- RColorBrewer::brewer.pal(n = length(levels(meta)), name = "Paired")
    }
    names(mycol) <- levels(meta)
    ta <- ComplexHeatmap::HeatmapAnnotation(df = meta, show_annotation_name=TRUE,
                                            col = list(df = mycol))
  } else {
    ta <- NULL
  }
  
  # rename genes
  rownames(donor_unfold_sub) <- sapply(rownames(donor_unfold_sub),function(x) {
    gn <- strsplit(x,split="_")[[1]][1]
    ct <- strsplit(x,split="_")[[1]][2]
    gn <- convert_gn(container,gn)
    return(paste0(gn,"_",ct))
  })
  
  rn_show <- sapply(rownames(donor_unfold_sub),function(x){
    strsplit(x,split="_")[[1]][[1]]
  })
  ct_show <- sapply(rownames(donor_unfold_sub),function(x){
    strsplit(x,split="_")[[1]][[2]]
  })
  ct_show <- factor(ct_show,levels=ctypes)
  
  set.seed(10)
  mycol <- RColorBrewer::brewer.pal(n = length(ctypes), name = "Accent")
  names(mycol) <- ctypes
  
  ct_annot <- ComplexHeatmap::rowAnnotation(cell_types=anno_simple(ct_show),
                                            show_annotation_name=FALSE,
                                            col = list(cell_types = mycol))
  
  # create the hmap
  col_fun = colorRamp2(c(min(donor_unfold_sub), 0, max(donor_unfold_sub)), c("blue", "white", "red"))
  names(rownames(donor_unfold_sub)) <- NULL
  if (!is.null(top_bot_n)) {
    donor_unfold_sub <- donor_unfold_sub[,c(1:top_bot_n,(ncol(donor_unfold_sub)-top_bot_n+1):ncol(donor_unfold_sub))]
  }
  myhmap <- Heatmap(donor_unfold_sub, name = "expr",
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    cluster_row_slices=FALSE,
                    column_names_gp = gpar(fontsize = 8),
                    row_names_gp = gpar(fontsize = 10),
                    col = col_fun, bottom_annotation=ha, row_split = ct_show,
                    row_labels=rn_show,border=TRUE, show_column_names=show_donor_labels,
                    left_annotation=ct_annot, show_row_dend = FALSE,
                    column_title = paste0('Factor ',as.character(factor_select)),
                    column_title_gp = gpar(fontsize = 20),
                    column_title_side = "top",
                    top_annotation=ta)
  
  container$plots$donor_sig_genes[[as.character(factor_select)]] <- myhmap
  return(container)
}

genes_plot <- c('IFI6_Th','ISG15_Th','MX1_Th','XAF1_Th','PTPN22_Th','IL18R1_Th','RTKN2_Th','FOXP3_Th','CALR_Th','IFI6_Tc','ISG15_Tc','MX1_Tc','XAF1_Tc','IFI6_cMono','ISG15_cMono','MX1_cMono','XAF1_cMono','JUP_cMono','BATF2_cMono','STOM_cMono','LILRB2_cMono')
names(genes_plot) <- sapply(genes_plot,function(x){
  strsplit(x,split='_')[[1]][[1]]
})
ct_in_hmap <- sapply(genes_plot,function(x){
  strsplit(x,split='_')[[1]][[2]]
})
names(ct_in_hmap) <- NULL
plot_donor_sig_genes_override(pbmc_container, factor_select=1, top_n_per_ctype=5,genes_plot=genes_plot,
                              ct_in_hmap=ct_in_hmap,    
                     ctypes_use=c('Th','Tc','NK','B','cMono','ncMono','cDC'), show_donor_labels=FALSE,
                     additional_meta=NULL, add_genes=NULL,top_bot_n=15)

p <- pbmc_container$plots$donor_sig_genes[[1]]

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/f1_top_bot_expression.pdf", useDingbats = FALSE,
#     width = 7, height = 5)
p
# dev.off()






















##### Plotting the loadings hmap for F2
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

t_all_spec_gns <- names(get_ct_specific_genes(pbmc_container,c('Th','Tc'),2,thresh = .05,signed = 'pos'))
cm_spec_gns <- names(get_ct_specific_genes(pbmc_container,c('cMono'),2,thresh = .01,signed = 'neg'))
cm_all_spec_gns <- names(get_ct_specific_genes(pbmc_container,c('cMono','ncMono','cDC'),2,thresh = .01,signed = 'neg'))
th_spec_gns <- names(get_ct_specific_genes(pbmc_container,c('Th'),2,thresh = .01,signed = 'neg'))
b_spec_gns <- names(get_ct_specific_genes(pbmc_container,c('B'),2,thresh = .01,signed = 'pos'))

spec_callouts <- c('CXCR5','NFKBIB','DUSP10','G0S2','MXD1','LAIR1','CD3D','TRIM21')

#### run gsea for a f2
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/F1_gsea_summary.pdf", useDingbats = FALSE,
#     width = 17, height = 10)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=2,direc='up',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))
# dev.off()
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=8)


### Figure 2f
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/F2_lds.pdf", useDingbats = FALSE,
#     width = 5, height = 5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts = TRUE, specific_callouts = spec_callouts,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                      clust_method='complete', h_w=c(5,6.5))
# dev.off()




##### computing activated B cell proportion association
## Computed my own subclusters to see if we could identify other canonical ones not in original publication
# add conos object generated in preprocessing/embedding_prep.R file
con <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_conos2.rds')
pbmc_container$embedding <- con

# large number of cores seems to hamper some stuff below
pbmc_container$embedding$n.cores <- 5

pbmc_container <- get_subtype_prop_associations(pbmc_container, max_res=.9, stat_type='adj_pval',
                                                min_cells_group=200)


# saveRDS(pbmc_container$subclusters,file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')
pbmc_container$subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')

# plot umap for the B cell subclusters
all_ctypes=c('B')
all_res=c(.8)

# get all subcluster umaps
pbmc_container <- get_subclust_umap(pbmc_container,all_ctypes=all_ctypes,
                                    all_res=all_res,n_col=1)

### 
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/B_subc_umap.pdf", useDingbats = FALSE,
#     width = 4, height = 3)
pbmc_container[["plots"]][["subc_umaps"]][["B:0.8"]] +   
  scale_y_reverse() +
  scale_x_reverse() 
# dev.off()

# get de genes
de_hmaps <- get_subclust_de_hmaps(pbmc_container,all_ctypes,all_res)
cowplot::plot_grid(plotlist=de_hmaps,nrow=1)

my_ctype <- 'B'
pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype=my_ctype,res=.8,n_col=2) # selecting resolution that contains the activated memory subtype
pbmc_container$plots$ctype_prop_factor_associations

dotplot <- get_subclust_enr_dotplot(pbmc_container,ctype='B',res=.8,subtype=4,factor_use=2) # get the association for the correct subtype

dotplot <- dotplot + ylim(0,.75) + ylab('Proportion Bact/B') + 
  theme(plot.title = element_blank())

### Figure S4f
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/F2_bact_props.pdf", useDingbats = FALSE,
#     width = 4, height = 3)
dotplot
# dev.off()

















## Plotting the loadings hmap for F3
t_all_spec_gns <- names(get_ct_specific_genes(pbmc_container,c('Th'),3,thresh = .05,signed = 'pos'))
t_all_spec_gns <- names(get_ct_specific_genes(pbmc_container,c('Th'),3,thresh = .05,signed = 'neg'))
B_all_spec_gns <- names(get_ct_specific_genes(pbmc_container,c('B'),3,thresh = .05,signed = 'neg'))
cM_all_spec_gns <- names(get_ct_specific_genes(pbmc_container,c('cMono'),3,thresh = .05,signed = 'neg'))


spec_callouts_ct <- c('ICA1','IL2RA','CD28','TCL1A','TLR7','IL1B')

spec_callouts <- c(spec_callouts_cort,spec_callouts_ct)

### Figure 2h
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/F3_lds.pdf", useDingbats = FALSE,
#     width = 5, height = 5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=3, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts = TRUE, specific_callouts = spec_callouts,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE,
                                      clust_method='complete', h_w=c(5,6.5))
# dev.off()





## showing Treg association with this factor as well
my_ctype <- 'Th'
subc <- pbmc@meta.data[pbmc@meta.data$cg_cov==my_ctype,'ct_cov']
subc_lev <- unique(subc)
print(subc_lev)
subc <- factor(subc,levels=subc_lev)
levels(subc) <- 1:length(subc_lev)
names(subc) <- rownames(pbmc@meta.data)[pbmc@meta.data$cg_cov==my_ctype]
# convert subtype names to numbers
pbmc_container$subclusters[[my_ctype]][['res:0.5']] <- subc
pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype=my_ctype,res=.5,n_col=2)
pbmc_container$plots$ctype_prop_factor_associations

dotplot <- get_subclust_enr_dotplot(pbmc_container,ctype='Th',res=.5,subtype=1,factor_use=3)

dotplot <- dotplot + ylim(0,.5) + ylab('Proportion Treg/Th') + 
  theme(plot.title = element_blank())

### Figure S5e
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/F3_Treg_props.pdf", useDingbats = FALSE,
#     width = 3.5, height = 3)
dotplot
# dev.off()












