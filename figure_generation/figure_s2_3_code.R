library(scITD)
library(Seurat)
library(ggplot2)
library(readxl)
library(MASS)

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


pbmc_container_full <- run_tucker_ica(pbmc_container_full, ranks=c(7,20),
                                 tucker_type = 'regular', rotation_type = 'hybrid')

# flip sign of F1 so high ISG expression is positive instead of negative (signs are arbitrary)
pbmc_container_full$tucker_results[[1]][,1] <- pbmc_container_full$tucker_results[[1]][,1] * -1
pbmc_container_full$tucker_results[[2]][1,] <- pbmc_container_full$tucker_results[[2]][1,] * -1




###### starting sle-only analysis
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

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Ethnicity'),
                                        stat_use='pval')

## plot donor score
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_only_dscores3.pdf", useDingbats = FALSE,
#     width = 5, height = 6)
pbmc_container$plots$donor_matrix
# dev.off()


# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)


## trying to draw loadings plots individually because I want them to have different params...
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(7,3.5))

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(7,3.5))

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=3, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(7,3.5))

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=4, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(7,3.5))

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=5, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(7,3.5))

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=6, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=9, h_w=c(7,3.5))

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=7, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(7,3.5))


# render them together because it's easier for the final figure
myfig <- render_multi_plots(pbmc_container,data_type='loadings',max_cols=4)

# pdf(file = "/home/jmitchel/figures/for_paper_v2/lds_all3.pdf", useDingbats = FALSE,
#     width = 15, height = 18)
myfig
dev.off()






##### comparing sle-only decomposition to full decomposition
tr_sle <- pbmc_container$tucker_results

tr_full <- pbmc_container_full$tucker_results # after re-running full decomp

pdf(file = "/home/jmitchel/figures/for_paper_v2/decomp_comparison3.pdf", useDingbats = FALSE,
    width = 8, height = 4.5)
compare_decompositions(tr_sle,tr_full,c('SLE only','Full dataset',use_text=TRUE))
dev.off()



##### running stability analysis
set.seed(1234)
pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(7,20),n_iterations=500,subset_type='subset', sub_prop=.85)

pdf(file = "/home/jmitchel/figures/for_paper_v2/stability_dsc2.pdf", useDingbats = FALSE,
    width = 4, height = 3.5)
pbmc_container$plots$stability_plot_dsc
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/stability_lds2.pdf", useDingbats = FALSE,
    width = 4, height = 3.5)
pbmc_container$plots$stability_plot_lds
dev.off()















##### now computing clinical associations
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
dev.off()



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


# co-occurrance of LN and dsdna with factor 2
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
pdf(file = "/home/jmitchel/figures/for_paper_v2/LN_antidsDNA_link2.pdf", useDingbats = FALSE,
    width = 4, height = 3.25)
ggplot(plot_df,aes(x=dscores,y=stored_counts)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  xlab('Factor 2 Donor Score (window center)') +
  ylab('Percent in window with LN') +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))
dev.off()


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




# plotting F1 against sledai score
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"sledaiscore"])]
tmp <- as.data.frame(cbind(dsc[d_keep,1], clin_vars[d_keep,"sledaiscore"]))
colnames(tmp) <- c('dscore','cvar')

# force cvar to be numeric so can calculate the lm line
tmp$cvar <- as.numeric(tmp$cvar)

# adding reg line
lmres <- lm(cvar~dscore,data=tmp)
line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_sledai3.pdf", useDingbats = FALSE,
#     width = 3.75, height = 3)
ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point(alpha = 0.25,pch=19,size=3) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Factor 1 Donor Score') +
  ylab('SLEDAI Score') +
  theme_classic()
dev.off()


##### getting av ifn associations for comparison to F1
rownames(clin_vars) <- old_names[rownames(clin_vars)]
for (ct in pbmc_container$experiment_params$ctypes_use) {
  pb <- pbmc_container[["scMinimal_ctype"]][[ct]][["pseudobulk"]]
  go_ifn <- read.csv(file='/home/jmitchel/IFN_gene_list.csv') # GO TI IFN set
  go_ifn <- go_ifn[go_ifn[,1] %in% colnames(pb),1]
  
  pb_ifn <- pb[,go_ifn] # subset expression to just the ifn genes
  trim_names <- sapply(rownames(pb_ifn), function(x) {
    strsplit(x,split='_')[[1]][[1]]
  })
  
  inf_egene <- rowMeans(pb_ifn) # average IFN gene expression
  names(inf_egene) <- trim_names
  
  dsc <- pbmc_container$tucker_results[[1]]
  dsc <- dsc[names(trim_names),]
  rownames(dsc) <- trim_names
  
  # evaluating association between IFN egene and sledai
  tmp <- cbind.data.frame(clin_vars[trim_names,'sledaiscore'],dsc[trim_names,1],inf_egene[trim_names])
  colnames(tmp) <- c('sledai','dsc','expres')
  
  print(ct)
  print(cor(tmp$sledai,tmp$expres))
}
print(cor(tmp$sledai,tmp$dsc))










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

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_prednisone_binary2.pdf", useDingbats = FALSE,
    width = 4.5, height = 3.5)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 3 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
dev.off()



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

# make plot with best fit line
pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_prednisone_dose2.pdf", useDingbats = FALSE,
    width = 6, height = 4)
ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point(alpha = 0.3,pch=19,size=3) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Factor 3 Donor Score') +
  ylab('Prednisone Dose (mg)') +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))
dev.off()

















##### getting enriched gene sets for factor 2
# run gsea for a f2
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))

plot_gsea_hmap_w_similarity(pbmc_container,factor_select=2,direc='up',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))

plot_gsea_sub(pbmc_container,thresh=.05,clust_select=1)

## f2 sets to show on loading hmap
gsets <- c('GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_STRESS',
           'GOBP_CELL_CYCLE',
           'GOBP_RESPONSE_TO_GROWTH_FACTOR',
           'GOBP_REGULATION_OF_CELL_DEATH',
           'GOBP_REGULATION_OF_CELL_DIFFERENTIATION',
           'GOBP_P38MAPK_CASCADE',
           'GOBP_CELL_MIGRATION',
           'GOBP_CELL_POPULATION_PROLIFERATION',
           'GOBP_POSITIVE_REGULATION_OF_CELL_ADHESION',
           'GOBP_POSITIVE_REGULATION_OF_CELL_GROWTH',
           'GOBP_LYMPHOCYTE_COSTIMULATION',
           'GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY')

gset_cmap <- rep('black',length(gsets))

names(gset_cmap) <- gsets

hm_list <- plot_select_sets(pbmc_container, 2, gsets, color_sets=gset_cmap, 
                            cl_rows=T, myfontsize=5, h_w=c(6,4.5))

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_f2_gsets2.pdf", useDingbats = FALSE,
    width = 6, height = 5)
hm_list
dev.off()








## getting enriched gene sets for factor 3
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=3, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=3,direc='up',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=8)

## sets to show on loading hmap
gsets <- c('GOBP_RESPONSE_TO_HORMONE',
           'GOBP_RESPONSE_TO_CORTICOSTEROID',
           'GOBP_CELL_MIGRATION',
           'GOBP_REGULATION_OF_T_CELL_ACTIVATION',
           'GOBP_POSITIVE_REGULATION_OF_CELL_DIFFERENTIATION',
           'GOBP_APOPTOTIC_PROCESS',
           'GOBP_LEUKOCYTE_PROLIFERATION',
           'GOBP_REGULATION_OF_CELL_ADHESION',
           'GOBP_NEGATIVE_REGULATION_OF_SIGNALING',
           'GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_OR_POLYSACCHARIDE_ANTIGEN_VIA_MHC_CLASS_II',
           'GOBP_POSITIVE_REGULATION_OF_LIPID_BIOSYNTHETIC_PROCESS')

gset_cmap <- rep('black',length(gsets))

names(gset_cmap) <- gsets

hm_list <- plot_select_sets(pbmc_container, 3, gsets, color_sets=gset_cmap, 
                            cl_rows=T, myfontsize=5, h_w=c(6,4.5))

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_f3_gsets2.pdf", useDingbats = FALSE,
    width = 6, height = 5)
hm_list
dev.off()
























