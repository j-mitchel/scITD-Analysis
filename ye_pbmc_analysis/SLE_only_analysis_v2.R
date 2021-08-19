


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


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20,7),
                                 tucker_type = 'regular', rotation_type = 'hybrid')

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Ethnicity'),
                                        stat_use='pval')

## plot donor score
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_only_dscores.pdf", useDingbats = FALSE,
    width = 5, height = 6)
pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 5, height = 6)
pbmc_container$plots$donor_matrix
dev.off()



## storing the container so I can compare it to full all donors decomp
pbmc_container_SLE <- pbmc_container
tr_sle <- pbmc_container_SLE$tucker_results

tr_full <- pbmc_container$tucker_results # after re-running full decomp

pdf(file = "/home/jmitchel/figures/for_paper_v2/decomp_comparison.pdf", useDingbats = FALSE,
    width = 8, height = 4.5)
compare_decompositions(tr_sle,tr_full,c('SLE only','Full dataset',use_text=TRUE))
dev.off()

pbmc_container <- pbmc_container_SLE # resetting pbmc_container to sle only decomp



# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)


## trying to draw loadings plots individually because I want them to have different params...
pdf(file = "/home/jmitchel/figures/for_paper_v2/lds_f1.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(7,3.5))
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/lds_f2.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(7,3.5))
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/lds_f3.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=3, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(7,3.5))
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/lds_f4.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=4, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(7,3.5))
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/lds_f5.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=5, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(7,3.5))
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/lds_f6.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=6, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=9, h_w=c(7,3.5))
dev.off()


pdf(file = "/home/jmitchel/figures/for_paper_v2/lds_f7.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=7, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5, h_w=c(7,3.5))
dev.off()


# will render them together because it's easier for the final figure
myfig <- render_multi_plots(pbmc_container,data_type='loadings',max_cols=4)

pdf(file = "/home/jmitchel/figures/for_paper_v2/lds_all.pdf", useDingbats = FALSE,
    width = 15, height = 18)
myfig
dev.off()



#### now plotting some clinical associations
# need to have run appropriate lines in clinical_tests.R first

### plot antibody hits for f1 sle only acrantidsdna
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

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_antidsdna.pdf", useDingbats = FALSE,
    width = 3.75, height = 2.75)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 1 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
dev.off()



### plot antibody hits for f1 sle only acrantismith
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

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_antismith.pdf", useDingbats = FALSE,
    width = 3.75, height = 2.75)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 1 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
dev.off()


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

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_sledai.pdf", useDingbats = FALSE,
    width = 3.75, height = 3)
ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point(alpha = 0.25,pch=19,size=3) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Factor 1 Donor Score') +
  ylab('SLEDAI Score') +
  theme_classic()
dev.off()



# co-occurrance of LN and dsdna with factor 2
# old_dsc <- dsc
dsc <- old_dsc
tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantidsdna']))
# tmp <- as.data.frame(cbind(dsc[,4],clin_vars[,'crflupusneph'],clin_vars[,'acrantismith']))
tmp <- tmp[order(tmp[,1],decreasing=TRUE),]
colnames(tmp) <-  c('dscore','ln','dsdna')

# remove rows that don't have dsdna==1
tmp <- tmp[tmp$dsdna==1,]
# window_size <- 21 #size I originally used
window_size <- 17
stored_counts <- c()
for (i in 1:(nrow(tmp)-window_size+1)) {
  tests <- tmp[i:(i+window_size-1),'ln']
  stored_counts <- c(stored_counts,sum(tests==1))
}
dscores <- tmp$dscore[(floor(window_size/2)+1):(nrow(tmp)-floor(window_size/2))]
plot(dscores,stored_counts)

lmres <- lm(stored_counts~dscores)
lmres <- summary(lmres)
myfstat <- lmres$fstatistic[[1]]
myfstat <- cor(stored_counts,dscores,method='spearman')
myfstat <- cor(stored_counts,dscores,method='pearson')
myfstat

plot_df <- cbind.data.frame(stored_counts,dscores)
pdf(file = "/home/jmitchel/figures/for_paper_v2/LN_antidsDNA_link.pdf", useDingbats = FALSE,
    width = 4, height = 3.25)
ggplot(plot_df,aes(x=dscores,y=stored_counts)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  xlab('Factor 2 Donor Score (window center)') +
  ylab('Number of Joint Occurrences\nof LN + anti_dsDNA') +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))
dev.off()


# ## seeing how convincing the age association looks
# dscores <- pbmc_container[["tucker_results"]][[1]]
# meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Status','Age')]
# meta <- unique(meta)
# rownames(meta) <- meta$donors
# meta$donors <- NULL
# 
# tmp <- cbind.data.frame(dscores[,4],meta[rownames(dscores),2])
# colnames(tmp) <- c('dsc','Age')
# f4_age <- ggplot(tmp,aes(x=dsc,y=Age)) +
#   geom_point(pch=19,alpha=.5) +
#   ylab('Age') +
#   xlab('Factor 4 Donor Score') +
#   theme_bw()
# 
# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_f4_age.pdf", useDingbats = FALSE,
#     width = 4.5, height = 2)
# f4_age
# dev.off()



## getting enriched gene sets for factor 2
# run gsea for a f1
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 14, height = 7)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=2,direc='up',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))
dev.off()

pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 14, height = 14)
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=10)
dev.off()

## f4 sets to show on loading hmap
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
           'GOBP_TUMOR_NECROSIS_FACTOR_MEDIATED_SIGNALING_PATHWAY',
           'GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY')

gset_cmap <- rep('black',length(gsets))

names(gset_cmap) <- gsets

hm_list <- plot_select_sets(pbmc_container, 2, gsets, color_sets=gset_cmap, 
                            cl_rows=T, myfontsize=5, h_w=c(6,4.5))

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_f2_gsets.pdf", useDingbats = FALSE,
    width = 6, height = 5)
hm_list
dev.off()




#### exploratory analysis of IFN associations with sledai score

# try plotting IFI6 levels vs sledai score
d_exp <- pbmc_container[["scMinimal_ctype"]][['Th']][["pseudobulk"]][,'MX2']

trim_names <- sapply(names(d_exp), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(d_exp) <- trim_names

tmp <- cbind.data.frame(clin_vars[trim_names,'sledaiscore'],dsc[trim_names,1],d_exp[trim_names])
colnames(tmp) <- c('sledai','dsc','expres')
lm1 <- lm(sledai~expres,data=tmp)
lm2 <- lm(sledai~dsc+expres,data=tmp)
anova(lm1,lm2)

# getting the eigengene for the GO set of IFN genes
pb <- pbmc_container[["scMinimal_ctype"]][['Th']][["pseudobulk"]]
go_ifn <- read.csv(file='/home/jmitchel/IFN_gene_list.csv')
go_ifn <- go_ifn[go_ifn[,1] %in% colnames(pb),1]

go_ifn <- c('HERC5', 'IFI27', 'IRF7', 'ISG15', 'LY6E', 'MX1', 'OAS2', 'OAS3', 'RSAD2', 'USP18', 'GBP5') # or use Rao set
go_ifn <- go_ifn[go_ifn %in% colnames(pb)]

pb_ifn <- pb[,go_ifn] # subset expression to just the ifn genes
inf_egene <- svd(pb_ifn)$u[,1] # get IFN eigengene
names(inf_egene) <- rownames(pb_ifn)
names(inf_egene) <- trim_names

# evaluating association between IFN egene and sledai
tmp <- cbind.data.frame(clin_vars[trim_names,'sledaiscore'],dsc[trim_names,1],inf_egene[names(trim_names)])
tmp <- cbind.data.frame(clin_vars[trim_names,'acrantidsdna'],dsc[trim_names,1],inf_egene[names(trim_names)])
# tmp <- cbind.data.frame(clin_vars[trim_names,'acrantidsdna'],dsc[trim_names,1],inf_egene[trim_names])
colnames(tmp) <- c('sledai','dsc','expres')
lm1 <- lm(sledai~expres,data=tmp)
lm2 <- lm(sledai~dsc+expres,data=tmp)
anova(lm1,lm2)

lm2 <- lm(sledai~dsc,data=tmp)

coxtest(lm1,lm2)
jtest(lm1,lm2)
encomptest(lm1,lm2)


# trying logistic regression method
tmp <- cbind.data.frame(clin_vars[trim_names,'acrantidsdna'],dsc[trim_names,1],inf_egene[trim_names])
tmp <- cbind.data.frame(clin_vars[trim_names,'acrantismith'],dsc[trim_names,1],inf_egene[trim_names])
colnames(tmp) <- c('dsdna','dsc','expres')
tmp$dsdna <- as.factor(tmp$dsdna)
m1 <- glm(dsdna~expres, data=tmp, family = "binomial") ##"full" mod
m2 <- glm(dsdna~dsc, data=tmp, family = "binomial") ##"full" mod

AIC(m1)
AIC(m2)


## trying it with the all cell types together pseudobulk
cdata <- t(pbmc_container$scMinimal_full$count_data)
donor_meta <- as.factor(pbmc_container$scMinimal_full$metadata$donors)
donor_sum_counts <- get_sums(cdata,donor_meta)
donor_sum_counts <- donor_sum_counts[2:nrow(donor_sum_counts),]
donor_sum_counts <- donor_sum_counts[rowSums(donor_sum_counts)>0,] # remove donors with no cells
donor_sum_counts <- donor_sum_counts[rownames(pbmc_container$tucker_results[[1]]),] # only keep donors tested
donor_sum_counts <- t(donor_sum_counts)
donor_sum_counts <- Matrix(donor_sum_counts, sparse = TRUE)
all_nf <- edgeR::calcNormFactors(donor_sum_counts)
lib_sizes <- Matrix::colSums(donor_sum_counts)
donor_sum_counts <- sweep(donor_sum_counts,MARGIN=2,lib_sizes*all_nf,FUN='/')
donor_sum_counts <- log1p(donor_sum_counts * 10000)
donor_sum_counts <- t(donor_sum_counts)
donor_sum_counts <- Matrix(donor_sum_counts, sparse = TRUE)

pb <- donor_sum_counts
trim_names <- sapply(names(d_exp), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(d_exp) <- trim_names








## getting enriched gene sets for factor 3
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=3, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 14, height = 7)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=3,direc='down',thresh=.05,
                            exclude_words=c('regulation','positive','negative'))
dev.off()

pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 14, height = 14)
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=8)
dev.off()

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

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_f3_gsets.pdf", useDingbats = FALSE,
    width = 6, height = 5)
hm_list
dev.off()




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

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_prednisone_binary.pdf", useDingbats = FALSE,
    width = 4.5, height = 3.5)
pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_prednisone_binary.pdf", useDingbats = FALSE,
    width = 4, height = 2.75)
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
pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_prednisone_dose.pdf", useDingbats = FALSE,
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











## running LR analysis with new functions
# prep for new LR analysis
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
  rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
  return(rname)
})
lr_pairs$interaction_name <- NULL


pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='pool')
sft_thresh <- c(12,14,12,10,12,9,12)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

lr_hmap <- compute_LR_interact_v4(pbmc_container, lr_pairs, sig_thresh=0.0000000001,
                                   percentile_exp_rec=0.75, add_ld_fact_sig=TRUE)

lr_hmap <- myhmap1
pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_new_lr.pdf", useDingbats = FALSE,
    width = 8.25, height = 8.25)
lr_hmap
dev.off()

# compute number of obtained interactions
sum(pbmc_container[["lr_res"]]<.05) # contains duplicates because may have multiple module per target ct

lr_res <- pbmc_container[["lr_res"]]
unique_channels <- c()
for (i in 1:nrow(lr_res)) {
  lig_ct_rec <- lr_res[i,]
  lig_ct_rec_name <- strsplit(rownames(lr_res)[i],split='_')[[1]]
  lig <- lig_ct_rec_name[[1]]
  source <- lig_ct_rec_name[[2]]
  n_rec_comps <- length(lig_ct_rec_name) - 2
  rec_elements <- lig_ct_rec_name[3:(2+n_rec_comps)]
  rec_nm <- rec_elements[1]
  if (n_rec_comps>1) {
    for (j in 2:length(rec_elements)) {
      rec_nm <- paste0(rec_nm,"_",rec_elements[j])
    }
  }

  for (j in 1:ncol(lr_res)) {
    pv <- lig_ct_rec[j]
    if (pv < 0.05) {
      target_ct <- strsplit(names(pv),split="_")[[1]][[1]]
      lig_source_rec_target <- paste0(lig,"_",source,"_",rec_nm,"_",target_ct)
      unique_channels <- c(unique_channels,lig_source_rec_target)
    }
  }
}
unique_channels <- unique(unique_channels)
length(unique_channels)



## compute enrichment of module with the nichenet scores
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
test = ligand_target_matrix[,'TNFSF13B']
test = ligand_target_matrix[,'ICOSLG']
b_mod = pbmc_container[["module_genes"]][["B"]]
b_mod = pbmc_container[["module_genes"]][["Th"]]
mymod <- c(5)
g_in_mod <- names(b_mod)[b_mod%in%mymod]
g_not_mod <- names(b_mod)[!(b_mod%in%mymod)]
tmp <- cbind.data.frame(c(g_in_mod,g_not_mod),
                        c(rep('B_m1',length(g_in_mod)),rep('other',length(g_not_mod))))
tmp <- cbind.data.frame(c(g_in_mod,g_not_mod),
                        c(rep('Th_m5',length(g_in_mod)),rep('other',length(g_not_mod))))
colnames(tmp) <- c('gn','in_mod')
tmp$in_mod <- factor(tmp$in_mod,levels=c('Th_m5','other'))
target_scores <- test[tmp$gn]
tmp$target_scores <- target_scores
tmp <- tmp[which(!is.na(target_scores)),]
ggplot(tmp,aes(x=as.factor(in_mod),y=target_scores)) +
  geom_boxplot(notch=TRUE) +
  ylab('NicheNet regulatory potential') +
  xlab('')
wilcox.test(target_scores~in_mod,data=tmp)


pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=1, 
                                      sig_thresh=0.001, percentile_exp_rec=.75,
                                      show_rec_sig=FALSE) # used for fig


pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_f1_lr.pdf", useDingbats = FALSE,
    width = 7.5, height = 9)
pbmc_container$plots$lr_analysis[['Factor1']]
dev.off()

# getting GO enrichment HMAPs for modules
ctypes <- c('B')
modules <- c(1)

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use='TF')
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.00001, db_use=c('GO'),max_plt_pval=.01)
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.01, db_use=c('Reactome'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.01, db_use=c('KEGG'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.01, db_use=c('BioCarta'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('Hallmark'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.01, db_use=c('BioCarta','KEGG','Reactome','TF','Hallmark','GO'))

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_TNFSF13B_gsets.pdf", useDingbats = FALSE,
    width = 7, height = 10)
mod_enr
dev.off()

## making simplify enrichmet hmap instead
mat = GO_similarity(rownames(myres2))
pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_TNFSF13B_gsets.pdf", useDingbats = FALSE,
    width = 6, height = 5)
df = simplifyGO(mat,word_cloud_grob_param = list(max_width = 80),fontsize_range=c(7,15),max_words=5)
dev.off()


lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,mod_ct='Th',mod=8,lig_ct='cMpno',lig='ICOSLG')
lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=5,mod_ct='cMono',mod=3,lig_ct='Th',lig='TNFSF8')
lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=1,mod_ct='B',mod=1,lig_ct='cMono',lig='TNFSF13B')
pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_TNFSF13B_trio.pdf", useDingbats = FALSE,
    width = 6, height = 5)
lig_mod_fact
dev.off()

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=1,mod_ct='Th',mod=2,lig_ct='ncMono',lig='CXCL10')
pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_CXCL10_trio.pdf", useDingbats = FALSE,
    width = 6, height = 5)
lig_mod_fact
dev.off()


pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=2, 
                                      sig_thresh=0.001, percentile_exp_rec=.75,
                                      show_rec_sig=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_f2_lr.pdf", useDingbats = FALSE,
    width = 6.25, height = 5.5)
pbmc_container$plots$lr_analysis[['Factor2']]
dev.off()


ctypes <- c('Th','Th')
modules <- c(6,5)
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.01, db_use=c('GO'),max_plt_pval=.01)
pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_ICOSLG_gsets.pdf", useDingbats = FALSE,
    width = 7.5, height = 9)
mod_enr
dev.off()

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,mod_ct='Th',mod=6,lig_ct='cMono',lig='ICOSLG')
pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_ICOSLG_trio.pdf", useDingbats = FALSE,
    width = 6, height = 5)
lig_mod_fact
dev.off()


## now comparing LR results to "standard" LR analysis
pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=1, 
                                      sig_thresh=0.05, percentile_exp_rec=.75,
                                      show_rec_sig=FALSE)
test <- pbmc_container$lr_res_raw[['Factor1']] 

# extract the comm channels
mythresh=.35
all_channels <- c()
for (i in 1:nrow(test)) {
  for (j in 1:ncol(test)) {
    if (abs(test[i,j]) > mythresh) {
      l_ct_r <- rownames(test)[i]
      target <- colnames(test)[j]
      target <- strsplit(target,split='_')[[1]][[1]]
      mychannel <- paste0(l_ct_r,"_",target)
      all_channels <- c(all_channels,mychannel)
    }
  }
}
all_channels <- unique(all_channels)
length(all_channels)

pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=2, 
                                      sig_thresh=0.05, percentile_exp_rec=.75,
                                      show_rec_sig=FALSE)
test <- pbmc_container$lr_res_raw[['Factor2']] 

# add channels found in factor 2
for (i in 1:nrow(test)) {
  for (j in 1:ncol(test)) {
    if (abs(test[i,j]) > mythresh) {
      l_ct_r <- rownames(test)[i]
      target <- colnames(test)[j]
      target <- strsplit(target,split='_')[[1]][[1]]
      mychannel <- paste0(l_ct_r,"_",target)
      all_channels <- c(all_channels,mychannel)
    }
  }
}
all_channels <- unique(all_channels)
length(all_channels)


pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=3, 
                                      sig_thresh=0.05, percentile_exp_rec=.75,
                                      show_rec_sig=FALSE)
test <- pbmc_container$lr_res_raw[['Factor3']] 

for (i in 1:nrow(test)) {
  for (j in 1:ncol(test)) {
    if (abs(test[i,j]) > mythresh) {
      l_ct_r <- rownames(test)[i]
      target <- colnames(test)[j]
      target <- strsplit(target,split='_')[[1]][[1]]
      mychannel <- paste0(l_ct_r,"_",target)
      all_channels <- c(all_channels,mychannel)
    }
  }
}
all_channels <- unique(all_channels)
length(all_channels)


pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=4, 
                                      sig_thresh=0.05, percentile_exp_rec=.75,
                                      show_rec_sig=FALSE)
test <- pbmc_container$lr_res_raw[['Factor4']]

for (i in 1:nrow(test)) {
  for (j in 1:ncol(test)) {
    if (abs(test[i,j]) > mythresh) {
      l_ct_r <- rownames(test)[i]
      target <- colnames(test)[j]
      target <- strsplit(target,split='_')[[1]][[1]]
      mychannel <- paste0(l_ct_r,"_",target)
      all_channels <- c(all_channels,mychannel)
    }
  }
}
all_channels <- unique(all_channels)
length(all_channels)


pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=5, 
                                      sig_thresh=0.05, percentile_exp_rec=.75,
                                      show_rec_sig=FALSE)
test <- pbmc_container$lr_res_raw[['Factor5']]

for (i in 1:nrow(test)) {
  for (j in 1:ncol(test)) {
    if (abs(test[i,j]) > mythresh) {
      l_ct_r <- rownames(test)[i]
      target <- colnames(test)[j]
      target <- strsplit(target,split='_')[[1]][[1]]
      mychannel <- paste0(l_ct_r,"_",target)
      all_channels <- c(all_channels,mychannel)
    }
  }
}
all_channels <- unique(all_channels)
length(all_channels)

# now comparing them to the full lr results
all_df_net <- readRDS(file='/home/jmitchel/data/lupus_data/LR_cell_chat_res_d_sep_list.rds')
all_df_net_comb <- readRDS(file='/home/jmitchel/data/lupus_data/LR_cell_chat_res_d_sep.rds')



all_df_net_comb <- all_df_net_comb[,c(1,2,3,4,7)]
## need to switch source/taget names to match mine since I changed them...
all_df_net_comb$source <- sapply(as.character(all_df_net_comb$source),function(x){
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
all_df_net_comb$target <- sapply(as.character(all_df_net_comb$target),function(x){
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
all_df_net_comb$channel <- sapply(1:nrow(all_df_net_comb),function(i){
  return(paste0(all_df_net_comb$ligand[i],'_',all_df_net_comb$source[i],
                '_',all_df_net_comb$receptor[i],'_',all_df_net_comb$target[i]))

})

dim(all_df_net_comb)

all_df_net_comb <- unique(all_df_net_comb)
dim(all_df_net_comb)

cchat_channels <- all_df_net_comb$channel

chan_both <- intersect(all_channels,cchat_channels)

scitd_only <- all_channels[!(all_channels %in% chan_both)]

cchat_only <- cchat_channels[!(cchat_channels %in% chan_both)]

length(chan_both)
length(cchat_only)
length(scitd_only)

chan_both
scitd_only
cchat_only





# add conos object for cell proportion analysis
con <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_conos2.rds')
pbmc_container$embedding <- con

# # in case below fn errors in the middle, need to save these objects
# orig_embed <- pbmc_container$embedding[["embedding"]]
# orig_clusts <- pbmc_container$embedding$clusters$leiden$groups


# to recover original embedding/cell assignments
pbmc_container$embedding[["embedding"]] <- orig_embed
pbmc_container$embedding$clusters$leiden$groups <- orig_clusts

pbmc_container$subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')


pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='Th',res=.6,n_col=2,alt_name='T4')
pdf(file = "/home/jmitchel/figures/for_paper_v2/Th_subc.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()

pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='cMono',res=.5,n_col=2,alt_name='cM')
pdf(file = "/home/jmitchel/figures/for_paper_v2/cMono_subc.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()

pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='NK',res=.6,n_col=2)
pdf(file = "/home/jmitchel/figures/for_paper_v2/NK_subc.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()

pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='Tc',res=.6,n_col=2,alt_name='T8')
pdf(file = "/home/jmitchel/figures/for_paper_v2/Tc_subc.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()


pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='ncMono',res=.6,n_col=2,alt_name='ncM')
pdf(file = "/home/jmitchel/figures/for_paper_v2/ncMono_subc.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()

pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='cDC',res=.5,n_col=2)
pdf(file = "/home/jmitchel/figures/for_paper_v2/cDC_subc.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()

pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='B',res=.8,n_col=2)
pdf(file = "/home/jmitchel/figures/for_paper_v2/B_subc.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()


## now getting full ctype prop associations
pbmc_container <- get_ctype_prop_associations(pbmc_container,'adj_pval',n_col=2)
pdf(file = "/home/jmitchel/figures/for_paper_v2/major_ctype_props.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()



## plotting F6 (sex-associated factor) associations with cell proportion shifts
f_6_subc_pv <- c(.82,.09,.075,1,.43,1,.054,.51)
f_6_subc_ct <- c('NK','Th','Tc','ncMono','cMono','cDC','B','Major ctypes')
f_6_ct_subc <- cbind.data.frame(-log10(f_6_subc_pv),f_6_subc_ct)
colnames(f_6_ct_subc) <- c('adj_pval','ctype')
f_6_ct_subc$ctype <- factor(f_6_ct_subc$ctype,levels=f_6_subc_ct)

p <- ggplot(f_6_ct_subc, aes(x = ctype, y = adj_pval)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(.05), 
             color = "red", size=1.5) +
  ylab("-log10(padj)") +
  xlab('Subclusters from this cell-type') +
  ggtitle('Factor 6\nsubcluster shift associations') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file = "/home/jmitchel/figures/for_paper_v2/f6_subc_all.pdf", useDingbats = FALSE,
    width = 5, height = 3.5)
p
dev.off()
  
  

# running stability analysis
pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(7,20,7),n_iterations=500,subset_type='subset', sub_prop=.75)

pdf(file = "/home/jmitchel/figures/for_paper_v2/stability_dsc.pdf", useDingbats = FALSE,
    width = 4, height = 3.5)
pbmc_container$plots$stability_plot_dsc
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/stability_lds.pdf", useDingbats = FALSE,
    width = 4, height = 3.5)
pbmc_container$plots$stability_plot_lds
dev.off()




## looking at location of FOXP3 expression in T4 cluster only
ct_res <- paste0('T4',':',as.character(0.6))
resolution_name <- paste0('res:',as.character(0.6))
subclusts <- container$subclusters[['T4']][[resolution_name]]

# append large cell type name to subclusters
subclusts <- sapply(subclusts,function(x){paste0('T4','_',x)})

# save original embedding
orig_embed <- con[["embedding"]]

# save original cluster labels
orig_clusts <- con$clusters$leiden$groups

# limit cells in subclusts to those that we actually have scores for
con[["embedding"]] <- orig_embed[names(subclusts),]


pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 4, height = 3.5)
con$plotGraph(gene='FOXP3')
dev.off()





## trying new lr method

























