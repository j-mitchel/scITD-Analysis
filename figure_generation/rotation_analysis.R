library(scITD)
library(plyr)
library(ComplexHeatmap)
library(circlize)

# comparison of three decompositions using the same data tensor

##### first need to generate the tensor
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

# using full dataset for this analysis (not SLE-only)
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


##### run and store results of the 3 different rotations
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20),
                                 tucker_type = 'regular', rotation_type = 'hybrid')
hybrid <- pbmc_container$tucker_results

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(9,20),
                                 tucker_type = 'regular', rotation_type = 'ica_dsc') 
ica_dmat <- pbmc_container$tucker_results





##### get factor associations with core set of IFN genes...
plot_ifn_associations <- function(container,tuck_res,ifn_genes) {
  dscores <- tuck_res[[1]]
  all_factors <- c()
  rsq <- c()
  for (i in 1:ncol(dscores)) {
    for (gn in ifn_genes) {
      for (ct in container$experiment_params$ctypes_use) {
        expr <- container[["scMinimal_ctype"]][[ct]][["pseudobulk"]][,gn]
        tmp <- cbind.data.frame(expr,dscores[names(expr),i])
        colnames(tmp) <- c('g_exp','dsc')
        lmres <- lm(dsc~g_exp,data=tmp)
        lmres <- summary(lmres)
        rsq <- c(rsq,lmres$r.squared)
        all_factors <- c(all_factors,paste0('Factor ',as.character(i)))
      }
    }
  }
  
  # plot the results
  tmp <- cbind.data.frame(rsq,all_factors)
  colnames(tmp) <- c('rsq','facts')
  tmp$facts <- as.factor(tmp$facts)
  df <- ddply(tmp, c("facts"), summarize, Mean = mean(rsq), SD = sd(rsq))
  p <- ggplot(df, aes(x = facts, y = Mean)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = Mean, ymax = Mean + SD), width = 0) +
    ylab("IFN-factor association (r-squared)") +
    xlab('') +
    ylim(0,1) +
    theme_bw()
  
  return(p)
  
}


ifn_core <- c('HERC5', 'IFI27', 'IRF7', 'ISG15', 'LY6E', 'MX1', 'OAS2', 'OAS3',
              'RSAD2', 'USP18', 'GBP5') # from Davenport et al. (2018) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6195724/
p = plot_ifn_associations(pbmc_container,hybrid,ifn_core)

### Supplementary Note 1 Figure B top
# pdf(file = "/home/jmitchel/figures/for_paper_v2/hybrid_ifn2.pdf", useDingbats = FALSE,
#     width = 4.75, height = 3)
p
# dev.off()

p = plot_ifn_associations(pbmc_container,ica_dmat,ifn_core)

### Supplementary Note 1 Figure B bottom
# pdf(file = "/home/jmitchel/figures/for_paper_v2/ica_dmat_ifn2.pdf", useDingbats = FALSE,
#     width = 5.75, height = 3)
p
# dev.off()








##### now showing difference in association strength for the ethnicity-associated factor
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Ethnicity','sex')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

dscores <- ica_dmat[[1]]
tmp <- cbind.data.frame(dscores[,4],meta[rownames(dscores),'Ethnicity'])
colnames(tmp) <- c('dsc','Ethnicity')
mf_association <- ggplot(tmp,aes(x=Ethnicity,y=dsc)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.6, binwidth = .008) +
  ylab('Factor 4 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()

### Supplementary Note 1 Figure C bottom
# pdf(file = "/home/jmitchel/figures/for_paper_v2/ica_dmat_eth2.pdf", useDingbats = FALSE,
#     width = 3.5, height = 3.25)
mf_association
# dev.off()


# now for hybrid
dscores <- hybrid[[1]]
tmp <- cbind.data.frame(dscores[,6],meta[rownames(dscores),'Ethnicity'])
colnames(tmp) <- c('dsc','Ethnicity')
mf_association <- ggplot(tmp,aes(x=Ethnicity,y=dsc)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.65, binwidth = .015) +
  ylab('Factor 6 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()

### Supplementary Note 1 Figure C top
# pdf(file = "/home/jmitchel/figures/for_paper_v2/hybrid_eth2.pdf", useDingbats = FALSE,
#     width = 3.5, height = 3.25)
mf_association
# dev.off()















