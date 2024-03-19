library(Seurat)
library(ggplot2)
library(coda.base)
library(dplyr)
library(devtools)
load_all('/home/jmitchel/scITD/')

### showing results stability for downsampling donors to x% of original

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


set.seed(0)
downsample_proportions <- c(.1,.2,.3,.4,.5,.6,.7,.8,.9)
num_donors <- length(pbmc_container$tensor_data[[1]])
mean_stability_dsc <- c()
mean_stability_lds <- c()
downsample_proportions_iter <- c() 
for (dprop in downsample_proportions) {
  pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(7,20),
                                           n_iterations=50,
                                           subset_type='subset', sub_prop=dprop)
  
  stability_res <- pbmc_container$stability_results
  iter_ndx <- sapply(1:50,function(x){
    return(rep(x,7))
  })
  stability_res$iteration <- c(iter_ndx)
  
  # compute average max correlations over factors
  stability_res_means <- stability_res %>%
    group_by(iteration) %>% 
    summarise(av_dsc=mean(dscores), av_lds=mean(ldngs))
  
  # store means
  mean_stability_dsc <- c(mean_stability_dsc,stability_res_means$av_dsc)
  mean_stability_lds <- c(mean_stability_lds,stability_res_means$av_lds)
  num_donors_down <- round(num_donors * dprop)
  downsample_proportions_iter <- c(downsample_proportions_iter,rep(num_donors_down,rep(length(stability_res_means$av_lds))))
}

# plot results
tmp <- cbind.data.frame(downsample_proportions_iter,mean_stability_dsc,mean_stability_lds)
colnames(tmp) <- c('prop','dsc_stability','lds_stability')
tmp2 <- tmp %>%
  group_by(prop) %>% 
  summarise(av_dsc=mean(dsc_stability), sd_dsc=sd(dsc_stability),
            av_lds=mean(lds_stability), sd_lds=sd(lds_stability))

all_means <- c(tmp2$av_dsc,tmp2$av_lds)
all_sds <- c(tmp2$sd_dsc,tmp2$sd_lds)
all_means <- as.data.frame(all_means)
all_sds <- as.data.frame(all_sds)
all_means$var_type <- c(rep('dsc',length(tmp2$av_dsc)),rep('lds',length(tmp2$av_lds)))
all_means$prop <- c(tmp2$prop,tmp2$prop)
tmp2 <- cbind.data.frame(all_means,all_sds)
colnames(tmp2) <- c('mean_cor','var_type','prop','sd_cor')
p <- ggplot(tmp2,aes(x=prop,y=mean_cor,color=var_type)) +
  geom_point() +
  geom_errorbar(aes(ymin=mean_cor-sd_cor, ymax=mean_cor+sd_cor), width=.2) +
  geom_line() +
  xlab('# of donors') +
  ylab('Correlation to original factors') +
  theme_bw()

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/downsample_analysis.pdf", useDingbats = FALSE,
#     width = 5, height = 3.5)
p
dev.off()






