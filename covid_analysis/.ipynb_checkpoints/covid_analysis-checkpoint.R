
pbmc <- readRDS(file="/home/jmitchel/data/covid_data/COVID19_subsetted_seurat.rds")

param_list <- initialize_params(ctypes_use = c("B","CD4","Mono"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(seurat_obj=pbmc,
                                     params=param_list,
                                     metadata_cols=c('PatientID',
                                                     "Sex",
                                                     "majorType"),
                                     metadata_col_nm=c('donors',
                                                       'sex',
                                                       'ctypes'))


pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=0,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5)


# checking for RPS4Y1 differences between M/F
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','sex')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

pb <- pbmc_container$scMinimal_ctype[['B']]$pseudobulk[,'XIST']
tmp <- cbind.data.frame(meta[names(pb),1],pb)
colnames(tmp) <- c('sex','expr')
ggplot(tmp,aes(x=as.factor(sex),y=expr)) +
  geom_boxplot()

















