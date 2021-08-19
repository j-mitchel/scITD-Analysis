library(Seurat)

pbmc <- readRDS(file="/home/jmitchel/data/covid_data/COVID19_subsetted_seurat.rds")

table(pbmc@meta.data$majorType)

colnames(pbmc@meta.data)

# make age to be numeric instead of character
pbmc@meta.data$Age[1000]
pbmc@meta.data$Age[10000]
pbmc@meta.data$Age <- as.numeric(pbmc@meta.data$Age)

param_list <- initialize_params(ctypes_use = c("B","CD4","CD8","NK","Mono","DC","Mega","Plasma"),
                                ncores = 30, rand_seed = 10)

param_list <- initialize_params(ctypes_use = c("B","CD4","CD8","NK","Mono","DC","Mega"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(seurat_obj=pbmc,
                                     params=param_list,
                                     metadata_cols=c('PatientID',
                                                     "Sex",
                                                     "Age",
                                                     "majorType",
                                                     'Sample.type',
                                                     'CoVID.19.severity',
                                                     'Sample.time',
                                                     'Single.cell.sequencing.platform',
                                                     'Outcome',
                                                     'Comorbidities',
                                                     'COVID.19.related.medication.and.anti.microbials'),
                                     metadata_col_nm=c('donors',
                                                       'sex',
                                                       'age',
                                                       'ctypes',
                                                       'fresh_frozen',
                                                       'covid_severity',
                                                       'sample_time',
                                                       'sequencing_platform',
                                                       'outcome',
                                                       'comorbidities',
                                                       'medications'))


pbmc_container <- form_tensor(pbmc_container, donor_min_cells=10, gene_min_cells=0,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.0001,
                              batch_var = 'sequencing_platform',
                              scale_var = TRUE, var_scale_power = .5)

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(10,25,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(25,50,7),
                                 tucker_type = 'regular', rotation_type = 'ica')


# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','age','fresh_frozen','covid_severity',
                                                                   'sample_time','sequencing_platform',
                                                                   'outcome','comorbidities','medications'),
                                        stat_use='pval')

# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    cluster_by_meta = 'sex',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

pbmc_container$plots$donor_matrix
dev.off()







# # checking for RPS4Y1 differences between M/F
# meta <- pbmc_container$scMinimal_full$metadata[,c('donors','sex')]
# meta <- unique(meta)
# rownames(meta) <- meta$donors
# meta$donors <- NULL
# 
# pb <- pbmc_container$scMinimal_ctype[['B']]$pseudobulk[,'XIST']
# tmp <- cbind.data.frame(meta[names(pb),1],pb)
# colnames(tmp) <- c('sex','expr')
# ggplot(tmp,aes(x=as.factor(sex),y=expr)) +
#   geom_boxplot()

















