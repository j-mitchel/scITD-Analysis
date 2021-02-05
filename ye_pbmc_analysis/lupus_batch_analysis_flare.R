


# load up the subsetted dataset
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_flare_v3.rds')

# set up project parameters
param_list <- initialize_params(ctypes_use = c("B","NK","T4","T8","cDC",
                                               "cM","ncM"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(seurat_obj=pbmc,
                                     params=param_list,
                                     metadata_cols=c('ind_cov_batch_cov',
                                                     "SLE_status",
                                                     "Status",
                                                     "cg_cov",
                                                     "sex",
                                                     "age",
                                                     "batch_cov",
                                                     "Processing_Cohort"),
                                     metadata_col_nm=c('donors',
                                                       'SLE_status',
                                                       'Status',
                                                       'ctypes',
                                                       'sex',
                                                       'age',
                                                       'pool',
                                                       'processing'))

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                              scale_var = TRUE, var_scale_power = 1.5)


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,12,6),
                                 tucker_type = 'regular', rotation_type = 'ica')
# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(8,16,6),
#                                  tucker_type = 'regular', rotation_type = 'ica')

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','pool','Status'))

# save decomposition in object for comparison later
tucker_res1 <- pbmc_container$tucker_results
meta_anno1 <- pbmc_container$meta_associations



# now do with combat batch correction
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                              scale_var = TRUE, var_scale_power = 1.5,
                              batch_var='pool')

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,10,6),
                                 tucker_type = 'regular', rotation_type = 'ica')
# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,12,6),
#                                  tucker_type = 'regular', rotation_type = 'ica')
# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,13,6),
#                                  tucker_type = 'regular', rotation_type = 'ica')

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','pool','Status'))


# save decomposition in object for comparison later
tucker_res2 <- pbmc_container$tucker_results
meta_anno2 <- pbmc_container$meta_associations





decomp_names <- c('no combat','with combat')
pdf(file = "/home/jmitchel/figures/for_paper/flare/flare_combat_compare.pdf", useDingbats = FALSE,
    width = 9, height = 5)
compare_decompositions(tucker_res1,tucker_res2,decomp_names,meta_anno1,meta_anno2)
dev.off()






