
library(scITD)
library(Seurat)

# load up the subsetted dataset
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v2.rds')

# prep the container
pbmc_scMinimal <- seurat_to_scMinimal(pbmc,normalize_counts=FALSE,
                                      metadata_cols=c('ind_cov_batch_cov',
                                                      "Genotype.ID",
                                                      "SLE_status",
                                                      "Status",
                                                      "cg_cov",
                                                      "sex",
                                                      "age",
                                                      "batch_cov",
                                                      "Processing_Cohort"),
                                      metadata_col_nm=c('donors',
                                                        'donor_genotype',
                                                        'SLE_status',
                                                        'Status',
                                                        'ctypes',
                                                        'sex',
                                                        'age',
                                                        'pool',
                                                        'processing'))

pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = c("B","NK","T4","T8","cDC","cM","ncM"),
                                     scale_var = TRUE,
                                     var_scale_power = 1.5,
                                     tucker_type = 'regular', rotation_type = 'ica',
                                     ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = c("NK","T4","T8","cM"),
                                     scale_var = TRUE,
                                     var_scale_power = 1.5,
                                     tucker_type = 'regular', rotation_type = 'ica',
                                     ncores = 30, rand_seed = 10)

# finish prepping the data
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=500)

# run ICA to large number of factors
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(10,70,7), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(20,140,7), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(30,210,7), shuffle=FALSE)

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(10,30,7), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(20,30,7), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(30,60,7), shuffle=FALSE)


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(20,30,4), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(20,40,4), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(25,50,4), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(30,50,4), shuffle=FALSE)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(30,120,4), shuffle=FALSE)
# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(15,25,4), shuffle=FALSE)

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','pool'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status','pool','processing'),
                                    show_donor_ids = FALSE,add_meta_associations=TRUE)
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','Status','pool','processing'),
                                    cluster_by_meta='Status',
                                    show_donor_ids = FALSE,add_meta_associations=TRUE)
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,add_meta_associations=TRUE)
pbmc_container$plots$donor_matrix

# save decomposition in object for comparison later
tucker_res1 <- pbmc_container$tucker_results
meta_anno1 <- pbmc_container$meta_associations

# now decompose using combat
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(8,16,4), shuffle=FALSE, batch_var='pool')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(12,18,4), shuffle=FALSE, batch_var='pool')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(12,30,4), shuffle=FALSE, batch_var='pool')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(12,50,4), shuffle=FALSE, batch_var='pool')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(15,30,4), shuffle=FALSE, batch_var='pool')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(20,80,4), shuffle=FALSE, batch_var='pool')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(20,140,7), shuffle=FALSE, batch_var='pool')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(13,20,7), shuffle=FALSE, batch_var='pool')

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','pool'))

pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,add_meta_associations=TRUE)
pbmc_container$plots$donor_matrix

# save decomposition in object for comparison later
tucker_res2 <- pbmc_container$tucker_results
meta_anno2 <- pbmc_container$meta_associations

decomp_names <- c('no combat','with combat')
compare_decompositions(tucker_res1,tucker_res2,decomp_names,meta_anno1,meta_anno2)





