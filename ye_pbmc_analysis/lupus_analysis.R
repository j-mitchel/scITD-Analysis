
library(scITD)

pbmc_scMinimal <- seurat_to_scMinimal(pbmc,normalize_counts=FALSE, 
                                      metadata_cols=c('ind_cov_batch_cov',
                                                      "Genotype.ID",
                                                      "SLE_status",
                                                      "Status",
                                                      "cg_cov",
                                                      "sex",
                                                      "age",
                                                      "batch_cov"),
                                      metadata_col_nm=c('donors',
                                                        'donor_genotype',
                                                        'SLE_status',
                                                        'Status',
                                                        'ctypes',
                                                        'sex',
                                                        'age',
                                                        'pool'))


# get median number cells per donor 
ctypes <- levels(pbmc_scMinimal$metadata$ctypes)
for (ctype in ctypes) {
  tmp <- pbmc_scMinimal$metadata[pbmc_scMinimal$metadata$ctypes==ctype,]
  print(ctype)
  print(median(table(tmp$donors)))
}

# remove original pbmc object because it takes up a lot of RAM
rm(pbmc)
gc()

pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = c("B","NK","T4","T8","cDC",
                                                    "cM","ncM"),
                                     scale_var = TRUE,
                                     var_scale_power = 1.25,
                                     tucker_type = 'sparse', rotation_type = 'ica',
                                     ncores = 30, rand_seed = 10)

rm(pbmc_scMinimal)
gc()


# finish prepping the data
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="empir", thresh=0.01)

# determine appropriate variance scaling parameter
pbmc_container <- optimize_var_scale_power(pbmc_container, min_ranks_test=c(5,8,5),
                                           max_ranks_test=c(10,15,5),
                                           min_power_test=1.5,
                                           max_power_test=2.25)
pbmc_container$plots$var_scale_plot

# determine appropriate ranks to use for decomposition
pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(6,10,5),
                                         method='svd', num_iter=5, shuffle_level='cells')
pbmc_container$plots$rank_determination_plot


