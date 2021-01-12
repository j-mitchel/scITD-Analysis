library(scITD)

# counts matrix
pbmc_counts <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_counts.rds')

# meta data matrix
pbmc_meta <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_meta.rds')

# ensembl to gene name conversions
feature.names <- readRDS('/home/jmitchel/data/van_der_wijst/genes.rds')

### first get results with standard pipeline (no batch added in yet)
# put data in project container
pbmc_scMinimal <- instantiate_scMinimal(count_data=pbmc_counts, meta_data=pbmc_meta)
pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = c("CD4+ T", "CD8+ T", "cMonocyte", "CD56(dim) NK", "B"),
                                     gn_convert = feature.names, scale_var = TRUE,
                                     var_scale_power = 2,
                                     tucker_type = 'sparse', rotation_type = 'ica',
                                     ncores = 30, rand_seed = 10)

# get sex meta data
pbmc_container <- identify_sex_metadata(pbmc_container)

# finish prepping the data
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=500)

# run tucker
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5), shuffle=FALSE)

# get factor associations with metadata
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes'))

# look at decomposition donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes'),
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)
pbmc_container$plots$donor_matrix

tucker_res1 <- pbmc_container$tucker_results
meta_anno1 <- pbmc_container$meta_associations

# now get decomposition with combat batch correction
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5), shuffle=FALSE, batch_var='lanes')
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes'))
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes'),
                                    show_donor_ids = FALSE,
                                    add_meta_associations=TRUE)
pbmc_container$plots$donor_matrix

tucker_res2 <- pbmc_container$tucker_results
meta_anno2 <- pbmc_container$meta_associations

# compare the two decompositions
decomp_names <- c('no combat','with combat')
compare_decompositions(tucker_res1,tucker_res2,decomp_names,meta_anno1,meta_anno2)

