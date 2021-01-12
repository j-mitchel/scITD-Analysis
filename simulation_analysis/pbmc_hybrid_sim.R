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

# save decomposition in object for comparison later
tucker_res1 <- pbmc_container$tucker_results
meta_anno1 <- pbmc_container$meta_associations



### now get results with added batch effect in log-space
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

# or make hybrid batch sim at log count level
pbmc_container <- get_hybrid_sim(pbmc_container,nbatches=2)

# finish prepping the data
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=500)

# get decomposition (to one more factor than needed)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,9,5), shuffle=FALSE)

# get meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes','batch'))

pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes','batch'),
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)
pbmc_container$plots$donor_matrix

# save decomposition results
tucker_res2 <- pbmc_container$tucker_results
meta_anno2 <- pbmc_container$meta_associations

# now get decomposition with combat batch correction
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5), shuffle=FALSE, batch_var='batch')
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes','batch'))
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes','batch'),
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)
pbmc_container$plots$donor_matrix

tucker_res3 <- pbmc_container$tucker_results
meta_anno3 <- pbmc_container$meta_associations

# compare decompositions
decomp_names <- c('no batch','with batch')
compare_decompositions(tucker_res1,tucker_res2,decomp_names,meta_anno1,meta_anno2)

decomp_names <- c('no batch','batch + combat')
compare_decompositions(tucker_res1,tucker_res3,decomp_names,meta_anno1,meta_anno3)

decomp_names <- c('with batch','batch + combat')
compare_decompositions(tucker_res2,tucker_res3,decomp_names,meta_anno2,meta_anno3)




### now get results with added batch effect in count-space
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

# make hybrid batch sim at count level
pbmc_container <- get_hybrid_sim(pbmc_container,counts=pbmc_counts,nbatches=2)

# finish prepping the data
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=500)

# get decomposition (to one more factor than needed)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,9,5), shuffle=FALSE)

# get meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes','batch'))

pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes','batch'),
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)
pbmc_container$plots$donor_matrix

tucker_res2 <- pbmc_container$tucker_results
meta_anno2 <- pbmc_container$meta_associations

# now get decomposition with combat batch correction
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5), shuffle=FALSE, batch_var='batch')
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes','batch'))
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes','batch'),
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)
pbmc_container$plots$donor_matrix

tucker_res3 <- pbmc_container$tucker_results
meta_anno3 <- pbmc_container$meta_associations


# compare decompositions
decomp_names <- c('no batch','with batch')
compare_decompositions(tucker_res1,tucker_res2,decomp_names,meta_anno1,meta_anno2)

decomp_names <- c('no batch','batch + combat')
compare_decompositions(tucker_res1,tucker_res3,decomp_names,meta_anno1,meta_anno3)

decomp_names <- c('with batch','batch + combat')
compare_decompositions(tucker_res2,tucker_res3,decomp_names,meta_anno2,meta_anno3)





### now get results with manually-specified batches in count-space
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

# make hybrid batch sim at count level
# complete batch confounding, equal batch sizes
mybatches <- list(
  sapply(c(40,38,45,7,1,2,3,8,9,12,13,14,18,19,23,24,25,29,30,34,35,36,41),function(x){paste0('s',as.character(x))}),
  sapply(c(4,5,6,10,11,15,16,17,20,21,22,26,27,28,31,32,33,37,39,42,43,44),function(x){paste0('s',as.character(x))})
)
pbmc_container <- get_hybrid_sim(pbmc_container,counts=pbmc_counts,nbatches=2,batch_override=mybatches)

# finish prepping the data
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="norm_var", thresh=500)

# get decomposition (to one more factor than needed)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,9,5), shuffle=FALSE)

# get meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes','batch'))

pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes','batch'),
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)
pbmc_container$plots$donor_matrix

tucker_res2 <- pbmc_container$tucker_results
meta_anno2 <- pbmc_container$meta_associations

# now get decomposition with combat batch correction
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5), shuffle=FALSE, batch_var='batch')
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes','batch'))
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes','batch'),
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)
pbmc_container$plots$donor_matrix

tucker_res3 <- pbmc_container$tucker_results
meta_anno3 <- pbmc_container$meta_associations


# compare decompositions
decomp_names <- c('no batch','with batch')
compare_decompositions(tucker_res1,tucker_res2,decomp_names,meta_anno1,meta_anno2)

decomp_names <- c('no batch','batch + combat')
compare_decompositions(tucker_res1,tucker_res3,decomp_names,meta_anno1,meta_anno3)

decomp_names <- c('with batch','batch + combat')
compare_decompositions(tucker_res2,tucker_res3,decomp_names,meta_anno2,meta_anno3)



