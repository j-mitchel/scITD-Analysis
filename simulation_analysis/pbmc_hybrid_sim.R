library(scITD)

# counts matrix
pbmc_counts <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_counts_v2.rds')

# meta data matrix
pbmc_meta <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_meta_v2.rds')

# ensembl to gene name conversions
feature.names <- readRDS('/home/jmitchel/data/van_der_wijst/genes.rds')



param_list <- initialize_params(ctypes_use = c("CD4+ T", "CD8+ T", "cMonocyte", "CD56(dim) NK", "B"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(count_data=pbmc_counts, meta_data=pbmc_meta,
                                     gn_convert = feature.names, params=param_list,
                                     label_donor_sex = TRUE)
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=5,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var', vargenes_thresh=1000,
                              scale_var = TRUE, var_scale_power = 2)


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5),
                                 tucker_type = 'regular', rotation_type = 'ica')


# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)

pdf(file = "/home/jmitchel/figures/for_paper/hybrid_original_dscores.pdf", useDingbats = FALSE,
    width = 5.5, height = 7)
pbmc_container$plots$donor_matrix
dev.off()

# save decomposition in object for comparison later
tucker_res1 <- pbmc_container$tucker_results
meta_anno1 <- pbmc_container$meta_associations




# now generate a hybrid simulation
pbmc_container <- get_hybrid_sim(pbmc_container,nbatches=2)

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=5,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var', vargenes_thresh=1000,
                              scale_var = TRUE, var_scale_power = 2)


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,8,5),
                                 tucker_type = 'regular', rotation_type = 'ica')


# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes','batch'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)

pdf(file = "/home/jmitchel/figures/for_paper/hybrid_w_batch_dscores.pdf", useDingbats = FALSE,
    width = 5.5, height = 7)
pbmc_container$plots$donor_matrix
dev.off()

# save decomposition in object for comparison later
tucker_res2 <- pbmc_container$tucker_results
meta_anno2 <- pbmc_container$meta_associations


# now use combat to remove batch effects
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=5,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var', vargenes_thresh=1000,
                              scale_var = TRUE, var_scale_power = 2,
                              batch_var='batch')


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5),
                                 tucker_type = 'regular', rotation_type = 'ica')


# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes','batch'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)
pdf(file = "/home/jmitchel/figures/for_paper/hybrid_w_batch_combat_dscores.pdf", useDingbats = FALSE,
    width = 5.5, height = 7)
pbmc_container$plots$donor_matrix
dev.off()


# save decomposition in object for comparison later
tucker_res3 <- pbmc_container$tucker_results
meta_anno3 <- pbmc_container$meta_associations



# compare decompositions
decomp_names <- c('no batch','with batch')
pdf(file = "/home/jmitchel/figures/for_paper/hybrid_lin_compare_orig_w_batch.pdf", useDingbats = FALSE,
    width = 8, height = 4.5)
compare_decompositions(tucker_res1,tucker_res2,decomp_names,meta_anno1,meta_anno2)
dev.off()

decomp_names <- c('no batch','batch + combat')
pdf(file = "/home/jmitchel/figures/for_paper/hybrid_lin_compare_orig_w_batch_combat.pdf", useDingbats = FALSE,
    width = 8, height = 4.5)
compare_decompositions(tucker_res1,tucker_res3,decomp_names,meta_anno1,meta_anno3)
dev.off()

decomp_names <- c('with batch','batch + combat')
pdf(file = "/home/jmitchel/figures/for_paper/hybrid_lin_compare_w_batch_w_batch_combat.pdf", useDingbats = FALSE,
    width = 8, height = 4.5)
compare_decompositions(tucker_res2,tucker_res3,decomp_names,meta_anno2,meta_anno3)
dev.off()












# testing out nonlinear hybrid batch simulation
pbmc_container <- make_new_container(count_data=pbmc_counts, meta_data=pbmc_meta,
                                     gn_convert = feature.names, params=param_list,
                                     label_donor_sex = TRUE)
pbmc_container <- parse_data_by_ctypes(pbmc_container)
pbmc_container <- clean_data(pbmc_container, donor_min_cells=5, gene_min_cells=5)
pbmc_container <- get_pseudobulk(pbmc_container)
pbmc_container <- normalize_pseudobulk(pbmc_container, method='trim', scale_factor=10000)
pbmc_container <- get_hybrid_sim_nonlin(pbmc_container,nbatches=2)
pbmc_container <- get_normalized_variance(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method='norm_var', thresh=1000)
pbmc_container <- scale_variance(pbmc_container,var_scale_power=2)
# pbmc_container <- apply_combat(pbmc_container,batch_var='batch')
pbmc_container <- stack_tensor(pbmc_container)


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,8,5),
                                 tucker_type = 'regular', rotation_type = 'ica')


# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes','batch'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)

pdf(file = "/home/jmitchel/figures/for_paper/hybrid_nonlin_w_batch_dscores.pdf", useDingbats = FALSE,
    width = 5.5, height = 7)
pbmc_container$plots$donor_matrix
dev.off()

# save decomposition in object for comparison later
tucker_res2 <- pbmc_container$tucker_results
meta_anno2 <- pbmc_container$meta_associations


pbmc_container <- apply_combat(pbmc_container,batch_var='batch')
pbmc_container <- stack_tensor(pbmc_container)


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5),
                                 tucker_type = 'regular', rotation_type = 'ica')


# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes','batch'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)

pdf(file = "/home/jmitchel/figures/for_paper/hybrid_nonlin_w_batch_combat_dscores.pdf", useDingbats = FALSE,
    width = 5.5, height = 7)
pbmc_container$plots$donor_matrix
dev.off()

# save decomposition in object for comparison later
tucker_res3 <- pbmc_container$tucker_results
meta_anno3 <- pbmc_container$meta_associations





# compare decompositions
decomp_names <- c('no batch','with batch')
pdf(file = "/home/jmitchel/figures/for_paper/hybrid_nonlin_compare_orig_w_batch.pdf", useDingbats = FALSE,
    width = 8, height = 4.5)
compare_decompositions(tucker_res1,tucker_res2,decomp_names,meta_anno1,meta_anno2)
dev.off()

decomp_names <- c('no batch','batch + combat')
pdf(file = "/home/jmitchel/figures/for_paper/hybrid_nonlin_compare_orig_w_batch_combat.pdf", useDingbats = FALSE,
    width = 8, height = 4.5)
compare_decompositions(tucker_res1,tucker_res3,decomp_names,meta_anno1,meta_anno3)
dev.off()

decomp_names <- c('with batch','batch + combat')
pdf(file = "/home/jmitchel/figures/for_paper/hybrid_nonlin_compare_w_batch_w_batch_combat.pdf", useDingbats = FALSE,
    width = 8, height = 4.5)
compare_decompositions(tucker_res2,tucker_res3,decomp_names,meta_anno2,meta_anno3)
dev.off()



