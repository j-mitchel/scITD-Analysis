


# read in the processed data
gtex_tpm_sub_transform <- readRDS(file='/home/jmitchel/data/gtex/7_tissues_counts.rds')
gtex_meta <- readRDS(file='/home/jmitchel/data/gtex/7_tissues_meta.rds')
feature.names.final <- readRDS(file='/home/jmitchel/data/gtex/genes.rds')

# make counts matrix into sparse matrix
gtex_tpm_sub_transform <- Matrix(gtex_tpm_sub_transform, sparse = TRUE)

ttypes_use <- unique(gtex_meta$ctypes)

# set up project parameters
param_list <- initialize_params(ctypes_use = ttypes_use,
                                ncores = 30, rand_seed = 10)


gtex_container <- make_new_container(count_data=gtex_tpm_sub_transform, meta_data=gtex_meta,
                                     gn_convert = feature.names.final, params=param_list,
                                     label_donor_sex = FALSE)

## need to do a modified pipeline for tensor formation because already 'pseudobulked' and normalized
gtex_container <- parse_data_by_ctypes(gtex_container)

# set donor_min_cells to 0 because each donor only has 1 sample per tissue type
gtex_container <- clean_data(gtex_container, donor_min_cells=0, gene_min_cells=20)

# still have to run this to get norm counts in the right spot
gtex_container <- get_pseudobulk(gtex_container)

for (ct in gtex_container$experiment_params$ctypes_use) {
  scMinimal <- gtex_container$scMinimal_ctype[[ct]]
  scMinimal$pseudobulk <- Matrix(scMinimal$pseudobulk,sparse=TRUE)
  scMinimal$pseudobulk <- t(scMinimal$pseudobulk)
}

# rest of pipeline is like normal
gtex_container <- get_normalized_variance(gtex_container)
gtex_container <- get_ctype_vargenes(gtex_container, method='norm_var', thresh=1000)
gtex_container <- scale_variance(gtex_container,var_scale_power=1.5)
gtex_container <- stack_tensor(gtex_container)


gtex_container <- run_tucker_ica(gtex_container, ranks=c(12,20,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
gtex_container <- run_tucker_ica(gtex_container, ranks=c(15,30,7),
                                 tucker_type = 'regular', rotation_type = 'ica')

# get factor-meta data associations
gtex_container <- get_meta_associations(gtex_container,vars_test=c('sex','age','dthhrdy'))
# include death type and age here as they are of interest...

# plot donor scores
# gtex_container <- plot_donor_matrix(gtex_container, meta_vars=c('sex'),
#                                     show_donor_ids = TRUE,
#                                     add_meta_associations=TRUE)
gtex_container <- plot_donor_matrix(gtex_container, meta_vars=c('sex','age','dthhrdy'),
                                    show_donor_ids = FALSE,
                                    cluster_by_meta='dthhrdy',
                                    add_meta_associations=TRUE)

gtex_container$plots$donor_matrix


# get assistance with rank determination
gtex_container <- determine_ranks_tucker(gtex_container, max_ranks_test=c(30,40,5),
                                         shuffle_level='tensor', shuffle_within=NULL,
                                         num_iter=10, batch_var=NULL,
                                         norm_method='trim',
                                         scale_factor=10000,
                                         scale_var=TRUE,
                                         var_scale_power=1.5)

gtex_container$plots$rank_determination_plot


gtex_container <- plot_scores_by_meta(gtex_container,'dthhrdy')

# pdf(file = "/home/jmitchel/figures/for_paper/lupus_status_associations.pdf", useDingbats = FALSE,
#     width = 11, height = 7)
gtex_container$plots$indv_meta_scores_associations
# dev.off()


# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(12,20,7), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')

# get loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=0.02,
                                           display_genes=FALSE,
                                           gene_callouts=TRUE,
                                           callout_n_gene_per_ctype=8)














