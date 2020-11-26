library(scITD)

sim_data <- get_sim_data(donors_total=50, cells_per_donor=250, n_processes=2,
                         donors_per_process=10, n_ctypes=2, n_genes=500, 
                         de_prob=0.07, de_strength=2, factor_overlap=TRUE, 
                         rseed=97854658)

sim_scMinimal <- instantiate_scMinimal(data_sparse = as.matrix(sim_data[[1]]), meta_data = sim_data[[2]])
sim_container <- make_new_container(sim_scMinimal,
                                    ctypes_use = c('ct1', 'ct2'),
                                    scale_var = T, var_scale_power = 1.25,
                                    tucker_type = 'sparse', 
                                    rotation_type = 'varimax',
                                    ncores = 30, rand_seed = 16)

# get ctype matrices
sim_container <- get_ctype_data(sim_container,make_clean=TRUE)


# show that rank determination method selects correct number of factors
sim_container <- determine_ranks_tucker(sim_container, max_ranks_test=c(6,10,5),
                                         method='svd', num_iter=5, shuffle_level='cells')
sim_container$plots$rank_determination_plot

# run tucker
sim_container <- run_tucker_ica(sim_container, ranks=c(2,4,2), shuffle=F)

# plot donor scores
sim_container <- plot_donor_matrix(sim_container, show_donor_ids = TRUE)
sim_container$plots$donor_matrix

# get significant genes
sim_container <- run_jackstraw(sim_container,n_fibers=20,n_iter=500)

# plot significant genes next to loadings and true de genes
sim_container <- get_all_lds_factor_plots(sim_container, use_sig_only=FALSE, 
                                           nonsig_to_zero=FALSE, 
                                           annot='sig_genes',
                                           sig_thresh=0.05, 
                                           sim_de_donor_group=c(2,1),
                                           display_genes=FALSE,
                                           gene_callouts=FALSE)
render_all_lds_plots(sim_container, n_rows=2)








