
library(scITD)
library(edgeR)

# counts matrix
pbmc_counts <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_counts.rds')

# meta data matrix
pbmc_meta <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_meta.rds')

# ensembl to gene name conversions
feature.names <- readRDS('/home/jmitchel/data/van_der_wijst/genes.rds')


### first generate plots using regular normalization method
# put data in project container
pbmc_scMinimal <- instantiate_scMinimal(count_data=pbmc_counts, meta_data=pbmc_meta)
pbmc_container <- make_new_container(pbmc_scMinimal,
                                     ctypes_use = c("CD4+ T", "CD8+ T", "cMonocyte", "CD56(dim) NK", "B"),
                                     gn_convert = feature.names, scale_var = TRUE,
                                     var_scale_power = 1.25,
                                     tucker_type = 'regular', rotation_type = 'ica',
                                     ncores = 30, rand_seed = 10)

# get sex meta data
pbmc_container <- identify_sex_metadata(pbmc_container)

# finish prepping the data
pbmc_container <- get_ctype_data(pbmc_container)
pbmc_container <- get_ctype_vargenes(pbmc_container, method="empir", thresh=0.01)

# run tucker
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5), shuffle=FALSE)

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes'),
                                    cluster_by_meta='sex', show_donor_ids = FALSE)
pbmc_container$plots$donor_matrix

# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, n_fibers=100, n_iter=500)

# show that significant genes correspond to large magnitude loadings
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE, 
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=0.05, 
                                           display_genes=FALSE,
                                           gene_callouts=FALSE)
render_all_lds_plots(pbmc_container, n_rows=2)


# run gsea
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=4, method="fgsea", thresh=0.05,
                                         db_use="GO", collapse_paths=TRUE)
pbmc_container$plots$gsea$Factor4


# get top n up loading genes in the IFN factor
tmp <- pbmc_container[["tucker_results"]][[2]][4,]
top_up_genes <- tmp[order(tmp,decreasing=TRUE)][1:100]
top_up_genes <- sapply(names(top_up_genes),function(x){
  strsplit(x,split=':')[[1]][[2]]
  })
top_up_genes <- convert_gn(pbmc_container,top_up_genes)
print(top_up_genes)






### now do a trimmed mean normalization preprocessing
all_nf <- calcNormFactors(pbmc_counts)

# divide by lib size and multiply by scale factor
lib_sizes <- Matrix::colSums(pbmc_counts)
pbmc_counts_tm <- sweep(pbmc_counts,MARGIN=2,lib_sizes*all_nf,FUN='/') * 10000

# log transform result
pbmc_counts_tm <- log1p(pbmc_counts_tm)

# put data in project container
pbmc_scMinimal_tm <- instantiate_scMinimal(data_sparse=pbmc_counts_tm, meta_data=pbmc_meta)
pbmc_container_tm <- make_new_container(pbmc_scMinimal_tm,
                                     ctypes_use = c("CD4+ T", "CD8+ T", "cMonocyte", "CD56(dim) NK", "B"),
                                     gn_convert = feature.names, scale_var = TRUE,
                                     var_scale_power = 1.25,
                                     tucker_type = 'regular', rotation_type = 'ica',
                                     ncores = 30, rand_seed = 10)

# get sex meta data
pbmc_container_tm <- identify_sex_metadata(pbmc_container_tm)

# finish prepping the data
pbmc_container_tm <- get_ctype_data(pbmc_container_tm)
pbmc_container_tm <- get_ctype_vargenes(pbmc_container_tm, method="empir", thresh=0.01)

# run tucker
pbmc_container_tm <- run_tucker_ica(pbmc_container_tm, ranks=c(5,8,5), shuffle=FALSE)

# plot donor scores
pbmc_container_tm <- plot_donor_matrix(pbmc_container_tm, meta_vars=c('sex','lanes'),
                                    cluster_by_meta='sex', show_donor_ids = FALSE)
pbmc_container_tm$plots$donor_matrix


# get significant genes
pbmc_container_tm <- run_jackstraw(pbmc_container_tm, n_fibers=100, n_iter=500)

# show that significant genes correspond to large magnitude loadings
pbmc_container_tm <- get_all_lds_factor_plots(pbmc_container_tm, use_sig_only=TRUE, 
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=0.05, 
                                           display_genes=FALSE,
                                           gene_callouts=FALSE)
render_all_lds_plots(pbmc_container_tm, n_rows=2)


# run gsea
pbmc_container_tm <- run_gsea_one_factor(pbmc_container_tm, factor_select=5, method="fgsea", thresh=0.05,
                                      db_use="GO", collapse_paths=TRUE)
pbmc_container_tm$plots$gsea$Factor5


# get top n up loading genes in the IFN factor
tmp <- pbmc_container_tm[["tucker_results"]][[2]][5,]
top_up_genes <- tmp[order(tmp,decreasing=TRUE)][1:100]
top_up_genes <- sapply(names(top_up_genes),function(x){
  strsplit(x,split=':')[[1]][[2]]
})
top_up_genes <- convert_gn(pbmc_container_tm,top_up_genes)
print(top_up_genes)


