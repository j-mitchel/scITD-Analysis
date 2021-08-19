
library(scITD)

# counts matrix
# pbmc_counts <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_counts_v2.rds')
pbmc_counts <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_counts_v2_winsorized_simple.rds')

# meta data matrix
pbmc_meta <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_meta_v2.rds')

# ensembl to gene name conversions
feature.names <- readRDS('/home/jmitchel/data/van_der_wijst/genes.rds')

# change names of ctypes to match those from sle dataset
pbmc_meta$ctypes <- sapply(as.character(pbmc_meta$ctypes),function(x){
  if (x=='CD4+ T') {
    return('Th')
  } else if (x=='cMonocyte') {
    return('cMono')
  } else if (x=='CD8+ T') {
    return('Tc')
  } else {
    return(x)
  }
})

pbmc_meta$ctypes <- as.factor(pbmc_meta$ctypes)

# set up project parameters
# param_list <- initialize_params(ctypes_use = c("CD4+ T", "CD8+ T", "cMonocyte", "CD56(dim) NK", "B"),
#                                 ncores = 30, rand_seed = 10)
param_list <- initialize_params(ctypes_use = c("Th", "Tc", "cMono", "CD56(dim) NK", "B"),
                                ncores = 30, rand_seed = 10)

# pbmc_container <- make_new_container(count_data=pbmc_counts, meta_data=pbmc_meta,
#                                      gn_convert = feature.names, params=param_list,
#                                      label_donor_sex = TRUE, winsorize_param=2)
# 
# saveRDS(pbmc_container$scMinimal_full$count_data,file='/home/jmitchel/data/van_der_wijst/pbmc_counts_v2_winsorized.rds')
# saveRDS(pbmc_counts,file='/home/jmitchel/data/van_der_wijst/pbmc_counts_v2_winsorized_simple.rds')

pbmc_container <- make_new_container(count_data=pbmc_counts, meta_data=pbmc_meta,
                                     gn_convert = feature.names, params=param_list,
                                     label_donor_sex = TRUE)

# pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=5,
#                               norm_method='trim', scale_factor=10000,
#                               vargenes_method='norm_var', vargenes_thresh=1000,
#                               scale_var = TRUE, var_scale_power = 2)
# 
# pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=5,
#                               norm_method='trim', scale_factor=10000,
#                               vargenes_method='norm_var', vargenes_thresh=500,
#                               scale_var = TRUE, var_scale_power = 2)
# 
# pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=5,
#                               norm_method='trim', scale_factor=10000,
#                               vargenes_method='anova', vargenes_thresh=.02,
#                               scale_var = TRUE, var_scale_power = 2)

# pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=5,
#                               norm_method='trim', scale_factor=10000,
#                               vargenes_method='norm_var_pvals', vargenes_thresh=.1,
#                               scale_var = TRUE, var_scale_power = 2)

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=5,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = 1.5)

# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,10,5),
#                                  tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,10,5),
                                 tucker_type = 'regular', rotation_type = 'ica')

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes'),
                                        stat_use='pval')
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex'),
                                        stat_use='pval')

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes'),
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes'),
                                    show_donor_ids = FALSE,
                                    cluster_by_meta='sex',
                                    add_meta_associations=TRUE)
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('lanes'),
                                    show_donor_ids = TRUE,
                                    cluster_by_meta='lanes',
                                    add_meta_associations=TRUE)

pbmc_container$plots$donor_matrix

pdf(file = "/home/jmitchel/figures/for_paper_v2/pbmc_dscores.pdf", useDingbats = FALSE,
    width = 6, height = 6.5)
pbmc_container$plots$donor_matrix
dev.off()


# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)


# get assistance with rank determination
pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(10,15,5),
  shuffle_level='tensor', shuffle_within=NULL,
  num_iter=10, batch_var=NULL,
  norm_method='trim',
  scale_factor=10000,
  scale_var=TRUE,
  var_scale_power=2)

pbmc_container$plots$rank_determination_plot

# pdf(file = "/home/jmitchel/figures/for_paper/pbmc_rank_determination_tensor.pdf", useDingbats = FALSE,
#     width = 7, height = 7)
# pbmc_container$plots$rank_determination_plot
# dev.off()

# optimize power for variance scaling
pbmc_container <- optimize_var_scale_power(pbmc_container,c(1,1,5),c(10,15,5),.5,2)
pbmc_container <- optimize_var_scale_power(pbmc_container,c(3,3,5),c(10,15,5),.5,2,tucker_type='sparse')
pbmc_container$plots$var_scale_plot


pbmc_container <- run_stability_analysis(pbmc_container, ranks=c(5,10,5), downsample_ratio=0.9, n_iter=1000,
                       norm_method='trim', scale_factor=10000,
                       batch_var=NULL, scale_var=TRUE,
                       var_scale_power=2, tucker_type='regular',
                       rotation_type='ica')

pbmc_container$plots$stability_plot

# pdf(file = "/home/jmitchel/figures/for_paper/pbmc_stability_analysis.pdf", useDingbats = FALSE,
#     width = 3.5, height = 4)
# pbmc_container$plots$stability_plot
# dev.off()


# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(5,10,5), n_fibers=25, n_iter=5000,
                                tucker_type='regular', rotation_type='ica')

# save significant genes
# saveRDS(pbmc_container[["gene_score_associations"]],file='/home/jmitchel/data/van_der_wijst/pbmc_jackstraw.rds')
pbmc_container[["gene_score_associations"]] <- readRDS(file='/home/jmitchel/data/van_der_wijst/pbmc_jackstraw.rds')


## get loadings plots (for paper)
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=.01,
                                           display_genes=F,
                                           gene_callouts = TRUE,
                                           callout_n_gene_per_ctype=5,
                                           show_var_explained = FALSE)


pdf(file = "/home/jmitchel/figures/for_paper_v2/pbmc_loadings_hbb_factor.pdf", useDingbats = FALSE,
    width = 4, height = 3.5)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["5"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["5"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/pbmc_loadings_xy2_factor.pdf", useDingbats = FALSE,
    width = 4, height = 4)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["3"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["3"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()


pdf(file = "/home/jmitchel/figures/for_paper/pbmc_loadings_xy2_factor.pdf", useDingbats = FALSE,
    width = 4, height = 4)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["2"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["2"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/pbmc_loadings_hbb_factor.pdf", useDingbats = FALSE,
    width = 4, height = 3.5)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["5"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["5"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()
##

pdf(file = "/home/jmitchel/figures/for_paper/pbmc_loadings_xy_oqe.pdf", useDingbats = FALSE,
    width = 4, height = 4)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["3"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["3"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/pbmc_loadings_hbb_oqe.pdf", useDingbats = FALSE,
    width = 4, height = 4)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["5"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["5"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/pbmc_loadings_IFN_oqe.pdf", useDingbats = FALSE,
    width = 4.5, height = 6)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["1"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["1"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()


pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=0.05,
                                           display_genes=FALSE,
                                           gene_callouts=TRUE,
                                           callout_n_gene_per_ctype=5)

myfig <- render_multi_plots(pbmc_container,data_type='loadings')
myfig

# pdf(file = "/home/jmitchel/figures/for_paper/pbmc_loadings.pdf", useDingbats = FALSE,
#     width = 18, height = 14)
# myfig
# dev.off()

# pdf(file = "/home/jmitchel/figures/for_paper/pbmc_f1_dsig_genes.pdf", useDingbats = FALSE,
#     width = 9, height = 13)
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=1,ctypes_use=c('cMonocyte','CD4+ T','CD8+ T','CD56(dim) NK'),
                                       top_n_per_ctype=c(20,10,10,20), show_donor_labels=TRUE)
dev.off()

# pdf(file = "/home/jmitchel/figures/for_paper/pbmc_f4_dsig_genes.pdf", useDingbats = FALSE,
#     width = 9, height = 13)
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=4,ctypes_use=c('CD8+ T','cMonocyte'),
                                       top_n_per_ctype=c(30,10), show_donor_labels=TRUE)
dev.off()


pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=1, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))

pdf(file = "/home/jmitchel/figures/for_paper/pbmc_f1_gsea_oqe.pdf", useDingbats = FALSE,
    width = 7, height = 8)
pbmc_container[["plots"]][["gsea"]][["1"]]
dev.off()

pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=4, method="fgsea", thresh=0.05,
                                      db_use=c("GO","Reactome","BioCarta"), collapse_paths=FALSE)

# pdf(file = "/home/jmitchel/figures/for_paper/pbmc_f4_gsea.pdf", useDingbats = FALSE,
#     width = 7, height = 8)
pbmc_container[["plots"]][["gsea"]][["4"]]
dev.off()

pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=5, method="fgsea", thresh=0.05,
                                      db_use=c("GO"), collapse_paths=FALSE)

pbmc_container[["plots"]][["gsea"]][["5"]]

plot_loadings_annot(pbmc_container, factor_select=5)

myfig <- render_multi_plots(pbmc_container,data_type='loadings')
myfig


pbmc_container <- get_subtype_prop_associations(pbmc_container,max_res=1,'adj_pval',integration_var='lanes')
pbmc_container <- get_subtype_prop_associations(pbmc_container,max_res=1,'adj_pval',n_col=3)

pdf(file = "/home/jmitchel/figures/for_paper/pbmc_subtype_associations.pdf", useDingbats = FALSE,
    width = 9, height = 7)
pbmc_container$plots$subtype_prop_factor_associations
dev.off()


# # save subclustering
# saveRDS(pbmc_container$subclusters,file='/home/jmitchel/data/van_der_wijst/pbmc_subcluster_data.rds')



pbmc_container <- get_ctype_prop_associations(pbmc_container,'adj_pval',n_col=3)
pdf(file = "/home/jmitchel/figures/for_paper/pbmc_major_ctype_associations.pdf", useDingbats = FALSE,
    width = 7, height = 7)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()

pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('CD4+ T','CD56(dim) NK',
                                                  'cMonocyte','B'),
                                         res=c(.6,.6,.6,.6),
                                         factors=c(1,1,1,1))

subc_fig <- render_subtype_plots_v2(pbmc_container)
pdf(file = "/home/jmitchel/figures/for_paper/pbmc_f1_subtypes.pdf", useDingbats = FALSE,
    width = 21, height = 16)
subc_fig
dev.off()


pbmc_container[["plots"]][["subc_plots"]] <- NULL
pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('CD8+ T'),
                                         res=c(.5),
                                         factors=c(4))

subc_fig <- render_subtype_plots_v2(pbmc_container)
pdf(file = "/home/jmitchel/figures/for_paper/pbmc_f4_subtypes.pdf", useDingbats = FALSE,
    width = 4, height = 11)
subc_fig
dev.off()


print('test')














### old pipeline

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

# run tucker
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5), shuffle=FALSE)

# or can first get rid of batch effect associated with 10x lanes
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,8,5), shuffle=FALSE, batch_var='lanes')

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex','lanes'),
                                    show_donor_ids = TRUE,
                                    add_meta_associations=TRUE)
pbmc_container$plots$donor_matrix

# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, n_fibers=100, n_iter=500)

# show that significant genes correspond to large magnitude loadings
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=FALSE,
                                           nonsig_to_zero=FALSE,
                                           annot='sig_genes',
                                           sig_thresh=0.05,
                                           display_genes=FALSE,
                                           gene_callouts=FALSE)
render_all_lds_plots(pbmc_container, n_rows=2)


# plot loadings with significant genes only
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE, annot='none',
                                           sig_thresh=0.05, display_genes=FALSE,
                                           gene_callouts=TRUE)
render_all_lds_plots(pbmc_container, n_rows=2)

# plot loadings with significant genes only and gene callouts
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE, annot='none',
                                           sig_thresh=0.05, display_genes=FALSE,
                                           gene_callouts=TRUE,
                                           callout_n_gene_per_ctype=5,
                                           callout_ctypes=list(c(NULL),
                                                               c('CD4+ T',
                                                                 'cMonocyte'),
                                                               c(NULL),
                                                               c('CD8+ T'),
                                                               c(NULL)))
render_all_lds_plots(pbmc_container, n_rows=2)

### plot_donor_sig_genes and run_gsea_one_factor should be made automatic for all
### all factors and resulting plots should be auto formatted...
# generate some plots of scaled expression for top loadings genes of a factor
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=2,
                                       top_n_per_ctype=c(30,30),
                                       ctypes_use=c('CD4+ T','cMonocyte'))
pbmc_container$plots$donor_sig_genes$Factor2

pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=3,
                                       top_n_per_ctype=c(5,30),
                                       ctypes_use=c('CD4+ T','CD8+ T'))
pbmc_container$plots$donor_sig_genes$Factor3

pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=4,
                                       top_n_per_ctype=c(5,20),
                                       ctypes_use=c('CD4+ T','CD8+ T'))
pbmc_container$plots$donor_sig_genes$Factor4


# run gsea for a few factors
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.01,
                    db_use="GO", collapse_paths=TRUE)
pbmc_container$plots$gsea$Factor2

pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.01,
                                      db_use="GO", collapse_paths=TRUE)
pbmc_container$plots$gsea$Factor4


# run cell subtype proportion analysis
pbmc_container <- get_subtype_prop_associations(pbmc_container,max_res=1.5,
                                                stat_type='adj_pval',
                                                integration_var='lanes')
container$plots$subtype_prop_factor_associations

# generate subcluster plots and association significance for individual subtypes
pbmc_container <- get_all_subclust_plots(pbmc_container,
                                         ctypes=c('CD4+ T','CD56(dim) NK',
                                                  'cMonocyte','CD8+ T'),
                                         res=c(.7,.6,.5,.5),
                                         factors=c(2,2,2,4))

render_subtype_plots(pbmc_container)


# get associations between proportions of major cell types and factors
pbmc_container <- get_ctype_prop_associations(pbmc_container,
                                              stat_type='adj_pval')
pbmc_container$plots$ctype_prop_factor_associations


# run stability analysis
pbmc_container <- run_stability_analysis(pbmc_container, downsample_ratio=0.9,
                                         n_iter=500)

# look at number of significant genes for additional factors
pbmc_container <- get_min_sig_genes(pbmc_container,donor_rank_range=c(4:9),thresh=0.05)
pbmc_container$plots$min_sig_genes





# testing out data reconstruction using new data
d1_dat <- pbmc_container[["tensor_data"]][[4]][3,,] #donor s40
rownames(d1_dat) <- pbmc_container$tensor_data[[2]]
colnames(d1_dat) <- pbmc_container$tensor_data[[3]]
f1_data <- get_one_factor(pbmc_container, factor_select=1)

head(f1_data[[2]])
head(d1_dat)

# order them the same way
d1_dat <- d1_dat[rownames(f1_data[[2]]),colnames(f1_data[[2]])]

# trying to normalize by frobenius norm
fn <- rTensor::fnorm(as.tensor(f1_data[[2]]))
f1_data[[2]] <- f1_data[[2]] / fn

# adjust donor scores as well
f1_data[[1]] <- f1_data[[1]] * fn

# multiply them together
new_score <- sum(d1_dat * f1_data[[2]])
print(new_score)



## see if I can get new scores doing pca
pb <- pbmc_container[["scMinimal_ctype"]][["CD4+ T"]][["pseudobulk"]]
mysvd <- svd(pb)
# u is left singular vector, d is singular values, v is right singular values

# calculate PCs
mypcs <- mysvd$u %*% diag(mysvd$d)
rownames(mypcs) <- rownames(pb)

# multiply expression times loadings for s1 pc1
expr <- pb[1,,drop=FALSE]
lds <- mysvd$v[,1,drop=FALSE]
lds <- t(lds)

sum(expr * lds)
# okay the logic seems to check out... 

# what happens if I multiply singular values by loadings, then try to norm with fnorm
mysvd$d[1]
lds_new <- diag(mysvd$d) %*% t(mysvd$v)
lds_new <- lds_new[1,,drop=FALSE]
# now what's the fnorm of this
rTensor::fnorm(as.tensor(lds_new))
# yup, the norm is exactly the singular value that I pulled out

# what if I multiply expression by the expanded right singular values (times the $d)
sum(expr * lds_new)
# nah this doesn't work as anticipated

## need to try an example with just regular tucker to see if it's ICA that's screwing things up
tnsr <- container$tensor_data[[4]]
tucker_decomp <- rTensor::tucker(rTensor::as.tensor(tnsr), ranks=c(5,10,5))
# create loadings matrix
gene_by_factors <- tucker_decomp$U[[2]]
ctype_by_factors <- tucker_decomp$U[[3]]
donor_mat <- tucker_decomp$U[[1]]
kron_prod <- kronecker(ctype_by_factors,gene_by_factors,make.dimnames = TRUE)
core_new <- t(as.matrix(donor_mat)) %*% rTensor::k_unfold(rTensor::as.tensor(tnsr),1)@data %*% kron_prod
ldngs <- core_new %*% t(kron_prod)

# I wonder if I actually need to do a per cell type normalization instead of using overall fnorm
# but then I'm not sure what I'd do with the donor scores...






## testing out if having more donors like S13 would make NK signature come out as separate factor
head(pbmc_meta)
s13_bcodes <- rownames(pbmc_meta)[pbmc_meta$donors=='s13']
s13_meta <- pbmc_meta[s13_bcodes,]
s13_counts <- pbmc_counts[,s13_bcodes]
new_names <- sapply(colnames(s13_counts),function(x){
  paste0(x,'_dup')
})
names(new_names) <- NULL
colnames(s13_counts) <- new_names
rownames(s13_meta) <- new_names

pbmc_counts2 <- cbind.data.frame(pbmc_counts,s13_counts)
pbmc_counts3 <- methods::as(as.matrix(pbmc_counts2),'sparseMatrix')

s13_meta$donors <- rep('s99',nrow(s13_meta))
pbmc_meta2 <- rbind.data.frame(pbmc_meta,s13_meta)













