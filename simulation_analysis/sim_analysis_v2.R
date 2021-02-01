
library(SymSim)

# first get ground truth counts
phyla2 <- read.tree('/home/jmitchel/scITD-Analysis/simulation_analysis/newick_tree_medium.txt')
# phyla2 <- read.tree('/home/jmitchel/scITD-Analysis/simulation_analysis/newick_tree_easy.txt')
# phyla2 <- read.tree('/home/jmitchel/scITD-Analysis/simulation_analysis/newick_tree_hard.txt')

pdf(file = "/home/jmitchel/figures/for_paper/sim_phyla.pdf", useDingbats = FALSE,
    width = 4, height = 4)
plot(phyla2)
dev.off()

ngenes <- 500
# # works well with hard but clusters are too separate
# true_counts_res <- SimulateTrueCounts(ncells_total=4000, min_popsize=600, 
#                                       ngenes=ngenes, nevf=60, 
#                                       evf_type="discrete", n_de_evf=60, 
#                                       vary="s", Sigma=0.05, phyla=phyla2, 
#                                       randseed=2)
# trying next with medium
true_counts_res <- SimulateTrueCounts(ncells_total=4000, min_popsize=600, 
                                      ngenes=ngenes, nevf=60, 
                                      evf_type="discrete", n_de_evf=45, 
                                      vary="s", Sigma=0.5, phyla=phyla2, 
                                      randseed=2)

true_counts_res_dis <- true_counts_res
tsne_true_counts <- PlotTsne(meta=true_counts_res[[3]], 
                             data=log2(true_counts_res[[1]]+1), 
                             evf_type="discrete", n_pc=20, label='pop', 
                             saving = F, 
                             plotname="discrete populations (true counts)")
tsne_true_counts[[2]]

# now get "observed" data
data(gene_len_pool)
gene_len <- sample(gene_len_pool, ngenes, replace = FALSE)
observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], 
                                       meta_cell=true_counts_res[[3]], 
                                       protocol="UMI", alpha_mean=0.25, 
                                       alpha_sd=0.02, gene_len=gene_len, 
                                       depth_mean=5e4, depth_sd=3e3)
tsne_UMI_counts <- PlotTsne(meta=observed_counts[[2]], 
                            data=normalize_counts(observed_counts[[1]]), 
                            evf_type="discrete", n_pc=20, label='pop', 
                            saving = F, plotname="observed counts UMI")
tsne_UMI_counts[[2]]


# plot as umap instead of tsne
library(umap)
mydat <- normalize_counts(observed_counts[[1]])
myumap <- umap(t(mydat))
myumap <- cbind(myumap$layout,observed_counts[[2]]$pop)
myumap <- as.data.frame(myumap)
colnames(myumap) <- c('UMAP1','UMAP2','pop')
p <- ggplot(myumap,aes(x=UMAP1,y=UMAP2,color=as.factor(pop))) +
  geom_point() +
  labs(color = "Cell Type")


## create heatmap of genes DE between the subtypes of ct1
comps <- list(c(1,3),c(2,3),c(1,2))
comps <- list(c(1,3),c(2,3))
ndx_keep <- c()
for (j in 1:length(comps)) {
  DEinfo <- getDEgenes(true_counts_res = true_counts_res_dis, popA = comps[[j]][1], popB = comps[[j]][2])
  ndx_de <- which(abs(DEinfo[["logFC_theoretical"]]) > 0.8 & DEinfo[["nDiffEVF"]] > 0)
  ndx_keep <- c(ndx_keep,ndx_de)
}
ndx_keep <- unique(ndx_keep)

# normalize counts
data_norm <- normalize_counts(observed_counts[[1]])
data_sub <- data_norm[ndx_keep,]

meta <- as.factor(observed_counts[[2]]$pop)

cells_keep <- which(meta %in%  c(1,2,3))

meta <- meta[cells_keep]

data_sub <- data_sub[,cells_keep]

data_sub <- t(scale(t(data_sub)))

ha = HeatmapAnnotation(bar = as.factor(meta))

hmap1 <- Heatmap(data_sub,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        top_annotation=ha,
        column_split = as.factor(meta))


# or show heatmap of all genes
meta <- as.factor(observed_counts[[2]]$pop)

data_sub <- normalize_counts(observed_counts[[1]])
data_sub <- t(scale(t(data_sub)))

ha = HeatmapAnnotation(bar = as.factor(meta))

hmap2 <- Heatmap(data_sub,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        top_annotation=ha,
        column_split = as.factor(meta))








# assign cell types to cells
meta <- observed_counts[[2]][,c('cellid','pop')]
rownames(meta) <- meta$cellid
meta$cellid <- NULL
meta$ctypes <- sapply(meta$pop,function(x) {
  if (x %in% c(1,2,3)) {
    return('ct1')
  } else {
    return('ct2')
  }
})

# assign donors to cells - will do roughly equal number donors per DE group

# meta$donors <- sapply(meta$pop,function(x) {
#   if (x %in% c(1,4)) {
#     tmp <- sample(1:10,1)
#   } else if (x %in% c(2,5)) {
#     tmp <- sample(11:20,1)
#   } else {
#     tmp <- sample(21:30,1)
#   }
#   return(paste0('s',as.character(tmp)))
# })

meta$donors <- sapply(meta$pop,function(x) {
  if (x %in% c(1,4)) {
    tmp <- sample(1:5,1)
  } else if (x %in% c(2,5)) {
    tmp <- sample(6:10,1)
  } else {
    tmp <- sample(11:25,1)
  }
  return(paste0('s',as.character(tmp)))
})

meta$dgroup <- sapply(meta$pop,function(x) {
  if (x %in% c(1,4)) {
    tmp <- 2
  } else if (x %in% c(2,5)) {
    tmp <- 3
  } else {
    tmp <- 1
  }
  return(as.character(tmp))
})

table(meta$donors)

meta$donors <- as.factor(meta$donors)
meta$ctypes <- as.factor(meta$ctypes)
meta$dgroup <- as.factor(meta$dgroup)

mycounts <- Matrix(observed_counts[[1]], sparse = TRUE)

colnames(mycounts) <- rownames(meta)
rownames(mycounts) <- sapply(1:nrow(mycounts),function(x){
  paste0('gene',as.character(x))
})

### remaking umap to specify ctypes and donor groups
mydat <- normalize_counts(mycounts)
myumap <- umap(as.matrix(t(mydat)))
myumap <- cbind(myumap$layout,meta$dgroup,meta$ctypes)
myumap <- as.data.frame(myumap)
colnames(myumap) <- c('UMAP1','UMAP2','dgroup','ctypes')
myumap$dgroup <- as.factor(myumap$dgroup)
myumap$ctypes <- as.factor(myumap$ctypes)
p <- ggplot(myumap,aes(x=UMAP1,y=UMAP2,shape=ctypes,color=dgroup)) +
  geom_point(size=1.5) +
  labs(shape = "Cell Type", color = "Donor Group")

# save plot
pdf(file = "/home/jmitchel/figures/for_paper/sim_UMAP.pdf", useDingbats = FALSE,
    width = 4, height = 3)
p
dev.off()
###




mycounts <- readRDS('/home/jmitchel/data/sim_data/sim1_counts.rds')
meta <- readRDS('/home/jmitchel/data/sim_data/sim1_meta.rds')

# set up project parameters
param_list <- initialize_params(ctypes_use = c("ct1", "ct2"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(count_data=mycounts, meta_data=meta,
                                     params=param_list)

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=5,
                              norm_method='trim', scale_factor=10000, 
                              vargenes_method='norm_var', vargenes_thresh=500,
                              scale_var = TRUE, var_scale_power = .5)

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(2,4,2),
                                 tucker_type = 'regular', rotation_type = 'ica')

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = TRUE,
                                    add_meta_associations=FALSE)

pdf(file = "/home/jmitchel/figures/for_paper/sim_dscores.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container$plots$donor_matrix
dev.off()



## need to create a df of fold-changes for DE genes (one for each process)

# for process 1 (donors 1-5)
DEinfo <- getDEgenes(true_counts_res = true_counts_res_dis, popA = 1, popB = 3)
p1_ct1 <- abs(DEinfo[["logFC_theoretical"]])
DEinfo <- getDEgenes(true_counts_res = true_counts_res_dis, popA = 4, popB = 6)
p1_ct2 <- abs(DEinfo[["logFC_theoretical"]])
df1 <- as.data.frame(cbind(p1_ct1,p1_ct2))
colnames(df1) <- c('ct1','ct2')
rownames(df1) <- sapply(1:nrow(df1),function(x){
  paste0('gene',as.character(x))
})

# for process 2 (donors 6-10)
DEinfo <- getDEgenes(true_counts_res = true_counts_res_dis, popA = 2, popB = 3)
p2_ct1 <- abs(DEinfo[["logFC_theoretical"]])
DEinfo <- getDEgenes(true_counts_res = true_counts_res_dis, popA = 5, popB = 6)
p2_ct2 <- abs(DEinfo[["logFC_theoretical"]])
df2 <- as.data.frame(cbind(p2_ct1,p2_ct2))
colnames(df2) <- c('ct1','ct2')
rownames(df2) <- sapply(1:nrow(df2),function(x){
  paste0('gene',as.character(x))
})


# make sure order of dfs matches donors in factor 1/2 respectively!
de_lfc <- list(df1,df2)
de_lfc <- list(df2,df1)

pbmc_container <- run_jackstraw(pbmc_container, ranks=c(2,4,2), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')

pdf(file = "/home/jmitchel/figures/for_paper/sim_ldngs_f1.pdf", useDingbats = FALSE,
    width = 8, height = 3.5)
plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=FALSE, nonsig_to_zero=FALSE, annot='sig_genes',
                    pathways=NULL, sim_de_donor_group=de_lfc, sig_thresh=0.0001, display_genes=FALSE, 
                    gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, show_xlab=TRUE,
                    show_var_explained=FALSE, reset_other_factor_plots=FALSE, draw_plot=TRUE)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/sim_ldngs_f2.pdf", useDingbats = FALSE,
    width = 8, height = 3.5)
plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=FALSE, nonsig_to_zero=FALSE, annot='sig_genes',
                    pathways=NULL, sim_de_donor_group=de_lfc, sig_thresh=0.0001, display_genes=FALSE, 
                    gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, show_xlab=TRUE,
                    show_var_explained=FALSE, reset_other_factor_plots=FALSE, draw_plot=TRUE)
dev.off()



# accuracy testing
library(pROC)
## Compute accuracy in predicting DE genes for factor 1
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=1, c('ct1','ct2'))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
pred <- c(sig_df)

# order de ground truth genes same way
de_true1 <- df1[rownames(sig_df),]
de_true1 <- unlist(de_true1)
de_true1 <- abs(de_true1) > .6 # make into logical mask

print(auc(de_true1, pred)) # response, predictor

pdf(file = "/home/jmitchel/figures/for_paper/sim_auc_f1.pdf", useDingbats = FALSE,
    width = 4, height = 4)
pROC_obj <- roc(de_true1,pred,
                smoothed = TRUE,
                ci=FALSE,
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,
                print.auc=TRUE, show.thres=TRUE)
dev.off()


# now compute AUC for factor 2
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=2, c('ct1','ct2'))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
pred <- c(sig_df)

# order de ground truth genes same way
de_true2 <- df2[rownames(sig_df),]
de_true2 <- unlist(de_true2)
de_true2 <- abs(de_true2) > .6 # make into logical mask

print(auc(de_true2, pred)) # response, predictor

pdf(file = "/home/jmitchel/figures/for_paper/sim_auc_f2.pdf", useDingbats = FALSE,
    width = 4, height = 4)
pROC_obj <- roc(de_true2,pred,
                smoothed = TRUE,
                ci=FALSE,
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,
                print.auc=TRUE, show.thres=TRUE)
dev.off()






# save raw sim data so can reproduce it
# parameters used were medium tree as ((A:2,B:2,C:0):4,(D:2,E:2,F:0):4);
# ncells_total=4000, min_popsize=600, ngenes=ngenes, nevf=60, evf_type="discrete", n_de_evf=45, vary="s", Sigma=0.5, phyla=phyla2, randseed=2
saveRDS(mycounts,'/home/jmitchel/data/sim_data/sim1_counts.rds')
saveRDS(meta,'/home/jmitchel/data/sim_data/sim1_meta.rds')
saveRDS(true_counts_res,'/home/jmitchel/data/sim_data/sim1_true_counts.rds')
saveRDS(list(tsne_true_counts[[2]],tsne_UMI_counts[[2]],p),'/home/jmitchel/data/sim_data/dim_reduction_plts.rds')


# test out rank determination accuracy
# get assistance with rank determination
pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(7,10,2),
                                         shuffle_level='cells', shuffle_within=NULL,
                                         num_iter=50, batch_var=NULL,
                                         norm_method='trim',
                                         scale_factor=10000,
                                         scale_var=TRUE,
                                         var_scale_power=.5)

pdf(file = "/home/jmitchel/figures/for_paper/sim_rank_determination.pdf", useDingbats = FALSE,
    width = 7, height = 7)
pbmc_container$plots$rank_determination_plot
dev.off()



