library(splatter)
library(PRROC)

## trying new approach by adding two DE datasets together
all_donors <- sapply(1:40,function(x) {
  paste0('s',x)
})

dg1 <- all_donors[1:20]
dg2 <- all_donors[21:40]
dg3 <- all_donors[c(1:10,21:30)]
dg4 <- all_donors[c(11:20,31:40)]

params <- newSplatParams(nGenes = 3000)
params <- setParams(params, update = list(batchCells = 8000, seed=10))

scsim_de1 <- splatSimulateGroups(params,
                                group.prob = c(.5,.5),
                                de.prob = c(0,.035),
                                de.downProb = c(0,.5),
                                de.facLoc = .3,
                                verbose = FALSE)

rd1 <- as.data.frame(rowData(scsim_de1)@listData)
rownames(rd1) <- rd1$Gene
cd1 <- as.data.frame(colData(scsim_de1)@listData)
rownames(cd1) <- cd1$Cell


# assign cells to donors
cells_de1 <- rownames(cd1)[cd1$Group=='Group1']
cells_de2 <- rownames(cd1)[cd1$Group=='Group2']
d_assigns1 <- sample(dg1,length(cells_de1),replace=TRUE)
d_assigns2 <- sample(dg2,length(cells_de2),replace=TRUE)
cd1$donors <- NULL
cd1[cells_de1,'donors'] <- d_assigns1
cd1[cells_de2,'donors'] <- d_assigns2


# normalize by lib size
sim_counts1 <- counts(scsim_de1)
lib_sizes1 <- colSums(sim_counts1)
sim_counts_norm1 <- sweep(sim_counts1,2,lib_sizes1,'/')




## now simulate a slightly bigger population of cells to add as the second factor
params <- newSplatParams(nGenes = 3000)
params <- setParams(params, update = list(batchCells = 12000, seed=11))

scsim_de2 <- splatSimulateGroups(params,
                                 group.prob = c(.5,.5),
                                 de.prob = c(0,.03),
                                 de.downProb = c(0,.5),
                                 de.facLoc = .25,
                                 verbose = FALSE)

rd2 <- as.data.frame(rowData(scsim_de2)@listData)
rownames(rd2) <- rd2$Gene
cd2 <- as.data.frame(colData(scsim_de2)@listData)
rownames(cd2) <- cd2$Cell


# assign cells to donors
cells_de1 <- rownames(cd2)[cd2$Group=='Group1']
cells_de2 <- rownames(cd2)[cd2$Group=='Group2']
d_assigns1 <- sample(dg3,length(cells_de1),replace=TRUE)
d_assigns2 <- sample(dg4,length(cells_de2),replace=TRUE)
cd2$donors <- NULL
cd2[cells_de1,'donors'] <- d_assigns1
cd2[cells_de2,'donors'] <- d_assigns2

# normalize by lib size
sim_counts2 <- counts(scsim_de2)
lib_sizes2 <- colSums(sim_counts2)
sim_counts_norm2 <- sweep(sim_counts2,2,lib_sizes2,'/')


# check that all donors have more cells in cd2 than in cd1
sum(table(cd2$donors)<table(cd1$donors)) # should be 0

# now link cells from second sim to those in first sim via same donor labels
for (d in all_donors) {
  d_cells1 <- rownames(cd1)[cd1$donors==d]

  d_cells2 <- rownames(cd2)[cd2$donors==d]
  d_counts_norm2 <- sim_counts_norm2[,d_cells2]
  
  # subset cells from second dataset to match number in first
  d_counts_norm2 <- d_counts_norm2[,1:length(d_cells1)]
  colnames(d_counts_norm2) <- d_cells1
  
  # add second dataset to first
  sim_counts_norm1[,d_cells1] <- sim_counts_norm1[,d_cells1] + d_counts_norm2
}

# check that colSums are each equal to 2
colSums(sim_counts_norm1)[1:100]

# bring it back to count space by multiplying each cell by half it's library size from first dataset
sim_counts_final <- sweep(sim_counts_norm1,2,lib_sizes1/2,'*')

# check we get the right library sizes
colSums(sim_counts_final)[1:10]

counts_ct1 <- sim_counts_final
meta_ct1 <- cd1
de1_ct1 <- rd1
de2_ct1 <- rd2
########################################### now repeat this for cell type 2



params <- newSplatParams(nGenes = 3000)
params <- setParams(params, update = list(batchCells = 8000, seed=12))

scsim_de1 <- splatSimulateGroups(params,
                                 group.prob = c(.5,.5),
                                 de.prob = c(0,.035),
                                 de.downProb = c(0,.5),
                                 de.facLoc = .3,
                                 verbose = FALSE)

rd1 <- as.data.frame(rowData(scsim_de1)@listData)
rownames(rd1) <- rd1$Gene
cd1 <- as.data.frame(colData(scsim_de1)@listData)
rownames(cd1) <- cd1$Cell


# assign cells to donors
cells_de1 <- rownames(cd1)[cd1$Group=='Group1']
cells_de2 <- rownames(cd1)[cd1$Group=='Group2']
d_assigns1 <- sample(dg1,length(cells_de1),replace=TRUE)
d_assigns2 <- sample(dg2,length(cells_de2),replace=TRUE)
cd1$donors <- NULL
cd1[cells_de1,'donors'] <- d_assigns1
cd1[cells_de2,'donors'] <- d_assigns2


# normalize by lib size
sim_counts1 <- counts(scsim_de1)
lib_sizes1 <- colSums(sim_counts1)
sim_counts_norm1 <- sweep(sim_counts1,2,lib_sizes1,'/')




## now simulate a slightly bigger population of cells to add as the second factor
params <- newSplatParams(nGenes = 3000)
params <- setParams(params, update = list(batchCells = 12000, seed=13))

scsim_de2 <- splatSimulateGroups(params,
                                 group.prob = c(.5,.5),
                                 de.prob = c(0,.03),
                                 de.downProb = c(0,.5),
                                 de.facLoc = .25,
                                 verbose = FALSE)

rd2 <- as.data.frame(rowData(scsim_de2)@listData)
rownames(rd2) <- rd2$Gene
cd2 <- as.data.frame(colData(scsim_de2)@listData)
rownames(cd2) <- cd2$Cell


# assign cells to donors
cells_de1 <- rownames(cd2)[cd2$Group=='Group1']
cells_de2 <- rownames(cd2)[cd2$Group=='Group2']
d_assigns1 <- sample(dg3,length(cells_de1),replace=TRUE)
d_assigns2 <- sample(dg4,length(cells_de2),replace=TRUE)
cd2$donors <- NULL
cd2[cells_de1,'donors'] <- d_assigns1
cd2[cells_de2,'donors'] <- d_assigns2

# normalize by lib size
sim_counts2 <- counts(scsim_de2)
lib_sizes2 <- colSums(sim_counts2)
sim_counts_norm2 <- sweep(sim_counts2,2,lib_sizes2,'/')


# check that all donors have more cells in cd2 than in cd1
sum(table(cd2$donors)<table(cd1$donors)) # should be 0

# now link cells from second sim to those in first sim via same donor labels
for (d in all_donors) {
  d_cells1 <- rownames(cd1)[cd1$donors==d]
  
  d_cells2 <- rownames(cd2)[cd2$donors==d]
  d_counts_norm2 <- sim_counts_norm2[,d_cells2]
  
  # subset cells from second dataset to match number in first
  d_counts_norm2 <- d_counts_norm2[,1:length(d_cells1)]
  colnames(d_counts_norm2) <- d_cells1
  
  # add second dataset to first
  sim_counts_norm1[,d_cells1] <- sim_counts_norm1[,d_cells1] + d_counts_norm2
}

# check that colSums are each equal to 2
colSums(sim_counts_norm1)[1:100]

# bring it back to count space by multiplying each cell by half it's library size from first dataset
sim_counts_final <- sweep(sim_counts_norm1,2,lib_sizes1/2,'*')

# check we get the right library sizes
colSums(sim_counts_final)[1:10]

counts_ct2 <- sim_counts_final
meta_ct2 <- cd1
de1_ct2 <- rd1
de2_ct2 <- rd2



# combine data from the two cell types and add ctype meta column
colnames(counts_ct1) <- sapply(colnames(counts_ct1),function(x) {
  paste0(x,"_ct1")
})
names(colnames(counts_ct1)) <- NULL

colnames(counts_ct2) <- sapply(colnames(counts_ct2),function(x) {
  paste0(x,"_ct2")
})
names(colnames(counts_ct2)) <- NULL

rownames(meta_ct1) <- sapply(rownames(meta_ct1),function(x) {
  paste0(x,"_ct1")
})

rownames(meta_ct2) <- sapply(rownames(meta_ct2),function(x) {
  paste0(x,"_ct2")
})

combined_counts <- cbind(counts_ct1,counts_ct2)

# make counts a sparse matrix
combined_counts_sparse <- Matrix(combined_counts, sparse = TRUE) 

meta_ct1$ctypes <- 'ct1'
meta_ct2$ctypes <- 'ct2'
combined_meta <- rbind(meta_ct1,meta_ct2)


## plot umap coloring by cells or donors in groups
pbmc <- CreateSeuratObject(counts = combined_counts_sparse,
                           meta.data = combined_meta,
                           min.cells = 0, min.features = 0)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:15)

pdf(file = "/home/jmitchel/figures/for_paper/sim_v3_umap_ctypes.pdf", useDingbats = FALSE,
    width = 5, height = 5)
jpeg("/home/jmitchel/figures/for_paper/sim_v3_umap_ctypes.jpg")
DimPlot(pbmc, reduction = "umap",group.by = 'ctypes')
dev.off()

pbmc@meta.data$process1 <- sapply(1:nrow(pbmc@meta.data),function(i){
  donor <- pbmc@meta.data$donors[i]
  if (donor %in% dg1) {
    return('on')
  } else {
    return('off')
  }
})
pbmc@meta.data$process2 <- sapply(1:nrow(pbmc@meta.data),function(i){
  donor <- pbmc@meta.data$donors[i]
  if (donor %in% dg3) {
    return('on')
  } else {
    return('off')
  }
})

pdf(file = "/home/jmitchel/figures/for_paper/sim_v3_umap_proc1.pdf", useDingbats = FALSE,
    width = 5, height = 5)
jpeg("/home/jmitchel/figures/for_paper/sim_v3_umap_proc1.jpg")
DimPlot(pbmc, reduction = "umap",group.by = 'process1')
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/sim_v3_umap_proc2.pdf", useDingbats = FALSE,
    width = 5, height = 5)
jpeg("/home/jmitchel/figures/for_paper/sim_v3_umap_proc2.jpg")
DimPlot(pbmc, reduction = "umap",group.by = 'process2')
dev.off()


## edit dimplots by hand
tmp <- DimPlot(pbmc, reduction = "umap",group.by = 'ctypes')
tmp2 <- tmp +
  ggtitle('Cell Types') +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  scale_color_brewer(palette="Dark2") +
  labs(color='') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=15))

tmp2


tmp <- DimPlot(pbmc, reduction = "umap",group.by = 'process1')
tmp2 <- tmp +
  ggtitle('Process 1') +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  labs(color='') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=15))

tmp2

tmp <- DimPlot(pbmc, reduction = "umap",group.by = 'process2')
tmp2 <- tmp +
  ggtitle('Process 2') +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  labs(color='') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size=15))

tmp2



# set up project parameters
param_list <- initialize_params(ctypes_use = c("ct1", "ct2"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(count_data=combined_counts_sparse, meta_data=combined_meta,
                                     params=param_list)

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=1,
                              scale_var = TRUE, var_scale_power = .5)
# pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=5,
#                               norm_method='trim', scale_factor=10000,
#                               vargenes_method='norm_var', vargenes_thresh=750,
#                               scale_var = TRUE, var_scale_power = .5)

# pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(2,4,2),
#                                  tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(2,4),
                                 tucker_type = 'regular', rotation_type = 'ica_dsc')

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = TRUE)

# pdf(file = "/home/jmitchel/figures/for_paper_v2/sim_v3_dscores.pdf", useDingbats = FALSE,
#     width = 6, height = 7)
pbmc_container$plots$donor_matrix
dev.off()

pbmc_container <- run_jackstraw(pbmc_container, ranks=c(2,4,2), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')
# or use the lm method
pbmc_container <- get_lm_pvals(pbmc_container)


## making loadings plots when limiting to just most variable genes so it's easier to see them but will
# use all genes when computing accuracy of jackstraw in determining ground truth DE genes
pdf(file = "/home/jmitchel/figures/for_paper_v2/sim_ldngs_v3_f1.pdf", useDingbats = FALSE,
    width = 8, height = 3.5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=FALSE, annot='sig_genes',
                    pathways=NULL, sim_de_donor_group=list(de1_ct1,de1_ct2), sig_thresh=0.1, display_genes=FALSE, 
                    gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, show_xlab=TRUE,
                    show_var_explained=FALSE, reset_other_factor_plots=FALSE, draw_plot=TRUE)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/sim_ldngs_v3_f2.pdf", useDingbats = FALSE,
    width = 8, height = 3.5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=FALSE, annot='sig_genes',
                    pathways=NULL, sim_de_donor_group=list(de2_ct1,de2_ct2), sig_thresh=0.1, display_genes=FALSE, 
                    gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, show_xlab=TRUE,
                    show_var_explained=FALSE, reset_other_factor_plots=FALSE, draw_plot=TRUE)
dev.off()


## computing AUC values and getting plots
# do factor 1 first
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=1, c('ct1','ct2'))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
pred1 <- c(sig_df)

## trying using loadings instead
sig_df <- get_one_factor(pbmc_container,1)[[2]]
pred1 <- abs(c(sig_df))
##

ct1_de <- de1_ct1
ct1_de <- ct1_de[rownames(sig_df),]
ct2_de <- de1_ct2
ct2_de <- ct2_de[rownames(sig_df),]
de_true1 <- c(ct1_de$DEFacGroup2!=1,ct2_de$DEFacGroup2!=1)

library(pROC)

# pdf(file = "/home/jmitchel/figures/for_paper_v2/sim_v3_auc_f1.pdf", useDingbats = FALSE,
#     width = 4, height = 4)
pROC_obj <- roc(de_true1,pred1,
                smoothed = TRUE,
                ci=FALSE,
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,
                print.auc=TRUE, show.thres=TRUE)
dev.off()

## now do factor 2
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=2, c('ct1','ct2'))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
pred2 <- c(sig_df)

ct1_de <- de2_ct1
ct1_de <- ct1_de[rownames(sig_df),]
ct2_de <- de2_ct2
ct2_de <- ct2_de[rownames(sig_df),]
de_true2 <- c(ct1_de$DEFacGroup2!=1,ct2_de$DEFacGroup2!=1)

pdf(file = "/home/jmitchel/figures/for_paper_v2/sim_v3_auc_f2.pdf", useDingbats = FALSE,
    width = 4, height = 4)
pROC_obj <- roc(de_true2,pred2,
                smoothed = TRUE,
                ci=FALSE,
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,
                print.auc=TRUE, show.thres=TRUE)
dev.off()



# show rank determination works
pbmc_container <- determine_ranks_tucker(pbmc_container, max_ranks_test=c(7,10,2),
                                         shuffle_level='cells', shuffle_within=NULL,
                                         num_iter=50, batch_var=NULL,
                                         norm_method='trim',
                                         scale_factor=10000,
                                         scale_var=TRUE,
                                         var_scale_power=.5)

pdf(file = "/home/jmitchel/figures/for_paper_v2/sim_v3_rank_determination.pdf", useDingbats = FALSE,
    width = 7, height = 7)
pbmc_container$plots$rank_determination_plot
dev.off()





### testing AUC with downsampling
# get number of donors
num_donors <- length(all_donors)

cells_per_donor <- table(combined_meta$donors)

sizes_test <- c(15,20,40,60,80,100,120,140,160)
# sizes_test <- c(100,120,140,160)
downsample_sizes <- num_donors * 2 * sizes_test

f1_aucs <- c()
f2_aucs <- c()
cpd <- c()
for (ds in downsample_sizes) {
  # subsample data to correct mean size 
  prev_subs <- list()
  prev_ident <- FALSE
  for (myiter in 1:5) {
    min_cpd <- 0
    while (min_cpd < 5 | prev_ident) {
      prev_ident <- FALSE
      cells_sampled <- sample(rownames(combined_meta),ds)
      meta_sub <- combined_meta[cells_sampled,]
      cells_per_donor <- table(meta_sub[,c('donors','ctypes')])
      min_cpd <- min(cells_per_donor)
      
      # determine whether subsampling was previously recorded
      if (length(prev_subs)>1) {
        for (j in 1:length(prev_subs)) {
          if (identical(prev_subs[[j]],cells_sampled)) {
            prev_ident <- TRUE
          }
        }
      }
    }
    prev_subs[[length(prev_subs)+1]] <- cells_sampled
    
    counts_sub <- combined_counts_sparse[,cells_sampled]
    
    # prep data
    param_list <- initialize_params(ctypes_use = c("ct1", "ct2"),
                                    ncores = 30, rand_seed = 10)
    
    pbmc_container <- make_new_container(count_data=counts_sub, meta_data=meta_sub,
                                         params=param_list)
    
    pbmc_container <- form_tensor(pbmc_container, donor_min_cells=0, 
                                  norm_method='trim', scale_factor=10000, 
                                  vargenes_method='norm_var_pvals', vargenes_thresh=1,
                                  scale_var = TRUE, var_scale_power = .5)
    
    # do decomposition
    pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(2,4),
                                     tucker_type = 'regular', rotation_type = 'ica_dsc')

    # determine which factor is process1
    dscores <- pbmc_container$tucker_results[[1]]
    m1 <- abs(mean(dscores[dg1,1]))
    m2 <- abs(mean(dscores[dg1,2]))
    ndx_max <- order(c(m1,m2),decreasing=TRUE)[1]
    ndx_other <- order(c(m1,m2),decreasing=TRUE)[2]
    
    # get significant genes
    # pbmc_container <- run_jackstraw(pbmc_container, ranks=c(2,4,2), n_fibers=100, n_iter=1000,
    #                                 tucker_type='regular', rotation_type='ica')
    pbmc_container <- get_lm_pvals(pbmc_container)
    
    # evaluate AUC
    # do factor 1 first
    sig_vectors <- get_significance_vectors(pbmc_container,
                                            factor_select=ndx_max, c('ct1','ct2'))
    # convert list to df
    sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
    pred1 <- c(sig_df)
    
    # ## trying using loadings instead
    # sig_df <- get_one_factor(pbmc_container,ndx_max)[[2]]
    # pred1 <- abs(c(sig_df))
    # ##
    
    ct1_de <- de1_ct1
    ct1_de <- ct1_de[rownames(sig_df),]
    ct2_de <- de1_ct2
    ct2_de <- ct2_de[rownames(sig_df),]
    de_true1 <- c(ct1_de$DEFacGroup2!=1,ct2_de$DEFacGroup2!=1)
    
    pROC_obj <- roc(de_true1,pred1,
                    smoothed = TRUE,
                    plot=FALSE, AUC=TRUE)
    auc1 <- pROC_obj[["auc"]]
    auc1
    ## now do factor 2
    sig_vectors <- get_significance_vectors(pbmc_container,
                                            factor_select=ndx_other, c('ct1','ct2'))
    # convert list to df
    sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
    pred2 <- c(sig_df)
    
    # ## trying using loadings instead
    # sig_df <- get_one_factor(pbmc_container,ndx_other)[[2]]
    # pred2 <- abs(c(sig_df))
    # ##
    
    ct1_de <- de2_ct1
    ct1_de <- ct1_de[rownames(sig_df),]
    ct2_de <- de2_ct2
    ct2_de <- ct2_de[rownames(sig_df),]
    de_true2 <- c(ct1_de$DEFacGroup2!=1,ct2_de$DEFacGroup2!=1)
    
    pROC_obj <- roc(de_true2,pred2,
                    smoothed = TRUE,
                    plot=FALSE, AUC=TRUE)
    auc2 <- pROC_obj[["auc"]]
    
    # store results
    if (ndx_max==1) {
      f1_aucs <- c(f1_aucs,auc1)
      f2_aucs <- c(f2_aucs,auc2)
    } else {
      f1_aucs <- c(f1_aucs,auc2)
      f2_aucs <- c(f2_aucs,auc1)
    }
    
    cpd <- c(cpd,ds)
  }
}


# plot results
tmp <- cbind.data.frame(c(f1_aucs,f2_aucs),c(cpd,cpd),c(rep('f1',length(f1_aucs)),rep('f2',length(f2_aucs))))
colnames(tmp) <- c('auc_val','cells_total','process')

tmp$cells_per <- sapply(tmp$cells_total,function(x) {
  return(x/(2*num_donors))
})

# get mean values for each factor at each value of cells_per
tmp2_means <- c()
tmp2_process <- c()
tmp2_cells_per <- c()
for (cp in unique(tmp$cells_per)) {
  tmp_sub <- tmp[tmp$cells_per==cp,]
  
  tmp_sub_sub <- tmp_sub[tmp_sub$process=='f1',]
  tmp2_means <- c(tmp2_means,mean(tmp_sub_sub$auc_val))
  tmp2_process <- c(tmp2_process,'f1')
  tmp2_cells_per <- c(tmp2_cells_per,cp)
  
  tmp_sub_sub <- tmp_sub[tmp_sub$process=='f2',]
  tmp2_means <- c(tmp2_means,mean(tmp_sub_sub$auc_val))
  tmp2_process <- c(tmp2_process,'f2')
  tmp2_cells_per <- c(tmp2_cells_per,cp)
}

tmp2 <- cbind.data.frame(tmp2_means,tmp2_process,tmp2_cells_per)
colnames(tmp2) <- c('auc_mean','pr','mcp')

pdf(file = "/home/jmitchel/figures/for_paper_v2/sim_AUC_vs_ncells.pdf", useDingbats = FALSE,
    width = 5, height = 5)
ggplot(tmp,aes(x=cells_per,y=auc_val,color=process)) +
  geom_point() +
  geom_line(data=tmp2,aes(x=mcp,y=auc_mean,color=pr)) +
  xlab('Av. cells per donor_cell type') +
  ylab('AUC') +
  theme_bw() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))
dev.off()


# ## save results so far
# save.image(file='/home/jmitchel/data/sim_data/sim_v3.RData')
load(file='/home/jmitchel/data/sim_data/sim_v3.RData')




#### comparing use of jackstraw to no jackstraw in controlling FDR
library(pROC)
# should see around 5% of p<0.05 as incorrect with jackstraw
# should see higher than 5% of p<0.05 as incorrect without jackstraw

# do factor 1 first
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=ndx_max, c('ct1','ct2'))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
pred1 <- c(sig_df)

ct1_de <- de1_ct1
ct1_de <- ct1_de[rownames(sig_df),]
ct2_de <- de1_ct2
ct2_de <- ct2_de[rownames(sig_df),]
de_true1 <- c(ct1_de$DEFacGroup2!=1,ct2_de$DEFacGroup2!=1)

ndx_pos <- which(pred1<0.05)
rv <- de_true1[ndx_pos]
tp <- sum(rv)
fp <- sum(!rv)
fp/length(rv) # 2% of padj<0.05 were false

# get associations as gene.ct.factor
pvals <- get_real_fstats(pbmc_container,ncores=4)
padj <- p.adjust(pvals,method='fdr')
names(padj) <- sapply(names(padj),function(x) {
  return(substr(x,1,nchar(x)-6))
})

pbmc_container[["gene_score_associations"]] <- padj
# recalculate the above...




