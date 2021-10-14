library("splatter")
library("scater")

vcf <- mockVCF(n.samples = 200)
gff <- mockGFF(n.genes = 2000)

# now trying for multi population two condition sim and adding more donors
params.cond <- newSplatPopParams(similarity.scale = 1.5,
                                 condition.prob = c(0.1, 0.9),
                                 cde.prob = c(0,.1),
                                 cde.downProb = c(0,.5),
                                 cde.facLoc = 1.05, 
                                 cde.facScale = 0.5,
                                 de.prob = 0.3,
                                 de.facLoc = 0.5, 
                                 de.facScale = 0.5,
                                 group.prob = c(0.5, 0.5),
                                 batchCells = c(10),
                                 batch.size = 200, nGenes=2000,
                                 eqtl.n = 0,
                                 seed = 98765)

sim.pop.cond1 <- splatPopSimulate(vcf = vcf, gff = gff,
                                 params = params.cond, 
                                 sparsify = FALSE)

# need to make cell names unique per donor and group
colnames(sim.pop.cond1) <- paste(sim.pop.cond1$Sample, sim.pop.cond1$Group, sim.pop.cond1$Cell, sep="_")

sim.pop.cond1 <- logNormCounts(sim.pop.cond1)
sim.pop.cond1 <- runPCA(sim.pop.cond1, ncomponents = 10)
plotPCA(sim.pop.cond1, colour_by = "Sample", shape_by = "Condition", point_size=3)

simcounts1 <- counts(sim.pop.cond1)
simcounts1 <- methods::as(as.matrix(simcounts1),'sparseMatrix')
simmeta1 <- as.data.frame(colData(sim.pop.cond1)@listData)
colnames(simmeta1)[4] <- 'donors'
colnames(simmeta1)[5] <- 'ctypes'

param_list <- initialize_params(ctypes_use = c("Group1", "Group2"),
                                ncores = 30, rand_seed = 10)

sim_container <- make_new_container(count_data=simcounts1, meta_data=simmeta1, params=param_list)

sim_container <- form_tensor(sim_container, donor_min_cells=3,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var', vargenes_thresh=1900,
                              scale_var = TRUE, var_scale_power = .5)

sim_container <- run_tucker_ica(sim_container, ranks=c(2,4),
                                 tucker_type = 'regular', rotation_type = 'ica_dsc')
sim_container <- run_tucker_ica(sim_container, ranks=c(2,4),
                                 tucker_type = 'regular', rotation_type = 'hybrid')
pca_unfolded(sim_container,2)


sim_container <- plot_donor_matrix(sim_container, 
                                    show_donor_ids = TRUE)

sim_container$plots$donor_matrix

sim_container <- get_lm_pvals(sim_container)

sim_container <- plot_loadings_annot(sim_container, factor_select=1, use_sig_only=FALSE, nonsig_to_zero=FALSE, annot='sig_genes',
                                      pathways=NULL, sig_thresh=0.02, display_genes=FALSE, 
                                      gene_callouts=FALSE, show_xlab=TRUE,
                                      show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)



## trying to add in another process by averaging normalized counts
params.cond <- newSplatPopParams(similarity.scale = 1.5,
                                 condition.prob = c(0.1, 0.9),
                                 cde.prob = c(0,.1),
                                 cde.downProb = c(0,.5),
                                 cde.facLoc = 1.05, 
                                 cde.facScale = 0.5,
                                 de.prob = 0.3,
                                 de.facLoc = 0.5, 
                                 de.facScale = 0.5,
                                 group.prob = c(0.5, 0.5),
                                 batchCells = c(10),
                                 batch.size = 200, nGenes=2000,
                                 eqtl.n = 0,
                                 seed = 12345)
sim.pop.cond2 <- splatPopSimulate(vcf = vcf, gff = gff,
                                 params = params.cond, 
                                 sparsify = FALSE)

# need to make cell names unique per donor and group
colnames(sim.pop.cond2) <- paste(sim.pop.cond2$Sample, sim.pop.cond2$Group, sim.pop.cond2$Cell, sep="_")

sim.pop.cond2 <- logNormCounts(sim.pop.cond2)
sim.pop.cond2 <- runPCA(sim.pop.cond2, ncomponents = 10)
plotPCA(sim.pop.cond2, colour_by = "Sample", shape_by = "Condition", point_size=3)


# simmeta2 <- as.data.frame(colData(sim.pop.cond2)@listData)
# simmeta2[1:50,]
# # appears that sample_x not in same condition groups as before
# dim(sim.pop.cond2)

simcounts2 <- counts(sim.pop.cond2)
simcounts2 <- methods::as(as.matrix(simcounts2),'sparseMatrix')
simmeta2 <- as.data.frame(colData(sim.pop.cond2)@listData)
colnames(simmeta2)[4] <- 'donors'
colnames(simmeta2)[5] <- 'ctypes'

# # seeing how the internal logcounts fn works
# test <- normalize_counts(simcounts2,mean(colSums(simcounts2))) # sce uses log base 2 instead of log base e
# test[1:10,1:10]

# # looking at process overlap between the two simulations
# simmeta1_2 <- unique(simmeta1[,c(4,6)])
# simmeta2_2 <- unique(simmeta2[,c(4,6)])

## loop through donors, match up cells, average log counts, revert to counts by first lib size
d_list <- unique(simmeta1$donors)
ctypes <- unique(simmeta1$ctypes)
new_log_norm <- matrix(ncol=ncol(logcounts(sim.pop.cond1)), nrow=nrow(logcounts(sim.pop.cond1)))
colnames(new_log_norm) <- colnames(logcounts(sim.pop.cond1))
rownames(new_log_norm) <- rownames(logcounts(sim.pop.cond1))
for (dnr in d_list) {
  for (ct in ctypes) {
    # get donor-ctype counts for first sim
    meta_sub1 <- simmeta1[simmeta1$donors==dnr,]
    d_cells1 <- rownames(meta_sub1)[meta_sub1$ctypes==ct]
    d_counts1 <- logcounts(sim.pop.cond1)[,d_cells1]
    
    # get donor-ctype counts for second sim
    meta_sub2 <- simmeta2[simmeta2$donors==dnr,]
    d_cells2 <- rownames(meta_sub2)[meta_sub2$ctypes==ct]
    d_counts2 <- logcounts(sim.pop.cond2)[,d_cells2]
    
    # take pairwise cell averages from the two simulations
    d_counts_av <- (d_counts1 + d_counts2) / 2
    
    # store results
    new_log_norm[rownames(d_counts_av),colnames(d_counts_av)] <- d_counts_av
  }
}

# check that there are no more NA values in the averaged results
sum(is.na(new_log_norm))

# exponentiate the data
new_log_norm2 <- 2**new_log_norm

# subtract 1
new_log_norm2 <- new_log_norm2 - 1

# divide by scale factor
lib_sizes <- colSums(simcounts1)
scale_factor <- mean(lib_sizes)
new_log_norm2 <- new_log_norm2 / scale_factor

# adjust fractional counts to sum to 1
adjustments <- 1/colSums(new_log_norm2)
new_log_norm2 <- sweep(new_log_norm2,MARGIN=2,adjustments,FUN='*')

# multiply to get correct original sizes
new_log_norm2 <- sweep(new_log_norm2,MARGIN=2,lib_sizes,FUN='*')
## these are honestly too large... should add sparsity or adjust the param

new_simcounts <- round(new_log_norm2)
new_simcounts <- methods::as(as.matrix(new_simcounts),'sparseMatrix')



# new_simcounts_store <- new_simcounts
# simmeta1_store <- simmeta1
# new_simcounts <- new_simcounts[,cells_sub]
# simmeta1 <- simmeta1[cells_sub,]



## should now get out two separate processes
param_list <- initialize_params(ctypes_use = c("Group1", "Group2"),
                                ncores = 30, rand_seed = 10)

sim_container <- make_new_container(count_data=new_simcounts, meta_data=simmeta1, params=param_list)

sim_container <- form_tensor(sim_container, donor_min_cells=3,
                             norm_method='trim', scale_factor=10000,
                             vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                             scale_var = TRUE, var_scale_power = .5)

sim_container <- run_tucker_ica(sim_container, ranks=c(2,4),
                                tucker_type = 'regular', rotation_type = 'ica_dsc')
sim_container <- run_tucker_ica(sim_container, ranks=c(2,4),
                                tucker_type = 'regular', rotation_type = 'hybrid')
sim_container <- run_tucker_ica(sim_container, ranks=c(5,12),
                                tucker_type = 'regular', rotation_type = 'hybrid')

pca_unfolded(sim_container,2)
nmf_unfolded(sim_container,2)

sim_container <- plot_donor_matrix(sim_container, 
                                   show_donor_ids = TRUE)

sim_container$plots$donor_matrix

sim_container <- get_lm_pvals(sim_container)

sim_container <- plot_loadings_annot(sim_container, factor_select=1, use_sig_only=FALSE, nonsig_to_zero=FALSE, annot='sig_genes',
                                     pathways=NULL, sig_thresh=0.02, display_genes=FALSE, 
                                     gene_callouts=FALSE, show_xlab=TRUE,
                                     show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)

cor(t(sim_container[["tucker_results"]][[2]]))

# rd1 <- as.data.frame(rowData(sim.pop.cond1)@listData)
# rownames(rd1) <- rd1$Gene
# 
# rd2 <- as.data.frame(rowData(sim.pop.cond2)@listData)
# rownames(rd2) <- rd2$Gene
# 
# sim_container <- plot_loadings_annot(sim_container, factor_select=1, use_sig_only=FALSE, nonsig_to_zero=FALSE, annot='sig_genes',
#                                      pathways=NULL, sig_thresh=0.02, display_genes=FALSE, 
#                                      gene_callouts=FALSE, show_xlab=TRUE, sim_de_donor_group=list(rd1,rd2),
#                                      show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)

# new_simcounts <- new_simcounts_store
# cells_sub <- sample(colnames(new_simcounts),3000)
# 
# simmeta1_sub <- simmeta1[cells_sub,]
# simmeta2_sub <- simmeta2[cells_sub,]
# coldat_lst <- list(simmeta1_sub,simmeta2_sub)


container <- sim_container
coldat_lst <- list(simmeta1,simmeta2)
rowdat_lst <- list(rowData(sim.pop.cond1),rowData(sim.pop.cond2))
get_sim_auc(sim_container,coldat_lst,rowdat_lst)

get_sim_auc <- function(container, coldat_lst, rowdat_lst) {
  all_auc <- c()
  all_pv <- c()
  for (i in 1:length(coldat_lst)) {
    # determine donors of condition 1 for the given process
    coldat <- coldat_lst[[i]]
    dg1 <- unique(coldat$donors[coldat$Condition=='Condition1'])
    
    # determine which factor best distinguishes these donors
    dscores <- container$tucker_results[[1]]
    fact_links <- abs(colMeans(dscores[dg1,]))
    ndx_max <- order(fact_links,decreasing=TRUE)[1]
    
    # evaluate AUC
    # # do factor 1 first
    # sig_vectors <- get_significance_vectors(container,
    #                                         factor_select=ndx_max, container$experiment_params$ctypes_use)
    # # convert list to df
    # sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
    # pred1 <- c(sig_df)
    
    ## trying using loadings instead
    f_res <- get_one_factor(container,ndx_max)
    sig_df <- f_res[[2]]
    pred1 <- abs(c(sig_df))
    ##
    
    # should make this work with more than two cell types
    rowdat <- rowdat_lst[[i]]
    ct1_de <- rowdat[rownames(sig_df),c('ConditionDE.Condition1','ConditionDE.Condition2')]
    ct1_de <- as.matrix(ct1_de)
    ct1_de_logical <- rowSums(ct1_de!=1)>0

    ## remove exceptions here
    
    ##
    
    # de_true1 <- c(ct1_de_logical1,ct2_de_logical)
    de_true1 <- rep(ct1_de_logical,ncol(sig_df)) # repeated for each ctype
    
    pROC_obj <- roc(de_true1,pred1,
                    smoothed = TRUE,
                    plot=FALSE, AUC=TRUE)
    auc1 <- pROC_obj[["auc"]]
    # print(auc1)
    
    # store AUC result
    all_auc <- c(all_auc,auc1)
    
    ## now compute donor process association p-values
    coldat <- unique(coldat[,c('donors','Condition')])
    rownames(coldat) <- coldat$donors
    dsc <- f_res[[1]]
    tmp <- cbind.data.frame(dsc,coldat[rownames(dsc),'Condition'])
    colnames(tmp) <- c('dscore','Condition')
    lmres <- summary(lm(as.formula('dscore~Condition'),data=tmp))
    pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
    all_pv <- c(all_pv,-log10(pval))
  }
  return(list(all_auc,all_pv))
}



generate_sim <- function(n_ctypes,ctype_de_prob,n_donors,cells_per_don,
                         n_genes,n_process,de_strength,condition.prob,
                         similarity.scale, sim_seed=1) {
  coldat_lst <- list()
  rowdat_lst <- list()
  sim_objs <- list()
  
  vcf <- mockVCF(n.samples = n_donors)
  gff <- mockGFF(n.genes = n_genes)
  
  for (i in 1:n_process) {
    seed <- i * sim_seed
    
    # make sure group.prob sums to 1
    gpr <- rep(1/n_ctypes, n_ctypes)
    gpr <- round(gpr,2)
    gpr[length(gpr)] <- 1 - sum(gpr[1:(length(gpr)-1)])
    
    # now trying for multi population two condition sim and adding more donors
    params.cond <- newSplatPopParams(similarity.scale = similarity.scale,
                                     condition.prob = condition.prob,
                                     cde.prob = c(0,.1),
                                     cde.downProb = c(0,.5),
                                     cde.facLoc = de_strength, 
                                     cde.facScale = 0.5,
                                     de.prob = ctype_de_prob,
                                     de.facLoc = 0.5, 
                                     de.facScale = 0.5,
                                     group.prob = gpr,
                                     batchCells = c(cells_per_don),
                                     batch.size = n_donors, nGenes=n_genes,
                                     eqtl.n = 0,
                                     seed = seed)
    
    sim.pop.cond1 <- splatPopSimulate(vcf = vcf, gff = gff,
                                      params = params.cond, 
                                      sparsify = FALSE)
    
    # need to make cell names unique per donor and group
    colnames(sim.pop.cond1) <- paste(sim.pop.cond1$Sample, sim.pop.cond1$Group, sim.pop.cond1$Cell, sep="_")
    
    sim.pop.cond1 <- logNormCounts(sim.pop.cond1)
    simcounts1 <- counts(sim.pop.cond1)
    simcounts1 <- methods::as(as.matrix(simcounts1),'sparseMatrix')
    simmeta1 <- as.data.frame(colData(sim.pop.cond1)@listData)
    colnames(simmeta1)[4] <- 'donors'
    colnames(simmeta1)[5] <- 'ctypes'
    
    coldat_lst[[i]] <- simmeta1
    rowdat_lst[[i]] <- rowData(sim.pop.cond1)
    sim_objs[[i]] <- sim.pop.cond1
  }
  
  new_simcounts <- combine_processes(sim_objs,coldat_lst,rowdat_lst)
  
  return(list(new_simcounts,coldat_lst,rowdat_lst,sim_objs))
}


combine_processes <- function(sim_objs,coldat_lst,rowdat_lst) {
  simmeta1 <- coldat_lst[[1]]
  sim.pop.cond1 <- sim_objs[[1]]
  d_list <- unique(simmeta1$donors)
  ctypes <- unique(simmeta1$ctypes)
  
  # matrix to store new combined data
  new_log_norm <- matrix(ncol=ncol(logcounts(sim.pop.cond1)), nrow=nrow(logcounts(sim.pop.cond1)))
  colnames(new_log_norm) <- colnames(logcounts(sim.pop.cond1))
  rownames(new_log_norm) <- rownames(logcounts(sim.pop.cond1))
  for (dnr in d_list) {
    for (ct in ctypes) {
      d_counts_list <- list()
      for (pr_ndx in 1:length(sim_objs)) {
        simmeta1 <- coldat_lst[[pr_ndx]]
        sim.pop.cond1 <- sim_objs[[pr_ndx]]
        meta_sub1 <- simmeta1[simmeta1$donors==dnr,]
        d_cells1 <- rownames(meta_sub1)[meta_sub1$ctypes==ct]
        d_counts1 <- logcounts(sim.pop.cond1)[,d_cells1]
        
        d_counts_list[[pr_ndx]] <- d_counts1
      }
      
      # compute average expression
      running_total <- d_counts_list[[1]]
      for (pr_ndx in 2:length(d_counts_list)) {
        running_total <- running_total + d_counts_list[[pr_ndx]]
      }
      d_counts_av <- running_total / length(d_counts_list)
      
      # store results
      new_log_norm[rownames(d_counts_av),colnames(d_counts_av)] <- d_counts_av
    }
  }
  
  # exponentiate the data
  new_log_norm2 <- 2**new_log_norm
  
  # subtract 1
  new_log_norm2 <- new_log_norm2 - 1
  
  # divide by scale factor
  lib_sizes <- colSums(counts(sim_objs[[1]])) # arbitrarily using the first sim for lib sizes
  scale_factor <- mean(lib_sizes)
  new_log_norm2 <- new_log_norm2 / scale_factor
  
  # adjust fractional counts to sum to 1
  adjustments <- 1/colSums(new_log_norm2)
  new_log_norm2 <- sweep(new_log_norm2,MARGIN=2,adjustments,FUN='*')
  
  # multiply to get correct original sizes
  new_log_norm2 <- sweep(new_log_norm2,MARGIN=2,lib_sizes,FUN='*')
  ## these are honestly too large... should add sparsity or adjust the param
  
  new_simcounts <- round(new_log_norm2)
  new_simcounts <- methods::as(as.matrix(new_simcounts),'sparseMatrix')
  
  return(new_simcounts)
}



sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=15,n_genes=2000,n_process=3,
                        de_strength=1.05,condition.prob=c(.25,.75),
                        similarity.scale=1.5) # works decently well

# trying to increase cells_per_donor
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=3,
                        de_strength=1.05,condition.prob=c(.25,.75),
                        similarity.scale=1.5) # gives roughly same results as above, maybe slightly better


# trying to increase de_strength
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=3,
                        de_strength=1.5,condition.prob=c(.25,.75),
                        similarity.scale=1.5) # this makes huge difference and Tucker much better than pca

# seeing what reducing de strength to less than 1 does 
# if it's essentially a mean FC, then a FC of 1.5 essentially equals FC of 1/1.5=.66 
# alternatively if it's a logFC, then FC of .66 is jsut less strong...
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=3,
                        de_strength=.666,condition.prob=c(.25,.75),
                        similarity.scale=1.5) # performance is much worse overall so it must be log2FC

# trying much higher de_strength
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=3,
                        de_strength=2.5,condition.prob=c(.25,.75),
                        similarity.scale=1.5) # performance better still, with Tucker doing best

# trying higher similarity.scale with middle range de_strength 
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=3,
                        de_strength=1.5,condition.prob=c(.25,.75),
                        similarity.scale=3) # performance increases overall, Tucker still better

# trying original similarity.scale, middle de_strength, and 50/50 condition.prob
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=3,
                        de_strength=1.5,condition.prob=c(.5,.5),
                        similarity.scale=1.5) # performance increases slightly, Tucker still better

# trying lower ctype_de_prob (more similar cells), middle de_strength, original everything else
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.1,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=3,
                        de_strength=1.5,condition.prob=c(.25,.75),
                        similarity.scale=1.5, sim_seed=120) # this doesn't really affect performance, as expected

# trying larger number of donors, middle de_strength, original everything else
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=200,
                        cells_per_don=50,n_genes=2000,n_process=3,
                        de_strength=1.5,condition.prob=c(.25,.75),
                        similarity.scale=1.5) # better overall performance, Tucker still far better


# going back to middle de_strength, original everything else - will compare to reducing number of processes
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=2,
                        de_strength=1.5,condition.prob=c(.25,.75),
                        similarity.scale=1.5) # with fewer processes performance is closer between Tucker/PCA, but Tucker still better

# trying same params as above but with lower de_strength
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=2,
                        de_strength=1.05,condition.prob=c(.25,.75),
                        similarity.scale=1.5) # with fewer processes + lower de_strength, performance basically exact same between Tucker/PCA

# trying same params as above but with equal condition.prob
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=2,
                        de_strength=1.05,condition.prob=c(.5,.5),
                        similarity.scale=1.5) # Tucker does much better than PCA here, strangely

# trying to change sim seed
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=2,
                        de_strength=1.5,condition.prob=c(.5,.5),
                        similarity.scale=1.5, sim_seed=12) # yields very similar results to one above

# trying to increase number of processes even more
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=4,
                        de_strength=2,condition.prob=c(.25,.75),
                        similarity.scale=1.5, sim_seed=12) # it does seem that more processes equals greater separation between the two methods

# keeping params same as above but increasing n_ctypes
sim_res <- generate_sim(n_ctypes=4,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=4,
                        de_strength=2,condition.prob=c(.25,.75),
                        similarity.scale=1.5, sim_seed=12) #

# trying to recapitulate result from below
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=3,
                        de_strength=1,condition.prob=c(.25,.75),
                        similarity.scale=1.5, sim_seed=193) #

# testing if n_genes makes a big difference: lower end
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=500,n_process=3,
                        de_strength=1.5,condition.prob=c(.25,.75),
                        similarity.scale=1.5, sim_seed=193) # around .05 diff

# testing if n_genes makes a big difference: higher end
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                        cells_per_don=50,n_genes=2000,n_process=3,
                        de_strength=1.5,condition.prob=c(.25,.75),
                        similarity.scale=1.5, sim_seed=193) # doesn't seem to be huge difference in the gap

# sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.3,n_donors=50,
#                         cells_per_don=20,n_genes=2000,n_process=4,
#                         de_strength=1.05,condition.prob=c(.25,.75),
#                         similarity.scale=1.5) # also works well!
# 
# sim_res <- generate_sim(n_ctypes=3,ctype_de_prob=.3,n_donors=50,
#                         cells_per_don=30,n_genes=2000,n_process=3,
#                         de_strength=1.05,condition.prob=c(.25,.75),
#                         similarity.scale=1.5) # works very well!

new_simcounts <- sim_res[[1]]
coldat_lst <- sim_res[[2]]
rowdat_lst <- sim_res[[3]]
sim_objs <- sim_res[[4]]


param_list <- initialize_params(ctypes_use = c("Group1", "Group2", "Group3", "Group4"),
                                ncores = 30, rand_seed = 10)

param_list <- initialize_params(ctypes_use = c("Group1", "Group2"),
                                ncores = 30, rand_seed = 10)

sim_container <- make_new_container(count_data=new_simcounts, meta_data=coldat_lst[[1]], params=param_list)

sim_container <- form_tensor(sim_container, donor_min_cells=3,
                             norm_method='trim', scale_factor=10000,
                             vargenes_method='norm_var_pvals', vargenes_thresh=1,
                             scale_var = TRUE, var_scale_power = .5)

# sim_container <- run_tucker_ica(sim_container, ranks=c(3,6),
#                                 tucker_type = 'regular', rotation_type = 'ica_dsc')
sim_container <- run_tucker_ica(sim_container, ranks=c(2,4),
                                tucker_type = 'regular', rotation_type = 'hybrid')
sim_container <- run_tucker_ica(sim_container, ranks=c(3,6),
                                tucker_type = 'regular', rotation_type = 'hybrid')

pca_unfolded(sim_container,3)
nmf_unfolded(sim_container,2)

sim_container <- plot_donor_matrix(sim_container, 
                                   show_donor_ids = TRUE)

sim_container$plots$donor_matrix

sim_container <- get_lm_pvals(sim_container)

get_sim_auc(sim_container,coldat_lst,rowdat_lst)


## now for systematic tests and plotting for paper...
# compute AUC values for each process as well as max dscore associations for each process
# will do comparisons of my method against both PCA and NMF on the unfolded tensor
# as a note: NNF only works with the gene significance AUC values because it won't correctly capture the loadings for negative DE genes!
# when showing factor associations NMF will do pretty well (but not as good as tucker), and when showing the loadings AUCs it will just do terribly.
##




get_performance_range <- function(param_test, param_range, n_iter, rand_seeds,
                                   n_ctypes=2,ctype_de_prob=.3,n_donors=50,
                                   cells_per_don=50,n_genes=2000,n_process=2,
                                   de_strength=1,condition.prob=c(.25,.75),
                                   similarity.scale=1.5) {
  
  # matrix to store results
  myres <- matrix(ncol=4,nrow=0)
  colnames(myres) <- c('mean_auc','mean_pv','method','param_val')
  myres_summary <- matrix(ncol=6,nrow=0)
  colnames(myres_summary) <- c('mean_auc','sd_auc','mean_pv','sd_pv','method','param_val')
  
  for (i in 1:length(param_range)) {
    
    # need to switch parameter value for param_test
    if (param_test=='n_ctypes') {
      n_ctypes <- param_range[i]
    } else if (param_test=='ctype_de_prob') {
      ctype_de_prob <- param_range[i]
    } else if (param_test=='n_donors') {
      n_donors <- param_range[i]
    } else if (param_test=='cells_per_don') {
      cells_per_don <- param_range[i]
    } else if (param_test=='n_genes') {
      n_genes <- param_range[i]
    } else if (param_test=='n_process') {
      n_process <- param_range[i]
    } else if (param_test=='de_strength') {
      de_strength <- param_range[i]
    } else if (param_test=='similarity.scale') {
      similarity.scale <- param_range[i]
    } else {
      stop('choose one of the available parameter options')
    }
    
    for (j in 1:n_iter) {
      print(i)
      print(j)
      sim_res <- generate_sim(n_ctypes=n_ctypes,ctype_de_prob=ctype_de_prob,
                              n_donors=n_donors,
                              cells_per_don=cells_per_don,n_genes=n_genes,n_process=n_process,
                              de_strength=de_strength,condition.prob=condition.prob,
                              similarity.scale=similarity.scale, sim_seed=(rand_seeds[i]+j))
      
      new_simcounts <- sim_res[[1]]
      coldat_lst <- sim_res[[2]]
      rowdat_lst <- sim_res[[3]]
      sim_objs <- sim_res[[4]]
      
      ctypes_use <- sapply(1:n_ctypes,function(x){
        paste0('Group',x)
      })
      
      param_list <- initialize_params(ctypes_use = ctypes_use,
                                      ncores = 30, rand_seed = 10)
      
      sim_container <- make_new_container(count_data=new_simcounts, meta_data=coldat_lst[[1]], params=param_list)
      
      sim_container <- form_tensor(sim_container, donor_min_cells=3,
                                   norm_method='trim', scale_factor=10000,
                                   vargenes_method='norm_var_pvals', vargenes_thresh=1,
                                   scale_var = TRUE, var_scale_power = .5)

      # get AUC and pvals with Tucker
      sim_container <- run_tucker_ica(sim_container, ranks=c(n_process,n_process*2),
                                      tucker_type = 'regular', rotation_type = 'hybrid')
      
      all_tucker_res <- get_sim_auc(sim_container,coldat_lst,rowdat_lst)
      tucker_auc <- mean(all_tucker_res[[1]])
      tucker_pv <- mean(all_tucker_res[[2]])
      row_add <- data.frame('mean_auc'=tucker_auc,'mean_pv'=tucker_pv,'method'='Tucker','param_val'=param_range[i])
      myres <- rbind.data.frame(myres,row_add) # store results
      
      # get AUC and pvals with PCA
      pca_unfolded(sim_container,n_process)
      all_pca_res <- get_sim_auc(sim_container,coldat_lst,rowdat_lst)
      pca_auc <- mean(all_pca_res[[1]])
      pca_pv <- mean(all_pca_res[[2]])
      row_add <- data.frame('mean_auc'=pca_auc,'mean_pv'=pca_pv,'method'='PCA','param_val'=param_range[i])
      myres <- rbind.data.frame(myres,row_add) # store results
      
      # get AUC and pvals with NMF
      nmf_unfolded(sim_container,n_process)
      all_nmf_res <- get_sim_auc(sim_container,coldat_lst,rowdat_lst)
      nmf_auc <- mean(all_nmf_res[[1]])
      nmf_pv <- mean(all_nmf_res[[2]])
      row_add <- data.frame('mean_auc'=nmf_auc,'mean_pv'=nmf_pv,'method'='NMF','param_val'=param_range[i])
      myres <- rbind.data.frame(myres,row_add) # store results
    }
    myres_sub <- myres[myres$param_val==param_range[i],]
    # auc_summary = tapply(myres_sub$mean_auc, myres_sub$method, mean) 
    # pv_summary = tapply(myres_sub$mean_pv, myres_sub$method, mean) 
    # pv_summary <- pv_summary[names(auc_summary)] # ordering them same way
    
    auc_mean = tapply(myres_sub$mean_auc, myres_sub$method, mean) 
    pv_mean = tapply(myres_sub$mean_pv, myres_sub$method, mean) 
    pv_mean <- pv_mean[names(auc_mean)] # ordering them same way
    
    auc_sd = tapply(myres_sub$mean_auc, myres_sub$method, sd) 
    pv_sd = tapply(myres_sub$mean_pv, myres_sub$method, sd) 
    pv_sd <- pv_sd[names(auc_sd)] # ordering them same way
    
    row_add <- data.frame('mean_auc'=auc_mean,'sd_auc'=auc_sd,'mean_pv'=pv_mean,'sd_pv'=pv_sd,'method'=names(pv_mean),'param_val'=param_range[i])
    myres_summary <- rbind.data.frame(myres_summary,row_add)
  }
  return(list(myres,myres_summary))
}

all_res <- get_performance_range(param_test='de_strength', param_range=c(1,1.5), n_iter=2, 
                      rand_seeds=c(192,511), n_ctypes=2, ctype_de_prob=.3, n_donors=50,
                      cells_per_don=20, n_genes=500, n_process=2,
                      de_strength=1, condition.prob=c(.25,.75),
                      similarity.scale=1.5)

all_res_saved <- all_res
all_res_saved2 <- all_res

r1 <- rbind.data.frame(all_res_saved[[1]],all_res_saved2[[1]])
r2 <- rbind.data.frame(all_res_saved[[2]],all_res_saved2[[2]])

r1 <- all_res[[1]]
r2 <- all_res[[2]]
p1 <- ggplot(r1,aes(x=param_val,y=mean_auc,color=method)) +
  geom_point() +
  geom_line(data=r2,aes(x=param_val,y=mean_auc,color=method)) +
  geom_errorbar(data=r2, aes(x=param_val, ymin=mean_auc-sd_auc, ymax=mean_auc+sd_auc, color=method),
                width = .05) +
  xlab('Process DE strength') +
  ylab('Mean AUC') +
  theme_bw()
p2 <- ggplot(r1,aes(x=param_val,y=mean_pv,color=method)) +
  geom_point() +
  geom_line(data=r2,aes(x=param_val,y=mean_pv,color=method)) +
  geom_errorbar(data=r2, aes(x=param_val, ymin=mean_pv-sd_pv, ymax=mean_pv+sd_pv, color=method),
                width = .05) +
  xlab('Process DE strength') +
  ylab('Mean -log10(pval)') +
  theme_bw()
combined_plt <- cowplot::plot_grid(p1,p2,nrow=2,align = 'v')










## post analysis of running the different iter sim tests
library(cowplot)

# n_proc data
all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_proc.rds')
all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_proc_2.rds')
all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_proc_3.rds')
all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_proc_4.rds')

# other
all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_donors.rds')
all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_donors_2.rds')
# all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_cells_per_don.rds')
# all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_cells_per_don_2.rds')
all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_de_strength.rds')
all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_de_strength_2.rds')
all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_ctypes.rds')
all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_ctypes_2.rds')
all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_similarity_scale.rds')
all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_similarity_scale_2.rds')

rname <- 'Number of processes'
rname <- 'Number of donors'
# rname <- 'Number of cells/donor/cell-type'
rname <- 'DE strength'
rname <- 'Number of cell types'
rname <- 'Donor-donor similarity'

r1 <- all_res[[1]]
r2 <- all_res[[2]]

# to use SEM instead of SD
r2[,'sd_auc'] <- r2[,'sd_auc'] / sqrt(10)
r2[,'sd_pv'] <- r2[,'sd_pv'] / sqrt(10)

p1 <- ggplot(r1,aes(x=param_val,y=mean_auc,color=method)) +
  # geom_point() +
  geom_line(data=r2,aes(x=param_val,y=mean_auc,color=method)) +
  geom_errorbar(data=r2, aes(x=param_val, ymin=mean_auc-sd_auc, ymax=mean_auc+sd_auc, color=method),
                width = .1) +
  xlab(rname) +
  ylab('Mean AUC') +
  theme_bw() +
  theme(legend.position="none")
p2 <- ggplot(r1,aes(x=param_val,y=mean_pv,color=method)) +
  # geom_point() +
  geom_line(data=r2,aes(x=param_val,y=mean_pv,color=method)) +
  geom_errorbar(data=r2, aes(x=param_val, ymin=mean_pv-sd_pv, ymax=mean_pv+sd_pv, color=method),
                width = .1) +
  xlab(rname) +
  ylab('Mean -log10(pval)') +
  theme_bw() +
  theme(legend.position="none")
combined_plt <- cowplot::plot_grid(p1,p2,nrow=2,align = 'v')


combined_plt1 <- combined_plt
combined_plt2 <- combined_plt
combined_plt3 <- combined_plt
combined_plt4 <- combined_plt
combined_plt5 <- combined_plt


my_fig <- cowplot::plot_grid(plotlist = list(combined_plt1,combined_plt2,
                                                   combined_plt3,combined_plt4,
                                                   combined_plt5),
                                   nrow=1,align = 'h')

pdf(file = "/home/jmitchel/figures/for_paper_v2/sim_method_compare.pdf", useDingbats = FALSE,
    width = 10, height = 5.5)
my_fig
dev.off()



# generating sim for my fig 1 umap and pattern
sim_res <- generate_sim(n_ctypes=2,ctype_de_prob=.2,n_donors=8,
                        cells_per_don=25,n_genes=5000,n_process=1,
                        de_strength=.5,condition.prob=c(.25,.75),
                        similarity.scale=1) # works decently well


generate_sim_simple <- function(n_ctypes,ctype_de_prob,n_donors,cells_per_don,
                         n_genes,n_process,de_strength,condition.prob,
                         similarity.scale, sim_seed=1) {
  coldat_lst <- list()
  rowdat_lst <- list()
  sim_objs <- list()
  
  vcf <- mockVCF(n.samples = n_donors)
  gff <- mockGFF(n.genes = n_genes)
  
  seed <- sim_seed
  
  # make sure group.prob sums to 1
  gpr <- rep(1/n_ctypes, n_ctypes)
  gpr <- round(gpr,2)
  gpr[length(gpr)] <- 1 - sum(gpr[1:(length(gpr)-1)])
  
  # now trying for multi population two condition sim and adding more donors
  params.cond <- newSplatPopParams(similarity.scale = similarity.scale,
                                   condition.prob = condition.prob,
                                   cde.prob = c(0,.07),
                                   cde.downProb = c(0,.5),
                                   cde.facLoc = de_strength, 
                                   cde.facScale = 0.5,
                                   de.prob = ctype_de_prob,
                                   de.facLoc = 2.5, 
                                   de.facScale = 0.5,
                                   group.prob = gpr,
                                   batchCells = c(cells_per_don),
                                   batch.size = n_donors, nGenes=n_genes,
                                   eqtl.n = 10,
                                   seed = seed)
  
  sim.pop.cond1 <- splatPopSimulate(vcf = vcf, gff = gff,
                                    params = params.cond, 
                                    sparsify = FALSE)
  
  # need to make cell names unique per donor and group
  colnames(sim.pop.cond1) <- paste(sim.pop.cond1$Sample, sim.pop.cond1$Group, sim.pop.cond1$Cell, sep="_")
  
  sim.pop.cond1 <- logNormCounts(sim.pop.cond1)
  simcounts1 <- counts(sim.pop.cond1)
  simcounts1 <- methods::as(as.matrix(simcounts1),'sparseMatrix')
  simmeta1 <- as.data.frame(colData(sim.pop.cond1)@listData)
  colnames(simmeta1)[4] <- 'donors'
  colnames(simmeta1)[5] <- 'ctypes'
  
  ## alter some genes for one of the cell types
  # get CDE genes vec
  de_info <- rowData(sim.pop.cond1)
  de1 <- rownames(de_info)[de_info$ConditionDE.Condition1!=1]
  de2 <- rownames(de_info)[de_info$ConditionDE.Condition2!=1]
  all_de <- unique(c(de1,de2))
  num_de <- length(all_de)
  num_change <- round(num_de-20)
  if (num_change%%2==1) {
    num_change <- num_change - 1
  }
  all_change <- sample(all_de,num_change)
  change_ct1 <- sample(all_change,num_change/2)
  change_ct2 <- all_change[!(all_change %in% change_ct1)]
  ## all cells of ct1 should get values for these genes from some non-de genes
  non_de_mask <- rownames(de_info)[de_info$ConditionDE.Condition1==1]
  de_info_sub <- de_info[non_de_mask,]
  non_de <- rownames(de_info_sub)[de_info_sub$ConditionDE.Condition2==1]
  change_ct1_to <- sample(num_change/2)
  change_ct2_to <- sample(non_de,num_change/2)
  ct1_cells <- rownames(simmeta1)[simmeta1$ctypes=='Group1']
  ct2_cells <- rownames(simmeta1)[simmeta1$ctypes=='Group2']
  # now actually change counts
  simcounts1[change_ct1,ct1_cells] <- simcounts1[change_ct1_to,ct1_cells]
  simcounts1[change_ct2,ct2_cells] <- simcounts1[change_ct2_to,ct2_cells]
  
  ## trying to invert some genes to have opp expression in diff ctypes
  # first get genes up in cond2 down cond1
  de_info_sub <- de_info[de_info$ConditionDE.Condition2>1,]
  de_up_down <- rownames(de_info_sub)[de_info_sub$ConditionDE.Condition1==1]
  g_flip1 <- all_de[!(all_de %in% change_ct2)]
  g_flip1 <- sample(g_flip1,100)
  for (g in g_flip1) {
    cond1 <- de_info[g,'ConditionDE.Condition1']
    cond2 <- de_info[g,'ConditionDE.Condition2']
    if (cond1>1 && cond2==1) {
      print('yup')
      g_use <- sample(de_up_down,1)
      simcounts1[g,ct1_cells] <- simcounts1[g_use,ct1_cells]
    }
  }
  
  return(list(simcounts1,simmeta1))
}

sim_res <- generate_sim_simple(n_ctypes=2,ctype_de_prob=.15,n_donors=8,
                        cells_per_don=50,n_genes=5000,n_process=1,
                        de_strength=.3,condition.prob=c(.25,.75),
                        similarity.scale=7)
simcounts1 <- sim_res[[1]]
simmeta1 <- sim_res[[2]]

param_list <- initialize_params(ctypes_use = c("Group1", "Group2"),
                                ncores = 30, rand_seed = 10)

sim_container <- make_new_container(count_data=simcounts1, meta_data=simmeta1, params=param_list)

sim_container <- form_tensor(sim_container, donor_min_cells=3,
                             norm_method='trim', scale_factor=10000,
                             vargenes_method='norm_var_pvals', vargenes_thresh=1,
                             scale_var = TRUE, var_scale_power = .5)

sim_container <- run_tucker_ica(sim_container, ranks=c(2,4),
                                tucker_type = 'regular', rotation_type = 'hybrid')

sim_container <- plot_donor_matrix(sim_container, 
                                   show_donor_ids = TRUE)

sim_container$plots$donor_matrix

sim_container <- get_lm_pvals(sim_container)

pdf(file = "/home/jmitchel/figures/for_paper_v2/sim_f1_lds.pdf", useDingbats = FALSE,
    width = 3.5, height = 4.25)
sim_container <- plot_loadings_annot(sim_container, factor_select=1, use_sig_only=T, nonsig_to_zero=T,
                                     pathways=NULL, sig_thresh=0.15, display_genes=FALSE, 
                                     gene_callouts=FALSE, show_xlab=TRUE,
                                     show_var_explained=FALSE, reset_other_factor_plots=FALSE, draw_plot=TRUE)
dev.off()


sim_seurat <- CreateSeuratObject(simcounts1, project = "SeuratProject", assay = "RNA",
                   min.cells = 0, min.features = 0, names.field = 1,
                   names.delim = "_", meta.data = simmeta1)
sim_seurat <- NormalizeData(sim_seurat)
# sim_seurat <- FindVariableFeatures(sim_seurat, selection.method = "vst", nfeatures = 220)
sim_seurat <- FindVariableFeatures(sim_seurat, selection.method = "vst", nfeatures = 3000)
sim_seurat <- ScaleData(sim_seurat)
sim_seurat <- RunPCA(sim_seurat, features = VariableFeatures(object = sim_seurat))
# all.genes <- rownames(sim_seurat)
# sim_seurat <- RunPCA(sim_seurat, features = all.genes)
DimPlot(sim_seurat, reduction = "pca", group.by = 'donors')
DimPlot(sim_seurat, reduction = "pca", group.by = 'donors') + xlim(-25,25) + ylim(-25,25)

sim_seurat <- RunUMAP(sim_seurat, dims = 1:10)
DimPlot(sim_seurat, reduction = "umap", group.by = 'donors')

sim_seurat <- RunTSNE(sim_seurat, dims = 1:17)
pdf(file = "/home/jmitchel/figures/for_paper_v2/sim_tsne_f1.pdf", useDingbats = FALSE,
    width = 5, height = 3.75)
DimPlot(sim_seurat, reduction = "tsne", group.by = 'donors_new')
dev.off()

DimPlot(sim_seurat, reduction = "tsne", group.by = 'ctypes')


# changing names from sample to donor
sim_seurat@meta.data$donors_new <- sapply(as.character(sim_seurat@meta.data$donors), function(x) {
  s_num <- as.character(strsplit(x,split="_")[[1]][[2]])
  return(paste0('donor_',s_num))
})

## saving in case need to make visual changes later
# save.image(file='/home/jmitchel/data/sim_data/sim_dat_fig1.RData')
# load(file='/home/jmitchel/data/sim_data/sim_dat_fig1.RData')
























# ## trying to get the pairwise differences to plot instead
# all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_proc.rds')
# all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_proc_2.rds')
# all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_proc_3.rds')
# all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_donors.rds')
# all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_cells_per_don.rds')
# all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_de_strength.rds')
# all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_n_ctypes.rds')
# all_res <- readRDS(file='/home/jmitchel/data/sim_data/iter_similarity_scale.rds')
# 
# r1 <- all_res[[1]]
# r2 <- all_res[[2]]
# 
# # matrix to store results
# myres <- matrix(ncol=4,nrow=0)
# colnames(myres) <- c('mean_auc','mean_pv','method','param_val')
# 
# myres_summary <- matrix(ncol=6,nrow=0)
# colnames(myres_summary) <- c('mean_auc','sd_auc','mean_pv','sd_pv','method','param_val')
# 
# for (i in 1:nrow(r1)) {
#   start_ndx <- i%%3
#   if (start_ndx==1) {
#     auc_diff1 <- r1[i,'mean_auc'] - r1[i+1,'mean_auc']
#     auc_diff2 <- r1[i,'mean_auc'] - r1[i+2,'mean_auc']
#     pv_diff1 <- r1[i,'mean_pv'] - r1[i+1,'mean_pv']
#     pv_diff2 <- r1[i,'mean_pv'] - r1[i+2,'mean_pv']
#     row_add1 <- data.frame('mean_auc'=auc_diff1,'mean_pv'=pv_diff1,'method'='Tucker_PCA','param_val'=r1[i,'param_val'])
#     row_add2 <- data.frame('mean_auc'=auc_diff2,'mean_pv'=pv_diff2,'method'='Tucker_NMF','param_val'=r1[i,'param_val'])
#     myres <- rbind.data.frame(myres,row_add1,row_add2) # store results
#   }
#   if (i%%30==0) {
#     # calculate summar stats
#     myres_sub <- myres[myres$param_val==r1[i,'param_val'],]
#     
#     auc_mean = tapply(myres_sub$mean_auc, myres_sub$method, mean) 
#     pv_mean = tapply(myres_sub$mean_pv, myres_sub$method, mean) 
#     pv_mean <- pv_mean[names(auc_mean)] # ordering them same way
#     
#     auc_sd = tapply(myres_sub$mean_auc, myres_sub$method, sd) 
#     pv_sd = tapply(myres_sub$mean_pv, myres_sub$method, sd) 
#     pv_sd <- pv_sd[names(auc_sd)] # ordering them same way
#     
#     row_add <- data.frame('mean_auc'=auc_mean,'sd_auc'=auc_sd,'mean_pv'=pv_mean,'sd_pv'=pv_sd,'method'=names(pv_mean),'param_val'=r1[i,'param_val'])
#     myres_summary <- rbind.data.frame(myres_summary,row_add)
#   }
# }
# 
# r1 <- myres
# r2 <- myres_summary
# p1 <- ggplot(r1,aes(x=param_val,y=mean_auc,color=method)) +
#   # geom_point() +
#   geom_line(data=r2,aes(x=param_val,y=mean_auc,color=method)) +
#   geom_errorbar(data=r2, aes(x=param_val, ymin=mean_auc-sd_auc, ymax=mean_auc+sd_auc, color=method),
#                 width = .1) +
#   xlab('iter val') +
#   ylab('Mean AUC') +
#   theme_bw()
# p2 <- ggplot(r1,aes(x=param_val,y=mean_pv,color=method)) +
#   # geom_point() +
#   geom_line(data=r2,aes(x=param_val,y=mean_pv,color=method)) +
#   geom_errorbar(data=r2, aes(x=param_val, ymin=mean_pv-sd_pv, ymax=mean_pv+sd_pv, color=method),
#                 width = .1) +
#   xlab('iter val') +
#   ylab('Mean -log10(pval)') +
#   theme_bw()
# combined_plt <- cowplot::plot_grid(p1,p2,nrow=2,align = 'v')
# combined_plt







