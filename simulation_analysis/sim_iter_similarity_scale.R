library("splatter")
library("scater")
library(pROC)
library(NMF)
library(devtools)
load_all('/home/jmitchel/scITD')

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
      
      # # get AUC and pvals with Tucker
      # sim_container <- run_tucker_ica(sim_container, ranks=c(n_process,n_process*2),
      #                                 tucker_type = 'regular', rotation_type = 'hybrid')
      # 
      # all_tucker_res <- get_sim_auc(sim_container,coldat_lst,rowdat_lst)
      # tucker_auc <- mean(all_tucker_res[[1]])
      # tucker_pv <- mean(all_tucker_res[[2]])
      
      # get AUC and pvals with Tucker
      sim_container <- run_tucker_ica(sim_container, ranks=c(n_process,n_process*2),
                                      tucker_type = 'regular', rotation_type = 'hybrid')
      
      all_tucker_res <- get_sim_auc(sim_container,coldat_lst,rowdat_lst)
      tucker_auc1 <- mean(all_tucker_res[[1]])
      tucker_pv1 <- mean(all_tucker_res[[2]])
      
      # trying Tucker with slightly higher ranks
      sim_container <- run_tucker_ica(sim_container, ranks=c(n_process+1,(n_process+1)*2),
                                      tucker_type = 'regular', rotation_type = 'hybrid')
      
      all_tucker_res <- get_sim_auc(sim_container,coldat_lst,rowdat_lst)
      tucker_auc2 <- mean(all_tucker_res[[1]])
      tucker_pv2 <- mean(all_tucker_res[[2]])
      
      tucker_best <- order(c(tucker_auc1,tucker_auc2),decreasing=TRUE)[1]
      tucker_auc <- c(tucker_auc1,tucker_auc2)[tucker_best]
      tucker_pv <- c(tucker_pv1,tucker_pv2)[tucker_best]
      
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

r_seeds <- c(901,3001,12345,55453,60032)

# all_res <- get_performance_range(param_test='similarity.scale', param_range=c(.25,.5,1,1.5,2.5), n_iter=10, 
#                                  rand_seeds=r_seeds, n_ctypes=2, ctype_de_prob=.3, n_donors=50,
#                                  cells_per_don=50, n_genes=500, n_process=3,
#                                  de_strength=1, condition.prob=c(.25,.75),
#                                  similarity.scale=1.5)
# 
# saveRDS(all_res,file='/home/jmitchel/data/sim_data/iter_similarity_scale.rds')


## trying with 1.5 de_strength as well as 100 cells per donor 100 donors and multiple tries for Tucker
all_res <- get_performance_range(param_test='similarity.scale', param_range=c(.25,.5,1,1.5,2.5), n_iter=10, 
                                 rand_seeds=r_seeds, n_ctypes=2, ctype_de_prob=.3, n_donors=100,
                                 cells_per_don=100, n_genes=500, n_process=3,
                                 de_strength=1.5, condition.prob=c(.25,.75),
                                 similarity.scale=1.5)

saveRDS(all_res,file='/home/jmitchel/data/sim_data/iter_similarity_scale_2.rds')




