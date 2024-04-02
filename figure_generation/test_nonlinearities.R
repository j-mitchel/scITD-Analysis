.libPaths(c("/home/jmitchel/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))
.libPaths(c("/home/jmitchel/R/x86_64-pc-linux-gnu-library/4.1", .libPaths()))
.libPaths(c("/home/jmitchel/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))

library(pROC)
library(devtools)
load_all('/home/jmitchel/splatter')

# n_samples <- 10000
# r_vals1 <- sample(10:100, n_samples, replace=TRUE)
# # r_vals1 <- sample(5:10, n_samples, replace=TRUE)
# r_vals2 <- r_vals1**(2) + rnorm(n_samples,sd = 1)
# counts1 <- rpois(n_samples,lambda=r_vals1)
# # counts2 <- rpois(n_samples,lambda=r_vals2)
# counts2 <- rpois(n_samples,lambda=r_vals2)
# 
# plot(counts1,counts2)
# plot(log(counts1),log(counts2))
# 
# n_samples <- 10000
# r_vals1 <- sample(10:15, n_samples, replace=TRUE) + rnorm(n_samples,sd = 1)
# r_vals1[1:(n_samples/2)] <- r_vals1[1:(n_samples/2)]**(2)
# samp_ndx <- sample(1:n_samples,n_samples/2)
# r_vals1[samp_ndx] <- r_vals1[samp_ndx]**(.5)
# r_vals2 <- sample(10:15, n_samples, replace=TRUE) + rnorm(n_samples,sd = 1)
# r_vals2[1:(n_samples/2)] <- r_vals2[1:(n_samples/2)]**(1)
# r_vals2[samp_ndx] <- r_vals2[samp_ndx]**(1.5)
# counts1 <- rpois(n_samples,lambda=r_vals1)
# counts2 <- rpois(n_samples,lambda=r_vals2)
# 
# plot(counts1,counts2)
# plot(r_vals1,r_vals2)





generate_sim <- function(param_obj,share_process_genes,n_processes,n_ctype,sep_ct_r_effects,nonlinear) {
  vcf <- mockVCF(n.samples = 50)
  gff <- mockGFF(n.genes = 500)
  
  sim.sc <- splatPopSimulate(n_processes=n_processes,
                             share_process_genes=share_process_genes,
                             sep_ct_r_effects=sep_ct_r_effects,
                             vcf = vcf,
                             gff = gff,
                             params = param_obj,
                             sparsify = FALSE,
                             eqtl = NULL,
                             means = NULL,
                             key = NULL,
                             counts.only = FALSE,
                             nonlinear=nonlinear,
                             verbose = TRUE)
  
  # rename genes
  new_gn_names <- sapply(rownames(sim.sc),function(x){
    mysplit=strsplit(x,split='_')[[1]]
    return(paste0(mysplit[1],'-',mysplit[2]))
  })
  rownames(sim.sc) <- new_gn_names
  
  # rename cells
  new_c_names <- sapply(1:ncol(sim.sc),function(x){
    return(paste0('cell',x))
  })
  colnames(sim.sc) <- new_c_names
  
  # prep ground truth de gene info
  de_gene_info <- as.data.frame(sim.sc@rowRanges@elementMetadata[,(ncol(sim.sc@rowRanges@elementMetadata)-(2*n_ctype*n_processes)-3):(ncol(sim.sc@rowRanges@elementMetadata)-4)])
  rownames(de_gene_info) <- rownames(sim.sc)
  
  meta_include <- as.data.frame(sim.sc@colData[,c(4:ncol(sim.sc@colData))])
  colnames(meta_include)[1] <- 'donors'
  colnames(meta_include)[ncol(meta_include)] <- 'ctypes'
  
  return(list(sim.sc,de_gene_info,meta_include))
}


sim_run_scITD <- function(sim.sc,de_gene_info,meta_include,n_processes,n_ctype,use_PCA) {
  ## now running scITD on the simulation
  ctypes_use <- unique(meta_include$ctypes)
  param_list <- initialize_params(ctypes_use = ctypes_use,
                                  ncores = 5, rand_seed = 10)
  
  sim_container <- make_new_container(count_data=counts(sim.sc), meta_data=meta_include, params=param_list)
  
  sim_container <- form_tensor(sim_container, donor_min_cells=3,
                               norm_method='trim', scale_factor=10000,
                               vargenes_method='norm_var_pvals', vargenes_thresh=1,
                               scale_var = TRUE, var_scale_power = 0)
  
  if (use_PCA) {
    # or use PCA on the unfolded tensor
    sim_container <- pca_unfolded(sim_container, n_processes)
  } else {
    # get AUC and pvals with Tucker
    sim_container <- run_tucker_ica(sim_container, ranks=c(n_processes,n_processes*n_ctype),
                                    tucker_type = 'regular', rotation_type = 'hybrid')
  }
  
  ### now to compute AUC of gene prediction using loadings and true DE genes
  pr_labs <- sapply(colnames(de_gene_info),function(x){
    strsplit(x,split='.',fixed = TRUE)[[1]][[4]]
  })
  
  max_aucs <- c()
  for (proc_ndx in 1:n_processes) {
    all_auc <- c()
    for (factor_select in 1:n_processes) {
      lds <- get_one_factor(sim_container,factor_select = factor_select)[[2]]
      ndx_keep <- which(pr_labs==paste0('process',proc_ndx))
      de_pr_sub <- de_gene_info[,ndx_keep]
      group_labs <- sapply(colnames(de_pr_sub),function(x){
        strsplit(x,split='.',fixed = TRUE)[[1]][[3]]
      })
      ct_true_all <- c()
      for (gr in 1:length(unique(group_labs))) {
        ndx_keep <- which(group_labs==paste0('Group',gr))
        de_pr_sub_ct <- de_pr_sub[,ndx_keep]
        ndx_de <- rowSums(de_pr_sub_ct!=1)>0
        ct_true_all <- c(ct_true_all,ndx_de)
      }
      
      pROC_obj <- roc(ct_true_all,abs(c(lds)),
                      smoothed = TRUE,
                      plot=FALSE, AUC=TRUE)
      auc <- pROC_obj[["auc"]]
      all_auc <- c(all_auc,auc)
    }
    max_aucs <- c(max_aucs,max(all_auc))
  }
  
  ## now compute max factor score-process associations
  max_rsq <- c()
  for (proc_ndx in 1:n_processes) {
    all_rsq <- c()
    for (factor_select in 1:n_processes) {
      dsc <- get_one_factor(sim_container,factor_select = factor_select)[[1]]
      sim_container <- get_donor_meta(sim_container,additional_meta = paste0('Process',proc_ndx))
      tmp <- cbind.data.frame(sim_container$donor_metadata[rownames(dsc),],dsc)
      colnames(tmp)[2:3] <- c('process','dscore')
      lmres <- summary(lm(dscore~process,data=tmp))
      rsq <- lmres$r.squared
      all_rsq <- c(all_rsq,rsq)
    }
    max_rsq <- c(max_rsq,max(all_rsq))
  }
  return(list(max_aucs,max_rsq))
}

### trying to use simple non-linearity in splat-Pop
n_processes <- 2 # keeping this constant for now
n_ctype <- 2
n_iter <- 20
share_process_genes <- FALSE
sep_ct_r_effects <- TRUE

#### in splatPopSimConditionalEffects, added an option to exponentiate effects instead of multiply


n_iter <- 20
rand_seeds <- sample(1:1000000000,n_iter,replace = FALSE)

res_auc <- data.frame(matrix(nrow=0,ncol=2),stringsAsFactors=FALSE)
colnames(res_auc) <- c('auc','method')
res_rsq <- data.frame(matrix(nrow=0,ncol=2),stringsAsFactors=FALSE)
colnames(res_rsq) <- c('rsq','method')

for (rand_iter in 1:n_iter) {
  print(rand_iter)
  ## make a parameters object
  sim_seed <- rand_seeds[rand_iter]
  params_obj <- newSplatPopParams(similarity.scale = 6, # used 20 with other changes
                                  lib.loc = 8, ### TRYING WITH SMALLER LIB SIZES
                                  condition.prob = c(.5,.5), # donors split equally between conditions
                                  cde.prob = c(0,.1), # .1 is my default
                                  cde.downProb = c(0,.5),
                                  cde.facLoc = .1,
                                  cde.facScale = .5,
                                  de.prob = .3,
                                  de.facLoc = 2, # cross cell type de is stronger than condition de
                                  de.facScale = 0.5,
                                  group.prob = c(.5,.5),
                                  batchCells = 100, # total cells per donor
                                  batch.size = 50, # 50 donors
                                  nGenes=500, # 500 genes
                                  eqtl.n = 0,
                                  seed=sim_seed)
  
  sim_dat_lin <- generate_sim(params_obj,share_process_genes,n_processes,n_ctype,sep_ct_r_effects,nonlinear = FALSE)
  sim_dat_nonlin <- generate_sim(params_obj,share_process_genes,n_processes,n_ctype,sep_ct_r_effects,nonlinear = TRUE)
  
  sim.sc <- sim_dat_lin[[1]]
  de_gene_info <- sim_dat_lin[[2]]
  meta_include <- sim_dat_lin[[3]]
  
  scITD_res <- tryCatch({
    sim_run_scITD(sim.sc,de_gene_info,meta_include,n_processes,n_ctype,use_PCA=FALSE)
  }, error=function(cond) {
    return(NA)
  })
  
  res_auc[nrow(res_auc)+1,] <- list(mean(scITD_res[[1]]),'linear')
  res_rsq[nrow(res_rsq)+1,] <- list(mean(scITD_res[[2]]),'linear')
  
  
  sim.sc <- sim_dat_nonlin[[1]]
  de_gene_info <- sim_dat_nonlin[[2]]
  meta_include <- sim_dat_nonlin[[3]]
  
  scITD_res <- tryCatch({
    sim_run_scITD(sim.sc,de_gene_info,meta_include,n_processes,n_ctype,use_PCA=FALSE)
  }, error=function(cond) {
    return(NA)
  })
  
  res_auc[nrow(res_auc)+1,] <- list(mean(scITD_res[[1]]),'nonlinear')
  res_rsq[nrow(res_rsq)+1,] <- list(mean(scITD_res[[2]]),'nonlinear')
}

# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_nonlinear.rds')



myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_nonlinear.rds') # changing % genes in pattern

res_rsq <- myres[[2]]
res_rsq$method <- as.factor(res_rsq$method)
p1 <- ggplot(res_rsq,aes(x=method,y=rsq)) +
  geom_boxplot() +
  xlab('Simulation type') +
  ylab('r-squared') +
  # ylim(0,1) +
  theme_bw()

res_auc <- myres[[1]]
res_auc$method <- as.factor(res_auc$method)
p2 <- ggplot(res_auc,aes(x=method,y=auc)) +
  geom_boxplot() +
  xlab('Simulation type') +
  ylab('AUC') +
  ylim(0,1) +
  theme_bw()

library(cowplot)
fig <- plot_grid(p1,p2,nrow=1,align = 'h')

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/sim_nonlinear.pdf", useDingbats = FALSE,
#     width = 5, height = 2.75)
fig
dev.off()





