.libPaths(c("/home/jmitchel/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))
.libPaths(c("/home/jmitchel/R/x86_64-pc-linux-gnu-library/4.1", .libPaths()))
.libPaths(c("/home/jmitchel/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))

library(MASS)
library(DIALOGUE)
library(scater)
library(pROC)
library(devtools)
library(Seurat)
library(scITD)
library(cowplot)
library(ggplot2)
library(dplyr)
library(MOFAcellulaR)
library(devtools)
load_all('/home/jmitchel/splatter')
library(reticulate)
reticulate::use_condaenv("sandbox", required=TRUE)

generate_sim <- function(param_obj,share_process_genes,n_processes,n_ctype,sep_ct_r_effects) {
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

sim_run_dialogue <- function(sim.sc,de_gene_info,meta_include,n_processes,n_ctype) {
  # normalize expression
  sim.sc <- logNormCounts(sim.sc)

  # prep each cell type separately
  gr_list <- list()
  for (gr_ndx in 1:n_ctype) {
    gr_nm <- paste0('Group',gr_ndx)
    cells_keep <- rownames(sim.sc@colData)[sim.sc@colData$Group==gr_nm]
    sim.sc_sub <- sim.sc[,cells_keep]
    samples <- sim.sc_sub@colData$Sample
    tpm <- sim.sc_sub@assays@data@listData[["logcounts"]]
    sim.sc_sub <- runPCA(sim.sc_sub, ncomponents = 5, scale=TRUE)
    # sim.sc_sub <- runPCA(sim.sc_sub, ncomponents = 10, scale=TRUE) # my default for new sims
    X <- reducedDims(sim.sc_sub)@listData[["PCA"]]
    metadata <- as.data.frame(sim.sc_sub@colData[,c("Sample",
                                                    "ExpLibSize")])
    colnames(metadata)[2] <- 'cellQ'
    ct_dat <- make.cell.type(name = gr_nm,tpm,samples,X,metadata,cellQ = metadata$cellQ)
    gr_list[[gr_nm]] <- ct_dat
  }

  param <- DLG.get.param(k = n_processes,
                         results.dir = "/home/jmitchel/data/dialogue_results/",
                         conf = c("cellQ"),
                         abn.c = 2,
                         p.anova = 1,
                         find.genes = F)
  
  # ### testing just getting out 10 factors with dialogue
  # param <- DLG.get.param(k = 10,
  #                        results.dir = "/home/jmitchel/data/dialogue_results/",
  #                        conf = c("cellQ"),
  #                        abn.c = 2,
  #                        p.anova = 1,
  #                        find.genes = F)
  

  dial_res <- DIALOGUE.run(rA = gr_list, # list of cell.type objects
                           main = "sim.data99",
                           param = param)

  ### now to compute AUC of gene prediction using loadings and true DE genes
  pr_labs <- sapply(colnames(de_gene_info),function(x){
    strsplit(x,split='.',fixed = TRUE)[[1]][[4]]
  })

  max_aucs <- c()
  for (proc_ndx in 1:n_processes) {
    all_auc <- c()
    # for (factor_select in 1:n_processes) {
    for (factor_select in 1:ncol(dial_res[["cca.gene.cor1"]][[1]])) {
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

      lds <- c()
      for (gr in 1:length(unique(group_labs))) {
        g_lds <- dial_res$cca.gene.cor1[[paste0('Group',gr)]]
        g_lds_vals <- g_lds[,paste0('MCP',factor_select)]
        lds <- c(lds,g_lds_vals)
      }

      names(ct_true_all) <- NULL
      names(lds) <- NULL
      pROC_obj <- roc(ct_true_all,abs(lds),
                      smoothed = TRUE,
                      plot=FALSE, AUC=TRUE)
      auc <- pROC_obj[["auc"]]
      all_auc <- c(all_auc,auc)
    }
    max_aucs <- c(max_aucs,max(all_auc))
  }

  dscores <- lapply(names(dial_res[["sample.PCs"]]), function(i) dial_res[["sample.PCs"]][[i]] %*% dial_res$cca$ws[[i]])
  dscores_cross_ct_av <- matrix(0,ncol=ncol(dscores[[1]]),nrow=nrow(dscores[[1]]))
  for (i in 1:length(dscores)) {
    dscores_cross_ct_av <- dscores_cross_ct_av + dscores[[i]]
  }
  dscores_cross_ct_av <- dscores_cross_ct_av / length(dscores)

  ## now compute max factor score-process associations
  d_meta <- unique(meta_include[,c(1:(ncol(meta_include)-1))])
  rownames(d_meta) <- d_meta$donors
  max_rsq <- c()
  for (proc_ndx in 1:n_processes) {
    all_rsq <- c()
    for (factor_select in 1:n_processes) {
      dsc <- dscores_cross_ct_av[,factor_select]
      tmp <- cbind.data.frame(d_meta[names(dsc),proc_ndx+1],dsc)
      colnames(tmp)[1:2] <- c('process','dscore')
      lmres <- summary(lm(dscore~process,data=tmp))
      rsq <- lmres$r.squared
      all_rsq <- c(all_rsq,rsq)
    }
    max_rsq <- c(max_rsq,max(all_rsq))
  }
  return(list(max_aucs,max_rsq))
}


sim_run_MOFA <- function(sim.sc,de_gene_info,meta_include,n_processes,n_ctype) {
  sim_pb <- scuttle::summarizeAssayByGroup(sim.sc,sim.sc@colData[,c('Sample','Group')],statistics='sum')
  sim_pb_meta <- as.data.frame(sim_pb@colData@listData)
  sim_pb_counts <- sim_pb@assays@data@listData[["sum"]]

  colnames(sim_pb_meta) <- c('donor_id','cell_type','cell_counts')
  concat_nms <- paste0(sim_pb_meta$cell_type,'_',sim_pb_meta$donor_id)
  rownames(sim_pb_meta) <- concat_nms
  colnames(sim_pb_counts) <- concat_nms

  sim_pb_obj <- MOFAcellulaR::create_init_exp(counts = sim_pb_counts,  coldata = sim_pb_meta)

  sim_ct_list <- MOFAcellulaR::filt_profiles(pb_dat = sim_pb_obj,
                                             cts = c("Group1","Group2"),
                                             ncells = 0, # don't remove any samples for not having enough cells
                                             counts_col = "cell_counts", # This refers to the column name in testcoldata where the number of cells per profile was stored
                                             ct_col = "cell_type") # This refers to the column name in testcoldata where the cell-type label was stored

  sim_ct_list <- MOFAcellulaR::tmm_trns(pb_dat_list = sim_ct_list,
                                        scale_factor = 1000000)

  ## convert the list to a MOFA object
  sim_multiview_dat <- pb_dat2MOFA(pb_dat_list = sim_ct_list,
                                   sample_column = "donor_id")




  ### running the MOFA model
  sim_MOFAobject <- MOFA2::create_mofa(sim_multiview_dat)

  data_opts <- MOFA2::get_default_data_options(sim_MOFAobject)
  train_opts <- MOFA2::get_default_training_options(sim_MOFAobject)
  model_opts <- MOFA2::get_default_model_options(sim_MOFAobject)

  # This avoids the regularization of multicellular programs per cell type.
  # This avoids less sparse gene weights
  model_opts$spikeslab_weights <- FALSE

  # Define the number of factors needed
  model_opts$num_factors <- n_processes

  # Prepare MOFA model:
  sim_MOFAobject <- MOFA2::prepare_mofa(object = sim_MOFAobject,
                                        data_options = data_opts,
                                        model_options = model_opts,
                                        training_options = train_opts)

  mofa_model <- MOFA2::run_mofa(sim_MOFAobject, outfile=NULL)


  ### now to compute AUC of gene prediction using loadings and true DE genes
  pr_labs <- sapply(colnames(de_gene_info),function(x){
    strsplit(x,split='.',fixed = TRUE)[[1]][[4]]
  })

  max_aucs <- c()
  for (proc_ndx in 1:n_processes) {
    all_auc <- c()
    for (factor_select in 1:n_processes) {
      mofa_lds <- MOFAcellulaR::get_geneweights(model = mofa_model, factor = paste0("Factor",factor_select))
      lds <- mofa_lds[,'value']
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

      pROC_obj <- roc(ct_true_all,abs(lds),
                      smoothed = TRUE,
                      plot=FALSE, AUC=TRUE)
      auc <- pROC_obj[["auc"]]
      all_auc <- c(all_auc,auc)
    }
    max_aucs <- c(max_aucs,max(all_auc))
  }

  ## now compute max factor score-process associations
  d_meta <- unique(meta_include[,c(1:(ncol(meta_include)-1))])
  rownames(d_meta) <- d_meta$donors
  colnames(d_meta)[1] <- c('sample')
  rownames(d_meta) <- NULL
  max_rsq <- c()
  for (proc_ndx in 1:n_processes) {
    all_rsq <- c()
    for (factor_select in 1:n_processes) {
      mofa_dsc <- MOFAcellulaR::get_tidy_factors(model = mofa_model,
                                                 metadata = d_meta,
                                                 factor = paste0("Factor",factor_select),
                                                 sample_id_column = "sample")
      tmp <- as.data.frame(mofa_dsc[,c(paste0('Process',proc_ndx),'value')])
      colnames(tmp)[1:2] <- c('process','dscore')
      lmres <- summary(lm(dscore~process,data=tmp))
      rsq <- lmres$r.squared
      all_rsq <- c(all_rsq,rsq)
    }
    max_rsq <- c(max_rsq,max(all_rsq))
  }

  return(list(max_aucs,max_rsq))
}

## pick some parameter to change
n_processes <- 2 # keeping this constant for now
n_ctype <- 2
n_iter <- 20
share_process_genes <- FALSE
# share_process_genes <- TRUE
sep_ct_r_effects <- TRUE
# sep_ct_r_effects <- FALSE
rand_seeds <- sample(1:1000000000,n_iter,replace = FALSE)
# param_range <- c(.001,.01,.1,1) # for de strength
# param_range <- c(.01,.05,.1,.15,.2) # for de %
# param_range <- c(2,3,4,5,6) # for variable numbers of processes
param_range <- c(1,6,11,16,21) # for different donor-donor similarity
res_auc <- data.frame(matrix(nrow=0,ncol=3),stringsAsFactors=FALSE)
colnames(res_auc) <- c('auc','method','param_val')
res_rsq <- data.frame(matrix(nrow=0,ncol=3),stringsAsFactors=FALSE)
colnames(res_rsq) <- c('rsq','method','param_val')
for (param_ndx in 1:length(param_range)) {
  param_val <- param_range[param_ndx]
  sims_save <- list()
  for (rand_iter in 1:n_iter) {
    print(param_ndx)
    print(rand_iter)
    sim_seed <- rand_seeds[rand_iter]

    # ### special for changing number of processes since it's used in multiple spots
    # n_processes <- param_val

    ## make a parameters object
    params_obj <- newSplatPopParams(similarity.scale = param_val, # used 20 with other changes
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

    # generate simulations
    sim_dat_all <- generate_sim(params_obj,share_process_genes,n_processes,n_ctype,sep_ct_r_effects)
    # sims_save[[rand_iter]] <- sim_dat_all
    sim.sc <- sim_dat_all[[1]]
    de_gene_info <- sim_dat_all[[2]]
    meta_include <- sim_dat_all[[3]]

    # ##### for downsampling cells
    # cells_keep <- sample(colnames(sim.sc),500)
    # sim.sc <- sim.sc[,cells_keep]
    # meta_include <- meta_include[cells_keep,]
    # #####

    # get results from each tool
    scITD_res <- tryCatch({
      sim_run_scITD(sim.sc,de_gene_info,meta_include,n_processes,n_ctype,use_PCA=FALSE)
    }, error=function(cond) {
      return(NA)
    })
    PCA_res <- tryCatch({
      sim_run_scITD(sim.sc,de_gene_info,meta_include,n_processes,n_ctype,use_PCA=TRUE)
    }, error=function(cond) {
      return(NA)
    })
    dialogue_res <- tryCatch({
      sim_run_dialogue(sim.sc,de_gene_info,meta_include,n_processes,n_ctype)
    }, error=function(cond) {
      return(NA)
    })
    MOFA_res <- tryCatch({
      sim_run_MOFA(sim.sc,de_gene_info,meta_include,n_processes,n_ctype)
    }, error=function(cond) {
      return(NA)
    })
    
    # invisible(utils::capture.output(
    #   scITD_res <- sim_run_scITD(sim.sc,de_gene_info,meta_include,n_processes,n_ctype,use_PCA=FALSE)
    # ))
    # invisible(utils::capture.output(
    #   PCA_res <- sim_run_scITD(sim.sc,de_gene_info,meta_include,n_processes,n_ctype,use_PCA=TRUE)
    # ))
    # invisible(utils::capture.output(
    #   dialogue_res <- sim_run_dialogue(sim.sc,de_gene_info,meta_include,n_processes,n_ctype)
    # ))
    # invisible(utils::capture.output(
    #   MOFA_res <- sim_run_MOFA(sim.sc,de_gene_info,meta_include,n_processes,n_ctype)
    # ))

    # store results
    if (is.na(scITD_res)[[1]][1]) {
      res_auc[nrow(res_auc)+1,] <- list(NA,'scITD',param_val)
      res_rsq[nrow(res_rsq)+1,] <- list(NA,'scITD',param_val)
    } else {
      res_auc[nrow(res_auc)+1,] <- list(mean(scITD_res[[1]]),'scITD',param_val)
      res_rsq[nrow(res_rsq)+1,] <- list(mean(scITD_res[[2]]),'scITD',param_val)
    }
    
    if (is.na(PCA_res)[[1]][1]) {
      res_auc[nrow(res_auc)+1,] <- list(NA,'PCA',param_val)
      res_rsq[nrow(res_rsq)+1,] <- list(NA,'PCA',param_val)
    } else {
      res_auc[nrow(res_auc)+1,] <- list(mean(PCA_res[[1]]),'PCA',param_val)
      res_rsq[nrow(res_rsq)+1,] <- list(mean(PCA_res[[2]]),'PCA',param_val)
    }
    
    if (is.na(dialogue_res)[[1]][1]) {
      res_auc[nrow(res_auc)+1,] <- list(NA,'dialogue',param_val)
      res_rsq[nrow(res_rsq)+1,] <- list(NA,'dialogue',param_val)
    } else {
      res_auc[nrow(res_auc)+1,] <- list(mean(dialogue_res[[1]]),'dialogue',param_val)
      res_rsq[nrow(res_rsq)+1,] <- list(mean(dialogue_res[[2]]),'dialogue',param_val)
    }
    
    if (is.na(MOFA_res)[[1]][1]) {
      res_auc[nrow(res_auc)+1,] <- list(NA,'MOFA',param_val)
      res_rsq[nrow(res_rsq)+1,] <- list(NA,'MOFA',param_val)
    } else {
      res_auc[nrow(res_auc)+1,] <- list(mean(MOFA_res[[1]]),'MOFA',param_val)
      res_rsq[nrow(res_rsq)+1,] <- list(mean(MOFA_res[[2]]),'MOFA',param_val)
    }




    # delete temporary file created by dialogue
    file.remove('/home/jmitchel/data/dialogue_results/DIALOGUE1_sim.data99.rds')
    unlink("/home/jmitchel/data/dialogue_results/DIALOGUE2_sim.data99/",recursive=TRUE)
  }
  # saveRDS(sims_save,file=paste0('/home/jmitchel/data/scITD_sim_res/sim_dat11_',param_ndx,'.rds'))
  # rm(sims_save)
  # gc()
}


# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res1.rds') # changing de strength param
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res2.rds') # changing de % param
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res3.rds') # changing num ctypes
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res4.rds') # changing de % param but with same genes per cell type
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res5.rds') # changing the similarity.scale parameter (~1/mean of the dist from which the cross-donor variance is samples)
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res6.rds') # changing de % param - using smaller library sizes this time
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res7.rds') # changing number of ground truth factors (dialogue w 10PCs now)
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res8.rds') # changing de % param, using different random effects per cell type
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res9.rds') # changing donor similarity param, using different random effects per cell type, realistic lib sizes
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res10.rds') # changing donor similarity param, using same random effects per cell type, realistic lib sizes
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res11.rds') # changing donor similarity param, using different random effects per cell type, realistic lib sizes, dialogue gets 10pcs
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res12.rds') # changing number of processes, using different random effects per cell type, high donor similarity realistic lib sizes
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res13.rds') # changing DE%, using different random effects per cell type, high donor similarity realistic lib sizes
# saveRDS(list(res_auc,res_rsq),file='/home/jmitchel/data/scITD_sim_res/sim_res14.rds') # changing donor similarity param, using different random effects per cell type, realistic lib sizes, using fewer PCs for dialoge



.libPaths(c("/home/jmitchel/R/x86_64-pc-linux-gnu-library/4.0", .libPaths()))
.libPaths(c("/home/jmitchel/R/x86_64-pc-linux-gnu-library/4.1", .libPaths()))
.libPaths(c("/home/jmitchel/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))

library(cowplot)
library(ggplot2)
library(dplyr)

# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res1.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res2.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res3.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res4.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res5.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res6.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res7.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res8.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res9.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res10.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res11.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res12.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res13.rds')
# myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res14.rds')

myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res13.rds') # changing % genes in pattern

res_rsq <- myres[[2]]

n_iter <- 20

res_rsq_sumstats <- res_rsq %>%
  group_by(method,param_val) %>%
  summarize(mean = mean(rsq,na.rm=TRUE),
            sd = sd(rsq,na.rm=TRUE),
            se = sd(rsq,na.rm=TRUE)/sqrt(n_iter))

res_rsq_sumstats$method <- as.factor(res_rsq_sumstats$method)
p1 <- ggplot(res_rsq_sumstats,aes(x=param_val,y=mean,color=method)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se,fill=method,alpha=0.5),show.legend=FALSE) +
  geom_point(size = 2) +
  xlab('% genes in multicellular pattern') +
  ylab('r-squared') +
  ylim(0,1) +
  theme_bw()



myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res12.rds') # changing number of patterns

res_rsq <- myres[[2]]

n_iter <- 20

res_rsq_sumstats <- res_rsq %>%
  group_by(method,param_val) %>%
  summarize(mean = mean(rsq,na.rm=TRUE),
            sd = sd(rsq,na.rm=TRUE),
            se = sd(rsq,na.rm=TRUE)/sqrt(n_iter))

res_rsq_sumstats$method <- as.factor(res_rsq_sumstats$method)
p2 <- ggplot(res_rsq_sumstats,aes(x=param_val,y=mean,color=method)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se,fill=method,alpha=0.5),show.legend=FALSE) +
  geom_point(size = 2) +
  xlab('Number of multicellular patterns') +
  ylab('r-squared') +
  ylim(0,1) +
  theme_bw()

fig1 <- plot_grid(p1,p2,nrow=1)

### Figure 3a
pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/sim_fig3.pdf", useDingbats = FALSE,
    width = 8, height = 2.75)
fig1
dev.off()


##### now to plot auc vals
myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res13.rds') # changing % genes in pattern

res_auc <- myres[[1]]
n_iter <- 20

res_auc_sumstats <- res_auc %>%
  group_by(method,param_val) %>%
  summarize(mean = mean(auc,na.rm=TRUE),
            sd = sd(auc,na.rm=TRUE),
            se = sd(auc,na.rm=TRUE)/sqrt(n_iter))

res_auc_sumstats$method <- as.factor(res_auc_sumstats$method)
p1 <- ggplot(res_auc_sumstats,aes(x=param_val,y=mean,color=method)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se,fill=method,alpha=0.5),show.legend=FALSE) +
  geom_point(size = 2) +
  xlab('% genes in multicellular pattern') +
  ylab('AUC') +
  ylim(.5,1)


myres <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res12.rds') # changing number of patterns

res_auc <- myres[[1]]
n_iter <- 20

res_auc_sumstats <- res_auc %>%
  group_by(method,param_val) %>%
  summarize(mean = mean(auc,na.rm=TRUE),
            sd = sd(auc,na.rm=TRUE),
            se = sd(auc,na.rm=TRUE)/sqrt(n_iter))

res_auc_sumstats$method <- as.factor(res_auc_sumstats$method)
p2 <- ggplot(res_auc_sumstats,aes(x=param_val,y=mean,color=method)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se,fill=method,alpha=0.5),show.legend=FALSE) +
  geom_point(size = 2) +
  xlab('Number of multicellular patterns') +
  ylab('AUC') +
  ylim(.5,1)


fig2 <- plot_grid(p1,p2,nrow=1)



### Figure S6b
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/sim_num_factors.pdf", useDingbats = FALSE,
#     width = 8, height = 2.75)
# fig1
# dev.off()
# 

### Figure S6c
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/sim_de_percent.pdf", useDingbats = FALSE,
#     width = 8, height = 2.75)
# fig1
# dev.off()




#### plotting dialogue versus and dialogue max of 10 factors on same plot versus donor-donor similarity - sims 9 and 11
myres1 <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res9.rds')
myres2 <- readRDS(file='/home/jmitchel/data/scITD_sim_res/sim_res11.rds')

myres1[[1]]$method <- as.factor(myres1[[1]]$method)
levels(myres1[[1]]$method) <- c("dialogue_2MCP","MOFA","PCA","scITD")
myres1[[2]]$method <- as.factor(myres1[[2]]$method)
levels(myres1[[2]]$method) <- c("dialogue_2MCP","MOFA","PCA","scITD")

myres2[[1]]$method <- as.factor(myres2[[1]]$method)
levels(myres2[[1]]$method) <- c("dialogue_10MCP","MOFA","PCA","scITD")
myres2[[2]]$method <- as.factor(myres2[[2]]$method)
levels(myres2[[2]]$method) <- c("dialogue_10MCP","MOFA","PCA","scITD")

auc_sub <- myres2[[1]][as.character(myres2[[1]]$method)=='dialogue_10MCP',]
rsq_sub <- myres2[[2]][as.character(myres2[[1]]$method)=='dialogue_10MCP',]

myres1[[1]] <- rbind.data.frame(myres1[[1]],auc_sub)
myres1[[2]] <- rbind.data.frame(myres1[[2]],rsq_sub)

res_auc <- myres1[[1]]
res_rsq <- myres1[[2]]

n_iter <- 20

res_rsq_sumstats <- res_rsq %>%
  group_by(method,param_val) %>%
  summarize(mean = mean(rsq,na.rm=TRUE),
            sd = sd(rsq,na.rm=TRUE),
            se = sd(rsq,na.rm=TRUE)/sqrt(n_iter))

res_auc_sumstats <- res_auc %>%
  group_by(method,param_val) %>%
  summarize(mean = mean(auc,na.rm=TRUE),
            sd = sd(auc,na.rm=TRUE),
            se = sd(auc,na.rm=TRUE)/sqrt(n_iter))


res_rsq_sumstats$method <- as.factor(res_rsq_sumstats$method)
p1 <- ggplot(res_rsq_sumstats,aes(x=param_val,y=mean,color=method)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se,fill=method,alpha=0.5),show.legend=FALSE) +
  geom_point(size = 2) +
  xlab('% genes in multicellular pattern') +
  # xlab('Donor-donor similarity (higher=more similar)') +
  # xlab('Number of multicellular patterns') +
  ylab('r-squared') +
  ylim(0,1)


res_auc_sumstats$method <- as.factor(res_auc_sumstats$method)
p2 <- ggplot(res_auc_sumstats,aes(x=param_val,y=mean,color=method)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se,fill=method,alpha=0.5),show.legend=FALSE) +
  geom_point(size = 2) +
  xlab('% genes in multicellular pattern') +
  # xlab('Donor-donor similarity (higher=more similar)') +
  # xlab('Number of multicellular patterns') +
  ylab('AUC') +
  ylim(.5,1)

fig1 <- plot_grid(p1,p2,nrow=1)

### Figure S6d
pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/sim_dialogue_10MCP.pdf", useDingbats = FALSE,
    width = 8, height = 2.75)
fig1
dev.off()





