library(Seurat)
library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(nichenetr)
library(dplyr)
library(scITD)


### This script is intended for comparing the results and performance of scITD LR analysis
# to that of Differential NicheNet

##### prepping the data the same way as in LR_analysis.R except that using NicheNet LR pairs
# load up the lupus dataset: see preprocessing/lupus_preprocessing.R 
# for code used to generate this object
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')

# converting shorthand cell type names to full names
new_names <- sapply(as.character(pbmc@meta.data$cg_cov), function(x){
  if (x=='cM') {
    return('cMono')
  } else if (x=='ncM') {
    return('ncMono')
  } else if (x=='T4') {
    return('Th')
  } else if (x=='T8') {
    return('Tc')
  } else {
    return(x)
  }
})
names(new_names) <- NULL
pbmc@meta.data$cg_cov <- factor(new_names,levels=unique(new_names))

# subset data to SLE patients only
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$Status=='Managed']
pbmc <- subset(pbmc,cells = cells_keep)

param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc","cDC",
                                               "cMono","ncMono"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(seurat_obj=pbmc,
                                     params=param_list,
                                     metadata_cols=c('ind_cov_batch_cov',
                                                     "SLE_status",
                                                     "Status",
                                                     "cg_cov",
                                                     "sex",
                                                     "Age",
                                                     "batch_cov",
                                                     "Processing_Cohort",
                                                     "Ethnicity"),
                                     metadata_col_nm=c('donors',
                                                       'SLE_status',
                                                       'Status',
                                                       'ctypes',
                                                       'sex',
                                                       'Age',
                                                       'pool',
                                                       'processing',
                                                       'Ethnicity'))


pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20),
                                 tucker_type = 'regular', rotation_type = 'hybrid')

pbmc_container$tucker_results[[1]][,1] <- pbmc_container$tucker_results[[1]][,1] * -1
pbmc_container$tucker_results[[2]][1,] <- pbmc_container$tucker_results[[2]][1,] * -1

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Ethnicity'),
                                        stat_use='pval')

## plot donor score
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# just to check that everything is as expected up to this point
pbmc_container$plots$donor_matrix



# function to get nichenet results for a factor
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
# ligand_target_matrix = readRDS("/home/jmitchel/data/NicheNet/ligand_target_matrix.rds")
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
# saveRDS(lr_network,file='/home/jmitchel/data/NicheNet/lr_network.rds')
# lr_network = readRDS("/home/jmitchel/data/NicheNet/lr_network.rds")
lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)

get_nichenet <- function(pbmc,pbmc_container,factor_test,ligand_target_matrix,lr_network) {
  dsc <- pbmc_container$tucker_results[[1]][,factor_test]
  
  # using top/bot quintiles
  f_qt <- quantile(dsc, probs = seq(0, 1, 0.2), names = TRUE)
  
  f_top <- names(dsc)[dsc>f_qt['80%']]
  f_bot <- names(dsc)[dsc<f_qt['20%']]
  
  f_cells <- sapply(1:nrow(pbmc@meta.data),function(i){
    if (pbmc@meta.data[i,'ind_cov_batch_cov'] %in% f_top) {
      return('high')
    } else if (pbmc@meta.data[i,'ind_cov_batch_cov'] %in% f_bot) {
      return('low')
    } else {
      return(NA)
    }
  })
  
  pbmc@meta.data$f_cells_test <- f_cells
  
  # now subsetting seurat object by those not NA for the factor
  cells_keep <- rownames(pbmc@meta.data)[!is.na(pbmc@meta.data$f_cells_test)]
  pbmc_sub <- subset(pbmc,cells=cells_keep)
  
  # removing unused cell types
  cells_keep <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$cg_cov %in% c('Th','Tc','cMono','ncMono','NK','B','cDC')]
  seurat_obj <- subset(pbmc_sub,cells=cells_keep)
  
  seurat_obj@meta.data$celltype_aggregate = paste(seurat_obj@meta.data$cg_cov, seurat_obj@meta.data$f_cells_test,sep = "_")
  
  celltype_id = "celltype_aggregate" # metadata column name of the cell type of interest
  seurat_obj = SetIdent(seurat_obj, value = seurat_obj[[celltype_id]])
  
  all_ctypes <- pbmc_container$experiment_params$ctypes_use
  res_lst <- list()
  for (ct in all_ctypes) {
    print(ct)
    # ct is the target ctype
    ct_other <- all_ctypes[all_ctypes!=ct]
    senders_high <- sapply(ct_other,function(x){
      paste0(x,"_high")
    })
    names(senders_high) <- NULL
    senders_low <- sapply(ct_other,function(x){
      paste0(x,"_low")
    })
    names(senders_low) <- NULL
    
    niches = list(
      "High_niche" = list(
        "sender" = senders_high,
        "receiver" = paste0(ct,"_high")),
      "Low_niche" = list(
        "sender" = senders_low,
        "receiver" = paste0(ct,"_low"))
    ) # user adaptation required on own dataset
    
    assay_oi = "RNA"
    DE_sender = nichenetr:::calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi) # only ligands important for sender cell types
    DE_receiver = nichenetr:::calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi) # only receptors now, later on: DE analysis to find targets
    
    expression_pct = 0.10
    DE_sender_processed = nichenetr:::process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
    DE_receiver_processed = nichenetr:::process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")
    
    specificity_score_LR_pairs = "min_lfc"
    DE_sender_receiver = nichenetr:::combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)
    
    include_spatial_info_sender = FALSE # if not spatial info to include: put this to false # user adaptation required on own dataset
    include_spatial_info_receiver = FALSE # if spatial info to include: put this to true # user adaptation required on own dataset
    
    spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
    
    sender_spatial_DE_processed = nichenetr:::get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "sender", lr_network = lr_network)
    sender_spatial_DE_processed = sender_spatial_DE_processed %>% mutate(scaled_ligand_score_spatial = nichenetr:::scale_quantile_adapted(ligand_score_spatial))
    receiver_spatial_DE_processed = get_non_spatial_de(niches = niches, spatial_info = spatial_info, type = "receiver", lr_network = lr_network)
    receiver_spatial_DE_processed = receiver_spatial_DE_processed %>% mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
    
    
    lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff. 
    specificity_score_targets = "min_lfc"
    
    DE_receiver_targets = nichenetr:::calculate_niche_de_targets(seurat_obj = seurat_obj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 
    DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = specificity_score_targets)
    
    background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
    geneset_niche1 = DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
    geneset_niche2 = DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
    
    
    top_n_target = 250
    
    niche_geneset_list = list(
      "High_niche" = list(
        "receiver" = niches[[1]]$receiver,
        "geneset" = geneset_niche1,
        "background" = background),
      "Low_nich" = list(
        "receiver" = niches[[2]]$receiver,
        "geneset" = geneset_niche2 ,
        "background" = background)
    )
    
    ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)
    
    features_oi = union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)
    
    dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
    exprs_tbl = dotplot$data %>% as_tibble()
    exprs_tbl = exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
      mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
    
    exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
    exprs_tbl_receptor = exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
    exprs_tbl_target = exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)
    
    exprs_tbl_ligand = exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
    
    exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
    
    exprs_sender_receiver = lr_network %>% 
      inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
      inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))
    
    ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% distinct() %>% ungroup() 
    
    prioritizing_weights = c("scaled_ligand_score" = 5,
                             "scaled_ligand_expression_scaled" = 1,
                             "ligand_fraction" = 1,
                             "scaled_ligand_score_spatial" = 0, 
                             "scaled_receptor_score" = 0.5,
                             "scaled_receptor_expression_scaled" = 0.5,
                             "receptor_fraction" = 1, 
                             "ligand_scaled_receptor_expression_fraction" = 1,
                             "scaled_receptor_score_spatial" = 0,
                             "scaled_activity" = 0,
                             "scaled_activity_normalized" = 1,
                             "bona_fide" = 1)
    
    output = list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
                  ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
    prioritization_tables = get_prioritization_tables(output, prioritizing_weights)
    
    res_lst[[ct]] <- prioritization_tables
  }
  return(res_lst)
}

# since I already have the data for factor 2, just need to run it for f1, f3
res_lst_f1 <- get_nichenet(pbmc,pbmc_container,factor_test=1,ligand_target_matrix,lr_network)
# saveRDS(res_lst_f1,file='/home/jmitchel/data/lupus_data/nichenet_f1.rds')

res_lst_f2 <- get_nichenet(pbmc,pbmc_container,factor_test=2,ligand_target_matrix,lr_network)
# saveRDS(res_lst_f2,file='/home/jmitchel/data/lupus_data/nichenet_f2.rds')

res_lst_f3 <- get_nichenet(pbmc,pbmc_container,factor_test=3,ligand_target_matrix,lr_network)
# saveRDS(res_lst_f3,file='/home/jmitchel/data/lupus_data/nichenet_f3.rds')

res_lst_f1 <- readRDS(file='/home/jmitchel/data/lupus_data/nichenet_f1.rds')
res_lst_f2 <- readRDS(file='/home/jmitchel/data/lupus_data/nichenet_f2.rds')
res_lst_f3 <- readRDS(file='/home/jmitchel/data/lupus_data/nichenet_f3.rds')




#### running scITD LR analysis using same LR pairs
lr_pairs <- as.matrix(lr_network)
lr_pairs <- lr_pairs[lr_pairs[,'bonafide']=='TRUE',]
lr_pairs <- lr_pairs[,c(1,2)]

# infer active LR interactions
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='pool')
sft_thresh <- c(12,14,12,10,12,9,12)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

lr_hmap <- compute_LR_interact(pbmc_container, lr_pairs, sig_thresh=.00000000005,
                               percentile_exp_rec=0.85, add_ld_fact_sig=TRUE)

lr_hmap




# for ligands with a significant interaction in at least 1 cell type, get the ones that
# are also associated with a factor of interest
get_lr_dat <- function(container,sig_thresh,factor_select) {
  myres_mat <- container[["lr_res"]]
  
  # remove rows where ligand not in the data
  lig_names <- sapply(rownames(myres_mat),function(x){
    strsplit(x,split='_')[[1]][[1]]
  })
  ndx_keep <- which(lig_names %in% colnames(container$scale_pb_extra[[1]]))
  myres_mat <- myres_mat[ndx_keep,]
  
  ##### from get_LR_interact
  # reduce to rows/columns with at least one significant hit
  myres_mat <- myres_mat[rowSums(myres_mat<sig_thresh)>0,]
  myres_mat <- myres_mat[,colSums(myres_mat<sig_thresh)>0]
  
  # log transform values
  myres_mat <- -log10(myres_mat)
  
  # put na values back where source ct == target ct
  rs <- sapply(rownames(myres_mat),function(x){
    strsplit(x,split='_')[[1]][[2]]
  })
  rs <- factor(rs,levels=unique(rs))
  cs <- sapply(colnames(myres_mat),function(x){
    strsplit(x,split='_')[[1]][[1]]
  })
  cs <- factor(cs,levels=levels(rs))
  myres_mat[outer(rs, cs, "==")] <- NA
  
  fact_res <- matrix(nrow=nrow(myres_mat),ncol=ncol(container$tucker_results[[1]]))
  colnames(fact_res) <- sapply(1:ncol(container$tucker_results[[1]]),function(x){
    paste0('Factor_',x)
  })
  for (i in 1:ncol(container$tucker_results[[1]])) {
    for (j in 1:nrow(myres_mat)) {
      lig <- strsplit(rownames(myres_mat)[j],split='_')[[1]][[1]]
      ct <- strsplit(rownames(myres_mat)[j],split='_')[[1]][[2]]
      
      lig_ct_exp <- container$scale_pb_extra[[ct]][,lig]
      
      if (sum(lig_ct_exp!=0)==0) {
        pval <- NA
      } else {
        tmp <- cbind.data.frame(container$tucker_results[[1]][names(lig_ct_exp),i],lig_ct_exp)
        colnames(tmp) <- c('dsc','l_exp')
        lmres <- lm(dsc~l_exp,data=tmp)
        lmres <- summary(lmres)
        pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
        fact_res[j,i] <- pval
      }
    }
  }
  # adjust p-values and log-transform
  fact_res2 <- matrix(p.adjust(fact_res,method='fdr'),ncol=ncol(fact_res),nrow=nrow(fact_res))
  colnames(fact_res2) <- colnames(fact_res)
  fact_res2[is.na(fact_res2)] <- 1
  fact_res2 <- fact_res2[,colSums(fact_res2<sig_thresh)>0,drop=FALSE]
  fact_res2 <- -log10(fact_res2)
  rownames(fact_res2) <- rownames(myres_mat)
  
  # limit both to only those that are significant in the selected factor
  res_keep <- rownames(fact_res2)[fact_res2[,paste0('Factor_',factor_select)]>(-log10(sig_thresh))]
  myres_mat_sub <- myres_mat[res_keep,]
  myres_mat_sub <- myres_mat_sub[,colSums(myres_mat_sub>-log10(sig_thresh),na.rm = TRUE)>0]
  fact_res2_sub <- fact_res2[res_keep,]
  
  return(myres_mat_sub)
  
}

# for the ligands associated with a factor of interest, get the ligand_target cell type
# channels where we inferred significant interactions
parse_scITD_LR_out <- function(myres_mat,sig_thresh) {
  ## get unique ligand_source_receptor_target channels from my analysis
  unique_channels <- c()
  for (i in 1:nrow(myres_mat)) {
    lig_ct_rec <- myres_mat[i,]
    lig_ct_rec_name <- strsplit(rownames(myres_mat)[i],split='_')[[1]]
    lig <- lig_ct_rec_name[[1]]
    source <- lig_ct_rec_name[[2]]
    rec <- lig_ct_rec_name[[3]]
    
    for (j in 1:ncol(myres_mat)) {
      pv <- lig_ct_rec[j]
      if (!is.na(pv)) {
        if (pv > (-log10(sig_thresh))) {
          target_ct <- strsplit(names(pv),split="_")[[1]][[1]]
          lig_source_rec_target <- paste0(lig,"_",target_ct)
          unique_channels <- c(unique_channels,lig_source_rec_target)
        }
      }
    }
  }
  scITD_channels <- unique(unique_channels)
  return(scITD_channels)
}





# get the top n ligand_target cell type interactions inferred by nichenet
parse_nnet_out <- function(res_lst,n_channels) {
  # concat all data together
  total_dat <- res_lst[[1]][[1]]
  for (i in 2:length(res_lst)) {
    tmp <- res_lst[[i]][[1]]
    total_dat <- rbind(total_dat,tmp)
  }
  
  # sort by highest prioritization_score
  total_dat <- total_dat[order(total_dat$prioritization_score,decreasing = TRUE),]
  
  # remove LR pairs not in bonafide
  total_dat_sub <- total_dat[total_dat$bonafide,]
  
  total_dat_sub$sender_ct <- sapply(total_dat_sub$sender,function(x) {
    strsplit(x,split="_")[[1]][[1]]
  })
  total_dat_sub$receiver_ct <- sapply(total_dat_sub$receiver,function(x) {
    strsplit(x,split="_")[[1]][[1]]
  })
  total_dat_sub_channels <- sapply(1:nrow(total_dat_sub),function(x) {
    paste0(total_dat_sub[x,'ligand'],"_",total_dat_sub[x,'receiver_ct'])
  })
  
  total_dat_sub_channels <- unique(total_dat_sub_channels)
  
  total_dat_sub_channels <- total_dat_sub_channels[1:n_channels]
  
  return(total_dat_sub_channels)
}


### pipeline for validation
# get significant factor associated genes (by LM) by testing all genes (not just those in tensor)
# get list of significant DE genes for each experiment
# loop through factors
# loop through methods
# for each lig_target hit compute fgsea enrichment using factor associated genes for correct cell type and experiment data for proper ligand
# do this for each available LR dataset and take the median result or bonferroni best result
# store results for each method as ligand_target_factor
# for each dataset, compute fraction of validation tests that were significant

get_sig_genes <- function(container, all_f_test, donor_min_cells=2,
                          norm_method='trim', scale_factor=10000,
                          batch_var = 'Batch', adjust_p=TRUE) {
  # Need to run a modified script to get pseudobulk expression without scaling, so can remove lowly expressed genes
  container <- parse_data_by_ctypes(container)
  container <- clean_data(container, donor_min_cells=donor_min_cells)
  container <- get_pseudobulk(container)
  container <- normalize_pseudobulk(container, method=norm_method, scale_factor=scale_factor)
  
  # remove lowly expressed genes here
  for (ct in container$experiment_params$ctypes_use) {
    # select only genes expressed to some amount in at least 5% of donors
    donor_thresh <- round(nrow(container$scMinimal_ctype[[ct]]$pseudobulk) * .05)
    g_keep <- colSums(container$scMinimal_ctype[[ct]]$pseudobulk>0) > donor_thresh
    container$scMinimal_ctype[[ct]]$pseudobulk <- container$scMinimal_ctype[[ct]]$pseudobulk[,g_keep]
  }
  
  container <- apply_combat(container,batch_var=batch_var)
  
  # unit scale each gene in each cell type
  for (ct in container$experiment_params$ctypes_use) {
    container$scMinimal_ctype[[ct]]$pseudobulk <- scale(container$scMinimal_ctype[[ct]]$pseudobulk)
  }
  
  # need to unit scale donor scores for each factor, so can compare betas
  container$tucker_results[[1]] <- scale(container$tucker_results[[1]])
  
  ## now compute expression-dsc associations for all factors and all cell types
  total_res <- list() # stores the results for all factors
  for (f_test in all_f_test) {
    print(f_test)
    dsc <- container$tucker_results[[1]][,f_test,drop=FALSE]
    
    f_res <- list() # stores the result from a single factor
    
    # loop through cell types
    for (ct in container$experiment_params$ctypes_use) {
      print(ct)
      ct_res <- c()
      pb <- container$scMinimal_ctype[[ct]]$pseudobulk
      # loop through genes
      for (g_ndx in 1:ncol(pb)) {
        tmp <- cbind.data.frame(dsc,pb[rownames(dsc),g_ndx])
        colnames(tmp) <- c('dscore','expr')
        x <- summary(lm(expr~dscore,data=tmp))
        beta <- x$coefficients['dscore','Estimate']
        ct_res <- c(ct_res,beta**2)
        
        # pval <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)
        # ct_res <- c(ct_res,pval)
        
        
        # # trying to use what I did in GSEA previously (val = sum(expression * donor score))
        # exp_transform <- tmp[,1] * tmp[,2]
        # de_val <- sum(exp_transform)
        
      }
      names(ct_res) <- colnames(pb)
      
      if (adjust_p) {
        ct_res <- p.adjust(ct_res,method='fdr')
      }
      f_res[[ct]] <- ct_res
    }
    total_res[[f_test]] <- f_res
  }
  
  return(total_res)
}


# helper function to keep code modular
# tests all significant channels for enrichment in experiment data where available
test_channels <- function(my_channels,sig_genes_factor,val_data_proc) {
  chan_res <- list()
  for (channel in my_channels) {
    lig <- strsplit(channel,split='_')[[1]][[1]]
    target <- strsplit(channel,split='_')[[1]][[2]]
    
    # only proceed if ligand is in the experiment data
    if (lig %in% names(val_data_proc)) {
      # get the significant genes for target ct of the specified factor
      assoc_padj <- sig_genes_factor[[target]]
      
      # get the ligand experiment sig genes for the ligand
      my_pathways <- val_data_proc[[lig]]
      
      # # run fgsea for each dataset sig genes
      # mystats <- -log10(assoc_padj)
      
      mystats <- assoc_padj
      
      # fgsea_res <- fgsea::fgsea(pathways = my_pathways,
      #                           stats = mystats,
      #                           minSize=0,
      #                           maxSize=1500,
      #                           eps=0,
      #                           gseaParam=0, # unweighted
      #                           scoreType = "pos",
      #                           nproc=4)
      
      fgsea_res <- fgsea::fgsea(pathways = my_pathways,
                                stats = mystats,
                                minSize=0,
                                maxSize=1500,
                                eps=0,
                                gseaParam=1, # weighted
                                scoreType = "pos",
                                nproc=4)
      
      val_pval <- min(fgsea_res$pval) * nrow(fgsea_res)
      chan_res[[channel]] <- min(val_pval,1)
    }
  }
  return(chan_res)
}


sig_genes <- get_sig_genes(pbmc_container, all_f_test=c(1,2,3), donor_min_cells=20,
                           norm_method='trim', scale_factor=10000,
                           batch_var='pool', adjust_p = FALSE)

# creating list of significant genes per ligand experiment
val_data = readRDS(url("https://zenodo.org/record/3260758/files/expression_settings.rds"))
# saveRDS(val_data2,file='/home/jmitchel/data/LR_datasets/NicheNet_validation_dataset.rds')
# val_data <- readRDS(file='/home/jmitchel/data/LR_datasets/NicheNet_validation_dataset.rds')
val_data_proc <- list()
for (i in 1:length(val_data)) {
  val_data_inst <- val_data[[i]]
  # only using ligand treatment datasets where 1 ligand was used
  if (length(val_data_inst$from)==1) {
    # using same thresholds as they use in the paper
    ndx_keep1 <- which(abs(val_data_inst$diffexp$lfc)>1)
    ndx_keep2 <- which(val_data_inst$diffexp$qval<.1)
    ndx_keep <- intersect(ndx_keep1,ndx_keep2)
    val_gene_sig <- val_data_inst$diffexp$gene[ndx_keep]
    val_data_proc[[val_data_inst$from]][[val_data_inst$name]] <- val_gene_sig
  }
}

## looping through factors
nnet_res_all <- list(res_lst_f1,res_lst_f2,res_lst_f3)
# sig_threshold <- .1 # used originally to get the most results
sig_threshold <- .05 # used originally to get the most results
final_res <- list()
for (f_test in 1:3) {
  print(f_test)
  for (method_test in c('scITD','nnet','rand_scITD','rand_nnet')) {
    print(method_test)
    
    if (method_test=='scITD') {
      myres_mat <- get_lr_dat(pbmc_container,sig_thresh=sig_threshold,factor_select=f_test)
      my_channels <- parse_scITD_LR_out(myres_mat,sig_thresh=sig_threshold)
      # extract all gene associations with the factor
      sig_genes_factor <- sig_genes[[f_test]]
      # run unsigned enrichment tests for each channel using validation data
      chan_res <- test_channels(my_channels,sig_genes_factor,val_data_proc)
      final_res[[method_test]][[f_test]] <- chan_res
      
    } else if (method_test=='nnet') {
      myres_mat <- get_lr_dat(pbmc_container,sig_thresh=sig_threshold,factor_select=f_test)
      my_channels <- parse_scITD_LR_out(myres_mat,sig_thresh=sig_threshold)
      my_channels <- parse_nnet_out(nnet_res_all[[f_test]],n_channels=length(my_channels))
      # extract all gene associations with the factor
      sig_genes_factor <- sig_genes[[f_test]]
      # run unsigned enrichment tests for each channel using validation data
      chan_res <- test_channels(my_channels,sig_genes_factor,val_data_proc)
      final_res[[method_test]][[f_test]] <- chan_res
      
    } else if (method_test=='rand_scITD' || method_test=='rand_nnet') {
      all_lig_ct_combos <- sapply(names(val_data_proc),function(x){
        sapply(names(nnet_res_all[[1]]),function(y){
          paste0(x,'_',y)
        })
      })
      all_lig_ct_combos <- c(all_lig_ct_combos)
      
      for (i in 1:100) {
        if (method_test=='rand_scITD') {
          my_channels_samp <- sample(all_lig_ct_combos,length(final_res[['scITD']][[f_test]]))
        } else if (method_test=='rand_nnet') {
          my_channels_samp <- sample(all_lig_ct_combos,length(final_res[['nnet']][[f_test]]))
        }
        # extract all gene associations with the factor
        sig_genes_factor <- sig_genes[[f_test]]
        # run unsigned enrichment tests for each channel using validation data
        chan_res <- test_channels(my_channels_samp,sig_genes_factor,val_data_proc)
        if (length(final_res[[method_test]])<i) {
          final_res[[method_test]][[i]] <- list()
        }
        final_res[[method_test]][[i]][[f_test]] <- chan_res
      }
    }
  }
}

# saveRDS(final_res,'/home/jmitchel/data/lupus_data/scITD_nichenet_compare_validation_res.rds')




fr1 <- unlist(final_res[['scITD']])
fr2 <- unlist(final_res[['nnet']])

fr1 <- p.adjust(fr1,method='fdr')
fr2 <- p.adjust(fr2,method='fdr')

p_thresh <- .001

# for rand_scITD
rand_sc_counts <- c()
rand_sc_fracs <- c()
for (i in 1:length(final_res[['rand_scITD']])) {
  iter_res <- unlist(final_res[['rand_scITD']][[i]])
  
  # p correct
  iter_res <- p.adjust(iter_res,method='fdr')
  
  # compute how many less than p_thresh
  count1 <- sum(iter_res<p_thresh)
  frac1 <- count1/length(iter_res)
  rand_sc_counts <- c(rand_sc_counts,count1)
  rand_sc_fracs <- c(rand_sc_fracs,frac1)
}

# for rand_nnet
rand_nn_counts <- c()
rand_nn_fracs <- c()
for (i in 1:length(final_res[['rand_nnet']])) {
  iter_res <- unlist(final_res[['rand_nnet']][[i]])
  
  # p correct
  iter_res <- p.adjust(iter_res,method='fdr')
  
  # compute how many less than p_thresh
  count1 <- sum(iter_res<p_thresh)
  frac1 <- count1/length(iter_res)
  rand_nn_counts <- c(rand_nn_counts,count1)
  rand_nn_fracs <- c(rand_nn_fracs,frac1)
}

r_nn_av_count <- mean(rand_nn_counts)
r_nn_av_frac <- mean(rand_nn_fracs)
r_nn_sd_count <- sd(rand_nn_counts)
r_nn_sd_count <- r_nn_sd_count/sqrt(length(rand_nn_counts))

r_sc_av_count <- mean(rand_sc_counts)
r_sc_av_frac <- mean(rand_sc_fracs)
r_sc_sd_count <- sd(rand_sc_counts)
r_sc_sd_count <- r_sc_sd_count/sqrt(length(rand_sc_counts))

count1 <- sum(fr1<p_thresh)
count2 <- sum(fr2<p_thresh)

frac1 <- count1/length(fr1)
frac2 <- count2/length(fr2)

tmp <- cbind.data.frame(c('scITD','nnet','scITD','nnet','rand_scITD','rand_nnet','rand_scITD','rand_nnet'),
                        c(count1,count2,length(fr1)-count1,length(fr2)-count2,r_sc_av_count,r_nn_av_count,length(fr1)-r_sc_av_count,length(fr2)-r_nn_av_count),
                        c('pos','pos','neg','neg','pos','pos','neg','neg'),
                        c(frac1,frac2,1-frac1,1-frac2,r_sc_av_frac,r_nn_av_frac,1-r_sc_av_frac,1-r_nn_av_frac),
                        c(NA,NA,NA,NA,r_sc_sd_count,r_nn_sd_count,NA,NA))
colnames(tmp) <- c('method','count','type','frac','count_sd')
tmp$frac <- round(tmp$frac,digits=2)

tmp$method <- factor(tmp$method,levels = c('nnet','rand_nnet','scITD','rand_scITD'))

p <- ggplot(tmp, aes(x = factor(method), y = count, fill = type)) +
  geom_bar(position = position_stack(), stat = "identity", width = .7) +
  geom_errorbar(aes(ymin = count-count_sd, ymax = count+count_sd),
                color="black", 
                width=0.1, position='identity') +
  geom_text(aes(label = frac), position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_manual(breaks = c("pos", "neg"), values=c("light blue", "dark gray")) +
  theme_bw()
p

### Figure S4G
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs2/scitd_vs_nnet_validation_v2.pdf", useDingbats = FALSE,
#     width = 5, height = 3.5) # using betase, .05 thresh, and weighted gsea
p
# dev.off()
# using scITD p_threshold of .05 with final p_threshold of .001 with betas and weighted version



length(intersect(names(final_res[['scITD']][[1]]),names(final_res[['nnet']][[1]])))
length(intersect(names(final_res[['scITD']][[2]]),names(final_res[['nnet']][[2]])))
length(intersect(names(final_res[['scITD']][[3]]),names(final_res[['nnet']][[3]])))

# only 4 channels were predicted by both methods











##### making comparison plots
#### trying to create a df of all top hits to sender ctypes for the "high" or "low" categories
res_lst <- readRDS(file='/home/jmitchel/data/lupus_data/nichenet_f2.rds')
all_ctypes <- names(res_lst)
prioritized_tbl_oi_all <- NULL
get_circos_plot <- function(res_lst,high_low) {
  for (ct in all_ctypes) {
    ct_other <- all_ctypes[all_ctypes!=ct]
    
    senders_high <- sapply(ct_other,function(x){
      paste0(x,"_high")
    })
    names(senders_high) <- NULL
    senders_low <- sapply(ct_other,function(x){
      paste0(x,"_low")
    })
    names(senders_low) <- NULL
    
    niches = list(
      "High_niche" = list(
        "sender" = senders_high,
        "receiver" = paste0(ct,"_high")),
      "Low_niche" = list(
        "sender" = senders_low,
        "receiver" = paste0(ct,"_low"))
    ) # user adaptation required on own dataset
    
    prioritization_tables <- res_lst[[ct]]
    
    top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
    top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% group_by(ligand, receptor) %>% top_n(1, prioritization_score) %>% ungroup() %>% select(ligand, receptor, niche) %>% rename(top_niche = niche)
    
    ligand_prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% select(niche, sender, receiver, ligand, prioritization_score) %>% group_by(ligand, niche) %>% top_n(1, prioritization_score) %>% ungroup() %>% distinct() %>% inner_join(top_ligand_niche_df) %>% filter(niche == top_niche) %>% group_by(niche) %>% top_n(50, prioritization_score) %>% ungroup() # get the top50 ligands per niche
    
    receiver_oi = paste0(ct,"_",high_low) 
    
    # just using top 5 ligands because likely going to be a lot
    filtered_ligands = ligand_prioritized_tbl_oi %>% filter(receiver == receiver_oi) %>% top_n(5, prioritization_score) %>% pull(ligand) %>% unique()
    
    prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>% filter(ligand %in% filtered_ligands) %>% select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% distinct() %>% inner_join(top_ligand_receptor_niche_df) %>% group_by(ligand) %>% filter(receiver == receiver_oi) %>% top_n(1, prioritization_score) %>% ungroup() 
    
    if (is.null(prioritized_tbl_oi_all)) {
      prioritized_tbl_oi_all <- prioritized_tbl_oi
    } else {
      prioritized_tbl_oi_all <- rbind(prioritized_tbl_oi_all,prioritized_tbl_oi)
    }
  }
  
  # dimensions of final table should be 5*2*7 = 70
  dim(prioritized_tbl_oi_all)
  
  # removing the CD99--CD99 ligand-receptor because it causes plot to error
  prioritized_tbl_oi_all <- prioritized_tbl_oi_all[prioritized_tbl_oi_all$ligand_receptor!="CD99--CD99",]
  
  colors_receiver = brewer.pal(n = prioritized_tbl_oi_all$receiver %>% unique() %>% sort() %>% length(), name = 'Spectral') %>% magrittr::set_names(prioritized_tbl_oi_all$receiver %>% unique() %>% sort())
  colors_sender = colors_receiver[unique(prioritized_tbl_oi_all$sender)]
  
  circos_output = make_circos_lr(prioritized_tbl_oi_all, colors_sender, colors_receiver)
  return(circos_output)
}

p1 <- get_circos_plot(res_lst,high_low='high')
p2 <- get_circos_plot(res_lst,high_low='low')

### Figure S4D right
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/F2_LR_high_circ.pdf", useDingbats = FALSE,
#     width = 8, height = 4)
p1
# dev.off()

### Figure S4D left
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/F2_LR_low_circ.pdf", useDingbats = FALSE,
#     width = 8, height = 4)
p2
dev.off()









#### now to do analysis on the full output table...
# should make sure to limit the results to just those that have TRUE for bonafide in the network

# concat all data together
total_dat <- res_lst[[1]][[1]]
for (i in 2:length(res_lst)) {
  tmp <- res_lst[[i]][[1]]
  total_dat <- rbind(total_dat,tmp)
}

# sort by highest prioritization_score
total_dat <- total_dat[order(total_dat$prioritization_score,decreasing = TRUE),]

# remove LR pairs not in bonafide
total_dat_sub <- total_dat[total_dat$bonafide,]

total_dat_sub$sender_ct <- sapply(total_dat_sub$sender,function(x) {
  strsplit(x,split="_")[[1]][[1]]
})
total_dat_sub$receiver_ct <- sapply(total_dat_sub$receiver,function(x) {
  strsplit(x,split="_")[[1]][[1]]
})
total_dat_sub_channels <- sapply(1:nrow(total_dat_sub),function(x) {
  paste0(total_dat_sub[x,'ligand'],"_",total_dat_sub[x,'sender_ct'],"_",
         total_dat_sub[x,'receptor'],"_",total_dat_sub[x,'receiver_ct'])
})
total_dat_sub_channels <- unique(total_dat_sub_channels)



col_select <- c('ligand','receptor','sender','receiver','prioritization_score')
n_select <- 10000
dat_sub <- total_dat_sub[1:n_select,col_select]
dat_sub$sender <- sapply(dat_sub$sender,function(x) {
  strsplit(x,split="_")[[1]][[1]]
})
dat_sub$receiver <- sapply(dat_sub$receiver,function(x) {
  strsplit(x,split="_")[[1]][[1]]
})
dat_sub_channels <- sapply(1:nrow(dat_sub),function(x) {
  paste0(dat_sub[x,'ligand'],"_",dat_sub[x,'sender'],"_",dat_sub[x,'receptor'],"_",dat_sub[x,'receiver'])
})
dat_sub_channels_unq <- unique(dat_sub_channels)
print(length(dat_sub_channels_unq))

## adjust channel extraction fn to give ligand_source_receptor_target
parse_scITD_LR_out_full <- function(myres_mat,sig_thresh) {
  ## get unique ligand_source_receptor_target channels from my analysis
  unique_channels <- c()
  for (i in 1:nrow(myres_mat)) {
    lig_ct_rec <- myres_mat[i,]
    lig_ct_rec_name <- strsplit(rownames(myres_mat)[i],split='_')[[1]]
    lig <- lig_ct_rec_name[[1]]
    source <- lig_ct_rec_name[[2]]
    rec <- lig_ct_rec_name[[3]]
    
    for (j in 1:ncol(myres_mat)) {
      pv <- lig_ct_rec[j]
      if (!is.na(pv)) {
        if (pv > (-log10(sig_thresh))) {
          target_ct <- strsplit(names(pv),split="_")[[1]][[1]]
          lig_source_rec_target <- paste0(lig,"_",source,"_",rec,"_",target_ct)
          unique_channels <- c(unique_channels,lig_source_rec_target)
        }
      }
    }
  }
  scITD_channels <- unique(unique_channels)
  return(scITD_channels)
}


res <- matrix(nrow=4,ncol=6)
colnames(res) <- c('padj_thresh','jaccard_coef','rand_jaccard_mean','rand_jaccard_sd','p_enr','num_res')
p_test <- c(.001,.01,.05,.1)
for (p in p_test) {
  lr_dat <- get_lr_dat(pbmc_container,sig_thresh=p,factor_select = 2)
  channels <- parse_scITD_LR_out_full(lr_dat,sig_thresh=p)
  dat_sub_channels_unq_top <- dat_sub_channels_unq[1:length(channels)]
  in_both <- intersect(dat_sub_channels_unq_top,channels)
  in_either <- union(dat_sub_channels_unq_top,channels)
  jaccard <- length(in_both)/length(in_either)
  
  # get jaccard from random channels
  rj_all <- sapply(1:10000,function(i) {
    ch_samp <- sample(total_dat_sub_channels,length(channels))
    in_both <- intersect(dat_sub_channels_unq_top,ch_samp)
    in_either <- union(dat_sub_channels_unq_top,ch_samp)
    rand_jaccard <- length(in_both)/length(in_either)
  })
  
  p_enr <- sum(rj_all>jaccard)/length(rj_all)
  r_ndx <- which(p_test==p)
  res[r_ndx,] <- c(p,jaccard,mean(rj_all),sd(rj_all),p_enr,length(channels))
}

res <- as.data.frame(res)
res1 <- res[,c('padj_thresh','jaccard_coef','rand_jaccard_sd')]
res1[,3] <- rep(0,4)
res2 <- res[,c('padj_thresh','rand_jaccard_mean','rand_jaccard_sd')]
colnames(res2) <- c('padj_thresh','jaccard_coef','rand_jaccard_sd')

res_transform <- rbind.data.frame(res1,res2)
res_transform$rand <- c(rep(FALSE,4),rep(TRUE,4))
res_transform$padj_thresh <- factor(res_transform$padj_thresh,
                                    levels=c(.1,.05,.01,.001))
p <- ggplot(res_transform,aes(x=padj_thresh,y=jaccard_coef,fill=rand)) +
  geom_bar(stat = 'identity',
           position = "dodge") +
  geom_errorbar(aes(ymin = jaccard_coef-rand_jaccard_sd, ymax = jaccard_coef+rand_jaccard_sd),
                color="black", 
                width=0.1, position=position_dodge(0.9)) +
  xlab('scITD padj threshold') +
  ylab('Jaccard coefficient') +
  theme_bw()
p

### Figure S4E
# pdf(file = "/home/jmitchel/figures/scITD_revision_figs/LR_jaccard.pdf", useDingbats = FALSE,
#     width = 4, height = 3)
p
# dev.off()






