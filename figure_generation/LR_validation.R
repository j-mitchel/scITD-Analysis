library(Seurat)
library(nichenetr)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(dplyr)
library(sccore)
library(devtools)
load_all('/home/jmitchel/scITD/')


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


pbmc_container <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_container_w_decomp.rds')


# function to get nichenet results for a factor
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
# ligand_target_matrix = readRDS("/home/jmitchel/data/NicheNet/ligand_target_matrix.rds")
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
# saveRDS(lr_network,file='/home/jmitchel/data/NicheNet/lr_network.rds')
# lr_network = readRDS("/home/jmitchel/data/NicheNet/lr_network.rds")
lr_network = lr_network %>% mutate(bonafide = ! database %in% c("ppi_prediction","ppi_prediction_go"))
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor, bonafide)




#### running scITD LR analysis using nichenet LR pairs
lr_pairs <- as.matrix(lr_network)
lr_pairs <- lr_pairs[lr_pairs[,'bonafide']=='TRUE',]
lr_pairs <- lr_pairs[,c(1,2)]

# infer active LR interactions
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='pool')
scITD_scale_pb_nnet <- pbmc_container$scale_pb_extra
# saveRDS(scITD_scale_pb_nnet,file='/home/jmitchel/data/lupus_data/scITD_nnet_pb_scaled.rds')
scITD_scale_pb_nnet <- readRDS(file='/home/jmitchel/data/lupus_data/scITD_nnet_pb_scaled.rds')

sft_thresh <- c(12,14,12,10,12,9,12)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

lr_hmap <- compute_LR_interact(pbmc_container, lr_pairs, sig_thresh=.00000000005,
                               percentile_exp_rec=0.85, add_ld_fact_sig=TRUE)

# saveRDS(pbmc_container$lr_res,file='/home/jmitchel/data/lupus_data/nnet_scITD_lr_res.rds')
scITD_nnet_dat <- pbmc_container$lr_res

scITD_nnet_dat <- readRDS(file='/home/jmitchel/data/lupus_data/nnet_scITD_lr_res.rds')





## function to run nichenet across a discretized scITD factor
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

# get the top n ligand_target cell type interactions inferred by nichenet for a given factor
parse_nnet_out <- function(res_lst,n_channels,channel_type,return_all_channels=FALSE) {
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
  
  if (channel_type=='ligand_rctype') {
    total_dat_sub_channels <- sapply(1:nrow(total_dat_sub),function(x) {
      paste0(total_dat_sub[x,'ligand'],"_",total_dat_sub[x,'receiver_ct'])
    })
  } else if (channel_type=='full') {
    total_dat_sub_channels <- sapply(1:nrow(total_dat_sub),function(x) {
      paste0(total_dat_sub[x,'ligand'],"_",total_dat_sub[x,'sender_ct'],"_",
             total_dat_sub[x,'receptor'],"_",total_dat_sub[x,'receiver_ct'])
    })
    
  }
  
  total_dat_sub_channels <- unique(total_dat_sub_channels)
  
  if (!return_all_channels) {
    total_dat_sub_channels <- total_dat_sub_channels[1:n_channels]
  }
  
  return(total_dat_sub_channels)
}

n_channels <- 400
nnet_f1 <- parse_nnet_out(res_lst_f1,n_channels=n_channels,channel_type='full')
nnet_f2 <- parse_nnet_out(res_lst_f2,n_channels=n_channels,channel_type='full')
nnet_f3 <- parse_nnet_out(res_lst_f3,n_channels=n_channels,channel_type='full')
nnet_sig_channels_full <- unique(c(nnet_f1,nnet_f2,nnet_f3))
print(length(nnet_sig_channels_full))


# selecting the same number of top pairs from scITD
get_matching_scITD_pairs <- function(lr_res,length_full) {
  mat <- lr_res
  mat_df <- as.data.frame(as.table(mat))
  colnames(mat_df) <- c("Row", "Column", "Value")
  
  # Sort the data frame by 'Value' column to get smallest elements
  sorted_df <- mat_df[order(mat_df$Value), ]
  
  # Select top N smallest elements
  top_N_smallest <- sorted_df
  top_N_smallest$Column <- as.character(top_N_smallest$Column)
  top_N_smallest$Row <- as.character(top_N_smallest$Row)
  top_N_smallest$rec_ct <- sapply(top_N_smallest$Column,function(x){
    strsplit(x,split='_')[[1]][[1]]
  })
  top_N_smallest$ligand <- sapply(top_N_smallest$Row,function(x){
    strsplit(x,split='_')[[1]][[1]]
  })
  top_N_smallest$channels_full <- paste0(top_N_smallest$Row,'_',top_N_smallest$rec_ct)
  top_N_smallest$channels_short <- paste0(top_N_smallest$ligand,'_',top_N_smallest$rec_ct)
  
  all_full <- unique(top_N_smallest$channels_full)
  all_short <- unique(top_N_smallest$channels_short)
  top_full <- all_full[1:length_full]
  # top_short <- all_short[1:length_short]
  return(top_full)
}
scITD_nnet_full <- get_matching_scITD_pairs(scITD_nnet_dat,length(nnet_sig_channels_full))



# loading up Tensor-C2C inferred channels
### loading cell2cell results
lr_ldngs <- read.csv('/home/jmitchel/data/lupus_data/tensor_c2c_lr_lds.csv',row.names = 1)
sample_ldngs <- read.csv('/home/jmitchel/data/lupus_data/tensor_c2c_donor_lds.csv',row.names = 1)
sender_ldngs <- read.csv('/home/jmitchel/data/lupus_data/tensor_c2c_ct_send_lds.csv',row.names = 1)
receiver_ldngs <- read.csv('/home/jmitchel/data/lupus_data/tensor_c2c_ct_rec_lds.csv',row.names = 1)
liana_res <- read.csv('/home/jmitchel/data/lupus_data/tensor_c2c_LR_tests.csv')

rownames(sample_ldngs) <- sapply(rownames(sample_ldngs),function(x){
  paste(strsplit(x,split=':')[[1]],collapse='')
})


### need to extract the significant channels associated with any factor in cell2cell
top_n_select <- 100
sig_thresh <- .05
sig_channels <- c()
for (f_ndx in 1:ncol(lr_ldngs)) {
  lds_fact <- lr_ldngs[,f_ndx,drop=FALSE]
  lds_fact <- lds_fact[order(lds_fact[,1],decreasing=TRUE),,drop=FALSE]
  top_lr <- rownames(lds_fact)[1:top_n_select]
  
  # get top from and to cell types, binarizing by one sd over median
  send_fact <- sender_ldngs[,f_ndx,drop=FALSE]
  rec_fact <- receiver_ldngs[,f_ndx,drop=FALSE]
  
  sd_send <- sd(send_fact[,1])
  sd_rec <- sd(rec_fact[,1])
  
  send_ct <- rownames(send_fact)[send_fact[,1]>(median(send_fact[,1]) + sd_send)]
  rec_ct <- rownames(rec_fact)[rec_fact[,1]>(median(rec_fact[,1]) + sd_rec)]
  
  if (length(send_ct)==0 | length(rec_ct)==0) {
    next
  }
  
  ## now selecting the subset of the lr ct combinations that were originally significant
  ct_combos <- expand.grid(send_ct, rec_ct)
  ct_combos[,1] <- as.character(ct_combos[,1])
  ct_combos[,2] <- as.character(ct_combos[,2])
  ct_combos <- ct_combos[ct_combos[,1] != ct_combos[,2],]
  
  for (lr in top_lr) {
    lig_rec <- strsplit(lr,split='^',fixed = TRUE)[[1]]
    lig <- lig_rec[[1]]
    rec <- lig_rec[[2]]
    for (ct_ndx in 1:nrow(ct_combos)) {
      ct_combo <- ct_combos[ct_ndx,]
      source_ct <- as.character(ct_combo[1])
      rec_ct <- as.character(ct_combo[2])
      lr_res_all_donors <- liana_res[liana_res$ligand_complex==lig & liana_res$receptor_complex==rec & liana_res$source==source_ct & liana_res$target==rec_ct,]
      # if significant in any donors, keep the interaction
      if (any(lr_res_all_donors$cellphone_pvals<sig_thresh)) {
        cur_channel <- paste(lig,source_ct,rec,rec_ct,sep = '_')
        sig_channels <- c(sig_channels,cur_channel)
      }
    }
  }
}
c2c_sig_channels_full <- unique(sig_channels)
print(length(c2c_sig_channels_full))


# loadings up scITD results run using tensor c2c lr pairs
# load up the lr pairs tested for with tensor cell2cell
lr_pairs <- read.csv('/home/jmitchel/data/lupus_data/tensor_c2c_LR_pairs.csv',row.names = 1)
rownames(lr_pairs) <- NULL
colnames(lr_pairs) <- c('ligand','receptor')

# infer active LR interactions
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='pool')
scITD_scale_pb_c2c <- pbmc_container$scale_pb_extra
# saveRDS(scITD_scale_pb_c2c,file='/home/jmitchel/data/lupus_data/scITD_c2c_pb_scaled.rds')
scITD_scale_pb_c2c <- readRDS(file='/home/jmitchel/data/lupus_data/scITD_c2c_pb_scaled.rds')

sft_thresh <- c(12,14,12,10,12,9,12)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

lr_hmap <- compute_LR_interact(pbmc_container, lr_pairs, sig_thresh=.05,
                               percentile_exp_rec=0.85, add_ld_fact_sig=TRUE)

# saveRDS(pbmc_container$lr_res,file='/home/jmitchel/data/lupus_data/C2C_scITD_lr_res.rds')
scITD_c2c_dat <- readRDS(file='/home/jmitchel/data/lupus_data/C2C_scITD_lr_res.rds')

# parsing scITD results with c2c lr pairs
scITD_c2c_full <- get_matching_scITD_pairs(scITD_c2c_dat,length(c2c_sig_channels_full))

# loading scTensor results from /scTensor_test.R file
gene_conv <- readRDS('/home/jmitchel/data/lupus_data/sle_genes_ncbi.rds') # to convert id's to gene symbols
sc_tensor_full <- readRDS(file='/home/jmitchel/data/lupus_data/scTensor_res.rds')
l_ct_f <- sc_tensor_full[[1]]
r_ct_f <- sc_tensor_full[[2]]
core_tensor <- sc_tensor_full[[3]]

# parse scTensor results into channels with just ligands and receptor cell types
# first, binarize scTensor factors
for (i in 1:nrow(l_ct_f)) {
  med_val <- median(l_ct_f[i,])
  sd_val <- sd(l_ct_f[i,])
  ct_rem_ndx <- which(l_ct_f[i,]<(med_val+sd_val))
  l_ct_f[i,ct_rem_ndx] <- 0
}

for (i in 1:nrow(r_ct_f)) {
  med_val <- median(r_ct_f[i,])
  sd_val <- sd(r_ct_f[i,])
  ct_rem_ndx <- which(r_ct_f[i,]<(med_val+sd_val))
  r_ct_f[i,ct_rem_ndx] <- 0
}

# next select the top 200 ligands with largest max values in core tensor
lr_pairs_in_tensor <- dimnames(core_tensor@data)[[3]]
lr_max_core <- c()
for (i in 1:length(lr_pairs_in_tensor)) {
  max_val <- max(core_tensor@data[,,i])
  lr_max_core <- c(lr_max_core,max_val)
}
ndx_top <- order(lr_max_core,decreasing=TRUE)[1:400] #indices of top 400 lr pairs


# for each lr pair, select biggest factor interaction from core tensor
sig_channels <- c()
for (i in ndx_top) {
  pair_current <- lr_pairs_in_tensor[i]
  
  lig <- strsplit(pair_current,split='_')[[1]][[1]]
  rec <- strsplit(pair_current,split='_')[[1]][[2]]
  
  # convert lig and rec back to gene symbols
  ndx_gene <- which(gene_conv[,'NCBI gene (formerly Entrezgene) ID']==lig)
  lig <- gene_conv[ndx_gene,'HGNC symbol']
  ndx_gene <- which(gene_conv[,'NCBI gene (formerly Entrezgene) ID']==rec)
  rec <- gene_conv[ndx_gene,'HGNC symbol']
  
  ndx_best <- which(core_tensor@data[,,i]==max(core_tensor@data[,,i]),arr.ind = TRUE) # indicates lig_factor, rec_factor
  
  ct_ndx <- which(l_ct_f[ndx_best[1],] != 0)
  lig_ct <- colnames(l_ct_f)[ct_ndx]
  
  ct_ndx <- which(r_ct_f[ndx_best[2],] != 0)
  rec_ct <- colnames(r_ct_f)[ct_ndx]
  
  # get combinations of ligands and receptors
  ct_combos <- expand.grid(lig_ct, rec_ct)
  ct_combos[,1] <- as.character(ct_combos[,1])
  ct_combos[,2] <- as.character(ct_combos[,2])
  ct_combos <- ct_combos[ct_combos[,1] != ct_combos[,2],]
  
  for (ct_ndx in 1:nrow(ct_combos)) {
    ct_combo <- ct_combos[ct_ndx,]
    source_ct <- as.character(ct_combo[1])
    rec_ct <- as.character(ct_combo[2])
    cur_channel <- paste(lig,source_ct,rec,rec_ct,sep = '_')
    sig_channels <- c(sig_channels,cur_channel)
  }
}

scT_sig_channels_full <- unique(sig_channels)
print(length(scT_sig_channels_full))


## loading precomputed scITD results using scTensor LR pairs
scITD_scT_dat <- readRDS(file='/home/jmitchel/data/lupus_data/scTensor_scITD_lr_res.rds')

scITD_scT_full <- get_matching_scITD_pairs(scITD_scT_dat,length(scT_sig_channels_full))







# creating list of significant genes per ligand experiment
val_data = readRDS(url("https://zenodo.org/record/3260758/files/expression_settings.rds"))
# saveRDS(val_data2,file='/home/jmitchel/data/LR_datasets/NicheNet_validation_dataset.rds')
# val_data <- readRDS(file='/home/jmitchel/data/LR_datasets/NicheNet_validation_dataset.rds')
val_data_proc <- list()
for (i in 1:length(val_data)) {
  val_data_inst <- val_data[[i]]
  # only using ligand treatment datasets where 1 ligand was used
  if (length(val_data_inst$from)==1) {
    val_data_proc[[val_data_inst$from]][[val_data_inst$name]] <- as.data.frame(val_data_inst$diffexp)
  }
}



scITD_scale_pb_scT <- readRDS(file='/home/jmitchel/data/lupus_data/scITD_scTensor_pb_scaled.rds')


get_target_genes <- function(full_channel,top_n_targets,expr_use) {
  chan_splt <- strsplit(full_channel,split='_')[[1]]
  lig <- chan_splt[[1]]
  lig_ct <- chan_splt[[2]]
  rec <- chan_splt[[3]]
  rec_ct <- chan_splt[[length(chan_splt)]]
  
  pb <- expr_use[[lig_ct]]
  lig_expr <- pb[,lig]
  
  pb <- expr_use[[rec_ct]]
  all_stats <- c()
  for (j in 1:ncol(pb)) {
    d_both <- intersect(names(lig_expr),rownames(pb))
    # don't test if the gene is the ligand
    if (colnames(pb)[j]==lig) {
      next
    }
    gene_expr <- pb[,j]
    tmp <- cbind.data.frame(lig_expr[d_both],gene_expr[d_both])
    colnames(tmp) <- c('ligand','target')
    lmres <- summary(lm(target~ligand,data=tmp))
    mystat <- lmres$r.squared
    # mystat <- abs(lmres$coefficients['ligand','t value'])
    names(mystat) <- colnames(pb)[j]
    all_stats <- c(all_stats,mystat)
  }
  all_stats <- all_stats[order(all_stats,decreasing=TRUE)][1:top_n_targets]
  return(all_stats)
}

enrichment_test <- function(channel_vec,top_n_targets,expr_use) {
  channel_enr_res <- list()
  
  # first get list of channels in the lig de datasets
  channels_test <- c()
  for (full_channel in channel_vec) {
    chan_splt <- strsplit(full_channel,split='_')[[1]]
    lig <- chan_splt[[1]]
    lig_ct <- chan_splt[[2]]
    rec <- chan_splt[[3]]
    rec_ct <- chan_splt[[4]]
    if (lig %in% names(val_data_proc)) {
      channels_test <- c(channels_test,full_channel)
    }
  }
  enr_pvals <- plapply(channels_test,function(full_channel) tryCatch({
    chan_splt <- strsplit(full_channel,split='_')[[1]]
    lig <- chan_splt[[1]]
    lig_ct <- chan_splt[[2]]
    rec <- chan_splt[[3]]
    rec_ct <- chan_splt[[length(chan_splt)]]
    
    # get top targets for the ligand in the target ct
    top_targets <- get_target_genes(full_channel, top_n_targets = top_n_targets, expr_use=expr_use)
    
    # get the ligand experiment sig genes for the ligand
    deg_results_all <- val_data_proc[[lig]]
    
    deg_enr_all <- c()
    for (i in 1:length(deg_results_all)) {
      print(i)
      deg_results <- deg_results_all[[i]]
      
      # to compute enrichment test using hypergeometric test
      g_both <- intersect(deg_results$gene,colnames(expr_use[[rec_ct]]))
      deg_results <- deg_results[deg_results$gene %in% g_both,]
      # using same thresholds as they use in the paper
      ndx_keep1 <- which(abs(deg_results$lfc)>1)
      ndx_keep2 <- which(deg_results$qval<.1)
      ndx_keep <- intersect(ndx_keep1,ndx_keep2)
      gene_sig <- deg_results$gene[ndx_keep]

      top_targets <- names(top_targets)[names(top_targets) %in% g_both]
      num_targets_sig <- sum(top_targets %in% gene_sig)
      pval <- phyper(num_targets_sig-1, length(top_targets), length(g_both)-length(top_targets), length(gene_sig), lower.tail=FALSE)
      
      deg_enr_all <- c(deg_enr_all,pval)
    }
    # take bonferroni corrected best pvalue
    pval2 <- min(deg_enr_all) * length(deg_enr_all)
    pval2 <- pmin(1,pval2)
    return(pval2)
  },error=function(e) paste0('error ',full_channel)), mc.preschedule=TRUE,n.cores=30,progress=TRUE)
  names(enr_pvals) <- channels_test
  return(enr_pvals)
}


## selecting random channels out of the ones tested for each dataset

get_rand_channels <- function(dat,expr_use) {
  mat_df <- as.data.frame(as.table(dat))
  colnames(mat_df) <- c("Row", "Column", "Value")
  mat_df$rec_ct <- sapply(as.character(mat_df$Column),function(x){
    strsplit(x,split='_')[[1]][[1]]
  })
  all_channels <- paste0(mat_df$Row,'_',mat_df$rec_ct)
  all_channels <- unique(all_channels)
  
  all_channels_in_val <- c()
  for (full_channel in all_channels) {
    chan_splt <- strsplit(full_channel,split='_')[[1]]
    lig <- chan_splt[[1]]
    lig_ct <- chan_splt[[2]]
    rec <- chan_splt[[3]]
    rec_ct <- chan_splt[[4]]
    if (lig %in% names(val_data_proc)) {
      # only add it if the ligand is expressed in the ligand ct
      if(!all(expr_use[[lig_ct]][,lig]==0)) {
        all_channels_in_val <- c(all_channels_in_val,full_channel)
      }
    }
  }
  
  return(all_channels_in_val)
}


## trying with full pb data
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=1,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')
dat <- list()
for (ct in names(pbmc_container$scMinimal_ctype)) {
  pb <- pbmc_container$scMinimal_ctype[[ct]]$pseudobulk
  dat[[ct]] <- pb
}
scITD_scale_pb_nnet <- dat
scITD_scale_pb_c2c <- dat
scITD_scale_pb_scT <- dat

# enrichment tests for all lr channel sets
top_n_targets <- 500
scITD_nnet_enr <- enrichment_test(scITD_nnet_full,top_n_targets=top_n_targets,scITD_scale_pb_nnet)
scITD_c2c_enr <- enrichment_test(scITD_c2c_full,top_n_targets=top_n_targets,scITD_scale_pb_c2c)
scITD_scT_enr <- enrichment_test(scITD_scT_full,top_n_targets=top_n_targets,scITD_scale_pb_scT)
nnet_enr <- enrichment_test(nnet_sig_channels_full,top_n_targets=top_n_targets,scITD_scale_pb_nnet)
c2c_enr <- enrichment_test(c2c_sig_channels_full,top_n_targets=top_n_targets,scITD_scale_pb_c2c)
scT_enr <- enrichment_test(scT_sig_channels_full,top_n_targets=top_n_targets,scITD_scale_pb_scT)

post_process_enr <- function(enr_res,thresh) {
  enr_res <- unlist(enr_res)

  num_tested <- length(enr_res)
  num_sig <- sum(enr_res<thresh)
  frac_sig <- num_sig/num_tested
  return(c(num_tested,num_sig,frac_sig))
}


thresh <- .10
methods_res <- list(scITD_nnet_enr,scITD_c2c_enr,
                    scITD_scT_enr,nnet_enr,
                    c2c_enr,scT_enr)
post_res_all <- lapply(methods_res,function(x) {
  post_process_enr(x,thresh)
})
post_res_all <- do.call(rbind.data.frame,post_res_all)
colnames(post_res_all) <- c('num_tested','num_sig','frac_sig')

# select equivalent numbers of randomly sampled lr channels
scITD_nnet_all_channels <- get_rand_channels(scITD_nnet_dat,scITD_scale_pb_nnet)
scITD_c2c_all_channels <- get_rand_channels(scITD_c2c_dat,scITD_scale_pb_c2c)
scITD_scT_all_channels <- get_rand_channels(scITD_scT_dat,scITD_scale_pb_scT)

rand_iter <- 10
rand_counts_all <- list()
rand_fracs_all <- list()
for (i in 1:rand_iter) {
  scITD_nnet_rand_chan <- sample(scITD_nnet_all_channels,length(scITD_nnet_enr))
  scITD_c2c_rand_chan <- sample(scITD_c2c_all_channels,length(scITD_c2c_enr))
  scITD_scT_rand_chan <- sample(scITD_scT_all_channels,length(scITD_scT_enr))
  nnet_rand_chan <- sample(scITD_nnet_all_channels,length(nnet_enr))
  c2c_rand_chan <- sample(scITD_c2c_all_channels,length(c2c_enr))
  scT_rand_chan <- sample(scITD_scT_all_channels,length(scT_enr))
  
  # enrichment tests for randomly sampled lr channels
  rand_scITD_nnet_enr <- enrichment_test(scITD_nnet_rand_chan,top_n_targets=top_n_targets,scITD_scale_pb_nnet)
  rand_scITD_c2c_enr <- enrichment_test(scITD_c2c_rand_chan,top_n_targets=top_n_targets,scITD_scale_pb_c2c)
  rand_scITD_scT_enr <- enrichment_test(scITD_scT_rand_chan,top_n_targets=top_n_targets,scITD_scale_pb_scT)
  rand_nnet_enr <- enrichment_test(nnet_rand_chan,top_n_targets=top_n_targets,scITD_scale_pb_nnet)
  rand_c2c_enr <- enrichment_test(c2c_rand_chan,top_n_targets=top_n_targets,scITD_scale_pb_c2c)
  rand_scT_enr <- enrichment_test(scT_rand_chan,top_n_targets=top_n_targets,scITD_scale_pb_scT)
  
  methods_res <- list(rand_scITD_nnet_enr,rand_scITD_c2c_enr,
                      rand_scITD_scT_enr,rand_nnet_enr,
                      rand_c2c_enr,rand_scT_enr)
  rand_post_res_all <- lapply(methods_res,function(x) {
    post_process_enr(x,thresh)
  })
  rand_post_res_all <- do.call(rbind.data.frame,rand_post_res_all)
  colnames(rand_post_res_all) <- c('num_tested','num_sig','frac_sig')
  rand_counts_all[[i]] <-  rand_post_res_all$num_sig
  rand_fracs_all[[i]] <- rand_post_res_all$frac_sig
}
rand_counts_all2 <- do.call(cbind.data.frame,rand_counts_all)
rand_fracs_all2 <- do.call(cbind.data.frame,rand_fracs_all)

av_counts_rand <- rowMeans(rand_counts_all2)
av_fracs_rand <- rowMeans(rand_fracs_all2)
sd_counts_rand <- transform(rand_counts_all2, SD=apply(rand_counts_all2,1, sd, na.rm = TRUE))$SD
se_counts_rand <- sd_counts_rand / sqrt(rand_iter)

post_res_all$posneg <- 'pos'
post_res_all_flip <- cbind.data.frame(post_res_all$num_tested,post_res_all$num_tested-post_res_all$num_sig,1-post_res_all$frac_sig,'neg')
colnames(post_res_all_flip) <- colnames(post_res_all)
post_res_all <- rbind.data.frame(post_res_all,post_res_all_flip)

rand_res_mat <- cbind.data.frame(post_res_all$num_tested[1:6],av_counts_rand,av_fracs_rand,se_counts_rand,'pos')
rand_res_mat_flip <- cbind.data.frame(post_res_all$num_tested[1:6],post_res_all$num_tested[1:6]-av_counts_rand,1-av_fracs_rand,NA,'neg')
colnames(rand_res_mat_flip) <- colnames(rand_res_mat) <- c('num_tested','num_sig','frac_sig','se_counts','posneg')
rand_res_mat <- rbind.data.frame(rand_res_mat,rand_res_mat_flip)

post_res_all$se_counts <- NA
rand_res_mat <- rand_res_mat[,colnames(post_res_all)]
full_res_mat <- rbind.data.frame(post_res_all,rand_res_mat)

full_res_mat$method <- rep(c('scITD_nnet','scITD_c2c','scITD_scT','NicheNet','TensorC2C','scTensor'),4)
full_res_mat$type <- c('real','real','real','real','real','real',
                       'real','real','real','real','real','real',
                       'random','random','random','random','random','random',
                       'random','random','random','random','random','random')

full_res_mat$frac <- round(full_res_mat$frac,digits=2)

full_res_mat$method <- factor(full_res_mat$method,levels = c('NicheNet','scITD_nnet','TensorC2C','scITD_c2c','scTensor','scITD_scT'))
full_res_mat$type <- factor(full_res_mat$type,levels = c('random','real'))

# saveRDS(full_res_mat,file='/home/jmitchel/data/lupus_data/full_lr_benchmarking_table6.rds')
full_res_mat <- readRDS(file='/home/jmitchel/data/lupus_data/full_lr_benchmarking_table6.rds')

p <- ggplot(full_res_mat,aes(y = type, x = num_sig, fill = posneg)) +
  geom_bar(position = "stack", width = .7,
           stat = "identity") +
  facet_grid(rows = vars(method), switch = "y") +
  geom_errorbar(aes(xmin = num_sig-se_counts, xmax = num_sig+se_counts),
                color="black",
                width=0.1, position='identity') +
  geom_text(aes(label = frac), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(breaks = c("pos", "neg"), values=c("light blue", "dark gray")) +
  theme_classic() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm"))
p

## Figure S7b
pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/lr_benchmarking_new2.pdf", useDingbats = FALSE,
    width = 5.5, height = 6.5)
p
dev.off()



## plotting jaccard overlap
jacc_nnet <- length(intersect(scITD_nnet_full,nnet_sig_channels_full)) / length(union(scITD_nnet_full,nnet_sig_channels_full))
jacc_c2c <- length(intersect(scITD_c2c_full,c2c_sig_channels_full)) / length(union(scITD_c2c_full,c2c_sig_channels_full))
jacc_scT <- length(intersect(scITD_scT_full,scT_sig_channels_full)) / length(union(scITD_scT_full,scT_sig_channels_full))

# get jaccard for randomly selected channels
rand_iter <- 10000
jacc_nnet_rand_all <- c()
jacc_c2c_rand_all <- c()
jacc_scT_rand_all <- c()
for (i in 1:rand_iter) {
  scITD_nnet_rand_chan <- sample(scITD_nnet_all_channels,length(scITD_nnet_enr))
  scITD_c2c_rand_chan <- sample(scITD_c2c_all_channels,length(scITD_c2c_enr))
  scITD_scT_rand_chan <- sample(scITD_scT_all_channels,length(scITD_scT_enr))
  
  jacc_nnet_rand <- length(intersect(scITD_nnet_rand_chan,nnet_sig_channels_full)) / length(union(scITD_nnet_rand_chan,nnet_sig_channels_full))
  jacc_c2c_rand <- length(intersect(scITD_c2c_rand_chan,c2c_sig_channels_full)) / length(union(scITD_c2c_rand_chan,c2c_sig_channels_full))
  jacc_scT_rand <- length(intersect(scITD_scT_rand_chan,scT_sig_channels_full)) / length(union(scITD_scT_rand_chan,scT_sig_channels_full))
  
  jacc_nnet_rand_all <- c(jacc_nnet_rand_all,jacc_nnet_rand)
  jacc_c2c_rand_all <- c(jacc_c2c_rand_all,jacc_c2c_rand)
  jacc_scT_rand_all <- c(jacc_scT_rand_all,jacc_scT_rand)
}

nnet_rand_mean <- mean(jacc_nnet_rand_all)
nnet_rand_sd <- sd(jacc_nnet_rand_all)

c2c_rand_mean <- mean(jacc_c2c_rand_all)
c2c_rand_sd <- sd(jacc_c2c_rand_all)

scT_rand_mean <- mean(jacc_scT_rand_all)
scT_rand_sd <- sd(jacc_scT_rand_all)

# calculating empirical pvalues
sum(jacc_nnet_rand_all > jacc_nnet) / length(jacc_nnet_rand_all)
sum(jacc_c2c_rand_all > jacc_c2c) / length(jacc_c2c_rand_all)
sum(jacc_scT_rand_all > jacc_scT) / length(jacc_scT_rand_all)

## columns are jaccard, sd, method
tmp <- cbind.data.frame(c(jacc_nnet, nnet_rand_mean, jacc_c2c, c2c_rand_mean, jacc_scT, scT_rand_mean),
                        c(NA, nnet_rand_sd, NA, c2c_rand_sd, NA, scT_rand_sd),
                        c('NicheNet','NicheNet','C2C','C2C','scTensor','scTensor'),
                        c('real','random','real','random','real','random'))
colnames(tmp) <- c('jaccard','std','method','type')
tmp$type <- factor(tmp$type,levels=c('real','random'))
tmp$method <- factor(tmp$method,levels=c('NicheNet','C2C','scTensor'))
p <- ggplot(tmp,aes(x=method,y=jaccard,fill=type)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  ggtitle('Overlap in inferred LR channels between methods') +
  geom_errorbar(aes(ymin=jaccard-std, ymax=jaccard+std), width=.2,
                position=position_dodge(.9)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

## Figure S7c
pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/scitd_nnet_c2c_jaccard3.pdf", useDingbats = FALSE,
    width = 4.5, height = 2.5)
p
dev.off()

# saveRDS(tmp,file='/home/jmitchel/data/lupus_data/lr_jaccard_table2.rds')


# things changed from previously
# type of enrichment test
# number of results included from the tests
# identification of target genes


