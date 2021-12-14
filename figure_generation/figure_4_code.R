library(scITD)
library(Seurat)
library(spqn) 
library(RColorBrewer)
library(iTALK)

# load up the subsetted dataset
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

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Ethnicity'),
                                        stat_use='pval')

## plot donor score
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# just to check that everything is as expected up to this point
pbmc_container$plots$donor_matrix


##### running LR analysis
# using cellchat database
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
  rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
  return(rname)
})
lr_pairs$interaction_name <- NULL

# infer active LR interactions
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='pool')
sft_thresh <- c(12,14,12,10,12,9,12)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

lr_hmap <- compute_LR_interact(pbmc_container, lr_pairs, sig_thresh=.00000000005,
                               percentile_exp_rec=0.85, add_ld_fact_sig=TRUE)

# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_new_lr8.pdf", useDingbats = FALSE,
#     width = 6, height = 7)
lr_hmap
# dev.off()


lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,mod_ct='Th',mod=5,lig_ct='cMono',lig='ICOSLG')
# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_ICOSLG_trio2.pdf", useDingbats = FALSE,
#     width = 6, height = 5)
lig_mod_fact
dev.off()

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=1,mod_ct='B',mod=1,lig_ct='cMono',lig='TNFSF13B')
# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_TNFSF13B_trio3.pdf", useDingbats = FALSE,
#     width = 6, height = 5)
lig_mod_fact
dev.off()


## extracting the pvalues for these two...
# running code within the compute_LR_interact fn
sig_thresh=.0000001
myres_mat <- pbmc_container$lr_res # at 285
container=pbmc_container

pbmc_container$lr_res['ICOSLG_cMono_ICOS','Th_m5']
which(rownames(myres_mat)=='ICOSLG_cMono_ICOS')
10**(-fact_res2[47,2])

pbmc_container$lr_res['TNFSF13B_cMono_TNFRSF13B','B_m1']
which(rownames(myres_mat)=='TNFSF13B_cMono_TNFRSF13B')
10**(-fact_res2[39,1])


# getting GO enrichment HMAPs for modules
ctypes <- c('Th')
modules <- c(5)

# using more stringent p-val threshold for figure plot
# mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.002, db_use=c('GO'),max_plt_pval=.002,h_w=c(7,3))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.005, 
                                 db_use=c('GO','BioCarta'),max_plt_pval=.005,h_w=c(14,2))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_ICOSLG_gsets2.pdf", useDingbats = FALSE,
#     width = 5, height = 8)
mod_enr
dev.off()

# also checking the enrichment of the TNFSF13B B module
ctypes <- c('B')
modules <- c(1)
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, 
                                 db_use=c('GO','Reactome','KEGG'),
                                 max_plt_pval=.05,h_w=c(20,5))
mod_enr



#### computing enrichment of modules with the nichenet scores
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

## ICOSLG target module first
test = ligand_target_matrix[,'ICOSLG']
t_mod = pbmc_container[["module_genes"]][["Th"]]
mymod <- c(5)
g_in_mod <- names(t_mod)[t_mod%in%mymod]
g_not_mod <- names(t_mod)[!(t_mod%in%mymod)]
tmp <- cbind.data.frame(c(g_in_mod,g_not_mod),
                        c(rep('Th_m5',length(g_in_mod)),rep('other',length(g_not_mod))))
colnames(tmp) <- c('gn','in_mod')
tmp$in_mod <- factor(tmp$in_mod,levels=c('Th_m5','other'))
target_scores <- test[tmp$gn]
tmp$target_scores <- target_scores
tmp <- tmp[which(!is.na(target_scores)),]
p <- ggplot(tmp,aes(x=as.factor(in_mod),y=target_scores)) +
  geom_boxplot(notch=TRUE) +
  ylab('ICOSLG NicheNet regulatory potential') +
  xlab('') +
  theme_bw()
test_res <- wilcox.test(target_scores~in_mod,data=tmp)
print(test_res$p.value)

pdf(file = "/home/jmitchel/figures/for_paper_v2/ICOSLG_NicheNet_enr2.pdf", useDingbats = FALSE,
    width = 4, height = 3.5)
p
dev.off()




## TNFSF13B target module now
test = ligand_target_matrix[,'TNFSF13B']
b_mod = pbmc_container[["module_genes"]][["B"]]
mymod <- c(1)
g_in_mod <- names(b_mod)[b_mod%in%mymod]
g_not_mod <- names(b_mod)[!(b_mod%in%mymod)]
tmp <- cbind.data.frame(c(g_in_mod,g_not_mod),
                        c(rep('B_m1',length(g_in_mod)),rep('other',length(g_not_mod))))
colnames(tmp) <- c('gn','in_mod')
tmp$in_mod <- factor(tmp$in_mod,levels=c('B_m1','other'))
target_scores <- test[tmp$gn]
tmp$target_scores <- target_scores
tmp <- tmp[which(!is.na(target_scores)),]
p <- ggplot(tmp,aes(x=as.factor(in_mod),y=target_scores)) +
  geom_boxplot(notch=TRUE) +
  ylab('TNFSF13B NicheNet regulatory potential') +
  xlab('') +
  theme_bw()
test_res <- wilcox.test(target_scores~in_mod,data=tmp)
print(test_res$p.value)

p















##### benchmarking analysis

## fns to use in benchmarking
test_lr_validity <- function(container,lr_test,ligand_target_matrix) {
  # all_pvals <- c()
  print(length(lr_test))
  all_pvals <- sccore::plapply(1:length(lr_test), function(lst_ndx) {
    lig_source_target <- lr_test[lst_ndx]
    lst_split <- strsplit(lig_source_target,split='_')[[1]]
    lig <- lst_split[[1]]
    source_ct <- lst_split[[2]]
    target_ct <- lst_split[[3]]
    
    lig_ct_exp <- container$scale_pb_extra[[source_ct]][,lig]
    
    all_cors <- cor_helper(container,lig_ct_exp,lig,source_ct,target_ct)
    
    # testing using wilcoxin test
    if (lig %in% colnames(ligand_target_matrix)) {
      gene_intersect <- intersect(names(all_cors),rownames(ligand_target_matrix))
      all_cors <- all_cors[gene_intersect]
      
      top_cors <- all_cors[order(abs(all_cors),decreasing=TRUE)][1:200]
      tmp <- cbind.data.frame(all_cors,ligand_target_matrix[gene_intersect,lig])
      colnames(tmp) <- c('mycors','target_scores')
      tmp$top <- sapply(1:nrow(tmp),function(i){
        if (rownames(tmp)[i] %in% names(top_cors)) {
          return('top_gn')
        } else {
          return('other')
        }
      })
      tmp$top <- as.factor(tmp$top)
      # my_res <- wilcox.test(target_scores~top,data=tmp)
      # pval <- my_res$p.value
      pval <- sig_test_helper(tmp,10000)
      return(pval)
    } else {
      return(NA)
    }
    
  },mc.preschedule=FALSE,n.cores=30,progress=TRUE)
  
  all_pvals2 <- all_pvals
  
  # add back names 
  names(all_pvals2) <- lr_test
  
  # remove na elements
  all_pvals2 <- all_pvals2[!is.na(all_pvals2)]
  
  # adjust pvalues
  all_pvals2 <- p.adjust(all_pvals2,method='fdr')
  
  # give warning if min non-zero adjusted pval not below .05
  if(min(all_pvals2[all_pvals2!=0])>.05) {
    print('Warning: minimum non-zero adjusted pvalue is > 0.05')
  }
  
  return(all_pvals2)
}


cor_helper <- function(container,lig_ct_exp,lig,source_ct,target_ct) {
  test <- cbind.data.frame(lig_ct_exp,container$scale_pb_extra[[target_ct]][names(lig_ct_exp),])
  cor_m <- cor(test)
  
  lig_ct_exp2 <- container$no_scale_pb_extra[[source_ct]][,lig]
  test2 <- cbind.data.frame(lig_ct_exp2,container$no_scale_pb_extra[[target_ct]][names(lig_ct_exp2),])
  ave_exp <- colMeans(test2)
  
  cor_m_spqn <- normalize_correlation(cor_m, ave_exp=ave_exp, ngrp=20, size_grp=300, ref_grp=5)
  all_cors <- cor_m_spqn[1,]
  names(all_cors) <- colnames(cor_m)
  all_cors <- all_cors[2:length(all_cors)]
  return(all_cors)
}

sig_test_helper <- function(tmp,n_iter) {
  
  ## to use with wilcoxon version
  my_res <- wilcox.test(target_scores~top,data=tmp)
  obs_pval <- my_res$p.value
  
  all_pv <- c()
  for (i in 1:n_iter) {
    tmp2 <- transform(tmp, top = sample(top) )
    my_res <- wilcox.test(target_scores~top,data=tmp2)
    pval <- my_res$p.value
    all_pv <- c(all_pv,pval)
  }
  pval <- sum(all_pv<=obs_pval) / length(all_pv)
  
  return(pval)
}



## for parsing my output to use in benchmarking... need to get LR channels out
parse_scITD_LR_out <- function(container) {
  myres_mat <- container$lr_res #includes NS results still but are FDR adjusted
  
  # reduce to rows/columns with at least one significant hit
  myres_mat <- myres_mat[rowSums(myres_mat<.0001)>0,]
  myres_mat <- myres_mat[,colSums(myres_mat<.0001)>0]
  
  ## get unique ligand_source_target channels from my analysis
  unique_channels <- c()
  for (i in 1:nrow(myres_mat)) {
    lig_ct_rec <- myres_mat[i,]
    lig_ct_rec_name <- strsplit(rownames(myres_mat)[i],split='_')[[1]]
    lig <- lig_ct_rec_name[[1]]
    source <- lig_ct_rec_name[[2]]
    
    for (j in 1:ncol(myres_mat)) {
      pv <- lig_ct_rec[j]
      if (pv < 0.0001) {
        target_ct <- strsplit(names(pv),split="_")[[1]][[1]]
        lig_source_target <- paste0(lig,"_",source,"_",target_ct)
        unique_channels <- c(unique_channels,lig_source_target)
      }
    }
  }
  scITD_channels <- unique(unique_channels)
  return(scITD_channels)
}


## for getting channels by the LR association method
# need to have run LR preprocessing with appropriate lr_pairs db first!
get_LR_association_chan <- function(container,lr_pairs) {
  my_channels <- c()
  all_pvals <- c()
  ctypes <- container$experiment_params$ctypes_use
  print(nrow(lr_pairs))
  for (i in 1:nrow(lr_pairs)) {
    if (i %% 100 == 0) {
      print(i)
    }
    ligand <- lr_pairs[i,1]
    rec_elements <- strsplit(lr_pairs[i,2],split='_')[[1]]
    for (source_ct in ctypes) {
      for (target_ct in ctypes) {
        if (source_ct != target_ct) {
          # loop through receptor components
          for (receptor in rec_elements) {
            if (ligand %in% colnames(container$scale_pb_extra[[source_ct]]) && receptor %in% colnames(container$scale_pb_extra[[target_ct]])) {
              lig_exp <- container$scale_pb_extra[[source_ct]][,ligand]
              rec_exp <- container$scale_pb_extra[[target_ct]][,receptor]
              
              # check that expression is not 0
              if (sum(lig_exp>.2)>(.01*length(lig_exp)) && sum(rec_exp>.2)>(.01*length(rec_exp))) {
                tmp <- cbind.data.frame(lig_exp,rec_exp[names(lig_exp)])
                colnames(tmp) <- c('l_exp','r_exp')
                lmres <- summary(lm(l_exp~r_exp,data=tmp))
                pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
                chan <- paste0(ligand,'_',source_ct,'_',target_ct)
                my_channels <- c(my_channels,chan)
                all_pvals <- c(all_pvals,pval)
              }
            }
          }
        }
      }
    }
  }
  
  # adjust pvals
  all_pvals <- p.adjust(all_pvals,method='fdr')
  
  ## get channels where there was significant associations with any receptor components
  my_channels <- my_channels[all_pvals<.05]
  my_channels <- unique(my_channels)
  return(my_channels)
}



## get n random LR channels
get_rand_chan <- function(container,lr_pairs,n_samp) {
  lr_pairs_sub <- lr_pairs[lr_pairs$ligand %in% colnames(container$scale_pb_extra[[1]]),]
  
  # get all two cell type combinations (no repeats)
  all_ctypes <- container$experiment_params$ctypes_use
  ct_combos <- combn(all_ctypes,2)
  ct_combos2 <- rbind(ct_combos[2,],ct_combos[1,])
  ct_combos <- cbind.data.frame(ct_combos,ct_combos2)
  
  unq_lig <- unique(lr_pairs_sub$ligand)
  
  all_chan <- sapply(1:length(unq_lig), function(i) {
    lig <- unq_lig[i]
    lig_source_target_combos <- sapply(1:ncol(ct_combos), function(j) {
      paste0(lig,"_",ct_combos[1,j],"_",ct_combos[2,j])
    })
  })
  my_samp <- sample(all_chan,n_samp)
  return(my_samp)
}




# parse scITD results for cell chat db which was used in above analysis
scITD_channels <- parse_scITD_LR_out(pbmc_container)

# compute NicheNet enrichment for each channel
scITD_res <- test_lr_validity(pbmc_container,scITD_channels,ligand_target_matrix)
scITD_frac <- sum(scITD_res<.05)/length(scITD_res)

## object to store results for plotting later
res_dat <- matrix(ncol=5,nrow=50)
colnames(res_dat) <- c('db','method','n_chan','n_chan_eval','percent_enr')
res_dat[1,] <-  c('CellChat','scITD',length(scITD_channels),length(scITD_res),scITD_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')

# making barplot for full scITD results on CellChat lr pairs
tmp <- cbind.data.frame(names(scITD_res),scITD_res)
colnames(tmp) <- c('channel','adj_pval')
tmp <- tmp[order(tmp$adj_pval,decreasing=FALSE),]
tmp$adj_pval[tmp$adj_pval==0] <- .0001
tmp$channel <- factor(tmp$channel,levels=tmp$channel)
tmp$adj_pval <- -log10(tmp$adj_pval)
mylabels <- tmp$channel
print(which(mylabels=='ICOSLG_cMono_Th'))
print(which(mylabels=='TNFSF13B_cMono_B'))
ndx_keep <- c(1,4,8,13,17,22,26)
mylabels2 <- rep('',length(mylabels))
mylabels2[ndx_keep] <- as.character(mylabels)[ndx_keep]
chan_bplot <- ggplot(tmp,aes(x=channel,y=adj_pval)) +
  geom_bar(stat='identity') +
  scale_x_discrete(labels= mylabels2) +
  scale_y_continuous(limits = c(0,4.25), expand = c(0, 0)) +
  ylab('Top target genes-NicheNet\nenrichment -log10(padj)') +
  xlab('scITD inferred channels (ligand_source_target)') +
  geom_hline(yintercept = -log10(.05), 
             color = "red", size=.75) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -30, vjust = 0.5, hjust=0),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/LR_chan_bplot3_2.pdf", useDingbats = FALSE,
#     width = 6, height = 4.5)
chan_bplot
dev.off()







## get cell chat inferred channels by LR association method
cchat_assoc_channels <- get_LR_association_chan(pbmc_container,lr_pairs)
cchat_assoc_res <- test_lr_validity(pbmc_container,cchat_assoc_channels,ligand_target_matrix)
cchat_assoc_frac <- sum(cchat_assoc_res<.05)/length(cchat_assoc_res)
res_dat[2,] <-  c('CellChat','LR_association',length(cchat_assoc_channels),length(cchat_assoc_res),cchat_assoc_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')

rand_channels <- get_rand_chan(pbmc_container,lr_pairs,n_samp=250)
rand_res <- test_lr_validity(pbmc_container,rand_channels,ligand_target_matrix)
rand_frac <- sum(rand_res<.05)/length(rand_res)
res_dat[3,] <-  c('CellChat','random',length(rand_channels),length(rand_res),rand_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')




### now for iTALK db analysis
# load LR pairs 
dim(iTALK::database)
head(database)
lr_pairs <- database[,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')]
colnames(lr_pairs) <- c('ligand','receptor')
lr_pairs <- unique(lr_pairs)

# run scITD analysis using these LR pairs
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='pool')
sft_thresh <- c(12,14,12,10,12,9,12)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)
lr_hmap <- compute_LR_interact(pbmc_container, lr_pairs, sig_thresh=.00000000005,
                               percentile_exp_rec=0.85, add_ld_fact_sig=TRUE)

# parse scITD inferred channels
scITD_channels <- parse_scITD_LR_out(pbmc_container)
scITD_res <- test_lr_validity(pbmc_container,scITD_channels,ligand_target_matrix)
scITD_frac <- sum(scITD_res<.05)/length(scITD_res)
res_dat[4,] <-  c('iTALK','scITD',length(scITD_channels),length(scITD_res),scITD_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')


italk_assoc_channels <- get_LR_association_chan(pbmc_container,lr_pairs)
italk_assoc_res <- test_lr_validity(pbmc_container,italk_assoc_channels,ligand_target_matrix)
italk_assoc_frac <- sum(italk_assoc_res<.05)/length(italk_assoc_res)
res_dat[5,] <-  c('iTALK','LR_association',length(italk_assoc_channels),
                  length(italk_assoc_res),italk_assoc_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')


rand_channels <- get_rand_chan(pbmc_container,lr_pairs,n_samp=250)
rand_res <- test_lr_validity(pbmc_container,rand_channels,ligand_target_matrix)
rand_frac <- sum(rand_res<.05)/length(rand_res)
res_dat[6,] <-  c('iTALK','random',length(rand_channels),length(rand_res),rand_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')


# res_dat <- readRDS(file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')

## infer LR channels using singlecellsig LR pair db
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/singlecellsignalr_LR.csv')

# infer active LR interactions
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='pool')
sft_thresh <- c(12,14,12,10,12,9,12)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)
lr_hmap <- compute_LR_interact(pbmc_container, lr_pairs, sig_thresh=.00000000005,
                               percentile_exp_rec=0.85, add_ld_fact_sig=TRUE)

scITD_channels <- parse_scITD_LR_out(pbmc_container)
scITD_res <- test_lr_validity(pbmc_container,scITD_channels,ligand_target_matrix)
scITD_frac <- sum(scITD_res<.05)/length(scITD_res)
res_dat[7,] <-  c('SingleCellSignalR','scITD',length(scITD_channels),length(scITD_res),scITD_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')


scs_assoc_channels <- get_LR_association_chan(pbmc_container,lr_pairs)
scs_assoc_res <- test_lr_validity(pbmc_container,scs_assoc_channels,ligand_target_matrix)
scs_assoc_frac <- sum(scs_assoc_res<.05)/length(scs_assoc_res)
res_dat[8,] <-  c('SingleCellSignalR','LR_association',length(scs_assoc_channels),
                  length(scs_assoc_res),scs_assoc_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')


rand_channels <- get_rand_chan(pbmc_container,lr_pairs,n_samp=250)
rand_res <- test_lr_validity(pbmc_container,rand_channels,ligand_target_matrix)
rand_frac <- sum(rand_res<.05)/length(rand_res)
res_dat[9,] <-  c('SingleCellSignalR','random',length(rand_channels),length(rand_res),rand_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')

res_dat <- readRDS(file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')

res_dat <- as.data.frame(res_dat)
res_dat$n_chan <- as.numeric(res_dat$n_chan)
res_dat$n_chan_eval <- as.numeric(res_dat$n_chan_eval)
res_dat$percent_enr <- as.numeric(res_dat$percent_enr)
res_dat$db <- as.factor(res_dat$db)
res_dat$method <- factor(res_dat$method,levels=c('random','LR_association','scITD'))

## adding in se error
res_dat$se_dat <- sqrt(res_dat$percent_enr*(1-res_dat$percent_enr)/res_dat$n_chan_eval)

mycol = brewer.pal(n = 8, name = "Dark2")
comp_bplot_size <- ggplot(res_dat,aes(x=db,y=n_chan,fill=method)) +
  geom_bar(stat='identity', position="dodge") +
  ylab('Number of inferred interactions') +
  xlab('LR database') +
  guides(fill=guide_legend(title="LR method")) +
  scale_fill_manual(values=c(mycol[8], mycol[1], mycol[6])) +
  theme_bw() +
  theme(legend.position = "none")

comp_bplot_frac <- ggplot(res_dat,aes(x=db,y=percent_enr,fill=method)) +
  geom_bar(stat='identity', position="dodge") +
  geom_errorbar(aes(ymin = percent_enr - se_dat, ymax = percent_enr + se_dat), 
                width=0.15, position=position_dodge(.9)) +
  ylab('Fraction of inferred interactions\nwith NicheNet enrichment') +
  xlab('LR database') +
  guides(fill=guide_legend(title="LR method")) +
  scale_fill_manual(values=c(mycol[8], mycol[1], mycol[6])) +
  theme_bw()

combined_plt <- cowplot::plot_grid(comp_bplot_size,comp_bplot_frac,nrow=1, rel_widths = c(.7,1.1))

# pdf(file = "/home/jmitchel/figures/for_paper_v2/LR_all_db_bplot3.pdf", useDingbats = FALSE,
#     width = 7.5, height = 3.5)
combined_plt
# dev.off()














##### testing for association between Treg proportions and F3
pbmc_container$subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')
t4_sub <- colnames(pbmc_container$scMinimal_ctype[['Th']]$count_data)

tmp <- pbmc@meta.data[t4_sub,'ct_cov']

ndx_mark <- which(tmp=='T4reg')
ndx_other <- which(tmp!='T4reg')

tmp <- as.character(tmp)
tmp[ndx_mark] <- 1
tmp[ndx_other] <- 2

pbmc_container[["subclusters"]][["T4"]][["res:0.6"]] <- as.numeric(tmp)
names(pbmc_container[["subclusters"]][["T4"]][["res:0.6"]]) <- t4_sub

ctype <- 'T4'
res <- 0.6
subtype=1
factor_use=3
ctype_cur='Th'

resolution_name <- paste0('res:',as.character(res))
subclusts <- pbmc_container$subclusters[[ctype]][[resolution_name]]

# append large cell type name to subclusters
subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})

# limit cells in subclusts to those that we actually have scores for
donor_scores <- pbmc_container$tucker_results[[1]]
cell_intersect <- intersect(names(subclusts),rownames(pbmc_container$scMinimal_full$metadata))
donor_vec <- pbmc_container$scMinimal_full$metadata[cell_intersect,'donors']
subclusts <- subclusts[cell_intersect]
subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]

# make subtype association plot
subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
scMinimal <- pbmc_container$scMinimal_ctype[[ctype_cur]]
sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

# get donor proportions of subclusters
donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)
donor_props <- donor_props[,subtype,drop=FALSE]
colnames(donor_props) <- 'prop'

# append dscores for factor 1
donor_props2 <- cbind(donor_props,donor_scores[rownames(donor_props),factor_use])
colnames(donor_props2)[ncol(donor_props2)] <- 'dsc'

donor_props2 <- as.data.frame(donor_props2)
donor_props2$dsc <- as.numeric(donor_props2$dsc)
donor_props2$prop <- as.numeric(donor_props2$prop)

lmres <- lm(prop~dsc,data=donor_props2)
line_range <- seq(min(donor_props2$dsc),max(donor_props2$dsc),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

p <- ggplot(donor_props2,aes(x=dsc,y=prop)) +
  geom_point(alpha = 0.75,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab(paste0('Factor ',as.character(factor_use),' Donor Score')) +
  ylab(paste0('Proportion Treg/Th')) +
  ylim(0,1) +
  ggtitle(paste0(ctype,'_',as.character(subtype),' Proportions')) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  ylim(0,.25) + ggtitle('')

p
dev.off()

donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)
donor_balances <- coda.base::coordinates(donor_props)
rownames(donor_balances) <- rownames(donor_props)

f1 <- get_one_factor(pbmc_container,3)
f1_dsc <- f1[[1]]
tmp <- cbind.data.frame(f1_dsc[rownames(donor_balances),1,drop=FALSE],donor_balances)
colnames(tmp) <- c('dsc','my_balance')
lmres <- summary(lm(my_balance~dsc,data=tmp))
pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
print(pval)



