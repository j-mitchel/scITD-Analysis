library(Seurat)

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

# set up project parameters
# param_list <- initialize_params(ctypes_use = c("B","NK","T4","T8","cDC",
#                                                "cM","ncM"),
#                                 ncores = 30, rand_seed = 10)
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
                                                       "Ethnicity"))


# pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
#                               norm_method='trim', scale_factor=10000,
#                               vargenes_method='norm_var_pvals', vargenes_thresh=.05,
#                               scale_var = TRUE, var_scale_power = 1.5)
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5)

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(25,37,7),
                                 tucker_type = 'regular', rotation_type = 'ica')


# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('pool','processing'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('pool','processing'),
                                    cluster_by_meta='pool',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='rsq')

pdf(file = "/home/jmitchel/figures/for_paper_v2/lupus_batch_dscores.pdf", useDingbats = FALSE,
    width = 7, height = 6.5)
pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 7, height = 6.5)
pbmc_container$plots$donor_matrix
dev.off()




# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)


## trying to draw loadings plots individually because I want them to have different params...
pdf(file = "/home/jmitchel/figures/for_paper_v2/batch_lds_f6.pdf", useDingbats = FALSE,
    width = 3.5, height = 4.5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=6, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=F)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/batch_lds_f2.pdf", useDingbats = FALSE,
    width = 3.5, height = 4.5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=F)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/batch_lds_f1.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=F)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/batch_lds_f10.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=10, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=F)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/batch_lds_f19.pdf", useDingbats = FALSE,
    width = 5, height = 5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=19, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=F)
dev.off()




# read in soup profile for a batch
soupProf <- readRDS(file="/home/jmitchel/data/lupus_data/SLE_droplets/YE_7-19/soupProf.rds")


test_soup_association <- function(container,my_factor,b_direc,comp_type='any_up',soupProf) {
  f_data <- get_one_factor(container, factor_select=my_factor)
  dscores <- f_data[[1]]
  lds <- f_data[[2]]
  
  sig_vectors <- get_significance_vectors(container,
                                          factor_select=my_factor, colnames(lds))
  # convert list to df
  sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
  
  # limit to just the genes in tmp_casted_num
  sig_df <- sig_df[rownames(lds),colnames(lds)]
  
  
  if (comp_type=="any_up") {
    g_keep <- c()
    for (i in 1:nrow(sig_df)) {
      for (j in 1:ncol(sig_df)) {
        sig_val <- sig_df[i,j]
        if (sig_val < 0.05) {
          if (b_direc=='up') {
            if (lds[i,j] > 0) {
              g_keep <- c(g_keep,rownames(sig_df)[i])
              break
            }
          } else if (b_direc=='down') {
            if (lds[i,j] < 0) {
              g_keep <- c(g_keep,rownames(sig_df)[i])
              break
            }
          }
        }
      }
    }
    
    g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
    g1 <- soupProf[g_keep,'est'] # set to length or gc
    g2 <- soupProf[g_not_keep,'est']
    g1 <- g1[!is.na(g1)]
    g1_nonzero <- g1[g1!=0]
    g1_nonzero_min <- min(g1_nonzero)
    g1[g1==0] <- g1_nonzero_min
    g2 <- g2[!is.na(g2)]
    g2_nonzero <- g2[g2!=0]
    g2_nonzero_min <- min(g2_nonzero)
    g2[g2==0] <- g2_nonzero_min
    g1 <- log10(g1)
    g2 <- log10(g2)
    tres <- t.test(g1, g2, alternative = "two.sided")
    
    print(tres$p.value)
    
    all_p <- c(g1,g2)
    all_t <- c(rep('Upreg genes',length(g1)),rep('NS genes',length(g2)))
    all_f <- rep(paste0('factor_',my_factor),length(all_p))
    mydf <- cbind.data.frame(all_p,all_t,all_f)
    mydf <- as.data.frame(mydf)
    
    return(list(tres$p.value,mydf))
  } else if (comp_type=='any') {
    g_keep <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
    g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
    g1 <- soupProf[g_keep,'est'] # set to length or gc
    g2 <- soupProf[g_not_keep,'est']
    g1 <- g1[!is.na(g1)]
    g1_nonzero <- g1[g1!=0]
    g1_nonzero_min <- min(g1_nonzero)
    g1[g1==0] <- g1_nonzero_min
    g2 <- g2[!is.na(g2)]
    g2_nonzero <- g2[g2!=0]
    g2_nonzero_min <- min(g2_nonzero)
    g2[g2==0] <- g2_nonzero_min
    g1 <- log10(g1)
    g2 <- log10(g2)
    tres <- t.test(g1, g2, alternative = "two.sided")
    
    print(tres$p.value)
    
    all_p <- c(g1,g2)
    all_t <- c(rep('sig_genes',length(g1)),rep('NS_genes',length(g2)))
    all_f <- rep(paste0('factor_',my_factor),length(all_p))
    mydf <- cbind.data.frame(all_p,all_t,all_f)
    mydf <- as.data.frame(mydf)
    
    return(list(tres$p.value,mydf))
  } else if (comp_type=='sig_up_consist') {
    all_sig <- rownames(sig_df)[rowSums(sig_df<0.05)==ncol(lds)]
    # keep only genes with consistent signs in the loadings in direction of upregulation
    tmp_lds <- lds[all_sig,]
    if (b_direc=='up') {
      g_keep <- rownames(tmp_lds)[rowSums(tmp_lds>0)==ncol(lds)]
    } else if (b_direc=='down') {
      g_keep <- rownames(tmp_lds)[rowSums(tmp_lds<0)==ncol(lds)]
    }
    
    g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
    
    g1 <- soupProf[g_keep,'est']
    g2 <- soupProf[g_not_keep,'est']
    tres <- t.test(g1, g2, alternative = "two.sided")
    
    print(tres$p.value)
    
    all_p <- c(g1,g2)
    all_t <- c(rep('g1',length(g1)),rep('g2',length(g2)))
    mydf <- cbind.data.frame(all_p,all_t)
    mydf <- as.data.frame(mydf)
    return(list(tres$p.value,mydf))
  } else if (comp_type=='all_v_some') {
    any_up <- c()
    for (i in 1:nrow(sig_df)) {
      for (j in 1:ncol(sig_df)) {
        sig_val <- sig_df[i,j]
        if (sig_val < 0.05) {
          if (b_direc=='up') {
            if (lds[i,j] > 0) {
              any_up <- c(any_up,rownames(sig_df)[i])
              break
            }
          } else if (b_direc=='down') {
            if (lds[i,j] < 0) {
              any_up <- c(any_up,rownames(sig_df)[i])
              break
            }
          }
        }
      }
    }
    
    
    all_sig <- rownames(sig_df)[rowSums(sig_df<0.05)==ncol(lds)]
    # keep only genes with consistent signs in the loadings in direction of upregulation
    tmp_lds <- lds[all_sig,]
    if (b_direc=='up') {
      g_keep <- rownames(tmp_lds)[rowSums(tmp_lds>0)==ncol(lds)]
    } else if (b_direc=='down') {
      g_keep <- rownames(tmp_lds)[rowSums(tmp_lds<0)==ncol(lds)]
    }
    
    g_not_keep <- any_up[!(any_up %in% g_keep)]
    g1 <- soupProf[g_keep,'est'] # set to length or gc
    g2 <- soupProf[g_not_keep,'est']
    g1 <- g1[!is.na(g1)]
    g1_nonzero <- g1[g1!=0]
    g1_nonzero_min <- min(g1_nonzero)
    g1[g1==0] <- g1_nonzero_min
    g2 <- g2[!is.na(g2)]
    g2_nonzero <- g2[g2!=0]
    g2_nonzero_min <- min(g2_nonzero)
    g2[g2==0] <- g2_nonzero_min
    g1 <- log10(g1)
    g2 <- log10(g2)
    tres <- t.test(g1, g2, alternative = "two.sided")
    
    print(tres$p.value)
    
    all_p <- c(g1,g2)
    all_t <- c(rep('Upreg all\ncell types',length(g1)),rep('Upreg some\ncell types',length(g2)))
    all_f <- rep(paste0('factor_',my_factor),length(all_p))
    mydf <- cbind.data.frame(all_p,all_t,all_f)
    mydf <- as.data.frame(mydf)
    
    return(list(tres$p.value,mydf))
  } else if (comp_type=='consist_up_v_down') {
    all_sig <- rownames(sig_df)[rowSums(sig_df<0.05)==ncol(lds)]
    # keep only genes with consistent signs in the loadings in direction of upregulation
    tmp_lds <- lds[all_sig,]
    if (b_direc=='up') {
      g_keep <- rownames(tmp_lds)[rowSums(tmp_lds>0)==ncol(lds)]
      g_not_keep <- rownames(tmp_lds)[rowSums(tmp_lds<0)==ncol(lds)]
    } else if (b_direc=='down') {
      g_keep <- rownames(tmp_lds)[rowSums(tmp_lds<0)==ncol(lds)]
      g_not_keep <- rownames(tmp_lds)[rowSums(tmp_lds>0)==ncol(lds)]
    }
    
    g1 <- soupProf[g_keep,'est']
    g2 <- soupProf[g_not_keep,'est']
    g1 <- g1[!is.na(g1)]
    g1_nonzero <- g1[g1!=0]
    g1_nonzero_min <- min(g1_nonzero)
    g1[g1==0] <- g1_nonzero_min
    g2 <- g2[!is.na(g2)]
    g2_nonzero <- g2[g2!=0]
    g2_nonzero_min <- min(g2_nonzero)
    g2[g2==0] <- g2_nonzero_min
    g1 <- log10(g1)
    g2 <- log10(g2)
    tres <- t.test(g1, g2, alternative = "two.sided")
    
    print(tres$p.value)
    
    all_p <- c(g1,g2)
    all_t <- c(rep('g1',length(g1)),rep('g2',length(g2)))
    mydf <- cbind.data.frame(all_p,all_t)
    mydf <- as.data.frame(mydf)
    return(list(tres$p.value,mydf))
  }  else if (comp_type=='any_up_v_down') {
    any_up <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
    lds_sub <- lds[any_up,]
    if (b_direc=='up') {
      any_up <- rownames(lds_sub)[rowSums(lds_sub>0)>0]
    } else if (b_direc=='down') {
      any_up <- rownames(lds_sub)[rowSums(lds_sub<0)>0]
    }
    
    any_down <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
    lds_sub <- lds[any_down,]
    if (b_direc=='up') {
      any_down <- rownames(lds_sub)[rowSums(lds_sub<0)>0]
    } else if (b_direc=='down') {
      any_down <- rownames(lds_sub)[rowSums(lds_sub>0)>0]
    }
    
    # only keep genes not in intersection of any_up and any_down
    int_genes <- intersect(any_up,any_down)
    any_up <- any_up[!(any_up %in% int_genes)]
    any_down <- any_down[!(any_down %in% int_genes)]
    
    g1 <- soupProf[any_up,'est']
    g2 <- soupProf[any_down,'est']
    tres <- t.test(log(g1), log(g2), alternative = "two.sided")
    
    print(tres$p.value)
    
    all_p <- c(g1,g2)
    all_t <- c(rep('g1',length(g1)),rep('g2',length(g2)))
    mydf <- cbind.data.frame(all_p,all_t)
    mydf <- as.data.frame(mydf)
    return(list(tres$p.value,mydf))
  }
}

# res1 <- test_soup_association(pbmc_container,6,'down','any',soupProf)
res1 <- test_soup_association(pbmc_container,6,'down','any_up',soupProf)
# res1 <- test_soup_association(pbmc_container,6,'down','consist_up_v_down',soupProf)

p <- ggplot(res1[[2]],aes(x = as.factor(all_t), y = all_p)) +
  geom_boxplot() +
  xlab('') +
  ylab('log10(soup fraction)') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))
pdf(file = "/home/jmitchel/figures/for_paper_v2/soup_any_up.pdf", useDingbats = FALSE,
    width = 3.5, height = 4)
pdf(file = "/home/jmitchel/figures/for_paper_v2/soup_any.pdf", useDingbats = FALSE,
    width = 3.5, height = 4)
p
dev.off()

res2 <- test_soup_association(pbmc_container,6,'down','all_v_some',soupProf)
p <- ggplot(res2[[2]],aes(x = as.factor(all_t), y = all_p)) +
  geom_boxplot() +
  xlab('') +
  ylab('log10(soup fraction)') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))
pdf(file = "/home/jmitchel/figures/for_paper_v2/soup_all_v_some.pdf", useDingbats = FALSE,
    width = 3.5, height = 4)
p
dev.off()



## now testing the batch factors for associations with GC content
test_gc_association <- function(container,my_factor,b_direc,comp_type='sig') {
  f_data <- get_one_factor(container, factor_select=my_factor)
  dscores <- f_data[[1]]
  lds <- f_data[[2]]
  
  if (is.null(container$gen_dat)) {
    library(biomaRt)
    library(EDASeq)
    
    ensembl = useMart("ensembl",
                      dataset="hsapiens_gene_ensembl")
    ens_id <- getBM(attributes=c('ensembl_gene_id',
                                 'hgnc_symbol'),
                    filters = 'hgnc_symbol',
                    mart = ensembl,
                    values = rownames(lds))
    
    gc_len <- getGeneLengthAndGCContent(ens_id[,1], 'hsa')
    
    gen_dat <- cbind.data.frame(ens_id,gc_len)
    
    gen_dat <- gen_dat[!(duplicated(gen_dat[,2]) | duplicated(gen_dat[,2], fromLast=TRUE)),]
    
    container$gen_dat <- gen_dat
  }
  
  gen_dat <- container$gen_dat
  rownames(gen_dat) <- gen_dat$hgnc_symbol
  
  sig_vectors <- get_significance_vectors(container,
                                          factor_select=my_factor, colnames(lds))
  # convert list to df
  sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
  
  # limit to just the genes in tmp_casted_num
  sig_df <- sig_df[rownames(lds),colnames(lds)]
  
  # loop through cell types to calculate 
  for (j in 1:ncol(lds)) {
    
    ct_ndx <- j
    ct_name <- colnames(lds)[j]
    
    if (comp_type=='sig') {
      lds_sig <- lds[sig_df[,ct_ndx]<0.05,ct_ndx,drop=FALSE]
      g_keep <- rownames(lds_sig)
      g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
      g1 <- gen_dat[g_keep,'gc'] # set to length or gc
      g2 <- gen_dat[g_not_keep,'gc']
      g1 <- g1[!is.na(g1)]
      g2 <- g2[!is.na(g2)]
      tres <- t.test(g1, g2, alternative = "two.sided")
      
      print(ct_name)
      print(tres$p.value)
    }
  }
  
  if (comp_type=="any_up") {
    # g_keep <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
    # # only keep ones that are upregulated in any cell type
    # lds_sub <- lds[g_keep,]
    # if (b_direc=='up') {
    #     g_keep <- rownames(lds_sub)[rowSums(lds_sub>0)>0]
    # } else if (b_direc=='down') {
    #     g_keep <- rownames(lds_sub)[rowSums(lds_sub<0)>0]
    # }
    
    g_keep <- c()
    for (i in 1:nrow(sig_df)) {
      for (j in 1:ncol(sig_df)) {
        sig_val <- sig_df[i,j]
        if (sig_val < 0.05) {
          if (b_direc=='up') {
            if (lds[i,j] > 0) {
              g_keep <- c(g_keep,rownames(sig_df)[i])
              break
            }
          } else if (b_direc=='down') {
            if (lds[i,j] < 0) {
              g_keep <- c(g_keep,rownames(sig_df)[i])
              break
            }
          }
        }
      }
    }
    
    g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
    g1 <- gen_dat[g_keep,'gc'] # set to length or gc
    g2 <- gen_dat[g_not_keep,'gc']
    g1 <- g1[!is.na(g1)]
    g2 <- g2[!is.na(g2)]
    tres <- t.test(g1, g2, alternative = "two.sided")
    
    print(tres$p.value)
    
    all_p <- c(g1,g2)
    all_t <- c(rep('sig_genes',length(g1)),rep('NS_genes',length(g2)))
    all_f <- rep(paste0('factor_',my_factor),length(all_p))
    mydf <- cbind.data.frame(all_p,all_t,all_f)
    mydf <- as.data.frame(mydf)
    
    return(list(tres$p.value,mydf))
  } else if (comp_type=='any') {
    g_keep <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
    g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
    g1 <- gen_dat[g_keep,'gc'] # set to length or gc
    g2 <- gen_dat[g_not_keep,'gc']
    g1 <- g1[!is.na(g1)]
    g2 <- g2[!is.na(g2)]
    tres <- t.test(g1, g2, alternative = "two.sided")
    
    print(tres$p.value)
    
    all_p <- c(g1,g2)
    all_t <- c(rep('sig_genes',length(g1)),rep('NS_genes',length(g2)))
    all_f <- rep(paste0('factor_',my_factor),length(all_p))
    mydf <- cbind.data.frame(all_p,all_t,all_f)
    mydf <- as.data.frame(mydf)
    
    return(list(tres$p.value,mydf))
  }
}

pv <- test_gc_association(pbmc_container,my_factor=14,b_direc='up',comp_type='any_up') #for testing the fn

all_pv <- c()
pv_list <- list()
my_direcs <- c('down','up','down','up','up','down','up','up','up','up',
               'down','up','up','up','up','up','up','up','up','up','up',
               'up','down','up','up')
# my_direcs <- rep('up',ncol(pbmc_container[["tucker_results"]][[1]]))
# my_direcs <- rep('down',ncol(pbmc_container[["tucker_results"]][[1]]))
for (i in 1:ncol(pbmc_container[["tucker_results"]][[1]])) {
  print(i)
  pv <- test_gc_association(pbmc_container,my_factor=i,b_direc=my_direcs[i],comp_type='any_up')
  # pv <- test_gc_association(pbmc_container,my_factor=i,b_direc='up',comp_type='any')
  all_pv <- c(all_pv,pv[[1]])
  pv_list[[i]] <- pv[[2]]
}
all_pv_adj <- p.adjust(all_pv,method='fdr')

# pv_no_combat <- all_pv
# pv_combat <- all_pv
# all_pv_adj <- p.adjust(c(pv_no_combat,pv_combat),method='fdr')
# 
# plt_no_combat <- pv_list
# plt_combat <- pv_list


all_pv_adj[c(1,2,3,4,5,6,7,10,11,19,23)]
full_df_batch <- rbind.data.frame(pv_list[[1]],pv_list[[2]],pv_list[[3]],
                                  pv_list[[4]],pv_list[[5]],pv_list[[6]],
                                  pv_list[[7]],pv_list[[10]],pv_list[[11]],
                                  pv_list[[19]],pv_list[[23]])
full_df_batch$all_t <- as.factor(full_df_batch$all_t)
# replace batch names with capitalized batch names
full_df_batch$all_f <- sapply(full_df_batch$all_f, function(x) {
  new_f_name <- paste0('Factor ',strsplit(x,split='_')[[1]][[2]])
  return(new_f_name)
})
f_names_batch <- sapply(c(1,2,3,4,5,6,7,10,11,19,23),function(x){paste0('Factor ',x)})
full_df_batch$all_f <- factor(full_df_batch$all_f, levels=f_names_batch)

p <- ggplot(full_df_batch,aes(x = all_f, y = all_p, fill = all_t)) +
  geom_boxplot() +
  xlab('') +
  ylab('GC %') +
  ylim(.3,.8) +
  ggtitle('Batch associated factors') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=16),
        axis.text=element_text(size=10),
        axis.title=element_text(size=14))

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_batch_gc.pdf", useDingbats = FALSE,
    width = 9.5, height = 4)
p
dev.off()


all_pv_adj[c(8,9,14,16,20,21,25)]
full_df_no_batch_small <- rbind.data.frame(pv_list[[8]],pv_list[[9]],pv_list[[14]],pv_list[[16]],
                                  pv_list[[20]],pv_list[[21]],pv_list[[25]])
full_df_no_batch_small$all_t <- as.factor(full_df_no_batch_small$all_t)
full_df_no_batch_small$all_f <- sapply(full_df_no_batch_small$all_f, function(x) {
  new_f_name <- paste0('Factor ',strsplit(x,split='_')[[1]][[2]])
  return(new_f_name)
})
f_names_no_batch_small <- sapply(c(8,9,14,16,20,21,25),function(x){paste0('Factor ',x)})
full_df_no_batch_small$all_f <- factor(full_df_no_batch_small$all_f, levels=f_names_no_batch_small)

p <- ggplot(full_df_no_batch_small,aes(x = all_f, y = all_p, fill = all_t)) +
  geom_boxplot() +
  xlab('') +
  ylab('GC %') +
  ylim(.3,.8) +
  ggtitle('Non-batch associated factors') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=16),
        axis.text=element_text(size=10),
        axis.title=element_text(size=14))

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_no_batch_gc.pdf", useDingbats = FALSE,
    width = 6.9, height = 4)
p
dev.off()








## computing loadings correlations and pairwise jaccard coefficients
cormat_dmat_rot <- cor(t(pbmc_container$tucker_results[[2]]))
cormat_dmat_rot
colnames(cormat_dmat_rot) <- sapply(1:ncol(cormat_dmat_rot),function(x){paste0('Factor ',x)})
rownames(cormat_dmat_rot) <- sapply(1:nrow(cormat_dmat_rot),function(x){paste0('Factor ',x)})


lds_hmap <- Heatmap(cormat_dmat_rot, name = "Pearson r",
                    cluster_columns = TRUE,
                    cluster_rows = TRUE,
                    column_names_gp = gpar(fontsize = 10),
                    row_names_gp = gpar(fontsize = 10),
                    col = col_fun,border=TRUE, show_column_names=TRUE,
                    show_row_names=TRUE,show_row_dend = FALSE,
                    show_column_dend = FALSE, row_names_side = 'left',
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid::grid.text(sprintf("%.2f", cormat_dmat_rot[i, j]), x, y, gp = gpar(fontsize = 10))
                    })

pdf(file = "/home/jmitchel/figures/for_paper_v2/batch_lds_cor.pdf", useDingbats = FALSE,
    width = 10, height = 10)
lds_hmap
dev.off()


# c(1,2,3,4,5,6,7,10,11,19,23) batch factors
myfactor=2
tmp <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- tmp[[1]]
lds <- tmp[[2]]

sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

# limit to just the genes in tmp_casted_num
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df1 <- c(sig_df)

myfactor=6
tmp <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- tmp[[1]]
lds <- tmp[[2]]

sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

# limit to just the genes in tmp_casted_num
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df2 <- c(sig_df)


# binarize significance
sbin1 <- sig_df1<0.05
sbin2 <- sig_df2<0.05

tmp <- cbind(sbin1,sbin2)
sum(rowSums(tmp)==2)


get_sig_vec <- function(container,myfactor) {
  tmp <- get_one_factor(container,factor_select=myfactor)
  dsc <- tmp[[1]]
  lds <- tmp[[2]]
  
  sig_vectors <- get_significance_vectors(container,
                                          factor_select=myfactor, colnames(lds))
  # convert list to df
  sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
  
  # limit to just the genes in tmp_casted_num
  sig_df <- sig_df[rownames(lds),colnames(lds)]
  sig_df1 <- c(sig_df)
  return(sig_df1)
}

fact_use <- c(1,2,3,4,5,6,7,10,11,19,23)
n_fact <- length(fact_use)
myres=matrix(ncol=n_fact,nrow=n_fact)
for (i in 1:n_fact) {
  f_i <- fact_use[i]
  sv1 <- get_sig_vec(pbmc_container,f_i)
  for (j in 1:n_fact) {
    f_j <- fact_use[j]
    
    print(f_i)
    print(f_j)
    
    sv2 <- get_sig_vec(pbmc_container,f_j)
    
    # binarize significance
    sbin1 <- sv1<0.01
    sbin2 <- sv2<0.01
    
    tmp <- cbind(sbin1,sbin2)
    num_both <- sum(rowSums(tmp)==2)
    jacc <- num_both/sum(rowSums(tmp)>0)
    print(jacc)
    print('')
    
    # store res
    myres[i,j] <- jacc
  }
}

colnames(myres) <- sapply(fact_use,function(x){paste0('Factor ',x)})
rownames(myres) <- sapply(fact_use,function(x){paste0('Factor ',x)})

lds_hmap <- Heatmap(myres, name = "Pearson r",
                    cluster_columns = TRUE,
                    cluster_rows = TRUE,
                    column_names_gp = gpar(fontsize = 10),
                    row_names_gp = gpar(fontsize = 10),
                    col = col_fun,border=TRUE, show_column_names=TRUE,
                    show_row_names=TRUE,show_row_dend = FALSE,
                    show_column_dend = FALSE, row_names_side = 'left',
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid::grid.text(sprintf("%.2f", myres[i, j]), x, y, gp = gpar(fontsize = 10))
                    })

pdf(file = "/home/jmitchel/figures/for_paper_v2/batch_jacc.pdf", useDingbats = FALSE,
    width = 8, height = 8)
lds_hmap
dev.off()


## computing loadings correlations for just the batch factors
cormat_dmat_rot <- cor(t(pbmc_container$tucker_results[[2]][fact_use,]))
cormat_dmat_rot
colnames(cormat_dmat_rot) <- sapply(fact_use,function(x){paste0('Factor ',x)})
rownames(cormat_dmat_rot) <- sapply(fact_use,function(x){paste0('Factor ',x)})


lds_hmap <- Heatmap(cormat_dmat_rot, name = "Pearson r",
                    cluster_columns = TRUE,
                    cluster_rows = TRUE,
                    column_names_gp = gpar(fontsize = 10),
                    row_names_gp = gpar(fontsize = 10),
                    col = col_fun,border=TRUE, show_column_names=TRUE,
                    show_row_names=TRUE,show_row_dend = FALSE,
                    show_column_dend = FALSE, row_names_side = 'left',
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid::grid.text(sprintf("%.2f", cormat_dmat_rot[i, j]), x, y, gp = gpar(fontsize = 10))
                    })

pdf(file = "/home/jmitchel/figures/for_paper_v2/batch_lds_cor.pdf", useDingbats = FALSE,
    width = 6.5, height = 5.5)
lds_hmap
dev.off()

