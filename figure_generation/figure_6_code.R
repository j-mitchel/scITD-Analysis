library(scITD)
library(ComplexHeatmap)

##### starting with PBMC dataset from van der Wijst et. al (2018)
### see preprocessing/vdw_preprocessing.R for code used to generate the follwing objects
# counts matrix
pbmc_counts <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_counts_v2.rds')

# meta data matrix
pbmc_meta <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_meta_v2.rds')

# ensembl to gene name conversions
feature.names <- readRDS('/home/jmitchel/data/van_der_wijst/genes.rds')

# change names of ctypes to match those from sle dataset
pbmc_meta$ctypes <- sapply(as.character(pbmc_meta$ctypes),function(x){
  if (x=='CD4+ T') {
    return('Th')
  } else if (x=='cMonocyte') {
    return('cMono')
  } else if (x=='CD8+ T') {
    return('Tc')
  } else {
    return(x)
  }
})

pbmc_meta$ctypes <- as.factor(pbmc_meta$ctypes)

# set up project parameters
param_list <- initialize_params(ctypes_use = c("Th", "Tc", "cMono", "CD56(dim) NK", "B"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(count_data=pbmc_counts, meta_data=pbmc_meta,
                                     gn_convert = feature.names, params=param_list,
                                     label_donor_sex = TRUE)

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = 1.5)

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,10),
                                 tucker_type = 'regular', rotation_type = 'ica_dsc')


# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','lanes'),
                                        stat_use='pval')

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('lanes'),
                                    show_donor_ids = TRUE,
                                    cluster_by_meta='lanes',
                                    add_meta_associations=TRUE)

# pdf(file = "/home/jmitchel/figures/for_paper_v2/pbmc_dscores2.pdf", useDingbats = FALSE,
#     width = 6, height = 6.5)
pbmc_container$plots$donor_matrix
# dev.off()


# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)

## get loadings plots (for paper)
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=.01,
                                           display_genes=F,
                                           gene_callouts = TRUE,
                                           callout_n_gene_per_ctype=5,
                                           show_var_explained = FALSE)


pdf(file = "/home/jmitchel/figures/for_paper_v2/pbmc_loadings_hbb_factor.pdf", useDingbats = FALSE,
    width = 4, height = 3.5)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["5"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["5"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()













##### now doing batch analysis with the SLE dataset
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

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5)

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(25,37),
                                 tucker_type = 'regular', rotation_type = 'ica_dsc')

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('pool'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('pool'),
                                    cluster_by_meta='pool',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='rsq')

# pdf(file = "/home/jmitchel/figures/for_paper_v2/lupus_batch_dscores.pdf", useDingbats = FALSE,
#     width = 7, height = 6.5)
pbmc_container$plots$donor_matrix
# dev.off()

# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)


# show loadings plots
pdf(file = "/home/jmitchel/figures/for_paper_v2/batch_lds_f1.pdf", useDingbats = FALSE,
    width = 3.5, height = 4.5)
# pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=6, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.05, display_genes=FALSE,
#                                       gene_callouts=F)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=F)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper_v2/batch_lds_f3.pdf", useDingbats = FALSE,
    width = 3.5, height = 4.5)
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=3, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=F)
dev.off()









##### computing loadings correlations
# plot is in supplemental figure s4 but is part of this analysis
cormat_dmat_rot <- cor(t(pbmc_container$tucker_results[[2]]))
# colnames(cormat_dmat_rot) <- sapply(1:ncol(cormat_dmat_rot),function(x){paste0('Factor ',x)})
# rownames(cormat_dmat_rot) <- sapply(1:nrow(cormat_dmat_rot),function(x){paste0('Factor ',x)})
colnames(cormat_hybrid_rot) <- as.character(1:ncol(cormat_hybrid_rot))
rownames(cormat_hybrid_rot) <- as.character(1:nrow(cormat_hybrid_rot))

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

# pdf(file = "/home/jmitchel/figures/for_paper_v2/batch_lds_cor.pdf", useDingbats = FALSE,
#     width = 10, height = 10)
lds_hmap
# dev.off()













##### soup analysis and gc content analysis
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

res1 <- test_soup_association(pbmc_container,1,'up','any_up',soupProf)


## plotting soup fractions for upregulated genes vs all others
p <- ggplot(res1[[2]],aes(x = as.factor(all_t), y = all_p)) +
  geom_boxplot(notch=TRUE) +
  xlab('') +
  ylab('log10(soup fraction)') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))
pdf(file = "/home/jmitchel/figures/for_paper_v2/soup_any_up.pdf", useDingbats = FALSE,
    width = 3.5, height = 4)
p
dev.off()

## plotting soup fraction for genes upreg in all ctypes versus those upreg only in some ctypes
res2 <- test_soup_association(pbmc_container,1,'up','all_v_some',soupProf)
p <- ggplot(res2[[2]],aes(x = as.factor(all_t), y = all_p)) +
  geom_boxplot(notch=TRUE) +
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
  } else if (comp_type=='lds') {
    g_keep <- intersect(rownames(lds),rownames(gen_dat))
    tmp <- cbind.data.frame(lds[g_keep,],gen_dat[g_keep,'gc'])
    colnames(tmp) <- c('lds_val','gc_content')
    lmres <- summary(lm(lds_val~gc_content,data=tmp))
    pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
    print(pval)
    
    return(list(pval,1))
  }
}

pv <- test_gc_association(pbmc_container,my_factor=1,b_direc='up',comp_type='any_up') #for testing the fn

all_pv <- c()
pv_list <- list()
my_direcs <- c('up','up','up','up','down','down','down','up','up','up',
               'up','up','up','up','up','up','down','up','up','up','up',
               'up','up','up','up')
# my_direcs <- rep('up',ncol(pbmc_container[["tucker_results"]][[1]]))
# my_direcs <- rep('down',ncol(pbmc_container[["tucker_results"]][[1]]))
for (i in 1:ncol(pbmc_container[["tucker_results"]][[1]])) {
  print(i)
  pv <- test_gc_association(pbmc_container,my_factor=i,b_direc=my_direcs[i],comp_type='any_up')
  all_pv <- c(all_pv,pv[[1]])
  pv_list[[i]] <- pv[[2]]
}
all_pv_adj <- p.adjust(all_pv,method='fdr')



# plotting gc associations for batch-associated factors
all_pv_adj[c(1,2,3,4,5,6,7,8,15,17,21)]
full_df_batch <- rbind.data.frame(pv_list[[1]],pv_list[[2]],pv_list[[3]],
                                  pv_list[[4]],pv_list[[5]],pv_list[[6]],
                                  pv_list[[7]],pv_list[[8]],pv_list[[15]],
                                  pv_list[[17]],pv_list[[21]])
full_df_batch$all_t <- as.factor(full_df_batch$all_t)
# replace batch names with capitalized batch names
full_df_batch$all_f <- sapply(full_df_batch$all_f, function(x) {
  new_f_name <- paste0('Factor ',strsplit(x,split='_')[[1]][[2]])
  return(new_f_name)
})
f_names_batch <- sapply(c(1,2,3,4,5,6,7,8,15,17,21),function(x){paste0('Factor ',x)})
full_df_batch$all_f <- factor(full_df_batch$all_f, levels=f_names_batch)

p <- ggplot(full_df_batch,aes(x = all_f, y = all_p, fill = all_t)) +
  geom_boxplot(notch=TRUE) +
  xlab('') +
  ylab('GC %') +
  ylim(.3,.9) +
  ggtitle('Batch associated factors') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=16),
        axis.text=element_text(size=10),
        axis.title=element_text(size=14))

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_batch_gc.pdf", useDingbats = FALSE,
    width = 14.5, height = 4)
p
dev.off()


# plotting gc associations for non-batch-associated factors
all_pv_adj[c(9,11,13,14,16,18,22,23,24)]
full_df_no_batch_small <- rbind.data.frame(pv_list[[9]],pv_list[[11]],pv_list[[13]],pv_list[[14]],
                                           pv_list[[16]],pv_list[[18]],pv_list[[22]],pv_list[[23]],pv_list[[24]])
full_df_no_batch_small$all_t <- as.factor(full_df_no_batch_small$all_t)
full_df_no_batch_small$all_f <- sapply(full_df_no_batch_small$all_f, function(x) {
  new_f_name <- paste0('Factor ',strsplit(x,split='_')[[1]][[2]])
  return(new_f_name)
})
f_names_no_batch_small <- sapply(c(9,11,13,14,16,18,22,23,24),function(x){paste0('Factor ',x)})
full_df_no_batch_small$all_f <- factor(full_df_no_batch_small$all_f, levels=f_names_no_batch_small)

p <- ggplot(full_df_no_batch_small,aes(x = all_f, y = all_p, fill = all_t)) +
  geom_boxplot(notch=TRUE) +
  xlab('') +
  ylab('GC %') +
  ylim(.3,.9) +
  ggtitle('Non-batch associated factors') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5,size=16),
        axis.text=element_text(size=10),
        axis.title=element_text(size=14))

pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_no_batch_gc.pdf", useDingbats = FALSE,
    width = 11.75, height = 4)
p
dev.off()









## now evaluating remaining gc-content associations after applying batch correction
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

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(25,37),
                                 tucker_type = 'regular', rotation_type = 'ica_dsc')

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('pool'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('pool'),
                                    cluster_by_meta='pool',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='rsq')

pbmc_container$plots$donor_matrix

# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)

## testing gc content associations for the factors in either direction
my_direcs <- rep('up',ncol(pbmc_container[["tucker_results"]][[1]]))
for (i in 1:ncol(pbmc_container[["tucker_results"]][[1]])) {
  print(i)
  pv <- test_gc_association(pbmc_container,my_factor=i,b_direc=my_direcs[i],comp_type='any_up')
  all_pv <- c(all_pv,pv[[1]])
  pv_list[[i]] <- pv[[2]]
}
all_pv_adj <- p.adjust(all_pv,method='fdr')
print(all_pv_adj)

my_direcs <- rep('down',ncol(pbmc_container[["tucker_results"]][[1]]))
for (i in 1:ncol(pbmc_container[["tucker_results"]][[1]])) {
  print(i)
  pv <- test_gc_association(pbmc_container,my_factor=i,b_direc=my_direcs[i],comp_type='any_up')
  all_pv <- c(all_pv,pv[[1]])
  pv_list[[i]] <- pv[[2]]
}
all_pv_adj <- p.adjust(all_pv,method='fdr')
print(all_pv_adj)











