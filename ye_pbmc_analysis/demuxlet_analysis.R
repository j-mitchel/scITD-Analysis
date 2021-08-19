library(Seurat)

# load cell meta data with cell types and donors
embed <- read.table("/home/jmitchel/data/lupus_data/demuxlet_tsne.csv",
                    row.names = 1, sep=',',header = TRUE)

# limit embed to just singlets
embed <- embed[embed$multiplets=='singlet',]

pbmc_stim <- Read10X(data.dir = "/home/jmitchel/data/demuxlet/stim")
pbmc_ctrl <- Read10X(data.dir = "/home/jmitchel/data/demuxlet/ctrl")

# rename barcodes for duplicates between stim and ctrl so that it matches the tsne file
bcode_dup_mask <- colnames(pbmc_stim) %in% colnames(pbmc_ctrl)
colnames(pbmc_stim)[bcode_dup_mask] <- sapply(colnames(pbmc_stim)[bcode_dup_mask], function(x) {
  paste0(x,'1')
})

# limit each counts matrix to only singlets kept in embed
stim_keep <- colnames(pbmc_stim) %in% rownames(embed)
pbmc_stim <- pbmc_stim[,stim_keep]

ctrl_keep <- colnames(pbmc_ctrl) %in% rownames(embed)
pbmc_ctrl <- pbmc_ctrl[,ctrl_keep]


# combine counts matrices
pbmc_all <- cbind(pbmc_ctrl,pbmc_stim)

pbmc_meta <- embed
colnames(pbmc_meta)[c(3,6)] <- c('donors','ctypes')
pbmc_meta$ctypes <- as.factor(pbmc_meta$ctypes)
pbmc_meta$donors <- as.factor(pbmc_meta$donors)
pbmc_meta$stim <- as.factor(pbmc_meta$stim)
class(pbmc_meta$donors)
class(pbmc_meta$ctypes)
class(pbmc_meta$stim)

# need to make each donor + stim combination a separate "donor"
pbmc_meta$donors <- sapply(1:nrow(pbmc_meta), function(i) {
  paste0(pbmc_meta[i,'donors'],"_",pbmc_meta[i,'stim'])
})
pbmc_meta$donors <- as.factor(pbmc_meta$donors)

# show cell counts of cell types available
print(table(pbmc_meta$ctypes))

# set up project parameters
param_list <- initialize_params(ctypes_use = c("B cells","CD14+ Monocytes",
                                               "CD4 T cells","CD8 T cells",
                                               "FCGR3A+ Monocytes","NK cells"),
                                ncores = 30, rand_seed = 10)
param_list <- initialize_params(ctypes_use = c("CD14+ Monocytes",
                                               "CD4 T cells"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(count_data=pbmc_all, meta_data=pbmc_meta,
                                     params=param_list,
                                     label_donor_sex = FALSE)

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=10,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var', vargenes_thresh=1500,
                              scale_var = TRUE, var_scale_power = .5)

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(3,6,6),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(1,6,6),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(2,4,2),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(2,8,2),
                                 tucker_type = 'regular', rotation_type = 'ica')

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('stim'),stat_use='pval')

pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('stim'),
                                    cluster_by_meta = 'stim',
                                    show_donor_ids = TRUE,
                                    add_meta_associations='pval')

pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 6, height = 7)
pbmc_container$plots$donor_matrix
dev.off()

# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(3,6,6), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(1,6,6), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(2,4,2), n_fibers=100, n_iter=2000,
                                tucker_type='regular', rotation_type='ica')
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(2,8,2), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')
pbmc_container <- run_jackstraw_v2(pbmc_container, ranks=c(2,8,2), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')
pbmc_container <- get_lm_pvals(pbmc_container)

# get loadings plots
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=1,
                                      use_sig_only=TRUE,
                                      nonsig_to_zero=TRUE,
                                      sig_thresh=.01,
                                      display_genes=FALSE,
                                      gene_callouts=TRUE,
                                      specific_callouts=c('MX1','CCL3',
                                                          'ANXA5','DDIT3',
                                                          'MYC'))

pdf(file = "/home/jmitchel/figures/for_paper_v2/demuxlet_loadings_v2.pdf", useDingbats = FALSE,
    width = 5.25, height = 6)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["1"]],
                         annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["1"]],
                         legend_grouping = "original",
                         newpage=TRUE)
dev.off()

# # try to reorganize the hmap by hand so it sort of matches the one from demuxlet paper
# lds <- pbmc_container[["plots"]][["lds_plots_data"]]




# compare results to de genes from demuxlet study

## first extract my loadings matrix and associated jackstraw p-values
ldngs <- pbmc_container$tucker_results[[2]]

# break down a factor from the loadings matrix
genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})

sr_col <- ldngs[1,]

tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)

sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=1, colnames(tmp_casted_num))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

# order df same way as in tmp_casted_num
sig_df <- sig_df[rownames(tmp_casted_num),colnames(tmp_casted_num)]
##




load("/home/jmitchel/data/lupus_data/cell.type.diffexp.RData")
# cd4.expressed.res
# cd14.expressed.res

# remove duplicate gene in cd4.expressed.res
ndx_check <- which(duplicated(cd4.expressed.res$featureData.symbol))
cd4.expressed.res[ndx_check,] # gene is NAA60
ndx_check <- which(cd4.expressed.res$featureData.symbol=='NAA60')
cd4.expressed.res[ndx_check[1],]
cd4.expressed.res[ndx_check[2],]
cd4.expressed.res <- cd4.expressed.res[-ndx_check[2],]
sum(duplicated(cd4.expressed.res$featureData.symbol))

# set rownames to be gene names, not ensemble symbols
rownames(cd4.expressed.res) <- cd4.expressed.res$featureData.symbol
rownames(cd14.expressed.res) <- cd14.expressed.res$featureData.symbol

sig_in_both <- 0
sig_tested_in_both <- 0
for (i in 1:nrow(sig_df)) {
  gene <- rownames(sig_df)[i]
  for (j in 1:ncol(sig_df)) {
    ctype <- colnames(sig_df)[j]
    pv <- sig_df[i,j]
    
    if (pv < 0.05) {
      if (ctype=='CD4 T Cells') {
        if (gene %in% rownames(cd4.expressed.res)) {
          demux_pv <- cd4.expressed.res[gene,'padj']
          if (demux_pv < 0.05) {
            sig_in_both <- sig_in_both + 1
          }
          sig_tested_in_both <- sig_tested_in_both + 1
        }
      } else if (ctype=='CD14+ Monocytes') {
        if (gene %in% rownames(cd14.expressed.res)) {
          demux_pv <- cd14.expressed.res[gene,'padj']
          if (demux_pv < 0.05) {
            sig_in_both <- sig_in_both + 1
          }
          sig_tested_in_both <- sig_tested_in_both + 1
        }
      }
    }
  }
}

sig_in_both
sig_tested_in_both


# save.image(file = "/home/jmitchel/data/lupus_data/demuxlet_scITD_comparison.RData")
# load("/home/jmitchel/data/lupus_data/demuxlet_scITD_comparison.RData")


# plotting just the first factor in the sample scores matrix
pbmc_container$tucker_results[[1]] <- pbmc_container$tucker_results[[1]][,1,drop=FALSE]
pbmc_container[["exp_var"]] <- pbmc_container[["exp_var"]][1]

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('stim'),stat_use='pval')

pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = TRUE)
pdf(file = "/home/jmitchel/figures/for_paper_v2/demuxlet_scores.pdf", useDingbats = FALSE,
    width = 3, height = 6)
pbmc_container$plots$donor_matrix
dev.off()



## calculating correlation between fold change and loadings for all gene_ctype combos present in both
lds_both <- c()
fc_both <- c()
sig_both <- c()
for (i in 1:nrow(tmp_casted_num)) {
  gene <- rownames(tmp_casted_num)[i]
  for (j in 1:ncol(tmp_casted_num)) {
    ctype <- colnames(tmp_casted_num)[j]
    if (ctype=='CD4 T cells') {
      if (gene %in% rownames(cd4.expressed.res)) {
        demux_fc <- cd4.expressed.res[gene,'log2FoldChange']
        fc_both <- c(fc_both,demux_fc)
        lds_both <- c(lds_both,tmp_casted_num[i,j])
        sig_both <- c(sig_both,sig_df[i,j])
      }
    } else if (ctype=='CD14+ Monocytes') {
      if (gene %in% rownames(cd14.expressed.res)) {
        demux_fc <- cd14.expressed.res[gene,'log2FoldChange']
        fc_both <- c(fc_both,demux_fc)
        lds_both <- c(lds_both,tmp_casted_num[i,j])
        sig_both <- c(sig_both,sig_df[i,j])
      }
    }
  }
}

cor(lds_both,fc_both,method = 'spearman')
tmp <- as.data.frame(cbind(lds_both,fc_both,sig_both))
colnames(tmp) <- c('loading','log2FC','sig_val')
tmp$sig_val <- tmp$sig_val<.01
pdf(file = "/home/jmitchel/figures/for_paper_v2/demuxlet_loading_vs_log2FC_v2.pdf", useDingbats = FALSE,
    width = 4, height = 1.5)
pdf(file = "/home/jmitchel/figures/for_paper_v2/demuxlet_loading_vs_log2FC_v2.pdf", useDingbats = FALSE,
    width = 4, height = 1.5)
ggplot(tmp,aes(x=log2FC,y=loading,color=sig_val)) +
  geom_point(alpha = .1, pch=19) +
  theme_bw() +
  scale_color_manual(values=c("#000000", "#FF0000")) +
  guides(color = FALSE)
dev.off()


# showing correlation of pvalues
lds_both <- c()
fc_both <- c()
gene_ct <- c()
for (i in 1:nrow(tmp_casted_num)) {
  gene <- rownames(tmp_casted_num)[i]
  for (j in 1:ncol(tmp_casted_num)) {
    ctype <- colnames(tmp_casted_num)[j]
    if (ctype=='CD4 T cells') {
      if (gene %in% rownames(cd4.expressed.res)) {
        demux_fc <- cd4.expressed.res[gene,'padj']
        fc_both <- c(fc_both,demux_fc)
        lds_both <- c(lds_both,sig_df[i,j])
        gene_ct <- c(gene_ct,paste0(gene,"_",ctype))
      }
    } else if (ctype=='CD14+ Monocytes') {
      if (gene %in% rownames(cd14.expressed.res)) {
        demux_fc <- cd14.expressed.res[gene,'padj']
        fc_both <- c(fc_both,demux_fc)
        lds_both <- c(lds_both,sig_df[i,j])
        gene_ct <- c(gene_ct,paste0(gene,"_",ctype))
      }
    }
  }
}

cor(lds_both,fc_both,method = 'spearman')
tmp <- as.data.frame(cbind(lds_both,fc_both))
colnames(tmp) <- c('scITD_padj','DE_padj')
rownames(tmp) <- gene_ct
tmp$scITD_padj[tmp$scITD_padj<.0001] <- .0001
tmp$DE_padj[tmp$DE_padj<.000000000000001] <- .000000000000001

tmp$scITD_padj[tmp$scITD_padj<.000001] <- .000001
tmp$DE_padj[tmp$DE_padj<.0000000000001] <- .0000000000001

pdf(file = "/home/jmitchel/figures/for_paper_v2/demuxlet_lm_vs_DEpval_v2.pdf", useDingbats = FALSE,
    width = 4.5, height = 3)
ggplot(tmp,aes(x=-log(DE_padj,base=10),y=-log(scITD_padj,base=10))) +
  geom_point(alpha = 0.3,pch=19) +
  theme_bw() +
  geom_hline(yintercept=-log10(.05), linetype="dashed", color = "red") +
  geom_vline(xintercept=-log10(.05), linetype="dashed", color = "red") +
  xlab('DE -log10(adj p-value)') +
  ylab('F1 gene significance\n-log10(adj p-value)')
dev.off()

# calculate number of genes above horizontal line and below vertical line
mask <- tmp
mask <- -log10(mask)
pv_thresh=.05
sum(mask$DE_padj<(-log10(pv_thresh)) & mask$scITD_padj>(-log10(pv_thresh)))

# now get vice versa
sum(mask$DE_padj>(-log10(pv_thresh)) & mask$scITD_padj<(-log10(pv_thresh)))

# get DE both
sum(mask$DE_padj>(-log10(pv_thresh)) & mask$scITD_padj>(-log10(pv_thresh)))


# #### trying with no jackstraw gene significance 
# # get associations as gene.ct.factor
# pvals <- get_real_fstats(pbmc_container,ncores=4) # using pvals for this fn
# padj <- p.adjust(pvals,method='fdr')
# names(padj) <- sapply(names(padj),function(x) {
#   return(substr(x,1,nchar(x)-6))
# })
# 
# pbmc_container[["gene_score_associations"]] <- padj
# # now rerun plot above
# 
# pdf(file = "/home/jmitchel/figures/for_paper/demuxlet_no_jackstraw_vs_DEpval_v2.pdf", useDingbats = FALSE,
#     width = 4.5, height = 3.5)
# ggplot(tmp,aes(x=-log(DE_padj,base=10),y=-log(scITD_padj,base=10))) +
#   geom_point(alpha = 0.3) +
#   theme_bw() +
#   geom_hline(yintercept=-log10(.05), linetype="dashed", color = "red") +
#   geom_vline(xintercept=-log10(.05), linetype="dashed", color = "red") +
#   xlab('DE -log10(adj p-value)') +
#   ylab('scITD Jackstraw\n-log10(adj p-value)')
# dev.off()

