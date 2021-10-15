library(scITD)
library(ComplexHeatmap)
library(ggplot2)
library(Seurat)
library(circlize)

##### plot the umap with just the two cell types

# data from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583 batch2
embed <- read.table("/home/jmitchel/data/lupus_data/demuxlet_tsne.csv",
                    row.names = 1, sep=',',header = TRUE)

embed <- embed[embed$multiplets=='singlet',]
embed <- embed[embed$cell %in% c("CD4 T cells","CD14+ Monocytes"),]

# create embedding plot
myembed <- ggplot(embed,aes(x=tsne1,y=tsne2,color=stim)) +
  geom_point() +
  labs(color = "Group") +
  theme_bw()

# pdf(file = "/home/jmitchel/figures/for_paper/cd4_cM_stim_embed.pdf", useDingbats = FALSE,
#     width = 6, height = 4)
myembed
# dev.off()

rm(list=ls())






##### prep counts data for additional analyses
# same embedding loaded as above
embed <- read.table("/home/jmitchel/data/lupus_data/demuxlet_tsne.csv",
                    row.names = 1, sep=',',header = TRUE)

# limit embed to just singlets
embed <- embed[embed$multiplets=='singlet',]

# data from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96583
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

# need to make each donor + stim combination a separate "donor"
pbmc_meta$donors <- sapply(1:nrow(pbmc_meta), function(i) {
  paste0(pbmc_meta[i,'donors'],"_",pbmc_meta[i,'stim'])
})
pbmc_meta$donors <- as.factor(pbmc_meta$donors)










##### generate DE heatmap

# data from: https://github.com/yelabucsf/demuxlet_paper_code/tree/master/s15_interferon_comparison
load("/home/jmitchel/data/lupus_data/cell.type.diffexp.RData")

cd4.expressed.res.sig <- rownames(cd4.expressed.res)[cd4.expressed.res$padj<.01]
cM.expressed.res.sig <- rownames(cd14.expressed.res)[cd14.expressed.res$padj<.01]
de_list <- list(cd4.expressed.res,
                cd14.expressed.res)

# get intersection of genes in present in all
gene_intersect <- intersect(rownames(de_list[[1]]),rownames(de_list[[2]]))

# now get union of DE genes from any cell type that are in gene_intersect
de_union <- unique(c(cd4.expressed.res.sig,
                     cM.expressed.res.sig))

de_union_in_all <- de_union[de_union %in% gene_intersect]

# convert to gene symbols
de_union_in_all <- de_list[[1]][de_union_in_all,'featureData.symbol']
de_union_in_all <- as.character(de_union_in_all)

# get pseudobulk data 
param_list <- initialize_params(ctypes_use = c("CD14+ Monocytes",
                                               "CD4 T cells"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(count_data=pbmc_all, meta_data=pbmc_meta,
                                     params=param_list,
                                     label_donor_sex = FALSE)

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=1,
                              scale_var = TRUE, var_scale_power = .85)

pb_cM <- pbmc_container[["scMinimal_ctype"]][["CD14+ Monocytes"]][["pseudobulk"]]
pb_cd4 <- pbmc_container[["scMinimal_ctype"]][["CD4 T cells"]][["pseudobulk"]]

# store old donor symbols
old_d_symbols <- c(rownames(pb_cM),rownames(pb_cd4))

# append ctype name to the donor symbols
rownames(pb_cM) <- sapply(rownames(pb_cM),function(x) {
  paste0(x,"_cM")
})
names(rownames(pb_cM)) <- NULL

rownames(pb_cd4) <- sapply(rownames(pb_cd4),function(x) {
  paste0(x,"_cd4")
})
names(rownames(pb_cd4)) <- NULL

# stack dataframes
pb_total <- rbind(pb_cM,pb_cd4)

pb_total <- t(pb_total)

pb_total <- as.matrix(pb_total)

# reduce it to just the significant DE genes (rows)
de_union_in_all <- de_union_in_all[de_union_in_all%in%rownames(pb_total)]
pb_total <- pb_total[de_union_in_all,]


## order rows as down in both, up in both, down cM only, up cM only, down T4 only, up T4 only
cd4.sig <- cd4.expressed.res[cd4.expressed.res.sig,]
cd4.sig.up <- cd4.sig[cd4.sig$log2FoldChange>0,]
cd4.sig.up.genes <- as.character(cd4.sig.up$featureData.symbol)
cd4.sig.down <- cd4.sig[cd4.sig$log2FoldChange<0,]
cd4.sig.down.genes <- as.character(cd4.sig.down$featureData.symbol)

cM.sig <- cd14.expressed.res[cM.expressed.res.sig,]
cM.sig.up <- cM.sig[cM.sig$log2FoldChange>0,]
cM.sig.up.genes <- as.character(cM.sig.up$featureData.symbol)
cM.sig.down <- cM.sig[cM.sig$log2FoldChange<0,]
cM.sig.down.genes <- as.character(cM.sig.down$featureData.symbol)

de_both_up <- intersect(cd4.sig.up.genes,cM.sig.up.genes)
de_both_down <- intersect(cd4.sig.down.genes,cM.sig.down.genes)
cM_only_up <- cM.sig.up.genes[!(cM.sig.up.genes %in% de_both_up)]
cM_only_down <- cM.sig.down.genes[!(cM.sig.down.genes %in% de_both_down)]
cd4_only_up <- cd4.sig.up.genes[!(cd4.sig.up.genes %in% de_both_up)]
cd4_only_down <- cd4.sig.down.genes[!(cd4.sig.down.genes %in% de_both_down)]

gene_order <- unique(c(de_both_up,de_both_down,cM_only_up,cM_only_down,cd4_only_up,cd4_only_down))
gene_order <- gene_order[gene_order %in% rownames(pb_total)]
pb_total <- pb_total[gene_order,]

# get gene callout annotation for hmap
gene_callouts <- c('MX1','CCL3',
                   'ANXA5','DDIT3',
                   'MYC')

ndx <- match(gene_callouts,rownames(pb_total))

myannot <- rowAnnotation(callouts = anno_mark(at = ndx, which='row',
                                              labels = gene_callouts))

# make heatmap
ta = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:3), labels = c('cMonocytes','CD4+ T')))
color_lim <- stats::quantile(as.matrix(abs(pb_total)), c(.95))
col_fun = colorRamp2(c(-color_lim, 0, color_lim), c("blue", "white", "red"))
myhmap <- Heatmap(pb_total,
                  name = 'expression',
                  show_row_names=FALSE,
                  show_column_names=TRUE,
                  show_row_dend=FALSE,
                  show_column_dend=FALSE,
                  cluster_rows=FALSE,
                  cluster_columns=FALSE,
                  column_title = "Stim vs Unstim",
                  column_title_side = "top",
                  column_title_gp = gpar(fontsize = 20),
                  border = TRUE,
                  column_split = factor(c(rep('cMonocytes',16),rep('CD4+ T',16)),levels=c('cMonocytes','CD4+ T')),
                  top_annotation = ta,
                  col = col_fun,
                  right_annotation = myannot)

# pdf(file = "/home/jmitchel/figures/for_paper_v2/demuxlet_de_hmap.pdf", useDingbats = FALSE,
#     width = 4.75, height = 5.5)
myhmap
# dev.off()











##### now running scITD to compare to DE results
param_list <- initialize_params(ctypes_use = c("CD14+ Monocytes",
                                               "CD4 T cells"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(count_data=pbmc_all, meta_data=pbmc_meta,
                                     params=param_list,
                                     label_donor_sex = FALSE)

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, 
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var', vargenes_thresh=1500,
                              scale_var = TRUE, var_scale_power = .5)

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(2,4),
                                 tucker_type = 'regular', rotation_type = 'ica_dsc')

# plotting just the first factor in the sample scores matrix
pbmc_container$tucker_results[[1]] <- pbmc_container$tucker_results[[1]][,1,drop=FALSE]
pbmc_container[["exp_var"]] <- pbmc_container[["exp_var"]][1]

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('stim'),stat_use='pval')

pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = TRUE)
# pdf(file = "/home/jmitchel/figures/for_paper_v2/demuxlet_scores.pdf", useDingbats = FALSE,
#     width = 3, height = 6)
pbmc_container$plots$donor_matrix
# dev.off()

# rerun tucker to keep both factors
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(2,4),
                                 tucker_type = 'regular', rotation_type = 'ica_dsc')

# get significant genes
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
                                                          'MYC'),
                                      show_var_explained = FALSE,
                                      clust_method = 'ward.D')

# pdf(file = "/home/jmitchel/figures/for_paper_v2/demuxlet_loadings_v2.pdf", useDingbats = FALSE,
#     width = 5.25, height = 6)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["1"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["1"]],
     legend_grouping = "original",
     newpage=TRUE)
# dev.off()







##### calculating correlation between fold change and loadings for all gene_ctype combos present in both
tmp_casted_num <- get_one_factor(pbmc_container,1)[[2]]

# get matrix of gene-association p-values
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=1, colnames(tmp_casted_num))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

# order df same way as in tmp_casted_num
sig_df <- sig_df[rownames(tmp_casted_num),colnames(tmp_casted_num)]
##

# remove duplicate gene in cd4.expressed.res
ndx_check <- which(duplicated(cd4.expressed.res$featureData.symbol))
ndx_check <- which(cd4.expressed.res$featureData.symbol=='NAA60')
cd4.expressed.res <- cd4.expressed.res[-ndx_check[2],]

# set rownames to be gene names, not ensemble symbols
rownames(cd4.expressed.res) <- cd4.expressed.res$featureData.symbol
rownames(cd14.expressed.res) <- cd14.expressed.res$featureData.symbol

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

# pdf(file = "/home/jmitchel/figures/for_paper_v2/demuxlet_loading_vs_log2FC_v2.pdf", useDingbats = FALSE,
#     width = 4, height = 1.5)
ggplot(tmp,aes(x=log2FC,y=loading,color=sig_val)) +
  geom_point(alpha = .1, pch=19) +
  theme_bw() +
  scale_color_manual(values=c("#000000", "#FF0000")) +
  guides(color = 'none')
dev.off()









##### calculating correlation of between scITD pvalues and DE pvalues
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

tmp <- as.data.frame(cbind(lds_both,fc_both))
colnames(tmp) <- c('scITD_padj','DE_padj')
rownames(tmp) <- gene_ct
tmp$scITD_padj[tmp$scITD_padj<.0001] <- .0001
tmp$DE_padj[tmp$DE_padj<.000000000000001] <- .000000000000001

# pdf(file = "/home/jmitchel/figures/for_paper_v2/demuxlet_lm_vs_DEpval_v2.pdf", useDingbats = FALSE,
#     width = 4.5, height = 3)
ggplot(tmp,aes(x=-log(DE_padj,base=10),y=-log(scITD_padj,base=10))) +
  geom_point(alpha = 0.3,pch=19) +
  theme_bw() +
  geom_hline(yintercept=-log10(.05), linetype="dashed", color = "red") +
  geom_vline(xintercept=-log10(.05), linetype="dashed", color = "red") +
  xlab('DE -log10(adj p-value)') +
  ylab('F1 gene significance\n-log10(adj p-value)')
# dev.off()

# calculate number of genes above horizontal line and below vertical line
mask <- tmp
mask <- -log10(mask)
pv_thresh=.05
sum(mask$DE_padj<(-log10(pv_thresh)) & mask$scITD_padj>(-log10(pv_thresh)))

# now get vice versa
sum(mask$DE_padj>(-log10(pv_thresh)) & mask$scITD_padj<(-log10(pv_thresh)))

# get DE both
sum(mask$DE_padj>(-log10(pv_thresh)) & mask$scITD_padj>(-log10(pv_thresh)))















##### doing subsampling and recalculation of loadings-FC correlation (part of figure s1 but uses above data)
cor_helper <- function(tmp_casted_num,cd4.expressed.res,cd14.expressed.res) {
  lds_both <- c()
  fc_both <- c()
  for (i in 1:nrow(tmp_casted_num)) {
    gene <- rownames(tmp_casted_num)[i]
    for (j in 1:ncol(tmp_casted_num)) {
      ctype <- colnames(tmp_casted_num)[j]
      if (ctype=='CD4 T cells') {
        if (gene %in% rownames(cd4.expressed.res)) {
          demux_fc <- cd4.expressed.res[gene,'log2FoldChange']
          fc_both <- c(fc_both,demux_fc)
          lds_both <- c(lds_both,tmp_casted_num[i,j])
        }
      } else if (ctype=='CD14+ Monocytes') {
        if (gene %in% rownames(cd14.expressed.res)) {
          demux_fc <- cd14.expressed.res[gene,'log2FoldChange']
          fc_both <- c(fc_both,demux_fc)
          lds_both <- c(lds_both,tmp_casted_num[i,j])
        }
      }
    }
  }
  
  mycor <- cor(lds_both,fc_both,method = 'spearman')
  return(mycor)
}


num_donors <- nrow(pbmc_container$scMinimal_ctype[[1]]$pseudobulk)

cells_per_donor <- table(pbmc_container$scMinimal_full$metadata[,c('donors','ctypes')])

combined_meta <- pbmc_container$scMinimal_full$metadata[,c('donors','ctypes')]
combined_meta <- combined_meta[combined_meta$ctypes %in% c('CD4 T cells','CD14+ Monocytes'),]
combined_meta$ctypes <- factor(combined_meta$ctypes,levels=unique(combined_meta$ctypes))

full_counts <- pbmc_container$scMinimal_full$count_data

sizes_test <- c(12,20,40,60,80,100,120)
downsample_sizes <- num_donors * 2 * sizes_test

all_cors <- c()
cpd <- c()
for (ds in downsample_sizes) {
  # subsample data to correct mean size 
  prev_subs <- list()
  prev_ident <- FALSE
  for (myiter in 1:5) {
    min_cpd <- 0
    while (min_cpd < 5 | prev_ident) {
      prev_ident <- FALSE
      cells_sampled <- sample(rownames(combined_meta),ds)
      meta_sub <- combined_meta[cells_sampled,]
      cells_per_donor <- table(meta_sub[,c('donors','ctypes')])
      min_cpd <- min(cells_per_donor)
      
      # determine whether subsampling was previously recorded
      if (length(prev_subs)>1) {
        for (j in 1:length(prev_subs)) {
          if (identical(prev_subs[[j]],cells_sampled)) {
            prev_ident <- TRUE
          }
        }
      }
    }
    prev_subs[[length(prev_subs)+1]] <- cells_sampled
    
    counts_sub <- full_counts[,cells_sampled]
    
    param_list <- initialize_params(ctypes_use = c("CD14+ Monocytes",
                                                   "CD4 T cells"),
                                    ncores = 30, rand_seed = 10)
    
    pbmc_container <- make_new_container(count_data=counts_sub, meta_data=meta_sub,
                                         params=param_list,
                                         label_donor_sex = FALSE)
    
    pbmc_container <- form_tensor(pbmc_container, donor_min_cells=0, 
                                  norm_method='trim', scale_factor=10000,
                                  vargenes_method='norm_var', vargenes_thresh=1500,
                                  scale_var = TRUE, var_scale_power = .5)
    
    pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(2,4,2),
                                     tucker_type = 'regular', rotation_type = 'ica_dsc')
    
    ## compute loadings-FC correlations
    f1 <- get_one_factor(pbmc_container,1)
    f1_lds <- f1[[2]]
    
    my_cor <- cor_helper(f1_lds,cd4.expressed.res,cd14.expressed.res)
    
    # store results
    all_cors <- c(all_cors,my_cor)
    cpd <- c(cpd,ds)
  }
}


# plot results
tmp <- cbind.data.frame(all_cors,cpd,cpd/(num_donors * 2))
colnames(tmp) <- c('spearman_cor','total_num_cells','cells_donor_ctype')
tmp$spearman_cor <- abs(tmp$spearman_cor)

# get mean values for each factor at each value of cells_per
tmp2_means <- c()
tmp2_cells_per <- c()
for (cp in unique(tmp$cells_donor_ctype)) {
  tmp_sub <- tmp[tmp$cells_donor_ctype==cp,]
  tmp2_means <- c(tmp2_means,mean(abs(tmp_sub$spearman_cor)))
  tmp2_cells_per <- c(tmp2_cells_per,cp)
}

tmp2 <- cbind.data.frame(tmp2_means,tmp2_cells_per)
colnames(tmp2) <- c('spearman_cor','cells_donor_ctype')

tmp <- tmp[tmp$cells_donor_ctype!=12,]
tmp2 <- tmp2[tmp2$cells_donor_ctype!=12,]

# pdf(file = "/home/jmitchel/figures/for_paper_v2/demux_loading_cor_subsamp3.pdf", useDingbats = FALSE,
#     width = 4, height = 2.25)
ggplot(tmp,aes(x=cells_donor_ctype,y=spearman_cor)) +
  geom_point() +
  geom_line(data=tmp2,aes(x=cells_donor_ctype,y=spearman_cor)) +
  xlab('Av. cells per donor_cell type') +
  ylab('Loadings-log2FC cor') +
  theme_bw()
# dev.off()





