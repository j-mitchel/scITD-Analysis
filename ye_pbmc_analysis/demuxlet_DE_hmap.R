
library(ComplexHeatmap)

# loading up the DE demuxlet data
load("/home/jmitchel/data/lupus_data/cell.type.diffexp.RData")

# get union of significant DE genes from any cell type
pbmc.expressed.res.sig <- rownames(pbmc.expressed.res)[pbmc.expressed.res$padj<.05]
nk.expressed.res.sig <- rownames(nk.expressed.res)[nk.expressed.res$padj<.05]
ncmono.expressed.res.sig <- rownames(ncmono.expressed.res)[ncmono.expressed.res$padj<.05]
cd8.expressed.res.sig <- rownames(cd8.expressed.res)[cd8.expressed.res$padj<.05]
cd4.expressed.res.sig <- rownames(cd4.expressed.res)[cd4.expressed.res$padj<.05]
B.expressed.res.sig <- rownames(cd19.expressed.res)[cd19.expressed.res$padj<.05]
cM.expressed.res.sig <- rownames(cd14.expressed.res)[cd14.expressed.res$padj<.05]
dc.expressed.res.sig <- rownames(dc.expressed.res)[dc.expressed.res$padj<.05]

# get genes in all comparisons
de_list <- list(pbmc.expressed.res, nk.expressed.res, ncmono.expressed.res,
                cd8.expressed.res, cd4.expressed.res, cd19.expressed.res,
                cd14.expressed.res, dc.expressed.res)

# get intersection of genes in present in all
gene_intersect <- rownames(de_list[[1]])
for (i in 2:length(de_list)) {
  gene_intersect <- intersect(gene_intersect,rownames(de_list[[i]]))
}


# now get union of DE genes from any cell type that are in gene_intersect
de_union <- unique(c(pbmc.expressed.res.sig,
         nk.expressed.res.sig,
         ncmono.expressed.res.sig,
         cd8.expressed.res.sig,
         cd4.expressed.res.sig,
         B.expressed.res.sig,
         cM.expressed.res.sig,
         dc.expressed.res.sig))

de_union_in_all <- de_union[de_union %in% gene_intersect]

tmp <- as.data.frame(cbind(
  de_list[[1]][de_union_in_all,'log2FoldChange'],
  de_list[[2]][de_union_in_all,'log2FoldChange'],
  de_list[[3]][de_union_in_all,'log2FoldChange'],
  de_list[[4]][de_union_in_all,'log2FoldChange'],
  de_list[[5]][de_union_in_all,'log2FoldChange'],
  de_list[[6]][de_union_in_all,'log2FoldChange'],
  de_list[[7]][de_union_in_all,'log2FoldChange'],
  de_list[[8]][de_union_in_all,'log2FoldChange']
))

colnames(tmp) <- c('pbmc', 'nk', 'ncmono',
                   'cd8', 'cd4', 'B',
                   'cM', 'dc')










###### for doing just t4 and cM
# get genes in all comparisons
de_list <- list(cd4.expressed.res,
                cd14.expressed.res)

# get intersection of genes in present in all
gene_intersect <- intersect(rownames(de_list[[1]]),rownames(de_list[[2]]))

# now get union of DE genes from any cell type that are in gene_intersect
de_union <- unique(c(cd4.expressed.res.sig,
                     cM.expressed.res.sig))

de_union_in_all <- de_union[de_union %in% gene_intersect]

tmp <- as.data.frame(cbind(
  de_list[[1]][de_union_in_all,'log2FoldChange'],
  de_list[[2]][de_union_in_all,'log2FoldChange']
))

colnames(tmp) <- c('cd4',
                   'cM')

Heatmap(as.matrix(tmp))



## should only show color for those genes that are significant! showing
# the logFC for all genes makes it look like there is more similarity than there actually is
rownames(tmp) <- de_union_in_all

genes_not_sig_t4 <- !(de_union_in_all %in% cd4.expressed.res.sig)
tmp[genes_not_sig_t4,'cd4'] <- 0

genes_not_sig_cM <- !(de_union_in_all %in% cM.expressed.res.sig)
tmp[genes_not_sig_cM,'cM'] <- 0

# order rows as down in both, up in both, down cM only, up cM only, down T4 only, up T4 only

down_both <- rownames(tmp)[rowSums(tmp < 0) == 2]
up_both <- rownames(tmp)[rowSums(tmp > 0) == 2]

down_cM <- rownames(tmp)[tmp[,'cM']<0]
down_t4 <- rownames(tmp)[tmp[,'cd4']<0]
up_cM <- rownames(tmp)[tmp[,'cM']>0]
up_t4 <- rownames(tmp)[tmp[,'cd4']>0]

down_cM_only <- down_cM[!(down_cM %in% down_t4)]
down_cM_only <- down_cM_only[!(down_cM_only %in% up_t4)]

up_cM_only <- up_cM[!(up_cM %in% down_t4)]
up_cM_only <- up_cM_only[!(up_cM_only %in% up_t4)]

down_t4_only <- down_t4[!(down_t4 %in% down_cM)]
down_t4_only <- down_t4_only[!(down_t4_only %in% up_cM)]

up_t4_only <- up_t4[!(up_t4 %in% down_cM)]
up_t4_only <- up_t4_only[!(up_t4_only %in% up_cM)]

new_row_order <- c(down_both,up_both,down_cM_only,up_cM_only,down_t4_only,up_t4_only)

# add the rest of the genes
new_row_order <- c(new_row_order,rownames(tmp)[!(rownames(tmp) %in% new_row_order)])

tmp <- tmp[new_row_order,]

tmp <- as.matrix(tmp)
colnames(tmp) <- c("CD4+ T","cMonocytes")
myhmap <- Heatmap(tmp,
                  name = 'log2FC',
                  show_row_names=FALSE,
                  show_row_dend=FALSE,
                  show_column_dend=FALSE,
                  cluster_rows=FALSE,
                  cluster_columns=FALSE,
                  column_title = "Stim vs Unstim",
                  column_title_side = "top",
                  column_title_gp = gpar(fontsize = 20),
                  border = TRUE)

pdf(file = "/home/jmitchel/figures/for_paper/cd4_cM_stim_DEG.pdf", useDingbats = FALSE,
    width = 5, height = 7)
myhmap
dev.off()

# plot the umap with just these two cell types
embed <- read.table("/home/jmitchel/data/lupus_data/demuxlet_tsne.csv",
                    row.names = 1, sep=',',header = TRUE)

embed <- embed[embed$multiplets=='singlet',]
dim(embed)
embed <- embed[embed$cell %in% c("CD4 T cells","CD14+ Monocytes"),]

# create embedding plot
library(ggplot2)
ggplot(embed,aes(x=tsne1,y=tsne2,color=cell)) +
  geom_point()

myembed <- ggplot(embed,aes(x=tsne1,y=tsne2,color=stim)) +
  geom_point() +
  labs(color = "Group") +
  theme_bw()

pdf(file = "/home/jmitchel/figures/for_paper/cd4_cM_stim_embed.pdf", useDingbats = FALSE,
    width = 6, height = 4)
myembed
dev.off()









### now creating a heatmap showing pseudobulked expression for each donor 
# in each cell type and in each condition

# first get genes tested for DE in both ctypes
load("/home/jmitchel/data/lupus_data/cell.type.diffexp.RData")

cd4.expressed.res_old <- cd4.expressed.res
cd14.expressed.res_old <- cd14.expressed.res

# cd4.expressed.res.sig <- rownames(cd4.expressed.res)[cd4.expressed.res$padj<.01]
# cM.expressed.res.sig <- rownames(cd14.expressed.res)[cd14.expressed.res$padj<.01]
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
load("/home/jmitchel/data/lupus_data/demuxlet_scITD_comparison.RData") # from demuxlet_analysis.R

param_list <- initialize_params(ctypes_use = c("CD14+ Monocytes",
                                               "CD4 T cells"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(count_data=pbmc_all, meta_data=pbmc_meta,
                                     params=param_list,
                                     label_donor_sex = FALSE)

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=10,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=1,
                              scale_var = TRUE, var_scale_power = .85)
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=5, gene_min_cells=10,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=1,
                              scale_var = TRUE, var_scale_power = .5) #testing this with new method

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

# # try scaling data before combining
# pb_cM <- scale(pb_cM)
# pb_cd4 <- scale(pb_cd4)

# stack dataframes
pb_total <- rbind(pb_cM,pb_cd4)

# # scale genes to unit variance
# pb_total <- scale(pb_total)

pb_total <- t(pb_total)

pb_total <- as.matrix(pb_total)

# reduce it to just the significant DE genes (rows)
de_union_in_all <- de_union_in_all[de_union_in_all%in%rownames(pb_total)]
pb_total <- pb_total[de_union_in_all,]


## order rows as down in both, up in both, down cM only, up cM only, down T4 only, up T4 only
# genes_plot <- rownames(pb_total)

# convert gene symbol names within ct de genes
cd4.expressed.res <- cd4.expressed.res_old
cd14.expressed.res <- cd14.expressed.res_old

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
# color_lim <- max(pb_total)
col_fun = colorRamp2(c(-color_lim, 0, color_lim), c("blue", "white", "red"))
# col_fun = colorRamp2(c(0, color_lim), c("white", "red"))
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

pdf(file = "/home/jmitchel/figures/for_paper/demuxlet_de_by_donors_v2.pdf", useDingbats = FALSE,
    width = 4.75, height = 5.5)
pdf(file = "/home/jmitchel/figures/for_paper_v2/demuxlet_de_hmap.pdf", useDingbats = FALSE,
    width = 4.75, height = 5.5)
myhmap
dev.off()














