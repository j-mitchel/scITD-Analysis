
library(ComplexHeatmap)

# loading up the DE demuxlet data
load("/home/jmitchel/data/lupus_data/cell.type.diffexp.RData")

# get union of significant DE genes from any cell type
pbmc.expressed.res.sig <- rownames(pbmc.expressed.res)[pbmc.expressed.res$padj<.001]
nk.expressed.res.sig <- rownames(nk.expressed.res)[nk.expressed.res$padj<.001]
ncmono.expressed.res.sig <- rownames(ncmono.expressed.res)[ncmono.expressed.res$padj<.001]
cd8.expressed.res.sig <- rownames(cd8.expressed.res)[cd8.expressed.res$padj<.001]
cd4.expressed.res.sig <- rownames(cd4.expressed.res)[cd4.expressed.res$padj<.001]
B.expressed.res.sig <- rownames(cd19.expressed.res)[cd19.expressed.res$padj<.001]
cM.expressed.res.sig <- rownames(cd14.expressed.res)[cd14.expressed.res$padj<.001]
dc.expressed.res.sig <- rownames(dc.expressed.res)[dc.expressed.res$padj<.001]

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





### for doing just t4 and cM
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




