library(plyr)

## will first get the different decompositions and they include:
# ICA on dscores
# ICA on loadings
# hybrid method on loadings

ica_dmat <- pbmc_container$tucker_results
ica_lds <- pbmc_container$tucker_results
hybrid <- pbmc_container$tucker_results



# get factor associations with core set of IFN genes...
ifn_core <- c('HERC5', 'IFI27', 'IRF7', 'ISG15', 'LY6E', 'MX1', 'OAS2', 'OAS3',
              'RSAD2', 'USP18', 'GBP5')
go_ifn1 <- read.csv(file='/home/jmitchel/IFN_gene_list.csv')[,1]
go_ifn1 <- go_ifn1[go_ifn1 %in% colnames(pbmc_container[["scMinimal_ctype"]][[1]][["pseudobulk"]])]

plot_ifn_associations <- function(container,tuck_res,ifn_genes) {
  dscores <- tuck_res[[1]]
  all_factors <- c()
  rsq <- c()
  for (i in 1:ncol(dscores)) {
    for (gn in ifn_genes) {
      for (ct in container$experiment_params$ctypes_use) {
        expr <- container[["scMinimal_ctype"]][[ct]][["pseudobulk"]][,gn]
        tmp <- cbind.data.frame(expr,dscores[names(expr),i])
        colnames(tmp) <- c('g_exp','dsc')
        lmres <- lm(dsc~g_exp,data=tmp)
        lmres <- summary(lmres)
        rsq <- c(rsq,lmres$r.squared)
        all_factors <- c(all_factors,paste0('Factor ',as.character(i)))
      }
    }
  }
  
  # plot the results
  tmp <- cbind.data.frame(rsq,all_factors)
  colnames(tmp) <- c('rsq','facts')
  tmp$facts <- as.factor(tmp$facts)
  df <- ddply(tmp, c("facts"), summarize, Mean = mean(rsq), SD = sd(rsq))
  p <- ggplot(df, aes(x = facts, y = Mean)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = Mean, ymax = Mean + SD), width = 0) +
    ylab("IFN-factor association (r-squared)") +
    xlab('') +
    ylim(0,1) +
    theme_bw()

  return(p)

}

p = plot_ifn_associations(pbmc_container,hybrid,ifn_core)
pdf(file = "/home/jmitchel/figures/for_paper_v2/hybrid_ifn.pdf", useDingbats = FALSE,
    width = 4.75, height = 3)
p
dev.off()

p = plot_ifn_associations(pbmc_container,ica_lds,ifn_core)
pdf(file = "/home/jmitchel/figures/for_paper_v2/ica_lds_ifn.pdf", useDingbats = FALSE,
    width = 4.75, height = 3)
p
dev.off()

p = plot_ifn_associations(pbmc_container,ica_dmat,ifn_core)
pdf(file = "/home/jmitchel/figures/for_paper_v2/ica_dmat_ifn.pdf", useDingbats = FALSE,
    width = 5.75, height = 3)
p
dev.off()


## show that batch factors have highly similar loadings patterns when do ICA on loadings
# do this by cor hmap

# looking at correlations of loadings for dmat rot
cormat_dmat_rot <- cor(t(ica_dmat[[2]]))
colnames(cormat_dmat_rot) <- sapply(1:ncol(cormat_dmat_rot),function(x){paste0('Factor ',x)})
rownames(cormat_dmat_rot) <- sapply(1:nrow(cormat_dmat_rot),function(x){paste0('Factor ',x)})

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
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

pdf(file = "/home/jmitchel/figures/for_paper_v2/ica_dmat_lds_cor.pdf", useDingbats = FALSE,
    width = 5, height = 4.5)
lds_hmap
dev.off()


cormat_dmat_rot <- cor(t(hybrid[[2]]))
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

pdf(file = "/home/jmitchel/figures/for_paper_v2/hybrid_lds_cor.pdf", useDingbats = FALSE,
    width = 5, height = 4.5)
lds_hmap
dev.off()



cormat_dmat_rot <- cor(t(ica_lds[[2]]))
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

pdf(file = "/home/jmitchel/figures/for_paper_v2/ica_lds_lds_cor.pdf", useDingbats = FALSE,
    width = 5, height = 4.5)
lds_hmap
dev.off()



# now showing better associations for sex-associated factor and ethnicity-associated factor
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Ethnicity','sex')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

# do the ica_dmat first
dscores <- ica_dmat[[1]]

tmp <- cbind.data.frame(dscores[,8],meta[rownames(dscores),'sex'])
colnames(tmp) <- c('dsc','sex')
mf_association <- ggplot(tmp,aes(x=sex,y=dsc)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.6, binwidth = .008) +
  ylab('Factor 8 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()

pdf(file = "/home/jmitchel/figures/for_paper_v2/ica_dmat_mf.pdf", useDingbats = FALSE,
    width = 4.5, height = 2.5)
mf_association
dev.off()


tmp <- cbind.data.frame(dscores[,3],meta[rownames(dscores),'Ethnicity'])
colnames(tmp) <- c('dsc','Ethnicity')
mf_association <- ggplot(tmp,aes(x=Ethnicity,y=dsc)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.6, binwidth = .008) +
  ylab('Factor 3 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()

pdf(file = "/home/jmitchel/figures/for_paper_v2/ica_dmat_eth.pdf", useDingbats = FALSE,
    width = 3.5, height = 3.25)
mf_association
dev.off()


# now for hybrid
dscores <- hybrid[[1]]

tmp <- cbind.data.frame(dscores[,7],meta[rownames(dscores),'sex'])
colnames(tmp) <- c('dsc','sex')
mf_association <- ggplot(tmp,aes(x=sex,y=dsc)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.6, binwidth = .008) +
  ylab('Factor 8 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()

pdf(file = "/home/jmitchel/figures/for_paper_v2/hybrid_mf.pdf", useDingbats = FALSE,
    width = 4.5, height = 2.5)
mf_association
dev.off()


tmp <- cbind.data.frame(dscores[,6],meta[rownames(dscores),'Ethnicity'])
colnames(tmp) <- c('dsc','Ethnicity')
mf_association <- ggplot(tmp,aes(x=Ethnicity,y=dsc)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.65, binwidth = .015) +
  ylab('Factor 3 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()

pdf(file = "/home/jmitchel/figures/for_paper_v2/hybrid_eth.pdf", useDingbats = FALSE,
    width = 3.5, height = 3.25)
mf_association
dev.off()













