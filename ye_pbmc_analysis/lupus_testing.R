



# generating all dsig genes plots
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=1,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['1']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=2,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['2']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=3,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['3']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=4,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['4']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=5,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['5']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=6,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['6']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=7,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=c(5,10,5,5,5,5,5),
                                       show_donor_labels=FALSE,
                                       additional_meta='Status',
                                       add_genes=c('IFI6','ISG15','MX1','XAF1'))
pbmc_container$plots$donor_sig_genes[['7']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=8,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=20,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['8']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=9,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=7,
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['9']]
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=10,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=c(5,5,5,50,5,5,5),
                                       show_donor_labels=FALSE,
                                       additional_meta='Status')
pbmc_container$plots$donor_sig_genes[['10']]

dsig_fig <- render_multi_plots(pbmc_container,data_type='dgenes',max_cols=5)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_dsig_genes.pdf", useDingbats = FALSE,
    width = 35, height = 35)
dsig_fig
dev.off()


# using a gene override
pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=5,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=c(5,5,5,5,5,5,5),
                                       show_donor_labels=FALSE,
                                       additional_meta='Status',
                                       add_genes=c('KNOP1','PPM1K'))
pbmc_container$plots$donor_sig_genes[['5']]


# generate dsig genes plot with donor scores heatmap appended so we can see that
# stripes correspond with donors from the other interferon pathways


pbmc_container <- plot_donor_sig_genes_v2(pbmc_container, factor_select=5,
                                          ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                          top_n_per_ctype=15, show_donor_labels=FALSE,
                                          additional_meta='Status')
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_dsig_f5_w_scores.pdf", useDingbats = FALSE,
    width = 20, height = 17)
draw(pbmc_container$plots$donor_sig_genes[['5']], column_title='Donors', column_title_side='bottom')
dev.off()

### for testing purposes only
# need to look at subtype proportions vs dscores
# Want to:
# double check associations
# see what the deal is with an overlapping top donors for contradictory factors
container <- pbmc_container
ctype <- 'cM'
res <- 0.5
resolution_name <- paste0('res:',as.character(res))
subclusts <- container$subclusters[[ctype]][[resolution_name]]

# append large cell type name to subclusters
subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})

# limit cells in subclusts to those that we actually have scores for
donor_scores <- container$tucker_results[[1]]
donor_vec <- container$scMinimal_full$metadata[names(subclusts),'donors']
subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]

# make subtype association plot
subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
scMinimal <- container$scMinimal_ctype[[ctype]]
sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

# get donor proportions of subclusters
donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)


subtype_associations <- get_indv_subtype_associations(pbmc_container,donor_props,5)




subtype <- donor_props[,5,drop=FALSE]

# add a second column with value of 1 - first column
subtype <- cbind(subtype,1-subtype)

# get balances
donor_balances <- coda.base::coordinates(subtype)



# for trying coordinates on counts
clusts <- subclusts_num
metadata <- sub_meta_tmp

# got donor_props count version from compute_donor_props
donor_balances <- coda.base::coordinates(donor_props)




# append dscores for factor 4
donor_props2 <- cbind(donor_props,donor_scores[rownames(donor_props),4])
colnames(donor_props2)[ncol(donor_props2)] <- 'dsc'

ggplot(as.data.frame(donor_props2),aes(x=as.numeric(dsc),y=as.numeric(K2))) +
  geom_point()

# append disease status
meta <- unique(container$scMinimal_full$metadata[,c('donors','Status')])
rownames(meta) <- meta$donors
donor_props2 <- cbind(donor_props2,as.character(meta[rownames(donor_props2),'Status']))
colnames(donor_props2)[ncol(donor_props2)] <- 'Status'

ggplot(as.data.frame(donor_props2),aes(x=as.numeric(dsc),y=as.numeric(K3),color=as.factor(Status))) +
  geom_point()



# see if any high loading donors from f1 are low loading in f3
top_f1 <- donor_scores[,1]
top_f1 <- top_f1[order(top_f1,decreasing=TRUE)][1:10]

top_f3 <- donor_scores[,3]
top_f3 <- top_f3[order(top_f3,decreasing=FALSE)][1:10]

sum(names(top_f3) %in% names(top_f1))

top_f3[names(top_f3) %in% names(top_f1)]







## trying to compute donor props using total cell numbers
clusts <- subclusts_num
metadata <- sub_meta_tmp

names(clusts) <- metadata[names(clusts),"donors"]
all_donors <- unique(as.character(metadata$donors))

# store results in df
donor_props <- data.frame(matrix(0,ncol=length(unique(clusts)),nrow = length(all_donors)))
colnames(donor_props) <- sapply(1:ncol(donor_props),function(x) {
  paste0('K',as.character(x))
})
rownames(donor_props) <- all_donors
for (d in all_donors) {
  tmp_clusts <- clusts[names(clusts)==d]
  counts <- table(tmp_clusts)
  names(counts) <- sapply(names(counts),function(x) {
    paste0('K',as.character(x))
  })
  for (j in 1:length(counts)) {
    donor_props[d,names(counts)[j]] <- counts[j]
  }
}
donor_props <- donor_props + 1 #adding pseudocount to avoid infinities when make balances
# donor_props <- t(apply(donor_props, 1, function(i) i/sum(i))) # counts -> props

new_totals <- table(container$scMinimal_full$metadata$donors)
donor_props <- t(sweep(t(donor_props),MARGIN=2,new_totals[rownames(donor_props)],FUN='/'))
subtype_associations <- get_indv_subtype_associations(container,donor_props,1)
subtype_associations





# I want to explore some of the differences a bit to see why some new things are significant and old are not...
container <- pbmc_container
ctype <- 'T8'
res <- 0.6
resolution_name <- paste0('res:',as.character(res))
subclusts <- container$subclusters[[ctype]][[resolution_name]]

# append large cell type name to subclusters
subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})

# limit cells in subclusts to those that we actually have scores for
donor_scores <- container$tucker_results[[1]]
donor_vec <- container$scMinimal_full$metadata[names(subclusts),'donors']
subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]

# make subtype association plot
subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
scMinimal <- container$scMinimal_ctype[[ctype]]
sub_meta_tmp <- scMinimal$metadata[names(subclusts),]

# get donor proportions of subclusters
donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)

# append dscores for factor 4
donor_props2 <- cbind(donor_props,donor_scores[rownames(donor_props),1])
colnames(donor_props2)[ncol(donor_props2)] <- 'dsc'

ggplot(as.data.frame(donor_props2),aes(x=as.numeric(dsc),y=as.numeric(K3))) +
  geom_point()

# append disease status
meta <- unique(container$scMinimal_full$metadata[,c('donors','Status')])
rownames(meta) <- meta$donors
donor_props2 <- cbind(donor_props2,as.character(meta[rownames(donor_props2),'Status']))
colnames(donor_props2)[ncol(donor_props2)] <- 'Status'

ggplot(as.data.frame(donor_props2),aes(x=as.numeric(dsc),y=as.numeric(K1),color=as.factor(Status))) +
  geom_point()


# trying different way to get balances
subtype_associations <- get_indv_subtype_associations(container,donor_props,5)

lmres <- lm(as.numeric(dsc)~as.numeric(K3),data=as.data.frame(donor_props2))
summary(lmres)

j <- 3
tmp <- donor_props[,j,drop=FALSE]
donor_props <- donor_props[,-j]
donor_props <- cbind(donor_props,tmp)
# donor_balances <- coda.base::coordinates(donor_props)
donor_balances <- compositions::ilr(donor_props)
rownames(donor_balances) <- rownames(donor_props)
donor_balances <- donor_balances[,ncol(donor_balances),drop=FALSE]
donor_props2 <- cbind(donor_balances,donor_scores[rownames(donor_balances),4])
colnames(donor_props2)[ncol(donor_props2)] <- 'dsc'
meta <- unique(container$scMinimal_full$metadata[,c('donors','Status')])
rownames(meta) <- meta$donors
donor_props2 <- cbind(donor_props2,as.character(meta[rownames(donor_props2),'Status']))
colnames(donor_props2)[ncol(donor_props2)] <- 'Status'
colnames(donor_props2)[1] <- 'ilr4'

ggplot(as.data.frame(donor_props2),aes(x=as.numeric(dsc),y=as.numeric(ilr4),color=as.factor(Status))) +
  geom_point()

head(donor_props2)

donor_props2 <- as.data.frame(donor_props2)
donor_props2$dsc <- as.numeric(donor_props2$dsc)
donor_props2$ilr4 <- as.numeric(donor_props2$ilr4)
rownames(donor_props2)[order(donor_props2[,'ilr4'],decreasing=F)][1:10]


container <- get_subclust_enr_hmap(container,all_ctypes,all_res,1:10)
container$plots$subc_enr_hmap


# getting donors with top scores for factor 7
tmp <- container$tucker_results[[1]][,7]
cool <- names(tmp)[tmp>.1]

test <- plot_dscore_enr(pbmc_container,factor_use=6,meta_var='Status')
test




# for testing enrichment of anergic pathway among factorss associated with naive cd4 depletion
my_pathways2 <- my_pathways[c('GSE46242_TH1_VS_ANERGIC_TH1_CD4_TCELL_DN','GSE46242_TH1_VS_ANERGIC_TH1_CD4_TCELL_UP','GSE5960_TH1_VS_ANERGIC_TH1_DN','GSE5960_TH1_VS_ANERGIC_TH1_UP')]

plt <- plotEnrichment(mypaths[[meta_vals[i]]],
                      myranks) + labs(title=paste0('',' - Factor ',as.character(6)))
plt <- plt +
  annotate(geom="text",  x=Inf, y=Inf, hjust=1,vjust=1, col="black",
           label=paste0('adj pval: ',
                        round(fgseaRes[fgseaRes$pathway==meta_vals[i],'padj'],digits=4)))


tmp <- fgsea_res[3,8][[1]][[1]]


pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=6,
                                       ctypes_use=c('cM','T8','NK','T4','ncM','cDC','B'),
                                       top_n_per_ctype=c(5,5,5,5,5,5,5),
                                       show_donor_labels=FALSE,
                                       additional_meta='Status',
                                       add_genes=tmp)

pbmc_container$plots$donor_sig_genes[['6']]

# a few caveats:
# -the relatively borderline pvalue
# -the gene set is from mouse
# -




# LR interaction analysis
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/NicheNet-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('from','to')]

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_f1.pdf", useDingbats = FALSE,
    width = 9, height = 9)
tmp <- get_LR_interact(container,lr_pairs,1)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_f2.pdf", useDingbats = FALSE,
    width = 9, height = 9)
tmp <- get_LR_interact(container,lr_pairs,2)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_f4.pdf", useDingbats = FALSE,
    width = 9, height = 9)
tmp <- get_LR_interact(container,lr_pairs,4)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_f5.pdf", useDingbats = FALSE,
    width = 9, height = 9)
tmp <- get_LR_interact(container,lr_pairs,5)
dev.off()











#### exploring new idea for LR interactions
container <- pbmc_container

factor_select <- 5

### prep stuff from lr fn
ctypes_use <- container$experiment_params$ctypes_use

# extract significance of all genes in all ctypes and put in list
sig_vectors <- get_significance_vectors(container,
                                        factor_select, ctypes_use)
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

# set 0 pvals to the min nonzero pval and take -log10
min_nonz <- min(sig_df[sig_df!=0])
sig_df[sig_df==0] <- min_nonz
sig_df <- -log10(sig_df)

# sign sig_df by loading
ldngs <- container$tucker_results[[2]]
genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})
sr_col <- ldngs[factor_select,]
tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)
tmp_casted_num <- tmp_casted_num[rownames(sig_df),colnames(sig_df)]
neg_mask <- tmp_casted_num < 0
sig_df[neg_mask] <- sig_df[neg_mask] * -1
###

lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')

# first identify ligands significantly associated with dscore for a factor
ligs <- lr_pairs[,'ligand']
ligs <- unique(ligs)

sig_ligs <- list()
for (lig in ligs) {
  if (lig %in% rownames(sig_df)) {
    for (ct in ctypes_use) {
      sig_val <- sig_df[lig,ct]
      if (abs(sig_val) > -log10(.01)) {
        sig_ligs[[paste0(lig,'_',ct)]] <- sig_val
      }
    } 
  }
}

print(sig_ligs)

# pick a ligand, see if receptor(s) present in any cell type
mylig <- 'RETN'
mylig <- 'LAIR1'
mylig <- 'THBS1'
mylig <- 'ADM'

myrec <- lr_pairs[lr_pairs$ligand==mylig,]
recs <- lapply(myrec$interaction_name,function(x) {
  tmp <- strsplit(x,split='_')[[1]]
  myrecs <- tmp[2:length(tmp)]
  return(myrecs)
})

print(recs)

# container <- get_pseudobulk(container)
# container <- normalize_pseudobulk(container, method='trim', scale_factor=10000)

lig_pos <- FALSE
r_thresh <- .01
for (j in 1:length(recs)) {
  rs <- recs[[j]]
  num_in_df <- sum(rs %in% rownames(sig_df))
  if (num_in_df == length(rs)) {
    for (ct in ctypes_use) {
      # need to use pseudobulked/normalized expression (but not scaled!)
      pb <- container$scMinimal_ctype[[ct]]$pseudobulk
      checks <- list()
      for (r in rs) {
        # determine direction of ligand expressing donors
        if (lig_pos) {
          dsc <- container$tucker_results[[1]]
          tmp <- dsc[,factor_select]
          tmp <- tmp[order(tmp,decreasing=TRUE)]
          top_n <- names(tmp)[1:10]
          d_exp <- sum(pb[top_n,r] > r_thresh)
          if (d_exp == 10) {
            checks[[r]] <- TRUE
          } else {
            checks[[r]] <- FALSE
          }
        } else {
          dsc <- container$tucker_results[[1]]
          tmp <- dsc[,factor_select]
          tmp <- tmp[order(tmp,decreasing=FALSE)]
          top_n <- names(tmp)[1:10]
          d_exp <- sum(pb[top_n,r] > r_thresh)
          if (d_exp == 10) {
            checks[[r]] <- TRUE
          } else {
            checks[[r]] <- FALSE
          }
        }
      }
      if (sum(unlist(checks))==length(checks)) {
        print(rs)
        print(ct)
      }
    }
  }
}


genes_test <- rownames(sig_df)[sig_df[,'ncM']<log10(.01)]
print(genes_test[1:30])

genes_test <- rownames(sig_df)[sig_df[,'cM']<log10(.01)]
print(genes_test[1:30])

genes_test <- rownames(sig_df)[sig_df[,'T4']<log10(.01)]
print(genes_test[1:30])

sig_df_tmp <- sig_df[,colnames(sig_df)!='cM']
sig_genes_other <- rownames(sig_df_tmp)[rowSums(sig_df_tmp<log10(.01)) > 0]
genes_test <- genes_test[!(genes_test %in% sig_genes_other)]

sig_df_tmp <- sig_df[,colnames(sig_df)!='T4']
sig_genes_other <- rownames(sig_df_tmp)[rowSums(sig_df_tmp<log10(.01)) > 0]
genes_test <- genes_test[!(genes_test %in% sig_genes_other)]

dsc <- container$tucker_results[[1]]
tmp <- dsc[,factor_select]

pb <- container$scMinimal_ctype[['ncM']]$pseudobulk
pb <- container$scMinimal_ctype[['cM']]$pseudobulk
pb <- container$scMinimal_ctype[['T4']]$pseudobulk

tmp2 <- cbind(tmp,pb[names(tmp),'HMGB2'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'IL1R2'])
plot(tmp2[,1],tmp2[,2])

# can look for genes with high spearman correlation to HMGB2
res <- list()
for (g in genes_test) {
  tmp2 <- cbind(pb[names(tmp),'UGCG'],pb[names(tmp),g])
  mycor <- cor(tmp2,method='spearman')[1,2]
  res[[g]] <- mycor
}
res <- unlist(res)
res <- res[order(res,decreasing=TRUE)]
print(res[1:10])

tmp2 <- cbind(tmp,pb[names(tmp),'FKBP5'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'ZFAND5'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'ETS2'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'CXCR4'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'ISG15'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'IFI6'])
plot(tmp2[,1],tmp2[,2])

# trying to normalize by receptor levels
rlevs <- pb[names(tmp),'CD36']
sum(rlevs==0)
tmp2 <- cbind(tmp,pb[names(tmp),'SESN1']/rlevs)
plot(tmp2[,1],tmp2[,2])



# negative controls
tmp2 <- cbind(tmp,pb[names(tmp),'ISG15'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'IFI6'])
plot(tmp2[,1],tmp2[,2])



# looking at ligand expression
pb <- container$scMinimal_ctype[['cDC']]$pseudobulk

tmp2 <- cbind(tmp,pb2[names(tmp),'RETN'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'THBS1'])
plot(tmp2[,1],tmp2[,2])

pb <- container$scMinimal_ctype[['cM']]$pseudobulk
tmp2 <- cbind(tmp,pb2[names(tmp),'ADM'])
plot(tmp2[,1],tmp2[,2])

pb <- container$scMinimal_ctype[['cM']]$pseudobulk
tmp2 <- cbind(tmp,pb[names(tmp),'ZBTB16'])
plot(tmp2[,1],tmp2[,2])


# look at combined ligand levels for all expressing ctypes
pb <- container$scMinimal_ctype[['cM']]$pseudobulk

pb1 <- container$scMinimal_ctype[['cDC']]$pseudobulk
pb2 <- container$scMinimal_ctype[['cM']]$pseudobulk
pb3 <- container$scMinimal_ctype[['ncM']]$pseudobulk

tmp2 <- cbind(tmp,pb1[names(tmp),'ADM']+pb2[names(tmp),'ADM']+pb3[names(tmp),'ADM'])
plot(tmp2[,1],tmp2[,2])


res <- list()
for (i in 1:ncol(pb)) {
  tmp3 <- cbind(pb[names(tmp),i],tmp2[,2])
  mycor <- cor(tmp3,method='pearson')[1,2]
  res[[colnames(pb)[i]]] <- mycor
}
res <- unlist(res)
res <- res[order(res,decreasing=TRUE)]
print(res[1:10])




pb1 <- container$scMinimal_ctype[['T8']]$pseudobulk
pb2 <- container$scMinimal_ctype[['T4']]$pseudobulk

tmp2 <- cbind(tmp,pb1[names(tmp),'CD70']+pb2[names(tmp),'CD70'])
plot(tmp2[,1],tmp2[,2])


# look at receptor levels
pb <- container$scMinimal_ctype[['T4']]$pseudobulk

tmp2 <- cbind(tmp,pb[names(tmp),'CD36'])
plot(tmp2[,1],tmp2[,2])



# incorporate ratio to receptor levels
# use total ligand expression levels from the different cell types (maybe)
# show plots for "negative controls" things like ISG15, which are associated with 
# score but shouldn't show the expected trend


# looking to see if CALCRL explains ZBTB16 levels

pb <- container$scMinimal_ctype[['cM']]$pseudobulk
tmp2 <- cbind(pb[names(tmp)[1:12],'CALCRL'],pb[names(tmp)[1:12],'ZBTB16'])
plot(tmp2[,1],tmp2[,2])

# doesnt really appear to explain it... 
# now I'll do unbiased search of receptors that might explain it
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/NicheNet-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('from','to')]
myres <- list()
rs_test <- unique(lr_pairs$to)
for (r in rs_test) {
  if (r %in% colnames(pb)) {
    tmp2 <- as.data.frame(cbind(pb[names(tmp)[1:12],r],pb[names(tmp)[1:12],'ZBTB16']))
    colnames(tmp2) <- c('rec','sig')
    lmres <- lm(sig~rec,data=tmp2)
    lmres <- summary(lmres)
    # myres[[r]] <-  lmres$fstatistic[[1]]
    myres[[r]] <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
  }
}

myres <- unlist(myres)
myres <- p.adjust(myres,method='fdr')
myres[order(myres,decreasing=F)][1:10]
# myres[order(myres,decreasing=TRUE)][1:10]

tmp2 <- cbind(pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'ZBTB16'])
plot(tmp2[,1],tmp2[,2])


# I also showed increased levels of THBS1 in these patients and this happens to be the
# ligand for the CD47 gene!
# I wonder if I can use this to estimate Kd values for signal transduction since I have
# relative ligand and receptor levels and saturation


# to gain further evidence, I should see if CD47 correlates with residuals of other genes with the same patterns
# will try FKBP5

tmp2 <- cbind(pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'FKBP5'])
plot(tmp2[,1],tmp2[,2])

# it looks pretty okay actually
# whats the pvalue though?

tmp2 <- as.data.frame(cbind(pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'FKBP5']))
colnames(tmp2) <- c('rec','sig')
lmres <- lm(sig~rec,data=tmp2)
lmres <- summary(lmres)


# trying also for CCND3
tmp2 <- cbind(pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'CCND3'])
plot(tmp2[,1],tmp2[,2])


# seeing how well THBS1 levels correlate without dscore
tmp2 <- cbind(tmp,pb[names(tmp),'THBS1'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb1[names(tmp),'CLSTN1'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'CD47'])
plot(tmp2[,1],tmp2[,2])

pb1 <- container$scMinimal_ctype[['cDC']]$pseudobulk
pb2 <- container$scMinimal_ctype[['cM']]$pseudobulk
pb3 <- container$scMinimal_ctype[['ncM']]$pseudobulk

pb4 <- container$scMinimal_ctype[['NK']]$pseudobulk
pb5 <- container$scMinimal_ctype[['B']]$pseudobulk
pb6 <- container$scMinimal_ctype[['T8']]$pseudobulk
pb7 <- container$scMinimal_ctype[['T4']]$pseudobulk


tmp2 <- cbind(tmp,pb1[names(tmp),'THBS1']+pb2[names(tmp),'THBS1']+pb3[names(tmp),'THBS1'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb3[names(tmp),'THBS1'])
plot(tmp2[,1],tmp2[,2])

sig_df['THBS1',]


# looking for additional genes that are similar to FKBP5 and ZBTB16
pb <- container$scMinimal_ctype[['cM']]$pseudobulk
res <- list()
for (i in 1:ncol(pb)) {
  tmp3 <- cbind(pb[names(tmp),i],pb[names(tmp),'ZBTB16'])
  mycor <- cor(tmp3,method='spearman')[1,2]
  res[[colnames(pb)[i]]] <- mycor
}
res <- unlist(res)
res <- res[order(res,decreasing=TRUE)]
print(res[1:10])


'ZBTB16' %in% rownames(sig_df)
'FKBP5' %in% rownames(sig_df)

tmp2 <- cbind(tmp,pb[names(tmp),'TSPAN14'])
plot(tmp2[,1],tmp2[,2])

# best candidates:
# ETS2, SMAP2
# second best:
# TSPAN14, IRS2, SLA, IRAK3, HMGB2

# now seeing if cd47 predicts expression at upper end for these hits
tmp2 <- as.data.frame(cbind(pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'SMAP2']))
colnames(tmp2) <- c('rec','sig')
lmres <- lm(sig~rec,data=tmp2)
lmres <- summary(lmres)
print(lmres)

# ones with variance explained by cd47
# SMAP2 and maybe TSPAN14

tmp2 <- cbind(pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'ZBTB16'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp)[1:15],'CD47'],pb[names(tmp)[1:15],'ZBTB16'])
plot(tmp2[,1],tmp2[,2])


# plot target expression vs ligand expression for top donors only
tmp2 <- cbind(pb[names(tmp)[1:12],'THBS1'],pb[names(tmp)[1:12],'ZBTB16'])
plot(tmp2[,1],tmp2[,2])

# trying it with all datapoints
tmp2 <- cbind(pb[names(tmp),'THBS1'],pb[names(tmp),'ZBTB16'])
plot(tmp2[,1],tmp2[,2])

# seeing if the lm is significant
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('THBS1','ZBTB16')
lmres <- lm(ZBTB16 ~ THBS1, data=tmp2)
summary(lmres)

# Looking for receptors linked to some of the other level-off genes
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
rs_test <- sapply(lr_pairs$interaction_name, function(x) {
  myspl <- strsplit(x,split='_')[[1]]
  return(myspl[[2]])
})

lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/NicheNet-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('from','to')]
rs_test <- unique(lr_pairs$to)
myres <- list()
for (r in rs_test) {
  if (r %in% colnames(pb)) {
    tmp2 <- as.data.frame(cbind(pb[names(tmp)[1:15],r],pb[names(tmp)[1:15],'SPTLC2']))
    colnames(tmp2) <- c('rec','sig')
    lmres <- lm(sig~rec,data=tmp2)
    lmres <- summary(lmres)
    # myres[[r]] <-  lmres$fstatistic[[1]]
    myres[[r]] <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
  }
}

myres <- unlist(myres)
myres <- p.adjust(myres,method='fdr')
myres[order(myres,decreasing=F)][1:10]


# for CD44, plotting ligand expression vs dscore
tmp2 <- cbind(tmp,pb[names(tmp),'LGALS9'])
tmp2 <- cbind(tmp,pb6[names(tmp),'LGALS9']+pb2[names(tmp),'LGALS9']+pb3[names(tmp),'LGALS9'])
plot(tmp2[,1],tmp2[,2])

# need to show ligand levels not associated with levels of SPTLC2 top genes
tmp2 <- cbind(pb[names(tmp),'LGALS9'],pb[names(tmp),'SPTLC2'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp),'CD44'],pb[names(tmp),'LGALS9'],pb[names(tmp),'SPTLC2'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('CD44','LGALS9','SPTLC2')
lmres <- lm(SPTLC2 ~ LGALS9 + CD44,data=tmp2)
summary(lmres)

fit1 <- lm(as.formula('SPTLC2 ~ LGALS9'), data=tmp2)
fit2 <- lm(as.formula('SPTLC2 ~ CD44 + LGALS9'), data=tmp2)
anova(fit1, fit2)

sum(is.na(tmp2))

# how well does the receptor expression correlate with dscore
tmp2 <- cbind(tmp,pb[names(tmp),'CD44'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp)[1:15],'CD44'],pb[names(tmp)[1:15],'SPTLC2'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp)[1:15],'LGALS9'],pb[names(tmp)[1:15],'SPTLC2'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp)[1:15],'THBS1'],pb[names(tmp)[1:15],'FKBP5'])
plot(tmp2[,1],tmp2[,2])


# trying to fit an exponential model
tmp2 <- cbind(tmp,pb[names(tmp),'FKBP5'])
plot(tmp2[,1],tmp2[,2])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('dsc','FKBP5')
mod <- nls(FKBP5 ~ -a**(b*dsc) + d, data = tmp2, start = list(a = 0, b = 0, d = 0))
mod <- nls(FKBP5 ~ a - dsc**2, data = tmp2, start = list(a = 0, b = 0),
           lower=c(-10,.1), upper=c(10,5))
lines(tmp2$dsc, predict(mod, list(x = tmp2$dsc)))

fm1 <- nls(FKBP5 ~ SSasymp( dsc, Asym, resp0, lrc), data = tmp2)

# fit3 <- lm(FKBP5~poly(dsc,3,raw=TRUE),data=tmp2)
# summary(fit3)
# xx <- seq(-.3,.1, length=50)
# cool <- predict(fit3, as.data.frame(as.matrix(xx)))
# lines(xx, , col="blue")

model <- drm(FKBP5 ~ dsc, fct = DRC.asymReg())


nlsfit <- nls(FKBP5 ~ SSasymp( dsc, Asym, resp0, lrc), data = tmp2)
p <- coef(nlsfit)

# Plot and add a curve. Note that 'p' is a named vector.
with(Loblolly, plot(age, height, xlim=c(0, 30), ylim=c(0,75)))
curve(fx(dsc, Asym=p["Asym"], R0=p["R0"], lrc=p["lrc"]), 
      add=T, col="red")

# Predict from the fitted model for the new dataframe.
newdat <- data.frame(age=dsc)
newdat$heightpred <- predict(nlsfit, newdata=newdat)
with(newdat, lines(tmp2$dsc, heightpred, col="blue"))

lines(tmp2$dsc, newdat$heightpred, col="blue")

library(easynls)
# trying a gompertzmodel
tmp2[,3] <- tmp2[,1]
tmp2[,1] <- tmp2[,2]
tmp2[,2] <- tmp2[,3]
tmp2[,3] <- NULL
colnames(tmp2) <- c('FKBP5','dsc')
tmp2$dsc <- tmp2$dsc*-1
tmp2$dsc <- tmp2$dsc + .2
head(tmp2)
model2 <- nlsfit(tmp2, model = 10, start = c(a = 1, b = .4, d = .5))

plot(tmp2)

nlsplot(tmp2, model = 10, start = c(a = 1, b = .4, c = .1), 
        xlab = "Days" , ylab = "Tumor Volume", position = 1)




# trying to plot ratio of ligand to receptor vs target levels
tmp2 <- cbind(pb[names(tmp)[1:12],'THBS1']/pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'ZBTB16'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp)[1:12],'THBS1']/pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'TSPAN14'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'FKBP5'])
plot(tmp2[,1],tmp2[,2])


# trying to cluster the jackstraw significant genes by umap
library(umap)
library(ggplot2)
sig_genes <- sig_df[,'cM']
sig_genes <- sig_genes[sig_genes<log10(.01)]
sig_exp <- pb[,names(sig_genes)]
sig_exp <- scale(sig_exp) # trying with scale
# ## trying with pca first
# sig_exp <- prcomp(t(sig_exp))
# sig_exp <- sig_exp$x[,1:10]
# sig_umap <- umap(as.matrix(sig_exp))
# ##
sig_umap <- umap(as.matrix(t(sig_exp)))
sig_umap <- sig_umap[["layout"]]
sig_umap <- as.data.frame(sig_umap)
colnames(sig_umap) <- c('UMAP1','UMAP2')
gene_search <- c('FKBP5','ZBTB16','TSPAN14','SMAP2')
gene_search <- c('FKBP5','ZBTB16')
gene_search <- c('SPTLC2')
sig_umap$test <- sapply(rownames(sig_umap),function(x) {
  if (x %in% gene_search) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})
ggplot(sig_umap,aes(x=UMAP1,y=UMAP2,color=test)) +
  geom_point()

# get knn genes to ZBTB16
library(FNN)
test <- get.knn(sig_umap, k=10)
zbt_ndx <- which(rownames(sig_umap)=='ZBTB16')
zbt_ndx <- which(rownames(sig_umap)=='CCND3')
rownames(sig_umap)[test[["nn.index"]][zbt_ndx,]]

# here is ZBTB16 again for reference
tmp2 <- cbind(tmp,pb[names(tmp),'ZBTB16'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,exp(pb[names(tmp),'SPTLC2']))
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'SPTLC2'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'CLEC4E'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp)[1:12],'CD47'],pb[names(tmp)[1:12],'CLEC4E'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'IRAK3'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp)[1:20],'CD47'],pb[names(tmp)[1:20],'IRAK3'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp)[1:20],'PLXNC1'],pb[names(tmp)[1:20],'SPRY1'])
plot(tmp2[,1],tmp2[,2])

lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/NicheNet-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('from','to')]
rs_test <- unique(lr_pairs$to)
myres <- list()
for (r in rs_test) {
  if (r %in% colnames(pb)) {
    tmp2 <- as.data.frame(cbind(pb[names(tmp)[1:20],r],pb[names(tmp)[1:20],'ZFP36L2']))
    colnames(tmp2) <- c('rec','sig')
    lmres <- lm(sig~rec,data=tmp2)
    lmres <- summary(lmres)
    # myres[[r]] <-  lmres$fstatistic[[1]]
    myres[[r]] <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
  }
}

myres <- unlist(myres)
myres <- p.adjust(myres,method='fdr')
myres[order(myres,decreasing=F)][1:10]



# can explicitly test target genes for links with CD47
pb <- container$scMinimal_ctype[['cM']]$pseudobulk
res <- list()
for (j in 1:nrow(sig_df)) {
  if (sig_df[j,'cM'] < log10(.001)) {
    i <- rownames(sig_df)[j]
    tmp3 <- cbind(pb[names(tmp)[1:15],'CD47'],pb[names(tmp)[1:15],i])
    mycor <- cor(tmp3,method='pearson')[1,2]
    res[[i]] <- mycor
  }
}
res <- unlist(res)
res <- res[order(res,decreasing=FALSE)]
print(res[1:10])

tmp2 <- cbind(tmp,pb[names(tmp),'ZFP36L2'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp)[1:15],'CD47'],pb[names(tmp)[1:15],'ZFP36L2'])
plot(tmp2[,1],tmp2[,2])

# ones that look the same
# 'CNIH4', 'PHF11'
# 'CHMP5', 'GCA', 'FPR2', 

tmp2 <- cbind(pb[names(tmp),'CD47'],pb[names(tmp),'ITSN1'])
plot(tmp2[,1],tmp2[,2])



# seeing if the CD47-FKBP5/ZBTB16 link can be found in other cell types too
pb1 <- container$scMinimal_ctype[['cDC']]$pseudobulk
pb2 <- container$scMinimal_ctype[['cM']]$pseudobulk
pb3 <- container$scMinimal_ctype[['ncM']]$pseudobulk

tmp2 <- cbind(tmp,pb1[names(tmp),'SMAP2'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb1[names(tmp)[1:10],'CD47'],pb1[names(tmp)[1:10],'SMAP2'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb1[names(tmp),'CD47'],pb1[names(tmp),'SH3TC1'])
plot(tmp2[,1],tmp2[,2])


# going to try making correlation hmap as umap wasn't very helpful...
library(ComplexHeatmap)
sig_genes <- sig_df[,'cM']
sig_genes <- sig_genes[sig_genes<log10(.001)]
sig_exp <- pb[,names(sig_genes)]
cool <- cor(as.matrix(sig_exp))
Heatmap(cool,name='gene cor',row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5))

tmp2 <- cbind(tmp,pb[names(tmp),'MYL12A'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb[names(tmp),'VPS29'])
plot(tmp2[,1],tmp2[,2])


lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/NicheNet-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('from','to')]
rs_test <- unique(lr_pairs$to)
myres <- list()
for (r in rs_test) {
  if (r %in% colnames(pb)) {
    tmp2 <- as.data.frame(cbind(pb[names(tmp)[1:20],r],pb[names(tmp)[1:20],'VPS29']))
    colnames(tmp2) <- c('rec','sig')
    lmres <- lm(sig~rec,data=tmp2)
    lmres <- summary(lmres)
    # myres[[r]] <-  lmres$fstatistic[[1]]
    myres[[r]] <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
  }
}

myres <- unlist(myres)
myres <- p.adjust(myres,method='fdr')
myres[order(myres,decreasing=F)][1:14]


tmp2 <- cbind(tmp,pb[names(tmp),'TNFRSF1B'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb3[names(tmp),'TNF'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp)[1:20],'TNFRSF1B'],pb[names(tmp)[1:20],'VPS29'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb[names(tmp),'TNFRSF1B'],pb[names(tmp),'VPS29'])
plot(tmp2[,1],tmp2[,2])


## looking at hmap for a different cell type
# going to try making correlation hmap as umap wasn't very helpful...
library(ComplexHeatmap)
pb3 <- container$scMinimal_ctype[['T8']]$pseudobulk
sig_genes <- sig_df[,'T8']
sig_genes <- sig_genes[sig_genes<log10(.01)]
sig_exp <- pb3[,names(sig_genes)]
cool <- cor(as.matrix(sig_exp),method='pearson')
Heatmap(cool,name='gene cor',row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7))

tmp2 <- cbind(tmp,scale(pb3[names(tmp),'IFITM1'])*(3.797986**1.5))
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb3[names(tmp),'TXNIP'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb3[names(tmp),'CCR5'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb3[names(tmp),'CCL4'],pb3[names(tmp),'KLRD1'],pb3[names(tmp),'CCR5'])
tmp2 <- cbind(pb3[names(tmp),'KLRD1'],pb3[names(tmp),'CCL4'],pb3[names(tmp),'CCR5'])
plot(tmp2[,3],tmp2[,2])

tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('KLRD1','CCL4','CCR5')
mod1 <- lm(KLRD1 ~ CCR5 + CCL4, data=tmp2)
mod2 <- lm(KLRD1 ~ CCR5, data=tmp2)
anova(mod2, mod1)
summary(mod2)

tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('dsc','expr')
lmres <- lm(dsc~expr,data=tmp2)
summary(lmres)

tmp2 <- cbind(tmp,pb3[names(tmp),'CCL4'])
plot(tmp2[,1],tmp2[,2])


tmp2 <- cbind(tmp,pb3[names(tmp),'FKBP5'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb3[names(tmp)[1:30],'CCR5'],pb3[names(tmp)[1:30],'KLRD1'])
plot(tmp2[,1],tmp2[,2])

# looking at total ccl4 levels across all cell types
pb1 <- container$scMinimal_ctype[['cDC']]$pseudobulk
pb2 <- container$scMinimal_ctype[['cM']]$pseudobulk
pb3 <- container$scMinimal_ctype[['T8']]$pseudobulk
pb4 <- container$scMinimal_ctype[['T4']]$pseudobulk
pb5 <- container$scMinimal_ctype[['ncM']]$pseudobulk
pb6 <- container$scMinimal_ctype[['B']]$pseudobulk
pb7 <- container$scMinimal_ctype[['NK']]$pseudobulk

tmp2 <- cbind(tmp,pb3[names(tmp),'CCL4']+pb6[names(tmp),'CCL4']+pb7[names(tmp),'CCL4'])
plot(tmp2[,1],tmp2[,2])

lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/NicheNet-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('from','to')]
rs_test <- unique(lr_pairs$to)
myres <- list()
for (r in rs_test) {
  if (r %in% colnames(pb3)) {
    tmp2 <- as.data.frame(cbind(pb3[names(tmp)[1:12],r],pb3[names(tmp)[1:12],'SRI']))
    colnames(tmp2) <- c('rec','sig')
    lmres <- lm(sig~rec,data=tmp2)
    lmres <- summary(lmres)
    # myres[[r]] <-  lmres$fstatistic[[1]]
    myres[[r]] <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
  }
}

myres <- unlist(myres)
myres <- p.adjust(myres,method='fdr')
myres[order(myres,decreasing=F)][1:14]

myres[order(myres,decreasing=T)][1:14]


tmp2 <- cbind(pb3[names(tmp)[1:10],'FAS'],pb3[names(tmp)[1:10],'PRDM1'])
plot(tmp2[,1],tmp2[,2])

tmp2 <- cbind(pb3[names(tmp),'FAS'],pb3[names(tmp),'PRDM1'])
plot(tmp2[,1],tmp2[,2])


tmp2 <- cbind(pb3[names(tmp),'TNF'],pb3[names(tmp),'PRDM1'])
plot(tmp2[,1],tmp2[,2])



# create plot with coloring by batch
meta <- container$scMinimal_full$metadata[,c('donors','pool')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

tmp2 <- cbind(tmp,pb3[names(tmp),'IFITM1'],meta[names(tmp),1,drop=F])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('dsc','eval','pool')
ggplot(tmp2,aes(x=dsc,y=eval,color=pool)) +
  geom_point()



### tests for my formal interaction model

# first want to show that for a negative control like ISG15, it shouldn't be associated with the ligand
# THBS1 and receptor CD47
tmp2 <- cbind(pb[names(tmp),'CD47'],pb[names(tmp),'THBS1'],pb[names(tmp),'ISG15'])
plot(tmp2[,1],tmp2[,2])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('CD47','THBS1','ISG15')
lmres <- lm(ISG15~CD47 * THBS1,data=tmp2)
summary(lmres)

# check again that THBS1 levels predict ZBTB16 levels
tmp2 <- cbind(pb[names(tmp),'CD47'],pb[names(tmp),'THBS1'],pb[names(tmp),'ZBTB16'])
plot(tmp2[,2],tmp2[,3])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('CD47','THBS1','ZBTB16')
# lmres <- lm(ZBTB16~CD47 * THBS1,data=tmp2)
lmres1 <- lm(ZBTB16~THBS1,data=tmp2)
lmres2 <- lm(ZBTB16~THBS1 + CD47,data=tmp2)
lmres3 <- lm(ZBTB16~THBS1 * CD47,data=tmp2)
anova(lmres2,lmres3)
summary(lmres)

tmp2 <- cbind(pb2[names(tmp)[1:15],'CD47'],pb2[names(tmp)[1:15],'FKBP5'])
plot(tmp2[,1],tmp2[,2])


# see if interaction improves model in T8 for CCL4 ligand, CCR5 receptor, and KLRD1 target
tmp2 <- cbind(pb3[names(tmp),'CCR5'],pb3[names(tmp),'CCL4'],pb3[names(tmp),'KLRD1'])
plot(tmp2[,1],tmp2[,3])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('CCR5','CCL4','KLRD1')
lmres1 <- lm(KLRD1~CCL4,data=tmp2)
lmres2 <- lm(KLRD1~CCL4 + CCR5,data=tmp2)
lmres3 <- lm(KLRD1~CCL4 * CCR5,data=tmp2)
lmres4 <- lm(KLRD1~CCR5,data=tmp2)
anova(lmres1,lmres2)
anova(lmres1,lmres3)
anova(lmres2,lmres3)
summary(lmres4)


# look in NK cells where CCR1 (receptor for CCL4) was significantly associated with factor 5
# first need to get a gene or module to test against
sig_genes <- sig_df[,'NK']
sig_genes <- sig_genes[sig_genes<log10(.05)]
sig_exp <- pb7[,names(sig_genes)]
cool <- cor(as.matrix(sig_exp),method='pearson')
Heatmap(cool,name='gene cor',row_names_gp = gpar(fontsize = 7),
        column_names_gp = gpar(fontsize = 7))
Heatmap(cool,name='gene cor',show_row_names=F,show_column_names=F)

tmp2 <- cbind(pb7[names(tmp),'CCR1'],pb3[names(tmp),'CCL4'],pb7[names(tmp),'ZBP1'])
plot(tmp2[,2],tmp2[,3])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('CCR1','CCL4','ZBP1')
lmres1 <- lm(ZBP1~CCL4,data=tmp2)
lmres2 <- lm(ZBP1~CCL4 + CCR1,data=tmp2)
lmres3 <- lm(ZBP1~CCL4 * CCR1,data=tmp2)
lmres4 <- lm(ZBP1~CCR1,data=tmp2)
anova(lmres1,lmres2)
anova(lmres1,lmres3)
anova(lmres2,lmres3)
summary(lmres3)

# interesting that interaction is significant but neither ligand nor receptor is alone...
tmp2 <- cbind(pb7[names(tmp),'CCR1'],pb7[names(tmp),'CCL4'],pb7[names(tmp),'TSC22D3'])
plot(tmp2[,2],tmp2[,3])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('CCR1','CCL4','TSC22D3')
lmres1 <- lm(TSC22D3~CCL4,data=tmp2)
lmres2 <- lm(TSC22D3~CCL4 + CCR1,data=tmp2)
lmres3 <- lm(TSC22D3~CCL4 * CCR1,data=tmp2)
lmres4 <- lm(TSC22D3~CCR1,data=tmp2)
anova(lmres1,lmres2)
anova(lmres1,lmres3)
anova(lmres2,lmres3)
summary(lmres1)

# this is a case of a bit stronger evidence where receptor and ligand both associated and each add to prediction of target levels

# now I'll try some other genes in the same module as TSC22D3
tmp2 <- cbind(pb7[names(tmp),'CCR1'],pb3[names(tmp),'CCL4'],pb7[names(tmp),'PIK3IP1'])
plot(tmp2[,2],tmp2[,3])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('CCR1','CCL4','PIK3IP1')
lmres1 <- lm(PIK3IP1~CCL4,data=tmp2)
lmres2 <- lm(PIK3IP1~CCL4 + CCR1,data=tmp2)
lmres3 <- lm(PIK3IP1~CCL4 * CCR1,data=tmp2)
lmres4 <- lm(PIK3IP1~CCR1,data=tmp2)
anova(lmres1,lmres2)
anova(lmres1,lmres3)
anova(lmres2,lmres3)
summary(lmres1)


tmp2 <- cbind(tmp,pb6[names(tmp),'TGFB1'])
plot(tmp2[,1],tmp2[,2])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('dsc','TGFB1')
lmres1 <- lm(TGFB1~dsc,data=tmp2)
summary(lmres1)

## need to ensure all ligands and receptors are in final dataset

# to get a measure of specificity I'm wondering I can just use cor coefficient or R-sq
# I'll do a comparison of SAP30/ADM which I consider speicific and ISG15/IFI6 which I don't
tmp2 <- cbind(tmp,pb2[names(tmp),'SAP30'])
plot(tmp2[,1],tmp2[,2])
cor(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb2[names(tmp),'ADM'])
plot(tmp2[,1],tmp2[,2])
cor(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb2[names(tmp),'IFI6'])
plot(tmp2[,1],tmp2[,2])
cor(tmp2[,1],tmp2[,2])

tmp2 <- cbind(tmp,pb2[names(tmp),'ISG15'])
plot(tmp2[,1],tmp2[,2])
cor(tmp2[,1],tmp2[,2])

# actually correlation coefficient seems to give a pretty solid sense of specificity to the factor...
# will look at rsq too
tmp2 <- cbind(tmp,pb2[names(tmp),'SAP30'])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('dsc','SAP30')
lmres1 <- lm(SAP30~dsc,data=tmp2)
summary(lmres1)

tmp2 <- cbind(tmp,pb2[names(tmp),'ADM'])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('dsc','ADM')
lmres1 <- lm(ADM~dsc,data=tmp2)
summary(lmres1)

tmp2 <- cbind(tmp,pb2[names(tmp),'IFI6'])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('dsc','IFI6')
lmres1 <- lm(IFI6~dsc,data=tmp2)
summary(lmres1)

tmp2 <- cbind(tmp,pb2[names(tmp),'ISG15'])
tmp2 <- as.data.frame(tmp2)
colnames(tmp2) <- c('dsc','ISG15')
lmres1 <- lm(ISG15~dsc,data=tmp2)
summary(lmres1)

# r-squared seems to do well too, and I like it better for it's interpretability...







# prep for new LR analysis
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
  rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
  return(rname)
})
lr_pairs$interaction_name <- NULL




# trying to cluster genes by gene factors
library(umap)
gene_by_factors
# gene_umap <- umap(abs(gene_by_factors))
gene_umap <- umap(gene_by_factors)
gene_umap <- gene_umap[["layout"]]
head(gene_umap)
gene_umap <- as.data.frame(gene_umap)
colnames(gene_umap) <- c('UMAP1','UMAP2')
g_color <- c('IFI6','ISG15','MX1','MX2','XAF1')
g_color <- c('HLA-A','HLA-B','HLA-C','HLA-D')
g_color <- c('C1QA','C1QB')
gene_umap$mycol <- sapply(rownames(gene_umap),function(x) {
  if (x %in% g_color) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})
ggplot(gene_umap,aes(x=UMAP1,y=UMAP2,color=gene_by_factors[,10])) +
  geom_point() +
  scale_color_gradient2(midpoint=0, low="blue", mid="white",
                        high="red" )
ggplot(gene_umap,aes(x=UMAP1,y=UMAP2,color=mycol)) +
  geom_point() 
ggplot(gene_umap,aes(x=UMAP1,y=UMAP2,color=gene_by_factors[,20])) +
  geom_point() +
  scale_color_gradient2(midpoint=0, low="blue", mid="white",
                        high="red" )
# I wonder if the high and low scoring genes for a factor should actually be considered together
# a gene factor represents a pattern across donors and cell types where some genes are highly expressed
# and some genes are lowely expressed.

# trying dynamic tree cut clustering on gene factors matrix
library(dynamicTreeCut)
library(cluster)
gf_dist <- dist(gene_by_factors)
g_dend <- hclust(gf_dist,method='average')
ct <- cutreeDynamic(g_dend,cutHeight=.90)
ct <- kmeans(gene_by_factors,centers=15)
ct <- pam(gene_umap,k=30)
ct <- pam(gene_by_factors,k=25)
ggplot(gene_umap,aes(x=UMAP1,y=UMAP2,color=as.factor(ct$cluster))) +
  geom_point() + 
  scale_fill_brewer(palette="Dark2")

ct$cluster['ISG15']
ct$cluster[ct$cluster==20]

# pb7 <- container$scMinimal_ctype[['NK']]$pseudobulk
pb_sub <- pb2[,names(ct$cluster[ct$cluster==20])]
# pb_sub <- pb4[,'XAF1']
# eigengene <- pb_sub

# dim(pb_sub)
eigengene <- prcomp(pb_sub)$x[,1]
dsc <- container$tucker_results[[1]][,5]
cor(dsc[names(eigengene)],eigengene)
plot(dsc[names(eigengene)],eigengene)

cool <- as.data.frame(cbind(dsc[names(eigengene)],eigengene))
colnames(cool) <- c('dsc','gn')
lmres <- lm(dsc~gn,data=cool)
summary(lmres)


ct$cluster[ct$cluster==15]

pb_sub <- pb2[,names(ct$cluster[ct$cluster==12])]
# pb_sub <- pb4[,'XAF1']
# eigengene <- pb_sub

# dim(pb_sub)
eigengene <- prcomp(pb_sub)$x[,1]
dsc <- container$tucker_results[[1]][,1]
cor(dsc[names(eigengene)],eigengene)
plot(dsc[names(eigengene)],eigengene)

cool <- as.data.frame(cbind(dsc[names(eigengene)],eigengene))
colnames(cool) <- c('dsc','gn')
lmres <- lm(dsc~gn,data=cool)
summary(lmres)

### making heatmap for correlation of gene clusters in a cell type with a factor
# remove clusters under certain size
myclusts <- ct$cluster
c_counts <- table(myclusts)
clust_keep <- names(c_counts)[c_counts>5]
myclusts <- myclusts[myclusts %in% clust_keep]
myclusts <- sapply(myclusts,function(x) {
  paste0('gclust_',x)
})

# set up results df
ctypes_use <- container$experiment_params$ctypes_use
myres<- data.frame(matrix(ncol=length(ctypes_use),nrow=length(unique(myclusts))))
colnames(myres) <- ctypes_use
rownames(myres) <- unique(myclusts)
factor_select <- 5
dsc <- container$tucker_results[[1]][,factor_select]

for (ct in ctypes_use) {
  pb <- container$scMinimal_ctype[[ct]]$pseudobulk
  
  for (i in 1:nrow(myres)) {
    gc <- rownames(myres)[i]
    pb_sub <- pb[,names(myclusts[myclusts==gc])]
    eigengene <- prcomp(pb_sub)$x[,1]
    
    # need to sign eigengene appropriately
    top_d <- names(eigengene)[order(eigengene,decreasing=TRUE)][1]
    pb_sub[top_d,]
    sum(pb_sub[top_d,]>0)
    eigengene <- pb[,'LL22NC03-2H8.5']
    
    # determine sign of donor score for top_d
    d_sign <- mean(pb_sub[top_d,]) > 0
    if (!d_sign) {
      eigengene <- eigengene * -1
    }
    
    cor(dsc[names(eigengene)],eigengene)
    plot(dsc[names(eigengene)],eigengene)
    
    
  }
  
}


# testing how pca can get you pc values for new samples
pb <- container$scMinimal_ctype[['NK']]$pseudobulk
pres <- prcomp(pb,center=F)
lds1 <- pres[["rotation"]][,1]
pc1 <- pres$x[,1]
samp1 <- pb[1,]
pc1_samp1 <- pc1[1]
print(pc1_samp1)
sum(samp1 * lds1)


### NEED TO KEEP THIS BIT OF CODE SOMEWHERE SO I STOP NEEDING TO RECHECK IT
# rechecking both ways of doing ttm operation
core_new_folded <- k_fold(core_new,m=2,modes=c(10,20,7))
tnsr <- rTensor::as.tensor(tnsr)
core_new_test <- ttl(tnsr,list(t(donor_mat),t(gene_by_factors),t(ctype_by_factors)),ms=c(1,2,3))
dim(core_new_test)
core_new_folded@data[1:10,1,1]
core_new_test@data[1:10,1,1]
all.equal(core_new_test@data,core_new_folded@data)
# this confirms the change that I made indeed does what I expected
###

## now altering the core tensor so I can visualize gene factor links to donor factors
# multiply by ctype factors to get ctypes
df_gf_ct <- ttm(core_new_test,ctype_by_factors,m=3)

# extract panel for donor factor 5
hmap_data <- df_gf_ct@data[10,,]
rownames(hmap_data) <- sapply(1:nrow(hmap_data),function(x) {
  paste0('gene_factor_',x)
}) 
colnames(hmap_data) <- rownames(ctype_by_factors)
# nintieth_per <- stats::quantile(as.matrix(abs(hmap_data)), c(.90))
# nintieth_per <- stats::quantile(as.matrix(abs(hmap_data)), c(.99))
# col_fun = colorRamp2(c(-nintieth_per, 0, nintieth_per), c("blue", "white", "red"))
col_fun = colorRamp2(c(-abs(max(hmap_data)), 0, abs(max(hmap_data))), c("blue", "white", "red"))
Heatmap(hmap_data,col=col_fun)
##


# need to confirm meaning of sign here. It's probably dependent on sign of genes in gene factors and donors in donor factor

## now I want to compute gene factor facing tensor to get corresponding loadings matrices
df_lds_tnsr <- ttl(core_new_test,list(donor_mat,ctype_by_factors),ms=c(1,3))
dim(df_lds_tnsr)
df_lds_tnsr_slice <- df_lds_tnsr@data[,12,] 
dim(df_lds_tnsr_slice)
rownames(df_lds_tnsr_slice) <- rownames(donor_mat)
colnames(df_lds_tnsr_slice) <- rownames(ctype_by_factors)
cor(df_lds_tnsr_slice[,1],donor_mat[,5])
# something seems wrong... I was expecting higher correlations for the things with 
# higher eigenvalues in donor factor gene factor intersection of core tensor

# it's possible that the labels no longer correspond (which is maybe why I used the kronecker method in the first place)
# need to compare df_lds_tnsr to lds tnsr unfolded from fn 
ldngs
ldngs_folded <- k_fold(ldngs,m=2,modes=c(171,20,7))
dim(ldngs_folded)
all.equal(df_lds_tnsr@data,ldngs_folded@data)

# break down a factor from the loadings matrix
genes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][2]})
ctypes <- sapply(colnames(ldngs),function(x){strsplit(x,split=":")[[1]][1]})

sr_col <- ldngs[1,]

tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)

# compare tmp_casted_num to loadings slice
head(ldngs_folded@data[,10,])
head(tmp_casted_num)
sum(round(tmp_casted_num[,5],4) %in% round(ldngs_folded@data[,16,2],4))
## tmp_casted_num shuffles the column order around

cor(donor_mat[rownames(tmp_casted_num),1],tmp_casted_num[,'T8'])
plot(donor_mat[rownames(tmp_casted_num),1],tmp_casted_num[,'T4'])

## why does donor_mat look so different from tucker_results[[1]]?? This should not be affected
## ohhhh it's because normally I reorder the factors by explained variance, so old factor 5 is not same in new

##

# i think it's actually the sex linked factor... let's check the gene factor for the correct genes
dim(gene_by_factors)
cool <- gene_by_factors[,20]
cool[order(cool,decreasing=F)][1:15]
cool[order(cool,decreasing=T)][1:15]
## my suspicions were confirmed! probably then I don't want to clip the eigenvalues quite so much as it looks like there are real processes here.
## though there may actually be some genuine signal since they're all lupus patients and male, so might have some unique signature...

# now see if a derived pc from a gene factor equals its respective column in its loadings matrix
# will use gene factor 12 in T8 to start
pb <- container$scMinimal_ctype[['cM']]$pseudobulk
cool <- gene_by_factors[,12]
cool <- cool[colnames(pb)]
new_pc <- apply(pb,MARGIN=1,function(x) {
  return(sum(x * cool))
})
new_pc[1:10]
# now compare to other one
sr_col <- ldngs[12,]
tmp_casted_num <- reshape_loadings(sr_col,genes,ctypes)
plot(donor_mat[rownames(tmp_casted_num),3],tmp_casted_num[,'NK'])
cor(donor_mat[names(new_pc),4],new_pc)

rssq <- sqrt(sum(new_pc**2))
new_pc2 <- new_pc / rssq
sum(new_pc2**2)
cor(donor_mat[names(new_pc),8],new_pc2)
var(new_pc2)
new_pc2 <- scale(new_pc2)

# am starting to think the gene factors are informative but using the correlations is spurious
# because they can easily be influenced by the many non-zero but non-significant genes in a
# gene factor. The sparse method seems to help a bit but takes forever


# I'm going to try my idea to use wgcna to get gene modules and then test for 
# significance of each module eigenvalue against each donor factor
library(WGCNA)
datExpr <- container$scMinimal_ctype[['cM']]$pseudobulk

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

net = blockwiseModules(datExpr, power = 2, maxBlockSize = 10000,
                       TOMType = "unsigned", minModuleSize = 15,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


factor_select <- 5
dsc <- container$tucker_results[[1]][,factor_select]

MEs <- net$MEs
col_ndx <- sapply(colnames(MEs),function(x) {
  as.numeric(strsplit(x,split='ME')[[1]][[2]])
})
MEs <- MEs[,order(col_ndx)]
head(MEs)
MEs[,1] <- NULL
head(MEs)

for (i in 1:ncol(MEs)) {
  ME <- MEs[,i]
  names(ME) <- rownames(MEs)
  myc <- cor(dsc[names(ME)],ME)
  print(i)
  print(myc)
}

net$colors[net$colors==4]


## more formally running LR analysis with new functions
# prep for new LR analysis
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
  rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
  return(rname)
})
lr_pairs$interaction_name <- NULL

container <- pbmc_container
container <- prep_LR_interact(container, lr_pairs, norm_method='trim', scale_factor=10000,
                              var_scale_power=1.5, batch_var='pool')
sft_thresh <- c(3,3,2,2,2,2,2)
container <- get_gene_modules(container,sft_thresh)

container <- compute_LR_interact(container, lr_pairs, factor_select=5, sig_thresh=0.05, percentile_exp_rec=.9)
container$plots$lr_analysis[['Factor5']]

container <- compute_LR_interact(container, lr_pairs, factor_select=4, sig_thresh=0.05, percentile_exp_rec=.9)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_f4.pdf", useDingbats = FALSE,
    width = 15, height = 13)
container$plots$lr_analysis[['Factor4']] # ICAM2 find demonstrates recapitulation of known biology
dev.off()

container <- compute_LR_interact(container, lr_pairs, factor_select=7, sig_thresh=0.05, percentile_exp_rec=.9)
container$plots$lr_analysis[['Factor7']]

container <- compute_LR_interact(container, lr_pairs, factor_select=2, sig_thresh=0.05, percentile_exp_rec=.9)
container$plots$lr_analysis[['Factor2']]

container <- compute_LR_interact(container, lr_pairs, factor_select=1, sig_thresh=0.05, percentile_exp_rec=.9)
container$plots$lr_analysis[['Factor1']]

container <- compute_LR_interact(container, lr_pairs, factor_select=3, sig_thresh=0.05, percentile_exp_rec=.9)
container$plots$lr_analysis[['Factor3']]

padj <- get_module_enr(container,ctype='cM',mod_select=4,db_use='GO')
padj <- get_module_enr(container,ctype='cM',mod_select=4,db_use='immuno')
padj <- get_module_enr(container,ctype='T4',mod_select=7,db_use='GO')
padj <- get_module_enr(container,ctype='T4',mod_select=7,db_use='immuno')
padj <- get_module_enr(container,ctype='T4',mod_select=2,db_use='TF')
padj <- get_module_enr(container,ctype='NK',mod_select=2,db_use='TF')

# looking at group of modules highly specifically correlated with the factor
padj <- get_module_enr(container,ctype='B',mod_select=7,db_use='TF')
padj <- get_module_enr(container,ctype='NK',mod_select=8,db_use='TF')
padj <- get_module_enr(container,ctype='cM',mod_select=4,db_use='TF')
padj <- get_module_enr(container,ctype='T4',mod_select=7,db_use='TF')
padj <- get_module_enr(container,ctype='cDC',mod_select=3,db_use='TF')
padj <- get_module_enr(container,ctype='ncM',mod_select=2,db_use='TF')
## we could be looking at a modified IFN program mediated by NCOA2
## I wonder if I would have caught this in my gsea on the loadings if I used TF db

padj[order(padj,decreasing=FALSE)][1:10]


ctypes <- c('NK','cM','T4','cDC','ncM')
modules <- c(8,4,7,3,2)

mod_enr <- plot_multi_module_enr(container, ctypes, modules, sig_thresh=.05, db_use='TF')
mod_enr

mod_enr <- plot_multi_module_enr(container, ctypes, modules, sig_thresh=.05, db_use=c('GO'))
mod_enr

mod_enr <- plot_multi_module_enr(container, ctypes, modules, sig_thresh=.05, db_use=c('Hallmark'))
mod_enr

ctypes <- c('NK','cDC','ncM','cM','T4','T8')
modules <- c(2,1,1,1,2,3)

mod_enr <- plot_multi_module_enr(container, ctypes, modules, sig_thresh=.1, db_use=c('TF'))
mod_enr

mod_enr <- plot_multi_module_enr(container, ctypes, modules, sig_thresh=.05, db_use=c('Hallmark'))
mod_enr



# dsc <- container$tucker_results[[1]][,5]
# tmp <- as.data.frame(cbind(dsc,container$scMinimal_ctype[['ncM']]$pseudobulk[,'TNF']))
# colnames(tmp) <- c('dsc','expr')
# ggplot(tmp,aes(x=dsc,y=expr)) +
#     geom_point()


# look at levels of ncM TNF in men vs women of the top donors for factor 5
dsc <- container$tucker_results[[1]][,5]
top_d <- dsc[order(dsc,decreasing=FALSE)][1:20]
meta <- container$scMinimal_ctype[['ncM']]$metadata
meta <- meta[,c('donors','sex','Status')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

male_don <- rownames(meta)[meta$sex=='Male']
female_don <- rownames(meta)[meta$sex=='Female']
male_don <- male_don[male_don %in% names(top_d)]
female_don <- female_don[female_don %in% names(top_d)]

m_exp <- container$scMinimal_ctype[['ncM']]$pseudobulk[male_don,'TNF']
f_exp <- container$scMinimal_ctype[['ncM']]$pseudobulk[female_don,'TNF']

tot_exp <- c(m_exp,f_exp)
add_on <- c(rep('M',length(m_exp)),rep('F',length(f_exp)))
m_f_exp <- cbind(tot_exp,add_on)
m_f_exp <- as.data.frame(m_f_exp)
m_f_exp$tot_exp <- as.numeric(m_f_exp$tot_exp)
# run t-test
t.test(tot_exp~add_on)


# checking TNFSF8 in T4
m_exp <- container$scMinimal_ctype[['T4']]$pseudobulk[male_don,'TNFSF8']
f_exp <- container$scMinimal_ctype[['T4']]$pseudobulk[female_don,'TNFSF8']

tot_exp <- c(m_exp,f_exp)
add_on <- c(rep('M',length(m_exp)),rep('F',length(f_exp)))
m_f_exp <- cbind(tot_exp,add_on)
m_f_exp <- as.data.frame(m_f_exp)
m_f_exp$tot_exp <- as.numeric(m_f_exp$tot_exp)
# run t-test
t.test(tot_exp~add_on)

cool <- container$scMinimal_ctype[['ncM']]$pseudobulk[names(top_d),'NCOA2']
'NCOA2' %in% rownames(container$scMinimal_full$count_data)
'NR3C1' %in% rownames(container$scMinimal_full$count_data)
'ESR1' %in% rownames(container$scMinimal_full$count_data)

# it wasnt in most variable genes so not in pb data, need to recompute it
container <- get_pseudobulk(container)
container <- normalize_pseudobulk(container, method='trim', scale_factor=10000)
"NCOA2" %in% colnames(container[["scMinimal_ctype"]][["cM"]]$pseudobulk)
tmp <- as.data.frame(cbind(container[["scMinimal_ctype"]][["cM"]]$pseudobulk[names(dsc),'NCOA2'],dsc))
colnames(tmp) <- c('expr','dsc')
plot(tmp$dsc,tmp$expr)

tmp <- as.data.frame(cbind(container[["scMinimal_ctype"]][["cM"]]$pseudobulk[names(dsc),'NR3C1'],dsc))
colnames(tmp) <- c('expr','dsc')
plot(tmp$dsc,tmp$expr)

# restrict it to just lupus patients
meta_l <- meta[meta$Status=='Managed',]
dsc_l <- dsc[rownames(meta_l)]

tmp <- as.data.frame(cbind(container[["scMinimal_ctype"]][["NK"]]$pseudobulk[names(dsc_l),'NCOA2'],dsc_l))
colnames(tmp) <- c('expr','dsc')
plot(tmp$dsc,tmp$expr)

tmp <- as.data.frame(cbind(container[["scMinimal_ctype"]][["T8"]]$pseudobulk[names(dsc_l),'NR3C1'],dsc_l))
colnames(tmp) <- c('expr','dsc')
plot(tmp$dsc,tmp$expr)

tmp <- as.data.frame(cbind(t(container[["scMinimal_ctype"]][["ncM"]]$pseudobulk)[names(dsc_l),'NR3C1'],dsc_l))
colnames(tmp) <- c('expr','dsc')
plot(tmp$dsc,tmp$expr)

# pretty tough to tell if there is any real signal...

# perhaps of the donors with high ifn signaling they have less?
isg_lev <- container$scale_pb_extra[["T4"]][names(dsc_l),'ISG15']
isg_lev <- isg_lev[order(isg_lev,decreasing=TRUE)][1:100]
don_lup_ifn <- names(isg_lev)

tmp <- as.data.frame(cbind(container[["scMinimal_ctype"]][["cDC"]]$pseudobulk[don_lup_ifn,'NCOA2'],dsc_l[don_lup_ifn]))
colnames(tmp) <- c('expr','dsc')
plot(tmp$dsc,tmp$expr)

# It seems like perhaps there is an enrichment of the factor 5 donors for having low NCOA2 NR3C1 but there is not a striking trend


# running gsea to see what functionally distinguishes these cd4 subtype 2 genes
# ITGB1 S100A4 S100A10 ANXA1 LGALS1 KLRB1 CRIP1 SH3BGRL3
mod_genes <- c('ITGB1', 'S100A4', 'S100A10', 'ANXA1', 'LGALS1', 'KLRB1', 'CRIP1', 'SH3BGRL3')
db_use='immuno'
pvals[order(pvals,decreasing=FALSE)][1:10]
## top 3 hits are all naive vs memory cd4 cell down, so these are likely markers of memory cd4 cells

# now will try for the T4_3 groups of genes, then all together
mod_genes <- c('JUN', 'JUNB', 'CD69', 'FOS', 'DUSP1', 'CXCR4', 'TSC22D3', 'LEPROTL1', 'ZFP36L2', 'PNRC1')
db_use='immuno'
db_use='ctype'
pvals[order(pvals,decreasing=FALSE)][1:30]

mod_genes <- c('ITGB1', 'S100A4', 'S100A10', 'ANXA1', 'LGALS1', 'KLRB1', 'CRIP1', 'SH3BGRL3','JUN', 'JUNB', 'CD69', 'FOS', 'DUSP1', 'CXCR4', 'TSC22D3', 'LEPROTL1', 'ZFP36L2', 'PNRC1')
db_use='immuno'
pvals[order(pvals,decreasing=FALSE)][1:30]


# seeing if age is correlated with dscore
dsc <- pbmc_container$tucker_results[[1]][,4]
meta <- pbmc_container$scMinimal_ctype[['ncM']]$metadata
meta <- meta[,c('donors','sex','Status','age')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL
meta$age <- as.numeric(as.character(meta$age))
cor(dsc,meta[names(dsc),'age'])
tmp <- as.data.frame(cbind(dsc,meta[names(dsc),c('age','Status')]))
colnames(tmp) <- c('dsc','age','Status')
lmres <- lm(age~dsc,data=tmp)
summary(lmres)
ggplot(tmp,aes(x=dsc,y=age,color=Status)) +
  geom_point()






### testing whether it's possible to plot gene callouts on loadings plots with different colors
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=5, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.05, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, show_xlab=TRUE,
                                      show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)

# extracting leading edge genes from GSEA
gsea_res <- pbmc_container[["gsea_res_full"]][["Factor5"]][["T4"]]
gsea_res$padj[1:100]
gsea_res$pathway[1:100]
callouts <- gsea_res$leadingEdge[gsea_res$pathway=="GO_RESPONSE_TO_TYPE_I_INTERFERON"][[1]][1:5]
callouts <- c(callouts,gsea_res$leadingEdge[gsea_res$pathway=="GO_APOPTOTIC_PROCESS"][[1]][1:5])
callouts <- c(callouts,gsea_res$leadingEdge[gsea_res$pathway=="GO_REGULATION_OF_CELL_ADHESION"][[1]][1:5])


gsea_res <- pbmc_container[["gsea_res_full"]][["Factor5"]][["cM"]]
gsea_res$pathway[1:60]
callouts <- c(callouts,gsea_res$leadingEdge[gsea_res$pathway=="GO_EXOCYTOSIS"][[1]][1:5])

tmp_casted_num <- pbmc_container$plots$lds_plots_data[['5']]
ndx <- match(callouts,rownames(tmp_casted_num))
# check that no genes appear in multple sets
if (length(unique(ndx))!=length(ndx)) {
  print('there are duplicates')
}

# if there are duplicates, need to get index of all instances and remove both them and their elements from color vector

ann_hmap <- pbmc_container[["plots"]][["all_lds_plots"]][["5"]] + 
  rowAnnotation(callouts = anno_mark(at = ndx, which='row',
                                     labels = callouts,
                                     labels_gp = gpar(col = c(rep('green',5),rep('orange',5)))))

# check that I can add legend back to heatmap
pd <- pbmc_container$plots$all_legends[['5']]

# optionally draw the plot
draw(ann_hmap,annotation_legend_list = pd,
     legend_grouping = "original",
     newpage=FALSE)

# make a legend to show gene set links to the thing of the thing of course
mycols <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_EXOCYTOSIS")
names(mycols) <- c('green','orange')
draw(Legend(labels = mycols, legend_gp = gpar(fill = names(mycols)), title = "gene_sets",
            grid_height = unit(1, "mm"), grid_width = unit(3, "mm")))











### very quickly testing my idea of gene set splitting using the IFN gene set enriched in all or most
# first, chekc that it's in all ctypes
gs <- "GO_RESPONSE_TO_TYPE_I_INTERFERON"
gs <- "GO_RESPONSE_TO_CYTOKINE"
for (ct in pbmc_container$experiment_params$ctypes_use) {
  gsea_res <- pbmc_container[["gsea_res_full"]][["Factor5"]][[ct]]
  print(gsea_res$padj[gsea_res$pathway==gs])
}

# get number of times each leading edge gene shows up
le_counts <- list()
for (ct in pbmc_container$experiment_params$ctypes_use) {
  gsea_res <- pbmc_container[["gsea_res_full"]][["Factor5"]][[ct]]
  le_genes <- gsea_res$leadingEdge[gsea_res$pathway==gs][[1]]
  for (g in le_genes) {
    if (g %in% names(le_counts)) {
      le_counts[[g]] <- le_counts[[g]] + 1
    } else {
      le_counts[[g]] <- 1
    }
  }
}

# compute fractions of unique genes for each ctype
for (ct in pbmc_container$experiment_params$ctypes_use) {
  gsea_res <- pbmc_container[["gsea_res_full"]][["Factor5"]][[ct]]
  le_genes <- gsea_res$leadingEdge[gsea_res$pathway==gs][[1]]
  num_unq <- 0
  for (g in le_genes) {
    if (le_counts[[g]] == 1) {
      num_unq <- num_unq + 1
    }
  }
  print(ct)
  print(num_unq / length(le_genes))
}

# most genes are found in the leading edge genes of at least one other cell type as expected




# testing whether some other exhaustion genes are enriched in F5
dsc <- container$tucker_results[[1]]
tmp <- dsc[,5]

pb <- container$scMinimal_ctype[['T4']]$pseudobulk

tmp2 <- cbind(tmp,pb[names(tmp),'ITGB1'])
plot(tmp2[,1],tmp2[,2])
# doesnt really appear to be significant










# compare similar factors to find out what's similar/different
pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,5), direction=c('up','down'),
                                  compare_type='same', sig_thresh=0.02)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_f2_f5_same_up.pdf", useDingbats = FALSE,
    width = 10, height = 14)
pbmc_container$plots$comparisons[['2_5']]
dev.off()

pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,5), direction=c('up','down'),
                                  compare_type='different', sig_thresh=0.02)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_f2_f5_different_up.pdf", useDingbats = FALSE,
    width = 10, height = 14)
pbmc_container$plots$comparisons[['2_5']]
dev.off()

pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,5), direction=c('down','up'),
                                  compare_type='same', sig_thresh=0.02)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_f2_f5_same_down.pdf", useDingbats = FALSE,
    width = 10, height = 14)
pbmc_container$plots$comparisons[['2_5']]
dev.off()

pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,5), direction=c('down','up'),
                                  compare_type='different', sig_thresh=0.02)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/lupus_f2_f5_different_down.pdf", useDingbats = FALSE,
    width = 10, height = 35)
pbmc_container$plots$comparisons[['2_5']]
dev.off()


pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,5), direction=c('down','up'),
                                  compare_type='different', sig_thresh=0.05)
pbmc_container$plots$comparisons[['2_5']]
diff_enr <- get_compare_go_enrich(pbmc_container,'cM',-1)
print(diff_enr[order(diff_enr,decreasing=F)][1:10])


pbmc_container <- compare_factors(pbmc_container, f_compare=c(2,4), direction=c('up','down'),
                                  compare_type='different', sig_thresh=0.05)
pbmc_container$plots$comparisons[['2_4']]
diff_enr <- get_compare_go_enrich(pbmc_container,'T8',-1)
print(diff_enr[order(diff_enr,decreasing=F)][1:15])




## testing out new function to get leading edge genes
gsets <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_EXOCYTOSIS")
cool <- get_leading_edge_genes(pbmc_container,factor_select=5,gsets,num_genes_per=5)




# this is very strange that HALLMARK_TNFA_SIGNALING_VIA_NFKB is up in F5 cM!! This is the exact opposite of what I saw in the module overrepresentation analysis...
# This could be explained two ways potentially:
# 1. maybe the specific cM module is actually of genes which are downregulated, like the signs are somehow mixed up, or maybe the module contains a mix of up AND down genes
# 2. maybe TNFA is sort of significant on both ends of F5 with the module contributing to the "wrong" end, but still significant in the module compared to rest...
# it's probably a good sign that we see cM behave separately from the rest of the ctypes with respect to NFKB signaling, but need to look into this as it's pretty critical
# for me making the argument later on that my LR analysis is working as intended.

# to check on these things I'll first look into the factor 5 specific cM module to check on directionality of the genes

### this code is copied from above
# prep for new LR analysis
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
  rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
  return(rname)
})
lr_pairs$interaction_name <- NULL

pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=1.5, batch_var='pool')
sft_thresh <- c(3,3,2,2,2,2,2)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=5, sig_thresh=0.05, percentile_exp_rec=.9)

pdf(file = "/home/jmitchel/figures/for_paper/lupus_v2/LR_f5.pdf", useDingbats = FALSE,
    width = 15, height = 13)
pbmc_container$plots$lr_analysis[['Factor5']]
dev.off()

cm_mod <- pbmc_container[["module_genes"]][["cM"]]
cool <- names(cm_mod)[cm_mod==4]

# seeing if TNFA leading edge genes in the module
test_g <- c('IL1B','PNRC1','CD83','PLEK','IER5','BTG2')
test_g %in% cool

# or test all le genes
le[[1]] %in% cool

dsc <- pbmc_container$tucker_results[[1]][,5]
MEs <- pbmc_container[["module_eigengenes"]][["cM"]]
head(MEs)
me <- MEs[,4]
names(me) <- rownames(MEs)
plot(dsc[names(me)],me)

# checking direction of expression for individual genes in module
pb <- pbmc_container$scMinimal_ctype[['cM']]$pseudobulk
plot(dsc[names(me)],pb[names(me),cool[13]])
plot(dsc[names(me)],pb[names(me),'KLF10'])
# module seems to be bidirectional

# see if TNF genes in all the up ones or not

cool[13]

m_df <- data.frame()
m_df <- rbind(m_df,msigdbr::msigdbr(species = "Homo sapiens",
                                    category = "H"))
my_pathways <- split(m_df$gene_symbol, f = m_df$gs_name)
cool %in% my_pathways[['HALLMARK_TNFA_SIGNALING_VIA_NFKB']]

print(cool[cool %in% my_pathways[['HALLMARK_TNFA_SIGNALING_VIA_NFKB']]])
plot(dsc[names(me)],pb[names(me),"RHOB"])

# need to make a NES plot to see that the genes fall close to either end of the spectrum...
fgsea::plotEnrichment(my_pathways[['HALLMARK_TNFA_SIGNALING_VIA_NFKB']],
                      exp_vals)
# this confirms my suspicion that there are TNF genes on either end, however there are indeed more
# on the positive end, meaning that ?
exp_vals['RHOB']

ctypes <- rep('cM',4)
modules <- c(1,3,4,5)
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('Hallmark'))
mod_enr

ctypes <- rep('cM',6)
modules <- c(1,3,4,5,6,7)
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('Hallmark'))
mod_enr


# trying regular gsea with hallmark sets
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=5, method="fgsea", thresh=0.05,
                                      db_use=c("Hallmark"), collapse_paths=FALSE)
tmp <- pbmc_container[["gsea_res_full"]][["Factor5"]][["cM"]]
le <- tmp[tmp$pathway=='HALLMARK_TNFA_SIGNALING_VIA_NFKB','leadingEdge'][[1]]


# look for modules with highest correlations to F2 dscores with just lupus patients
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Status')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL

dsc <- pbmc_container$tucker_results[[1]][,2]
tmp <- cbind.data.frame(dsc,meta[names(dsc),])
colnames(tmp) <- c('dsc','Status')
tmp <- tmp[tmp$Status=='Managed',]
dsc <- tmp[,1,drop=FALSE]

mod_fact_r_all <- list()
ME_pvals <- list()
for (ct in pbmc_container$experiment_params$ctypes_use) {
  MEs <- pbmc_container$module_eigengenes[[ct]]
  for (j in 1:ncol(MEs)) {
    ME <- MEs[,j,drop=FALSE]
    
    # calculate significance of association with factor as well as Rsq
    tmp <- as.data.frame(cbind(dsc,ME[rownames(dsc),]))
    colnames(tmp) <- c('dsc','eg')
    lmres <- lm(dsc~eg,data=tmp)
    lmres <- summary(lmres)
    pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
    mod_fact_r <- cor(tmp[,1],tmp[,2])
    mod_fact_r_all[[paste0(ct,'_',j)]] <- mod_fact_r
    ME_pvals[[paste0(ct,"_",as.character(j))]] <- pval
  }
}

