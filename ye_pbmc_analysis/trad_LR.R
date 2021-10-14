

con <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_conos2.rds')

con$plotGraph()

### can try reducing number of cells by randomly downsampling
# # save original embedding
# orig_embed <- con[["embedding"]]
#
# # save original cluster labels
# orig_clusts <- con$clusters$leiden$groups

# downsample cells
cells_keep <- sample(rownames(orig_embed),300000)

con$clusters$leiden$groups <- orig_clusts[cells_keep]
con[["embedding"]] <- orig_embed[cells_keep,]

# get subtype DE results heamap
myde <- con$getDifferentialGenes(groups=con$clusters$leiden$groups,append.auc=TRUE,z.threshold=0,upregulated.only=TRUE)




#### now identify number of interactions by co-upregulated LR pairs

# load up the LR database
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
  rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
  return(rname)
})
lr_pairs$interaction_name <- NULL

## loop through cell types to identify source ligands
all_lig_hits <- list()
for (ct in pbmc_container$experiment_params$ctypes_use) {
  de_res <- myde[[ct]]
  de_res_sig <- de_res[de_res$PAdj<.05,]
  de_res_sig <- de_res_sig[de_res_sig$AUC>.55,]

  lig_hits <- rownames(de_res_sig)[rownames(de_res_sig) %in% lr_pairs[,1]]
  all_lig_hits[[ct]] <- lig_hits
}



# maybe I should do this per individual

# loop through individuals

# subset to their cells

# run DE

# get putative interactions


# now get donor status meta data
meta <- pbmc_container$scMinimal_full$metadata[,c('donors','Status','Age')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL
head(meta)

donors <- rownames(meta)
for (d in donors) {
  cells_use <- rownames(pbmc_container$scMinimal_full$metadata)[pbmc_container$scMinimal_full$metadata$donors==d]
}









## trying with CellChat to see if it gives me more results
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')

# reduce to cells of the cell types that I'm using
# cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$cg_cov %in% c("B","NK","T4","T8","cDC","cM","ncM")]
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$cg_cov %in% c("B","NK","Th","Tc","cDC","cMono","ncMono")]

cells_keep <- sample(cells_keep,50000)
pbmc <- subset(pbmc,cells = cells_keep)
# pbmc@meta.data$cg_cov <- factor(pbmc@meta.data$cg_cov,levels=c("B","NK","T4","T8","cDC","cM","ncM"))
pbmc@meta.data$cg_cov <- factor(pbmc@meta.data$cg_cov,levels=c("B","NK","Th","Tc","cDC","cMono","ncMono"))
pbmc@meta.data$idents <- pbmc@meta.data$cg_cov
Idents(pbmc) <- pbmc@meta.data$idents

cellchat <- createCellChat(object = pbmc, group.by = "cg_cov")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)

# future::plan("multiprocess", workers = 15) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)


cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

df.net <- subsetCommunication(cellchat)
df.net <- subsetCommunication(cellchat,slot.name = "netP")
df.net[1:50,]
unique(df.net$ligand)
max(df.net$pval)

df.net[df.net$ligand=='IL16',]

# saveRDS(df.net,file='/home/jmitchel/data/lupus_data/LR_cell_chat_res.rds')
# saveRDS(df.net,file='/home/jmitchel/data/lupus_data/LR_cell_chat_res_all_together.rds')




## will try to run it on each donor since I didn't get results like TNFSF13B which are to be expected
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')

CellChatDB <- CellChatDB.human # the LR database to use

# reduce to cells of the cell types that I'm using
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$cg_cov %in% c("B","NK","T4","T8","cDC","cM","ncM")]
pbmc <- subset(pbmc,cells = cells_keep)
pbmc@meta.data$cg_cov <- factor(pbmc@meta.data$cg_cov,levels=c("B","NK","T4","T8","cDC","cM","ncM"))

# get vec of donors
donors <- rownames(pbmc_container$tucker_results[[1]])

# loop through donors and get LR interactions and store in list
all_df_net <- mclapply(1:length(donors), function(i) {
  d <- donors[i]
  cells_use <- rownames(pbmc@meta.data)[pbmc@meta.data$ind_cov_batch_cov==d]
  pbmc_sub <- subset(pbmc,cells = cells_use)
  
  cellchat <- createCellChat(object = pbmc_sub, group.by = "cg_cov")
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  return(df.net)
}, mc.cores=20)


# combine lists
all_df_net_comb <- do.call("rbind", all_df_net)

# add donor labels to the list version
names(all_df_net) <- donors

# saveRDS(all_df_net,file='/home/jmitchel/data/lupus_data/LR_cell_chat_res_d_sep_list.rds')
# saveRDS(all_df_net_comb,file='/home/jmitchel/data/lupus_data/LR_cell_chat_res_d_sep.rds')
all_df_net_comb <- readRDS(file='/home/jmitchel/data/lupus_data/LR_cell_chat_res_d_sep.rds')


head(all_df_net_comb)

unique(all_df_net_comb$ligand)

all_df_net_comb_sub <- all_df_net_comb[,c(1,2,3,4,6)]

all_df_net_comb_sub <- unique(all_df_net_comb_sub)
dim(all_df_net_comb_sub)

head(all_df_net_comb_sub)
all_df_net_comb_sub[all_df_net_comb_sub$ligand=='TNFSF13B',]

tnf_yes <- c()
dsc <- c()
for (i in 1:length(all_df_net)) {
  tmp <- all_df_net[[i]]
  # print('TNFSF13B' %in% tmp$ligand)
  # print(pbmc_container$tucker_results[[1]][names(all_df_net)[i],1])
  
  tnf_yes <- c(tnf_yes,'TNFSF13B' %in% tmp$ligand)
  dsc <- c(dsc,pbmc_container$tucker_results[[1]][names(all_df_net)[i],1])
  
}

my_ord <- order(dsc)
dsc[my_ord]
tnf_yes[my_ord]



## trying out nichenet ligand to target stuff
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

# seeing if ICOSLG results are enriched for cell cycle and T signaling genes
test = ligand_target_matrix[,'ICOSLG']
test = ligand_target_matrix[,'TNFSF13B']
test = ligand_target_matrix[,'TNF']
test2 <- test[order(test,decreasing=TRUE)][1:50]

library(enrichR)
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
enriched <- enrichr(names(test2),dbs)




### trying to scale up my nichenet analysis
###
myres_mat <- pbmc_container$lr_res #includes NS results still but are FDR adjusted

# reduce to rows/columns with at least one significant hit
myres_mat <- myres_mat[rowSums(myres_mat<.05)>0,]
myres_mat <- myres_mat[,colSums(myres_mat<.05)>0]

myres_mat <- myres_mat[rowSums(myres_mat<.04)>0,]
myres_mat <- myres_mat[,colSums(myres_mat<.04)>0]

myres_mat <- myres_mat[rowSums(myres_mat<.0001)>0,]
myres_mat <- myres_mat[,colSums(myres_mat<.0001)>0]

myres_mat <- myres_mat[rowSums(myres_mat<.001)>0,]
myres_mat <- myres_mat[,colSums(myres_mat<.001)>0]

myres_mat <- myres_mat[rowSums(myres_mat<.000000000001)>0,]
myres_mat <- myres_mat[,colSums(myres_mat<.000000000001)>0]

## get unique ligand_source_target channels from my analysis
unique_channels <- c()
for (i in 1:nrow(myres_mat)) {
  lig_ct_rec <- myres_mat[i,]
  lig_ct_rec_name <- strsplit(rownames(myres_mat)[i],split='_')[[1]]
  lig <- lig_ct_rec_name[[1]]
  source <- lig_ct_rec_name[[2]]
  
  for (j in 1:ncol(myres_mat)) {
    pv <- lig_ct_rec[j]
    if (pv < 0.001) {
      target_ct <- strsplit(names(pv),split="_")[[1]][[1]]
      lig_source_target <- paste0(lig,"_",source,"_",target_ct)
      unique_channels <- c(unique_channels,lig_source_target)
    }
  }
}
unique_channels <- unique(unique_channels)
my_channels <- unique_channels

# get unique channels for the cellchat results
all_df_net_comb <- readRDS(file='/home/jmitchel/data/lupus_data/LR_cell_chat_res_all_together.rds')
# all_df_net_comb <- readRDS(file='/home/jmitchel/data/lupus_data/LR_cell_chat_res_d_sep.rds')
all_df_net_comb_sub <- all_df_net_comb[,c(1,2,3,6)]
sum(all_df_net_comb_sub$pval>0.04)
all_df_net_comb_sub$pval <- NULL
all_df_net_comb_sub <- unique(all_df_net_comb_sub)
dim(all_df_net_comb_sub)

# update source and target names to match ones I use
unique(all_df_net_comb_sub$source)
unique(all_df_net_comb_sub$target)

all_df_net_comb_sub$source <- as.character(all_df_net_comb_sub$source)
new_source <- sapply(all_df_net_comb_sub$source,function(x) {
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
names(new_source) <- NULL
all_df_net_comb_sub$source <- new_source

all_df_net_comb_sub$target <- as.character(all_df_net_comb_sub$target)
new_target <- sapply(all_df_net_comb_sub$target,function(x) {
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
names(new_target) <- NULL
all_df_net_comb_sub$target <- new_target

cc_channels <- sapply(1:nrow(all_df_net_comb_sub),function(i){
  source <- all_df_net_comb_sub[i,1]
  target <- all_df_net_comb_sub[i,2]
  ligand <- all_df_net_comb_sub[i,3]
  if (source != target) {
    return(paste0(ligand, "_",
                  source, "_",
                  target))
  } else {
    return(NA)
  }
  
})
cc_channels <- cc_channels[!is.na(cc_channels)]
length(cc_channels)
cc_channels[1:10]

# look at overlap between the two sets
length(cc_channels)
length(my_channels)
c_both <- intersect(cc_channels,my_channels)
length(c_both)


'TNFSF13B_cMono_B' %in% my_channels
'ICOSLG_cMono_Th' %in% my_channels
'TNFSF13B_cMono_B' %in% cc_channels
'ICOSLG_cMono_Th' %in% cc_channels


# library(spqn) # consider using this to improve gene associations...

# now loop through channels get associated genes and test 
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

test_lr_validity <- function(container,lr_test,ligand_target_matrix) {
  # all_pvals <- c()
  print(length(lr_test))
  all_pvals <- sccore::plapply(1:length(lr_test), function(lst_ndx) {
  # all_pvals <- lapply(1:length(lr_test), function(lst_ndx) {
    lig_source_target <- lr_test[lst_ndx]
    lst_split <- strsplit(lig_source_target,split='_')[[1]]
    lig <- lst_split[[1]]
    source_ct <- lst_split[[2]]
    target_ct <- lst_split[[3]]
    
    lig_ct_exp <- container$scale_pb_extra[[source_ct]][,lig]
    
    all_cors <- cor_helper(container,lig_ct_exp,lig,source_ct,target_ct)

    # # testing using wilcoxin test
    # if (lig %in% colnames(ligand_target_matrix)) {
    #   gene_intersect <- intersect(names(all_cors),rownames(ligand_target_matrix))
    #   all_cors <- all_cors[gene_intersect]
    # 
    #   top_cors <- all_cors[order(abs(all_cors),decreasing=TRUE)][1:200]
    #   # top_cors <- all_cors[order(all_cors,decreasing=TRUE)][1:200]
    #   tmp <- cbind.data.frame(all_cors,ligand_target_matrix[gene_intersect,lig])
    #   colnames(tmp) <- c('mycors','target_scores')
    #   tmp$top <- sapply(1:nrow(tmp),function(i){
    #     if (rownames(tmp)[i] %in% names(top_cors)) {
    #       return('top_gn')
    #     } else {
    #       return('other')
    #     }
    #   })
    #   tmp$top <- as.factor(tmp$top)
    #   my_res <- wilcox.test(target_scores~top,data=tmp)
    #   pval <- my_res$p.value
    #   # pval <- sig_test_helper(tmp,1000)
    #   return(pval)
    # } else {
    #   return(NA)
    # }
    
    # testing using linear model
    if (lig %in% colnames(ligand_target_matrix)) {
      gene_intersect <- intersect(names(all_cors),rownames(ligand_target_matrix))
      # tmp <- cbind.data.frame(all_cors[gene_intersect],ligand_target_matrix[gene_intersect,lig])
      tmp <- cbind.data.frame(abs(all_cors[gene_intersect]),ligand_target_matrix[gene_intersect,lig])
      colnames(tmp) <- c('mycors','target_scores')
      lmres <- summary(lm(mycors~target_scores,data=tmp))
      pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
      # pval <- sig_test_helper(tmp,10000)
      return(pval)
    } else {
      return(NA)
    }
    
    # # testing using an enrichment test
    # if (lig %in% colnames(ligand_target_matrix)) {
    #   gene_intersect <- intersect(names(all_cors),rownames(ligand_target_matrix))
    #   all_cors <- all_cors[gene_intersect]
    #   tmp <- cbind.data.frame(all_cors,ligand_target_matrix[gene_intersect,lig])
    #   colnames(tmp) <- c('mycors','target_scores')
    # 
    #   # top 100 genes by predicted downstream perturbation
    #   mypaths <- list()
    #   mypaths[['signal_genes']] <- rownames(tmp)[order(tmp$target_scores,decreasing=TRUE)[1:100]]
    #   myranks <- all_cors
    #   fgseaRes <- fgsea::fgseaSimple(pathways = mypaths,
    #                            stats    = myranks,
    #                            minSize  = 0,
    #                            maxSize  = 5000,
    #                            nperm = 10000)
    #   return(fgseaRes$pval[1])
    # } else {
    #   return(NA)
    # }
  # })
  },mc.preschedule=TRUE,n.cores=30,progress=TRUE)
  
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



library(spqn) 
cor_helper <- function(container,lig_ct_exp,lig,source_ct,target_ct) {
  test <- cbind.data.frame(lig_ct_exp,container$scale_pb_extra[[target_ct]][names(lig_ct_exp),])
  cor_m <- cor(test)
  
  lig_ct_exp2 <- container$no_scale_pb_extra[[source_ct]][,lig]
  test2 <- cbind.data.frame(lig_ct_exp2,container$no_scale_pb_extra[[target_ct]][names(lig_ct_exp2),])
  ave_exp <- colMeans(test2)
  
  cor_m_spqn <- normalize_correlation(cor_m, ave_exp=ave_exp, ngrp=20, size_grp=300, ref_grp=5)
  # cor_m_spqn <- normalize_correlation(cor_m, ave_exp=ave_exp, ngrp=20, size_grp=300, ref_grp=18)
  all_cors <- cor_m_spqn[1,]
  names(all_cors) <- colnames(cor_m)
  all_cors <- all_cors[2:length(all_cors)]
  return(all_cors)
}

sig_test_helper <- function(tmp,n_iter) {
  ## to use with lm version
  # lmres <- summary(lm(mycors~target_scores,data=tmp))
  # obs_fstat <- lmres$fstatistic[1]
  # 
  # all_fstats <- c()
  # for (i in 1:n_iter) {
  #   tmp2 <- transform(tmp, target_scores = sample(target_scores) )
  #   lmres <- summary(lm(mycors~target_scores,data=tmp2))
  #   fstat <- lmres$fstatistic[1]
  #   all_fstats <- c(all_fstats,fstat)
  # }
  # pval <- sum(all_fstats>=obs_fstat) / length(all_fstats)
  
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

cc_channels <- readRDS(file='/home/jmitchel/data/lupus_data/LR_cchat_res_orig.rds')

cc_res <- test_lr_validity(pbmc_container,cc_channels,ligand_target_matrix)
scITD_res <- test_lr_validity(pbmc_container,my_channels,ligand_target_matrix)
scITD_res2 <- test_lr_validity(pbmc_container,my_channels,ligand_target_matrix)
italk_res <- test_lr_validity(pbmc_container,italk_channels,ligand_target_matrix)
rand_res <- test_lr_validity(pbmc_container,my_samp,ligand_target_matrix)

sum(cc_res<.05)/length(cc_res)
sum(scITD_res<.05)/length(scITD_res)
sum(scITD_res2<.05)/length(scITD_res2)
sum(italk_res<.05)/length(italk_res)
sum(rand_res<.05)/length(rand_res)

is_hla <- sapply(names(cc_res),function(x){
  my_h <- strsplit(x,split='-')[[1]][[1]]
  if (my_h=='HLA') {
    return(TRUE)
  } else {
    return(FALSE)
  }
})
cc_res <- cc_res[!is_hla]

is_hla <- sapply(names(scITD_res),function(x){
  my_h <- strsplit(x,split='-')[[1]][[1]]
  if (my_h=='HLA') {
    return(TRUE)
  } else {
    return(FALSE)
  }
})
scITD_res <- scITD_res[!is_hla]

## store results for plotting later
res_dat <- matrix(ncol=4,nrow=6)
colnames(res_dat) <- c('db','method','n_chan','percent_enr')
res_dat[1,] <-  c('cellchat','scITD',length(my_channels),sum(scITD_res<.05)/length(scITD_res))
res_dat[2,] <-  c('cellchat','LR_association',length(italk_channels),sum(italk_res<.05)/length(italk_res))
res_dat[3,] <-  c('italk','scITD',length(my_channels),sum(scITD_res<.05)/length(scITD_res))
res_dat[4,] <-  c('italk','LR_association',length(italk_channels),sum(italk_res<.05)/length(italk_res))
res_dat[5,] <-  c('SingleCellSignalR','scITD',length(my_channels),sum(scITD_res<.05)/length(scITD_res))
res_dat[6,] <-  c('SingleCellSignalR','LR_association',length(italk_channels),sum(italk_res<.05)/length(italk_res))

# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison.rds')
res_dat <- readRDS(file='/home/jmitchel/data/lupus_data/LR_method_comparison.rds')










## making visualization for the results
myres <- scITD_res
tmp <- cbind.data.frame(names(myres),myres)
colnames(tmp) <- c('channel','adj_pval')
tmp <- tmp[order(tmp$adj_pval,decreasing=FALSE),]
tmp$channel <- factor(tmp$channel,levels=tmp$channel)
tmp$adj_pval <- -log10(tmp$adj_pval)
mylabels <- tmp$channel
ndx_keep <- c(1,7,15,25,38,46,60,80)
mylabels2 <- rep('',length(mylabels))
mylabels2[ndx_keep] <- as.character(mylabels)[ndx_keep]
chan_bplot <- ggplot(tmp,aes(x=channel,y=adj_pval)) +
  geom_bar(stat='identity') +
  scale_x_discrete(labels= mylabels2) +
  scale_y_continuous(limits = c(0,20), expand = c(0, 0)) +
  ylab('NicheNet target gene \nenrichment -log10(padj)') +
  xlab('ligand_source_target') +
  ggtitle('scITD inferred channels') +
  geom_hline(yintercept = -log10(.05), 
             color = "red", size=.75) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -30, vjust = 0.5, hjust=0),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

pdf(file = "/home/jmitchel/figures/for_paper_v2/LR_chan_bplot2.pdf", useDingbats = FALSE,
    width = 8.5, height = 3.5)
chan_bplot
dev.off()


## make comparison bplot
# first go to 624

# change names italk to iTALK and cellchat to CellChat
res_dat[c(1,2),1] <- 'CellChat'
res_dat[c(3,4),1] <- 'iTALK'

res_dat <- as.data.frame(res_dat)
res_dat$n_chan <- as.numeric(res_dat$n_chan)
res_dat$percent_enr <- as.numeric(res_dat$percent_enr)
res_dat$db <- as.factor(res_dat$db)
res_dat$method <- as.factor(res_dat$method)

library(RColorBrewer)
mycol = brewer.pal(n = 8, name = "Dark2")

comp_bplot_size <- ggplot(res_dat,aes(x=db,y=n_chan,fill=method)) +
  geom_bar(stat='identity', position="dodge") +
  ylab('Number of inferred interactions') +
  xlab('LR database') +
  guides(fill=guide_legend(title="LR method")) +
  scale_fill_manual(values=c(mycol[8], mycol[1])) +
  theme_bw() +
  theme(legend.position = "none")

comp_bplot_size

comp_bplot_frac <- ggplot(res_dat,aes(x=db,y=percent_enr,fill=method)) +
  geom_bar(stat='identity', position="dodge") +
  ylab('Fraction of inferred interactions\nwith NicheNet enrichment') +
  xlab('LR database') +
  guides(fill=guide_legend(title="LR method")) +
  scale_fill_manual(values=c(mycol[8], mycol[1])) +
  theme_bw()

comp_bplot_frac

combined_plt <- cowplot::plot_grid(comp_bplot_size,comp_bplot_frac,nrow=1, rel_widths = c(.7,1))

pdf(file = "/home/jmitchel/figures/for_paper_v2/LR_all_db_bplot.pdf", useDingbats = FALSE,
    width = 8.5, height = 4)
combined_plt
dev.off()


## trying iTALK for LR inference
library(iTALK)
dim(iTALK::database)
head(database)

lr_pairs <- database[,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')]
colnames(lr_pairs) <- c('ligand','receptor')
lr_pairs <- unique(lr_pairs)

## finding significantly associated ligands and receptors
italk_channels <- c()
italk_pvals <- c()
ctypes <- pbmc_container$experiment_params$ctypes_use
print(nrow(lr_pairs))
for (i in 1:nrow(lr_pairs)) {
  if (i %% 100 == 0) {
    print(i)
  }
  ligand <- lr_pairs[i,1]
  receptor <- lr_pairs[i,2]
  for (source_ct in ctypes) {
    for (target_ct in ctypes) {
      if (source_ct != target_ct) {
        if (ligand %in% colnames(pbmc_container$scale_pb_extra[[source_ct]]) && receptor %in% colnames(pbmc_container$scale_pb_extra[[target_ct]])) {
          lig_exp <- pbmc_container$scale_pb_extra[[source_ct]][,ligand]
          rec_exp <- pbmc_container$scale_pb_extra[[target_ct]][,receptor]
          
          # check that expression is above certain threshold
          # if (sum(lig_exp>.2)>(.05*length(lig_exp)) && sum(rec_exp>.2)>(.05*length(rec_exp))) {
          if (sum(lig_exp>.2)>(.01*length(lig_exp)) && sum(rec_exp>.2)>(.01*length(rec_exp))) {
            tmp <- cbind.data.frame(lig_exp,rec_exp[names(lig_exp)])
            colnames(tmp) <- c('l_exp','r_exp')
            lmres <- summary(lm(l_exp~r_exp,data=tmp))
            pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
            chan <- paste0(ligand,'_',source_ct,'_',target_ct)
            italk_channels <- c(italk_channels,chan)
            italk_pvals <- c(italk_pvals,pval)
          }
          # if (sum(lig_exp!=0)>0 && sum(rec_exp!=0)>0) {
          #   tmp <- cbind.data.frame(lig_exp,rec_exp[names(lig_exp)])
          #   colnames(tmp) <- c('l_exp','r_exp')
          #   lmres <- summary(lm(l_exp~r_exp,data=tmp))
          #   pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
          #   chan <- paste0(ligand,'_',source_ct,'_',target_ct)
          #   italk_channels <- c(italk_channels,chan)
          #   italk_pvals <- c(italk_pvals,pval)
          # }
        }
      }
    }
  }
}
names(italk_pvals) <- italk_channels
italk_pvals <- p.adjust(italk_pvals,method='fdr')

# saveRDS(italk_pvals,file='/home/jmitchel/data/lupus_data/LR_italk_res.rds')

sum(italk_pvals<.05)
italk_channels <- names(italk_pvals)[italk_pvals<.05]
italk_channels <- unique(italk_channels)
length(italk_channels)







## trying the LR association method with multi receptor db like cellchat
italk_channels <- c()
italk_pvals <- c()
all_full_chan <- c()
ctypes <- pbmc_container$experiment_params$ctypes_use
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
          if (ligand %in% colnames(pbmc_container$scale_pb_extra[[source_ct]]) && receptor %in% colnames(pbmc_container$scale_pb_extra[[target_ct]])) {
            lig_exp <- pbmc_container$scale_pb_extra[[source_ct]][,ligand]
            rec_exp <- pbmc_container$scale_pb_extra[[target_ct]][,receptor]
            
            # check that expression is not 0
            if (sum(lig_exp>.2)>(.01*length(lig_exp)) && sum(rec_exp>.2)>(.01*length(rec_exp))) {
              tmp <- cbind.data.frame(lig_exp,rec_exp[names(lig_exp)])
              colnames(tmp) <- c('l_exp','r_exp')
              lmres <- summary(lm(l_exp~r_exp,data=tmp))
              pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
              chan <- paste0(ligand,'_',source_ct,'_',target_ct)
              full_chan <- paste0(ligand,'_',source_ct,'_',target_ct,'_',lr_pairs[i,2])
              italk_channels <- c(italk_channels,chan)
              italk_pvals <- c(italk_pvals,pval)
              all_full_chan <- c(all_full_chan,full_chan)
            }
          }
        }
      }
    }
  }
}

# adjust pvals
italk_pvals <- p.adjust(italk_pvals,method='fdr')

## get channels where there was significant associations with any receptor components
italk_channels2 <- italk_channels[italk_pvals<.05]
italk_channels <- italk_channels2
italk_channels <- unique(italk_channels)
length(italk_channels)

## for parsing full chan results
new_italk_channels <- c()
all_full_chan_unq <- unique(all_full_chan)
for (i in 1:length(all_full_chan_unq)) {
  afq <- all_full_chan_unq[i]
  ndx_test <- which(all_full_chan==afq)
  pv_test <- italk_pvals[ndx_test]
  # check all pv less than thresh
  if (sum(pv_test<.05)==length(pv_test)) {
    new_italk_channels <- c(new_italk_channels, italk_channels[ndx_test[1]])
  }
}
length(new_italk_channels)
new_italk_channels <- unique(new_italk_channels)
length(new_italk_channels)


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
}




