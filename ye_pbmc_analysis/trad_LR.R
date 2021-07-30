

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
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$cg_cov %in% c("B","NK","T4","T8","cDC","cM","ncM")]

cells_keep <- sample(cells_keep,50000)
pbmc <- subset(pbmc,cells = cells_keep)
pbmc@meta.data$cg_cov <- factor(pbmc@meta.data$cg_cov,levels=c("B","NK","T4","T8","cDC","cM","ncM"))

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


































