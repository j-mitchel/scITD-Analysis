



## preprocessing my scITD output to use
myres_mat <- pbmc_container$lr_res #includes NS results still but are FDR adjusted

# reduce to rows/columns with at least one significant hit
myres_mat <- myres_mat[rowSums(myres_mat<.0001)>0,]
myres_mat <- myres_mat[,colSums(myres_mat<.0001)>0]

## get unique ligand_source_target channels from my analysis
unique_channels <- c()
for (i in 1:nrow(myres_mat)) {
  lig_ct_rec <- myres_mat[i,]
  lig_ct_rec_name <- strsplit(rownames(myres_mat)[i],split='_')[[1]]
  lig <- lig_ct_rec_name[[1]]
  source <- lig_ct_rec_name[[2]]
  
  for (j in 1:ncol(myres_mat)) {
    pv <- lig_ct_rec[j]
    if (pv < 0.0001) {
      target_ct <- strsplit(names(pv),split="_")[[1]][[1]]
      lig_source_target <- paste0(lig,"_",source,"_",target_ct)
      unique_channels <- c(unique_channels,lig_source_target)
    }
  }
}
unique_channels <- unique(unique_channels)
my_channels <- unique_channels
# saveRDS(my_channels,file='/home/jmitchel/data/lupus_data/LR_scITD_res.rds')
# saveRDS(my_channels,file='/home/jmitchel/data/lupus_data/LR_iTALK_scITD_res.rds')
# saveRDS(my_channels,file='/home/jmitchel/data/lupus_data/LR_singlecellsig_scITD_res.rds')












## for getting channels by the LR association method
# need to have run LR preprocessing with appropriate lr_pairs db first!
get_LR_association_chan <- function(container,lr_pairs) {
  my_channels <- c()
  all_pvals <- c()
  ctypes <- container$experiment_params$ctypes_use
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
            if (ligand %in% colnames(container$scale_pb_extra[[source_ct]]) && receptor %in% colnames(container$scale_pb_extra[[target_ct]])) {
              lig_exp <- container$scale_pb_extra[[source_ct]][,ligand]
              rec_exp <- container$scale_pb_extra[[target_ct]][,receptor]
              
              # check that expression is not 0
              if (sum(lig_exp>.2)>(.01*length(lig_exp)) && sum(rec_exp>.2)>(.01*length(rec_exp))) {
                tmp <- cbind.data.frame(lig_exp,rec_exp[names(lig_exp)])
                colnames(tmp) <- c('l_exp','r_exp')
                lmres <- summary(lm(l_exp~r_exp,data=tmp))
                pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
                chan <- paste0(ligand,'_',source_ct,'_',target_ct)
                my_channels <- c(my_channels,chan)
                all_pvals <- c(all_pvals,pval)
              }
            }
          }
        }
      }
    }
  }
  
  # adjust pvals
  all_pvals <- p.adjust(all_pvals,method='fdr')
  
  ## get channels where there was significant associations with any receptor components
  my_channels <- my_channels[all_pvals<.05]
  my_channels <- unique(my_channels)
  return(my_channels)
}









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
  return(my_samp)
}














##### generating data to use with cellphonedb
library(Seurat)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

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


# subset data to SLE patients only
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$Status=='Managed']
pbmc <- subset(pbmc,cells = cells_keep)

# reduce to cells of the cell types that I'm using
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$cg_cov %in% c("B","NK","Th","Tc","cDC","cMono","ncMono")]

cells_keep <- sample(cells_keep,50000)
pbmc <- subset(pbmc,cells = cells_keep)
pbmc@meta.data$cg_cov <- factor(pbmc@meta.data$cg_cov,levels=c("B","NK","Th","Tc","cDC","cMono","ncMono"))
pbmc@meta.data$idents <- pbmc@meta.data$cg_cov
Idents(pbmc) <- pbmc@meta.data$idents
pbmc@meta.data$cell_type <- pbmc@meta.data$cg_cov
write.csv(as.matrix(pbmc@assays$RNA@data), file = '/home/jmitchel/data/LR_datasets/cellphone/lupus_subsetted_cpdb_counts.csv')
write.csv(pbmc@meta.data, file = '/home/jmitchel/data/LR_datasets/cellphone/lupus_subsetted_cpdb_meta.csv')

## read in cpdb results and parse them


