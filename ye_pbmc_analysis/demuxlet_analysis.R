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

pbmc_container <- make_new_container(count_data=pbmc_all, meta_data=pbmc_meta,
                                     params=param_list,
                                     label_donor_sex = TRUE)

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=10, gene_min_cells=10,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                              scale_var = TRUE, var_scale_power = 1.5)

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(3,6,6),
                                 tucker_type = 'regular', rotation_type = 'ica')

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('stim'),stat_use='pval')

pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('stim'),
                                    cluster_by_meta = 'stim',
                                    show_donor_ids = TRUE,
                                    add_meta_associations='pval')
pbmc_container$plots$donor_matrix



