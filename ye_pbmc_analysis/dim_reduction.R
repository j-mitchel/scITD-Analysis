library(scITD)
library(Seurat)

# load up the subsetted dataset
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v2.rds')

# set up project parameters
param_list <- initialize_params(ctypes_use = c("B","NK","T4","T8","cDC",
                                               "cM","ncM"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(seurat_obj=pbmc,
                                     params=param_list,
                                     metadata_cols=c('ind_cov_batch_cov',
                                                     "SLE_status",
                                                     "Status",
                                                     "cg_cov",
                                                     "sex",
                                                     "age",
                                                     "batch_cov",
                                                     "Processing_Cohort"),
                                     metadata_col_nm=c('donors',
                                                       'SLE_status',
                                                       'Status',
                                                       'ctypes',
                                                       'sex',
                                                       'age',
                                                       'pool',
                                                       'processing'))

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=10, gene_min_cells=10,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var', vargenes_thresh=500,
                              scale_var = TRUE, var_scale_power = 1.5,
                              batch_var='pool')




## doing dim reduction outside of function for fear it will crash after 2 hours
integration_var <- 'pool'
container <- pbmc_container
rm(pbmc_container)
ncores <- container$experiment_params$ncores

# some cells have been removed because donors had too few cells per ctype
# need to make sure the full data is limited to the cells used in analysis
all_cells <- c()
for (ct in container$experiment_params$ctypes_use) {
  cells_in_ctype <- rownames(container$scMinimal_ctype[[ct]]$metadata)
  all_cells <- c(all_cells,cells_in_ctype)
}

container$scMinimal_full$metadata <- container$scMinimal_full$metadata[all_cells,]
container$scMinimal_full$count_data <- container$scMinimal_full$count_data[,all_cells]

# create a list of subsetted data matrices (one per var value)
panel <- list()
meta <- as.character(container$scMinimal_full$metadata[,integration_var])
var_vals <- unique(meta)
for (v in var_vals) {
  cell_ndx <- which(meta == v)
  panel[[v]] <- container$scMinimal_full$count_data[,cell_ndx]
}

# turn the list of matrices to list of pagoda2 objects
panel.preprocessed <- lapply(panel, pagoda2::basicP2proc, n.cores=ncores,
                             min.cells.per.gene=0, n.odgenes=2e3,
                             get.largevis=FALSE, make.geneknn=FALSE)

saveRDS(panel.preprocessed,file='/home/jmitchel/data/lupus_data/lupus_pagoda.rds')

con <- conos::Conos$new(panel.preprocessed, n.cores=ncores)

# build graph
con$buildGraph()

saveRDS(con,file='/home/jmitchel/data/lupus_data/lupus_conos.rds')

# make umap embedding
con$embedGraph(method="UMAP", min.dist=0.01, spread=15, n.cores=ncores, min.prob.lower=1e-3)

# assign ctype names to the cells
con$findCommunities(method=conos::leiden.community, resolution=1)
cell_assigns <- container$scMinimal_full$metadata[,"ctypes"]
names(cell_assigns) <- rownames(container$scMinimal_full$metadata)
con$clusters$leiden$groups <- cell_assigns[names(con$clusters$leiden$groups)]

saveRDS(con,file='/home/jmitchel/data/lupus_data/lupus_conos.rds')