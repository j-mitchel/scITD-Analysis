library(Seurat)
library(ggplot2)
library(devtools)
load_all('/home/jmitchel/scITD/')

# load up the lupus dataset: see preprocessing/lupus_preprocessing.R 
# for code used to generate this object
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

param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc","cDC",
                                               "cMono","ncMono"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(seurat_obj=pbmc,
                                     params=param_list,
                                     metadata_cols=c('ind_cov_batch_cov',
                                                     "SLE_status",
                                                     "Status",
                                                     "cg_cov",
                                                     "sex",
                                                     "Age",
                                                     "batch_cov",
                                                     "Processing_Cohort",
                                                     "Ethnicity"),
                                     metadata_col_nm=c('donors',
                                                       'SLE_status',
                                                       'Status',
                                                       'ctypes',
                                                       'sex',
                                                       'Age',
                                                       'pool',
                                                       'processing',
                                                       'Ethnicity'))


pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')



# extract pseudobulk data for a couple cell types
pb1 <- pbmc_container$scMinimal_ctype$cMono$pseudobulk
pb2 <- pbmc_container$scMinimal_ctype$Th$pseudobulk

pcs_ct1 <- prcomp(pb1)$x[,1:5]
pcs_ct2 <- prcomp(pb2)$x[,1:5]

cormat <- abs(cor(pcs_ct1,pcs_ct2[rownames(pcs_ct1),]))
colnames(cormat) <- paste0('cMono_PC_',1:ncol(cormat))
rownames(cormat) <- paste0('Th_PC_',1:nrow(cormat))

col_fun = colorRamp2(c(0, 1), c("white", "red"))
hmap <- Heatmap(cormat,name = "pearson r",
                cluster_columns = FALSE,
                col = col_fun,
                cluster_rows = FALSE,
                show_row_names = TRUE,
                show_column_names = TRUE,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                column_names_side = 'top',
                row_names_side = 'left',
                column_names_rot = 30,
                column_names_gp = grid::gpar(fontsize = 8),
                row_names_gp = grid::gpar(fontsize = 8),
                border = TRUE,
                column_title = 'PCA on cell types separately\nPC correlations',
                column_title_gp = gpar(fontsize = 14),
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid::grid.text(sprintf("%.2f", cormat[i, j]), x, y, gp = gpar(fontsize = 10))
                })


# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/two_ct_pca_cors.pdf", useDingbats = FALSE,
#     width = 4.5, height = 3.5)
hmap
dev.off()









