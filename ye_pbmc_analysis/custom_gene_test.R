

# load up the subsetted dataset
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



##### running edgeR DE analysis for LN covariate
library(edgeR)
# need to get pseudobulk counts to start
pbmc_container <- parse_data_by_ctypes(pbmc_container)
pbmc_container <- clean_data(pbmc_container, donor_min_cells=20)
pbmc_container <- get_pseudobulk(pbmc_container)

# using Th counts only to start off with
pb <- pbmc_container[["scMinimal_ctype"]][["Th"]][["pseudobulk"]]

# getting clinical annotations for LN
library(readxl)
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_categorical.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL

## get donors in both dsc and in clin_vars
# trim donor IDs in dsc
trim_names <- sapply(colnames(pb), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})

# subset clin_vars to be same patients
clin_vars <- clin_vars[trim_names,]
rownames(clin_vars) <- names(trim_names)

LN <- as.factor(clin_vars[,'crflupusneph']) # should be factor of 1s and 0s
names(LN) <- rownames(clin_vars)
y <- DGEList(counts=pb, genes=rownames(pb))
isexpr <- filterByExpr(y, group=LN)
table(isexpr)

hasannot <- rowSums(is.na(y$genes))==0
y <- y[isexpr & hasannot, , keep.lib.sizes=FALSE]
dim(y)
y <- calcNormFactors(y)
head(y$samples)

## getting batch variable to regress out
pbmc_container <- get_donor_meta(pbmc_container,additional_meta = 'pool')
head(pbmc_container$donor_metadata)
mybatch <- pbmc_container$donor_metadata[names(LN),'pool']
names(mybatch) <- names(LN)

# design <- model.matrix(~LN)
design <- model.matrix(~mybatch+LN)
rownames(design) <- colnames(y)
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)
fit <- glmQLFit(y, design, robust=TRUE)
plotQLDisp(fit)
qlf <- glmQLFTest(fit)
top_de <- rownames(topTags(qlf,n=1000))
#####

## genes to use in tensor are a combo of top genes and random selected genes
mygenes_rand <- sample(all_genes,1000)
mygenes <- unique(c(mygenes_rand,top_de))

# add LN factor to meta_data
pbmc@meta.data$LN <- sapply(pbmc@meta.data$ind_cov_batch_cov,function(x) {
  return(LN[as.character(x)])
})


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
                                                     "Ethnicity",
                                                     "LN"),
                                     metadata_col_nm=c('donors',
                                                       'SLE_status',
                                                       'Status',
                                                       'ctypes',
                                                       'sex',
                                                       'Age',
                                                       'pool',
                                                       'processing',
                                                       'Ethnicity',
                                                       "LN"))

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool', custom_genes = mygenes)


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(10,30),
                                 tucker_type = 'regular', rotation_type = 'hybrid')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(18,50),
                                 tucker_type = 'regular', rotation_type = 'hybrid')


pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Ethnicity','LN'),
                                        stat_use='pval')

## plot donor score
pbmc_container <- plot_donor_matrix(pbmc_container,meta_vars = c('LN','Ethnicity'),
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

pbmc_container$plots$donor_matrix
