library(Seurat)

# load up the subsetted dataset
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')

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


# form tensor with batch correction applied
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                              scale_var = TRUE, var_scale_power = 1.5,
                              batch_var='pool')
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                              scale_var = TRUE, var_scale_power = 0,
                              batch_var='pool')

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,20,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,16,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,24,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(14,24,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(16,40,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(3,8,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,16,7),
                                 tucker_type = 'sparse', rotation_type = 'ica',
                                 sparsity=2)
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,16,7),
                                 tucker_type = 'sparse', rotation_type = 'ica',
                                 sparsity=sqrt(2))
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(10,20,7),
                                 tucker_type = 'sparse', rotation_type = 'ica')

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Status','Ethnicity'),
                                        stat_use='pval')

# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    cluster_by_meta = 'sex',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# pdf(file = "/home/jmitchel/figures/for_paper/lupus_dscores_v2.pdf", useDingbats = FALSE,
#     width = 6, height = 7)
pbmc_container$plots$donor_matrix
dev.off()

# 
# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(10,20,7), n_fibers=100, n_iter=1000,
                                tucker_type='regular', rotation_type='ica')

# saveRDS(pbmc_container[["gene_score_associations"]],file='/home/jmitchel/data/lupus_data/lupus_jackstraw_lds_ica.rds')
pbmc_container[["gene_score_associations"]] <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_jackstraw_lds_ica.rds')


# get loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=0.01,
                                           display_genes=FALSE,
                                           gene_callouts=TRUE,
                                           callout_n_gene_per_ctype=5)

# render loadings figure
myfig <- render_multi_plots(pbmc_container,data_type='loadings',max_cols=3)
myfig

pdf(file = "/home/jmitchel/figures/for_paper/lds_ica_lds.pdf", useDingbats = FALSE,
    width = 18, height = 23)
myfig
dev.off()





# plotting IFN response gene expression against f1
myfactor <- 3
mygene <- 'IFI6'
mygene <- 'IFI44L'
# mygene <- 'G0S2'
# mygene <- 'ZNF331'
mygene <- 'CD69'
d_exp <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,mygene]
# d_exp <- pbmc_container[["scMinimal_ctype"]][['NK']][["pseudobulk"]][,mygene]
dsc <- pbmc_container$tucker_results[[1]][,myfactor]
tmp <- cbind.data.frame(dsc[names(d_exp)],d_exp)
colnames(tmp) <- c('dscore','expres')

# add regression line
lmres <- lm(expres~dscore,data=tmp)
line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

ggplot(tmp,aes(x=dscore,y=expres)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  ylab(paste0(mygene,' expression (CD4+ T)')) +
  xlab(paste0('Factor ',myfactor,' donor scores')) +
  theme_bw()
dev.off()


pbmc_container <- plot_donor_sig_genes(pbmc_container, factor_select=2,
                                       top_n_per_ctype=10, show_donor_labels=F)

pbmc_container[["plots"]][["donor_sig_genes"]][["2"]]

pbmc_container <- get_lm_pvals(pbmc_container)
myfactor <- 6
f1_data <- get_one_factor(pbmc_container, factor_select=myfactor)
dscores <- f1_data[[1]]
lds <- f1_data[[2]]

lds['ISG15',]
lds['IFI6',]
lds['IFIT1',]
lds['IFI27',]
lds['IFITM3',]
lds['IFITM3',]
lds['G0S2',]
lds['ZNF331',]
lds['GIMAP7',]

sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

# limit to just the genes in tmp_casted_num
sig_df <- sig_df[rownames(lds),colnames(lds)]
# sum(sig_df<0.05)
sig_df['ISG15',]
sig_df['IFI6',]
sig_df['IFIT1',]
sig_df['IFI27',]
sig_df['IFITM3',]
sig_df['G0S2',]
sig_df['ZNF331',]

sig_df['CCL5',]
sig_df['EGR1',]
sig_df['GBP7',]
sig_df['IFI30',]
sig_df['IFNG',]
sig_df['IRF8',]
sig_df['PTAFR',]
sig_df['SETD2',]
sig_df['TLR4',]
sig_df['XCL2',]

## compute factor associations for IFI6 manually since it takes long time tcko do jackstraw
container <- pbmc_container
tensor_data <- container$tensor_data
tucker_results <- container$tucker_results

n_genes <- length(tensor_data[[2]])
n_ctypes <- length(tensor_data[[3]])
fstats_real <- data.frame(matrix(ncol=3,nrow=0))
ncores <- 20
fstats_real <- mclapply(1:n_genes, function(i) {
  gene <- tensor_data[[2]][i]
  gene_res <- list()
  for (j in 1:n_ctypes) {
    ctype <- tensor_data[[3]][j]
    gene_res[[ctype]] <- list()
    for (k in 1:ncol(tucker_results[[1]])) {
      tmp_fiber <- tensor_data[[4]][,i,j]
      df_test <- as.data.frame(cbind(tmp_fiber, tucker_results[[1]][,k]))
      colnames(df_test) <- c('fiber','factor')
      lmres <- lm(factor~fiber,df_test)

      # testing getting pval instead of fstat
      x <- summary(lmres)
      fstat <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)
      
      gene_res[[ctype]][[as.character(k)]] <- fstat
    }
  }
  return(gene_res)
}, mc.cores = ncores)

names(fstats_real) <- tensor_data[[2]]

# unpack the list
fstats_real <- unlist(fstats_real)

fstats_real <- p.adjust(fstats_real,method='fdr')

# fstats_real['IFI6.T4.1.value']

new_names <- sapply(names(fstats_real),function(x) {
  tmp <- strsplit(x,split = '.', fixed = TRUE)[[1]]
  if (length(tmp)==4) {
    return(paste0(tmp[[1]],'.',tmp[[2]],'.',tmp[[3]]))
  } else if (length(tmp)==5) {
    return(paste0(tmp[[1]],'.',tmp[[2]],'.',tmp[[3]],'.',tmp[[4]]))
  } else if (length(tmp)==6) {
    return(paste0(tmp[[1]],'.',tmp[[2]],'.',tmp[[3]],'.',tmp[[4]],'.',tmp[[5]]))
  }
})
names(new_names) <- NULL
names(fstats_real) <- new_names
pbmc_container[["gene_score_associations"]] <- fstats_real
identical(container,pbmc_container)


# get loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=0.01,
                                           display_genes=FALSE,
                                           gene_callouts=TRUE,
                                           callout_n_gene_per_ctype=5)

# render loadings figure
myfig <- render_multi_plots(pbmc_container,data_type='loadings',max_cols=3)

# pdf(file = "/home/jmitchel/figures/for_paper/lds_varimax_lds.pdf", useDingbats = FALSE,
#     width = 18, height = 23)
myfig
dev.off()

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, sig_thresh=0.01, display_genes=FALSE,
                                      gene_callouts=TRUE, callout_n_gene_per_ctype=5)


pbmc_container <- get_lm_pvals(pbmc_container)

myfactor <- 5
f1_data <- get_one_factor(pbmc_container, factor_select=myfactor)
dscores <- f1_data[[1]]
lds <- f1_data[[2]]
lds['ISG15',]
lds['IFI6',]
lds['IFITM3',]
lds['IFIT3',]
lds['IFI27',]
lds['ISG20',]
lds['CD69',]

sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

# limit to just the genes in tmp_casted_num
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df['ISG15',]
sig_df['IFI6',]
sig_df['IFITM3',]
sig_df['IFIT3',]
sig_df['IFI27',]
sig_df['ISG20',]
sig_df['SOCS1',]
sig_df['CX3CR1',]
sig_df['LGALS1',]
sig_df['CD69',]


## trying tucker version from package CA3variants: I think this uses tuckals3 but not sure
tucker_decomp <- CA3variants::tucker(tnsr,10,20,7)


# trying an ica rotation of just the core
core_un <- rTensor::k_unfold(tucker_decomp$Z,1)@data
ica_res <- ica::icafast(t(core_un),ranks[1],center=FALSE,alg='def')
ldngs <- t(ica_res$S) %*% t(kron_prod)
donor_mat <- donor_mat %*% solve(ica_res$W)

# or with varimax
core_un <- rTensor::k_unfold(tucker_decomp$Z,1)@data
ica_res <- varimax(t(core_un))
core_un <- t(t(core_un) %*% ica_res$rotmat)
ldngs <- core_un %*% t(kron_prod)
donor_mat <- donor_mat %*% solve(t(ica_res$rotmat))



# trying version of tucker from package rrcov3way
library(rrcov3way)
tucker_decomp <- Tucker3(tnsr,10,20,7,robust=TRUE)
gene_by_factors <- tucker_decomp$B
rownames(gene_by_factors) <- gene_nm
ctype_by_factors <- tucker_decomp$C
rownames(ctype_by_factors) <- ctype_nm
donor_mat <- tucker_decomp$A
rownames(donor_mat) <- donor_nm

kron_prod <- kronecker(ctype_by_factors,gene_by_factors,make.dimnames = TRUE)

ldngs <- tucker_decomp$GA %*% t(kron_prod)

# lets see if the core is actually orthogonal
dim(tucker_decomp$GA)
core_tmp <- t(tucker_decomp$GA)
t(core_tmp[,1]) %*% core_tmp[,3]
# doesnt seem to be

# trying ica rotation of loadings with it
ica_res <- ica::icafast(t(ldngs),ranks[1],center=FALSE,alg='def')
ica_res <- ica::icafast(t(ldngs),ranks[1],center=FALSE,alg='par')
# ica_res <- ica::icafast(t(ldngs),ranks[1],center=FALSE,alg='def',maxit = 1000,tol = 1e-15)
ldngs <- t(ica_res$S)
donor_mat <- donor_mat %*% solve(ica_res$W)
# seeing if rotation matrix is orthonormal
t(solve(ica_res$W)[,1]) %*% solve(ica_res$W)[,2]
sqrt(sum(donor_mat[,1]**2))
# nope not orthonormal again...

# I almost wonder if there is a bug in this ICA code
# like it seems as if they shouldn't be transposing the rotation matrix
# because the rotation matrix is orthogonal when not transposed
# could try a different ICA package to see if this "fixes" the issue




# trying the decomposition with fastICA instead of icafast
library(fastICA)
ica_res <- fastICA(t(ldngs),ranks[1])
ica_res$S
mysource <- ica_res$X %*% ica_res$K %*% ica_res$W
all.equal(ica_res$S,mysource)

rotmat <- ica_res$K %*% ica_res$W
# ldngs <- t(ica_res$S)
ldngs <- t(rotmat) %*% ldngs
donor_mat <- donor_mat %*% solve(t(rotmat))






## checking to see if loadings matrices seem sparse without reducing them to sig genes only
pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=6, use_sig_only=F, nonsig_to_zero=F)


# get loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=F,
                                           nonsig_to_zero=F)

# render loadings figure
myfig <- render_multi_plots(pbmc_container,data_type='loadings',max_cols=3)
myfig



## trying to see if "wrong" loadings signs are because it is accounting for only some residuals

tnsr <- rTensor::as.tensor(container$tensor_data[[4]])
donor_mat <- container$tucker_results[[1]]
ldngs <- container$tucker_results[[2]]

factor_use <- c(1,4)
gene1 <- 'T4:IFI6'
gene2 <- 'T4:ISG15'
gene1 <- 'T4:IFI27'
gene1 <- 'T4:IFI44L'

factor_use <- c(1)
gene1 <- 'T4:IFI6'
gene1 <- 'T8:IFI6'
gene2 <- 'T4:IFI27'
gene1 <- 'cDC:IFI27'
gene1 <- 'T4:CD69'

recon <- donor_mat[,factor_use,drop=FALSE] %*% ldngs[factor_use,,drop=FALSE]
# recon_ndx1 <- which(colnames(recon)==gene1) %% tnsr@modes[2]
# recon_ndx2 <- which(colnames(recon)==gene2) %% tnsr@modes[2]
r1 <- recon[,gene1]
r2 <- recon[,gene2]

# recon_tnsr <- rTensor::k_fold(recon,m=1,modes=tnsr@modes)

# test1 <- as.matrix(recon[1:5,(2*tnsr@modes[2]+recon_ndx):(2*tnsr@modes[2]+recon_ndx+2)])
# test2 <- as.matrix(recon_tnsr@data[1:5,recon_ndx:(recon_ndx+2),3])
# all.equal(test1,test2,check.attributes=FALSE)

# tmp <- tnsr[,recon_ndx,3]@data
# names(tmp) <- container$tensor_data[[1]]

myf <- get_one_factor(container,factor_select=3)
# dsc <- myf[[1]]
dsc <- container$tucker_results[[1]][,3]
lds <- myf[[2]]

lds['IFI6','T4']
lds['IFI6','T8']
lds['IFI44L','T4']
lds['IFI27','cDC']
lds['CD69','T4']

d_exp1 <- pbmc_container[["scMinimal_ctype"]][['T8']][["pseudobulk"]][,'IFI6']
d_exp2 <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,'ISG15']
d_exp1 <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,'IFI27']
d_exp1 <- pbmc_container[["scMinimal_ctype"]][['cDC']][["pseudobulk"]][,'IFI27']
d_exp1 <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,'IFI44L']
d_exp1 <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,'CD69']

d_exp1 <- d_exp1[names(r1)]
d_exp2 <- d_exp2[names(r2)]

resid1 <- d_exp1[names(r1)] - r1
# plot(resid1,dsc[names(resid1)])
plot(dsc[names(resid1)],resid1)
plot(dsc_saved[names(resid1)],resid1)

resid2 <- d_exp2[names(r2)] - r2
# plot(resid2,dsc[names(resid2)])
plot(dsc[names(resid2)],resid2)

## seeing if the pred donors with positive residuals have low F1 scores
pred_don <- names(resid2)[dsc[names(resid2)] < -.15]
resid2[pred_don]

dsc2 <- container$tucker_results[[1]][,1]
dsc2[pred_don]
dsc2[order(dsc2,decreasing=T)][1:10]

# need to explicity look at the pred donors WITH known IFN signaling, then 
# see their F1 scores
myfactor <- 4
mygene <- 'IFI27'
d_exp <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]][,mygene]
dsc <- pbmc_container$tucker_results[[1]][,myfactor]
tmp <- cbind.data.frame(dsc[names(d_exp)],d_exp)
colnames(tmp) <- c('dscore','expres')

# add regression line
lmres <- lm(expres~dscore,data=tmp)
line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

ggplot(tmp,aes(x=dscore,y=expres)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  ylab(paste0(mygene,' expression (CD4+ T)')) +
  xlab(paste0('Factor ',myfactor,' donor scores')) +
  theme_bw()

tmp_pred <- tmp[tmp$dscore<(-.1),]
tmp_pred_ifn <- tmp_pred[tmp_pred$expres>4,]
head(tmp_pred_ifn)

# see if the tmp_pred_ifn donors have large or medium size residuals
plot(dsc[rownames(tmp_pred_ifn)],resid2[rownames(tmp_pred_ifn)])
# yup they basically have some amount of under-prediction for IFI27
# I wonder if degree of under-prediction is correlated with F1 score
dsc2 <- container$tucker_results[[1]][,1]
dsc2[rownames(tmp_pred_ifn)]
resid2[rownames(tmp_pred_ifn)]



# seeing if the wrong loading signs






myf <- get_one_factor(container,factor_select=4)
dsc <- myf[[1]]
lds <- myf[[2]]

dsc_ica <- dsc
dsc_vari <- dsc

# get IFN genes from vari F1
pbmc_container <- get_lm_pvals(pbmc_container)
myfactor <- 1
f1_data <- get_one_factor(pbmc_container, factor_select=myfactor)
dscores <- f1_data[[1]]
lds <- f1_data[[2]]
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
IFN_genes <- rownames(sig_df)[rowSums(sig_df<0.01)>0]
IFN_genes2 <- rownames(sig_df)[rowSums(sig_df<0.000000000000000000001)>0]
IFN_genes[1:5]

# use GO IFN gene set
go_ifn1 <- read.csv(file='/home/jmitchel/IFN_gene_list.csv')
go_ifn2 <- read.csv(file='/home/jmitchel/IFN_gamma_gene_list.csv')
go_ifn_all <- unique(c(go_ifn1[,1],go_ifn2[,1]))
sum(IFN_genes %in% go_ifn_all)
sum(rownames(sig_df) %in% go_ifn_all)
ifn_pres <- rownames(sig_df)[rownames(sig_df) %in% go_ifn_all]
f1_ifn <- IFN_genes[IFN_genes %in% go_ifn_all]
f_oth_ifn <- ifn_pres[!(ifn_pres %in% f1_ifn)]

# now removing IFN genes from the tensor...
ndx_keep <- which(!(pbmc_container[["tensor_data"]][[2]] %in% IFN_genes))
ndx_keep <- which(!(pbmc_container[["tensor_data"]][[2]] %in% f1_ifn))
ndx_keep <- which(!(pbmc_container[["tensor_data"]][[2]] %in% IFN_genes2))
pbmc_container[["tensor_data"]][[3]][1]
test1 <- pbmc_container[["scMinimal_ctype"]][["B"]][["pseudobulk"]][,ndx_keep]
test2 <- pbmc_container[["tensor_data"]][[4]][,ndx_keep,1]
test1[1:5,1:5]
test2[1:5,1:5]

# old_tensor_dat <- pbmc_container$tensor_data
pbmc_container$tensor_data <- old_tensor_dat
pbmc_container[["tensor_data"]][[4]] <- pbmc_container[["tensor_data"]][[4]][,ndx_keep,]
pbmc_container[["tensor_data"]][[2]] <- pbmc_container[["tensor_data"]][[2]][ndx_keep]

# now rerunning decomposition using ICA method
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(8,22,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,16,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(6,16,7),
                                 tucker_type = 'regular', rotation_type = 'ica')
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,14,7),
                                 tucker_type = 'regular', rotation_type = 'ica')

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Status'),
                                        stat_use='pval')

pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    cluster_by_meta = 'sex',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

pbmc_container$plots$donor_matrix


# calculating whether the donors scores are still associated with IFN genes after removal
myf <- get_one_factor(pbmc_container,factor_select=2)
dsc <- myf[[1]]
lds <- myf[[2]]

pb <- pbmc_container[["scMinimal_ctype"]][['T4']][["pseudobulk"]]
all_pv <- c()
# for (gn in core_set) {
for (gn in IFN_genes) {
  tmp <- cbind.data.frame(dsc,pb[rownames(dsc),gn])
  # tmp <- cbind.data.frame(dsc_pred_no_IFN,pb[rownames(dsc_pred_no_IFN),gn])
  # tmp <- cbind.data.frame(dsc_t4_no_IFN,pb[rownames(dsc_t4_no_IFN),gn])
  colnames(tmp) <- c('dscores','expr')
  lmres <- lm(expr~dscores,data=tmp)
  x <- summary(lmres)
  pval <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)
  all_pv <- c(all_pv,pval)
}
# names(all_pv) <- core_set
names(all_pv) <- IFN_genes
all_pv2 <- p.adjust(all_pv,method='fdr')
all_pv2[order(all_pv2)][1:100]
all_pv2[core_set]
all_pv2['ISG20']
all_pv2['CX3CR1']

# do hypergeometric test among IFN_genes and associated to show no IFN enrichment

# now extracting the new prednisone factor
myf <- get_one_factor(container,factor_select=4)
dsc <- myf[[1]]
lds <- myf[[2]]

identical(rownames(dsc_vari),rownames(dsc))
identical(rownames(dsc_ica),rownames(dsc))

plot(dsc[,1],dsc_ica[,1])
plot(dsc[,1],dsc_vari[,1])
plot(dsc_ica[,1],dsc_vari[,1])
cor(dsc[,1],dsc_ica[,1])
cor(dsc[,1],dsc_vari[,1])


## now using var scale power of .5 and saving the decomps that include IFN genes
myf <- get_one_factor(pbmc_container,factor_select=4)
dsc_vari_4 <- myf[[1]]
lds_vari_4 <- myf[[2]]

myf <- get_one_factor(pbmc_container,factor_select=5)
dsc_vari_5 <- myf[[1]]
lds_vari_5 <- myf[[2]]

myf <- get_one_factor(pbmc_container,factor_select=4)
dsc_ica_4 <- myf[[1]]
lds_ica_4 <- myf[[2]]

myf <- get_one_factor(pbmc_container,factor_select=6)
dsc_ica_6 <- myf[[1]]
lds_ica_6 <- myf[[2]]

identical(rownames(dsc_vari),rownames(dsc_ica))
plot(dsc_vari[,1],dsc_ica[,1])
identical(rownames(lds_vari),rownames(lds_ica))
plot(lds_vari[,'T4'],lds_ica[,'T4'])

myf <- get_one_factor(pbmc_container,factor_select=4)
dsc_ica_no_ifn <- myf[[1]]
lds_ica_no_ifn <- myf[[2]]

myf <- get_one_factor(pbmc_container,factor_select=3)
dsc_ica_no_ifn_3 <- myf[[1]]
lds_ica_no_ifn_3 <- myf[[2]]

myf <- get_one_factor(pbmc_container,factor_select=7)
dsc_ica_no_ifn_7 <- myf[[1]]
lds_ica_no_ifn_7 <- myf[[2]]

plot(dsc_ica_no_ifn_3[,1],dsc_ica_4[,1])
plot(dsc_ica_no_ifn_3[,1],dsc_vari_4[,1])
cor(dsc_ica_no_ifn_3[,1],dsc_ica_4[,1])
cor(dsc_ica_no_ifn_3[,1],dsc_vari_4[,1])

plot(dsc_ica_no_ifn_7[,1],-dsc_ica_6[,1])
plot(dsc_ica_no_ifn_7[,1],dsc_vari_5[,1])
cor(dsc_ica_no_ifn_7[,1],dsc_ica_6[,1])
cor(dsc_ica_no_ifn_7[,1],dsc_vari_5[,1])

plot(dsc_ica_no_ifn[,1],-dsc_ica[,1])
plot(dsc_ica_no_ifn[,1],dsc_vari[,1])

plot(lds_ica_no_ifn[,1],lds_vari[rownames(lds_ica_no_ifn),1])
plot(lds_ica_no_ifn[,'T4'],lds_ica[rownames(lds_ica_no_ifn),'T4'])

cor(dsc_ica_no_ifn[,1],dsc_ica[,1])
cor(dsc_ica_no_ifn[,1],dsc_vari[,1])








## I want to see cases where ICA correctly predicts an association but varimax misses it
# x-axis can be ICA significance values and y-axis varimax significance values
pbmc_container <- get_lm_pvals(pbmc_container)

fstats_real_ica <- pbmc_container[["gene_score_associations"]]
fstats_real_vari <- pbmc_container[["gene_score_associations"]]
fstats_real_combo <- pbmc_container[["gene_score_associations"]]
pbmc_container[["gene_score_associations"]] <- fstats_real_ica
pbmc_container[["gene_score_associations"]] <- fstats_real_vari
pbmc_container[["gene_score_associations"]] <- fstats_real_combo

myfactor <- 6
myf <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- myf[[1]]
lds <- myf[[2]]
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df_ica <- sig_df
sig_df_vari <- sig_df
sig_df_combo <- sig_df

identical(rownames(sig_df_vari),rownames(sig_df_ica))
plot(-log10(sig_df_vari[,2]),-log10(sig_df_ica[,2]))

sig_df_ica_IFN <- sig_df_ica[rownames(sig_df_ica) %in% IFN_genes,]
sig_df_vari_IFN <- sig_df_vari[rownames(sig_df_vari) %in% IFN_genes,]
sig_df_combo_IFN <- sig_df_combo[rownames(sig_df_combo) %in% IFN_genes,]

sig_df_ica_IFN <- sig_df_ica[rownames(sig_df_ica) %in% IFN_genes2,]
sig_df_vari_IFN <- sig_df_vari[rownames(sig_df_vari) %in% IFN_genes2,]
sig_df_combo_IFN <- sig_df_combo[rownames(sig_df_combo) %in% IFN_genes2,]

plot(-log10(sig_df_vari_IFN[,6]),-log10(sig_df_ica_IFN[,6]))

tmp <- cbind.data.frame(-log10(sig_df_vari[,2]),-log10(sig_df_ica[,2]))
tmp <- cbind.data.frame(-log10(sig_df_vari_IFN[,2]),-log10(sig_df_ica_IFN[,2]))
tmp <- cbind.data.frame(-log10(sig_df_ica_IFN[names(all_pv2),3]),-log10(all_pv2))
tmp <- cbind.data.frame(-log10(sig_df_vari_IFN[names(all_pv2),3]),-log10(all_pv2))
tmp <- cbind.data.frame(-log10(sig_df_combo_IFN[names(all_pv2),3]),-log10(all_pv2))

colnames(tmp) <- c('vari','ica')
colnames(tmp) <- c('ica','true')
colnames(tmp) <- c('vari','true')
colnames(tmp) <- c('combo','true')
# ggplot(tmp,aes(x=vari,y=ica)) +
ggplot(tmp,aes(x=true,y=combo)) +
  geom_point(alpha=.4) +
  # xlim(c(0,9)) +
  # ylim(c(0,9)) +
  geom_vline(xintercept = -log10(.05), 
             color = "blue", size=1.5) +
  geom_hline(yintercept = -log10(.05), 
             color = "blue", size=1.5)
  

# maybe I should have done this with var power of 1.5 because it had more IFN associations for ICA comparitively






# instead of using F1 from either varimax or ICA, I will use eigengene method to get all IFN genes
core_set <- c('MX1','ISG15','IFI27','OAS1','XAF1','IFITM3')
# check that core set is highly upregulated in F1
myf <- get_one_factor(pbmc_container,factor_select=1)
dsc <- myf[[1]]
lds <- myf[[2]]
lds[core_set,]

sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=1, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df[core_set,]

all_g_lists <- c()
g_list <- c()
pv_list <- c()
for (ct in container$experiment_params$ctypes_use) {
  pb <- pbmc_container[["scMinimal_ctype"]][[ct]][["pseudobulk"]]
  pb_exp <- pb[,core_set]
  svd_res <- svd(pb_exp)
  egene <- svd_res$u[,1]
  for (i in 1:ncol(pb)) {
    gn_name <- colnames(pb)[i]
    tmp <- cbind.data.frame(egene,pb[,i])
    colnames(tmp) <- c('egene','gene')
    lmres <- lm(gene~egene,data=tmp)
    x <- summary(lmres)
    pval <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)
    g_list <- c(g_list,gn_name)
    pv_list <- c(pv_list,pval)
  }
  pv_list <- p.adjust(pv_list,method='fdr')
  g_keep <- g_list[pv_list<0.01]
  all_g_lists <- c(all_g_lists,g_keep)
}
IFN_genes <- unique(all_g_lists)
length(IFN_genes)




## checking number of significant genes for vari F7 that are in my IFN list
myf <- get_one_factor(pbmc_container,factor_select=4)
dsc <- myf[[1]]
lds <- myf[[2]]
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=4, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
f7_sig <- rownames(sig_df)[rowSums(sig_df<0.01)>0]
f7_sig <- rownames(sig_df)[sig_df[,3]<0.01]
sum(f7_sig %in% IFN_genes)
f7_sig[(f7_sig %in% IFN_genes)]
length(f7_sig)

IFN_genes <- IFN_genes[!(IFN_genes %in% f7_sig)]

myf <- get_one_factor(pbmc_container,factor_select=1)
dsc <- myf[[1]]
lds <- myf[[2]]
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=1, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
f1_sig <- rownames(sig_df)[rowSums(sig_df<0.01)>0]
f1_sig <- rownames(sig_df)[sig_df[,3]<0.01]
length(intersect(f1_sig,f7_sig))

lds[f7_sig[f7_sig %in% f1_sig],6]



myf <- get_one_factor(container,factor_select=4)
dsc_ica <- myf[[1]]
lds_ica <- myf[[2]]











## here I'm going to plot loadings values versus signed association p-values
pbmc_container <- get_lm_pvals(pbmc_container)

myfactor <- 9
myf <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- myf[[1]]
lds <- myf[[2]]
# pv_saved <- pbmc_container[["gene_score_associations"]]
pbmc_container[["gene_score_associations"]] <- pv_saved
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]

# next I get signs of associations using code at bottom
# pv2_saved <- pbmc_container[["gene_score_associations"]]
pbmc_container[["gene_score_associations"]] <- pv2_saved
sign_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sign_df <- t(as.data.frame(do.call(rbind, sign_vectors)))
sign_df <- sign_df[rownames(lds),colnames(lds)]

sig_df <- -log10(sig_df)
sign_mask <- sign_df>0
sig_df[!sign_mask] <- sig_df[!sign_mask]*-1

# now plot the results
sig_df2 <- c(sig_df)
lds2 <- c(lds)
# ## to get out specific genes
# sig_df2 <- sig_df[,6]
# lds2 <- lds[,6]
tmp <- cbind.data.frame(sig_df2,lds2)
# rownames(tmp) <- rownames(lds)
# ##
colnames(tmp) <- c('sig_val','lds_val')
ggplot(tmp,aes(x=sig_val,y=lds_val)) +
  geom_point(alpha=.3) +
  geom_vline(xintercept = -log10(.01), linetype="dotted", 
               color = "blue", size=1.5) +
  geom_vline(xintercept = log10(.01), linetype="dotted", 
             color = "blue", size=1.5) +
  geom_vline(xintercept = 0, 
             color = "black", size=1.5) +
  geom_hline(yintercept = 0, 
             color = "black", size=1.5) +
  ggtitle(myfactor)

## count num in bad vs good quadrants
sig_up <- tmp[tmp$sig_val>(-log10(0.01)),]
sig_down <- tmp[tmp$sig_val<(log10(0.01)),]
sum(sig_up$lds_val<0)
sum(sig_up$lds_val>0)
sum(sig_down$lds_val>0)
sum(sig_down$lds_val<0)

gtest <- rownames(sig_down[sig_down$lds_val>0,])
length(gtest)
sum(gtest %in% IFN_genes)
gtest
## seeing if gtest have higher than average norm_var pvals
mean(pbmc_container[["scMinimal_ctype"]][["T4"]][["norm_variances"]][gtest])
mean(pbmc_container[["scMinimal_ctype"]][["T4"]][["norm_variances"]][rownames(tmp)])

## stuff below is to get the sign of gene-factor associations
container <- pbmc_container
tensor_data <- container$tensor_data
tucker_results <- container$tucker_results

n_genes <- length(tensor_data[[2]])
n_ctypes <- length(tensor_data[[3]])
fstats_real <- data.frame(matrix(ncol=3,nrow=0))
ncores <- 20
fstats_real <- mclapply(1:n_genes, function(i) {
  gene <- tensor_data[[2]][i]
  gene_res <- list()
  for (j in 1:n_ctypes) {
    ctype <- tensor_data[[3]][j]
    gene_res[[ctype]] <- list()
    for (k in 1:ncol(tucker_results[[1]])) {
      tmp_fiber <- tensor_data[[4]][,i,j]
      df_test <- as.data.frame(cbind(tmp_fiber, tucker_results[[1]][,k]))
      colnames(df_test) <- c('fiber','factor')
      lmres <- lm(factor~fiber,df_test)
      
      # testing getting pval instead of fstat
      x <- summary(lmres)
      gene_res[[ctype]][[as.character(k)]] <- x$coefficients['fiber','Estimate']
    }
  }
  return(gene_res)
}, mc.cores = ncores)

names(fstats_real) <- tensor_data[[2]]

# unpack the list
fstats_real <- unlist(fstats_real)

pbmc_container[["gene_score_associations"]] <- fstats_real















##### making ICA/varimax comparison plot
# with ICA first
pbmc_container <- get_lm_pvals(pbmc_container)

myfactor <- 6
myf <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- myf[[1]]
lds <- myf[[2]]
dsc_ica_6 <- dsc
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df_ica_t4t8 <- sig_df

myfactor <- 4
myf <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- myf[[1]]
lds <- myf[[2]]
dsc_ica_4 <- dsc
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df_ica_pred <- sig_df


# now varimax
pbmc_container <- get_lm_pvals(pbmc_container)

myfactor <- 5
myf <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- myf[[1]]
lds <- myf[[2]]
dsc_vari_5 <- dsc
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df_vari_t4t8 <- sig_df

myfactor <- 4
myf <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- myf[[1]]
lds <- myf[[2]]
dsc_vari_4 <- dsc
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df_vari_pred <- sig_df


## for pred factor
# to use all genes
sig_df_ica_IFN <- sig_df_ica_pred
sig_df_vari_IFN <- sig_df_vari_pred

# using all F1 significant genes
sig_df_ica_IFN <- sig_df_ica_pred[rownames(sig_df_ica_pred) %in% IFN_genes,]
sig_df_vari_IFN <- sig_df_vari_pred[rownames(sig_df_vari_pred) %in% IFN_genes,]

# using only F1 significant genes in GO IFN sets
sig_df_ica_IFN <- sig_df_ica_pred[rownames(sig_df_ica_pred) %in% f1_ifn,]
sig_df_vari_IFN <- sig_df_vari_pred[rownames(sig_df_vari_pred) %in% f1_ifn,]

## for t4t8 factor
# to use all genes
sig_df_ica_IFN <- sig_df_ica_t4t8
sig_df_vari_IFN <- sig_df_vari_t4t8

# using all F1 significant genes
sig_df_ica_IFN <- sig_df_ica_t4t8[rownames(sig_df_ica_t4t8) %in% IFN_genes,]
sig_df_vari_IFN <- sig_df_vari_t4t8[rownames(sig_df_vari_t4t8) %in% IFN_genes,]

# using only F1 significant genes in GO IFN sets
sig_df_ica_IFN <- sig_df_ica_t4t8[rownames(sig_df_ica_t4t8) %in% f1_ifn,]
sig_df_vari_IFN <- sig_df_vari_t4t8[rownames(sig_df_vari_t4t8) %in% f1_ifn,]

# ctype <- 'T8'
# ct_ndx <- which(colnames(sig_df_vari_IFN)==ctype)
# tmp <- cbind.data.frame(-log10(sig_df_vari_IFN[,ct_ndx]),-log10(sig_df_ica_IFN[,ct_ndx]))
# 
# colnames(tmp) <- c('vari','ica')

# for comparing combo to vari
sig_df_vari_IFN <- sig_df_vari_t4t8
sig_df_combo_IFN <- sig_df_combo_t4t8
ctype <- 'NK'
ct_ndx <- which(colnames(sig_df_vari_IFN)==ctype)
tmp <- cbind.data.frame(-log10(sig_df_vari_IFN[,ct_ndx]),-log10(sig_df_combo_IFN[,ct_ndx]))
colnames(tmp) <- c('vari','combo')

ggplot(tmp,aes(x=vari,y=combo)) +
  geom_point(alpha=.2) +
  geom_vline(xintercept = -log10(.05), 
             color = "blue", size=1.5) +
  geom_hline(yintercept = -log10(.05), 
             color = "blue", size=1.5) + 
  ggtitle(ctype) +
  theme(plot.title = element_text(hjust = 0.5))

sig_ica <- rownames(sig_df_ica_IFN)[sig_df_ica_IFN[,ct_ndx]<.05]
sig_vari_sub <- sig_df_vari_IFN[sig_ica,]
opp_genes <- rownames(sig_vari_sub)[sig_vari_sub[,ct_ndx]>.05]

pval <- stats::phyper(sum(opp_genes %in% f1_ifn)-1, length(f1_ifn), nrow(sig_df_ica_IFN) - length(f1_ifn),
                      length(opp_genes), lower.tail = FALSE)

pval <- stats::phyper(sum(opp_genes %in% IFN_genes)-1, length(IFN_genes), nrow(sig_df_ica_IFN) - length(IFN_genes),
                      length(opp_genes), lower.tail = FALSE)


# comparing dscores for ica vs vari
identical(rownames(dsc_ica_4),rownames(dsc_vari_4))
tmp <- cbind.data.frame(dsc_ica_4,dsc_vari_4)
colnames(tmp) <- c('ica','vari')
ggplot(tmp,aes(x=vari,y=ica)) +
  geom_point(alpha=.4)

identical(rownames(dsc_ica_6),rownames(dsc_vari_5))
tmp <- cbind.data.frame(dsc_ica_6,dsc_vari_5)
colnames(tmp) <- c('ica','vari')
ggplot(tmp,aes(x=vari,y=ica)) +
  geom_point(alpha=.4)


## now running results for combo method and for IFN-negative dataset for other comparisons
# a question I still have is whether or not to use the whole F1 when I compare these distributions
# or whether I should just use the GO f1 genes. 

# with combo method
pbmc_container <- get_lm_pvals(pbmc_container)

myfactor <- 6
myf <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- myf[[1]]
lds <- myf[[2]]
dsc_combo_6 <- dsc
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df_combo_t4t8 <- sig_df

sig_df_combo_IFN <- sig_df_combo_t4t8[rownames(sig_df_combo_t4t8) %in% IFN_genes,]
sig_df_combo_IFN <- sig_df_combo_t4t8[rownames(sig_df_combo_t4t8) %in% f1_ifn,]

myfactor <- 4
myf <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- myf[[1]]
lds <- myf[[2]]
dsc_combo_4 <- dsc
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df_combo_pred <- sig_df

sig_df_combo_IFN <- sig_df_combo_pred[rownames(sig_df_combo_pred) %in% IFN_genes,]
sig_df_combo_IFN <- sig_df_combo_pred[rownames(sig_df_combo_pred) %in% f1_ifn,]


# now getting the decomposition without any IFN genes present
pbmc_container <- get_lm_pvals(pbmc_container)

myfactor <- 7
myf <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- myf[[1]]
lds <- myf[[2]]
dsc_no_ifn_7 <- dsc

myfactor <- 3
myf <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- myf[[1]]
lds <- myf[[2]]
dsc_no_ifn_3 <- dsc


myf <- get_one_factor(pbmc_container,factor_select=5)
dsc <- myf[[1]]
lds <- myf[[2]]

ctype <- 'T4'
pb <- pbmc_container[["scMinimal_ctype"]][[ctype]][["pseudobulk"]]
all_pv <- c()
for (gn in IFN_genes) {
# for (gn in f1_ifn) {
  tmp <- cbind.data.frame(dsc,pb[rownames(dsc),gn])
  # tmp <- cbind.data.frame(dsc_for_pred,pb[rownames(dsc_for_pred),gn])
  # tmp <- cbind.data.frame(dsc_for_t4,pb[rownames(dsc_for_t4),gn])
  colnames(tmp) <- c('dscores','expr')
  lmres <- lm(expr~dscores,data=tmp)
  x <- summary(lmres)
  pval <- stats::pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3],lower.tail=FALSE)
  all_pv <- c(all_pv,pval)
}
names(all_pv) <- IFN_genes
# names(all_pv) <- f1_ifn
all_pv2 <- p.adjust(all_pv,method='fdr')

# tmp <- cbind.data.frame(-log10(all_pv2_ica) ,-log10(all_pv2_vari))
tmp <- cbind.data.frame(-log10(sig_df_combo[names(all_pv2),7]) ,-log10(all_pv2))
colnames(tmp) <- c('combo','true')
ggplot(tmp,aes(x=true,y=combo)) + geom_point(alpha=.2) + xlim(c(0,5)) + ylim(c(0,5)) +
  geom_vline(xintercept = -log10(.05),
             color = "blue", size=1.5) +
  geom_hline(yintercept = -log10(.05),
             color = "blue", size=1.5)

ct_ndx <- which(colnames(sig_df_vari_IFN)==ctype)

tmp <- cbind.data.frame(-log10(sig_df_ica_IFN[names(all_pv2),ct_ndx]),-log10(all_pv2))
tmp <- cbind.data.frame(-log10(sig_df_vari_IFN[names(all_pv2),ct_ndx]),-log10(all_pv2))
tmp <- cbind.data.frame(-log10(sig_df_combo_IFN[names(all_pv2),ct_ndx]),-log10(all_pv2))

tmp <- cbind.data.frame(-log10(sig_df_combo[names(all_pv2),6]),-log10(all_pv2))


colnames(tmp) <- c('ica','true')
colnames(tmp) <- c('vari','true')
colnames(tmp) <- c('combo','true')
ggplot(tmp,aes(x=true,y=combo)) +
  geom_point(alpha=.2) +
  # xlim(c(0,9)) +
  # ylim(c(0,9)) +
  geom_vline(xintercept = -log10(.05), 
             color = "blue", size=1.5) +
  geom_hline(yintercept = -log10(.05), 
             color = "blue", size=1.5)


# calculate number of false positives
# for ica
tmp_sub <- tmp[tmp$ica>(-log10(.05)),]
sum(tmp_sub$true<(-log10(.05)))
# for varimax
tmp_sub <- tmp[tmp$vari>(-log10(.05)),]
sum(tmp_sub$true<(-log10(.05)))
# for combo
tmp_sub <- tmp[tmp$combo>(-log10(.05)),]
sum(tmp_sub$true<(-log10(.05)))

# calculate number of false negatives
# for ica
tmp_sub <- tmp[tmp$ica<(-log10(.05)),]
sum(tmp_sub$true>(-log10(.05)))
# for varimax
tmp_sub <- tmp[tmp$vari<(-log10(.05)),]
sum(tmp_sub$true>(-log10(.05)))
# for combo
tmp_sub <- tmp[tmp$combo<(-log10(.05)),]
sum(tmp_sub$true>(-log10(.05)))

# showing dscore correlations between combo method and others
identical(rownames(dsc_vari_4),rownames(dsc_combo_4))
tmp <- cbind.data.frame(dsc_vari_4,dsc_combo_4)
colnames(tmp) <- c('vari','combo')
ggplot(tmp,aes(x=vari,y=combo)) +
  geom_point(alpha=.4)

tmp <- cbind.data.frame(dsc_vari_5,dsc_combo_6)
colnames(tmp) <- c('vari','combo')
ggplot(tmp,aes(x=vari,y=combo)) +
  geom_point(alpha=.4)





# trying GSEA with the different methods to see how the resuls differ
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,16,7),
                                 tucker_type = 'regular', rotation_type = 'ica')

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Status'),
                                        stat_use='pval')

pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    cluster_by_meta = 'sex',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

pbmc_container$plots$donor_matrix

pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=2,direc='up',thresh=.05)
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=2,direc='down',thresh=.05)

# investigate several clusters of enriched sets
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=3)


pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=6, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=6,direc='down',thresh=.05)


pbmc_container <- run_cv_cors(pbmc_container, ranks=c(7,16,7), n_splits=4)

pbmc_container$plots$stability_plot_dsc
pbmc_container$plots$stability_plot_lds








## for comparing significant genes from two decompositions
pbmc_container <- get_lm_pvals(pbmc_container)

## need to save gene score associations and tucker results
# fstats_real_ica <- pbmc_container[["gene_score_associations"]]
# fstats_real_vari <- pbmc_container[["gene_score_associations"]]
# fstats_real_combo <- pbmc_container[["gene_score_associations"]]
pbmc_container[["gene_score_associations"]] <- fstats_real_ica
pbmc_container[["gene_score_associations"]] <- fstats_real_vari
pbmc_container[["gene_score_associations"]] <- fstats_real_combo

# tucker_ica <- pbmc_container[["tucker_results"]]
# tucker_vari <- pbmc_container[["tucker_results"]]
# tucker_combo <- pbmc_container[["tucker_results"]]
pbmc_container[["tucker_results"]] <- tucker_ica
pbmc_container[["tucker_results"]] <- tucker_vari
pbmc_container[["tucker_results"]] <- tucker_combo


myfactor <- 4
pbmc_container[["gene_score_associations"]] <- fstats_real_ica
pbmc_container[["tucker_results"]] <- tucker_ica
myf <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- myf[[1]]
lds <- myf[[2]]
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df_ica <- sig_df
# IFN_genes <- rownames(sig_df)[rowSums(sig_df<0.01)>0]


myfactor <- 4
pbmc_container[["gene_score_associations"]] <- fstats_real_vari
pbmc_container[["tucker_results"]] <- tucker_vari
myf <- get_one_factor(pbmc_container,factor_select=myfactor)
dsc <- myf[[1]]
lds <- myf[[2]]
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df_vari <- sig_df

ctype <- 'T8'
ct_ndx <- which(colnames(sig_df_ica)==ctype)
tmp <- cbind.data.frame(-log10(sig_df_vari[,ct_ndx]),-log10(sig_df_ica[,ct_ndx]))
colnames(tmp) <- c('hybrid','ica')

ggplot(tmp,aes(x=hybrid,y=ica)) +
  geom_point(alpha=.2) +
  geom_vline(xintercept = -log10(.05), 
             color = "blue", size=1.5) +
  geom_hline(yintercept = -log10(.05), 
             color = "blue", size=1.5) + 
  ggtitle(ctype) +
  theme(plot.title = element_text(hjust = 0.5))

# see what the genes are that vari "underpredicts"
tmp_sub <- tmp[tmp$ica>(-log10(.05)),]
opp_genes <- rownames(tmp_sub)[tmp_sub$hybrid<(-log10(.05))]
length(opp_genes)
sum(opp_genes %in% IFN_genes)

## seeing if these genes are enriched for pred gene sets
all_genes <- rownames(tmp)
total_num_genes <- length(all_genes)
sig_genes <- opp_genes
my_pathways <- my_pathways[c('GO_RESPONSE_TO_STEROID_HORMONE',"GO_RESPONSE_TO_CORTICOSTEROID","GO_RESPONSE_TO_LIPID","GO_RESPONSE_TO_ESTRADIOL")]
my_pathways <- my_pathways[c('GO_RESPONSE_TO_TYPE_I_INTERFERON',"GO_RESPONSE_TO_INTERFERON_GAMMA")]
padj[order(padj)][1:30]

# testing for enrichment of F1 genes in opp_genes
pval <- stats::phyper(sum(opp_genes %in% IFN_genes)-1, length(IFN_genes), total_num_genes - length(IFN_genes),
                      length(opp_genes), lower.tail = FALSE)














myfactor <- 1
f1_data <- get_one_factor(pbmc_container, factor_select=myfactor)
dscores <- f1_data[[1]]
lds <- f1_data[[2]]
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
IFN_genes <- rownames(sig_df)[rowSums(sig_df<0.1)>0]







# run gsea
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=5, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=5,direc='up',thresh=.05)

# investigate several clusters of enriched sets
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=12)






# running stability analysis
pbmc_container <- run_stability_analysis(pbmc_container, ranks=c(7,16,7), downsample_ratio=0.9, n_iter=5,
                                         norm_method='trim', scale_factor=10000,
                                         batch_var='pool', scale_var=TRUE,
                                         var_scale_power=1.5, tucker_type='regular',
                                         rotation_type='ica')
pbmc_container$plots$stability_plot_dsc
pbmc_container$plots$stability_plot_lds





# get factor-meta data associations
container_sub <- get_meta_associations(container_sub,vars_test=c('sex','Age','pool','processing','Status'),
                                        stat_use='pval')

# plot donor scores by status
container_sub <- plot_donor_matrix(container_sub, meta_vars=c('sex'),
                                    cluster_by_meta = 'sex',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

container_sub$plots$donor_matrix


d_stored <- pbmc_container$tensor_data
pbmc_container <- run_cv_cors(pbmc_container,ranks=c(7,16,7),n_splits=4)
pbmc_container <- run_cv_cors(pbmc_container,ranks=c(10,20,7),n_splits=4)
pbmc_container <- run_cv_cors(pbmc_container,ranks=c(12,20,7),n_splits=4)
container$plots$stability_plot_dsc
container$plots$stability_plot_lds

# mydsc <- pbmc_container$tucker_results[[1]]
# cv_res <- run_cv(pbmc_container,ranks=c(7,16,7),n_splits=4)
# f_sel=4
# cor(cv_res[[f_sel]],mydsc[rownames(cv_res[[f_sel]]),])



pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(7,20,7),n_iterations=50,subset_type='bootstrap')
pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(7,20,7),n_iterations=500,subset_type='subset', sub_prop=.75)
pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(9,20,7),n_iterations=10,subset_type='subset', sub_prop=.75)

pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(10,20,7),n_iterations=50,subset_type='bootstrap')
pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(10,20,7),n_iterations=50,subset_type='subset', sub_prop=.75)

pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(7,16,7),n_iterations=50,subset_type='bootstrap',
                                         tucker_type='sparse',sparsity=2)
pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(11,20,7),n_iterations=50,subset_type='bootstrap',
                                         tucker_type='sparse',sparsity=sqrt(2))
pbmc_container <- run_stability_analysis(pbmc_container,ranks=c(4,10,5),n_iterations=50,subset_type='subset', sub_prop=.75,)

pdf(file = "/home/jmitchel/figures/test.pdf", useDingbats = FALSE,
    width = 5, height = 6)
pbmc_container$plots$stability_plot_dsc
dev.off()

pbmc_container$plots$stability_plot_lds





















## running LR analysis with new functions
# prep for new LR analysis
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
  rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
  return(rname)
})
lr_pairs$interaction_name <- NULL


pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='pool')
sft_thresh <- c(3,3,5,2,2,2,2) #unsigned
sft_thresh <- c(9,9,6,5,8,5,6) #signed
sft_thresh <- c(12,12,12,10,12,9,10) #signed new
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

pbmc_container <- compute_LR_interact(pbmc_container, lr_pairs, factor_select=2, 
                                      sig_thresh=0.01, percentile_exp_rec=.85,
                                      show_rec_sig=FALSE)

pbmc_container$plots$lr_analysis[['Factor2']]


# getting GO enrichment HMAPs for modules
ctypes <- c('T4','T8','T4')
modules <- c(7,3,5)
# modules <- c(3,5)
ctypes <- c('T8')
modules <- c(8)

mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use='TF')
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('GO'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('Reactome'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('KEGG'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('BioCarta'))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.05, db_use=c('Hallmark'))
mod_enr

lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,mod_ct='T4',mod=8,lig_ct='cM',lig='ICOSLG')
lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=5,mod_ct='cM',mod=3,lig_ct='T4',lig='TNFSF8')
lig_mod_fact







## looking at correlation of lm vs jackstraw gene pvalues
pbmc_container <- get_lm_pvals(pbmc_container)
lm_only <- container[["gene_score_associations"]]

pbmc_container <- run_jackstraw_v2(pbmc_container, ranks=c(7,20,7), n_fibers=250, n_iter=1000,
                                   tucker_type='regular', rotation_type='ica')
w_jack <- pbmc_container[["gene_score_associations"]]

# saveRDS(w_jack,file='/home/jmitchel/data/lupus_data/lupus_jackstraw_lds_rot_hybrid.rds')


tmp <- cbind.data.frame(w_jack,lm_only)
colnames(tmp) <- c('jackstraw','lm_simple')
ggplot(tmp,aes(x=-log10(lm_simple),y=-log10(jackstraw))) +
  geom_point() +
  xlim(0,4) +
  ylim(0,4)

# looking at just the F1 sig genes because not sure if the other ones matched up

container[["gene_score_associations"]] <- w_jack

myfactor <- 1
f1_data <- get_one_factor(pbmc_container, factor_select=myfactor)
dscores <- f1_data[[1]]
lds <- f1_data[[2]]
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df_jack <- c(sig_df)

container[["gene_score_associations"]] <- lm_only

myfactor <- 1
f1_data <- get_one_factor(pbmc_container, factor_select=myfactor)
dscores <- f1_data[[1]]
lds <- f1_data[[2]]
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=myfactor, colnames(lds))
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
sig_df <- sig_df[rownames(lds),colnames(lds)]
sig_df_lm <- c(sig_df)





tmp <- cbind.data.frame(sig_df_jack,sig_df_lm)
colnames(tmp) <- c('jackstraw','lm_simple')
ggplot(tmp,aes(x=-log10(lm_simple),y=-log10(jackstraw))) +
  geom_point() +
  xlim(0,6) +
  ylim(0,6)





