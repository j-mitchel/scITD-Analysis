


# counts matrix
pbmc_counts <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_counts_v2.rds')

# meta data matrix
pbmc_meta <- readRDS('/home/jmitchel/data/van_der_wijst/pbmc_meta_v2.rds')

# ensembl to gene name conversions
feature.names <- readRDS('/home/jmitchel/data/van_der_wijst/genes.rds')

# set up project parameters
param_list <- initialize_params(ctypes_use = c("CD4+ T", "CD8+ T", "cMonocyte", "CD56(dim) NK", "B"),
                                ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(count_data=pbmc_counts, meta_data=pbmc_meta,
                                     gn_convert = feature.names, params=param_list,
                                     label_donor_sex = TRUE)
pbmc_container <- parse_data_by_ctypes(pbmc_container)
pbmc_container <- clean_data(pbmc_container, donor_min_cells=5, gene_min_cells=5)
pbmc_container <- get_pseudobulk(pbmc_container)
pbmc_container <- normalize_pseudobulk(pbmc_container, method='trim', scale_factor=1000)

scMinimal <- pbmc_container$scMinimal_ctype[['CD4+ T']]
donor_sum_counts <- scMinimal$pseudobulk

df <- colMeanVars(donor_sum_counts, rowSel = NULL)
df$m <- log(df$m); df$v <- log(df$v);
rownames(df) <- colnames(donor_sum_counts);

gam.k <- 5
min.gene.cells <- 0
vi <- which(is.finite(df$v) & df$nobs>=min.gene.cells);
if(length(vi)<gam.k*1.5) { gam.k=1 };# too few genes
if(gam.k<2) {
  m <- lm(v ~ m, data = df[vi,])
} else {
  m <- mgcv::gam(stats::as.formula(paste0('v ~ s(m, k = ',gam.k,')')), data = df[vi,])
}

df$res <- -Inf 
df$res <- stats::resid(m,type='response')
n.obs <- df$nobs;
suppressWarnings(df$lp <- as.numeric(stats::pf(exp(df$res),n.obs,n.obs,lower.tail=FALSE,log.p=FALSE)))
df$lpa <- log(p.adjust(df$lp,method='fdr'))
df$lp <- log(df$lp)
n.cells <- nrow(donor_sum_counts)
scaled_var <- as.numeric(stats::qchisq(df$lp, n.cells-1, lower.tail = FALSE,log.p=TRUE)/n.cells)
names(scaled_var) <- colnames(donor_sum_counts)

df$qv <- as.numeric(stats::qchisq(df$lp, n.cells-1, lower.tail = FALSE,log.p=TRUE)/n.cells)

ods <- which(df$lpa<log(5e-2))

par(mfrow=c(1,2), mar = c(3.5,3.5,2.0,0.5), mgp = c(2,0.65,0), cex = 1.0)
suppressWarnings(smoothScatter(log10(exp(1))*df$m, log10(exp(1))*df$v, main='', xlab='log10[ magnitude ]',ylab='log10[ variance ]'))
vi <- which(is.finite(log10(exp(1))*df$v) & df$nobs>=min.gene.cells)
grid <- seq(min(log10(exp(1))*df$m[vi]), max(log10(exp(1))*df$m[vi]), length.out=1000)
## re-calculate m
if (gam.k < 2) {
  if (verbose) message(" using lm ")
  m <- lm(v ~ m, data = log10(exp(1))*df[vi,])
} else {
  if (verbose) message(" using gam ")
  m <- mgcv::gam(as.formula(paste0('v ~ s(m, k = ',gam.k,')')), data = log10(exp(1))*df[vi,])
}
lines(grid,predict(m, newdata=data.frame(m=grid)), col="blue")
# points(log10(exp(1))*df$m[1:1000], log10(exp(1))*df$v[1:1000], pch='.',col=2,cex=1)
points(log10(exp(1))*df$m[ods], log10(exp(1))*df$v[ods], pch='.',col=2,cex=1)
# suppressWarnings(smoothScatter(log10(exp(1))*df$m[vi], log10(exp(1))*df$qv[vi], xlab='log10[ magnitude ]',ylab='',main='adjusted'))
suppressWarnings(smoothScatter(log10(exp(1))*df$m[vi], df$qv[vi], xlab='log10[ magnitude ]',ylab='',main='adjusted'))
abline(h=1,lty=2,col=8)
abline(h=1e3, lty=2, col=1)
# points(log10(exp(1))*df$m[1:1000], log10(exp(1))*df$qv[1:1000], col=2, pch='.')
# points(log10(exp(1))*df$m, log10(exp(1))*df$qv, col=2, pch='.')
points(log10(exp(1))*df$m[ods], df$qv[ods], col=2, pch='.')




# trying it with pagoda2 to see what is going on
scMinimal <- pbmc_container$scMinimal_ctype[['CD4+ T']]
counts <- scMinimal$pseudobulk

counts <- pbmc_container$scMinimal_full$count_data[,1:10000]

library(pagoda2)
r <- Pagoda2$new(counts,log.scale=TRUE)
r$adjustVariance(plot=T,gam.k=10)


library(pagoda2)
cd <- read.10x.matrices(list(PBMC8K='/home/jmitchel/GRCh38/'))
counts <- gene.vs.molecule.cell.filter(cd,min.cell.size=500)
counts <- counts[rowSums(counts)>=10,]
rownames(counts) <- make.unique(rownames(counts))
r <- Pagoda2$new(counts,log.scale=T, n.cores=20)
r$adjustVariance(plot=T,gam.k=10)



















