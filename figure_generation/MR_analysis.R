
library(Seurat)
library(MatrixEQTL)
library(sccore)
library(ggplot2)
library(locuszoomr)
library(EnsDb.Hsapiens.v75)
library(ivreg)
library(OpenImageR)
library(preprocessCore)
library(devtools)
load_all('/home/jmitchel/scITD/')
source('/home/jmitchel/smr/plot/plot_SMR.r')

# load the SNP data
vcf_mat <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_total_vcf.rds')

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


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20),
                                 tucker_type = 'regular', rotation_type = 'hybrid')


# flip sign of F1 so high ISG expression is positive instead of negative (signs are arbitrary)
pbmc_container$tucker_results[[1]][,1] <- pbmc_container$tucker_results[[1]][,1] * -1
pbmc_container$tucker_results[[2]][1,] <- pbmc_container$tucker_results[[2]][1,] * -1
pbmc_container$projection_data[[1]][1,] <- pbmc_container$projection_data[[1]][1,] * -1


# get matrix of covariates to include in the model
pbmc_container <- get_donor_meta(pbmc_container,additional_meta = c('sex','Age',
                                                                    'Ethnicity','processing'))

d_meta <- pbmc_container$donor_metadata
new_names <- sapply(rownames(d_meta),function(x){
  strsplit(x,split='dmx')[[1]][[1]]
})
rownames(d_meta) <- new_names


## ICOSLG chr21:45,642,874 (start)
gene_start <- 45642874
window_start <- gene_start-250000
window_end <- gene_start+250000

# picking genetic variants with 200kb of the gene start
chr_labs <- sapply(rownames(vcf_mat),function(x){
  strsplit(x,split=':')[[1]][[1]]
})

# subset to chr21
vcf_mat_sub <- vcf_mat[chr_labs=='21',]

pos_labs <- sapply(rownames(vcf_mat_sub),function(x){
  strsplit(x,split=':')[[1]][[2]]
})
pos_labs <- as.numeric(pos_labs)

ndx_keep1 <- which(pos_labs > window_start)
ndx_keep2 <- which(pos_labs < window_end)
ndx_keep <- intersect(ndx_keep1,ndx_keep2)
snps_keep <- rownames(vcf_mat_sub)[ndx_keep]
pos_keep <- pos_labs[ndx_keep]
vcf_mat_sub <- vcf_mat_sub[snps_keep,]


# get gene expression to use in eQTL mapping
pb <- pbmc_container[["scMinimal_ctype"]][["cMono"]][["pseudobulk"]]
new_names <- sapply(rownames(pb),function(x){
  strsplit(x,split='dmx')[[1]][[1]]
})
rownames(pb) <- new_names

d_both <- intersect(rownames(pb),colnames(vcf_mat_sub))
pb <- scale(pb)
pb <- t(pb)
snps1 <- SlicedData$new(as.matrix(vcf_mat_sub[,d_both]))
gene1_input <- pb['ICOSLG',d_both,drop=FALSE]
gene1 <- SlicedData$new(gene1_input)

# create donor level covariates matrix
d_meta <- d_meta[d_both,]
d_meta$donors <- NULL
d_meta$Ethnicity <- factor(d_meta$Ethnicity,levels=c('European','Asian'))
d_meta$processing <- factor(d_meta$processing,levels=unique(d_meta$processing))

# ##### testing adding the lead snp for first peak as a covariate!
# # using just the first peak
# d_meta$g_covar <- t(vcf_mat_sub[c('21:45624136:C:T'),rownames(d_meta)])
# # using both lead eQTL snps
# d_meta$g_covar <- t(vcf_mat_sub[c('21:45624136:C:T',"21:45660835:C:A"),rownames(d_meta)])
# #####

mm <- model.matrix(~., d_meta)
mm <- mm[,colSums(mm)!=0]
mm <- mm[,2:ncol(mm)] # remove intercept term
mm <- t(mm)
cvrt1 <- SlicedData$new(mm)

filename <- NULL
me <- Matrix_eQTL_main(
  snps = snps1,
  gene = gene1,
  cvrt = cvrt1,
  output_file_name = filename,
  pvOutputThreshold = 1,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = FALSE,
  pvalue.hist = FALSE )

head(me[["all"]][["eqtls"]])
e_res <- me[["all"]][["eqtls"]]
rownames(e_res) <- e_res$snps
e_res$pos <- sapply(e_res$snps,function(x){
  strsplit(x,split=':')[[1]][[2]]
})
e_res$pos <- as.numeric(e_res$pos)
e_res$chr <- sapply(e_res$snps,function(x){
  strsplit(x,split=':')[[1]][[1]]
})
e_res$chr <- as.numeric(e_res$chr)
ggplot(e_res,aes(x=pos,y=-log10(pvalue))) +
  geom_point()







## get the factor to test MR with
# first recompute factor without ICOSLG included
umi_mat <- pbmc_container[["scMinimal_full"]][["count_data"]]
g_keep <- rownames(umi_mat)[rownames(umi_mat)!='ICOSLG']
umi_mat <- umi_mat[g_keep,]
pbmc_container[["scMinimal_full"]][["count_data"]] <- umi_mat
tucker_stored <- pbmc_container$tucker_results
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20),
                                 tucker_type = 'regular', rotation_type = 'hybrid')

print(cor(tucker_stored[[1]],pbmc_container$tucker_results[[1]])) # factor 2 is almost identical

# flipping sign of F2 so it's the same as before
pbmc_container$tucker_results[[1]][,2] <- pbmc_container$tucker_results[[1]][,2] * -1
pbmc_container$tucker_results[[2]][2,] <- pbmc_container$tucker_results[[2]][2,] * -1

f_test <- get_one_factor(pbmc_container,2)[[1]] # factor 2 dscores
new_names <- sapply(rownames(f_test),function(x){
  strsplit(x,split='dmx')[[1]][[1]]
})
rownames(f_test) <- new_names
f_test <- scale(f_test)
f_test <- t(f_test)
dscores_input <- f_test[1,d_both,drop=FALSE]
dscores <- SlicedData$new(dscores_input)

pr_eQTL <- Matrix_eQTL_main(
  snps = snps1,
  gene = dscores,
  cvrt = cvrt1,
  output_file_name = filename,
  pvOutputThreshold = 1,
  useModel = modelLINEAR,
  errorCovariance = numeric(),
  verbose = FALSE,
  pvalue.hist = FALSE )

head(pr_eQTL[["all"]][["eqtls"]])
pr_res <- pr_eQTL[["all"]][["eqtls"]]
rownames(pr_res) <- pr_res$snps
pr_res$pos <- sapply(pr_res$snps,function(x){
  strsplit(x,split=':')[[1]][[2]]
})
pr_res$pos <- as.numeric(pr_res$pos)
pr_res$chr <- sapply(pr_res$snps,function(x){
  strsplit(x,split=':')[[1]][[1]]
})
pr_res$chr <- as.numeric(pr_res$chr)

ggplot(pr_res,aes(x=pos,y=-log10(pvalue))) +
  geom_point()




## getting LD with lead SNP to show on locus zoom plot
# index_snp_show <- e_res$snps[1] # choosing the lead snp of the second peak 21:45660835:C:A
index_snp_show <- e_res$snps[2] # choosing the lead snp of the first peak 21:45624136:C:T

snp_cormat <- cor(t(vcf_mat_sub[e_res$snps,]))**2
e_lead_ld <- snp_cormat[index_snp_show,]
e_res$r2 <- e_lead_ld
snp_cormat <- cor(t(vcf_mat_sub[pr_res$snps,]))**2
pr_lead_ld <- snp_cormat[index_snp_show,]
pr_res$r2 <- pr_lead_ld

loc1 <- locus(e_res, gene = 'ICOSLG', flank = 5e4,LD = "r2",
              ens_db = "EnsDb.Hsapiens.v75")
# p1 <- locus_plot(loc1,filter_gene_biotype = "protein_coding")


loc2 <- locus(pr_res, gene = 'ICOSLG', flank = 5e4,LD = "r2",
              ens_db = "EnsDb.Hsapiens.v75")
# p2 <- locus_plot(loc2,filter_gene_biotype = "protein_coding")

### Figure 4b
pdf("/home/jmitchel/figures/scITD_revision_figs2/MR_plot_tmp.pdf", width = 5, height = 3)
oldpar <- set_layers(2)
scatter_plot(loc1, xticks = FALSE,index_snp=index_snp_show)
scatter_plot(loc2, xticks = FALSE,index_snp=index_snp_show)
genetracks(loc2,filter_gene_biotype = "protein_coding")
par(oldpar)  # revert par() settings
dev.off()







########## running jlim colocalization test for each peak
library(jlimR)
source('/home/jmitchel/jlim/jlimR/R/jlim2.R')

main_tr <- pr_res
rownames(main_tr) <- main_tr$snps
main_tr <- main_tr[,c('chr','pos','pvalue')]
colnames(main_tr) <- c('CHR','BP','P')
main_tr <- main_tr[match(sort(main_tr$BP),main_tr$BP),] # just sorts variants by pos
ref.LD <- '/home/jmitchel/jlim/refld.1kg.nfe.b37'
start.bp <- min(main_tr$BP)
end.bp <- max(main_tr$BP)
CHR <- main_tr$CHR[1]
min.MAF <- 0.05 # variants below this MAF are excluded
refld.file <- find.panel(ref.LD, start.bp=start.bp,
                         end.bp=end.bp, CHR=CHR)
refgt <- read.delim(file=refld.file, header=FALSE, sep="\t", stringsAsFactors=FALSE)

verbose <- F
min.SNPs.count <- 15 # only runs jlim if there are at least this many snps in the window
min.pvalue <- 1 # only runs jlim if there is an eQTL pval below this. Set to 1 so all loci are tested.
max.perm <-  100000


## make sec_tr df
sec_tr <- cbind.data.frame(main_tr$CHR,main_tr$BP,e_res[rownames(main_tr),'pvalue'])
colnames(sec_tr) <- c('CHR','BP','P')

# running the analysis
sectr.ref.db <- c("")
secld.file <- NULL
sectr.sample.size <- length(gene1_input)

sec_tr_saved <- sec_tr
main_tr_saved <- main_tr

# testing specifically the first peak
sec_tr <- sec_tr[sec_tr$BP>45590000,]
sec_tr <- sec_tr[sec_tr$BP<45640000,]
main_tr <- main_tr[main_tr$BP>45590000,]
main_tr <- main_tr[main_tr$BP<45640000,]
indexSNP <- e_res$pos[2]

results.allgenes <- jlim.test2(maintr=main_tr, sectr=sec_tr, refgt=refgt,
                               secld.file=secld.file,CHR = CHR,indSNP=indexSNP, sectr.ref.db=sectr.ref.db,
                               start.bp=start.bp, end.bp=end.bp, min.MAF=min.MAF,
                               sectr.sample.size=sectr.sample.size,
                               min.SNPs.count=min.SNPs.count, min.pvalue=min.pvalue, perm.count=max.perm,
                               resultFileName=NULL,withPerm = FALSE)


jlim_pval1 <- as.numeric(results.allgenes[1,'pvalue'])

# testing specifically the second peak
sec_tr <- sec_tr_saved
main_tr <- main_tr_saved
sec_tr <- sec_tr[sec_tr$BP>45645000,]
sec_tr <- sec_tr[sec_tr$BP<45700000,]
main_tr <- main_tr[main_tr$BP>45645000,]
main_tr <- main_tr[main_tr$BP<45700000,]
indexSNP <- pr_res$pos[1]

results.allgenes <- jlim.test2(maintr=main_tr, sectr=sec_tr, refgt=refgt,
                               secld.file=secld.file,CHR = CHR,indSNP=indexSNP, sectr.ref.db=sectr.ref.db,
                               start.bp=start.bp, end.bp=end.bp, min.MAF=min.MAF,
                               sectr.sample.size=sectr.sample.size,
                               min.SNPs.count=min.SNPs.count, min.pvalue=min.pvalue, perm.count=max.perm,
                               resultFileName=NULL,withPerm = FALSE)


jlim_pval2 <- as.numeric(results.allgenes[1,'pvalue'])


plot(sec_tr$BP,-log10(sec_tr$P))
plot(main_tr$BP,-log10(main_tr$P))


#### now running the ivreg analysis for the top snp in the first peak
tmp <- cbind.data.frame(t(dscores_input),t(vcf_mat_sub[e_res$snps[2],d_both]),t(gene1_input),d_meta)
colnames(tmp)[1:2] <- c('dscore','gt')
fit_iv <- ivreg(dscore~sex+Age+Ethnicity+processing | ICOSLG | gt, data=tmp)
iv_res <- summary(fit_iv)
pval <- iv_res$coefficients['ICOSLG','Pr(>|t|)']


#### now running the ivreg analysis for the top eQTL snp in the second peak
tmp <- cbind.data.frame(t(dscores_input),t(vcf_mat_sub[e_res$snps[1],d_both]),t(gene1_input),d_meta)
colnames(tmp)[1:2] <- c('dscore','gt')
fit_iv <- ivreg(dscore~sex+Age+Ethnicity+processing | ICOSLG | gt, data=tmp)
iv_res <- summary(fit_iv)
pval <- iv_res$coefficients['ICOSLG','Pr(>|t|)']

#### now running the ivreg analysis for the top prQTL snp in the second peak
tmp <- cbind.data.frame(t(dscores_input),t(vcf_mat_sub[pr_res$snps[1],d_both]),t(gene1_input),d_meta)
colnames(tmp)[1:2] <- c('dscore','gt')
fit_iv <- ivreg(dscore~sex+Age+Ethnicity+processing | ICOSLG | gt, data=tmp)
iv_res <- summary(fit_iv)
pval <- iv_res$coefficients['ICOSLG','Pr(>|t|)']









### getting rsid of lead snps
variants <- read.table("/home/jmitchel/data/GWAS_data/variants.tsv",
                       sep="\t", header=TRUE)

# peak 1 "21:45624136:C:T"
# peak 2 "21:45660835:C:A"
v1 <- "21:45624136:C:T"
v2 <- "21:45660835:C:A"

ndx_get <- which(variants$variant==v1)
variants[ndx_get,] # rs8126495

ndx_get <- which(variants$variant==v2)
variants[ndx_get,] # rs2847224






### running command line SMR
# need to add MAF and beta sd for input to SMR
vcf_mat_sub <- vcf_mat[e_res$snps,]
af <- rowSums(vcf_mat_sub) / (ncol(vcf_mat_sub)*2)
e_res$af <- af
pr_res <- pr_res[rownames(e_res),]
pr_res$af <- af

a1 <- sapply(rownames(e_res),function(x){
  strsplit(x,split=':')[[1]][[3]]
})
a2 <- sapply(rownames(e_res),function(x){
  strsplit(x,split=':')[[1]][[4]]
})

e_res$A1 <- a1
e_res$A2 <- a2
pr_res$A1 <- a1
pr_res$A2 <- a2

# need to calculate se
e_res$se <- e_res$beta / e_res$statistic
pr_res$se <- pr_res$beta / pr_res$statistic

# add n stats
e_res$n <- ncol(pb)
pr_res$n <- ncol(pb)

# Columns needed for SMR: SNP    A1  A2  freq    b   se  p   n
e_res <- e_res[,c('snps','A1','A2','af','beta','se','pvalue','n')]
pr_res <- pr_res[,c('snps','A1','A2','af','beta','se','pvalue','n')]
colnames(e_res) <- c('SNP','A1','A2','freq','b','se','p','n')
colnames(pr_res) <- c('SNP','A1','A2','freq','b','se','p','n')
pr_res <- pr_res[,c('SNP','A2','A1','freq','b','se','p','n')]
colnames(pr_res) <- c('SNP','A1','A2','freq','b','se','p','n')
# write.table(pr_res,file='/home/jmitchel/data/lupus_data/for_smr/pr_res.ma',sep = ' ',row.names = FALSE,quote = FALSE)

## eQTL data needs to have slightly different format
# Columns are chromosome, SNP, genetic distance (can be any arbitary value since it will not be used in the SMR analysis),
# basepair position, the effect (coded) allele, the other allele and frequency of the effect allele
e_res$chromosome <- sapply(rownames(e_res),function(x){
  strsplit(x,split=':')[[1]][[1]]
})
e_res$pos <- sapply(rownames(e_res),function(x){
  strsplit(x,split=':')[[1]][[2]]
})
e_res$genetic_distance <- 0
e_res_sub <- e_res[,c('chromosome','SNP','genetic_distance','pos','A2','A1','freq')]
colnames(e_res_sub) <- NULL
# write.table(e_res_sub,file='/home/jmitchel/data/lupus_data/for_smr/e_res.esi',sep = ' ',row.names = FALSE,quote = FALSE)

gene_pos_dat <- as.data.frame(t(c(21,'ICOSLG',0,45642874,'ICOSLG','-')))
colnames(gene_pos_dat) <- NULL
# write.table(gene_pos_dat,file='/home/jmitchel/data/lupus_data/for_smr/e_res.epi',sep = ' ',row.names = FALSE,quote = FALSE)

# SNP    gene    beta    t-stat  p-value FDR
e_res2 <- me[["all"]][["eqtls"]]
rownames(e_res2) <- e_res2$snps
e_res2 <- e_res2[rownames(e_res),]
e_res2 <- e_res2[,c('snps','gene','beta','statistic','pvalue','FDR')]
colnames(e_res2) <- c('SNP','gene','beta','t-stat','p-value','FDR')
# write.table(e_res2,file='/home/jmitchel/data/lupus_data/for_smr/e_res_ss.txt',sep = ' ',row.names = FALSE,quote = FALSE)

# run in command line
# ./smr --eqtl-summary /home/jmitchel/data/lupus_data/for_smr/e_res_ss.txt --matrix-eqtl-format --make-besd --out /home/jmitchel/data/lupus_data/for_smr/make_besd/e_res
# then, moved the besd file to /for/smr/
# ./smr --bfile /home/jmitchel/data/lupus_data/for_smr/plink --gwas-summary /home/jmitchel/data/lupus_data/for_smr/pr_res.ma --beqtl-summary /home/jmitchel/data/lupus_data/for_smr/e_res --peqtl-smr .05 --out /home/jmitchel/data/lupus_data/for_smr/smr_out --thread-num 1

## now trying this but with subsets of SNPs for each peak
snps_pos_peak1 <- c(45590000,45640000)
snps_pos_peak2 <- c(45645000,45700000)

ndx_keep <- which(e_res$pos>snps_pos_peak1[1] & e_res$pos<snps_pos_peak1[2])
pr_res_p1 <- pr_res[ndx_keep,]
e_res_sub_p1 <- e_res_sub[ndx_keep,]
e_res2_p1 <- e_res2[ndx_keep,]
# write.table(pr_res_p1,file='/home/jmitchel/data/lupus_data/for_smr/pr_res_p1.ma',sep = ' ',row.names = FALSE,quote = FALSE)
# write.table(e_res_sub_p1,file='/home/jmitchel/data/lupus_data/for_smr/e_res_p1.esi',sep = ' ',row.names = FALSE,quote = FALSE)
# write.table(e_res2_p1,file='/home/jmitchel/data/lupus_data/for_smr/e_res_ss_p1.txt',sep = ' ',row.names = FALSE,quote = FALSE)

# ./smr --eqtl-summary /home/jmitchel/data/lupus_data/for_smr/e_res_ss_p1.txt --matrix-eqtl-format --make-besd --out /home/jmitchel/data/lupus_data/for_smr/make_besd/e_res_p1
# then, moved the besd file to /for/smr/. Also needed to make copy of e_res.epi to e_res_p1.epi
# ./smr --bfile /home/jmitchel/data/lupus_data/for_smr/plink --gwas-summary /home/jmitchel/data/lupus_data/for_smr/pr_res_p1.ma --beqtl-summary /home/jmitchel/data/lupus_data/for_smr/e_res_p1 --peqtl-smr .05 --out /home/jmitchel/data/lupus_data/for_smr/smr_out_p1 --thread-num 1


ndx_keep <- which(e_res$pos>snps_pos_peak2[1] & e_res$pos<snps_pos_peak2[2])
pr_res_p2 <- pr_res[ndx_keep,]
e_res_sub_p2 <- e_res_sub[ndx_keep,]
e_res2_p2 <- e_res2[ndx_keep,]
# write.table(pr_res_p2,file='/home/jmitchel/data/lupus_data/for_smr/pr_res_p2.ma',sep = ' ',row.names = FALSE,quote = FALSE)
# write.table(e_res_sub_p2,file='/home/jmitchel/data/lupus_data/for_smr/e_res_p2.esi',sep = ' ',row.names = FALSE,quote = FALSE)
# write.table(e_res2_p2,file='/home/jmitchel/data/lupus_data/for_smr/e_res_ss_p2.txt',sep = ' ',row.names = FALSE,quote = FALSE)

# ./smr --eqtl-summary /home/jmitchel/data/lupus_data/for_smr/e_res_ss_p2.txt --matrix-eqtl-format --make-besd --out /home/jmitchel/data/lupus_data/for_smr/make_besd/e_res_p2
# then, moved the besd file to /for/smr/. Also needed to make copy of e_res.epi to e_res_p1.epi
# ./smr --bfile /home/jmitchel/data/lupus_data/for_smr/plink --gwas-summary /home/jmitchel/data/lupus_data/for_smr/pr_res_p2.ma --beqtl-summary /home/jmitchel/data/lupus_data/for_smr/e_res_p2 --peqtl-smr .05 --out /home/jmitchel/data/lupus_data/for_smr/smr_out_p2 --thread-num 1


### for plotting results
image.scale <- function(z, zlim, col = heat.colors(12),
                        breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
  if(!missing(breaks)){
    if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
  }
  if(missing(breaks) & !missing(zlim)){
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  if(missing(breaks) & missing(zlim)){
    zlim <- range(z, na.rm=TRUE)
    zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
    zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
    breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
  }
  poly <- vector(mode="list", length(col))
  for(i in seq(poly)){
    poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
  }
  xaxt <- ifelse(horiz, "s", "n")
  yaxt <- ifelse(horiz, "n", "s")
  if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
  if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
  if(missing(xlim)) xlim=XLIM
  if(missing(ylim)) ylim=YLIM
  plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt, xaxs="i", yaxs="i", ...)
  for(i in seq(poly)){
    if(horiz){
      polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
    }
    if(!horiz){
      polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
    }
  }
}


# chromosome code, start of gene, end of gene and gene ID
gene_pos_dat2 <- as.data.frame(t(c(21,45642874,45660826,'ICOSLG')))
colnames(gene_pos_dat2) <- NULL
# write.table(gene_pos_dat2,file='/home/jmitchel/data/lupus_data/for_smr/glist-hg19',sep = ' ',row.names = FALSE,quote = FALSE)

# ./smr --bfile /home/jmitchel/data/lupus_data/for_smr/plink --gwas-summary /home/jmitchel/data/lupus_data/for_smr/pr_res_p1.ma --beqtl-summary /home/jmitchel/data/lupus_data/for_smr/e_res_p1 --peqtl-smr .001 --ld-upper-limit .95 --out /home/jmitchel/data/lupus_data/for_smr/smr_plot_p1 --plot --probe ICOSLG --probe-wind 500 --gene-list /home/jmitchel/data/lupus_data/for_smr/glist-hg19
SMRData = ReadSMRData("/home/jmitchel/data/lupus_data/for_smr/plot/smr_plot_p1.ICOSLG.txt")

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/smr_p1.pdf", useDingbats = FALSE,
#     width = 7, height = 5)
SMREffectPlot(data=SMRData, trait_name="aSLE_F2")
dev.off()


# ./smr --bfile /home/jmitchel/data/lupus_data/for_smr/plink --gwas-summary /home/jmitchel/data/lupus_data/for_smr/pr_res_p2.ma --beqtl-summary /home/jmitchel/data/lupus_data/for_smr/e_res_p2 --peqtl-smr .001 --ld-upper-limit .95 --out /home/jmitchel/data/lupus_data/for_smr/smr_plot_p2 --plot --probe ICOSLG --probe-wind 500 --gene-list /home/jmitchel/data/lupus_data/for_smr/glist-hg19
SMRData = ReadSMRData("/home/jmitchel/data/lupus_data/for_smr/plot/smr_plot_p2.ICOSLG.txt")

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/smr_p2.pdf", useDingbats = FALSE,
#     width = 7, height = 6)
SMREffectPlot(data=SMRData, trait_name="aSLE_F2")
dev.off()


# ./smr --bfile /home/jmitchel/data/lupus_data/for_smr/plink --gwas-summary /home/jmitchel/data/lupus_data/for_smr/pr_res.ma --beqtl-summary /home/jmitchel/data/lupus_data/for_smr/e_res --peqtl-smr .001 --ld-upper-limit .95 --ld-lower-limit 0 --out /home/jmitchel/data/lupus_data/for_smr/smr_plot --plot --probe ICOSLG --probe-wind 500 --gene-list /home/jmitchel/data/lupus_data/for_smr/glist-hg19
SMRData = ReadSMRData("/home/jmitchel/data/lupus_data/for_smr/plot/smr_plot.ICOSLG.txt")

# pdf(file = "/home/jmitchel/figures/scITD_revision_figs3/smr_full.pdf", useDingbats = FALSE,
#     width = 7, height = 5)
SMREffectPlot(data=SMRData, trait_name="aSLE_F2")
dev.off()




