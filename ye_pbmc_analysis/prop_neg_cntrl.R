
library("lmtest")

## assuming already built pbmc_container and run ICA up to this point
# also have subclusters loaded
pbmc_container$subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')

# loop through major cell types
ctypes_use <- pbmc_container$experiment_params$ctypes_use
print(ctypes_use)
all_res <- c(.8,.6,.6,.6,.5,.5,.6)
names(all_res) <- ctypes_use
mypvals <- c()
ctype_vals <- c()
subtype_vals <- c()
for (ctype in ctypes_use) {
  # get subtype proportions
  res <- all_res[ctype]
  resolution_name <- paste0('res:',as.character(res))
  subclusts <- pbmc_container$subclusters[[ctype]][[resolution_name]]
  
  # append large cell type name to subclusters
  subclusts <- sapply(subclusts,function(x){paste0(ctype,'_',x)})
  
  # limit cells in subclusts to those that we actually have scores for
  donor_scores <- pbmc_container$tucker_results[[1]]
  donor_vec <- pbmc_container$scMinimal_full$metadata[names(subclusts),'donors']
  subclusts <- subclusts[donor_vec %in% rownames(donor_scores)]
  
  # make subtype association plot
  subclusts_num <- sapply(subclusts,function(x){as.numeric(strsplit(x,split="_")[[1]][[2]])})
  scMinimal <- pbmc_container$scMinimal_ctype[[ctype]]
  sub_meta_tmp <- scMinimal$metadata[names(subclusts),]
  
  # get donor proportions of subclusters
  donor_props <- compute_donor_props(subclusts_num,sub_meta_tmp)
  
  # get batch each donor is in
  meta <- pbmc_container$scMinimal_ctype[[1]]$metadata[,c('donors','pool','Status')]
  meta <- unique(meta)
  rownames(meta) <- meta$donors
  meta$donors <- NULL
  
  ### to optionally shuffle to batches
  meta$pool <- sample(meta$pool)
  ###
  
  # # trying by getting balances first
  # saved_names <- rownames(donor_props)
  # donor_props <- coda.base::coordinates(donor_props)
  # rownames(donor_props) <- saved_names
  
  # loop through subtypes
  for (j in 1:ncol(donor_props)) {
    # first test if healthy donors have the cell type in sufficient proportions
    # tmp <- as.data.frame(cbind(donor_props[rownames(meta),j],as.character(meta[,2])))
    # colnames(tmp) <- c('prop','Status')
    # tmp$prop <- as.numeric(tmp$prop)
    # tmp$Status <- as.factor(tmp$Status)
    # tmp_sub <- tmp[tmp$Status=='Healthy',]
    # frac_with_subtype <- sum(tmp_sub$prop>.15)/nrow(tmp_sub)
    # if (frac_with_subtype < .8) {
    #   next
    # } 
    
    # test whether batch variable significantly predicts subtype proportion
    tmp <- as.data.frame(cbind(donor_props[rownames(meta),j],as.character(meta[,1])))
    colnames(tmp) <- c('prop','batch')
    tmp$prop <- as.numeric(tmp$prop)
    tmp$batch <- as.factor(tmp$batch)
    
    # # using linear model
    # lmres <- lm(prop~batch,data=tmp) #should use the beta model instead
    # lmres <- summary(lmres)
    # pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
    
    # breg <- betareg::betareg(prop~batch, data = tmp)
    # tmp <- summary(breg)
    # reg_stat <- tmp$coefficients$mean['dscore','Pr(>|z|)']
    
    ## tyring liklihood ratio test for beta regression overall significance
    breg1 <- betareg::betareg(prop~1, data = tmp)
    breg2 <- betareg::betareg(prop~batch, data = tmp)
    pval <- lrtest(breg1, breg2)$`Pr(>Chisq)`[2]

    # store pvalue
    mypvals <- c(mypvals,pval)
    ctype_vals <- c(ctype_vals,ctype)
    subtype_vals <- c(subtype_vals,paste0(ctype,'_',j))
  }
}

tmp <- as.data.frame(cbind(mypvals,ctype_vals,subtype_vals))
colnames(tmp) <- c('pvals','ctypes','subtype')
tmp$pvals <- as.numeric(tmp$pvals)
tmp$pvals <- p.adjust(tmp$pvals,method='fdr')
tmp$pvals <- -log10(tmp$pvals)

# # see if it's the really small subclusters that have super high associations
# ctype <- 'NK'
# res <- all_res[ctype]
# resolution_name <- paste0('res:',as.character(res))
# subclusts <- pbmc_container$subclusters[[ctype]][[resolution_name]]
# subc_select <- subclusts[subclusts==4]

# get fractional sizes each population
frac_pops <- c()
subtype_track <- c()
for (ctype in ctypes_use) {
  res <- all_res[ctype]
  resolution_name <- paste0('res:',as.character(res))
  subclusts <- pbmc_container$subclusters[[ctype]][[resolution_name]]
  subc_dif <- unique(subclusts)
  for (sc in 1:length(subc_dif)) {
    subc_select <- subclusts[subclusts==sc]
    frac_pop <- round(length(subc_select) / length(subclusts),3)
    frac_pops <- c(frac_pops,frac_pop)
    subtype_track <- c(subtype_track,paste0(ctype,"_",sc))
  }
}
names(frac_pops) <- subtype_track
frac_pops <- frac_pops[tmp$subtype]

# get relative entropy for each subcluster
# made compute donor props return counts instead
xt <- donor_props
# need to collapse these by pool
pools <- as.character(unique(meta$pool))
pool_counts <- as.data.frame(matrix(0,ncol=ncol(xt),nrow=length(pools)))
rownames(pool_counts) <- pools
for (i in 1:nrow(pool_counts)) {
  mypool <- rownames(pool_counts)[i]
  meta_sub <- meta[meta$pool==mypool,,drop=FALSE]
  d_select <- rownames(meta_sub)
  for (d in d_select) {
    pool_counts[i,] <- pool_counts[i,] + xt[d,]
  }
}
xt <- pool_counts
n.samples <- nrow(xt)
ne <- 1-apply(xt, 2, entropy::KL.empirical, y2=rowSums(xt), unit=c('log2')) / log2(n.samples) # relative entropy

# count number of donors with at least 10% of NK cells as NK4
# made a tmp with prop and Status
colnames(tmp) <- c('prop','Status')
tmp_sub <- tmp[tmp$Status=='Healthy',]
sum(tmp_sub$prop>.1)
dim(tmp_sub)

# plot results
ggplot(tmp, aes(x=as.factor(subtype),y=pvals)) +
  geom_col() +
  geom_text(aes(label = frac_pops), vjust = -0.2) +
  geom_hline(yintercept=-log10(.01),linetype="dashed",color = "red") +
  xlab('Cell subtype') +
  ylab('-log10(padj)')

ggplot(tmp, aes(x=as.factor(subtype),y=pvals)) +
  geom_col() +
  geom_hline(yintercept=-log10(.01),linetype="dashed",color = "red") +
  xlab('Cell subtype') +
  ylab('-log10(padj)') +
  ggtitle("Shuffled Batch Structure") +
  theme(plot.title = element_text(hjust = 0.5))


## need to see if those first two batches have more managed donors
meta <- pbmc_container$scMinimal_ctype[[1]]$metadata[,c('donors','pool','Status')]
meta <- unique(meta)
rownames(meta) <- meta$donors
meta$donors <- NULL
meta_sub <- meta[meta$pool=='dmx_YE_8-16',]
table(meta_sub$Status)

# could also manually check to see if pool is associated with F5
dsc <- pbmc_container$tucker_results[[1]][,5]
tmp <- as.data.frame(cbind(dsc,as.character(meta[names(dsc),1])))
colnames(tmp) <- c('dsc','pool')
tmp$pool <- as.factor(tmp$pool)
tmp$dsc <- as.numeric(tmp$dsc)
head(tmp)
class(tmp$pool)
class(tmp$dsc)
lmres <- lm(dsc~pool,data=tmp)
summary(lmres)
# confirms that there is no association with the factor


# I wonder if pool if confounded with study because this could explain the systematic
# differences in proportions between pools
meta_sub <- unique(pbmc@meta.data[,c("batch_cov","Study",'Status','ind_cov_batch_cov','Processing_Cohort')])
unique(meta_sub$batch_cov)
sub_sub <- meta_sub[meta_sub$batch_cov=='dmx_YS-JY-22_pool6',]
rownames(sub_sub) <- NULL
sub_sub


















