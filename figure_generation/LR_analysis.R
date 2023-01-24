
library(Seurat)
library(spqn) 
library(RColorBrewer)
library(iTALK)
library(scITD)

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

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Ethnicity'),
                                        stat_use='pval')

## plot donor score
pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# just to check that everything is as expected up to this point
pbmc_container$plots$donor_matrix


##### running LR analysis
# using cellchat database
lr_pairs <- read.csv(file='/home/jmitchel/data/LR_datasets/Human-2020-Jin-LR-pairs.csv')
lr_pairs <- lr_pairs[,c('ligand','interaction_name')]
lr_pairs$receptor <- sapply(lr_pairs$interaction_name,function(x) {
  rname <- regmatches(x, regexpr("_", x), invert = TRUE)[[1]][[2]]
  return(rname)
})
lr_pairs$interaction_name <- NULL

# infer active LR interactions
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='pool')
sft_thresh <- c(12,14,12,10,12,9,12)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

lr_hmap <- compute_LR_interact(pbmc_container, lr_pairs, sig_thresh=.00000000005,
                               percentile_exp_rec=0.85, add_ld_fact_sig=TRUE)

### Figure S3B
# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_new_lr8.pdf", useDingbats = FALSE,
#     width = 6, height = 7)
lr_hmap
# dev.off()


lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,mod_ct='Th',mod=5,lig_ct='cMono',lig='ICOSLG')

### Figure 3E
# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_ICOSLG_trio2.pdf", useDingbats = FALSE,
#     width = 6, height = 5)
lig_mod_fact
# dev.off()


# getting GO enrichment HMAPs for modules
ctypes <- c('Th')
modules <- c(5)

# using more stringent p-val threshold for figure plot
# mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.002, db_use=c('GO'),max_plt_pval=.002,h_w=c(7,3))
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.005, 
                                 db_use=c('GO','BioCarta'),max_plt_pval=.005,h_w=c(14,2))

### Figure 3F
# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_ICOSLG_gsets2.pdf", useDingbats = FALSE,
#     width = 5, height = 8)
mod_enr
# dev.off()


#### computing enrichment of modules with the nichenet scores
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

## ICOSLG target module first
test = ligand_target_matrix[,'ICOSLG']
t_mod = pbmc_container[["module_genes"]][["Th"]]
mymod <- c(5)
g_in_mod <- names(t_mod)[t_mod%in%mymod]
g_not_mod <- names(t_mod)[!(t_mod%in%mymod)]
tmp <- cbind.data.frame(c(g_in_mod,g_not_mod),
                        c(rep('Th_m5',length(g_in_mod)),rep('other',length(g_not_mod))))
colnames(tmp) <- c('gn','in_mod')
tmp$in_mod <- factor(tmp$in_mod,levels=c('Th_m5','other'))
target_scores <- test[tmp$gn]
tmp$target_scores <- target_scores
tmp <- tmp[which(!is.na(target_scores)),]
p <- ggplot(tmp,aes(x=as.factor(in_mod),y=target_scores)) +
  geom_boxplot(notch=TRUE) +
  ylab('ICOSLG NicheNet regulatory potential') +
  xlab('') +
  theme_bw()
test_res <- wilcox.test(target_scores~in_mod,data=tmp)
print(test_res$p.value)

### Figure S3C
# pdf(file = "/home/jmitchel/figures/for_paper_v2/ICOSLG_NicheNet_enr2.pdf", useDingbats = FALSE,
#     width = 4, height = 3.5)
p
# dev.off()



## showing THBS1 LR interaction relevant for factor 3
lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=3,mod_ct='Th',mod=9,lig_ct='cDC',lig='THBS1')

### Figure 4E
# pdf(file = "/home/jmitchel/figures/for_paper_v2/sle_THBS1_trio.pdf", useDingbats = FALSE,
#     width = 6, height = 5)
lig_mod_fact
# dev.off()


## extracting the pvalues for these two examples
# running code within the compute_LR_interact fn
sig_thresh=.0000001
myres_mat <- pbmc_container$lr_res # at 285
container=pbmc_container

pbmc_container$lr_res['ICOSLG_cMono_ICOS','Th_m5']
which(rownames(myres_mat)=='ICOSLG_cMono_ICOS')
10**(-fact_res2[47,2])

pbmc_container$lr_res['THBS1_cDC_CD47','Th_m9']
which(rownames(myres_mat)=='THBS1_cDC_CD47')
10**(-fact_res2[35,3])


