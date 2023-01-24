library(Seurat)
library(cowplot)
library(scITD)

# load up the covid dataset: see preprocessing/covid_uk_preprocessing.ipynb
# for code used to generate this object
pbmc_covid <- readRDS(file="/home/jmitchel/data/covid_data_uk/haniffa21_subset_no_lps.rds")

# collapsing cell subtypes and converting cell type names 
new_names <- sapply(as.character(pbmc_covid@meta.data$full_clustering), function(x){
  if (x=='B_exhausted' || x=='B_immature' || x=='B_naive' || x=='B_non-switched_memory' || x=='B_switched_memory') {
    return('B')
  } else if (x=='CD14_mono' || x=='CD83_CD14_mono') {
    return('cMono')
  } else if (x=='C1_CD16_mono' || x=='CD16_mono') {
    return('ncMono')
  } else if (x=='CD4.CM' || x=='CD4.EM' || x=='CD4.IL22' || x=='CD4.Naive' || x=='CD4.Prolif' || x=='CD4.Tfh' || x=='CD4.Th1') {
    return('Th')
  } else if (x=='CD8.EM' || x=='CD8.Naive' || x=='CD8.Prolif' || x=='CD8.TE') {
    return('Tc')
  } else if (x=='NK_16hi' || x=='NK_56hi' || x=='NK_prolif') {
    return('NK')
  } else {
    return(x)
  }
})

names(new_names) <- NULL
pbmc_covid@meta.data$initial_clustering <- factor(new_names,levels=unique(new_names))

# remove plasmablast-like B cells
cells_keep1 <- rownames(pbmc_covid@reductions[["umap"]]@cell.embeddings)[pbmc_covid@reductions[["umap"]]@cell.embeddings[,1]>-1]
cells_keep2 <- rownames(pbmc_covid@reductions[["umap"]]@cell.embeddings)[pbmc_covid@reductions[["umap"]]@cell.embeddings[,2]>.5]
cells_keep <- unique(c(cells_keep1,cells_keep2))

pbmc_covid <- subset(pbmc_covid,cells=cells_keep)



# remove antibody "genes"
g_ndx_rem <- which(pbmc_covid@assays[["RNA"]]@meta.features=='Antibody Capture')
g_ndx_keep <- c(1:nrow(pbmc_covid@assays[["raw"]]@counts))[!(c(1:nrow(pbmc_covid@assays[["raw"]]@counts)) %in% g_ndx_rem)]
dim(pbmc_covid)
pbmc_covid <- CreateSeuratObject(pbmc_covid@assays[["raw"]]@counts[g_ndx_keep,], project = "SeuratProject", assay = "raw",
                                 meta.data = pbmc_covid@meta.data)
dim(pbmc_covid)


# starting analysis
param_list <- initialize_params(ctypes_use = c("B","Tc","Th","NK","cMono"), ncores = 30, rand_seed = 10)

pbmc_container <- make_new_container(count_data = pbmc_covid@assays$raw@counts,
                                           meta_data = pbmc_covid@meta.data,
                                           params=param_list,
                                           metadata_cols=c('patient_id',
                                                           "Sex",
                                                           "Age_interval",
                                                           "initial_clustering",
                                                           "Status",
                                                           "Status_on_day_collection",
                                                           "Status_on_day_collection_summary",
                                                           'Days_from_onset',
                                                           "Site",
                                                           "Smoker",
                                                           'Worst_Clinical_Status',
                                                           'Outcome',
                                                           'Swab_result',
                                                           'time_after_LPS'),
                                           metadata_col_nm=c('donors',
                                                             'sex',
                                                             'age',
                                                             'ctypes',
                                                             'status',
                                                             'status_on_day_collection',
                                                             'status_on_day_collection_summary',
                                                             'days_from_onset',
                                                             'site',
                                                             "smoker",
                                                             'worst_Clinical_Status',
                                                             'outcome',
                                                             'swab_result',
                                                             'time_after_LPS'))

pbmc_container <- form_tensor(pbmc_container, donor_min_cells=2, 
                                    norm_method='trim', scale_factor=10000,
                                    vargenes_method='norm_var_pvals', vargenes_thresh=.0001,
                                    batch_var = 'site',
                                    scale_var = TRUE, var_scale_power = .5) 

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(9,38), 
                                       tucker_type = 'regular', rotation_type = 'hybrid') # best with vargenes_thresh=.0001 also .01

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','age',
                                                                               'status_on_day_collection_summary',
                                                                               'site',"smoker"),stat_use='pval')
# plot donor scores by status
pbmc_container <- plot_donor_matrix(pbmc_container,
                                          show_donor_ids = FALSE,
                                          add_meta_associations='pval',h_w=c(10,8))


pbmc_container$plots$donor_matrix


## flipping signs of f3 and f5
pbmc_container$tucker_results[[1]][,3] <- pbmc_container$tucker_results[[1]][,3] * -1
pbmc_container$tucker_results[[2]][3,] <- pbmc_container$tucker_results[[2]][3,] * -1
pbmc_container$tucker_results[[1]][,5] <- pbmc_container$tucker_results[[1]][,5] * -1
pbmc_container$tucker_results[[2]][5,] <- pbmc_container$tucker_results[[2]][5,] * -1




# load up the covid validation dataset: see preprocessing/covid_ye_preprocessing.R
# for code used to generate this object
pbmc <- readRDS(file='/home/jmitchel/data/covid_data_ye/covid_subset.rds')
# making pool variable not have batches without representation
pbmc@meta.data$pool <- factor(pbmc@meta.data$pool,levels=unique(pbmc@meta.data$pool))

# changing cell type names to match others used
new_names <- sapply(as.character(pbmc@meta.data$ct1), function(x){
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
pbmc@meta.data$ct1 <- factor(new_names,levels=unique(new_names))

param_list <- initialize_params(ctypes_use = c("B","NK","Th","Tc",
                                               "cMono"),
                                ncores = 30, rand_seed = 10)
pbmc_ye <- make_new_container(seurat_obj=pbmc,
                              params=param_list,
                              metadata_cols=c('donor',
                                              "COVID_status",
                                              'COVID_severity',
                                              "COVID_severity_merged",
                                              "ct1",
                                              "sex",
                                              "pool",
                                              "ethnicity",
                                              "respiratory_support_D0",
                                              'onset_to_D0_days',
                                              'intubated_days',
                                              'admission_to_discharge',
                                              'death',
                                              'pulmonary_infection',
                                              'non_pulmonary_infection',
                                              'WBC_count1',
                                              'WBC_count2',
                                              'WBC_count3',
                                              'Monocyte_count',
                                              'respiratory_support',
                                              'NIH_clinical',
                                              'NIH_ordinal',
                                              'admission_level'),
                              metadata_col_nm=c('donors',
                                                'covid_status',
                                                'covid_severity',
                                                'covid_severity_merged',
                                                'ctypes',
                                                'sex',
                                                'pool',
                                                'Ethnicity',
                                                'respiratory_support_D0',
                                                'onset_to_D0_days',
                                                'intubated_days',
                                                'admission_to_discharge',
                                                'death',
                                                'pulmonary_infection',
                                                'non_pulmonary_infection',
                                                'WBC_count1',
                                                'WBC_count2',
                                                'WBC_count3',
                                                'Monocyte_count',
                                                'respiratory_support',
                                                'NIH_clinical',
                                                'NIH_ordinal',
                                                'admission_level'))

pbmc_ye <- form_tensor(pbmc_ye, donor_min_cells=2, 
                       norm_method='trim', scale_factor=10000,
                       vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                       scale_var = TRUE, var_scale_power = .5,
                       batch_var='pool')

pbmc_ye <- run_tucker_ica(pbmc_ye, ranks=c(10,30), # tested well with comparison below
                          tucker_type = 'regular', rotation_type = 'hybrid')

# get factor-meta data associations
pbmc_ye <- get_meta_associations(pbmc_ye,vars_test=c('sex','pool',
                                                     'covid_status',
                                                     'covid_severity',
                                                     'covid_severity_merged',
                                                     'death'),
                                 stat_use='pval')

# plot donor scores by status
pbmc_ye <- plot_donor_matrix(pbmc_ye, meta_vars=c('sex'),
                             cluster_by_meta = 'sex',
                             show_donor_ids = FALSE,
                             add_meta_associations='pval')

pbmc_ye$plots$donor_matrix

pbmc_ye <- get_lm_pvals(pbmc_ye)

# getting count of covid and healthy patients
pbmc_ye <- get_donor_meta(pbmc_ye,'covid_status')
table(pbmc_ye$donor_metadata$covid_status)

# getting # cells used from the major clusters
d_kept <- rownames(pbmc_ye$scMinimal_ctype[[1]]$pseudobulk)
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$donor %in% d_kept]
pbmc <- subset(pbmc,cells=cells_keep)

ct_kept <- pbmc_ye$experiment_params$ctypes_use
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$ct1 %in% ct_kept]
pbmc <- subset(pbmc,cells=cells_keep)
dim(pbmc)


##### get median number cells per donor and total number of cells used
d_keep <- rownames(pbmc_ye$scMinimal_ctype[[1]]$pseudobulk)
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$donor %in% d_keep]
pbmc_sub <- subset(pbmc,cells=cells_keep)
ctypes <- c("B","Tc","Th","NK","cMono")
cells_keep <- rownames(pbmc_sub@meta.data)[pbmc_sub@meta.data$ct1 %in% ctypes]
pbmc_sub2 <- subset(pbmc_sub,cells=cells_keep)
pbmc_sub2@meta.data$donor <- factor(as.character(pbmc_sub2@meta.data$donor),levels=unique(as.character(pbmc_sub2@meta.data$donor)))
for (ctype in ctypes) {
  tmp <- pbmc_sub2@meta.data[pbmc@meta.data$ct1==ctype,]
  print(ctype)
  print(median(table(tmp$donor)))
}





## projecting UK factors onto Jimmie covid data
pbmc_ye <- project_new_data(pbmc_ye,pbmc_container)
cor(pbmc_ye$projected_scores,pbmc_ye$tucker_results[[1]])

# now plotting projected scores versus metadata
dsc <- pbmc_ye$projected_scores[,2,drop=FALSE]
pbmc_ye <- get_donor_meta(pbmc_ye,additional_meta = c('covid_status',
                                                      'covid_severity_merged'),only_analyzed = TRUE)
tmp <- cbind.data.frame(dsc,pbmc_ye$donor_metadata[rownames(dsc),'covid_severity_merged'],
                        pbmc_ye$donor_metadata[rownames(dsc),'covid_status'])
colnames(tmp) <- c('dscore','severity','status')
levels(tmp$severity)
levels(tmp$severity) <- c('Critical','Healthy','Moderate','Severe','Critical_anti-IFN','Negative')
tmp$severity <- factor(tmp$severity,levels=c('Healthy','Moderate','Severe','Critical','Critical_anti-IFN'))

mycol <- RColorBrewer::brewer.pal(n = 7, name = "Dark2")

tmp2 <- tmp
tmp2$status <- rep('violin',nrow(tmp2))
p <- ggplot(tmp,aes(x=severity,y=dscore,fill=status)) +
  geom_violin(data=tmp2) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=1.85, binwidth = .01) +
  ylab('Projected Factor 2 Scores') +
  xlab('Severity on collection day') +
  coord_flip() +
  scale_fill_manual(values=c(mycol[2], mycol[5], 'light gray')) +
  theme_bw() +
  theme(axis.text=element_text(size=24),
        axis.title=element_text(size=26))

### Figure 7C
# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_f2_proj.pdf", useDingbats = FALSE,
#     width = 12.5, height = 8.5)
p
# dev.off()

## calculate significance of the association
dsc <- pbmc_ye$projected_scores[,2,drop=FALSE]
tmp <- cbind.data.frame(dsc,
                        pbmc_ye$donor_metadata[rownames(dsc),'covid_severity_merged'])
colnames(tmp) <- c('proj_dsc','status')
lmres <- summary(lm(proj_dsc~status,data=tmp))
lmres
pv1 <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)

pv2 <- pbmc_container[["meta_associations"]]['status_on_day_collection_summary',2]

print(pv1)
print(pv2)

## combining severity association p-values for meta analysis
library(metap)
sumlog(c(pv1,pv2))
# p =  4.342494e-17 








##### running LR analysis - first with main covid dataset
# iTalk db
library(iTALK)
lr_pairs <- database[,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')]
colnames(lr_pairs) <- c('ligand','receptor')
lr_pairs <- unique(lr_pairs)

# infer active LR interactions
pbmc_container <- prep_LR_interact(pbmc_container, lr_pairs, norm_method='trim', scale_factor=10000,
                                   var_scale_power=.5, batch_var='site')
sft_thresh <- c(12,12,12,12,10)
pbmc_container <- get_gene_modules(pbmc_container,sft_thresh)

lr_hmap <- compute_LR_interact(pbmc_container, lr_pairs, sig_thresh=.0000001,
                               percentile_exp_rec=0.85, add_ld_fact_sig=TRUE)

lr_hmap


### looking at IL16 hit
lig_mod_fact <- plot_mod_and_lig(pbmc_container,factor_select=2,
                                 mod_ct='cMono',mod=14,lig_ct='Th',lig='IL16')
lig_mod_fact

### Figure 7E
# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_IL16_trio3.pdf", useDingbats = FALSE,
#     width = 5, height = 5.5)
lig_mod_fact 
# dev.off()


ctypes <- c('cMono')
modules <- c(14)
mod_enr <- plot_multi_module_enr(pbmc_container, ctypes, modules, sig_thresh=.009, 
                                 db_use=c('GO'),
                                 max_plt_pval=.009,h_w=c(12,6))

### Figure 7F
# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_IL16_enr.pdf", useDingbats = FALSE,
#     width = 6, height = 7)
mod_enr
# dev.off()

# getting p-value 
pbmc_container$lr_res['IL16_Th_CD4','cMono_m14']











## running the full lr analysis for jimmie's covid data
pbmc_ye <- prep_LR_interact(pbmc_ye, lr_pairs, norm_method='trim', scale_factor=10000,
                            var_scale_power=.5, batch_var='pool')
sft_thresh <- c(12,14,10,14,12)
pbmc_ye <- get_gene_modules(pbmc_ye,sft_thresh)

lr_hmap <- compute_LR_interact(pbmc_ye, lr_pairs, sig_thresh=.0001,
                               percentile_exp_rec=0.85, add_ld_fact_sig=TRUE)

lr_hmap

# make sure that dscores flipped to match those of Stephenson et al. dataset
pbmc_ye$tucker_results[[1]][,5] <- pbmc_ye$tucker_results[[1]][,5]*-1
lig_mod_fact <- plot_mod_and_lig(pbmc_ye,factor_select=5,
                                 mod_ct='cMono',mod=4,lig_ct='Th',lig='IL16')

### Figure S6C
# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_jimmie_IL16_trio.pdf", useDingbats = FALSE,
#     width = 5, height = 5.5)
lig_mod_fact
# dev.off()

ctypes <- c('cMono')
modules <- c(4)
mod_enr <- plot_multi_module_enr(pbmc_ye, ctypes, modules, sig_thresh=.005, 
                                 db_use=c('GO'),
                                 max_plt_pval=.005,h_w=c(12,6))

### Figure S6D
# pdf(file = "/home/jmitchel/figures/for_paper_v2/covid_jimmie_IL16_GO.pdf", useDingbats = FALSE,
#     width = 6, height = 7)
mod_enr
# dev.off()

pbmc_ye$lr_res['IL16_Th_CD4','cMono_m4']












