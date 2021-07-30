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
                                                     "Processing_Cohort"),
                                     metadata_col_nm=c('donors',
                                                       'SLE_status',
                                                       'Status',
                                                       'ctypes',
                                                       'sex',
                                                       'Age',
                                                       'pool',
                                                       'processing'))


pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20, gene_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.05,
                              scale_var = TRUE, var_scale_power = 1.5)

pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(20,30,7),
                                 tucker_type = 'regular', rotation_type = 'ica')

# get factor-meta data associations
pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('pool','processing'))

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('pool','processing'),
                                    cluster_by_meta='pool',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='rsq')

# pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_dscores.pdf", useDingbats = FALSE,
#     width = 7, height = 6.5)
pbmc_container$plots$donor_matrix
dev.off()

pbmc_container <- get_all_lds_factor_plots(pbmc_container, use_sig_only=FALSE,
                                           nonsig_to_zero=FALSE,
                                           display_genes=FALSE,
                                           gene_callouts=FALSE,
                                           show_var_explained = FALSE)

# 2,4,6,15
pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_f2_lds.pdf", useDingbats = FALSE,
    width = 4, height = 5)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["2"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["2"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_f4_lds.pdf", useDingbats = FALSE,
    width = 4, height = 5)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["4"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["4"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_f6_lds.pdf", useDingbats = FALSE,
    width = 4, height = 5)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["6"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["6"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()

pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_f15_lds.pdf", useDingbats = FALSE,
    width = 4, height = 5)
draw(pbmc_container[["plots"]][["all_lds_plots"]][["15"]],
     annotation_legend_list = pbmc_container[["plots"]][["all_legends"]][["15"]],
     legend_grouping = "original",
     newpage=TRUE)
dev.off()



# get significant genes
pbmc_container <- run_jackstraw(pbmc_container, ranks=c(20,30,7), n_fibers=100, n_iter=500,
                                tucker_type='regular', rotation_type='ica')

# saveRDS(pbmc_container[["gene_score_associations"]],file='/home/jmitchel/data/lupus_data/lupus_no_combat_jackstraw.rds')
pbmc_container[["gene_score_associations"]] <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_no_combat_jackstraw.rds')

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=15, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.15, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)

pbmc_container <- plot_loadings_annot(pbmc_container, factor_select=2, use_sig_only=TRUE, nonsig_to_zero=TRUE, annot='none',
                                      pathways=NULL, sim_de_donor_group=NULL, sig_thresh=0.15, display_genes=FALSE,
                                      gene_callouts=FALSE, callout_n_gene_per_ctype=5, callout_ctypes=NULL, show_le_legend=FALSE,
                                      show_xlab=TRUE, show_var_explained=TRUE, reset_other_factor_plots=FALSE, draw_plot=TRUE)

# run gsea
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=15, method="hypergeometric", thresh=0.05,
                                      db_use=c("GO"))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=15,direc='up',thresh=.05)


pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=2, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=2,direc='up',thresh=.05)
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=1)

gsets <- c("GO_NEGATIVE_REGULATION_OF_BIOSYNTHETIC_PROCESS",
           "GO_NEGATIVE_REGULATION_OF_RNA_BIOSYNTHETIC_PROCESS",
           "GO_POSITIVE_REGULATION_OF_RNA_BIOSYNTHETIC_PROCESS",
           "GO_POSITIVE_REGULATION_OF_PROTEIN_METABOLIC_PROCESS",
           "GO_NEGATIVE_REGULATION_OF_PHOSPHORUS_METABOLIC_PROCESS",
           "GO_MICROTUBULE_CYTOSKELETON_ORGANIZATION",
           "GO_CHROMOSOME_ORGANIZATION",
           "GO_CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS",
           "GO_MUSCLE_STRUCTURE_DEVELOPMENT",
           "GO_NEUROGENESIS",
           "GO_POSITIVE_REGULATION_OF_CELL_DIFFERENTIATION",
           "GO_BIOLOGICAL_ADHESION",
           "GO_NON_CANONICAL_WNT_SIGNALING_PATHWAY",
           "GO_DETOXIFICATION",
           "GO_ESTABLISHMENT_OF_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM",
           "GO_OXIDATIVE_PHOSPHORYLATION",
           "GO_MORPHOGENESIS_OF_A_POLARIZED_EPITHELIUM")

gset_cmap <- rep('black',length(gsets))

names(gset_cmap) <- gsets

pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_f2_gsea.pdf", useDingbats = FALSE,
    width = 9, height = 4)
hm_list <- plot_select_sets(pbmc_container, 2, gsets, color_sets=gset_cmap, 
                            cl_rows=TRUE, h_w=c(8.15,9), myfontsize=7.15)
dev.off()




## looking into a factor to see how the genes are perturbed
f1_data <- get_one_factor(pbmc_container, factor_select=4)
dscores <- f1_data[[1]]
lds <- f1_data[[2]]

# get a gene with largest magnitude loadings
which(lds[,1]==lds[,1][order(lds[,1],decreasing=FALSE)][5])
lds['COX6A1P2',]

# plot expression of this gene
pb <- pbmc_container[["scMinimal_ctype"]][["B"]][["pseudobulk"]]
expr <- pb[,'RPSAP58']
tmp <- cbind.data.frame(expr[rownames(dscores)],dscores)
colnames(tmp) <- c('expr','dsc')
ggplot(tmp,aes(x=dsc,y=expr)) +
    geom_point()

expr <- pb[,c('RPS10','TPT1')]
tmp <- cbind.data.frame(expr[rownames(dscores),1],expr[rownames(dscores),2],dscores)
colnames(tmp) <- c('expr1','expr2','dsc')
ggplot(tmp,aes(x=expr1,y=expr2,color=dsc)) +
    geom_point()




## trying to run decontX
library(SingleCellExperiment)
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$batch_cov=='dmx_YS-JY-20_pool3']
pbmc_sub <- subset(pbmc,cells = cells_keep)
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
pbmc_sub <- as.SingleCellExperiment(pbmc_sub) # convert to single-cell experiment
saveRDS(pbmc_sub,file='/home/jmitchel/data/lupus_data/lupus_test_batch.rds')
# dc <- decontX(pbmc_sub)

# from decontX done locally for the above batch
eta <- readRDS('/home/jmitchel/decontX_eta2.RDS')

f_data <- get_one_factor(pbmc_container, factor_select=8)
f_data <- get_one_factor(pbmc_container, factor_select=15)
dscores <- f_data[[1]]
lds <- f_data[[2]]

head(lds)

eta <- eta[rownames(lds),colnames(lds)]

plot(lds[,1],eta[,1])
plot(lds[,7],log(eta[,7]))

cor(lds[,7],eta[,7],method='spearman')


# trying to see if there is a more clear trend just among the significant genes
sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=15, colnames(lds))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

# limit to just the genes in tmp_casted_num
sig_df <- sig_df[rownames(lds),colnames(lds)]

lds <- lds[sig_df[,3]<0.05,3,drop=FALSE]


plot(lds,log(eta[rownames(lds),3]))
cor(lds,log(eta[rownames(lds),3]),method='spearman')





# trying same thing but with soupx
library(SoupX)
toc = counts(pbmc_sub)
scNoDrops = SoupChannel(toc, toc, calcSoupProfile = FALSE)
# Calculate soup profile
soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))

sum(scNoDrops$soupProfile$est)

f_data <- get_one_factor(pbmc_container, factor_select=15)
f_data <- get_one_factor(pbmc_container, factor_select=7)
dscores <- f_data[[1]]
lds <- f_data[[2]]

soupProf <- soupProf[rownames(lds),]

plot(lds[,1],soupProf$est)
plot(lds[,1],log(soupProf$est))

cor(lds[,1],soupProf$est,method='spearman')



# get gene properties from biomart
library(biomaRt)
f_data <- get_one_factor(pbmc_container, factor_select=1)
dscores <- f_data[[1]]
lds <- f_data[[2]]

ensembl = useMart("ensembl",
                  dataset="hsapiens_gene_ensembl")


ens_id <- getBM(attributes=c('ensembl_gene_id',
                          'hgnc_symbol'),
             filters = 'hgnc_symbol',
             mart = ensembl,
             values = rownames(lds))

# ens_id <- getBM(attributes=c('percentage_gene_gc_content',
#                              'hgnc_symbol'),
#                 filters = 'hgnc_symbol',
#                 mart = ensembl,
#                 values = rownames(lds))


library(EDASeq)

gc_len <- getGeneLengthAndGCContent(ens_id[,1], 'hsa')

gen_dat <- cbind.data.frame(ens_id,gc_len)

# reduce to non-repeating genes only
# first get genes with multiple ensg ID
gen_dat <- gen_dat[!(duplicated(gen_dat[,2]) | duplicated(gen_dat[,2], fromLast=TRUE)),]


# now add on some other features
listAttributes(ensembl, page = "feature_page")
library(GenomicFeatures)
library(dplyr)
refSeq <- makeTxDbFromUCSC(genom="hg38",tablename="refGene")   
threeUTRs <- threeUTRsByTranscript(refSeq, use.names=TRUE)
length_threeUTRs <- width(ranges(threeUTRs))
the_lengths <- as.data.frame(length_threeUTRs)
the_lengths <- the_lengths %>% group_by(group, group_name) %>% summarise(sum(value))
the_lengths <- unique(the_lengths[,c("group_name", "sum(value)")])
colnames(the_lengths) <- c("RefSeq Transcript", "3' UTR Length")

ref_id <- getBM(attributes=c('refseq_mrna',
                             'hgnc_symbol'),
                filters = 'hgnc_symbol',
                mart = ensembl,
                values = rownames(lds))

# remove ones where no ref_id is present
ref_id <- ref_id[ref_id[,1]!='',]

# remove duplicates again (seems like most genes have multiple refseq transcripts...)
ref_id <- ref_id[!(duplicated(ref_id[,2]) | duplicated(ref_id[,2], fromLast=TRUE)),]

# add on the 3' lengths
the_lengths <- as.data.frame(the_lengths)
the_lengths2 <- the_lengths[the_lengths[,1] %in% ref_id[,1],]
# rownames(the_lengths2) <- the_lengths2[,1]
ref_id_full <- cbind.data.frame(ref_id,the_lengths2[match(ref_id[,1],the_lengths2[,1]),2])



### formalizing tests of batch factors
library(SoupX)

# specify main batch associated with the factor
b_name <- 'dmx_YS-JY-21_pool1'
b_name <- 'dmx_YS-JY-20_pool3'
b_direc <- 'down'
my_factor <- 1

# subset counts to just the one pool
cells_keep <- rownames(pbmc@meta.data)[pbmc@meta.data$batch_cov==b_name]
pbmc_sub <- subset(pbmc,cells = cells_keep)
pbmc_sub <- pbmc_sub@assays$RNA@counts

# get soup
toc = pbmc_sub
# Calculate soup profile
soupProf = data.frame(row.names = rownames(toc), est = rowSums(toc)/sum(toc), counts = rowSums(toc))

# get significant gene loadings
f_data <- get_one_factor(pbmc_container, factor_select=my_factor)
dscores <- f_data[[1]]
lds <- f_data[[2]]

sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=my_factor, colnames(lds))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

# limit to just the genes in tmp_casted_num
sig_df <- sig_df[rownames(lds),colnames(lds)]

# first test is to see if genes upregulated for batch are enriched for soup genes compared to downregulated genes

# second test is to see if among the upregulated genes, cross-cell type significant genes are enriched for soup genes
# compared to non-cross cell type significant genes

## alternatively could test cross-upregulated genes against all else
# get genes significant in all ctypes
all_sig <- rownames(sig_df)[rowSums(sig_df<0.05)==ncol(lds)]
# keep only genes with consistent signs in the loadings in direction of upregulation
tmp_lds <- lds[all_sig,]
if (b_direc=='up') {
    g_keep <- rownames(tmp_lds)[rowSums(tmp_lds>0)==ncol(lds)]
} else if (b_direc=='down') {
    g_keep <- rownames(tmp_lds)[rowSums(tmp_lds<0)==ncol(lds)]
}

g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]

g1 <- soupProf[g_keep,'est']
g2 <- soupProf[g_not_keep,'est']
all_p <- c(g1,g2)
all_t <- c(rep('g1',length(g1)),rep('g2',length(g2)))
mydf <- cbind.data.frame(all_p,all_t)
mydf <- as.data.frame(mydf)
ggplot(mydf,aes(x=as.factor(all_t),y=log(all_p))) +
    geom_boxplot()

wilcox.test(g1, g2, alternative = "two.sided")

# trying just significant consistent up vs significant consistent down
g_not_keep <- rownames(tmp_lds)[rowSums(tmp_lds>0)==ncol(lds)]

# trying just one cell type all sig up vs all sig down
ct_ndx <- 3
lds_sig <- lds[sig_df[,ct_ndx]<0.05,ct_ndx,drop=FALSE]
lds_up <- rownames(lds_sig)[lds_sig>0]
lds_down <- rownames(lds_sig)[lds_sig<0]
g1 <- soupProf[lds_up,'est']
g2 <- soupProf[lds_down,'est']
all_p <- c(g1,g2)
all_t <- c(rep('g1',length(g1)),rep('g2',length(g2)))
mydf <- cbind.data.frame(all_p,all_t)
mydf <- as.data.frame(mydf)
ggplot(mydf,aes(x=as.factor(all_t),y=log(all_p))) +
    geom_boxplot()
wilcox.test(g1, g2, alternative = "two.sided")

# trying just one cell type all sig vs not sig
g_keep <- rownames(lds_sig)
g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
g1 <- soupProf[g_keep,'est']
g2 <- soupProf[g_not_keep,'est']
all_p <- c(g1,g2)
all_t <- c(rep('g1',length(g1)),rep('g2',length(g2)))
mydf <- cbind.data.frame(all_p,all_t)
mydf <- as.data.frame(mydf)
ggplot(mydf,aes(x=as.factor(all_t),y=log(all_p))) +
    geom_boxplot()
wilcox.test(g1, g2, alternative = "two.sided")


# comparing sig up across all ctypes to sig up some ctypes
# any_sig <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
any_sig <- c()
for (i in 1:nrow(sig_df)) {
    for (j in 1:ncol(sig_df)) {
        sig_val <- sig_df[i,j]
        if (sig_val < 0.05) {
            if (lds[i,j] < 0) {
                any_sig <- c(any_sig,rownames(sig_df)[i])
                break
            }
        }
    }
}
g_keep # from before when did all consistent upreg
g_not_keep <- any_sig[!(any_sig %in% g_keep)]
g1 <- soupProf[g_keep,'est']
g2 <- soupProf[g_not_keep,'est']
all_p <- c(g1,g2)
all_t <- c(rep('g1',length(g1)),rep('g2',length(g2)))
mydf <- cbind.data.frame(all_p,all_t)
mydf <- as.data.frame(mydf)
ggplot(mydf,aes(x=as.factor(all_t),y=log(all_p))) +
    geom_boxplot()
wilcox.test(g1, g2, alternative = "two.sided")

# can probably improve this further by comparing ones significantly up in at 
# least 4 cell types vs those only in 1-3 ctypes
# Could also plot x=soup prob y=#ct sig up in



## seeing if there is association with gene length
my_factor <- 1
f_data <- get_one_factor(pbmc_container, factor_select=my_factor)
dscores <- f_data[[1]]
lds <- f_data[[2]]

sig_vectors <- get_significance_vectors(pbmc_container,
                                        factor_select=my_factor, colnames(lds))
# convert list to df
sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))

# limit to just the genes in tmp_casted_num
sig_df <- sig_df[rownames(lds),colnames(lds)]
ct_ndx <- 3
lds_sig <- lds[sig_df[,ct_ndx]<0.05,ct_ndx,drop=FALSE]

rownames(gen_dat) <- gen_dat$hgnc_symbol
# rownames(ref_id_full) <- ref_id_full[,2]


g_keep <- rownames(lds_sig)
g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
g1 <- gen_dat[g_keep,'gc'] # set to length or gc
g2 <- gen_dat[g_not_keep,'gc']
# g1 <- ref_id_full[g_keep,3] # set to length or gc
# g2 <- ref_id_full[g_not_keep,3]
g1 <- g1[!is.na(g1)]
g2 <- g2[!is.na(g2)]

all_p <- c(g1,g2)
all_t <- c(rep('significant_across',length(g1)),rep('not_significant',length(g2)))
mydf <- cbind.data.frame(all_p,all_t)
mydf <- as.data.frame(mydf)
ggplot(mydf,aes(x=as.factor(all_t),y=all_p)) +
    geom_boxplot() +
    ylab('GC%') +
    xlab('')
wilcox.test(g1, g2, alternative = "two.sided")


# trying just one cell type all sig up vs all sig down
ct_ndx <- 1
lds_sig <- lds[sig_df[,ct_ndx]<0.05,ct_ndx,drop=FALSE]
lds_up <- rownames(lds_sig)[lds_sig>0]
lds_down <- rownames(lds_sig)[lds_sig<0]
g1 <- gen_dat[lds_up,'gc'] # set to length or gc
g2 <- gen_dat[lds_down,'gc']
# g1 <- ref_id_full[g_keep,3] # set to length or gc
# g2 <- ref_id_full[g_not_keep,3]
g1 <- g1[!is.na(g1)]
g2 <- g2[!is.na(g2)]
all_p <- c(g1,g2)
all_t <- c(rep('significant_across',length(g1)),rep('not_significant',length(g2)))
mydf <- cbind.data.frame(all_p,all_t)
mydf <- as.data.frame(mydf)
ggplot(mydf,aes(x=as.factor(all_t),y=all_p)) +
    geom_boxplot() +
    ylab('GC%') +
    xlab('')
wilcox.test(g1, g2, alternative = "two.sided")
length(g1)
length(g2)



test_gc_association(pbmc_container,my_factor=2,b_direc='down',comp_type='sig')
test_gc_association(pbmc_container,my_factor=16,b_direc='down',comp_type='any')

test_gc_association <- function(container,my_factor,b_direc,comp_type='sig') {
    f_data <- get_one_factor(container, factor_select=my_factor)
    dscores <- f_data[[1]]
    lds <- f_data[[2]]
    
    if (is.null(container$gen_dat)) {
        library(biomaRt)
        library(EDASeq)
        
        ensembl = useMart("ensembl",
                          dataset="hsapiens_gene_ensembl")
        ens_id <- getBM(attributes=c('ensembl_gene_id',
                                     'hgnc_symbol'),
                        filters = 'hgnc_symbol',
                        mart = ensembl,
                        values = rownames(lds))
        
        gc_len <- getGeneLengthAndGCContent(ens_id[,1], 'hsa')
        
        gen_dat <- cbind.data.frame(ens_id,gc_len)
        
        gen_dat <- gen_dat[!(duplicated(gen_dat[,2]) | duplicated(gen_dat[,2], fromLast=TRUE)),]
        
        container$gen_dat <- gen_dat
    }
    
    gen_dat <- container$gen_dat
    rownames(gen_dat) <- gen_dat$hgnc_symbol
    
    sig_vectors <- get_significance_vectors(container,
                                            factor_select=my_factor, colnames(lds))
    # convert list to df
    sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
    
    # limit to just the genes in tmp_casted_num
    sig_df <- sig_df[rownames(lds),colnames(lds)]
    
    # loop through cell types to calculate 
    for (j in 1:ncol(lds)) {
        
        ct_ndx <- j
        ct_name <- colnames(lds)[j]
        
        if (comp_type=='sig') {
            lds_sig <- lds[sig_df[,ct_ndx]<0.05,ct_ndx,drop=FALSE]
            g_keep <- rownames(lds_sig)
            g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
            g1 <- gen_dat[g_keep,'gc'] # set to length or gc
            g2 <- gen_dat[g_not_keep,'gc']
            g1 <- g1[!is.na(g1)]
            g2 <- g2[!is.na(g2)]
            tres <- t.test(g1, g2, alternative = "two.sided")
            
            print(ct_name)
            print(tres$p.value)
        }
    }
    
    if (comp_type=="any_up") {
        # g_keep <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
        # # only keep ones that are upregulated in any cell type
        # lds_sub <- lds[g_keep,]
        # if (b_direc=='up') {
        #     g_keep <- rownames(lds_sub)[rowSums(lds_sub>0)>0]
        # } else if (b_direc=='down') {
        #     g_keep <- rownames(lds_sub)[rowSums(lds_sub<0)>0]
        # }
        
        g_keep <- c()
        for (i in 1:nrow(sig_df)) {
            for (j in 1:ncol(sig_df)) {
                sig_val <- sig_df[i,j]
                if (sig_val < 0.05) {
                    if (b_direc=='up') {
                        if (lds[i,j] > 0) {
                            g_keep <- c(g_keep,rownames(sig_df)[i])
                            break
                        }
                    } else if (b_direc=='down') {
                        if (lds[i,j] < 0) {
                            g_keep <- c(g_keep,rownames(sig_df)[i])
                            break
                        }
                    }
                }
            }
        }
        
        g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
        g1 <- gen_dat[g_keep,'gc'] # set to length or gc
        g2 <- gen_dat[g_not_keep,'gc']
        g1 <- g1[!is.na(g1)]
        g2 <- g2[!is.na(g2)]
        tres <- t.test(g1, g2, alternative = "two.sided")
        
        print(tres$p.value)

        all_p <- c(g1,g2)
        all_t <- c(rep('sig_genes',length(g1)),rep('NS_genes',length(g2)))
        all_f <- rep(paste0('factor_',my_factor),length(all_p))
        mydf <- cbind.data.frame(all_p,all_t,all_f)
        mydf <- as.data.frame(mydf)

        return(list(tres$p.value,mydf))
    } else if (comp_type=='any') {
        g_keep <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
        g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
        g1 <- gen_dat[g_keep,'gc'] # set to length or gc
        g2 <- gen_dat[g_not_keep,'gc']
        g1 <- g1[!is.na(g1)]
        g2 <- g2[!is.na(g2)]
        tres <- t.test(g1, g2, alternative = "two.sided")
        
        print(tres$p.value)
        
        all_p <- c(g1,g2)
        all_t <- c(rep('sig_genes',length(g1)),rep('NS_genes',length(g2)))
        all_f <- rep(paste0('factor_',my_factor),length(all_p))
        mydf <- cbind.data.frame(all_p,all_t,all_f)
        mydf <- as.data.frame(mydf)
        
        return(list(tres$p.value,mydf))
    }
}

all_pv <- c()
pv_list <- list()
my_direcs <- c('down','up','up','up','up','down','up','down','up','up','up','up','up','up','down','up','up','up','up','up')
my_direcs <- rep('up',ncol(pbmc_container[["tucker_results"]][[1]]))
my_direcs <- rep('down',ncol(pbmc_container[["tucker_results"]][[1]]))
for (i in 1:ncol(pbmc_container[["tucker_results"]][[1]])) {
    print(i)
    pv <- test_gc_association(pbmc_container,my_factor=i,b_direc=my_direcs[i],comp_type='any_up')
    # pv <- test_gc_association(pbmc_container,my_factor=i,b_direc='up',comp_type='any')
    all_pv <- c(all_pv,pv[[1]])
    pv_list[[i]] <- pv[[2]]
}
all_pv_adj <- p.adjust(all_pv,method='fdr')

pv_no_combat <- all_pv
pv_combat <- all_pv
all_pv_adj <- p.adjust(c(pv_no_combat,pv_combat),method='fdr')

plt_no_combat <- pv_list
plt_combat <- pv_list


# make boxplots for select factors where significant genes associated with gene sets. Include some control factors
# full_df <- rbind.data.frame(pv_list[[1]],pv_list[[2]],pv_list[[4]],pv_list[[8]],pv_list[[10]],pv_list[[15]],pv_list[[17]],pv_list[[3]],
#                             pv_list[[7]],pv_list[[13]],pv_list[[14]],pv_list[[16]],pv_list[[18]],pv_list[[19]],pv_list[[20]])
full_df <- rbind.data.frame(plt_combat[[1]],plt_combat[[4]],plt_combat[[6]],plt_combat[[8]],plt_combat[[10]],plt_combat[[17]])
full_df <- rbind.data.frame(plt_no_combat[[1]],plt_no_combat[[2]],plt_no_combat[[3]],
                            plt_no_combat[[4]],plt_no_combat[[5]],plt_no_combat[[6]],
                            plt_no_combat[[7]],plt_no_combat[[10]])

full_df$all_f <- as.factor(full_df$all_f)
full_df$all_t <- as.factor(full_df$all_t)
# full_df$all_f <- factor(full_df$all_f,levels=c('factor_1','factor_2','factor_4','factor_8','factor_10','factor_15','factor_17',
#                                                'factor_3','factor_7','factor_13','factor_14','factor_16','factor_18','factor_19','factor_20'))
full_df$all_f <- factor(full_df$all_f,levels=c('factor_1','factor_4','factor_6','factor_8','factor_10','factor_17'))
full_df$all_f <- factor(full_df$all_f,levels=c('factor_1','factor_2','factor_3','factor_4','factor_5','factor_6','factor_7','factor_10'))
p <- ggplot(full_df,aes(x = all_f, y = all_p, fill = all_t)) +
    geom_boxplot() +
    xlab('') +
    ylab('GC %') +
    ylim(.3,.8) +
    ggtitle('With ComBat batch correction') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5,size=18),
          axis.text=element_text(size=14),
          axis.title=element_text(size=16))
pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_gc_no_combat.pdf", useDingbats = FALSE,
    width = 9, height = 4.5)
pdf(file = "/home/jmitchel/figures/for_paper/lupus_batch_gc_with_combat.pdf", useDingbats = FALSE,
    width = 10.75, height = 4.5)
p
dev.off()



# # trying to add pval text on top of boxes
# pvs1 <- all_pv_adj[1:10]
# all_f_ndx <- sapply(full_df$all_f,function(x){
#     return(as.numeric(strsplit(as.character(x),split='_')[[1]][[2]]))
# })
# full_df$pvs <- sapply(all_f_ndx,function(x) {
#     return(pvs1[x])
# })
# 
# 
# p + geom_text(data=full_df,aes(x=all_f,y=all_p,label=pvs))
# p
# 
# ## trying different way of adding labels
# f_ndx <- c(1,2,3,4,5,6,7,10)
# f_names <- c('factor_1','factor_2','factor_3','factor_4','factor_5','factor_6','factor_7','factor_10')
# pvs1 <- all_pv_adj[1:10]
# pvs1 <- pvs1[f_ndx]
# max_p <- sapply(f_names,function(x) {
#     max(full_df$all_p[full_df$all_f==x])
# })
# mydf <- cbind.data.frame(f_names,
#                          max_p,
#                          pvs1)
# colnames(mydf) <- c('factor_name','max_y','pvalues')
# mydf$factor_name <- factor(mydf$factor_name,levels=f_names)
# p + geom_text(data=data.frame(), aes(x = mydf$factor_name, y = mydf$max_y, label = mydf$pvalues))




## test whether batch genes have longer poly-a tails
polya <- read.csv(file='/home/jmitchel/poly_a_hela.csv')
rownames(polya) <- polya$gene
polya <- polya[,2,drop=FALSE]
pbmc_container$gen_dat <- polya

pv <- test_polya_association(pbmc_container,my_factor=9,b_direc='down',comp_type='any_up')

all_pv <- c()
pv_list <- list()
my_direcs <- c('down','up','up','up','up','down','up','down','up','up','up','up','up','up','down','up','up','up','up','up')
for (i in 1:ncol(pbmc_container[["tucker_results"]][[1]])) {
    print(i)
    pv <- test_polya_association(pbmc_container,my_factor=i,b_direc=my_direcs[i],comp_type='any')
    # pv <- test_gc_association(pbmc_container,my_factor=i,b_direc='up',comp_type='any')
    all_pv <- c(all_pv,pv[[1]])
    pv_list[[i]] <- pv[[2]]
}
all_pv_adj <- p.adjust(all_pv,method='fdr')

test_polya_association <- function(container,my_factor,b_direc,comp_type='sig') {
    f_data <- get_one_factor(container, factor_select=my_factor)
    dscores <- f_data[[1]]
    lds <- f_data[[2]]
    
    gen_dat <- container$gen_dat

    sig_vectors <- get_significance_vectors(container,
                                            factor_select=my_factor, colnames(lds))
    # convert list to df
    sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
    
    # limit to just the genes in tmp_casted_num
    sig_df <- sig_df[rownames(lds),colnames(lds)]
    
    
    if (comp_type=="any_up") {
        # g_keep <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
        # # only keep ones that are upregulated in any cell type
        # lds_sub <- lds[g_keep,]
        # if (b_direc=='up') {
        #     g_keep <- rownames(lds_sub)[rowSums(lds_sub>0)>0]
        # } else if (b_direc=='down') {
        #     g_keep <- rownames(lds_sub)[rowSums(lds_sub<0)>0]
        # }
        
        g_keep <- c()
        for (i in 1:nrow(sig_df)) {
            for (j in 1:ncol(sig_df)) {
                sig_val <- sig_df[i,j]
                if (sig_val < 0.05) {
                    if (b_direc=='up') {
                        if (lds[i,j] > 0) {
                            g_keep <- c(g_keep,rownames(sig_df)[i])
                            break
                        }
                    } else if (b_direc=='down') {
                        if (lds[i,j] < 0) {
                            g_keep <- c(g_keep,rownames(sig_df)[i])
                            break
                        }
                    }
                }
            }
        }
        
        g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
        g1 <- gen_dat[g_keep,'med_len'] # set to length or gc
        g2 <- gen_dat[g_not_keep,'med_len']
        g1 <- g1[!is.na(g1)]
        g2 <- g2[!is.na(g2)]
        tres <- t.test(g1, g2, alternative = "two.sided")
        
        print(tres$p.value)
        
        all_p <- c(g1,g2)
        all_t <- c(rep('sig_genes',length(g1)),rep('NS_genes',length(g2)))
        all_f <- rep(paste0('factor_',my_factor),length(all_p))
        mydf <- cbind.data.frame(all_p,all_t,all_f)
        mydf <- as.data.frame(mydf)
        
        return(list(tres$p.value,mydf))
    } else if (comp_type=='any') {
        g_keep <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
        g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
        g1 <- gen_dat[g_keep,'med_len'] # set to length or gc
        g2 <- gen_dat[g_not_keep,'med_len']
        g1 <- g1[!is.na(g1)]
        g2 <- g2[!is.na(g2)]
        tres <- t.test(g1, g2, alternative = "two.sided")
        
        print(tres$p.value)
        
        all_p <- c(g1,g2)
        all_t <- c(rep('sig_genes',length(g1)),rep('NS_genes',length(g2)))
        all_f <- rep(paste0('factor_',my_factor),length(all_p))
        mydf <- cbind.data.frame(all_p,all_t,all_f)
        mydf <- as.data.frame(mydf)
        
        return(list(tres$p.value,mydf))
    }
}


# ## testing soupX example
# library(SoupX)
# 
# tmpDir = tempdir(check = TRUE)
# download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz", 
#               destfile = file.path(tmpDir, "tod.tar.gz"))
# download.file("https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_filtered_gene_bc_matrices.tar.gz", 
#               destfile = file.path(tmpDir, "toc.tar.gz"))
# untar(file.path(tmpDir, "tod.tar.gz"), exdir = tmpDir)
# untar(file.path(tmpDir, "toc.tar.gz"), exdir = tmpDir)
# sc = load10X(tmpDir)
# 
# ## seeing if i can get some multi batch data other ways...
# # code from cellmixs work on kang dataset
# library(ExperimentHub)
# sc = ExperimentHub()
# 
# sce <- sc[["EH2259"]]
# ## Filter out genes that are not expressed in any cell
# sce <- sce[which(rowSums(counts(sce) > 0) > 0), ]
# sce$patient <- sce$ind
# 
# sce$patient <- factor(sce$patient)
# 
# table(sce$patient)
# 
# sce <- sce[,!sce$stim %in% "stim"]
# dim(sce)
# sce <- sce[,!sce$multiplets %in% "doublet"]
# 
# suppressPackageStartupMessages({
#     library(scran)
#     library(magrittr)
#     library(dplyr)
#     library(CellBench)
#     library(plyr)
#     library(EnsDb.Hsapiens.v86)
#     library(AnnotationDbi)
# })
# 
# seed <- 1000
# sc_data <- load_sc_data()
# 
# colData(sc_data[[1]])$protocol <- rep(names(sc_data)[1], ncol(sc_data[[1]]))
# sce <- sc_data[[1]]
# 
# for(i in 2:length(sc_data)){
#     colData(sc_data[[i]])$protocol <- rep(names(sc_data)[i], ncol(sc_data[[i]]))
#     gene_overlap <- intersect(rownames(sce), rownames(sc_data[[i]]))
#     coldata_overlap <- intersect(names(colData(sce)), names(colData(sc_data[[i]])))
#     sc_data[[i]] <- sc_data[[i]][gene_overlap,]
#     colData(sc_data[[i]]) <- colData(sc_data[[i]])[, coldata_overlap]
#     colData(sce) <- colData(sce)[, coldata_overlap]
#     sce <- sce[gene_overlap,]
#     sce <- cbind(sce, sc_data[[i]])
# }
# colnames(sce) <- paste0(colnames(sce), "_", sce$protocol)
# dim(sce)









## run gsea on factor 9 to see if it's enriched for cell viability markers and/or mitochondrial genes
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=9, method="fgsea", thresh=0.05,
                                      db_use=c("GO"))
plot_gsea_hmap_w_similarity(pbmc_container,factor_select=9,direc='down',thresh=.05)
plot_gsea_sub(pbmc_container,thresh=.05,clust_select=8)
# look at genes to see if they're mitochondrial
f1_data <- get_one_factor(pbmc_container, factor_select=9)
dscores <- f1_data[[1]]
lds <- f1_data[[2]]
lds_ct <- lds[,4]
lds_ct[order(lds_ct,decreasing=F)][1:50]





## now trying soupx with one test batch
pbmc_batch <- pbmc_container
# cool <- test_soup_association(pbmc_batch,6,'down','all_v_some',soupProf)
res1 <- test_soup_association(pbmc_batch,6,'down','any_up',soupProf)
# cool <- test_soup_association(pbmc_batch,6,'down','any',soupProf)
# cool <- test_soup_association(pbmc_batch,6,'down','sig_up_consist',soupProf)
# cool <- test_soup_association(pbmc_batch,6,'down','consist_up_v_down',soupProf)
# cool <- test_soup_association(pbmc_batch,6,'down','any_up_v_down',soupProf)


# cool <- test_soup_association(pbmc_container,10,'up','all_v_some',soupProf)

test_soup_association <- function(container,my_factor,b_direc,comp_type='any_up',soupProf) {
    f_data <- get_one_factor(container, factor_select=my_factor)
    dscores <- f_data[[1]]
    lds <- f_data[[2]]
    
    sig_vectors <- get_significance_vectors(container,
                                            factor_select=my_factor, colnames(lds))
    # convert list to df
    sig_df <- t(as.data.frame(do.call(rbind, sig_vectors)))
    
    # limit to just the genes in tmp_casted_num
    sig_df <- sig_df[rownames(lds),colnames(lds)]
    
    
    if (comp_type=="any_up") {
        g_keep <- c()
        for (i in 1:nrow(sig_df)) {
            for (j in 1:ncol(sig_df)) {
                sig_val <- sig_df[i,j]
                if (sig_val < 0.05) {
                    if (b_direc=='up') {
                        if (lds[i,j] > 0) {
                            g_keep <- c(g_keep,rownames(sig_df)[i])
                            break
                        }
                    } else if (b_direc=='down') {
                        if (lds[i,j] < 0) {
                            g_keep <- c(g_keep,rownames(sig_df)[i])
                            break
                        }
                    }
                }
            }
        }
        
        g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
        g1 <- soupProf[g_keep,'est'] # set to length or gc
        g2 <- soupProf[g_not_keep,'est']
        g1 <- g1[!is.na(g1)]
        g1_nonzero <- g1[g1!=0]
        g1_nonzero_min <- min(g1_nonzero)
        g1[g1==0] <- g1_nonzero_min
        g2 <- g2[!is.na(g2)]
        g2_nonzero <- g2[g2!=0]
        g2_nonzero_min <- min(g2_nonzero)
        g2[g2==0] <- g2_nonzero_min
        g1 <- log10(g1)
        g2 <- log10(g2)
        tres <- t.test(g1, g2, alternative = "two.sided")
        
        print(tres$p.value)
        
        all_p <- c(g1,g2)
        all_t <- c(rep('Upreg genes',length(g1)),rep('NS genes',length(g2)))
        all_f <- rep(paste0('factor_',my_factor),length(all_p))
        mydf <- cbind.data.frame(all_p,all_t,all_f)
        mydf <- as.data.frame(mydf)
        
        return(list(tres$p.value,mydf))
    } else if (comp_type=='any') {
        g_keep <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
        g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
        g1 <- soupProf[g_keep,'est'] # set to length or gc
        g2 <- soupProf[g_not_keep,'est']
        g1 <- g1[!is.na(g1)]
        g1_nonzero <- g1[g1!=0]
        g1_nonzero_min <- min(g1_nonzero)
        g1[g1==0] <- g1_nonzero_min
        g2 <- g2[!is.na(g2)]
        g2_nonzero <- g2[g2!=0]
        g2_nonzero_min <- min(g2_nonzero)
        g2[g2==0] <- g2_nonzero_min
        g1 <- log10(g1)
        g2 <- log10(g2)
        tres <- t.test(g1, g2, alternative = "two.sided")
        
        print(tres$p.value)
        
        all_p <- c(g1,g2)
        all_t <- c(rep('sig_genes',length(g1)),rep('NS_genes',length(g2)))
        all_f <- rep(paste0('factor_',my_factor),length(all_p))
        mydf <- cbind.data.frame(all_p,all_t,all_f)
        mydf <- as.data.frame(mydf)
        
        return(list(tres$p.value,mydf))
    } else if (comp_type=='sig_up_consist') {
        all_sig <- rownames(sig_df)[rowSums(sig_df<0.05)==ncol(lds)]
        # keep only genes with consistent signs in the loadings in direction of upregulation
        tmp_lds <- lds[all_sig,]
        if (b_direc=='up') {
            g_keep <- rownames(tmp_lds)[rowSums(tmp_lds>0)==ncol(lds)]
        } else if (b_direc=='down') {
            g_keep <- rownames(tmp_lds)[rowSums(tmp_lds<0)==ncol(lds)]
        }
        
        g_not_keep <- rownames(lds)[!(rownames(lds) %in% g_keep)]
        
        g1 <- soupProf[g_keep,'est']
        g2 <- soupProf[g_not_keep,'est']
        tres <- t.test(g1, g2, alternative = "two.sided")
        
        print(tres$p.value)
        
        all_p <- c(g1,g2)
        all_t <- c(rep('g1',length(g1)),rep('g2',length(g2)))
        mydf <- cbind.data.frame(all_p,all_t)
        mydf <- as.data.frame(mydf)
        return(list(tres$p.value,mydf))
    } else if (comp_type=='all_v_some') {
        any_up <- c()
        for (i in 1:nrow(sig_df)) {
            for (j in 1:ncol(sig_df)) {
                sig_val <- sig_df[i,j]
                if (sig_val < 0.05) {
                    if (b_direc=='up') {
                        if (lds[i,j] > 0) {
                            any_up <- c(any_up,rownames(sig_df)[i])
                            break
                        }
                    } else if (b_direc=='down') {
                        if (lds[i,j] < 0) {
                            any_up <- c(any_up,rownames(sig_df)[i])
                            break
                        }
                    }
                }
            }
        }
        
        
        all_sig <- rownames(sig_df)[rowSums(sig_df<0.05)==ncol(lds)]
        # keep only genes with consistent signs in the loadings in direction of upregulation
        tmp_lds <- lds[all_sig,]
        if (b_direc=='up') {
            g_keep <- rownames(tmp_lds)[rowSums(tmp_lds>0)==ncol(lds)]
        } else if (b_direc=='down') {
            g_keep <- rownames(tmp_lds)[rowSums(tmp_lds<0)==ncol(lds)]
        }
        
        g_not_keep <- any_up[!(any_up %in% g_keep)]
        g1 <- soupProf[g_keep,'est'] # set to length or gc
        g2 <- soupProf[g_not_keep,'est']
        g1 <- g1[!is.na(g1)]
        g1_nonzero <- g1[g1!=0]
        g1_nonzero_min <- min(g1_nonzero)
        g1[g1==0] <- g1_nonzero_min
        g2 <- g2[!is.na(g2)]
        g2_nonzero <- g2[g2!=0]
        g2_nonzero_min <- min(g2_nonzero)
        g2[g2==0] <- g2_nonzero_min
        g1 <- log10(g1)
        g2 <- log10(g2)
        tres <- t.test(g1, g2, alternative = "two.sided")
        
        print(tres$p.value)
        
        all_p <- c(g1,g2)
        all_t <- c(rep('Upreg all\ncell types',length(g1)),rep('Upreg some\ncell types',length(g2)))
        all_f <- rep(paste0('factor_',my_factor),length(all_p))
        mydf <- cbind.data.frame(all_p,all_t,all_f)
        mydf <- as.data.frame(mydf)
        
        return(list(tres$p.value,mydf))
    } else if (comp_type=='consist_up_v_down') {
        all_sig <- rownames(sig_df)[rowSums(sig_df<0.05)==ncol(lds)]
        # keep only genes with consistent signs in the loadings in direction of upregulation
        tmp_lds <- lds[all_sig,]
        if (b_direc=='up') {
            g_keep <- rownames(tmp_lds)[rowSums(tmp_lds>0)==ncol(lds)]
            g_not_keep <- rownames(tmp_lds)[rowSums(tmp_lds<0)==ncol(lds)]
        } else if (b_direc=='down') {
            g_keep <- rownames(tmp_lds)[rowSums(tmp_lds<0)==ncol(lds)]
            g_not_keep <- rownames(tmp_lds)[rowSums(tmp_lds>0)==ncol(lds)]
        }
        
        g1 <- soupProf[g_keep,'est']
        g2 <- soupProf[g_not_keep,'est']
        g1 <- g1[!is.na(g1)]
        g1_nonzero <- g1[g1!=0]
        g1_nonzero_min <- min(g1_nonzero)
        g1[g1==0] <- g1_nonzero_min
        g2 <- g2[!is.na(g2)]
        g2_nonzero <- g2[g2!=0]
        g2_nonzero_min <- min(g2_nonzero)
        g2[g2==0] <- g2_nonzero_min
        g1 <- log10(g1)
        g2 <- log10(g2)
        tres <- t.test(g1, g2, alternative = "two.sided")
        
        print(tres$p.value)
        
        all_p <- c(g1,g2)
        all_t <- c(rep('g1',length(g1)),rep('g2',length(g2)))
        mydf <- cbind.data.frame(all_p,all_t)
        mydf <- as.data.frame(mydf)
        return(list(tres$p.value,mydf))
    }  else if (comp_type=='any_up_v_down') {
        any_up <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
        lds_sub <- lds[any_up,]
        if (b_direc=='up') {
            any_up <- rownames(lds_sub)[rowSums(lds_sub>0)>0]
        } else if (b_direc=='down') {
            any_up <- rownames(lds_sub)[rowSums(lds_sub<0)>0]
        }
        
        any_down <- rownames(sig_df)[rowSums(sig_df<0.05)>0]
        lds_sub <- lds[any_down,]
        if (b_direc=='up') {
            any_down <- rownames(lds_sub)[rowSums(lds_sub<0)>0]
        } else if (b_direc=='down') {
            any_down <- rownames(lds_sub)[rowSums(lds_sub>0)>0]
        }
        
        # only keep genes not in intersection of any_up and any_down
        int_genes <- intersect(any_up,any_down)
        any_up <- any_up[!(any_up %in% int_genes)]
        any_down <- any_down[!(any_down %in% int_genes)]
        
        g1 <- soupProf[any_up,'est']
        g2 <- soupProf[any_down,'est']
        tres <- t.test(log(g1), log(g2), alternative = "two.sided")
        
        print(tres$p.value)
        
        all_p <- c(g1,g2)
        all_t <- c(rep('g1',length(g1)),rep('g2',length(g2)))
        mydf <- cbind.data.frame(all_p,all_t)
        mydf <- as.data.frame(mydf)
        return(list(tres$p.value,mydf))
    }
}


# all_pv <- c()
# pv_list <- list()
# my_direcs <- rep('down',10)
# for (i in 3:10) {
#     print(i)
#     # pv <- test_soup_association(pbmc_container,my_factor=i,b_direc=my_direcs[i],comp_type='all_v_some',soupProf)
#     pv <- test_soup_association(pbmc_container,my_factor=i,b_direc=my_direcs[i],comp_type='any_up_v_down',soupProf)
#     all_pv <- c(all_pv,pv[[1]])
#     pv_list[[i]] <- pv[[2]]
# }
# all_pv_adj <- p.adjust(all_pv,method='fdr')

# res1 <- test_soup_association(pbmc_batch,6,'down','any',soupProf)
res1 <- test_soup_association(pbmc_batch,6,'down','any_up',soupProf)
# res1 <- test_soup_association(pbmc_batch,6,'down','consist_up_v_down',soupProf)

p <- ggplot(res1[[2]],aes(x = as.factor(all_t), y = all_p)) +
    geom_boxplot() +
    xlab('') +
    ylab('log10(soup fraction)') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5,size=18),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
pdf(file = "/home/jmitchel/figures/for_paper/soup_any_up.pdf", useDingbats = FALSE,
    width = 3.5, height = 4)
p
dev.off()

res2 <- test_soup_association(pbmc_batch,6,'down','all_v_some',soupProf)
p <- ggplot(res2[[2]],aes(x = as.factor(all_t), y = all_p)) +
    geom_boxplot() +
    xlab('') +
    ylab('log10(soup fraction)') +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5,size=18),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14))
pdf(file = "/home/jmitchel/figures/for_paper/soup_all_v_some.pdf", useDingbats = FALSE,
    width = 3.5, height = 4)
p
dev.off()










