
library(scITD)
library(Seurat)
library(ggplot2)

# load up the dataset
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')

# using full dataset to generate subclusters, also using old ct names since this was
# used to get subcluster embeddings
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


pbmc_container <- form_tensor(pbmc_container, donor_min_cells=20,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.15,
                              scale_var = TRUE, var_scale_power = .5,
                              batch_var='pool')


pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(7,20),
                                 tucker_type = 'regular', rotation_type = 'hybrid')

pbmc_container <- get_meta_associations(pbmc_container,vars_test=c('sex','Age','pool','processing','Status','Ethnicity'),
                                        stat_use='pval')

## plot donor score
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('sex'),
                                    cluster_by_meta = 'sex',
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

# just to check that everything is as expected up to this point
pbmc_container$plots$donor_matrix














##### getting subclusters
# add conos object generated in get_embedding file
con <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_conos2.rds')
pbmc_container$embedding <- con

# in case below fn errors in the middle, need to save these objects
orig_embed <- pbmc_container$embedding[["embedding"]]
orig_clusts <- pbmc_container$embedding$clusters$leiden$groups

# # to recover original embedding/cell assignments
# pbmc_container$embedding[["embedding"]] <- orig_embed
# pbmc_container$embedding$clusters$leiden$groups <- orig_clusts

# large number of cores seems to hamper some stuff below
pbmc_container$embedding$n.cores <- 5

pbmc_container <- get_subtype_prop_associations(pbmc_container, max_res=.9, stat_type='adj_pval',
                                                min_cells_group=200)


# saveRDS(pbmc_container$subclusters,file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')
pbmc_container$subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')


### generate figure with all ctype information for all ctypes/factors
# first determine what resolution of each ctype to choose
all_ctypes=c('cM','cM','cM','cM',
             'ncM','ncM','ncM','ncM',
             'cDC','cDC','cDC','cDC',
             'B','B','B','B',
             'T4','T4','T4','T4',
             'T8','T8','T8','T8',
             'NK','NK','NK','NK')
all_res=c(.5,.7,.8,.9,
          .5,.6,.7,.8,
          .5,.7,.8,.9,
          .6,.7,.8,.9,
          .5,.6,.8,.9,
          .5,.6,.8,.9,
          .6,.7,.8,.9)


# get all subcluster umaps
pbmc_container <- get_subclust_umap(pbmc_container,all_ctypes=all_ctypes,
                                    all_res=all_res,n_col=4)

## subc plots at resolutions used for downstream analysis
pbmc_container[["plots"]][["subc_umaps"]][["B:0.8"]] +  scale_y_reverse() +
  scale_x_reverse()
pbmc_container[["plots"]][["subc_umaps"]][["T4:0.6"]] +  scale_y_reverse() +
  scale_x_reverse()
pbmc_container[["plots"]][["subc_umaps"]][["T8:0.6"]] +  scale_y_reverse() +
  scale_x_reverse()
pbmc_container[["plots"]][["subc_umaps"]][["cM:0.5"]] +  scale_y_reverse() +
  scale_x_reverse()
















##### now do the proportion analysis using the SLE-only portion of the dataset
## also using the new cell type naming scheme
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

pbmc_container <- plot_donor_matrix(pbmc_container,
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

pbmc_container$plots$donor_matrix


# add conos object generated in embedding prep file
con <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_conos2.rds')
pbmc_container$embedding <- con

# add above subclustering to the container
pbmc_container$subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')

## since I forgot to have the fn output the pvalues I'm going to rerun a modified version below
# it's already updated in dev version of the package
get_subtype_prop_associations_tmp <- function(container, max_res, stat_type,
                                          integration_var=NULL, min_cells_group=50,
                                          use_existing_subc=FALSE,
                                          alt_ct_names=NULL,n_col=2) {
  if (!(stat_type %in% c("fstat","adj_rsq","adj_pval"))) {
    stop("stat_type parameter is not one of the three options")
  }
  
  if (is.null(integration_var)) {
    if (!use_existing_subc) {
      if (is.null(container$embedding)) {
        stop("need to set integration_var parameter to get an embedding")
      }
    }
  } else {
    container <- reduce_dimensions(container,integration_var)
  }
  
  # make sure that groups doesn't contain cell types not present
  container$embedding$clusters$leiden$groups <- factor(container$embedding$clusters$leiden$groups,
                                                       levels=unique(container$embedding$clusters$leiden$groups))
  
  donor_scores <- container$tucker_results[[1]]
  
  # create dataframe to store association results
  res <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(res) <- c(stat_type,'resolution','factor','ctype')
  
  # make list to store subclustering results
  if (use_existing_subc) {
    subc_all <- container$subclusters
  } else {
    subc_all <- list()
  }
  # loop through cell types
  for (ct in container$experiment_params$ctypes_use) {
    print(ct)
    scMinimal <- container[["scMinimal_ctype"]][[ct]]
    
    # loop through increasing clustering resolutions
    cluster_res <- seq(.5,max_res,by=.1)
    for (r in cluster_res) {
      if (!use_existing_subc) {
        print(r)
        # run clustering
        # subclusts <- get_subclusters(container,ct,r,min_cells_group=min_cells_group,
        #                              small_clust_action='merge')
        subclusts <- get_subclusters(container,ct,r,min_cells_group=min_cells_group,
                                     small_clust_action='remove')
        subclusts <- subclusts + 1 # moves subcluster index from 0 to 1
        subc_all[[ct]][[paste0('res:',as.character(r))]] <- subclusts
      } else {
        if (!is.null(alt_ct_names)) {
          ct_ndx <- which(container$experiment_params$ctypes_use==ct)
          ct_new <- alt_ct_names[ct_ndx]
          subclusts <- container$subclusters[[ct_new]][[paste0('res:',as.character(r))]]
        } else {
          subclusts <- container$subclusters[[ct]][[paste0('res:',as.character(r))]]
        }
      }
      
      num_subclusts <- length(unique(subclusts))
      
      if (num_subclusts > 1) {
        # get cells in both metadata and subclusts
        cell_intersect <- intersect(names(subclusts),rownames(scMinimal$metadata))
        
        sub_meta_tmp <- scMinimal$metadata[cell_intersect,]
        
        # get donor proportions of subclusters
        donor_props <- compute_donor_props(subclusts,sub_meta_tmp)
        
        # transform from proportions to balances
        donor_balances <- coda.base::coordinates(donor_props)
        rownames(donor_balances) <- rownames(donor_props)
        
        # compute regression statistics
        reg_stats <- compute_associations(donor_balances,donor_scores,stat_type)
        
        # rename donor_props columns for generating plot of donor proportions and scores
        colnames(donor_props) <- sapply(1:ncol(donor_props),function(x){paste0(ct,'_',x)})
        
      } else {
        if (stat_type=='fstat' || stat_type=='adj_rsq') {
          reg_stats <- rep(0,ncol(container$tucker_results[[1]]))
        } else if (stat_type=='adj_pval') {
          reg_stats <- rep(1,ncol(container$tucker_results[[1]]))
        }
      }
      
      # store association results
      for (i in 1:length(reg_stats)) {
        new_row <- as.data.frame(list(reg_stats[i], r, paste0("Factor ", as.character(i)), ct),stringsAsFactors = F)
        colnames(new_row) <- colnames(res)
        res <- rbind(res,new_row)
      }
    }
  }
  
  # adjust p-values if using adj_pval stat_type
  if (stat_type=='adj_pval') {
    res$adj_pval <- p.adjust(res$adj_pval,method = 'fdr')
  }
  
  # generate plot of associations
  reg_stat_plots <- plot_subclust_associations(res,n_col=n_col)
  
  # save results
  container$plots$subtype_prop_factor_associations <- reg_stat_plots
  container$subclusters <- subc_all
  container$subc_factor_association_res <- res
  
  return(container)
}


# generate all subcluster-factor associations using the above subclusterings
pbmc_container <- get_subtype_prop_associations_tmp(pbmc_container, max_res=.9,
                                                use_existing_subc=TRUE,
                                                stat_type='adj_pval',
                                                min_cells_group=200,
                                                alt_ct_names=c('B','NK','T4','T8','cDC','cM','ncM'))


pv_table <- pbmc_container$subc_factor_association_res



# also need this tmp fig generation fn as it was updated in develop branch
plot_subclust_associations_tmp <- function(res,n_col=2) {
  
  stat_type <- colnames(res)[1]
  
  if (stat_type == 'adj_pval') {
    res[,stat_type] <- -log10(res[,stat_type])
  }
  
  if (stat_type=='fstat') {
    y_axis_name <- 'F-Statistic'
  } else if (stat_type=='adj_rsq') {
    y_axis_name <- 'adj r-sq'
  } else if (stat_type == 'adj_pval') {
    y_axis_name <- '-log10(adj p-val)'
  }
  
  num_factors <- length(unique(res$factor))
  ctypes <- unique(res$ctype)
  plot_list <- list()
  
  for (f in 1:num_factors) {
    factor_name <- paste0("Factor ",as.character(f))
    res_factor <- res[res$factor==factor_name,]
    
    p <- ggplot(res_factor,aes_string(x='resolution',y=stat_type,color='ctype')) +
      geom_line() +
      xlab("Leiden Resolution") +
      ylab(y_axis_name) +
      labs(color = "Cell Type") +
      ggtitle(factor_name) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position="bottom")
    
    # if plotting r-squared change y-limits to 0-1
    if (stat_type == 'adj_rsq') {
      p <- p + ylim(c(-.1,1))
    }
    
    # if plotting -log10 pvals draw significance line
    if (stat_type == 'adj_pval') {
      p <- p + geom_hline(yintercept=-log10(.01), linetype="dashed", color = "red")
    }
    
    legend <- cowplot::get_legend(
      p + theme(legend.box.margin = margin(0, 0, 30, 0))
    )
    
    p <- p + theme(legend.position="none")
    
    plot_list[[factor_name]] <- p
    
  }
  
  fig <- cowplot::plot_grid(plotlist=plot_list, ncol=n_col)
  
  fig <- cowplot::plot_grid(fig, legend, ncol = 1, rel_heights = c(1, .1))
  
  return(fig)
}

pv_table <- pv_table[pv_table$ctype=='Th',]
reg_stat_plots <- plot_subclust_associations(pv_table,n_col=2) # to generate the plot

pdf(file = "/home/jmitchel/figures/for_paper_v2/Th_multi_res2.pdf", useDingbats = FALSE,
    width = 5, height = 6.25)
reg_stat_plots
dev.off()











##### now plotting some individual ctype associations
# Note that the p-values on these plots are using fdr correction for each ctype batch of plots.
# The better way to do it is to use the fdr correction that was done when all association tests were initially run.
# These pvalues are found in pbmc_container$subc_factor_association_res and are the ones that appear in the final figure.

pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='Th',res=.6,n_col=2,alt_name='T4')
pdf(file = "/home/jmitchel/figures/for_paper_v2/Th_subc2.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()

pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='cMono',res=.5,n_col=2,alt_name='cM')
pdf(file = "/home/jmitchel/figures/for_paper_v2/cMono_subc2.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()

pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='NK',res=.6,n_col=2)
pdf(file = "/home/jmitchel/figures/for_paper_v2/NK_subc2.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()

pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='Tc',res=.6,n_col=2,alt_name='T8')
pdf(file = "/home/jmitchel/figures/for_paper_v2/Tc_subc2.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()


pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='ncMono',res=.6,n_col=2,alt_name='ncM')
pdf(file = "/home/jmitchel/figures/for_paper_v2/ncMono_subc2.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()

pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='cDC',res=.5,n_col=2)
pdf(file = "/home/jmitchel/figures/for_paper_v2/cDC_subc2.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()

pbmc_container <- get_ctype_subc_prop_associations(pbmc_container,ctype='B',res=.8,n_col=2)
pdf(file = "/home/jmitchel/figures/for_paper_v2/B_subc2.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()


## now getting full ctype prop associations
pbmc_container <- get_ctype_prop_associations(pbmc_container,'adj_pval',n_col=2)
pdf(file = "/home/jmitchel/figures/for_paper_v2/major_ctype_props2.pdf", useDingbats = FALSE,
    width = 8, height = 8)
pbmc_container$plots$ctype_prop_factor_associations
dev.off()

# getting an example dotplot to show for demonstrating the process
myplot <- get_subclust_enr_dotplot(pbmc_container,'T4',0.6,subtype=1,factor_use=1,ctype_cur='Th')
pdf(file = "/home/jmitchel/figures/for_paper_v2/CD4_f1_sub1_dot2.pdf", useDingbats = FALSE,
    width = 5.25, height = 3.5)
myplot
dev.off()






## plotting F6 (sex-associated factor) associations with cell proportion shifts
pv_table <- pbmc_container$subc_factor_association_res
pv_table <- pv_table[pv_table$factor=='Factor 6',]
print(pv_table)
f_6_subc_pv <- c(.175,.112,.008,.794,.157,.72,.040,.45)
f_6_subc_ct <- c('NK','Th','Tc','ncMono','cMono','cDC','B','Major ctypes')
f_6_ct_subc <- cbind.data.frame(-log10(f_6_subc_pv),f_6_subc_ct)
colnames(f_6_ct_subc) <- c('adj_pval','ctype')
f_6_ct_subc$ctype <- factor(f_6_ct_subc$ctype,levels=f_6_subc_ct)

p <- ggplot(f_6_ct_subc, aes(x = ctype, y = adj_pval)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(.001), 
             color = "red", size=1.5) +
  ylab("-log10(padj)") +
  xlab('Subclusters from this cell-type') +
  ggtitle('Factor 6\nsubcluster shift associations') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

pdf(file = "/home/jmitchel/figures/for_paper_v2/f6_subc_all2.pdf", useDingbats = FALSE,
    width = 5, height = 3.5)
p
dev.off()











##### render cell subtype marker dotplots
## using full dataset for these since clusters were generated using all donors
pbmc <- readRDS('/home/jmitchel/data/lupus_data/lupus_subsetted_seurat_v3.rds')

# read in the subclusters identified above
subclusters <- readRDS(file='/home/jmitchel/data/lupus_data/lupus_subcluster_data.rds')


# T8
subclusts <- subclusters[["T8"]][["res:0.6"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('T8','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('GZMH','FGFBP2','GNLY','GZMB','NKG7',
            'GZMA','CCL5','CST7','LGALS1','KLRD1',
            'DUSP2','GZMK','LYAR','CMC1','PIK3R1',
            'CD69','ZFP36L2','CD8B','NOSIP','CCR7',
            'EIF3E','C6orf48','PRKCQ-AS1','LTB','LDHB',
            'SELL')
t8_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()

# T4
subclusts <- subclusters[["T4"]][["res:0.6"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('T4','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('CCR7','FHIT','AIF1','NUCB2','SELL','PRKCQ-AS1','LINC00861','CD7',
            'RGS10','LEF1','ITGB1','S100A4','ANXA1','S100A10','S100A11','LGALS1',
            'KLRB1','EMP3','SH3BGRL3','CRIP1','JUN','JUNB','CD69','FOS','DUSP1',
            'CXCR4','TSC22D3','LEPROTL1','ZFP36L2','PNRC1','IFI6','ISG15','IFI44L',
            'LY6E','MT2A','MX1','IFIT3','IFITM1','EIF2AK2','ISG20','NEAT1',
            'JUND','MT-ND4L','SYNE2','ANKRD11','ETS1','DDX17','NKTR','MT-ND5','C1QA',
            'C1QB','COX5B','PLAC8','PSAP')
t4_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()


# NK
subclusts <- subclusters[["NK"]][["res:0.6"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('NK','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('FGFBP2','SPON2','GZMH','FCGR3A','LGALS1','PRF1','LAIR2','CCL4',
            'GZMB','IGFBP7','CMC1','CXCR4','DUSP2','CD160','PIK3R1','XCL2',
            'MAP3K8','CCL3','ZFP36','JUNB','GZMK','XCL1','SELL','KLRC1','CD44',
            'NFKBIA','COTL1','C1orf56','CDC42SE1','HNRNPH1','APOBEC3C','MDM4',
            'B4GALT1','CDC42','CTNNB1','SAR1A','SET')
nk_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()



# cM
subclusts <- subclusters[["res:0.5"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('cM','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('IFI6','ISG15','LY6E','IFI44L','IFITM3','IFI44','MX1','MT2A',
            'S100A12','EPSTI1','HLA-DPB1','HLA-DPA1','HLA-DMA','HLA-DQB1',
            'HLA-DQA1','LGALS2','CPVL','HLA-DRB1','SLC25A5','EEF1B2','ALOX5AP',
            'MGST1','RBP7','VCAN','PLBD1','CDA','METTL9','RETN','VNN2',
            'APOBEC3A','MARCKS','PSME2','IL8','IL1B','G0S2','CCL3','EREG',
            'TNFAIP3','NFKBIZ','SOD2','IER3','NAMPT')
cM_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()


# ncM
subclusts <- subclusters[["ncM"]][["res:0.6"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('ncM','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('TESC','VMO1','LYPD2','PPM1N','CKB','CDKN1C','ICAM4','SOD1','MEG3',
            'IFI6','ISG15','IFI44L','APOBEC3A','EPSTI1','IFIT3','TNFSF10',
            'PLSCR1','PLAC8','MX1','C1QA','C1QB','C1QC','HLA-DQA1','HLA-DMA',
            'VAMP8','HLA-DQB1','VAMP5','HLA-DMB','MS4A6A','JUND','LGALS2','G0S2',
            'VCAN','GPX1','MT-ND4L')
ncM_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()



# cDC
subclusts <- subclusters[["cDC"]][["res:0.5"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('cDC','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('VCAN','FCN1','S100A12','CD14','RNASE2','MS4A6A','CFD','CD36',
            'S100A8','SERPINA1','CD1C','FCER1A','ENHO','IL2RG','NDRG2','BASP1',
            'FCGR2B','CD1E','ARL4C','CLEC10A','CLEC9A','C1orf54','IRF8',
            'DNASE1L3','CPNE3','WDFY4','HLA-DOB','IDO1','GYPC')
cDC_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()


# B
subclusts <- subclusters[["B"]][["res:0.8"]]
pbmc_sub <- subset(pbmc, cells=names(subclusts))
identical(rownames(pbmc_sub@meta.data),names(subclusts))
new_cg <- sapply(subclusts,function(x) {
  return(paste0('B','_',x))
})
pbmc_sub@meta.data$cg_cov <- new_cg
Idents(pbmc_sub) <- pbmc_sub@meta.data$cg_cov
sub_lev <- levels(Idents(pbmc_sub))
Idents(pbmc_sub) <- factor(Idents(pbmc_sub),levels=sub_lev[order(sub_lev)])
myfeat <- c('GPR183','CRIP1','COTL1','EMP3','CLECL1','LGALS1','GAPDH','S100A6',
            'S100A10','IL4R','PPAPDC1B','FCER2','MEF2C','ADAM28','TCL1A','BIRC3',
            'STAG3','HVCN1','VPREB3','CD79B','IGLL5','SNX29P2','FCRLA','PPP1R14A',
            'MZB1','CD72','IER2','CD69','FOS','FOSB','JUN','NFKBIA','ZFP36','JUNB',
            'YPEL5','CD83','C1orf56','HNRNPH1','CDC42SE1','MDM4','B4GALT1',
            'PPP3CA','APOBEC3C','CDC42','CTNNB1','FGD2')
b_dot <- DotPlot(pbmc_sub, features = myfeat, cluster.idents = FALSE) + coord_flip()



























