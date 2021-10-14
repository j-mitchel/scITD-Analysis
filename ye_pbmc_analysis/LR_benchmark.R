library(spqn) 

test_lr_validity <- function(container,lr_test,ligand_target_matrix) {
  # all_pvals <- c()
  print(length(lr_test))
  all_pvals <- sccore::plapply(1:length(lr_test), function(lst_ndx) {
    lig_source_target <- lr_test[lst_ndx]
    lst_split <- strsplit(lig_source_target,split='_')[[1]]
    lig <- lst_split[[1]]
    source_ct <- lst_split[[2]]
    target_ct <- lst_split[[3]]
    
    lig_ct_exp <- container$scale_pb_extra[[source_ct]][,lig]
    
    all_cors <- cor_helper(container,lig_ct_exp,lig,source_ct,target_ct)
    
    # testing using wilcoxin test
    if (lig %in% colnames(ligand_target_matrix)) {
      gene_intersect <- intersect(names(all_cors),rownames(ligand_target_matrix))
      all_cors <- all_cors[gene_intersect]
      
      top_cors <- all_cors[order(abs(all_cors),decreasing=TRUE)][1:200]
      tmp <- cbind.data.frame(all_cors,ligand_target_matrix[gene_intersect,lig])
      colnames(tmp) <- c('mycors','target_scores')
      tmp$top <- sapply(1:nrow(tmp),function(i){
        if (rownames(tmp)[i] %in% names(top_cors)) {
          return('top_gn')
        } else {
          return('other')
        }
      })
      tmp$top <- as.factor(tmp$top)
      # my_res <- wilcox.test(target_scores~top,data=tmp)
      # pval <- my_res$p.value
      pval <- sig_test_helper(tmp,10000)
      return(pval)
    } else {
      return(NA)
    }
    
  },mc.preschedule=FALSE,n.cores=30,progress=TRUE)
  
  all_pvals2 <- all_pvals
  
  # add back names 
  names(all_pvals2) <- lr_test
  
  # remove na elements
  all_pvals2 <- all_pvals2[!is.na(all_pvals2)]
  
  # adjust pvalues
  all_pvals2 <- p.adjust(all_pvals2,method='fdr')
  
  # give warning if min non-zero adjusted pval not below .05
  if(min(all_pvals2[all_pvals2!=0])>.05) {
    print('Warning: minimum non-zero adjusted pvalue is > 0.05')
  }
  
  return(all_pvals2)
}


cor_helper <- function(container,lig_ct_exp,lig,source_ct,target_ct) {
  test <- cbind.data.frame(lig_ct_exp,container$scale_pb_extra[[target_ct]][names(lig_ct_exp),])
  cor_m <- cor(test)
  
  lig_ct_exp2 <- container$no_scale_pb_extra[[source_ct]][,lig]
  test2 <- cbind.data.frame(lig_ct_exp2,container$no_scale_pb_extra[[target_ct]][names(lig_ct_exp2),])
  ave_exp <- colMeans(test2)
  
  cor_m_spqn <- normalize_correlation(cor_m, ave_exp=ave_exp, ngrp=20, size_grp=300, ref_grp=5)
  all_cors <- cor_m_spqn[1,]
  names(all_cors) <- colnames(cor_m)
  all_cors <- all_cors[2:length(all_cors)]
  return(all_cors)
}

sig_test_helper <- function(tmp,n_iter) {
  
  ## to use with wilcoxon version
  my_res <- wilcox.test(target_scores~top,data=tmp)
  obs_pval <- my_res$p.value
  
  all_pv <- c()
  for (i in 1:n_iter) {
    tmp2 <- transform(tmp, top = sample(top) )
    my_res <- wilcox.test(target_scores~top,data=tmp2)
    pval <- my_res$p.value
    all_pv <- c(all_pv,pval)
  }
  pval <- sum(all_pv<=obs_pval) / length(all_pv)
  
  return(pval)
}

ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds")) # ligand-target NicheNet regulatory potential scores

scITD_channels <- readRDS(file='/home/jmitchel/data/lupus_data/LR_scITD_res.rds')

scITD_res <- test_lr_validity(pbmc_container,scITD_channels,ligand_target_matrix)
scITD_frac <- sum(scITD_res<.05)/length(scITD_res)

## store results for plotting later
res_dat <- matrix(ncol=5,nrow=50)
colnames(res_dat) <- c('db','method','n_chan','n_chan_eval','percent_enr')
res_dat[1,] <-  c('CellChat','scITD',length(scITD_channels),length(scITD_res),scITD_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')

# making barplot for full scITD results on CellChat lr pairs
tmp <- cbind.data.frame(names(scITD_res),scITD_res)
colnames(tmp) <- c('channel','adj_pval')
tmp <- tmp[order(tmp$adj_pval,decreasing=FALSE),]
tmp$adj_pval[tmp$adj_pva==0] <- .0001
tmp$channel <- factor(tmp$channel,levels=tmp$channel)
tmp$adj_pval <- -log10(tmp$adj_pval)
mylabels <- tmp$channel
print(which(mylabels=='ICOSLG_cMono_Th'))
print(which(mylabels=='TNFSF13B_cMono_B'))
ndx_keep <- c(1,4,8,12,16,21,26)
mylabels2 <- rep('',length(mylabels))
mylabels2[ndx_keep] <- as.character(mylabels)[ndx_keep]
chan_bplot <- ggplot(tmp,aes(x=channel,y=adj_pval)) +
  geom_bar(stat='identity') +
  scale_x_discrete(labels= mylabels2) +
  scale_y_continuous(limits = c(0,4.25), expand = c(0, 0)) +
  ylab('Top target genes-NicheNet\nenrichment -log10(padj)') +
  xlab('scITD inferred channels (ligand_source_target)') +
  geom_hline(yintercept = -log10(.05), 
             color = "red", size=.75) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = -30, vjust = 0.5, hjust=0),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

pdf(file = "/home/jmitchel/figures/for_paper_v2/LR_chan_bplot3.pdf", useDingbats = FALSE,
    width = 9, height = 4)
chan_bplot
dev.off()



cc_channels <- readRDS(file='/home/jmitchel/data/lupus_data/LR_cchat_res_orig.rds')
cc_res <- test_lr_validity(pbmc_container,cc_channels,ligand_target_matrix)
cc_frac <- sum(cc_res<.05)/length(cc_res)
res_dat[2,] <-  c('CellChat','CellChat',length(cc_channels),length(cc_res),cc_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')

cchat_assoc_channels <- get_LR_association_chan(pbmc_container,lr_pairs)
cchat_assoc_res <- test_lr_validity(pbmc_container,cchat_assoc_channels,ligand_target_matrix)
cchat_assoc_frac <- sum(cchat_assoc_res<.05)/length(cchat_assoc_res)
res_dat[3,] <-  c('CellChat','LR_association',length(cchat_assoc_channels),length(cchat_assoc_res),cchat_assoc_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')

rand_channels <- get_rand_chan(pbmc_container,lr_pairs,n_samp=250)
rand_res <- test_lr_validity(pbmc_container,rand_channels,ligand_target_matrix)
rand_frac <- sum(rand_res<.05)/length(rand_res)
res_dat[4,] <-  c('CellChat','random',length(rand_channels),length(rand_res),rand_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')




### now for iTALK db analysis
scITD_channels <- readRDS(file='/home/jmitchel/data/lupus_data/LR_iTALK_scITD_res.rds')
scITD_res <- test_lr_validity(pbmc_container,scITD_channels,ligand_target_matrix)
scITD_frac <- sum(scITD_res<.05)/length(scITD_res)
res_dat[5,] <-  c('iTALK','scITD',length(scITD_channels),length(scITD_res),scITD_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')





italk_channels <- readRDS(file='/home/jmitchel/data/lupus_data/LR_iTALK_res_parsed.rds')
italk_channels <- italk_channels[1:600]
italk_channels_sub <- italk_channels[1:200]
iTALK_res <- test_lr_validity(pbmc_container,italk_channels,ligand_target_matrix)
iTALK_frac <- sum(iTALK_res<.05)/length(iTALK_res)
# res_dat[6,] <-  c('iTALK','iTALK',length(italk_channels),length(iTALK_res),iTALK_frac)


# ## looking at HLA downstream targets in NicheNet
# test1 <- ligand_target_matrix[,'HLA-C']
# test2 <- ligand_target_matrix[,'HLA-B']
# plot(test1,test2)
# 
# ## remove the promiscuous interactions
# ligs <- sapply(names(iTALK_res), function(x) {
#   my_h <- strsplit(x,split='_')[[1]][[1]]
# })
# table(ligs)
# ligs_keep <- names(table(ligs))[table(ligs)<=3]
# chan_keep <- names(ligs)[ligs %in% ligs_keep]
# iTALK_res2 <- iTALK_res[chan_keep]
# sum(iTALK_res2<.05)/length(iTALK_res2)
# 
# ligs <- sapply(names(scITD_res), function(x) {
#   my_h <- strsplit(x,split='_')[[1]][[1]]
# })
# table(ligs)
# ligs_keep <- names(table(ligs))[table(ligs)<=3]
# chan_keep <- names(ligs)[ligs %in% ligs_keep]
# scITD_res2 <- scITD_res[chan_keep]
# sum(scITD_res2<.05)/length(scITD_res2)
# 
# is_hla <- sapply(names(iTALK_res_3),function(x){
#   my_h <- strsplit(x,split='-')[[1]][[1]]
#   if (my_h=='HLA') {
#     return(TRUE)
#   } else {
#     return(FALSE)
#   }
# })
# iTALK_res2 <- iTALK_res_3[!is_hla]
# sum(iTALK_res2<.05)/length(iTALK_res2)



italk_assoc_channels <- get_LR_association_chan(pbmc_container,lr_pairs)
italk_assoc_res <- test_lr_validity(pbmc_container,italk_assoc_channels,ligand_target_matrix)
italk_assoc_frac <- sum(italk_assoc_res<.05)/length(italk_assoc_res)
res_dat[6,] <-  c('iTALK','LR_association',length(italk_assoc_channels),
                  length(italk_assoc_res),italk_assoc_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')


rand_channels <- get_rand_chan(pbmc_container,lr_pairs,n_samp=250)
rand_res <- test_lr_validity(pbmc_container,rand_channels,ligand_target_matrix)
rand_frac <- sum(rand_res<.05)/length(rand_res)
res_dat[7,] <-  c('iTALK','random',length(rand_channels),length(rand_res),rand_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')


res_dat <- readRDS(file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')


scITD_channels <- readRDS(file='/home/jmitchel/data/lupus_data/LR_singlecellsig_scITD_res.rds')
scITD_res <- test_lr_validity(pbmc_container,scITD_channels,ligand_target_matrix)
scITD_frac <- sum(scITD_res<.05)/length(scITD_res)
res_dat[8,] <-  c('SingleCellSignalR','scITD',length(scITD_channels),length(scITD_res),scITD_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')


scs_assoc_channels <- get_LR_association_chan(pbmc_container,lr_pairs)
scs_assoc_res <- test_lr_validity(pbmc_container,scs_assoc_channels,ligand_target_matrix)
scs_assoc_frac <- sum(scs_assoc_res<.05)/length(scs_assoc_res)
res_dat[9,] <-  c('SingleCellSignalR','LR_association',length(scs_assoc_channels),
                  length(scs_assoc_res),scs_assoc_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')


rand_channels <- get_rand_chan(pbmc_container,lr_pairs,n_samp=250)
rand_res <- test_lr_validity(pbmc_container,rand_channels,ligand_target_matrix)
rand_frac <- sum(rand_res<.05)/length(rand_res)
res_dat[10,] <-  c('SingleCellSignalR','random',length(rand_channels),length(rand_res),rand_frac)
# saveRDS(res_dat,file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')


res_dat <- readRDS(file='/home/jmitchel/data/lupus_data/LR_method_comparison_v2.rds')
res_dat <- res_dat[c(1,3,4,5,6,7,8,9,10),]

res_dat <- as.data.frame(res_dat)
res_dat$n_chan <- as.numeric(res_dat$n_chan)
res_dat$n_chan_eval <- as.numeric(res_dat$n_chan_eval)
res_dat$percent_enr <- as.numeric(res_dat$percent_enr)
res_dat$db <- as.factor(res_dat$db)
res_dat$method <- factor(res_dat$method,levels=c('random','LR_association','scITD'))

## adding in se error
res_dat$se_dat <- sqrt(res_dat$percent_enr*(1-res_dat$percent_enr)/res_dat$n_chan_eval)

library(RColorBrewer)
mycol = brewer.pal(n = 8, name = "Dark2")

comp_bplot_size <- ggplot(res_dat,aes(x=db,y=n_chan,fill=method)) +
  geom_bar(stat='identity', position="dodge") +
  ylab('Number of inferred interactions') +
  xlab('LR database') +
  guides(fill=guide_legend(title="LR method")) +
  scale_fill_manual(values=c(mycol[8], mycol[1], mycol[6])) +
  theme_bw() +
  theme(legend.position = "none")

comp_bplot_size

comp_bplot_frac <- ggplot(res_dat,aes(x=db,y=percent_enr,fill=method)) +
  geom_bar(stat='identity', position="dodge") +
  geom_errorbar(aes(ymin = percent_enr - se_dat, ymax = percent_enr + se_dat), 
                width=0.15, position=position_dodge(.9)) +
  ylab('Fraction of inferred interactions\nwith NicheNet enrichment') +
  xlab('LR database') +
  guides(fill=guide_legend(title="LR method")) +
  scale_fill_manual(values=c(mycol[8], mycol[1], mycol[6])) +
  theme_bw()

comp_bplot_frac

combined_plt <- cowplot::plot_grid(comp_bplot_size,comp_bplot_frac,nrow=1, rel_widths = c(.7,1.1))

pdf(file = "/home/jmitchel/figures/for_paper_v2/LR_all_db_bplot3.pdf", useDingbats = FALSE,
    width = 7.5, height = 3.5)
combined_plt
dev.off()




















