library(scITD)

sim_data <- get_sim_data(donors_total=50, cells_per_donor=250, n_processes=2,
                         donors_per_process=10, n_ctypes=2, n_genes=500, 
                         de_prob=0.07, de_strength=2, factor_overlap=TRUE, 
                         rseed=97854658)

sim_scMinimal <- instantiate_scMinimal(data_sparse = as.matrix(sim_data[[1]]), meta_data = sim_data[[2]])
sim_container <- make_new_container(sim_scMinimal,
                                    ctypes_use = c('ct1', 'ct2'),
                                    scale_var = T, var_scale_power = 1.25,
                                    tucker_type = 'sparse', 
                                    rotation_type = 'varimax',
                                    ncores = 30, rand_seed = 16)

# get ctype matrices
sim_container <- get_ctype_data(sim_container,make_clean=TRUE)


# show that rank determination method selects correct number of factors
sim_container <- determine_ranks_tucker(sim_container, max_ranks_test=c(6,10,5),
                                         method='svd', num_iter=5, shuffle_level='cells')
sim_container$plots$rank_determination_plot

# run tucker
sim_container <- run_tucker_ica(sim_container, ranks=c(2,4,2), shuffle=F)

# plot donor scores
sim_container <- plot_donor_matrix(sim_container, show_donor_ids = TRUE)
sim_container$plots$donor_matrix

# get significant genes
sim_container <- run_jackstraw(sim_container,n_fibers=20,n_iter=500)

# plot significant genes next to loadings and true de genes
sim_container <- get_all_lds_factor_plots(sim_container, use_sig_only=FALSE, 
                                           nonsig_to_zero=FALSE, 
                                           annot='sig_genes',
                                           sig_thresh=0.05, 
                                           sim_de_donor_group=c(2,1),
                                           display_genes=FALSE,
                                           gene_callouts=FALSE)
render_all_lds_plots(sim_container, n_rows=2)







# batch effect analysis
batch_donors1 <- c(1:5,11:15,21:35)
batch_donors1 <- sapply(batch_donors1,function(x){paste0('s',as.character(x))})
batch_donors2 <- c(6:10,16:20,36:50)
batch_donors2 <- sapply(batch_donors2,function(x){paste0('s',as.character(x))})

# batch effect analysis for unbalanced batches
batch_donors1 <- c(1:8,11,12,21:25)
batch_donors1 <- sapply(batch_donors1,function(x){paste0('s',as.character(x))})
batch_donors2 <- c(9:10,13:20,26:50)
batch_donors2 <- sapply(batch_donors2,function(x){paste0('s',as.character(x))})

# batch effect analysis for slightly less unbalanced batches
batch_donors1 <- c(1:10,21:35)
batch_donors1 <- sapply(batch_donors1,function(x){paste0('s',as.character(x))})
batch_donors2 <- c(11:20,36:50)
batch_donors2 <- sapply(batch_donors2,function(x){paste0('s',as.character(x))})

# batch effect analysis where a process is only present in two batches
batch_donors1 <- c(1:5,11:13,21:30)
batch_donors1 <- sapply(batch_donors1,function(x){paste0('s',as.character(x))})
batch_donors2 <- c(5:10,14:16,31:40)
batch_donors2 <- sapply(batch_donors2,function(x){paste0('s',as.character(x))})
batch_donors3 <- c(17:20,41:50)
batch_donors3 <- sapply(batch_donors3,function(x){paste0('s',as.character(x))})

# batch effect analysis where a process is only present in two batches
batch_donors1 <- c(1:5,21:30)
batch_donors1 <- sapply(batch_donors1,function(x){paste0('s',as.character(x))})
batch_donors2 <- c(5:10,11:15,31:40)
batch_donors2 <- sapply(batch_donors2,function(x){paste0('s',as.character(x))})
batch_donors3 <- c(16:20,41:50)
batch_donors3 <- sapply(batch_donors3,function(x){paste0('s',as.character(x))})

# what if we have one extra small factor thats all donors of one process
batch_donors1 <- c(1:4)
batch_donors1 <- sapply(batch_donors1,function(x){paste0('s',as.character(x))})
batch_donors2 <- c(5:7,11:15,21:35)
batch_donors2 <- sapply(batch_donors2,function(x){paste0('s',as.character(x))})
batch_donors3 <- c(8:10,16:20,36:50)
batch_donors3 <- sapply(batch_donors3,function(x){paste0('s',as.character(x))})

# what if process only in one batch but that batch also has donors from other processes
batch_donors1 <- c(1:10,11:14,21:30)
batch_donors1 <- sapply(batch_donors1,function(x){paste0('s',as.character(x))})
batch_donors2 <- c(15:17,31:40)
batch_donors2 <- sapply(batch_donors2,function(x){paste0('s',as.character(x))})
batch_donors3 <- c(18:20,41:50)
batch_donors3 <- sapply(batch_donors3,function(x){paste0('s',as.character(x))})

sim_data <- get_sim_data(donors_total=50, cells_per_donor=250, n_processes=2,
                         donors_per_process=10, n_ctypes=2, n_genes=500, 
                         de_prob=0.07, de_strength=2, factor_overlap=FALSE, 
                         add_batch=TRUE, 
                         batch_donors=list(batch_donors1,batch_donors2,batch_donors3),
                         rseed=50)
donors_total=50
cells_per_donor=250
n_processes=2
donors_per_process=10
n_ctypes=2
n_genes=500
de_prob=0.07
de_strength=2
factor_overlap=FALSE
add_batch=TRUE
batch_donors=list(batch_donors1,batch_donors2)
rseed=10

# add batch metadata
sim_data[[2]]$batch <- sapply(sim_data[[2]]$donors, function(x) {
  if (x %in% batch_donors1) {
    return('Batch1')
  } else {
    return('Batch2')
  }
})

# sim_scMinimal <- instantiate_scMinimal(data_sparse = as.matrix(sim_counts), meta_data = sim_meta)
sim_scMinimal <- instantiate_scMinimal(data_sparse = as.matrix(sim_data[[1]]), meta_data = sim_data[[2]])
sim_container <- make_new_container(sim_scMinimal,
                                    ctypes_use = c('ct1', 'ct2'),
                                    scale_var = T, var_scale_power = .5,
                                    tucker_type = 'sparse', 
                                    rotation_type = 'ica',
                                    ncores = 30, rand_seed = 16)


sim_container <- get_ctype_data(sim_container,make_clean=TRUE)

sim_container <- get_ctype_vargenes(sim_container, method="norm_var", thresh=250)
# sim_container <- collapse_by_donors(sim_container, shuffle=FALSE)

sim_container <- run_tucker_ica(sim_container, ranks=c(2,8,2), shuffle=F)
sim_container <- run_tucker_ica(sim_container, ranks=c(3,6,2), shuffle=F,
                                batch_var='batch')

sim_container <- plot_donor_matrix(sim_container, show_donor_ids = TRUE)
sim_container$plots$donor_matrix




# checking that donors 6-10 have different expression compared to 1-5 for same gene
test <- sim_container[["scMinimal_ctype"]][["ct1"]][["data_means"]]
testm <- sim_container[["scMinimal_ctype"]][["ct1"]][["metadata"]]
testm <- unique(testm)
rownames(testm) <- testm$donors
testm <- testm[rownames(test),]

tmp <- as.data.frame(cbind(test[,2],rownames(testm)))
colnames(tmp) <- c('expr','dnr')
as.numeric(tmp[tmp$dnr%in%c('s1','s2','s3','s4','s5'),]$expr)
as.numeric(tmp[tmp$dnr%in%c('s6','s7','s8','s9','s10'),]$expr)


# 2 and 13 have pretty big effect sizes
gene <- 2
tmp <- as.data.frame(cbind(as.numeric(batch_counts[gene,]),batch_meta$Batch))
colnames(tmp) <- c('expr','batch')
tmp$expr <- as.numeric(tmp$expr)
tmp$batch <- as.factor(tmp$batch)

# ggplot(tmp, aes(x=batch, y=expr)) + 
#   geom_boxplot()



tmp$dg <- sapply(batch_meta$donors, function(x) {
  if (x %in% c('s1','s2','s3','s4','s5')) {
    return('g1')
  } else if ((x %in% c('s6','s7','s8','s9','s10'))) {
    return('g2')
  } else {
    return('g3')
  }
})
tmp$dg <- as.factor(tmp$dg)

ggplot(tmp, aes(x=dg, y=expr)) + 
  geom_boxplot()







# need to compute by hand to see how mean values change when add in batch cells
test1 <- sim_container[["scMinimal_ctype"]][["ct1"]][["data_sparse"]]
m1 <- sim_container[["scMinimal_ctype"]][["ct1"]][["metadata"]]

dim(test1)
dim(m1)

batch_mask <- sapply(colnames(test1),function(x){
  print(x)
  tmp <- strsplit(x,split="_")
  if (length(tmp[[1]])==2) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})
batch_data <- test1[,batch_mask]
batch_meta <- m1[batch_mask,]
nonbatch_data <- test1[,!batch_mask]
nonbatch_meta <- m1[batch_mask,]

# hist(batch_data[2,batch_meta$donors=='s1'])
# hist(nonbatch_data[2,nonbatch_meta$donors=='s1'])
mean(batch_data[2,batch_meta$donors=='s1'])
mean(nonbatch_data[2,nonbatch_meta$donors=='s1'])
mean(test1[2,m1$donors=='s1'])

mean(batch_data[2,batch_meta$donors=='s9'])
mean(nonbatch_data[2,nonbatch_meta$donors=='s9'])
mean(test1[2,m1$donors=='s9'])

# what if center the batch cell expression around the non batch expression
new_batch <- scale(batch_data[2,],center=TRUE,scale=FALSE) + mean(nonbatch_data[2,])

bd1 <- batch_data[2,batch_meta$donors=='s1']
bd1 <- scale(bd1,center=TRUE,scale=FALSE)
bd1 <- bd1 + mean(nonbatch_data[2,nonbatch_meta$donors=='s1'])
mean(bd1)
mean(c(bd1,nonbatch_data[2,nonbatch_meta$donors=='s1']))

bd1 <- new_batch[batch_meta$donors=='s1']
mean(bd1)
bd2 <- new_batch[batch_meta$donors=='s9']
mean(bd2)










