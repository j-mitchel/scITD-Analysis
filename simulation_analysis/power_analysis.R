library(ggplot2)

num_variants <- c(200000,500000,1000000,5000000,10000000)
MAFs <- c(.01,.02,.03,.04,.05,.075,.10,.12,.15,.20,.30,.40)
power_vals <- c()
MAF_vals <- c()
num_var_vals <- c()
for (num_var in num_variants) {
  print(num_var)
  for (MAF in MAFs) {
    print(MAF)
    sig_thresh <- .05/num_var
    samp_size <- 171
    num_alleles <- samp_size * 2
    num_minor <- ceiling(num_alleles * MAF)
    num_major <- num_alleles - num_minor
    
    trials <- 1000
    num_sig <- 0
    for (k in 1:trials) {
      all_alleles <- c(rep(1,num_major),rep(2,num_minor))
      AA <- 0
      AB <- 0
      BB <- 0
      for (i in 1:samp_size) {
        # select alleles for the individual
        d_allele <- sample(all_alleles,2,replace = FALSE)
        
        # add to count of genotype
        if (d_allele[1]==1 && d_allele[2]==1) {
          AA <- AA + 1
        } else if (d_allele[1]==2 && d_allele[2]==2){
          BB <- BB + 1
        } else {
          AB <- AB + 1
        }
        
        # remove those alleles from the population
        for (j in d_allele) {
          ndx_match <- which(all_alleles==j)
          ndx_rem <- ndx_match[1]
          ndx_keep <- c(1:length(all_alleles)!=ndx_rem)
          all_alleles <- all_alleles[ndx_keep]
        }
      }
      
      y1 <- rnorm(AA,.1,.075)
      y2 <- rnorm(AB,0,.075)
      y3 <- rnorm(BB,-.1,.075)
      
      dscores <- c(y1,y2,y3)
      x_val <- c(rep(0,length(y1)),rep(1,length(y2)),rep(2,length(y3)))
      tmp <- as.data.frame(cbind(dscores,x_val))
      colnames(tmp) <- c('dscores','minor_allele')
      
      lmres <- lm(dscores~minor_allele,data=tmp)
      lmres <- summary(lmres)
      pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
      
      sig <- pval < sig_thresh
      if (sig) {
        num_sig <- num_sig + 1
      }
    }
    
    power <- num_sig/trials
    
    # store results
    power_vals <- c(power_vals,power)
    MAF_vals <- c(MAF_vals,MAF)
    num_var_vals <- c(num_var_vals,num_var)
  }
}

tmp <- as.data.frame(cbind(power_vals,MAF_vals,num_var_vals))
tmp$MAF_vals <- tmp$MAF_vals * 100
ggplot(tmp,aes(x=MAF_vals,y=power_vals,color=as.factor(num_var_vals))) +
  geom_point() +
  geom_line() +
  xlab('MAF (percent)') +
  ylab('Power') +
  labs(color = "Num variants tested") +
  scale_x_continuous(trans='log2',breaks = c(1,2,seq(0, 30, by = 5)))
  

# get best fit line for single example plot
lmres <- lm(dscores~minor_allele,data=tmp)
lmres <- summary(lmres)
xrange <- seq(0,2,.1)
y_line <- (lmres$coefficients[2,1] * xrange) + lmres$coefficients[1,1]
myline <- as.data.frame(cbind(y_line,xrange))
colnames(myline) <- c('yval','xval')

## for single plot
# mysummary <- as.data.frame(cbind(c(mean(y1),mean(y2),mean(y3)),c(0,1,2)))
# colnames(mysummary) <- c('mean_val','x_val')
ggplot(tmp,aes(x=minor_allele,y=dscores)) +
  geom_point() +
  xlab('Minor allele count') +
  ylab('Factor X donor score') +
  geom_line(data=myline, aes(x=xval,y=yval)) +
  ggtitle('Example with MAF = 15%') +
  theme(plot.title = element_text(hjust = 0.5))
  # geom_line(data=mysummary, aes(x=x_val,y=mean_val,color='red'))


