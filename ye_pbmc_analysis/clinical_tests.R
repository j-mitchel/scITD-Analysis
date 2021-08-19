

library("readxl")

# get clinical data for eth info
clin_vars <- as.data.frame(read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data.xlsx'))
rownames(clin_vars) <- clin_vars$subjectid
eth = clin_vars[,'race',drop=FALSE]
eth = clin_vars[,'raceeth',drop=FALSE]
# eth = clin_vars[,'age',drop=FALSE]

# load data of categorical variables
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_categorical.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL

# get tucker donor scores to test
dsc <- pbmc_container$tucker_results[[1]]

## get donors in both dsc and in clin_vars
# trim donor IDs in dsc
trim_names <- sapply(rownames(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
old_names <- rownames(dsc)
names(old_names) <- trim_names
rownames(dsc) <- trim_names

# get donors in both dataframes
d_both <- rownames(clin_vars)[rownames(clin_vars) %in% rownames(dsc)]
old_id_both <- old_names[d_both]

# # get other meta data for this intersection of donors
# meta <- pbmc_container$scMinimal_ctype[[1]]$metadata
# meta <- meta[meta$donors%in%old_id_both,c('donors','Status')]
# meta <- unique(meta)
# rownames(meta) <- meta$donors

# ### limiting it to just patients with lupus nephritis
# clin_vars <- clin_vars[d_both,]
# d_both <- rownames(clin_vars)[clin_vars$crflupusarth==1]
# ###

# limit both dataframes to just the intersection of donors and in the same order
dsc <- dsc[d_both,]
clin_vars <- clin_vars[d_both,]
eth <- eth[d_both,,drop=FALSE]

#### apply other checks to clin vars before testing!!
ndx_rem <- c()
for (j in 1:ncol(clin_vars)) {
  d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]
  tmp <- clin_vars[d_keep,j]
  num_levs <- length(unique(tmp))
  if (num_levs==1) {
    ndx_rem <- c(ndx_rem,j)
  }
}
clin_vars <- clin_vars[,-ndx_rem]
####


# # do enrichment tests
# library(fgsea)
# mypaths <- list()
# for (j in 1:ncol(clin_vars)) {
#   d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]
#   d_in_set <- rownames(clin_vars)[clin_vars[d_keep,j]!=0]
#   mypaths[[colnames(clin_vars)[j]]] <- d_in_set
# }
# myranks <- dsc[,12]
# 
# # tmp_names <- names(mypaths)[c(1:20)]
# # mypaths <- mypaths[1:50]
# # names(mypaths) <- tmp_names
# 
# fgseaRes <- fgseaSimple(pathways = mypaths,
#                   stats    = myranks,
#                   minSize  = 3,
#                   maxSize  = 80,
#                   nperm=10000)
# fgseaRes <- fgseaRes[order(fgseaRes$padj,decreasing=FALSE),]
# head(fgseaRes)
# ##


all_pvals <- c()
f_tested <- c()
c_tested <- c()
# loop through the variables to test
for (j in 1:ncol(clin_vars)) {
  print(j)
  # loop through factors
  for (f in 1:ncol(dsc)) {
  # for (f in 1:4) {
    # get donors in clin var that don't have an NA value
    d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]

    tmp <- as.data.frame(cbind(dsc[d_keep,f], clin_vars[d_keep,j]))
    colnames(tmp) <- c('dscore','cvar')
    
    # # to try with eth covar
    # tmp <- as.data.frame(cbind(dsc[d_keep,f], clin_vars[d_keep,j], eth[d_keep,1]))
    # colnames(tmp) <- c('dscore','cvar','eth')
    # tmp$eth <- as.factor(tmp$eth)
    
    # force cvar to be factor
    tmp$cvar <- as.factor(tmp$cvar)

    # if smallest level has less thatn n donors skip this one
    if (min(table(tmp$cvar)) < 15) {
      next
    }

    # # compute linear model
    # lmres <- summary(lm(dscore ~ cvar, data=tmp))
    # pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)

    # # trying lm with regress out covariates
    # lm1 <- lm(dscore ~ eth, data=tmp)
    # lm2 <- lm(dscore ~ cvar + eth, data=tmp)
    # anova_res <- anova(lm1,lm2)
    # pval <- anova_res$`Pr(>F)`[2]
    
    
    # # trying with two sample t-test
    # t_res <- try(stats::t.test(dscore ~ cvar, data=tmp,
    #                        alternative = "two.sided",var.equal = FALSE))
    # if (class(t_res) == 'try-error') {
    #   next
    # }
    # pval <- t_res$p.value


    # trying with logistic regression model
    # fmod <- glm(cvar~dscore+eth, data=tmp, family = "binomial") ##"full" mod
    # nmod <- glm(cvar~1+eth, data=tmp, family = 'binomial') ##"null" mod
    fmod <- glm(cvar~dscore, data=tmp, family = "binomial") ##"full" mod
    nmod <- glm(cvar~1, data=tmp, family = 'binomial') ##"null" mod
    a_res <- anova(nmod, fmod, test = 'Chisq')
    pval <- a_res$`Pr(>Chi)`[2]

    all_pvals <- c(all_pvals,pval)
    f_tested <- c(f_tested,f)
    c_tested <- c(c_tested,colnames(clin_vars)[j])

  }
}
all_pvals <- p.adjust(all_pvals,method='fdr')
all_pvals[order(all_pvals,decreasing=FALSE)][1:40]
f_tested[order(all_pvals,decreasing=FALSE)][1:40]
c_tested[order(all_pvals,decreasing=FALSE)][1:40]


### plot antibody hits for f1 sle only acrantidsdna
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"acrantidsdna"])]
tmp <- as.data.frame(cbind(dsc[d_keep,1], clin_vars[d_keep,"acrantidsdna"]))
colnames(tmp) <- c('dscore','cvar')

# force cvar to be factor
tmp$cvar <- as.factor(tmp$cvar)

tmp$cvar_word <- sapply(tmp$cvar,function(x){
  if (x==1) {
    return('anti-dsDNA')
  } else {
    return('No anti-dsDNA')
  }
})
tmp$cvar_word <- as.factor(tmp$cvar_word)

# pdf(file = "/home/jmitchel/figures/for_paper/sle_only_antidsdna.pdf", useDingbats = FALSE,
#     width = 3.75, height = 2.75)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 1 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
dev.off()

### plot antibody hits for f1 sle only acrantismith
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"acrantismith"])]
tmp <- as.data.frame(cbind(dsc[d_keep,1], clin_vars[d_keep,"acrantismith"]))
colnames(tmp) <- c('dscore','cvar')

# force cvar to be factor
tmp$cvar <- as.factor(tmp$cvar)

tmp$cvar_word <- sapply(tmp$cvar,function(x){
  if (x==1) {
    return('anti-smith')
  } else {
    return('No anti-smith')
  }
})
tmp$cvar_word <- as.factor(tmp$cvar_word)

pdf(file = "/home/jmitchel/figures/for_paper/sle_only_antismith.pdf", useDingbats = FALSE,
    width = 3.75, height = 2.75)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 1 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
dev.off()



# plotting dscores against ethnicity
myfactor=2
# tmp <- as.data.frame(cbind(dsc[d_keep,6], eth[d_keep,1]))
tmp <- as.data.frame(cbind(dsc[d_keep,myfactor], eth[d_keep,1]))
colnames(tmp) <- c('dscore','eth')
tmp$eth <- as.factor(tmp$eth)

fmod <- glm(eth~dscore, data=tmp, family = "binomial") ##"full" mod
nmod <- glm(eth~1, data=tmp, family = 'binomial') ##"null" mod
a_res <- anova(nmod, fmod, test = 'Chisq')
pval <- a_res$`Pr(>Chi)`[2]
pval

ggplot(tmp,aes(x=eth,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 1 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()



d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"crfmucuulcers"])]
tmp <- as.data.frame(cbind(dsc[d_keep,5], clin_vars[d_keep,"crfmucuulcers"]))
colnames(tmp) <- c('dscore','cvar')

d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"crfothrashsle"])]
tmp <- as.data.frame(cbind(dsc[d_keep,5], clin_vars[d_keep,"crfothrashsle"]))
colnames(tmp) <- c('dscore','cvar')

d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"crflymphadeno"])]
tmp <- as.data.frame(cbind(dsc[d_keep,6], clin_vars[d_keep,"crflymphadeno"]))
colnames(tmp) <- c('dscore','cvar')

### plot lupusneph result
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"crflupusneph"])]
tmp <- as.data.frame(cbind(dsc[d_keep,2], clin_vars[d_keep,"crflupusneph"]))
colnames(tmp) <- c('dscore','cvar')


# force cvar to be factor
tmp$cvar <- as.factor(tmp$cvar)

tmp$cvar_word <- sapply(tmp$cvar,function(x){
  if (x==1) {
    return('Lupus Nephritis')
  } else {
    return('No Lupus Nephritis')
  }
})
tmp$cvar_word <- as.factor(tmp$cvar_word)

# pdf(file = "/home/jmitchel/figures/for_paper/LN_factor2.pdf", useDingbats = FALSE,
#     width = 4.5, height = 3.5)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 2 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
dev.off()

fmod <- glm(cvar~dscore, data=tmp, family = "binomial") ##"full" mod
nmod <- glm(cvar~1, data=tmp, family = 'binomial') ##"null" mod
a_res <- anova(nmod, fmod, test = 'Chisq')
pval <- a_res$`Pr(>Chi)`[2]


### plot acrantidsdna result
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"acrantidsdna"])]
tmp <- as.data.frame(cbind(dsc[d_keep,2], clin_vars[d_keep,"acrantidsdna"]))
colnames(tmp) <- c('dscore','cvar')

# force cvar to be factor
tmp$cvar <- as.factor(tmp$cvar)

tmp$cvar_word <- sapply(tmp$cvar,function(x){
  if (x==1) {
    return('Anti-dsDNA Present')
  } else {
    return('No Anti-dsDNA Present')
  }
})
tmp$cvar_word <- as.factor(tmp$cvar_word)

pdf(file = "/home/jmitchel/figures/for_paper/dsDNA_factor2.pdf", useDingbats = FALSE,
    width = 4.5, height = 3.5)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 2 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
dev.off()

fmod <- glm(cvar~dscore, data=tmp, family = "binomial") ##"full" mod
nmod <- glm(cvar~1, data=tmp, family = 'binomial') ##"null" mod
a_res <- anova(nmod, fmod, test = 'Chisq')
pval <- a_res$`Pr(>Chi)`[2]



### plot crflymphadeno result
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"crflymphadeno"])]
tmp <- as.data.frame(cbind(dsc[d_keep,2], clin_vars[d_keep,"crflymphadeno"]))
colnames(tmp) <- c('dscore','cvar')

# force cvar to be factor
tmp$cvar <- as.factor(tmp$cvar)

tmp$cvar_word <- sapply(tmp$cvar,function(x){
  if (x==1) {
    return('Lymphadenopathy')
  } else {
    return('No Lymphadenopathy')
  }
})
tmp$cvar_word <- as.factor(tmp$cvar_word)

pdf(file = "/home/jmitchel/figures/for_paper/lymph_factor2.pdf", useDingbats = FALSE,
    width = 4.5, height = 3.5)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 2 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
dev.off()

fmod <- glm(cvar~dscore, data=tmp, family = "binomial") ##"full" mod
nmod <- glm(cvar~1, data=tmp, family = 'binomial') ##"null" mod
a_res <- anova(nmod, fmod, test = 'Chisq')
pval <- a_res$`Pr(>Chi)`[2]











# load data of ordinal variables
library(MASS)
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data_ordinal.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL

clin_vars$sliccmalignancy <- NULL
clin_vars$lupusseverityindex <- NULL
clin_vars$smokestat <- NULL
clin_vars$acrcsum <- NULL
clin_vars$sliccavasnec <- NULL
clin_vars$slicccva <- NULL

clin_vars <- clin_vars[d_both,]

# all_pvals <- c()
# f_tested <- c()
# c_tested <- c()
# # loop through the variables to test
# for (j in 1:ncol(clin_vars)) {
#   print(j)
#   # loop through factors
#   for (f in 1:ncol(dsc)) {
#     # get donors in clin var that don't have an NA value
#     d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]
# 
#     # tmp <- as.data.frame(cbind(dsc[d_keep,f], clin_vars[d_keep,j]))
#     # colnames(tmp) <- c('dscore','cvar')
#     
#     # to try with eth covar
#     tmp <- as.data.frame(cbind(dsc[d_keep,f], clin_vars[d_keep,j], eth[d_keep,1]))
#     colnames(tmp) <- c('dscore','cvar','eth')
#     tmp$eth <- as.factor(tmp$eth)
# 
#     # compute linear model
#     # force cvar to be numeric
#     tmp$cvar <- as.numeric(tmp$cvar)
# 
#     # lmres <- summary(lm(dscore ~ cvar, data=tmp))
#     # pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
# 
#     fmod <- lm(cvar~dscore+eth, data=tmp) ##"full" mod
#     nmod <- lm(cvar~1+eth, data=tmp) ##"null" mod
#     a_res <- anova(nmod, fmod, test = 'Chisq')
#     pval <- a_res$`Pr(>Chi)`[2]
#     
#     # # using ordinal logistic regression
#     # tmp$cvar <- as.factor(tmp$cvar)
#     # m <- polr(cvar ~ dscore, data = tmp, Hess=TRUE, method='probit')
#     #
#     # ## view a summary of the model
#     # ctable <- coef(summary(m))
#     # ## calculate and store p values
#     # p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
#     # pval <- p[1][[1]]
# 
#     all_pvals <- c(all_pvals,pval)
#     f_tested <- c(f_tested,f)
#     c_tested <- c(c_tested,colnames(clin_vars)[j])
#   }
# }

all_pvals <- c()
f_tested <- c()
c_tested <- c()
for (j in 1:ncol(clin_vars)) {
  print(j)

  # get donors in clin var that don't have an NA value
  d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]
  
  tmp <- as.data.frame(cbind(dsc[d_keep,], clin_vars[d_keep,j]))
  colnames(tmp)[1:ncol(dsc)] <- sapply(1:ncol(dsc),function(x){paste0('Factor',x)})
  colnames(tmp)[ncol(dsc)+1] <- 'cvar'

  # using ordinal logistic regression
  tmp$cvar <- as.factor(tmp$cvar)
  m <- polr(cvar ~ ., data = tmp, Hess=TRUE, method='probit')
    
  ## view a summary of the model
  ctable <- coef(summary(m))
  ## calculate and store p values
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  pval <- p[1:ncol(dsc)]
    
  all_pvals <- c(all_pvals,pval)
  f_tested <- c(f_tested,1:ncol(dsc))
  c_tested <- c(c_tested,rep(colnames(clin_vars)[j],ncol(dsc)))
}
all_pvals <- p.adjust(all_pvals,method='fdr')
all_pvals[order(all_pvals,decreasing=FALSE)][1:10]
f_tested[order(all_pvals,decreasing=FALSE)][1:10]
c_tested[order(all_pvals,decreasing=FALSE)][1:10]


# plotting f1 sle only sledaiscore and sliccscore associations
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"sliccscore"])]
tmp <- as.data.frame(cbind(dsc[d_keep,1], clin_vars[d_keep,"sliccscore"]))
colnames(tmp) <- c('dscore','cvar')

# # force cvar to be factor
# tmp$cvar <- as.factor(tmp$cvar)
tmp$cvar <- as.numeric(tmp$cvar)


lmres <- lm(cvar~dscore,data=tmp)
line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

pdf(file = "/home/jmitchel/figures/for_paper/sle_only_f1_slicc.pdf", useDingbats = FALSE,
    width = 4.5, height = 3.5)
ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point(alpha = 0.3,pch=19,size=3) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Factor 1 Donor Score') +
  ylab('SLICC Score') +
  theme_classic()
dev.off()

d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"sledaiscore"])]
tmp <- as.data.frame(cbind(dsc[d_keep,1], clin_vars[d_keep,"sledaiscore"]))
colnames(tmp) <- c('dscore','cvar')

# # force cvar to be factor
# tmp$cvar <- as.factor(tmp$cvar)
tmp$cvar <- as.numeric(tmp$cvar)

# adding reg line
lmres <- lm(cvar~dscore,data=tmp)
line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

pdf(file = "/home/jmitchel/figures/for_paper/sle_only_f1_sledai.pdf", useDingbats = FALSE,
    width = 4.5, height = 3.5)
ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point(alpha = 0.3,pch=19,size=3) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Factor 1 Donor Score') +
  ylab('SLEDAI Score') +
  theme_classic()
dev.off()


d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"sliccscore"])]
tmp <- as.data.frame(cbind(dsc[d_keep,7], clin_vars[d_keep,"sliccscore"]))
colnames(tmp) <- c('dscore','cvar')

# # force cvar to be factor
# tmp$cvar <- as.factor(tmp$cvar)

lmres <- lm(cvar~dscore,data=tmp)
line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

pdf(file = "/home/jmitchel/figures/for_paper/sliccscore_factor7.pdf", useDingbats = FALSE,
    width = 4.5, height = 3.5)
pdf(file = "/home/jmitchel/figures/for_paper/sliccscore_factor7.pdf", useDingbats = FALSE,
    width = 6, height = 3.5)
ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point(alpha = 0.3,pch=19,size=3) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Factor 7 Donor Score') +
  ylab('SLICC Score') +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))
dev.off()




d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"sledaiscore"])]
tmp <- as.data.frame(cbind(dsc[d_keep,3], clin_vars[d_keep,"sledaiscore"]))
colnames(tmp) <- c('dscore','cvar')

# force cvar to be factor
tmp$cvar <- as.factor(tmp$cvar)

ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point() +
  xlab('Factor 3 Donor Score') +
  ylab('SLEDAI score') +
  theme_classic()
dev.off()


# trying out the ordinal approach
require(MASS)
colnames(clin_vars)
var_test <- 6
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,var_test])]
tmp <- as.data.frame(cbind(dsc[d_keep,7], clin_vars[d_keep,var_test]))
colnames(tmp) <- c('dscore','cvar')

# force cvar to be factor
tmp$cvar <- as.factor(round(tmp$cvar))

m <- polr(cvar ~ dscore, data = tmp, Hess=TRUE, method='logistic')

## view a summary of the model
summary(m)
(ctable <- coef(summary(m)))
## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
(ctable <- cbind(ctable, "p value" = p))
ctable



# trying to use all variables in prediction
require(MASS)
colnames(clin_vars)
var_test <- 'crflymphadeno'
var_test <- 'crflupusneph'
# var_test <- 'sliccesrd'
# var_test <- 'sledaiantidsdna'
var_test <- 'acrantidsdna'
var_test <- 'sledaiscore'
var_test <- 'sliccscore'
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,var_test])]
tmp <- as.data.frame(cbind(dsc[d_keep,], clin_vars[d_keep,var_test]))
colnames(tmp)[ncol(tmp)] <- c('cvar')
colnames(tmp)[1:(ncol(tmp)-1)] <- sapply(1:(ncol(tmp)-1),function(x){
  paste0('factor',x)
})

# # force cvar to be factor
# tmp$cvar <- as.factor(round(tmp$cvar))
tmp$cvar <- as.factor(tmp$cvar)
# tmp$cvar <- as.numeric(tmp$cvar)

# lmres <- lm(cvar ~ ., data = tmp)
# summary(lmres)

# plot(tmp$cvar,tmp[,2])

# trying with logistic regression model
fmod <- glm(cvar~ ., data=tmp, family = "binomial") ##"full" mod
nmod <- glm(cvar~1, data=tmp, family = 'binomial') ##"null" mod
a_res <- anova(nmod, fmod, test = 'Chisq')
pval <- a_res$`Pr(>Chi)`[2]
summary(fmod)

# m <- polr(cvar ~ dscore1 + dscore2 + dscore1*dscore2, data = tmp, Hess=TRUE)
m <- polr(cvar ~ ., data = tmp, Hess=TRUE, method='probit')

## view a summary of the model
summary(m)
ctable <- coef(summary(m))
## calculate and store p values
p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
ctable <- cbind(ctable, "p value" = p)
ctable
ctable[1:10,'p value']
## run all the tests with the ordinal regression!!!!


ggplot(tmp,aes(x=cvar,y=factor2)) +
  geom_point()
ggplot(tmp,aes(x=cvar,y=factor2)) +
  geom_violin()


old_dsc <- dsc
# for shuffling testing
for (j in 1:ncol(dsc)) {
  dsc[,j] <- sample(dsc[,j])
}



# testing co-occurrance of LN and dsdna with factor 2
tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantidsdna']))
# tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantismith']))
tmp <- tmp[order(tmp[,1],decreasing=TRUE),]
colnames(tmp) <-  c('dscore','ln','dsdna')

# tmp_sub <- tmp[tmp$dscore > .1,]
# ln_count <- sum(tmp_sub$ln==1)
# ds_count <- sum(tmp_sub$dsdna==1)
# both_count <- sum(tmp_sub$ln==1 & tmp_sub$dsdna==1)
# # print(both_count/ln_count)
# print(both_count/ds_count)
#
# tmp_sub <- tmp[tmp$dscore > .05,]
# ln_count <- sum(tmp_sub$ln==1)
# ds_count <- sum(tmp_sub$dsdna==1)
# both_count <- sum(tmp_sub$ln==1 & tmp_sub$dsdna==1)
# # print(both_count/ln_count)
# print(both_count/ds_count)
#
# tmp_sub <- tmp
# ln_count <- sum(tmp_sub$ln==1)
# ds_count <- sum(tmp_sub$dsdna==1)
# both_count <- sum(tmp_sub$ln==1 & tmp_sub$dsdna==1)
# # print(both_count/ln_count)
# print(both_count/ds_count)

run_ds <- c()
run_both <- c()
dscores <- c()
ds_count <- 0
both_count <- 0
for (i in 1:nrow(tmp)) {
  if (tmp$dsdna[i]==1) {
    ds_count <- ds_count + 1
    if (tmp$ln[i]==1) {
      both_count <- both_count + 1
    }
    run_ds <- c(run_ds,ds_count)
    run_both <- c(run_both,both_count)
    dscores <- c(dscores,tmp$dscore[i])
  }
}
fracs <- c()
for (i in 1:length(run_ds)) {
  fracs <- c(fracs,run_both[i]/run_ds[i])
  # fracs <- c(fracs,log(run_both[i]/(run_ds[i]-run_both[i])))
  
}
plot(dscores,fracs)
AUC = trapz(rev(dscores),rev(fracs))

lmres <- lm(fracs~dscores)
summary(lmres)




## significance testing by shuffling and getting AUC values
require(pracma)
shuff_AUC <- c()
all_dsc <- list()
all_fracs <- list()
for (myiter in 1:10000) {
  # shuffle donor scores
  for (j in 1:ncol(dsc)) {
    dsc[,j] <- sample(dsc[,j])
  }

  # testing co-occurrance of LN and dsdna with factor 2
  tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantidsdna']))
  tmp <- tmp[order(tmp[,1],decreasing=TRUE),]
  colnames(tmp) <-  c('dscore','ln','dsdna')

  run_ds <- c()
  run_both <- c()
  dscores <- c()
  ds_count <- 0
  both_count <- 0
  for (i in 1:nrow(tmp)) {
    if (tmp$dsdna[i]==1) {
      ds_count <- ds_count + 1
      if (tmp$ln[i]==1) {
        both_count <- both_count + 1
      }
      run_ds <- c(run_ds,ds_count)
      run_both <- c(run_both,both_count)
      dscores <- c(dscores,tmp$dscore[i])
    }
  }
  fracs <- c()
  for (i in 1:length(run_ds)) {
    fracs <- c(fracs,run_both[i]/run_ds[i])
  }

  AUC = trapz(rev(dscores),rev(fracs))
  shuff_AUC <- c(shuff_AUC,AUC)
  all_dsc[[myiter]] <- dscores
  all_fracs[[myiter]] <- fracs
}
hist(shuff_AUC)
sum(shuff_AUC>.179)/10000
max_ndx <- order(shuff_AUC,decreasing=TRUE)[1]
plot(all_dsc[[max_ndx]],all_fracs[[max_ndx]])






# ## trying out a more 'enrichment' like version of this
# dsc <- dsc_old
# tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantidsdna']))
# # tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantismith']))
# tmp <- tmp[order(tmp[,1],decreasing=FALSE),]
# colnames(tmp) <-  c('dscore','ln','dsdna')

running_tab <- 0
myscore <- 0
running_track <- c(0)
dscores <- c()
for (i in 1:nrow(tmp)) {
  if (tmp$dsdna[i]==1) {
    if (tmp$ln[i]==1) {
      running_tab <- running_tab + 1
    } else {
      running_tab <- running_tab - 1
    }
    delta <- running_tab - running_track[length(running_track)]
    running_track <- c(running_track,running_tab)
    myscore <- myscore + delta*tmp$dscore[i]
    dscores <- c(dscores,tmp$dscore[i])
  }
}

plot(dscores,running_track)
mycor <- cor(dscores,running_track,method='spearman')
myAUC <- trapz(dscores,running_track)
myscore <- mycor * myAUC



all_scores <- c()
for (testndx in 1:1000) {

  # shuffle donor scores
  for (j in 1:ncol(dsc)) {
    dsc[,j] <- sample(dsc[,j])
  }

  # testing co-occurrance of LN and dsdna with factor 2
  tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantidsdna']))
  tmp <- tmp[order(tmp[,1],decreasing=FALSE),]
  colnames(tmp) <-  c('dscore','ln','dsdna')

  running_tab <- 0
  this_sc <- 0
  running_track <- c(0)
  dscores <- c()
  for (i in 1:nrow(tmp)) {
    if (tmp$dsdna[i]==1) {
      if (tmp$ln[i]==1) {
        running_tab <- running_tab + 1
      } else {
        running_tab <- running_tab - 1
      }
      delta <- running_tab - running_track[length(running_track)]
      running_track <- c(running_track,running_tab)
      this_sc <- myscore + delta*tmp$dscore[i]
      dscores <- c(dscores,tmp$dscore[i])
    }
  }
  all_scores <- c(all_scores,this_sc)

}

sum(all_scores>myscore)/1000






## trying a sliding window approach
# testing co-occurrance of LN and dsdna with factor 2
# old_dsc <- dsc
dsc <- old_dsc
tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantidsdna']))
# tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantismith']))
tmp <- tmp[order(tmp[,1],decreasing=TRUE),]
colnames(tmp) <-  c('dscore','ln','dsdna')

# remove rows that don't have dsdna==1
tmp <- tmp[tmp$dsdna==1,]
# window_size <- 21 #size I originally used
window_size <- 17
stored_counts <- c()
for (i in 1:(nrow(tmp)-window_size+1)) {
  tests <- tmp[i:(i+window_size-1),'ln']
  stored_counts <- c(stored_counts,sum(tests==1))
}
dscores <- tmp$dscore[(floor(window_size/2)+1):(nrow(tmp)-floor(window_size/2))]
plot(dscores,stored_counts)

lmres <- lm(stored_counts~dscores)
lmres <- summary(lmres)
myfstat <- lmres$fstatistic[[1]]
myfstat <- cor(stored_counts,dscores,method='spearman')
myfstat <- cor(stored_counts,dscores,method='pearson')
myfstat

plot_df <- cbind.data.frame(stored_counts,dscores)
pdf(file = "/home/jmitchel/figures/for_paper/LN_antidsDNA_link.pdf", useDingbats = FALSE,
    width = 4, height = 3.25)
ggplot(plot_df,aes(x=dscores,y=stored_counts)) +
  geom_point(alpha = 0.3,pch=19,size=2) +
  xlab('Factor 2 Donor Score (window center)') +
  ylab('Number of Joint Occurrences\nof LN + anti_dsDNA') +
  theme_classic() +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14))
dev.off()


all_scores <- c()
for (testndx in 1:1000) {
  
  # shuffle donor scores
  for (j in 1:ncol(dsc)) {
    dsc[,j] <- sample(dsc[,j])
  }
  
  # testing co-occurrance of LN and dsdna with factor 2
  # tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantidsdna']))
  # tmp <- as.data.frame(cbind(dsc[,4],clin_vars[,'crflupusneph'],clin_vars[,'acrantidsdna']))
  tmp <- as.data.frame(cbind(dsc[,2],clin_vars[,'crflupusneph'],clin_vars[,'acrantismith']))
  tmp <- tmp[order(tmp[,1],decreasing=TRUE),]
  colnames(tmp) <-  c('dscore','ln','dsdna')
  
  tmp <- tmp[tmp$dsdna==1,]
  # tmp <- tmp[tmp$dsdna==0,]
  
  stored_counts <- c()
  for (i in 1:(nrow(tmp)-window_size+1)) {
    tests <- tmp[i:(i+window_size-1),'ln']
    stored_counts <- c(stored_counts,sum(tests==1))
  }
  dscores <- tmp$dscore[(floor(window_size/2)+1):(nrow(tmp)-floor(window_size/2))]
  lmres <- lm(stored_counts~dscores)
  lmres <- summary(lmres)
  # fstat <- lmres$fstatistic[[1]]
  fstat <- cor(stored_counts,dscores,method='spearman')
  all_scores <- c(all_scores,fstat)
}

sum(all_scores>myfstat)/1000
sum(all_scores<myfstat)/1000

sum(all_scores>myfstat)/10000
sum(all_scores<myfstat)/10000




## plotting LR curves for some predictors
dsc <- old_dsc
var_test <- 'crflymphadeno'
var_test <- 'crflupusneph'
# var_test <- 'sliccesrd'
# var_test <- 'sledaiantidsdna'
var_test <- 'acrantidsdna'
var_test <- 'sledaiscore'
var_test <- 'sliccscore'
d_keep <- rownames(clin_vars)[!is.na(clin_vars[,var_test])]
tmp <- as.data.frame(cbind(dsc[d_keep,], clin_vars[d_keep,var_test]))
colnames(tmp)[ncol(tmp)] <- c('cvar')
colnames(tmp)[1:(ncol(tmp)-1)] <- sapply(1:(ncol(tmp)-1),function(x){
  paste0('factor',x)
})

tmp$cvar <- as.factor(tmp$cvar)
# ggplot(tmp, aes(x=factor2, y=cvar)) + geom_point() + 
#   stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)

ggplot(tmp,aes(x=cvar,y=factor2)) +
  geom_violin(scale='count')
ggplot(tmp,aes(x=cvar,y=factor2)) +
  geom_point()


# trying with logistic regression model
fmod <- glm(cvar~factor2, data=tmp, family = "binomial") ##"full" mod
nmod <- glm(cvar~1, data=tmp, family = 'binomial') ##"null" mod
a_res <- anova(nmod, fmod, test = 'Chisq')
pval <- a_res$`Pr(>Chi)`[2]
summary(fmod)




## getting meta data statistics information
meta <- pbmc_container$scMinimal_ctype[[1]]$metadata
meta <- unique(meta[,c('donors','Status','pool','Age')])
rownames(meta) <- meta$donors
meta$donors <- NULL
table(meta$pool)
table(meta$Status)
table(meta$Age)
dim(pbmc_container$scMinimal_full$count_data)
##


## compute full model AUC for predicting disease status
library(pROC)
dsc <- pbmc_container$tucker_results[[1]]
tmp <- cbind.data.frame(dsc, as.character(meta[rownames(dsc),'Status']))
colnames(tmp)[ncol(tmp)] <- c('cvar')
colnames(tmp)[1:(ncol(tmp)-1)] <- sapply(1:(ncol(tmp)-1),function(x){
  paste0('factor',x)
})

tmp$cvar <- as.factor(tmp$cvar)

sle_mod <- glm(cvar~factor1 + factor2 + factor3 + factor4 + factor5 +
                 factor6 + factor7 + factor8 + factor10, data=tmp, family = "binomial")
train_pred <- predict(sle_mod, tmp[,1:10], type="response")

pdf(file = "/home/jmitchel/figures/for_paper/status_auc_plot.pdf", useDingbats = FALSE,
    width = 4, height = 4)
pROC_obj <- roc(tmp$cvar,train_pred,
                smoothed = TRUE,
                ci=FALSE,
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=FALSE,
                print.auc=TRUE, show.thres=TRUE)
dev.off()

##




### testing factors against meds
 
# load data of categorical variables
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_meds_cleaned.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars[,'Sample ID']
clin_vars[,'Sample ID'] <- NULL

# make all NA into zeros, since no 0 are put in the table
clin_vars[is.na(clin_vars)] <- 0

# separate out pred dose as it's the only continuous variable here
pred_dose <- clin_vars[,'pred_dose',drop=FALSE]
clin_vars[,'pred_dose'] <- NULL

# make sure there are no columns of all zeros
colSums(clin_vars)

# need to remove a few columns that have only 1 or 0 donors on the med
clin_vars[,c('solumedrol','rx_abatacept','rx_cyclophosphamide','rx_etanercept',
             'rx_IGG','rx_leflunomide','rx_rituximab','rx_sulfasalazine')] <- NULL


# get tucker donor scores to test
dsc <- pbmc_container$tucker_results[[1]]

## get donors in both dsc and in clin_vars
# trim donor IDs in dsc
trim_names <- sapply(rownames(dsc), function(x) {
  strsplit(x,split='_')[[1]][[1]]
})
names(trim_names) <- c()
old_names <- rownames(dsc)
names(old_names) <- trim_names
rownames(dsc) <- trim_names

# get donors in both dataframes
d_both <- rownames(clin_vars)[rownames(clin_vars) %in% rownames(dsc)]

# limit both dataframes to just the intersection of donors and in the same order
dsc <- dsc[d_both,]
clin_vars <- clin_vars[d_both,]



all_pvals <- c()
f_tested <- c()
c_tested <- c()
# loop through the variables to test
for (j in 1:ncol(clin_vars)) {
  print(j)
  # loop through factors
  for (f in 1:ncol(dsc)) {
    # get donors in clin var that don't have an NA value
    d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]
    
    tmp <- as.data.frame(cbind(dsc[d_keep,f], clin_vars[d_keep,j]))
    colnames(tmp) <- c('dscore','cvar')
    
    # force cvar to be factor
    tmp$cvar <- as.factor(tmp$cvar)
    
    # # compute linear model
    # lmres <- summary(lm(dscore ~ cvar, data=tmp))
    # pval <- stats::pf(lmres$fstatistic[1],lmres$fstatistic[2],lmres$fstatistic[3],lower.tail=FALSE)
    
    # # trying with two sample t-test
    # t_res <- try(stats::t.test(dscore ~ cvar, data=tmp,
    #                        alternative = "two.sided",var.equal = FALSE))
    # if (class(t_res) == 'try-error') {
    #   next
    # }
    # pval <- t_res$p.value
    
    
    # trying with logistic regression model
    fmod <- glm(cvar~dscore, data=tmp, family = "binomial") ##"full" mod
    nmod <- glm(cvar~1, data=tmp, family = 'binomial') ##"null" mod
    a_res <- anova(nmod, fmod, test = 'Chisq')
    pval <- a_res$`Pr(>Chi)`[2]
    
    all_pvals <- c(all_pvals,pval)
    f_tested <- c(f_tested,f)
    c_tested <- c(c_tested,colnames(clin_vars)[j])
    
  }
}
all_pvals <- p.adjust(all_pvals,method='fdr')
all_pvals[order(all_pvals,decreasing=FALSE)]
f_tested[order(all_pvals,decreasing=FALSE)]
c_tested[order(all_pvals,decreasing=FALSE)]




tmp <- cbind.data.frame(dsc[,5],clin_vars[,'prednisone'])
colnames(tmp) <- c('dscore','cvar')
head(tmp)
ggplot(tmp,aes(x=cvar,y=dscore)) +
  geom_point()

# testing factor 5 against prednisone level
tmp <- cbind.data.frame(dsc[,5],pred_dose[d_both,1])
colnames(tmp) <- c('dscore','cvar')
# head(tmp)
ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point()
# looks okay but might need to adjust for weight...

tmp <- cbind.data.frame(dsc[,5],pred_dose[d_both,1])
colnames(tmp) <- c('dscore','cvar')
head(tmp)
ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point()


# get column of totalweight or bmi
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars$subjectid
clin_vars$subjectid <- NULL

tmp <- cbind.data.frame(dsc[,5],pred_dose[d_both,1]/clin_vars[d_both,'totalweight'])
colnames(tmp) <- c('dscore','cvar')
head(tmp)
ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point()


## remove outlier for lm
tmp <- cbind.data.frame(dsc[,5],pred_dose[d_both,1])
colnames(tmp) <- c('dscore','cvar')
row_max <- which(tmp$cvar==max(tmp$cvar))
tmp2 <- tmp[-row_max,]
lmres <- lm(cvar~dscore,data=tmp2)
summary(lmres)

line_range <- seq(min(tmp$dscore),max(tmp$dscore),.001)
line_dat <- c(line_range*lmres$coefficients[[2]] + lmres$coefficients[[1]])
line_df <- cbind.data.frame(line_range,line_dat)
colnames(line_df) <- c('myx','myy')

# make plot with best fit line
pdf(file = "/home/jmitchel/figures/for_paper/prednisone_dose_reg.pdf", useDingbats = FALSE,
    width = 6, height = 4)
ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point(alpha = 0.3,pch=19,size=3) +
  geom_line(data=line_df,aes(x=myx,y=myy)) +
  xlab('Factor 5 Donor Score') +
  ylab('Prednisone Dose') +
  theme_classic() +
  theme(axis.text=element_text(size=10),
      axis.title=element_text(size=14))
dev.off()

# make plot for binary classification
# need to reload med clinical data
clin_vars <- read_excel('/home/jmitchel/data/lupus_data/SLE_meds_cleaned.xlsx')
clin_vars <- as.data.frame(clin_vars)
rownames(clin_vars) <- clin_vars[,'Sample ID']
clin_vars[,'Sample ID'] <- NULL
clin_vars[is.na(clin_vars)] <- 0
clin_vars <- clin_vars[d_both,]

tmp <- cbind.data.frame(dsc[,5],clin_vars[,'prednisone'])
colnames(tmp) <- c('dscore','cvar')
tmp$cvar_word <- sapply(tmp$cvar,function(x){
  if (x==1) {
    return('On Prednisone')
  } else {
    return('Not on Prednisone')
  }
})
tmp$cvar_word <- as.factor(tmp$cvar_word)

pdf(file = "/home/jmitchel/figures/for_paper/prednisone_LR.pdf", useDingbats = FALSE,
    width = 4.5, height = 3.5)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 5 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
dev.off()






#crferosivearth sliccdefarthritis
tmp <- cbind.data.frame(dsc[,3],clin_vars[,'sliccdefarthritis'])
colnames(tmp) <- c('dscore','cvar')
tmp$cvar_word <- sapply(tmp$cvar,function(x){
  if (x==1) {
    return('deforming arthritis')
  } else {
    return('no deforming arthritis')
  }
})
tmp$cvar_word <- as.factor(tmp$cvar_word)
pdf(file = "/home/jmitchel/figures/for_paper/SLE_only_arthritis.pdf", useDingbats = FALSE,
    width = 4.5, height = 3)
ggplot(tmp,aes(x=cvar_word,y=dscore)) +
  geom_violin() +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.75, binwidth = .01) +
  ylab('Factor 3 Donor Score') +
  xlab('') +
  coord_flip() +
  theme_bw()
dev.off()

tmp3 <- tmp[tmp$cvar==1,]








d_keep <- rownames(clin_vars)[!is.na(clin_vars[,"sliccscore"])]
tmp <- as.data.frame(cbind(dsc[d_keep,3], clin_vars[d_keep,"sliccscore"]))
colnames(tmp) <- c('dscore','cvar')

# force cvar to be factor
tmp$cvar <- as.factor(tmp$cvar)
pdf(file = "/home/jmitchel/figures/for_paper/SLE_only_f3_sliccscore.pdf", useDingbats = FALSE,
    width = 4.5, height = 3)
ggplot(tmp,aes(x=dscore,y=cvar)) +
  geom_point() +
  xlab('Factor 3 Donor Score') +
  ylab('SLICC Score') +
  theme_classic()
dev.off()


tmp2 <- tmp[tmp$dscore>.05,]






# interaction tests between ifn factor and other factors in determining symptoms
all_pvals <- c()
f_tested <- c()
c_tested <- c()
# loop through the variables to test
for (j in 1:ncol(clin_vars)) {
  print(j)
  # loop through factors
  # for (f in 1:ncol(dsc)) {
  for (f in 2:4) {
    # get donors in clin var that don't have an NA value
    d_keep <- rownames(clin_vars)[!is.na(clin_vars[,j])]
    
    tmp <- as.data.frame(cbind(dsc[d_keep,1],dsc[d_keep,f], clin_vars[d_keep,j]))
    colnames(tmp) <- c('f1_dsc','f2_dsc','cvar')
    
    # # force cvar to be factor
    # tmp$cvar <- as.factor(tmp$cvar)
    
    # # if smallest level has less thatn n donors skip this one
    # if (min(table(tmp$cvar)) < 10) {
    #   next
    # }
    
    # # trying with logistic regression model
    # fmod <- glm(cvar ~ f1_dsc + f2_dsc + f1_dsc*f2_dsc, data=tmp, family = "binomial") ##"full" mod
    # nmod <- glm(cvar ~ f1_dsc + f2_dsc, data=tmp, family = 'binomial') ##"null" mod
    # a_res <- anova(nmod, fmod, test = 'Chisq')
    # pval <- a_res$`Pr(>Chi)`[2]
    
    fmod <- lm(cvar ~ f1_dsc + f2_dsc + f1_dsc*f2_dsc, data=tmp)
    nmod <- lm(cvar ~ f1_dsc + f2_dsc, data=tmp)
    a_res <- anova(nmod, fmod, test = 'Chisq')
    pval <- a_res$`Pr(>Chi)`[2]
    
    all_pvals <- c(all_pvals,pval)
    f_tested <- c(f_tested,f)
    c_tested <- c(c_tested,colnames(clin_vars)[j])
    
  }
}
all_pvals <- p.adjust(all_pvals,method='fdr')
all_pvals[order(all_pvals,decreasing=FALSE)][1:40]
f_tested[order(all_pvals,decreasing=FALSE)][1:40]
c_tested[order(all_pvals,decreasing=FALSE)][1:40]








## trying to limit by ethnicity
cv <- read_excel('/home/jmitchel/data/lupus_data/SLE_clinical_data.xlsx')
cv <- as.data.frame(cv)
rownames(cv) <- cv$subjectid
rasian <- cv[,'raceasian',drop=F]
head(rasian)
rasian <- rasian[d_both,,drop=FALSE]
d_both <- d_both[rasian==1]
d_both <- d_both[!is.na(d_both)]





## testing LN dsdna w linear model
head(tmp)
tmp$ln <- as.factor(tmp$ln)
t.test(dscore~ln,data=tmp)

# with logistic regression
fmod <- glm(ln~dscore, data=tmp, family = "binomial") ##"full" mod
nmod <- glm(ln~1, data=tmp, family = 'binomial') ##"null" mod
fmod <- glm(dsdna~dscore, data=tmp, family = "binomial") ##"full" mod
nmod <- glm(dsdna~1, data=tmp, family = 'binomial') ##"null" mod
a_res <- anova(nmod, fmod, test = 'Chisq')
pval <- a_res$`Pr(>Chi)`[2]
pval

tmp$ln <- as.factor(tmp$ln)
ggplot(tmp,aes(x=ln,y=dscore)) +
  geom_violin(scale='count')
ggplot(tmp,aes(x=ln,y=dscore)) +
  geom_violin(scale='area')

tmp <- tmp[tmp$ln==1,]
tmp$dsdna <- as.factor(tmp$dsdna)
ggplot(tmp,aes(x=dsdna,y=dscore)) +
  geom_violin(scale='count')

## trying with enrichment
mypaths=list()
mypaths[['ln']] <- rownames(tmp)[tmp[,'ln']==1]

myranks <- tmp[,'dscore']
names(myranks) <- rownames(tmp)

fgseaRes <- fgsea::fgseaSimple(pathways = mypaths,
                         stats    = myranks,
                         minSize  = 0,
                         nperm=10000,
                         maxSize  = 5000,
                         gseaParam=2)

print(fgseaRes)



## trying interaction model
tmp$ln <- as.factor(tmp$ln)
tmp$dsdna <- as.factor(tmp$dsdna)

nmod <- lm(dscore~ln+dsdna,data=tmp)
fmod <- lm(dscore~ln+dsdna+ln*dsdna,data=tmp)
a_res <- anova(nmod, fmod, test = 'Chisq')
pval <- a_res$`Pr(>Chi)`[2]
pval

nmod <- glm(dsdna~ln, data=tmp, family = 'binomial') ##"null" mod
fmod <- glm(dsdna~ln+dscore, data=tmp, family = "binomial") ##"full" mod
a_res <- anova(nmod, fmod, test = 'Chisq')
pval <- a_res$`Pr(>Chi)`[2]
pval



