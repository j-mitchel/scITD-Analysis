
library(ggplot2)

# say we have two genes g1 and g2
g1 <- c(1:10)
g2 <- g1*2

dat <- as.data.frame(cbind(g1,g2))

# here is the data plotted as centered, not scaled
ggplot(dat,aes(x=g1,y=g2)) +
  geom_point()

# now select points (samples) to be in each batch
b1 <- sample(1:10,5)
b2 <- c(1:10)[!(c(1:10) %in% b1)]

# first try multiplicative batch effect.
dat_mult <- dat
dat_mult[b1,1] <- dat_mult[b1,1] * 1.5
dat_mult[b1,2] <- dat_mult[b1,2] * .9

dat_mult[b2,1] <- dat_mult[b2,1] * .81
dat_mult[b2,2] <- dat_mult[b2,2] * 1.15

# plot data as is
ggplot(dat_mult,aes(x=g1,y=g2)) +
  geom_point()

# run pca
mult_pca <- prcomp(dat_mult,scale.=FALSE,center=TRUE)
mult_pc1 <- mult_pca[["rotation"]][,1]

# center data
dat_mult_cent <- as.data.frame(scale(dat_mult,scale=FALSE))

# plot centered, showing pc1
line_x_min <- min(dat_mult_cent[,1])
line_x_max <- max(dat_mult_cent[,1])

ggplot(dat_mult_cent,aes(x=g1,y=g2)) +
  geom_point() +
  geom_segment(aes(x=0,y=0,
                   xend=line_x_max,
                   yend=(line_x_max*mult_pc1[2]/mult_pc1[1]))) +
  geom_segment(aes(x=0,y=0,
                   xend=line_x_min,
                   yend=(line_x_min*mult_pc1[2]/mult_pc1[1])))
  
# now do it with an additive batch effect
dat_add <- dat
dat_add[b1,1] <- dat_add[b1,1] + .6
dat_add[b1,2] <- dat_add[b1,2] - .5

dat_add[b2,1] <- dat_add[b2,1] - .45
dat_add[b2,2] <- dat_add[b2,2] + .65

# plot data as is
ggplot(dat_add,aes(x=g1,y=g2)) +
  geom_point()

# run pca
add_pca <- prcomp(dat_add,scale.=FALSE,center=TRUE)
add_pc1 <- add_pca[["rotation"]][,1]

# center data
dat_add_cent <- as.data.frame(scale(dat_add,scale=FALSE))

# plot centered, showing pc1
line_x_min <- min(dat_add_cent[,1])
line_x_max <- max(dat_add_cent[,1])

ggplot(dat_add_cent,aes(x=g1,y=g2)) +
  geom_point() +
  geom_segment(aes(x=0,y=0,
                   xend=line_x_max,
                   yend=(line_x_max*add_pc1[2]/add_pc1[1]))) +
  geom_segment(aes(x=0,y=0,
                   xend=line_x_min,
                   yend=(line_x_min*add_pc1[2]/add_pc1[1])))

ggplot(dat_add_cent,aes(x=g1,y=g2)) +
  geom_point() +
  geom_segment(aes(x=line_x_min,
                   y=line_x_min * add_pc1[2]/add_pc1[1],
                   xend=line_x_max,
                   yend=line_x_max*add_pc1[2]/add_pc1[1]))






