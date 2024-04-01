

n_samples <- 10000
r_vals1 <- sample(10:100, n_samples, replace=TRUE)
# r_vals1 <- sample(5:10, n_samples, replace=TRUE)
r_vals2 <- r_vals1**(3) + rnorm(n_samples,sd = 1)
counts1 <- rpois(n_samples,lambda=r_vals1)
# counts2 <- rpois(n_samples,lambda=r_vals2)
counts2 <- rpois(n_samples,lambda=r_vals2)

plot(counts1,counts2)
plot(log(counts1),log(counts2))










