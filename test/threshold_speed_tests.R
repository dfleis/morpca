library(Rcpp)
sourceCpp("src/helpers.cpp")
library(morpca)
library(microbenchmark)

#==========================#
#===== set parameters =====#
#==========================#
n1 <- 500
n2 <- 600
gamma <- 0.2
sparsity_mat <- matrix(1, nrow = n1, ncol = n2)

#========================#
#===== generate data ====#
#========================#
X <- round(matrix(rnorm(n1 * n2), nrow = n1), 3)

Xthresh1 <- threshold1(X, gamma, sparsity_mat)
Xthresh2 <- threshold2(X, gamma)
Xthresh3 <- threshold3(X, gamma)
Xthresh4 <- threshold4(X, gamma)

norm(Xthresh1 - Xthresh2, "f")
norm(Xthresh1 - Xthresh3, "f")
norm(Xthresh1 - Xthresh4, "f")
norm(Xthresh3 - Xthresh4, "f")

#======================#
#===== speed tests ====#
#======================#
pt <- proc.time()
mb <- microbenchmark(threshold1(X, gamma, sparsity_mat),
                     threshold2(X, gamma),
                     threshold3(X, gamma),
                     threshold4(X, gamma), times = 10, unit = "s")
proc.time() - pt
summary(mb)
boxplot(mb, outline = F)

pt <- proc.time()
row_pctls <- apply(abs(X), 1, quantile, probs = 1 - gamma)
col_pctls <- apply(abs(X), 2, quantile, probs = 1 - gamma)
mb <- microbenchmark(percentile_threshold_cpp(X, row_pctls, col_pctls),
                     times = 100, unit = "s")
proc.time() - pt
summary(mb)


x <- rnorm(5)
probs <- 0.2
n <- length(x)


qs == quantile(x, probs = 0.2)

quantile(x, probs = 0.2)












