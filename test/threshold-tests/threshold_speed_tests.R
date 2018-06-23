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
Xthresh5 <- threshold5(X, gamma)

norm(Xthresh1 - Xthresh2, "f")
norm(Xthresh1 - Xthresh3, "f")
norm(Xthresh1 - Xthresh4, "f")
norm(Xthresh1 - Xthresh5, "f")

#======================#
#===== speed tests ====#
#======================#
pt <- proc.time()
mb <- microbenchmark(threshold1(X, gamma, sparsity_mat),
                     threshold2(X, gamma),
                     threshold3(X, gamma),
                     threshold4(X, gamma),
                     threshold5(X, gamma), times = 10, unit = "s")
proc.time() - pt
summary(mb)
boxplot(mb, outline = F)


pt <- proc.time()
mb <- microbenchmark(apply(X, 2, percentile_cpp, prob = gamma),
                     col_pctls_cpp(X, gamma), times = 100, unit = "ms")
proc.time() - pt
summary(mb)
boxplot(mb, outline = F)

pt <- proc.time()
mb <- microbenchmark(apply(X, 1, percentile_cpp, prob = gamma),
                     row_pctls_cpp(X, gamma), times = 100, unit = "ms")
proc.time() - pt
summary(mb)
boxplot(mb, outline = F)






