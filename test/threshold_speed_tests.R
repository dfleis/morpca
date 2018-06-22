library(Rcpp)
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
X <- matrix(rnorm(n1 * n2), nrow = n1)

#======================#
#===== speed tests ====#
#======================#
pt <- proc.time()
Y <- threshold(X, gamma, sparsity_mat)
proc.time() - pt










