#
# Compares the time taken for the entire orthographic and projective
# retraction functions vs. only the thresholding function
#

library(Rcpp)
sourceCpp("src/helpers.cpp")
library(morpca)
library(microbenchmark)

set.seed(124)
n1 <- 500
n2 <- 600
r <- 5
gamma <- 0.2
eta <- 0.7

L <- matrix(rnorm(n1 * n2), nrow = n1)
Y <- matrix(rnorm(n1 * n2), nrow = n1)
sparsity <- matrix(1, nrow = n1, ncol = n2)

pt <- proc.time()
mb <- microbenchmark(orthographic_retraction(L, Y, r, gamma, eta, sparsity),
                     threshold(L - Y, gamma, sparsity), times = 100, unit = "s")
proc.time() - pt
summary(mb)
boxplot(mb, outline = F)



