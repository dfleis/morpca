library(microbenchmark)
library(Rcpp)
source("test/yi/code/threshold.R")
source("R/percentile_threshold.R")
sourceCpp("src/helpers.cpp")


n1 <- 50
n2 <- 60
gamma <- alpha <- 0.7

X <- matrix(rnorm(n1 * n2), nrow = n1)
sparsity <- matrix(rep(1, nrow(X) * ncol(X)), nrow = nrow(X))

T1 <- threshold(X, alpha, sparsity)
T2 <- percentile_threshold(X, gamma)

pt <- proc.time()
mb <- microbenchmark(threshold(X, alpha, sparsity), percentile_threshold(X, gamma), times = 100, unit = "ms")
tm <- proc.time() - pt

boxplot(mb, outline = F)
summary(mb)
tm

sum(abs(T1 - T2))
