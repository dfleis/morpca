library(morpca)
gamma <- 0.2
n1 <- 3
n2 <- 2

A <- matrix(rnorm(n1 * n2), nrow = n1)

row_pctls <- apply(abs(A), 1, quantile, 1 - gamma)
col_pctls <- apply(abs(A), 2, quantile, 1 - gamma)

B <- percentile_threshold(A, gamma)
sum(B == 0)/(n1 * n2)
