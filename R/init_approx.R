init_approx <- function(Y, r) {
  # computes the best rank-r approximation of percentile_threshold(Y)
  # via the Eckhard-Young Theorem

  # To do:
  #   * Check if `r` is an integer.

  svd_out <- svd(Y)

  U <- svd_out$u[,1:r]
  SIGMA <- diag(svd_out$d[1:r])
  V <- svd_out$v[,1:r]

  U %*% tcrossprod(SIGMA, V)
}
