my_quantile <- function(x, prob) {
  n <- length(x)
  index <- 1 + (n - 1) * prob
  lo <- floor(index)
  hi <- ceiling(index)
  x <- sort(x, partial = unique(c(lo, hi)))
  h <- index - lo

  (1 - h) * x[lo] + h * x[hi]
}
threshold1 <- function(X, gamma, sparsity) {
  # Thresholding
  #
  # TO DO:
  #   * Documentation
  #   * C++ implementation
  n1 <- nrow(X); n2 <- ncol(X)

  t1 <- rep(1, n1); t2 <- rep(1, n2)
  X <- X * sparsity

  for (i in 1:n1) {
    tt <- sort(abs(X[i, ]), decreasing = TRUE)
    t1[i] <- tt[floor(gamma * sum(sparsity[i, ])) + 1]
  }
  for (j in 1:n2) {
    tt <- sort(abs(X[, j]), decreasing = TRUE)
    t2[j] <- tt[floor(gamma * sum(sparsity[, j])) + 1]
  }

  threshold1 <- abs(X) <= matrix(rep(t1, each = n2), ncol = n2, byrow = TRUE)
  threshold2 <- abs(X) <= matrix(rep(t2, each = n1), nrow = n1)

  X_thresholded <- X * (as.double(threshold1 + threshold2) >= 1)

  X_thresholded
}
threshold2 <- function(X, gamma) {
  row_pctls <- apply(abs(X), 1, quantile, probs = 1 - gamma)
  col_pctls <- apply(abs(X), 2, quantile, probs = 1 - gamma)

  percentile_threshold_cpp(X, row_pctls, col_pctls)
}
threshold3 <- function(X, gamma) {
  row_pctls <- apply(abs(X), 1, my_quantile, prob = 1 - gamma)
  col_pctls <- apply(abs(X), 2, my_quantile, prob = 1 - gamma)

  percentile_threshold_cpp(X, row_pctls, col_pctls)
}
threshold4 <- function(X, gamma) {
  row_pctls <- apply(abs(X), 1, percentile_cpp, prob = 1 - gamma)
  col_pctls <- apply(abs(X), 2, percentile_cpp, prob = 1 - gamma)

  percentile_threshold_cpp(X, row_pctls, col_pctls)
}





