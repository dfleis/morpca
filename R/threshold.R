threshold <- function(X, gamma, sparsity) {
  # Thresholding
  #
  # TO DO:
  #   * Documentation
  #   * Sparsity
  #   * Possibly implement the apply() functionality within C++?
  X_abs <- abs(X)
  row_pctls <- apply(X_abs, 1, percentile_cpp, prob = 1 - gamma)
  col_pctls <- apply(X_abs, 2, percentile_cpp, prob = 1 - gamma)

  percentile_threshold_cpp(X, X_abs, row_pctls, col_pctls)
}

