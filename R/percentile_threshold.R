percentile_threshold <- function(A, gamma) {
  # Thresholding function F in our objective function
  # 0.5 ||F(L - Y)||^2_2 and gradient F(L - Y)
  #
  # TO DO:
  #   * Error handling

  row_pctls <- apply(abs(A), 1, quantile, gamma)
  col_pctls <- apply(abs(A), 2, quantile, gamma)

  percentile_threshold_cpp(A, row_pctls, col_pctls)
}
