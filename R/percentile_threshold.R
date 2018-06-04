percentile_threshold <- function(A, gamma) {
  row_pctls <- apply(abs(A), 1, quantile, gamma)
  col_pctls <- apply(abs(A), 2, quantile, gamma)

  percentile_threshold_cpp(A, row_pctls, col_pctls)
}
