percentile_threshold <- function(A, gamma) {
  # Thresholding function F in the objective function
  # 0.5 ||F(L - Y)||^2_2
  #
  # TO DO:
  #   * Implement in C++
  #   * Error handling

  row_pctl <- apply(abs(A), 1, quantile, gamma)
  col_pctl <- apply(abs(A), 2, quantile, gamma)

  # to do: think of a faster way to do this...
  A_out <- A
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A)) {
      if ((abs(A[i,j]) <= row_pctl[i]) & (abs(A[i,j]) <= col_pctl[j])) {
        A_out[i,j] <- 0
      }
    }
  }
  A_out
}
