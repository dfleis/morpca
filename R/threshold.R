#' Percentile thresholding of matrix entries
#'
#' Hard-thresholds the entries of a given real-valued matrix \eqn{X}
#' according to some thresholding percentile \eqn{\gamma} for which
#' the \eqn{(i, j)}-th entry of of \eqn{X} is set to zero if it
#' simultaneously exceeds the \eqn{(1-\gamma)}-th percentile of the
#' absolute values in row \eqn{i} and column \eqn{j}. That is,
#' \eqn{X_{i,j}} is set to zero if it is among the most extreme
#' \eqn{\gamma\times 100}% of values within its row and column.
#'
#' @details The hard-thresholding procedure will set the \eqn{(i,j)}-th
#' entry of X to 0 if
#' \deqn{|X_{i,j}| > |X^{[\gamma]}_{i,.}|}
#' and
#' \deqn{|X_{i,j}| > |X^{[\gamma]}_{.,j}|,}
#' where \eqn{X_{i,.}} denotes the \eqn{i}-th row of \eqn{X}, and
#' \eqn{X_{.,j}} represents the \eqn{j}-th column of \eqn{X}. Additionally,
#' \eqn{X^{[\gamma]}_{i,.}} and \eqn{X^{[\gamma]}_{.,j}} represent the
#' \eqn{(1-\gamma)}-th percentile of the absolute values of the entries of the
#' \eqn{i}-th row and the \eqn{j}-th column, respectively.
#'
#' @param X Input real-valued matrix.
#' @param gamma Thresholding percentile (in \eqn{[0, 1]}).
#'
#' @return Returns a real-valued matrix with \eqn{\gamma}-thresholded
#' entries.
#'
#' @export
threshold <- function(X, gamma) {
  # TO DO:
  #   * Sparsity?
  #   * Make more gooder.
  n1 <- nrow(X); n2 <- ncol(X)

  t1 <- rep(1, n1); t2 <- rep(1, n2)
  X <- X #* sparsity

  for (i in 1:n1) {
    tt <- sort(abs(X[i, ]), decreasing = TRUE)
    t1[i] <- tt[floor(gamma * n2) + 1] # why use this method instead of quantile()? speed?
    #t1[i] <- tt[floor(gamma * sum(sparsity[i, ])) + 1]
  }
  for (j in 1:n2) {
    tt <- sort(abs(X[, j]), decreasing = TRUE)
    t2[j] <- tt[floor(gamma * n1) + 1]
    #t2[j] <- tt[floor(gamma * sum(sparsity[, j])) + 1]
  }

  threshold1 <- abs(X) <= matrix(rep(t1, each = n2), ncol = n2, byrow = TRUE)
  threshold2 <- abs(X) <= matrix(rep(t2, each = n1), nrow = n1)

  X_thresholded <- X * (as.double(threshold1 + threshold2) >= 1)

  X_thresholded
}
