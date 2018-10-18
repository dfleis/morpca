#' Hard-Thresholding
#'
#' Computes the hard-thresholding of a given real-valued matrix \eqn{X}
#' according to thresholding percentile \eqn{\gamma}.
#'
#' @details The hard-thresholding procedure will set the \eqn{(i,j)}-th
#' entry of X to 0 if
#' \deqn{|X_{i,j}| > |X^{[\gamma]}_{i,.}|}
#' and
#' \deqn{|X_{i,j}| > |X^{[\gamma]}_{.,j}|},
#' where \eqn{X_{i,.}} denotes the \eqn{i}-th row of \eqn{X}, and
#' \eqn{X_{.,j}} represents the \eqn{j}-th column of \eqn{X}. Additionally,
#' \eqn{X^{[\gamma]}_{i,.}} and \eqn{X^{[\gamma]}_{.,j}} represent the
#' \eqn{\gamma}-th percentile of the absolute values of the entries of the
#' \eqn{i}-th row and the \eqn{j}-th column, respectively.
#'
#' That is, the elements of \eqn{X} which simultaneously exceed the
#' \eqn{\gamma}-th most extreme fraction of entries (in terms of absolute
#' value) in the corresponding rows and columns of \eqn{X} will be removed
#' and replaced with 0.
#'
#' @param X Input real-valued data matrix.
#' @param gamma Thresholding percentile.
#' @param sparsity TO DO...
#'
#' @return Returns a real-valued matrix with thresholded entries.
#'
#' @export
threshold <- function(X, gamma, sparsity) {
  # TO DO:
  #   * Sparsity
  #   * Possibly implement the apply() functionality within C++?
  # n1 <- nrow(X); n2 <- ncol(X)
  #
  # t1 <- rep(1, n1); t2 <- rep(1, n2)
  # X <- X * sparsity
  #
  # for (i in 1:n1) {
  #   tt <- sort(abs(X[i, ]), decreasing = TRUE)
  #   t1[i] <- tt[floor(gamma * sum(sparsity[i, ])) + 1]
  # }
  # for (j in 1:n2) {
  #   tt <- sort(abs(X[, j]), decreasing = TRUE)
  #   t2[j] <- tt[floor(gamma * sum(sparsity[, j])) + 1]
  # }
  #
  # threshold1 <- abs(X) <= matrix(rep(t1, each = n2), ncol = n2, byrow = TRUE)
  # threshold2 <- abs(X) <= matrix(rep(t2, each = n1), nrow = n1)
  #
  # X_thresholded <- X * (as.double(threshold1 + threshold2) >= 1)
  #
  # X_thresholded
  X_abs <- abs(X)
  row_pctls <- row_pctls_cpp(X_abs, 1 - gamma)
  col_pctls <- col_pctls_cpp(X_abs, 1 - gamma)

  percentile_threshold_cpp(X, X_abs, row_pctls, col_pctls)
}

