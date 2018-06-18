orthographic_retraction <- function(L, Y, eta, gamma, sparsity) {
  # Computes the orthographic retraction of L (with corresponding
  # observations Y and threshold gamma) onto the manifold of
  # low-rank matrices via gradient descent with step size eta.
  #
  # approximation L
  # observations Y
  # step size eta
  # threshold gamma
  # sparsity matrix sparsity
  #
  # TO DO:
  #   * Implement matrix calculations in C++

  n1 <- nrow(Y); n2 <- ncol(Y)
  gradient <- threshold(L - Y, gamma, sparsity)

  Q <- L[,sample(n2, r)] # any r independent columns of L
  R <- t(L[sample(n1, r),]) # any r independent rows of L

  # descent step
  L_tmp <- L - eta * gradient
  L_out <- (L_tmp %*% R) %*% solve(crossprod(Q, L_tmp) %*% R) %*% crossprod(Q, L_tmp)

  # return solution and corresponding gradient
  list("L" = L_out, "gradient" = gradient)
}
