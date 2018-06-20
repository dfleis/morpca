orthographic.retraction <- function(L, Y, r, gamma, eta, n1, n2, sparsity) {
  # Computes the orthographic retraction of L (with corresponding
  # observations Y and threshold gamma) onto the manifold of
  # low-rank matrices via gradient descent with step size eta.
  #
  # approximation L
  # observations Y
  # low rank target r
  # threshold gamma
  # step size eta
  # nb of rows of the observation matrix n1
  # nb of cols of the observation matrix n2
  # sparsity matrix sparsity
  #
  # TO DO:
  #   * Implement matrix calculations in C++

  gradient <- threshold(L - Y, gamma, n1, n2, sparsity)

  Q <- L[,sample(n2, r)] # any r independent columns of L
  R <- t(L[sample(n1, r),]) # any r independent rows of L

  # descent step
  L_tmp <- L - eta * gradient
  L_out <- (L_tmp %*% R) %*% solve(crossprod(Q, L_tmp) %*% R) %*% crossprod(Q, L_tmp)

  # return solution and corresponding gradient
  list("L" = L_out, "gradient" = gradient)
}
