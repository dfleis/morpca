orthographic_retraction <- function(L, Y, r, gamma, eta, n1, n2, sparsity) {
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
  #   * Check if n1 > n2 or n2 < n1 since we can speed up
  #     multiplications by knowing which case we're in

  gradient <- threshold(L - Y, gamma, n1, n2, sparsity)

  Q <- L[,sample(n2, r)] # any r independent columns of L
  R <- L[sample(n1, r),] # any r independent rows of L

  # descent step
  L_tmp <- L - eta * gradient
  QtL_tmp <- crossprod(Q, L_tmp)
  L_out <- tcrossprod(L_tmp, R) %*% solve(tcrossprod(QtL_tmp, R)) %*% QtL_tmp

  #L_out <- orthographic_retraction_cpp(L_tmp, Q, R)


  # return solution and corresponding gradient
  list("L" = L_out, "gradient" = gradient)
}
