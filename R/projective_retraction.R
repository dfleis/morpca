projective_retraction <- function(L, Y, gamma, eta) {
  # Computes the projective retraction of Y (?... check this)
  #
  # approximation L
  # observations Y
  # threshold gamma
  # stepsize eta

  # handle partial (NA) observations
  #Y_tmp <- Y
  #Y[is.na(Y_tmp)] <- 0

  return (L - eta * riemann_gradient(L, threshold(L - Y, gamma)))
}
riemann_gradient <- function(L, D) {
  # Calculates the Riemannian gradient of L with Euclidean gradient D
  #
  # To do:
  #   * Implement in C++
  #   * Error handling

  ### PRETTY SLOW!!!
  svd_out <- svd(L)
  U <- svd_out$u
  V <- svd_out$v

  ### VERY SLOW!!!!
  return (tcrossprod(U) %*% D + D %*% tcrossprod(V) - tcrossprod(U) %*% D %*% tcrossprod(V))
}
