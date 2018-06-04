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
  tcrossprod(U) %*% D + D %*% tcrossprod(V) - tcrossprod(U) %*% D %*% tcrossprod(V)

}
