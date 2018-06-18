projective_retraction <- function(L, Y, r, gamma, eta, n1, n2, sparsity) {
  # Computes the projective retraction of L (with corresponding
  # observations Y and threshold gamma) onto the manifold of
  # low-rank matrices via gradient descent with step size eta.
  #
  # approximation L
  # observations Y
  # threshold gamma
  # step size eta
  # nb of rows of the observation matrix n1
  # nb of cols of the observation matrix n2
  # sparsity matrix sparsity

  # compute descent step
  gradient <- threshold(L - Y, gamma, n1, n2, sparsity)
  L_tmp <- L - eta * riemann_gradient_cpp(L, gradient)
  L_out <- rank_r_approx_cpp(L_tmp, r)

  # return solution and corresponding gradient
  list("L" = L_out, "gradient" = gradient)
}
