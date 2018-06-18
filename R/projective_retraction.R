projective_retraction <- function(L, Y, eta, gamma, sparsity) {
  # Computes the projective retraction of L (with corresponding
  # observations Y and threshold gamma) onto the manifold of
  # low-rank matrices via gradient descent with step size eta.
  #
  # approximation L
  # observations Y
  # step size eta
  # threshold gamma
  # sparsity matrix sparsity

  # compute descent step
  gradient <- threshold(L - Y, gamma, sparsity)
  L_tmp <- L - eta * riemann_gradient_cpp(L, gradient)
  L_out <- rank_r_approx_cpp(L_tmp, r)

  # return solution and corresponding gradient
  list("L" = L_out, "gradient" = gradient)
}
