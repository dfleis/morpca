projective_retraction <- function(L, Y, eta, gamma) {
  # Computes the projective retraction of L (with corresponding
  # observations Y and threshold gamma) onto the manifold of
  # low-rank matrices via gradient descent with step size eta.
  #
  # approximation L
  # observations Y
  # step size eta
  # threshold gamma

  # compute descent step
  L_tmp <- L - eta * riemann_gradient_cpp(L, percentile_threshold(L - Y, gamma))

  # return the best rank-r approximation of the descent step
  rank_r_approx_cpp(L_tmp, r)
}
