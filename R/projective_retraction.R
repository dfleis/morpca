projective_retraction <- function(L, Y, eta, gamma) {
  # Computes the projective retraction of Y (?... check this)
  #
  # approximation L
  # observations Y
  # stepsize eta
  # threshold gamma

  # handle partial (NA) observations
  Y_tmp <- Y
  Y_tmp[is.na(Y_tmp)] <- 0

  L - eta * riemann_gradient(L, percentile_threshold(L - Y, gamma))
}
