morpca <- function(Y, r, gamma,
                   retraction = c("projective", "orthogonal"),
                   step_size, step_max, tol = .Machine$double.eps,
                   steps_out = T) {
  # TO DO:
  #   * Set escape condition under a sufficient tolerance (i.e. 10^-10 or something)

  # Initialize best rank-r approximation of Y
  # (? check this since the paper says to intialize with the best rank-r
  # approx. of f(Y) which appears to be a typo ?)
  L0 <- rank_r_approx_cpp(Y, r)
  L_list <- vector(mode = 'list')
  L_list[[1]] <- L0

  if (retraction[1] == "projective") {
    for (k in 1:(step_max - 1)) {
      L_tmp   <- projective_retraction(L_list[[k]], Y, eta = step_size, gamma)
      L_list[[k + 1]] <- rank_r_approx_cpp(L_tmp, r)
    }
  } else if (retraction[1] == "orthogonal") {
    ### TO DO
  } else {
    ### handle error
  }

  if (steps_out) {
    L_list
  } else {
    L_list[[length(L_list)]]
  }
}
