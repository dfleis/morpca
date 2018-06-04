morpca <- function(Y, r, gamma,
                   retraction = c("projective", "orthogonal"),
                   step_size, step_max, tol = .Machine$double.eps) {
  # TO DO:
  #   * Set escape condition under a sufficient tolerance (i.e. 10^-10 or something)

  # Initialize best rank-r approximation of F(L - Y)
  # (? check this since I'm doing the best rank-r approx of Y)
  L0 <- init_approx(Y, r)
  L_list <- vector(mode = 'list')
  L_list[[1]] <- L0

  if (retraction[1] == "projective") {
    for (k in 1:(step_max - 1)) {
      L_list[[k + 1]] <- projective_retraction(L_list[[k]], Y, eta = step_size, gamma)
    }
  } else if (retraction[1] == "orthogonal") {
    ### TO DO
  } else {
    ### handle error
  }

  L_list
}
