morpca <- function(Y, r, gamma, sparsity,
                   retraction = c("projective", "orthographic"),
                   stepsize,
                   maxiter,
                   tol       = .Machine$double.eps, #to do...
                   stepsout  = F,
                   verbose   = F) {
  # TO DO:
  #   * Set escape condition under a sufficient tolerance (i.e. 10^-10 or something)
  #   * Handle partial observations (NA values)
  #   * Handle missing and invalid inputs
  #   * Make 'verbose' work more elegantly
  #   * Set default inputs

  # set up data structures
  L_list <- gradient_list <- vector(mode = 'list', length = maxiter + 1)
  n1 <- nrow(Y); n2 <- ncol(Y)

  ##################
  ### INITIALIZE ###
  ##################

  ### OLD:
  # initialize best rank-r approximation of Y
  # (? check this since the paper says to intialize with the best rank-r
  # approx. of f(Y) which appears to be a typo ?)
  #
  #L_list[[1]] <- rank_r_approx_cpp(Y, r)
  #gradient_list[[1]] <- percentile_threshold(L_list[[1]] - Y, gamma)

  ### NEW: (not sure what this starting point is?)
  L <- threshold(Y, gamma, n1, n2, sparsity)
  SVD <- svd(tcrossprod(L))
  U <- SVD$u

  L_list[[1]] <- U[,1:r] %*% crossprod(U[,1:r], L)
  gradient_list[[1]] <- L

  if (verbose) {
    print(
      paste0("k = 0. Gradient/objective value = ",
             sprintf("%.5g", norm(gradient_list[[1]], "f"), 5)))
  }

  if (retraction[1] == "projective" |
      retraction[1] == "proj" |
      retraction[1] == "p") {
    #########################
    # PROJECTIVE RETRACTION #
    #########################

    for (k in 1:maxiter) {
      out <- projective_retraction(L = L_list[[k]], Y = Y, r = r, gamma = gamma,
                                   eta = stepsize, n1 = n1, n2 = n2,
                                   sparsity = sparsity)
      L_list[[k + 1]]        <- out$L
      gradient_list[[k + 1]] <- out$gradient

      if (verbose) {
        print(
          paste0("k = ", k, ". Gradient/objective value = ",
                 sprintf("%.5g", norm(gradient_list[[k + 1]], "f"), 5)))
      }
    }

  } else if (retraction[1] == "orthographic" |
             retraction[1] == "orth" |
             retraction[1] == "o") {
    ###########################
    # ORTHOGRAPHIC RETRACTION #
    ###########################

    for (k in 1:maxiter) {
      out <- orthographic_retraction(L = L_list[[k]], Y = Y, r = r, gamma = gamma,
                                     eta = stepsize, n1 = n1, n2 = n2,
                                     sparsity = sparsity)
      L_list[[k + 1]]        <- out$L
      gradient_list[[k + 1]] <- out$gradient

      if (verbose) {
        print(
          paste0("k = ", k, ". Gradient/objective value = ",
                 sprintf("%.5g", norm(gradient_list[[k + 1]], "f"), 5)))
      }
    }

  } else {
    ############################
    # INVALID RETRACTION INPUT #
    ############################

    stop("Error in morpca(): Argument 'retraction' must be specified as either 'projective' or 'orthographic'.")
  }

  # return outputs
  # check if the user requests the full descent path output
  # or only the final low rank matrix
  if (stepsout) {
    list("solution"   = L_list,
         "gradient"   = gradient_list,
         "retraction" = retraction)
  } else {
    list("solution"   = L_list[[length(L_list)]],
         "gradient"   = gradient_list[[length(gradient_list)]],
         "retraction" = retraction)
  }
}
