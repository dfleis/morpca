morpca <- function(Y, r, gamma, sparsity,
                   retraction = c("projective", "orthographic"),
                   stepsize,
                   maxiter,
                   tol       = .Machine$double.eps,
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
  L <- threshold(Y, gamma, sparsity)
  SVD <- svd(tcrossprod(L))
  U <- SVD$u

  L_list[[1]] <- U[,1:r] %*% crossprod(U[,1:r], L)
  gradient_list[[1]] <- L

  if (verbose) {
    print(
      paste0("k = 0. Gradient/objective value = ",
             round(norm(gradient_list[[1]], "f"), 5)))
  }

  if (retraction[1] == "projective") {
    #########################
    # PROJECTIVE RETRACTION #
    #########################

    for (k in 1:maxiter) {
      out <- projective_retraction(L_list[[k]], Y, stepsize, gamma, sparsity)
      L_list[[k + 1]]        <- out$L
      gradient_list[[k + 1]] <- out$gradient

      if (verbose) {
        print(
          paste0("k = ", k, ". Gradient/objective value = ",
                 round(norm(gradient_list[[k + 1]], "f"), 5)))
      }
    }

  } else if (retraction[1] == "orthographic") {
    ###########################
    # ORTHOGRAPHIC RETRACTION #
    ###########################

    # for (k in 1:maxiter) {
    #   L_list[[k + 1]] <- orthographic_retraction(L_list[[k]], Y, stepsize, gamma)
    #   gradient_list[[k + 1]] <- percentile_threshold(L_list[[k + 1]] - Y, gamma)
    #
    #   if (verbose) {
    #     print(
    #       paste0("k = ", k, ". Gradient/objective value = ",
    #              round(norm(gradient_list[[k + 1]], "f"), 5)))
    #   }
    # }

    for (k in 1:maxiter) {
      out <- orthographic_retraction(L_list[[k]], Y, stepsize, gamma, sparsity)
      L_list[[k + 1]]        <- out$L
      gradient_list[[k + 1]] <- out$gradient

      # L <- L_list[[k]]
      # U <- L[,sample(n2,r)]
      # V <- t(L[sample(n1,r),])
      # #gradient <- threshold(L - Y, gamma, sparsity)
      # gradient <- percentile_threshold(L - Y, gamma)
      # L1 <- L - stepsize * gradient
      #
      # L_list[[k + 1]] <- (L1 %*% V) %*% solve((t(U) %*% L1 %*% V)) %*% (t(U) %*% L1)
      # gradient_list[[k + 1]] <- gradient

      if (verbose) {
        print(
          paste0("k = ", k, ". Gradient/objective value = ",
                 round(norm(gradient_list[[k + 1]], "f"), 5)))
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
