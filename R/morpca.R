morpca <- function(Y, r, gamma, sparsity,
                   retraction = c("projective", "orthographic"),
                   stepsize  = NULL,
                   maxiter   = 100,
                   tol       = .Machine$double.eps, #to do...
                   stepsout  = F,
                   verbose   = F) {
  # TO DO:
  #   * Set escape condition under a sufficient tolerance (i.e. 10^-10 or something)
  #   * Handle partial observations (NA values)
  #   * Handle missing and invalid inputs
  #   * Make 'verbose' work more elegantly
  #   * Set default behavior for inputs
  #   * Create morpca type object in a list and output this object
  #   * Rename input 'gamma' to 'threshold' NOTE: WE MUST RENAME THE CORRESPONDING
  #     threshold() FUNCTION!

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

  #========================#
  # BEGIN GRADIENT DESCENT #
  #========================#
  if (retraction[1] == "projective" |
      retraction[1] == "proj" |
      retraction[1] == "p") {

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

  } else { # invalid 'retraction' input specified
    stop("Error in morpca(): Argument 'retraction' must be specified as either 'projective' or 'orthographic'.")
  }

  #================#
  # RETURN OUTPTUS #
  #================#

  out <- list("data"       = Y,
              "rank"       = r,
              "gamma"      = gamma,
              "sparsity"   = sparsity,
              "retraction" = retraction,
              "stepsize"   = stepsize,
              "maxiter"    = maxiter,
              "tol"        = tol,
              "solution"   = L_list,
              "gradient"   = gradient_list)

  # check if the user requests the full descent path output
  # or only the final low rank matrix
  if (stepsout) {
    out
  } else {
    out$solution <- out$solution[[length(L_list)]]
    out$gradient <- out$gradient[[length(gradient_list)]]
    out
  }
}
