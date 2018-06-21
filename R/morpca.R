#' Manifold Optimization for Robust PCA
#'
#' @details Implementation of the robust PCA algorithms outlined in
#' Zhang, T. and Yang, Y. (forthcoming) to recover the underlying low rank
#' matrix \eqn{L^*} given an input matrix \eqn{Y} assumed to be of the form
#' \deqn{Y = L^* + S^*,} where \eqn{S^*} is a sparse matrix representing the
#' noise and \eqn{L^*} representing the signal. This implementation offers
#' both the projective and orthographic retractions as methods to map from
#' the input matrix's tangent space back to the manifold of low rank matrices.
#'
#' @param Y Input matrix composed of the sum of matrices \eqn{L^*}
#'          (signal) and \eqn{S^*} (noise).
#' @param r Rank of the underlying matrix \eqn{L^*} and its estimate.
#' @param gamma Value between 0 and 1 corresponding to percentile
#'              for the hard thresholding procedure.
#' @param sparsity
#' @param retraction String specifying which retraction technique
#'                   should be applied. Currently implemented are the
#'                   \code{"projective"} and \code{"orthographic"}
#'                   retractions.
#' @param stepsize Positive nonzero number corresponding to the step size \eqn{\eta}
#'                 in the gradient descent algorithm.
#' @param maxiter Positive integer specifying the maximum number of steps
#'                compute for the gradient descent algorithm.
#' @param tol
#' @param stepsout Boolean value. If \code{stepsout = T} then the function
#'                 returns the output low rank matrix estimate and corresponding
#'                 gradient for every step in the gradient descent algorithm.
#' @param verbose Boolean value. If \code{verbose = T} then the gradient
#'                descent algorithm reports the current step number and
#'                value of the objective function at each iteration.
#'
#' @return Returns an object of class \code{morpca} with elements: TO DO...
#'
#' @export
morpca <- function(Y = NULL, r = NULL, gamma = NULL, sparsity = NULL,
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
  #   * Determine what the default retraction ought to be
  #   * Make 'verbose' work more elegantly
  #   * Set default behavior for inputs
  #   * Rename input 'gamma' to 'threshold' NOTE: WE MUST RENAME THE CORRESPONDING
  #     threshold() FUNCTION!

  if (is.null(Y)) {
    # return warning/set default? return error?
  }
  if (is.null(r)) {
    # return warning/set default? return error?
  }
  if (is.null(gamma)) {
    # return warning/set default? return error?
  }
  if (is.null(sparsity)) {
    # return warning/set default? return error?
  }
  if (is.null(stepsize)) {
    # return warning/set default? return error?
    # check if stepsize <= 0
  }

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
  out <- list()
  out <- list("data"       = Y,
              "rank"       = r,
              "gamma"      = gamma,
              "sparsity"   = sparsity,
              "retraction" = retraction,
              "stepsize"   = stepsize,
              "maxiter"    = maxiter,
              "tol"        = tol,
              "solution"   = L_list,
              "gradient"   = gradient_list,
              "call"       = match.call())
  # check if the user requests the full descent path output
  # or only the solution for the final low rank matrix
  if (!stepsout) {
    out$solution <- out$solution[[length(L_list)]]
    out$gradient <- out$gradient[[length(gradient_list)]]
  }
  class(out) <- "morpca"
  out
}
