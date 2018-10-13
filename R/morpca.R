#' Manifold Optimization for Robust PCA
#'
#' TO DO...
#'
#' @details Implementation of the robust PCA algorithms outlined in
#' Zhang, T. and Yang, Y. (forthcoming) to recover the underlying low rank
#' matrix \eqn{L^*} given an input matrix \eqn{Y} assumed to be of the form
#' \deqn{Y = L^* + S^*,} where \eqn{S^*} is a sparse matrix representing the
#' noise and \eqn{L^*} representing the signal. This implementation offers
#' both the projective and orthographic retractions as methods to map from
#' the input matrix's tangent space back to the manifold of low rank matrices.
#' Robustness is achieved via a hard thresholding function which sets
#' entry \eqn{(i,j)} to 0 if it exceeds the \eqn{\gamma}-th percentile of the
#' \eqn{i}-th row and \eqn{j}-th column (see Zhang, T. and Yang, Y. (forthcoming)
#' for details).
#'
#' TO DO... Talk about the objective function and a brief overview of
#' the gradient descent process.
#'
#' @param Y Input matrix composed of the sum of matrices \eqn{L^*}
#'          (signal) and \eqn{S^*} (noise).
#' @param r Rank of the underlying matrix \eqn{L^*} and its estimate.
#' @param gamma Value between 0 and 1 corresponding to percentile
#'              for the hard thresholding procedure.
#' @param sparsity TO DO...
#' @param retraction String specifying which retraction technique
#'                   should be applied. Currently implemented are the
#'                   \code{"projective"} and \code{"orthographic"}
#'                   retractions.
#' @param stepsize Positive nonzero number corresponding to the step size
#'                 \eqn{\eta} in the gradient descent algorithm.
#' @param maxiter Positive nonzero integer specifying the maximum number
#'                of steps to compute for the gradient descent algorithm.
#' @param stepsout Boolean value. If \code{stepsout = T} then the function
#'                 returns the output low rank matrix estimate and corresponding
#'                 gradient for every step in the gradient descent algorithm.
#' @param verbose Boolean value. If \code{verbose = T} then the gradient
#'                descent algorithm reports the current step number and
#'                value of the objective function at each iteration.
#'
#' @return Returns an object of class \code{morpca} with elements:
#'
#' \item{Y}{The original observations used as the input matrix for which we
#' seek an underlying low rank matrix.}
#' \item{rank}{Rank of the estimated target underlying matrix.}
#' \item{gamma}{Thresholding percentile.}
#' \item{sparsity}{TO DO...}
#' \item{retraction}{TO DO...}
#' \item{stepsize}{TO DO...}
#' \item{maxiter}{TO DO...}
#' \item{solution}{Estimated underlying low rank matrix. List of matrices if \code{stepsout = T}.}
#' \item{gradient}{Corresponding gradient matrix of the objective function. List of matrices if \code{stepsout = T}.}
#' \item{objective}{Corresponding value of the objective function (Frobenius norm of the gradient matrix).}
#'
#' @export
morpca <- function(Y = NULL, r = NULL, gamma = NULL, sparsity = NULL,
                   retraction = c("projective", "orthographic"),
                   stepsize  = NULL,
                   maxiter   = 100,
                   stepsout  = F,
                   verbose   = F) {
  # TO DO:
  #   * Set escape condition under a sufficient tolerance (i.e. 10^-10 or something)
  #     NOTE: tol should depend on the size of the input matrix,
  #           for example something like n1 * n2 * .Machine$double.eps
  #   * Handle partial observations (NA values)
  #   * Handle missing and invalid inputs
  #   * Determine what the default retraction ought to be
  #   * Make 'verbose' work more elegantly
  #   * Set default behavior for inputs
  #   * Rename input 'gamma' to 'threshold' NOTE: WE MUST RENAME THE CORRESPONDING
  #     threshold() FUNCTION!
  #   * Check if gamma = 0 or gamma = 1 and apply a special case/warning if
  #     such a scenario is specified?
  #   * Re-implement sparsity

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
  # convergence tolerance
  #if (is.null(tol)) { # (not actually implemented yet)
  #  # there's probably some theoretical guideline in the paper
  # # look through it to set something up...
  #  tol <- ncol(Y) * nrow(Y) * .Machine$double.eps
  #}

  # set up data structures
  L_list <- gradient_list <- vector(mode = 'list', length = maxiter + 1)
  objective <- vector(mode = 'numeric', length = maxiter + 1)
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

  ### NEW: (not sure what this starting point is/what is its theoretical
  # justification?)
  L <- threshold(Y, gamma, sparsity)
  SVD <- svd(tcrossprod(L))
  U <- SVD$u

  L_list[[1]] <- U[,1:r] %*% crossprod(U[,1:r], L)
  gradient_list[[1]] <- L
  objective[1] <- norm(gradient_list[[1]], "f")

  if (verbose) {
    print(paste0("k = 0. Objective value = ", sprintf("%.5g", objective[1], 5)))
  }

  #========================#
  # BEGIN GRADIENT DESCENT #
  #========================#
  out <- vector(mode = "list", length = 2)

  for (k in 1:maxiter) {

    if (retraction[1] == "projective" |
        retraction[1] == "proj" |
        retraction[1] == "p") {

      out <- projective_retraction(L = L_list[[k]], Y = Y, r = r, gamma = gamma,
                                   eta = stepsize, sparsity = sparsity)

    } else if (retraction[1] == "orthographic" |
               retraction[1] == "orth" |
               retraction[1] == "o") {

      out <- orthographic_retraction(L = L_list[[k]], Y = Y, r = r, gamma = gamma,
                                     eta = stepsize, sparsity = sparsity)

    } else { # invalid 'retraction' input specified
      stop("Error in morpca(): Argument 'retraction' must be specified as either 'projective' or 'orthographic'.")
    }

    L_list[[k + 1]]        <- out$L
    gradient_list[[k + 1]] <- out$gradient
    objective[k + 1]       <- norm(gradient_list[[k + 1]], "f")

    if (verbose) {
     print(
       paste0("k = ", k, ". Objective value = ",
              sprintf("%.5g", objective[k + 1], 5)))
    }

    if (!stepsout) { # clear memory
       L_list[[k]]        <- NA
       gradient_list[[k]] <- NA
    }
  }

  # remove erased elements in output if stepsout = T
  L_list        <- L_list[!is.na(L_list)]
  gradient_list <- gradient_list[!is.na(gradient_list)]

  #================#
  # RETURN OUTPTUS #
  #================#
  out <- list("Y"          = Y,
              "rank"       = r,
              "gamma"      = gamma,
              "sparsity"   = sparsity,
              "retraction" = retraction,
              "stepsize"   = stepsize,
              "maxiter"    = maxiter,
              "solution"   = L_list[!sapply(L_list, is.null)],
              "gradient"   = gradient_list[!sapply(gradient_list, is.null)],
              "objective"  = objective,
              "call"       = match.call())
  # check if the user requests the full descent path output
  # or only the solution for the final low rank matrix
  if (!stepsout) {
    out$solution <- out$solution[[length(out$solution)]]
    out$gradient <- out$gradient[[length(out$gradient)]]
  }
  class(out) <- "morpca"
  out
}
