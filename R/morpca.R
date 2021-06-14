#' Manifold Optimization for Robust PCA
#'
#' A robust PCA algorithm designed to extract an underlying low-rank matrix
#' with grossly corrupted observations.
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
#' To recover \eqn{L^*} we solve the following optimization problem:
#' \deqn{argmin_{rank(L) = r} f(L),}
#' where
#' \deqn{f(L) = 1/2 || F(L - Y) ||^2_F,}
#' such that \eqn{F} is the hard-thresholding function described above.
#'
#' @param Y Input matrix composed of the sum of matrices \eqn{L^*}
#'          (signal) and \eqn{S^*} (noise).
#' @param r Rank of the underlying matrix \eqn{L^*} and its estimate.
#' @param gamma Value between 0 and 1 corresponding to percentile
#'              for the hard thresholding procedure.
#' @param retraction String specifying which retraction technique
#'                   should be applied. Currently implemented are the
#'                   \code{"projective"} and \code{"orthographic"}
#'                   retractions.
#' @param stepsize Positive nonzero number corresponding to the step size
#'                 \eqn{\eta} in the gradient descent algorithm.
#' @param maxiter Positive nonzero integer specifying the maximum number
#'                of steps to compute for the gradient descent algorithm.
#'
#' @return To do...
#'
#' @export
morpca <- function(Y, r, gamma, stepsize, retraction = c("orthographic", "projective"), maxiter = 100) {
  # Data Y \in R^{n \times D} with additive decomposition Y = L^* + S^* for rank-d signal L^*
  # and sparse noise S^*
  retraction <- match.arg(retraction)

  n <- nrow(Y) # observations
  D <- ncol(Y) # variables/features

  # if missing L0 then initialize L...
  #L_list <- vector(mode = 'list', length = maxiter + 1)
  #L_list[[1]] <- low_rank_approx(threshold(Y, gamma), r)
  L <- low_rank_approx(threshold(Y, gamma), r)

  #err <- objective <- rep(NA, maxiter + 1)
  #err[1] <- 0.5 * norm(threshold(L_list[[1]] - Y, gamma), "f")^2
  for (k in 1:maxiter) {
    #L <- L_list[[k]]
    svdL <- svd(L)
    U <- svdL$u; SIGMA <- diag(svdL$d); V <- svdL$v
    D <- threshold(L - Y, gamma)
    if (retraction == "orthographic") {
      ### "OPTION 2"
      L <- (L - stepsize * D) %*% V %&% solve(t(U) %*% (L - stepsize * D) %*% V) %*% t(U) %*% (L - stepsize * D)
    } else if (retraction == "projective") {
      ### "OPTION 1"

    }


    #err[k + 1] <- 0.5 * norm(threshold(L - Y, gamma), "F")^2
    # print(proc.time() - pt)
  }

  return (list(L = L_list[[k]], S = Y - L_list[[k]], err = err, r = r, gamma = gamma, stepsize = stepsize))
}


Y <- matrix(rnorm(n * p), nrow = n)
L_list <- gradient_list <- vector(mode = 'list', length = maxiter + 1)
objective <- vector(mode = 'numeric', length = maxiter + 1)
n1 <- nrow(Y); n2 <- ncol(Y)


L <- threshold(Y, gamma)
SVD <- svd(tcrossprod(L))
U <- SVD$u

L_list[[1]] <- U[,1:r] %*% crossprod(U[,1:r], L)
gradient_list[[1]] <- L
objective[1] <- norm(gradient_list[[1]], "f")


for (k in 1:maxiter) {
  print(paste0("step = ", k))
  #tm <- system.time({
  if (retraction == "o") {
    out <- orthographic_retraction(L = L_list[[k]], Y = Y, r = r, gamma = gamma, eta = stepsize)
  } else if (retraction == "p") {
    out <- projective_retraction(L = L_list[[k]], Y = Y, r = r, gamma = gamma, eta = stepsize)
  } else {
    stop("Something went wrong... Invalid retraction.")
  }
  #})
  #print(tm)

  L_list[[k + 1]]        <- out$L
  gradient_list[[k + 1]] <- out$gradient
  objective[k + 1]       <- norm(gradient_list[[k + 1]], "f")
}

