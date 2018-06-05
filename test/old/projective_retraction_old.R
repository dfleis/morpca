#===========================================================#
#
# Recreating the projective retractive algorithm from
# 'Robust PCA by Manifold Optimization'
#
#===========================================================#

#===== functions =====#
F_thresh <- function(A, gamma) {
  # thresholding function in the objective

  row_pctl <- apply(abs(A), 1, quantile, gamma)
  col_pctl <- apply(abs(A), 2, quantile, gamma)

  # to do: think of a faster way to do this...
  A_out <- A
  for (i in 1:nrow(A)) {
    for (j in 1:ncol(A)) {
      if ((abs(A[i,j]) <= row_pctl[i]) & (abs(A[i,j]) <= col_pctl[j])) {
        A_out[i,j] <- 0
      }
    }
  }
  A_out
}
P <- function(L_k, D) {
  # calculates the Riemannian gradient of L given gradient D

  ### PRETTY SLOW!!!
  svd_out <- svd(L_k)
  U <- svd_out$u
  V <- svd_out$v

  ### VERY SLOW!!!!
  tcrossprod(U) %*% D + D %*% tcrossprod(V) - tcrossprod(U) %*% D %*% tcrossprod(V)

}
R_proj <- function(L_k, Y, eta, gamma) {
  # projective retraction
  #
  # approximation L_k
  # observations Y
  # stepsize eta
  # threshold gamma

  # handle partial observations
  Y_tmp <- Y
  Y_tmp[is.na(Y_tmp)] <- 0

  L_k - eta * P(L_k, F_thresh(L_k - Y, gamma))
}

init_approx <- function(Y, r) {
  # computes the best rank-r approximation of F_thresh(Y)
  # via the Eckhard-Young Theorem

  # to do... check if r is an integer
  svd_out <- svd(Y)

  U <- svd_out$u[,1:r]
  SIG <- diag(svd_out$d[1:r])
  V <- svd_out$v[,1:r]

  U %*% tcrossprod(SIG, V)
}

#===== parameters =====#
n1 <- 500 # rows
n2 <- 600 # columns
r <- 5 # r must be <= min(n1, n2)

SIGMA <- diag(rep(1, r))
gamma <- 0.2

eta <- 0.5 # step size
step_max <- 10 # max nb of steps

#===== generate data ====#
X <- matrix(rnorm(n1 * n2), nrow = n1)
svd_X <- svd(X)
U <- svd_X$u[,1:r]
V <- svd_X$v[,1:r]

Y0 <- U %*% tcrossprod(SIGMA, V)
Y <- apply(Y0, 2, function(y) {y[sample(n1, 25)] <- rnorm(25); y})

#==========================================#
#=============== START WORK ===============#
#==========================================#
### optimize
L0 <- init_approx(Y, r)
L_list <- vector(mode = 'list')
L_list[[1]] <- L0

for (k in 1:(step_max - 1)) {

  L_list[[k + 1]] <- R_proj(L_list[[k]], Y, eta, gamma)

}


### compute errors
err <- sapply(L_list, function(L) sum((L - Y)^2))


#===== plots =====#
plot(err, log = 'y', type = 'l')








