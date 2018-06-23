#
# Initial exploration/replication of the stochastic (online) 
# optimization method for PCA described by Arora, R., et al. (2012).
#
# NOTE: Both the 'expm' and 'Matrix' libraries contain functions named
#       expm(), so be sure to specify we wish to use expm::expm() since
#       it appears to be considerably faster.
# NOTE: Recall that we can manually perform PCA by computing the centered
#       eigendecomposition, i.e. (up to a sign difference +/-1)
# Xc <- scale(X, scale = F)
# n <- nrow(Xc)
# V <- eigen(1/d * crossprod(Xc))$vectors
# scores <- Xc %*% V
# pca_X <- princomp(X)
# scores 
# pca_X$scores
#

library(microbenchmark)
library(expm)

#===========================#
#===== define functions ====#
#===========================#
update_step <- function(X, M, eta) {
  P_re(expm::expm(expm::logm(M) - eta * crossprod(X)))
}
P_re <- function(A) {
  # projects M onto the desired contraints with respect to the 
  # relative entropy  (see: Arora, R., et al. (2012). 
  # "Stochastic Optimization for PCA and PLS" for details)
  
}

#==========================#
#===== set parameters =====#
#==========================#
set.seed(125)
d <- 5 # number of variables
k <- 3  # number of observations

maxiter <- 100
tol <- .Machine$double.eps
eta <- 0.7

#=========================#
#===== generate data =====#
#=========================#
X <- matrix(rnorm(k * d), nrow = k)

#=============================#
#===== begin computations ====#
#=============================#
Xc <- scale(X, scale = F)
eigen_Xc <- eigen(crossprod(1/d * Xc))
U <- eigen_Xc$vectors

M0 <- diag(1/(d - k), d)/sum(diag(M0))
norm(M0, "2") <= 1/(d - k)

M.list <- vector(mode = "list", length = maxiter + 1)
M.list[[1]] <- M0

pt <- proc.time()
for (i in 1:(maxiter - 1)) {
  print(i)
  M.list[[i + 1]] <- proj_re(X, M.list[[i]], eta)
}
proc.time() - pt

M.list <- M.list[!sapply(M.list, is.null)]
x <- sapply(M.list, function(m) norm(m - U, "f"))
plot(x, type = 'l', log = 'y')


lapply(M.list, function(m) sum(diag(m)))

