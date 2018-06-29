#library(expm)
library(matrixcalc) # for is.positive.semi.definite() to check PSD while testing
library(Rcpp)
library(microbenchmark)
sourceCpp("test/logmat.cpp")
sourceCpp("test/expmat.cpp")

#===== functions =====#
project_re <- function(M, d, k) {
  eig <- eigen(M)
  V <- eig$vectors
  s <- eig$values
  
  S_new <- sapply(0:(d - 1), function(k_prime) {
    s[1:k_prime] <- 1/(d - k_prime) 
    s[(k_prime + 1):d] <- 
      (1 - sum(s[1:k_prime])) * s[(k_prime + 1):d]/sum(s[(k_prime + 1):d])
    s
  })
  
  S_all_idx <- which(apply(S_new, 2, function(x) sum(x <= 1/(d - k))) == d)
  S_idx <- max(S_all_idx)
  s_new <- S_new[,S_idx]
  
  V %*% tcrossprod(diag(s_new), V)
}

#===== parameters =====#
set.seed(124)
n <- 100
d <- 5
k <- 3

maxiter <- n

#===== generate data =====#
X <- matrix(rnorm(n * d), nrow = n)
I <- diag(d)

#===== initialize =====#
# generate PSD matrix M0 such that the largest eigenvalue 
# is <= 1/(d - k) and tr(M0) = 1
INIT_TMP <- matrix(rnorm(d^2), nrow = d, ncol = d)
SVD <- svd(INIT_TMP)
U0 <- SVD$u

#maxsig <- 1/(d - k)
#nb_maxsig <- 1/maxsig
#SIG0 <- diag(c(rep(maxsig, nb_maxsig), rep(0, d - nb_maxsig)))

SIG0 <- diag(abs(rnorm(d)))
M0 <- round(U0 %*% SIG0 %*% t(U0), 6)

M_list <- vector(mode = 'list', length = maxiter + 1)
M_list[[1]] <- M0

#===== start calculations =====#
pt <- proc.time()
for (i in 1:maxiter) {
  eta <- sqrt(1/i)
  x <- X[i,]
  
  xxt <- tcrossprod(x)
  #M_tmp <- expmat_cpp(logmat_cpp(M_list[[i]]) - eta * xxt)
  M_tmp <- expmat_cpp(expm::logm(M_list[[i]], method = "Eigen") - eta * xxt)
  
  M_list[[i + 1]] <- project_re(M_tmp, d, k)
}
proc.time() - pt

# calculate objective in (5)
obj <- sapply(M_list, function(M) {
  obj_all <- apply(X, 1, function(x) {
    sum(diag(M %*% tcrossprod(x)))
  })
  mean(obj_all)
})


M <- M_list[[maxiter + 1]]
Q <- qr.Q(qr(M))
R <- qr.R(qr(M))
chol(I - M)
rankMatrix(M)


#===== figures =====#
plot(obj[2:length(obj)], type = 'o', pch = 21, cex = 0.75, bg = 'white')

sapply(M_list, rankMatrix)

