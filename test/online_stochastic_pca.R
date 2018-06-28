library(expm)

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
n <- 9
d <- 4
k <- 2
eta <- 1
maxiter <- 10

#===== generate data =====#
X <- matrix(rnorm(n * d), nrow = n)

#===== initialize =====#
# generate PSD matrix M0 such that the largest eigenvalue 
# is <= 1/(d - k) and tr(M0) = 1
INIT_TMP <- matrix(rnorm(d^2), nrow = d, ncol = d)
SVD <- svd(INIT_TMP)
U0 <- SVD$u

maxsig <- 1/(d - k)
nb_maxsig <- 1/maxsig
SIG0 <- diag(c(rep(maxsig, nb_maxsig), rep(0, d - nb_maxsig)))

M0 <- U0 %*% SIG0 %*% t(U0)

M_list <- vector(mode = 'list', length = maxiter + 1)
M_list[[1]] <- M0

#===== start calculations =====#


pt <- proc.time()
for (i in 1:maxiter) {
  x <- X[sample(1:n, 1),]
  xxt <- tcrossprod(x)
  M_tmp <- expm(logm(M_list[[i]]) - eta * xxt)
  
  M_list[[i + 1]] <- project_re(M_tmp, d, k)
}
proc.time() - pt

# calculate objective in (5)
obj <- sapply(M_list, function(M) {
  obj_all <- apply(X, 1, function(x) {
    sum(diag(crossprod(M, tcrossprod(x))))
  })
  mean(obj_all)
})

#===== figures =====#
plot(obj, type = 'o', pch = 21, cex = 0.75, bg = 'white')





