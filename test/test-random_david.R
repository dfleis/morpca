#=====================#
#===== LIBRARIES =====#
#=====================#
library(Rcpp) # loading helpers.cpp
library(Matrix)

library(morpca)
sourceCpp("src/helpers.cpp")
source("test/yi/code/threshold.R")
source("test/yi/code/gradient_descent.R")

#==========================#
#===== SET PARAMETERS =====#
#==========================#
set.seed(1)
n1 <- 500
n2 <- 600
r <- 5

alpha_bnd <- 0.2 # percentile threshold
p <- 0.5 # pct of corrupted/missing data

SIGMA <- diag(r)

#=========================#
#===== GENERATE DATA =====#
#=========================#

X <- matrix(rnorm(n1 * n2), nrow = n1)

# generate orthogonal matrices U (n1 x r) and V (n2 x r)
svd_X <- svd(X)
U <- svd_X$u[,1:r]
V <- svd_X$v[,1:r]

Lstar <- U %*% SIGMA %*% t(V)

Y <- apply(Lstar, 2, function(l) {l[sample(n1, 25)] <- rnorm(25); l})

#====================#
#===== OPTIMIZE =====#
#====================#
pt <- proc.time()
L.opt <- gradient_descent(Y = Y, r = r, alpha = alpha_bnd,
                          stepsize = 0.7, opt = 1, maxiter = 100,
                          sparsity = matrix(rep(1, nrow(Y) * ncol(Y)), nrow = nrow(Y)))
proc.time() - pt

plot(unlist(L.opt[[2]]), type = 'o', pch = 21, bg = 'white', cex = 0.5, log = 'y')

#=====================#
#===== VISUALIZE =====#
#=====================#

Lstar.plot <- t(apply(Lstar, 2, rev))
image(Lstar.plot, col = colorRampPalette(c("black", "white"))(32),
      xaxt = "n", yaxt = "n", main = expression(L^"*"))


Y.plot <- t(apply(Y, 2, rev))
image(Y.plot, col = colorRampPalette(c("black", "white"))(32),
      xaxt = "n", yaxt = "n", main = "Y")






