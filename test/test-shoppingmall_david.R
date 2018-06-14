#=====================#
#===== LIBRARIES =====#
#=====================#
library(Rcpp) # loading helpers.cpp
library(Matrix)
library(ggplot2)
library(R.matlab) # reading .mat files

library(morpca)
sourceCpp("src/helpers.cpp")
source("test/yi/code/threshold.R")
source("test/yi/code/gradient_descent.R")

#=====================#
#===== LOAD DATA =====#
#=====================#
#shoppingmall <- read.csv("test/yi/data/shoppingmall.csv", header = F, check.names = F)
shoppingmall_mat <- readMat("test/yi/data/shoppingmall.mat")

Y <- as.matrix(shoppingmall_mat$X)

#==========================#
#===== SET PARAMETERS =====#
#==========================#
set.seed(1)
r <- 3
alpha_bnd <- 0.1
p <- 0.5

#====================================#
#===== GENERATE DATA STRUCTURES =====#
#====================================#
n1 <- nrow(Y)
n2 <- ncol(Y)
ncol <- floor(p * n1)
n <- n2 * ncol

I0 <- rep(0, n)
J0 <- rep(0, n)
X0 <- rep(0, n)
sparsity <- matrix(rep(0, n1 * n2), nrow = n1)

for (j in 1:n2) {
  I0[((j - 1) * ncol + 1):(j * ncol)] <- sort(sample(n1, ncol))
  J0[((j - 1) * ncol + 1):(j * ncol)] <- j

  X0[((j - 1) * ncol + 1):(j * ncol)] <- Y[I0[((j - 1) * ncol + 1):(j * ncol)], j]
  sparsity[I0[((j - 1) * ncol + 1):(j * ncol)], j] <- 1
}

Y0 <- sparseMatrix(I0, J0, x = X0)
Y1 <- as.matrix(Y0)

#===============================#
#===== VISUALIZE THE INPUT =====#
#===============================#
# visualizing the shoppingmall data
idx <- 1
Mtmp1 <- Y1[idx, 1:ncol(Y)]
Mtmp <- matrix(Mtmp1, nrow = 256, byrow = F)
M <- t(apply(Mtmp, 2, rev))
image(M, col = colorRampPalette(c("black", "white"))(64),
      xaxt = 'n', yaxt = 'n')

#====================#
#===== OPTIMIZE =====#
#====================#
pt <- proc.time()
L <- gradient_descent(Y = Y, r = r, alpha = alpha_bnd,
                      stepsize = 0.7, opt = 1, maxiter = 100,
                      sparsity = matrix(rep(1, nrow(Y) * ncol(Y)), nrow = nrow(Y)))[[1]]
proc.time() - pt





