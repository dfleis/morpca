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

#====================#
#===== OPTIMIZE =====#
#====================#
pt <- proc.time()
L <- gradient_descent(Y = Y, r = r, alpha = alpha_bnd,
                      stepsize = 0.7, opt = 1, maxiter = 100,
                      sparsity = matrix(rep(1, nrow(Y) * ncol(Y)), nrow = nrow(Y)))[[1]]
proc.time() - pt


#================================#
#===== VISUALIZE IN/OUTPUTS =====#
#================================#
frame.idx <- 99

# inputs
input.tmp0 <- Y[frame.idx, 1:ncol(Y)]
input.tmp1 <- matrix(input.tmp0, nrow = 256, byrow = F)
input.plot <- t(apply(input.tmp1, 2, rev))
image(input.plot, col = colorRampPalette(c("black", "white"))(64),
      xaxt = 'n', yaxt = 'n')

# outputs
output.tmp0 <- L[frame.idx, 1:ncol(L)]
output.tmp1 <- matrix(output.tmp0, nrow = 256, byrow = F)
output.plot <- t(apply(output.tmp1, 2, rev))
image(output.plot, col = colorRampPalette(c("black", "white"))(64),
      xaxt = 'n', yaxt = 'n')



