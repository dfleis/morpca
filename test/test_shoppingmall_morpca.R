#=====================#
#===== LIBRARIES =====#
#=====================#
#library(devtools)
#install_github("dfleis/morpca")
library(morpca)

library(Matrix)
library(ggplot2)
library(R.matlab) # reading .mat files

#=====================#
#===== LOAD DATA =====#
#=====================#
#shoppingmall <- read.csv("test/yi/data/shoppingmall.csv", header = F, check.names = F)
shoppingmall_mat <- readMat("test/yi/data/shoppingmall.mat")

Y <- as.matrix(shoppingmall_mat$X)

#==========================#
#===== SET PARAMETERS =====#
#==========================#
set.seed(124)

r <- 3
gamma <- 0.1
retraction <- "o"
stepsize <- 0.7
maxiter <- 100

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

sparsity <- matrix(1, nrow = nrow(Y), ncol = ncol(Y))

#====================#
#===== OPTIMIZE =====#
#====================#
pt <- proc.time()
L.opt <- morpca(Y = Y, r = r, gamma = gamma, sparsity = sparsity,
                retraction = retraction, stepsize = stepsize,
                maxiter = maxiter, stepsout = T, verbose = T)
proc.time() - pt



#================================#
#===== VISUALIZE IN/OUTPUTS =====#
#================================#
opt.idx <- 2
frame.idx <- 2

L <- L.opt$solution[[opt.idx]]

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
      xaxt = 'n', yaxt = 'n', main = opt.idx)

# difference
diff.plot <- input.plot - output.plot
image(diff.plot, col = colorRampPalette(c("black", "white"))(64),
      xaxt = 'n', yaxt = 'n', main = opt.idx)


grad_norm <- sapply(L.opt$gradient, function(D) norm(D, "f"))
err_norm <- sapply(L.opt$solution, function(l)
  norm(threshold(l - Y, gamma = gamma, sparsity = sparsity), "f"))
plot(grad_norm, type = 'o', pch = 21, bg = 'white',
     cex = 0.75, main = "Gradient Norm")
plot(err_norm, type = 'o', pch = 21, bg = 'white',
     cex = 0.75, main = "Objective")
