#library(devtools)
#install_github("dfleis/morpca")
library(morpca)

#==========================#
#===== set parameters =====#
#==========================#
retraction <- "o"
n1 <- 500 # rows
n2 <- 600 # columns
r <- 5 # r must be 1 <= r <= min(n1, n2)

SIGMA <- diag(r)
gamma <- 0.2 # threshold

stepsize <- 0.5 # step size
maxiter  <- 100 # max nb of steps

#========================#
#===== generate data ====#
#========================#
X <- matrix(rnorm(n1 * n2), nrow = n1)

# generate orthogonal matrices U (n1 x r) and V (n2 x r)
svd_X <- svd(X)
U <- svd_X$u[,1:r]
V <- svd_X$v[,1:r]

Lstar <- U %*% SIGMA %*% t(V)

Y <- apply(Lstar, 2, function(l) {l[sample(n1, 25)] <- rnorm(25); l})
Sstar <- Y - Lstar

sparsity_mat <- matrix(1, nrow = n1, ncol = n2)

#=====================================#
#=============== TESTS ===============#
#=====================================#
pt <- proc.time()
L.opt <- morpca(Y = Y, r = r, gamma = gamma,
                sparsity   = sparsity_mat,
                retraction = retraction,
                stepsize   = stepsize,
                maxiter    = maxiter,
                stepsout   = T,
                verbose    = T)
proc.time() - pt

# objective function (Frobenius norm of the gradient matrix)
obj <- L.opt$objective
# error from true (low-rank) signal
err <- sapply(L.opt$solution, function(L) norm(L - Lstar, "f"))

#===================#
#===== FIGURES =====#
#===================#
#image(Y, col = colorRampPalette(c("white", "black"))(64))
#image(Lstar, col = colorRampPalette(c("white", "black"))(64))

plot(obj, type = 'l', log = 'y', main = "Objective Value",
     sub = "Frob. Norm of the Gradient Matrix")
plot(err, type = 'l', log = 'y', main = "Error from True Signal (Frob. Norm)")
