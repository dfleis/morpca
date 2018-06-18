#library(Rcpp)
#library(devtools)
#install_github("dfleis/morpca")
library(morpca)

#==========================#
#===== set parameters =====#
#==========================#
retraction <- "orthographic"
n1 <- 500 # rows
n2 <- 600 # columns
r <- 5 # r must be 0 < r <= min(n1, n2)

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

grad1 <- sapply(L.opt$gradient, function(grad) norm(grad, "f"))
grad2 <- grad1/norm(Y, "f")
err <- sapply(L.opt$solution, function(L) norm(L - Lstar, "f"))

plot(grad1, type = 'l', log = 'y')
plot(grad2, type = 'l', log = 'y')
plot(err, type = 'l', log = 'y')

#===================#
#===== FIGURES =====#
#===================#
#image(Y, col = colorRampPalette(c("white", "black"))(64))
#image(Lstar, col = colorRampPalette(c("white", "black"))(64))

err <- sapply(L, function(L_i) sqrt(sum((L_i - Lstar)^2)))
plot(err, type = 'l')
