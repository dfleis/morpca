library(devtools)
install_github("dfleis/morpca")

library(morpca)

#==========================#
#===== set parameters =====#
#==========================#
n1 <- 100 # rows
n2 <- 100 # columns
r <- 5 # r must be <= min(n1, n2)

SIGMA <- diag(rep(1, r))
gamma <- 0.2

step_size <- 0.5 # step size
step_max <- 10 # max nb of steps

#========================#
#===== generate data ====#
#========================#
X <- matrix(rnorm(n1 * n2), nrow = n1)
svd_X <- svd(X)
U <- svd_X$u[,1:r]
V <- svd_X$v[,1:r]

Y0 <- U %*% tcrossprod(SIGMA, V)
Y <- apply(Y0, 2, function(y) {y[sample(n1, 25)] <- rnorm(25); y})

#=====================================#
#=============== TESTS ===============#
#=====================================#
L <- morpca(Y = Y, r = r, gamma = gamma,
            retraction = "projective",
            step_size  = step_size,
            step_max   = step_max)


