library(devtools)
install_github("dfleis/morpca")

library(morpca)

#==========================#
#===== set parameters =====#
#==========================#
n1 <- 500 # rows
n2 <- 600 # columns
r <- 5 # r must be 0 < r <= min(n1, n2)

SIGMA <- diag(rep(1, r))
gamma <- 0.2

step_size <- 1 # step size
step_max <- 500 # max nb of steps

#========================#
#===== generate data ====#
#========================#
X <- matrix(rnorm(n1 * n2), nrow = n1)
svd_X <- svd(X)
U <- svd_X$u[,1:r]
V <- svd_X$v[,1:r]

Lstar <- U %*% SIGMA %*% t(V)
Y <- apply(Lstar, 2, function(y) {y[sample(n1, 25)] <- rnorm(25); y})

#=====================================#
#=============== TESTS ===============#
#=====================================#
pt <- proc.time()
L <- morpca(Y = Y, r = r, gamma = gamma,
            retraction = "projective",
            step_size  = step_size,
            step_max   = step_max,
            steps_out  = T)
proc.time() - pt

err <- sapply(L, function(l) sqrt(sum((l - Lstar)^2)))
plot(err, log = 'y', type = 'l')


