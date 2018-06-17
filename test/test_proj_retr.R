library(Rcpp)
library(devtools)
#install_github("dfleis/morpca")
library(morpca)


#==========================#
#===== set parameters =====#
#==========================#
n1 <- 50 # rows
n2 <- 60 # columns
r <- 5 # r must be 0 < r <= min(n1, n2)

SIGMA <- diag(rep(1, r))
gamma <- 0.2

step_size <- 0.05 # step size
step_max <- 10 # max nb of steps

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

#=====================================#
#=============== TESTS ===============#
#=====================================#
pt <- proc.time()
L <- morpca(Y = Y, r = r, gamma = gamma,
            retraction = "projective",
            step_size  = step_size,
            step_max   = step_max,
            steps_out  = T,
            verbose    = T)
proc.time() - pt


#===================#
#===== FIGURES =====#
#===================#
#image(Y, col = colorRampPalette(c("white", "black"))(64))
#image(Lstar, col = colorRampPalette(c("white", "black"))(64))

err <- sapply(L, function(L_i) sqrt(sum((L_i - Lstar)^2)))
plot(err, type = 'l')
