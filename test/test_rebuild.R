#===== libraries ======#
library(morpca)
library(matrixvis)

#===== parameters =====#
# data parameters
n1 <- 50 # rows
n2 <- 60 # columns
r <- 5 # underlying signal matrix rank
nb_noise <- 5 # number of elements to be replaced with noise (columnwise)

# algorithm parameters
retraction <- "o" # orthographic retraction
gamma <- 0.2 # hard-thresholding percentile
stepsize <- 0.5 # step size
maxiter  <- 100 # max nb of steps

#===== generate data =====#
X <- matrix(rnorm(n1 * n2), nrow = n1, ncol = n2)

# generate rank-r data matrix Lstar
SIGMA <- diag(r)
svd_X <- svd(X)
U <- svd_X$u[,1:r]
V <- svd_X$v[,1:r]
Lstar <- U %*% SIGMA %*% t(V)

# generate input data matrix Y = Lstar + Sstar
# columnwise replace data with random noise to 25 random elements
Y <- apply(Lstar, 2, function(l)
  {l[sample(n1, nb_noise)] <- rnorm(nb_noise); l})

# sparse noise matrix
Sstar <- Y - Lstar

#===== solve =====#
pt <- proc.time()
L.opt <- morpca(Y = Y, r = r, gamma = gamma,
                sparsity   = sparsity_mat,
                retraction = retraction,
                stepsize   = stepsize,
                maxiter    = maxiter,
                stepsout   = T,
                verbose    = T)
proc.time() - pt


L.est <- L.opt$solution

x <- sapply(L.est, function(L) norm(L - Lstar, "f"))

plot(x, type = 'l', log = 'y')


#===== visualization =====#
bwr_col <- colorRampPalette(c("blue", "white", "red"))
image_matrix(Lstar, main = "Signal", col = bwr_col(32), legend = T)
image_matrix(Sstar, main = "Sparse Noise", col = bwr_col(32), legend = T)
image_matrix(Y, main = "Input Data", col = bwr_col(32), legend = T)
