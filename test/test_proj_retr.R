library(pracma) # generating random orthogonal matrices via randortho()

library(devtools)
#install_github("dfleis/morpca")
library(morpca)

#==========================#
#===== set parameters =====#
#==========================#
n1 <- 500 # rows
n2 <- 600 # columns
r <- 5 # r must be 0 < r <= min(n1, n2)

SIGMA <- diag(rep(1, r))
gamma <- 0.2

step_size <- 0.05 # step size
step_max <- 10 # max nb of steps

#========================#
#===== generate data ====#
#========================#
U <- randortho(n1)[,1:r]
V <- randortho(n2)[,1:r]

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
