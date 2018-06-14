library(Matrix)
library(ggplot2)
set.seed(1)

shoppingmall <- read.csv("shoppingmall.csv",header=FALSE,check.names=FALSE)
Y <- as.matrix(shoppingmall)



r <- 3
alpha_bnd <- 0.1
d1 <- nrow(Y)
d2 <- ncol(Y)

p <- 0.5
ncol <- floor(p*d1)
n <- d2*ncol
I0 <- rep(0,n)
J0 <- rep(0,n)
X0 <- rep(0,n)
sparsity <- matrix(rep(0,d1*d2),nrow=d1)

for (j in 1:d2){

   I0[((j-1)*ncol+1):(j*ncol)] <- sort(sample(d1,ncol))
   J0[((j-1)*ncol+1):(j*ncol)] <- j

   X0[((j-1)*ncol+1):(j*ncol)] <- Y[I0[((j-1)*ncol+1):(j*ncol)],j]
   sparsity[I0[((j-1)*ncol+1):(j*ncol)],j] <- 1
}
Y0 <- sparseMatrix(I0,J0,x=X0)
Y1 <- as.matrix(Y0)

#allstepsizes <- c(0.05,0.1,0.2,0.4,0.7,1,1.5,2.5)


L <- gradient_descent(Y,r,alpha_bnd,0.7,1,100,matrix(rep(1,nrow(Y)*ncol(Y)),nrow=nrow(Y)))[[1]]

  
Lp <- gradient_descent( Y1,r,alpha_bnd,0.7/p,1,100,sparsity)[[1]]




allcases <- c(50,100,700)



images <- matrix(rep(1,(266*length(allcases)-10)*(330*5-10)),nrow=266*length(allcases)-10)
for (i in 1:length(allcases)){
  ind <- allcases[i]
  images[(1:256)+266*(i-1),(1:320)+330*(0)] <- matrix(Y[ind,],nrow=256,ncol=320)/255
  images[(1:256)+266*(i-1),(1:320)+330*(1)] <- matrix(L[ind,],nrow=256,ncol=320)/255
  images[(1:256)+266*(i-1),(1:320)+330*(2)] <- matrix(Lp[ind,],nrow=256,ncol=320)/255
  images[(1:256)+266*(i-1),(1:320)+330*(3)] <- matrix(Y[ind,]-L[ind,],nrow=256,ncol=320)/255
  images[(1:256)+266*(i-1),(1:320)+330*(4)] <- matrix(Y[ind,]-Lp[ind,],nrow=256,ncol=320)/255
}



rotate <- function(x) t(apply(x, 2, rev))
image(rotate(images),axes = FALSE, col = grey(seq(0, 1, length = 256)))



