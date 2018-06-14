gradient_descent <- function(Y, r, alpha, stepsize, opt, maxiter, sparsity) {

    L <- threshold(Y, alpha, sparsity)
    SVD <- svd(L %*% t(L))

    U <- SVD$u
    S <- diag(SVD$d)
    V <- SVD$v

    L <- U[,1:r]%*%(t(U[,1:r])%*%L)
    iter <- 1
    normY <- norm(Y,"f")
    n1 <- nrow(L)
    n2 <- ncol(L)
    #all_L <- list()
    alle <- list()
    #all_L[iter] <- L

    if(opt==1){
        iter <- 1

         while (iter<maxiter){
          print(iter)
          U <- L[,sample(n2,r)]
          V <- t(L[sample(n1,r),])
          gradient <- threshold( L-Y,alpha,sparsity )
          L1 <- L-stepsize*gradient
          L <- L1%*%V%*%solve((t(U)%*%L1%*%V))%*%(t(U)%*%L1)
          iter <- iter+1
          #all_L[iter] <- L
          alle[iter] <- norm(gradient,"f")/normY
       }
        }
        else{
          iter <- 1

          while (iter<maxiter){
              gradient <- threshold( L-Y,alpha,sparsity )
              projected_gradient <- gradient%*%V[,1:r]%*%t(V[,1:r])+U[,1:r]%*%t(U[,1:r])%*%gradient-U[,1:r]%*%t(U[,1:r])%*%gradient%*%V[,1:r]%*%t(V[,1:r])
              L1 <- L-stepsize*projected_gradient
              SVDL <- svd(L1)
              U <- SVDL$u
              S <- diag(SVDL$d)
              V <- SVDL$v

              L <- U[,1:r]%*%S[1:r,1:r]%*%t(V[,1:r])
              iter <- iter+1
              #all_L[iter] <- L
              alle[iter] <- norm(L-truth,"f")

          }
       }

    return(list(L,alle))
}
