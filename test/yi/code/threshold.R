threshold <- function(X,alpha,sparsity){
    N <- nrow(X)
    D <- ncol(X)
    X <- X*sparsity
    t1 <- rep(1,N)
    t2 <- rep(1,D)
    for (i in 1:N){   
         tt <- sort(abs(X[i,]),decreasing=TRUE)
         t1[i] <- tt[floor(alpha*sum(sparsity[i,]))+1]
     }
  
    for (j in 1:D){  
         tt <- sort(abs(X[,j]),decreasing=TRUE)
         t2[j] <- tt[floor(alpha*sum(sparsity[,j]))+1]
         
     }
  threshold1 <- matrix(rep(1,N*D),ncol=D)
  threshold2 <- matrix(rep(1,N*D),nrow=N)
  threshold1 <- abs(X)<=matrix(rep(t1,each=D),ncol=D,byrow = TRUE)
  threshold2 <- abs(X)<=matrix(rep(t2,each=N),nrow=N)
X_thresholded <- X*(as.double(threshold1+threshold2)>=1)
return(X_thresholded)
}

