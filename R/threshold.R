threshold <- function(X, alpha, sparsity) {
  n1 <- nrow(X); n2 <- ncol(X)
  t1 <- rep(1, n1); t2 <- rep(1, n2)

  X <- X * sparsity

  for (i in 1:n1) {
    tt <- sort(abs(X[i,]), decreasing=TRUE)
    t1[i] <- tt[floor(alpha * sum(sparsity[i,])) + 1]
  }
  for (j in 1:n2) {
    tt <- sort(abs(X[,j]), decreasing = TRUE)
    t2[j] <- tt[floor(alpha * sum(sparsity[, j])) + 1]
  }

  threshold1 <- matrix(rep(1, n1 * n2), ncol = n2)
  threshold2 <- matrix(rep(1, n1 * n2), nrow = n1)
  threshold1 <- abs(X) <= matrix(rep(t1, each = n2), ncol = n2,byrow = TRUE)
  threshold2 <- abs(X) <= matrix(rep(t2, each = n1), nrow = n1)

  X_thresholded <- X * (as.double(threshold1 + threshold2) >= 1)

  X_thresholded
}

