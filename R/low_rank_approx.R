low_rank_approx <- function(M, r) {
  # returns the best rank-r matrix to M (wrt the Frob. norm)
  # as determined by Eckart-Young theorem

  # in the case of our morpca implementation, this procedure
  # can likely be optimized in such a way that doesn't require
  # computing the svd each time (unless we only use it during
  # the initialization phase, then it doesn't matter)
  svdM <- svd(M, nu = r, nv = r)
  return (svdM$u %*% diag(svdM$d[1:r]) %*% t(svdM$v))
}
