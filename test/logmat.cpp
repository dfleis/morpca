#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat logmat_cpp(arma::mat M) {
  return arma::logmat_sympd(M);
}

