#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat expmat_cpp(arma::mat M) {
  return arma::expmat_sym(M);
}

