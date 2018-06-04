// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "helpers.h"

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat riemann_gradient_cpp(arma::mat L, 
							   arma::mat D) {
	// Computes the Riemann gradient of a matrix L 
	// given Euclidean gradient D
	
	arma::mat U;
	arma::vec s;
	arma::mat V;
	
	svd_econ(U, s, V, L);
	
	// is there a faster way to do this?
	return U * U.t() * D + D * V * V.t() - U * U.t() * D * V * V.t();
}