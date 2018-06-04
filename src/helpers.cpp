// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath> // std::abs


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

// [[Rcpp::export]]
arma::mat percentile_threshold_cpp(arma::mat A,
								   arma::vec row_pctls,
								   arma::vec col_pctls) {
	int nrow = A.n_rows;
	int ncol = A.n_cols;							   
	
	arma::mat A_out = A;	
	
	for (int i = 0; i <	nrow; i++) {
		for (int j = 0; j < ncol; j++) {
			double Aij_abs = std::abs(A(i,j));		
		
			if ((Aij_abs <= row_pctls(i)) & (Aij_abs <= col_pctls(j))) {
				A_out(i,j) = 0;
			}
		}
	}
	return A_out;							   
}								   
