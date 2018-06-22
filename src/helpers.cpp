// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat riemann_gradient_cpp(arma::mat L, arma::mat D) {
	// Computes the Riemann gradient of a matrix L 
	// given Euclidean gradient D
	
	arma::mat U, UUt, V, VVt;
	arma::vec s;
	
	svd_econ(U, s, V, L);
	
	// is there a faster way to do this?
	UUt = U * U.t();
	VVt = V * V.t();
	return UUt * D + D * VVt - UUt * D * VVt;
}

// [[Rcpp::export]]
arma::mat orthographic_retraction_cpp(arma::mat L_tmp, arma::mat Q, arma::mat R) {
	// computes the descent step for the orthographic retraction given
	// L_tmp = L^{(k)} - eta * gradient (Euclidean descent)
	// Q = any r independent columns of L^{(k)} (dimension n1 x r)
	// R = any r independent rows of L^{(k)} (dimension n2 x r)
	
	// is there a faster way of doing this?
	arma::mat QtL_tmp, QtL_tmpR_inv;
	
	QtL_tmp = Q.t() * L_tmp;
	QtL_tmpR_inv = inv(QtL_tmp * R);
	
	return L_tmp * R * QtL_tmpR_inv * QtL_tmp;
}
							   
// [[Rcpp::export]]
arma::mat rank_r_approx_cpp(arma::mat Y, int r) {
	// Computes the best rank-r approximation of matrix Y (via the
	// Eckheart-Young Theorem?)
	
	arma::mat U; arma::mat U_r;
	arma::vec s; arma::vec s_r; arma::mat S_r;
	arma::mat V; arma::mat V_r;
	
	svd_econ(U, s, V, Y);
	
	U_r = U.submat(0, 0, U.n_rows - 1, r - 1);	
	s_r = s.subvec(0, r - 1);
	V_r = V.submat(0, 0, V.n_rows - 1, r - 1);
	
	S_r = diagmat(s_r);
	
	return U_r * S_r * V_r.t();
}









