// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <cmath> // ceil() and floor()

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

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
arma::mat percentile_threshold_cpp(arma::mat A, arma::mat A_abs, arma::vec row_pctls, arma::vec col_pctls) {
    // gamma-th percentile thresholding of matrix A given row and 
    // column percentiles								   
    int nrow = A.n_rows;
    int ncol = A.n_cols;							   
    
    arma::mat A_out = A;
    double Aij_abs;
    
    for (int i = 0; i <	nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            Aij_abs = A_abs(i,j);		
            
            if ((Aij_abs > row_pctls(i)) & (Aij_abs > col_pctls(j))) {
                A_out(i,j) = 0;
            }
        }
    }
    return A_out;							   
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

// [[Rcpp::export]]
double percentile_cpp(arma::vec x, double prob) { 
    // Computes the prob-th percentile of a vector x
    // analogously to quantile() in R
    
    int n = x.size();
    double index = 1 + (n - 1) * prob;
    int lo = std::floor(index);
    int hi = std::ceil(index);
    x = arma::sort(x);
    //std::partial_sort(x.begin(), x.begin() + hi, x.end());
    double h = index - lo;
    
    return (1 - h) * x(lo - 1) + h * x(hi - 1);
}

// [[Rcpp::export]]
arma::vec row_pctls_cpp(arma::mat A, double prob) {
    // computes the rowwise prob-th percentiles of a matrix A
    
    int nrow = A.n_rows;
    
    arma::vec row_pctls;
    row_pctls.zeros(nrow); // look for more elegant way to initialize
    
    arma::mat At = A.t(); // can we do this without taking the transpose?
    
    for (int i = 0; i < nrow; i++) {
        row_pctls(i) = percentile_cpp(At.col(i), prob);
    }
    
    return row_pctls;
}

// [[Rcpp::export]]
arma::vec col_pctls_cpp(arma::mat A, double prob) {
    // computes the columnwise prob-th percentiles of a matrix A
    int ncol = A.n_cols;
    
    arma::vec col_pctls;
    col_pctls.zeros(ncol); // look for more elegant way to initialize
    
    for (int j = 0; j < ncol; j++) {
        col_pctls(j) = percentile_cpp(A.col(j), prob);
    }
    
    return col_pctls;
}