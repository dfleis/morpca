// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "helpers.h"

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec GetProxOne(arma::vec y,
                     arma::vec weights) {
  // This function evaluates the prox given by
  // 1/2 * || y - beta ||_2^2 + P(\beta), where
  // P(\beta) = \sum_{j=1}^{p} weights[j] * norm(beta[j:p]).
  //
  // NOTE: This is an internal function used for the additivehierbasis function
  //      only. The difference between this and GetProx is that this function
  //      solves the prox for a SINGLE lambda value.
  //
  // Args:
  //    y: The given input vector y.
  //    weights: The vector of weights used in the penalty function.
  // Returns:
  //    beta: The solution vector of the optimization problem.
  //
  // Initialize the dimension of vector y.
  int p = y.size();

  // Initialize the beta vector which will be returned.
  arma::vec beta(y.begin(), y.size(), true);


  // main loop for solving the proximal problem
  for(int j = p - 1; j >= 0; --j) {
    // the first term corresponds to the scaling factor
    beta.subvec(j, p - 1) =
      max(1 - weights(j) / norm(beta.subvec(j, p - 1)), 0) *
      beta.subvec(j, p - 1);
  }

  return beta;
}

// [[Rcpp::export]]
List FitAdditive(arma::vec y,
                 arma::mat weights, arma::vec ak,
                 NumericVector x,
                 arma::mat beta,
                 double max_lambda, double lam_min_ratio,
                 double alpha,
                 double tol, int p, int J, int n,
                 int nlam, double max_iter,
                 bool beta_is_zero, arma::vec active_set) {

  // Initialize some objects.
  IntegerVector dimX = x.attr("dim");
  arma::cube X(x.begin(), dimX[0], dimX[1], dimX[2]);
  arma::cube x_mats(n, J, p);
  arma::cube r_mats(J, J, p);
  arma::vec max_lam_values(p);

  // This loop does the QR decompositon and generates the Q, R matrices.
  // It also helps us find the maximum lambda value when it is not specified.
  for(int i = 0; i < p; ++i) {
    arma::mat temp_x_mat, temp_r_mat;

    // Perform an 'economial' QR decomposition.
    arma::qr_econ(temp_x_mat, temp_r_mat, X.slice(i));

    // Generate the x_mat and the r_mat.
    temp_x_mat = temp_x_mat * sqrt(n);
    temp_r_mat = temp_r_mat / sqrt(n);

    x_mats.slice(i) = temp_x_mat;
    r_mats.slice(i) = temp_r_mat;

    // If max_lambda = NULL, then we select the maximum lambda value ourselves.
    if(R_IsNA(max_lambda)) {
      arma::vec v_temp = temp_x_mat.t() * (y/n);
      if(R_IsNA(alpha)) {
        arma::vec temp_lam_max = sqrt(abs(v_temp)/ ak);

        // This is obtained by solving the inequality
        // lambda^2 + lambda >= |v_1|.
        temp_lam_max(0) = 0.5 * (sqrt(4 * fabs(v_temp(0)) + 1) - 1);
        max_lam_values(i) = max(temp_lam_max);
      } else {
        arma::vec temp_lam_max =  abs(v_temp)/ (ak * alpha);

        temp_lam_max(0) = fabs(v_temp(0));
        max_lam_values(i) = max(temp_lam_max);
      }

    }
  }

  if(R_IsNA(max_lambda)) {
    max_lambda = max(max_lam_values);
  }

  // Generate the full lambda sequence.
  arma::vec lambdas = linspace<vec>(log10(max_lambda),
                                    log10(max_lambda * lam_min_ratio),
                                    nlam);
  lambdas = exp10(lambdas);

  // Generate matrix of weights.
  if(!R_IsNA(alpha)) {
    weights.each_row() %= alpha * lambdas.t();
    weights.row(0) = weights.row(0) + (1 - alpha) * lambdas.t();
  } else {
    weights.each_row() %= pow(lambdas.t(), 2);
    weights.row(0) = weights.row(0) + lambdas.t();
  }


  // If the user left initial beta == NULL, then we don't need to
  // do all the matrix multiplication.
  arma::mat x_beta(n, p, fill::zeros);
  if(!beta_is_zero) {
    for(int i = 0; i < p; ++i) {
      x_beta.col(i) = x_mats.slice(i) * beta.col(i);
    }
  }

  // Turns out WE CAN VECTORIZE AND MAKE
  // ONE BIG SPARSE MATRIX.
  arma::cube beta_ans(J, p, nlam);

  // Initialize some vectors and matrices.
  arma::vec temp_weights;
  arma::vec temp_y;
  arma::vec temp_v;
  arma::vec temp_beta_j;

  arma::vec temp_vec_beta;

  double temp_norm_old;
  double change;



  // Begin main loop for each value of lambda.
  for(int i = 0; i < nlam; i++) {
    // Rcout << "nlam: " << i<<"\n";
    temp_weights = weights.col(i) ;
    int  counter = 0;
    bool converged_final = false;
    while(counter < max_iter && !converged_final) {

      // We will use this to check convergence.
      arma::mat old_beta(beta.begin(), J, p, true);

      // One loop of the block coordinate descent algorithm.
      for(int j = 0; j < p; j++) {
        if(active_set(j) != 0) {
          temp_y = y - sum(x_beta, 1) + (x_mats.slice(j) * beta.col(j));
          temp_v = trans(x_mats.slice(j)) *  (temp_y / n);
          temp_beta_j = GetProxOne(temp_v, temp_weights);
          // Update the vector x_beta (X_j\beta_j).
          x_beta.col(j) = x_mats.slice(j) * temp_beta_j;
        } else {
          temp_beta_j = zeros(J);
          // Update the vector x_beta (X_j\beta_j).
          x_beta.col(j) = zeros(n);
        }
        // Update the matrix beta.
        beta.col(j) = temp_beta_j;
      }

      temp_vec_beta = vectorise(old_beta);

      // Obtain the value of the relative change.
      temp_norm_old = norm(temp_vec_beta);
      change = norm(vectorise(beta)) - temp_norm_old;
      // Rcout << fabs(change) << "\n";
      if(fabs(change) < tol) {
        beta_ans.slice(i) = beta;
        converged_final = true;

        // One loop of the block coordinate descent algorithm.
        // To update the active set and check for final convergence.
        for(int j = 0; j < p; j++) {
          temp_y = y - sum(x_beta, 1) + (x_mats.slice(j) * beta.col(j));
          temp_v = trans(x_mats.slice(j)) *  (temp_y / n);
          temp_beta_j = GetProxOne(temp_v, temp_weights);
          // Update the vector x_beta (X_j\beta_j).
          x_beta.col(j) = x_mats.slice(j) * temp_beta_j;

          // Update the matrix beta.
          beta.col(j) = temp_beta_j;

          if(any(temp_beta_j) !=0 && active_set(j) == 0) {
            active_set(j) = 1;
            converged_final = false;
          }
        }

        // If we have not converged but only updated the active set then we increment
        // the counter.
        if(!converged_final){
          counter = counter + 1;
        }

      } else {
        counter = counter + 1;

        if(counter == max_iter) {
          Function warning("warning");
          warning("Function did not converge");
        }

      }

    }
  }

  arma::sp_mat beta_final(p * J, nlam);
  for(int i = 0; i < p; i++) {
    arma::mat temp_slice = beta_ans.tube(0, i, J-1, i);
    beta_ans.tube(0, i, J-1, i) = solve(trimatu(r_mats.slice(i)), temp_slice);
  }

  for(int i = 0; i < nlam; i++) {
    beta_final.col(i) = vectorise(beta_ans.slice(i));
  }

  return List::create(Named("beta") = beta_final,
                      Named("lambdas") = lambdas);
}
