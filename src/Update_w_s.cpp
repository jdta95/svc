#include <RcppArmadillo.h>
#include "common_fxns.h" 

using namespace Rcpp;

// Function to update w_s using Gibbs sampling
// [[Rcpp::export]]
double update_w_s(double tau_square, double sigma_w_square, double phi_w, const arma::mat& Knots) {
  // Calculate the covariance matrix K based on the locations and phi_w
  arma::mat C = calc_C(Knots, phi_w);
  
  // Calculate the inverse of the covariance matrix using Cholesky decomposition
  arma::mat C_inv = inv_Chol(C);

  // Calculate the variance parameter for the posterior distribution
  arma::mat post_var = inv_Chol(inv_Chol(tau_square*arma::eye(C.n_rows, C.n_rows)) + inv_Chol(sigma_w_square*C));

  // Return beta_r
  return arma::mvnrnd(arma::zeros(C.n_rows), post_var);
}
