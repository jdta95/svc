#include <RcppArmadillo.h>
#include "common_fxns.h" 

using namespace Rcpp;

// Function to update beta_r using Gibbs sampling
// [[Rcpp::export]]
double update_beta_r(double tau_square, double sigma_r_square, double phi_r, const arma::mat& Knots, const arma::mat& s) {
  // Calculate the covariance matrix K based on the locations and phi_r
  arma::mat C = calc_C(Knots, phi_r);
  
  // Calculate the inverse of the covariance matrix using Cholesky decomposition
  arma::mat C_inv = inv_Chol(C);

  // Calculate the variance parameter for the posterior distribution
  arma::mat post_var = inv_Chol(inv_Chol(tau_square*arma::eye(C.n_rows, C.n_rows)) + inv_Chol(sigma_r_square*C));

  // Return beta_r
  return calc_c(s,Knots, phi)*C_inv*arma::mvnrnd(arma::zeros(C.n_rows), post_var);
}
