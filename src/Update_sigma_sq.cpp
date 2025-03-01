#include <RcppArmadillo.h>
#include "common_fxns.h" 

using namespace Rcpp;

// Function to update sigma_r^2 using Gibbs sampling
// [[Rcpp::export]]
double update_sigma2_r(const arma::vec& beta_r, double a_r, double b_r, double phi_r, const arma::mat& Knots) {
  // Calculate the covariance matrix K based on the locations and phi_r
  arma::mat C = calc_C(Knots, phi_r);
  
  // Calculate the inverse of the covariance matrix using Cholesky decomposition
  arma::mat C_inv = inv_Chol(C);
  
  // Calculate (beta_r' C^{-1} beta_r)
  double beta_r_sum_square = arma::as_scalar(beta_r.t() * C_inv * beta_r);
  
  // Calculate the shape and scale parameters for the inverse gamma distribution
  double a_r_post = a_r + beta_r.n_elem / 2.0;
  double b_r_post = b_r + beta_r_sum_square / 2.0;
  
  // Sample sigma_r^2 from the inverse gamma distribution
  return 1.0 / R::rgamma(a_r_post, 1.0 / b_r_post);
}