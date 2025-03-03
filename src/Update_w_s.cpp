//#include <RcppArmadillo.h>
//#include "common_fxns.h" 

//// [[Rcpp::depends(RcppArmadillo)]]
//using namespace Rcpp;

// Function to update w_s using Gibbs sampling
// // [[Rcpp::export]]
//arma::vec update_w_s(const arma::vec& Y, const arma::mat& X, const arma::mat& beta_knots, double sigmasq_w, double phi_w, double tausq, const arma::mat& knots) {
  // Calculate the covariance matrix K based on the locations and phi_w
  //arma::mat C = calc_C(knots, phi_w);
  
  // Calculate the inverse of the covariance matrix using Cholesky decomposition
  //arma::mat C_inv = inv_Chol(C);
  
  // Calculate the variance parameter for the posterior distribution

  //arma::mat post_var = inv_Chol(inv_Chol(tausq * arma::eye(C.n_rows, C.n_rows)) + inv_Chol(sigmasq_w * C));
  
  // Calculate the mean parameter for the posterior distribution
  //arma::vec post_mean = post_var * (C_inv * (Y - X * beta_knots));
  
  // Sample w_s from the posterior distribution
  //return post_mean + arma::chol(post_var) * arma::randn(C.n_rows);
//}

