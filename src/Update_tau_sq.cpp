#include <RcppArmadillo.h>
#include "common_fxns.h" 

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Function to update tau^2 using Gibbs sampling
// [[Rcpp::export]]
double update_tau2_r(const arma::vec& Y, const arma::mat& X, const arma::vec& beta, const arma::vec& w, double a_t, double b_t) {
  arma::vec residual = Y - X * beta - w;
  double residual_sum_square = arma::dot(residual, residual);
  double a_t_post = a_t + Y.n_elem / 2.0;
  double b_t_post = b_t + residual_sum_square / 2.0;
 
   //Sample tau^2 from the inverse gamma distribution
  return 1.0 / arma::randg(arma::distr_param(a_t_post, 1.0 / b_t_post));
}