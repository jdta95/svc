#include <RcppArmadillo.h>
#include "common_fxns.h"

// [[Rcpp::export]]
Rcpp::List GP_Gibbs(
    arma::vec Y,               // Response vector
    arma::mat X,               // Design matrix
    arma::mat s,               // Spatial locations (n x 2 matrix)
    arma::mat knots,           // Knot locations (m x 2 matrix)
    double sigmasq_start,      // Initial value for sigma^2
    double tausq_start,        // Initial value for tau^2
    arma::vec phi_start,       // Initial value for phi
    arma::vec w_start,         // Initial value for w
    arma::vec beta_start,      // Initial value for beta
    double a_r,                // Prior shape parameter for sigma^2
    double b_r,                // Prior scale parameter for sigma^2
    double a_t,                // Prior shape parameter for tau^2
    double b_t,                // Prior scale parameter for tau^2
    arma::vec lower,           // Lower bound for phi
    arma::vec upper,           // Upper bound for phi
    int mcmc = 1000            // Number of MCMC iterations
) {
  int n = Y.n_elem;            // Number of observations
  int p = X.n_cols;            // Number of covariates
  int m = knots.n_rows;        // Number of knots
  
  // Initialize storage for posterior samples
  arma::mat beta_samples = arma::zeros(mcmc, p);
  arma::mat w_samples = arma::zeros(mcmc, n);
  arma::vec sigmasq_samples = arma::zeros(mcmc);
  arma::vec tausq_samples = arma::zeros(mcmc);
  arma::mat phi_samples = arma::zeros(mcmc, phi_start.n_elem);
  
  // Initialize current values
  arma::vec beta_cur = beta_start;
  arma::vec w_cur = w_start;
  double sigmasq_cur = sigmasq_start;
  double tausq_cur = tausq_start;
  arma::vec phi_cur = phi_start;
  
  // Adaptive MCMC for phi
  arma::mat phi_metrop_sd = 0.05 * arma::eye(phi_start.n_elem, phi_start.n_elem);
  RAMAdapt phi_adapt(phi_start.n_elem, phi_metrop_sd, 0.234);
  
  // Gibbs sampling loop
  for (int i = 0; i < mcmc; i++) {
    // Update phi 
    phi_cur = phi_RW(knots, w_cur, sigmasq_cur, phi_cur(0), phi_adapt.paramsd, arma::join_horiz(lower, upper));
    phi_samples.row(i) = phi_cur.t();
    
    // Update beta 
   // beta_cur = beta_update(Y, X, tausq_cur, w_cur, V);
  //  beta_samples.row(i) = beta_cur.t();
    
    // Update w (placeholder for group member's function)
  //  w_cur = w_update(Y, X, beta_cur, sigmasq_cur, tausq_cur, calc_C_tilde(s, knots, phi_cur(0)), n);
  //  w_samples.row(i) = w_cur.t();
    
    // Update sigma^2 using your function
    sigmasq_cur = update_sigma2_r(beta_cur, a_r, b_r, phi_cur, knots);
    sigmasq_samples(i) = sigmasq_cur;
    
    // Update tau^2 using your function
    tausq_cur = update_tau2_r(Y, X, beta_cur, w_cur, a_t, b_t);
    tausq_samples(i) = tausq_cur;
    
    // Adapt phi proposal variance
    phi_adapt.update_ratios();
    if (i % 10 == 0) {
      phi_adapt.adapt(arma::randn(phi_start.n_elem), 0.234, i);
    }
  }
  
  // Return posterior samples
  Rcpp::List output;
  output["beta"] = beta_samples;
  output["w"] = w_samples;
  output["sigmasq"] = sigmasq_samples;
  output["tausq"] = tausq_samples;
  output["phi"] = phi_samples;
  
  return output;
} 