#include <RcppArmadillo.h>
#include "common_fxns.h"

// [[Rcpp::export]]
Rcpp::List GP_Gibbs(
    const arma::vec& Y,               // Response vector
    const arma::mat& X,               // Design matrix
    const arma::mat& s,               // Spatial locations (n x 2 matrix)
    const arma::mat& knots,           // Knot locations (m x 2 matrix)
    arma::mat beta_knots_start,         // Initial value for coefficients: beta 1, ..., beta p
    arma::vec w_knots_start,          // Initial value for w
    arma::vec phi_beta_start,       // Initial value for phi: beta 1, ..., beta p
    double phi_w_start,         // Initial value for phi: w
    arma::vec sigmasq_beta_start,      // Initial value for sigma^2: beta 1, ..., beta p
    double sigmsq_w_start,       // Initial value for sigma^2: w
    double tausq_start,        // Initial value for tau^2
    const arma::vec& phi_beta_proposal_sd,       // Proposal standard deviation for phi: beta 1, ..., beta p
    double phi_w_proposal_sd,        // Proposal standard deviation for phi: w
    const arma::vec& a_beta,                // Prior shape parameter for sigma^2: beta 1, ..., beta p
    const arma::vec& b_beta,                // Prior scale parameter for sigma^2: beta 1, ..., beta p
    double a_w,                // Prior shape parameter for sigma^2: w
    double b_w,               // Prior scale parameter for sigma^2: w
    double a_t,                // Prior shape parameter for tau^2
    double b_t,                // Prior scale parameter for tau^2
    const arma::vec& lower_beta,           // Lower bound for phi: beta 1, ..., beta p
    const arma::vec& upper_beta,           // Upper bound for phi: beta 1, ..., beta p
    const arma::vec& lower_w,             // Lower bound for phi: w
    const arma::vec& upper_w,             // Upper bound for phi: w
    int mcmc = 1000            // Number of MCMC iterations
) {
  int n = Y.n_elem;            // Number of observations
  int p = X.n_cols;            // Number of covariates
  int m = knots.n_rows;        // Number of knots
  
  // Initialize storage for posterior samples
  arma::cube beta_samples = arma::zeros(mcmc, n, p);
  arma::mat w_samples = arma::zeros(mcmc, n);
  arma::mat phi_beta_samples = arma::zeros(mcmc, p); // phi: beta 1, ..., beta p
  arma::vec phi_w_samples = arma::zeros(mcmc); // phi: w
  arma::mat sigmasq_beta_samples = arma::zeros(mcmc, p); // sigma^2: beta 1, ..., beta p
  arma::vec sigmasq_w_samples = arma::zeros(mcmc); // sigma^2: w
  arma::vec tausq_samples = arma::zeros(mcmc); // tau^2
  
  // Initialize current values
  arma::mat beta_knots_cur = beta_knots_start;
  arma::vec w_knots_cur = w_knots_start;
  arma::vec phi_beta_cur = phi_beta_start;
  double phi_w_cur = phi_w_start;
  arma::vec sigmasq_beta_cur = sigmasq_beta_start;
  double sigmasq_w_cur = sigmsq_w_start;
  double tausq_cur = tausq_start;
  
  // Create lower and upper bounds matrices
  arma::mat phi_beta_bounds = arma::join_horiz(lower_beta, upper_beta);
  arma::mat phi_w_bounds = arma::join_horiz(lower_w, upper_w);
  
  // Gibbs sampling loop
  for (int i = 0; i < mcmc; i++) {
    // Update phi 
    phi_w_cur = phi_RW(
      knots, w_knots_cur, sigmasq_w_cur, phi_w_cur, phi_w_proposal_sd, phi_w_bounds
    );
    phi_w_samples(i) = phi_w_cur;
    
    for (int j = 0; j < p; j++) {
      phi_beta_cur(j) = phi_RW(
        knots, beta_knots_cur.col(j), sigmasq_beta_cur(j), phi_beta_cur(j), phi_beta_proposal_sd(j), phi_beta_bounds.row(j));
    }
    phi_beta_samples.row(i) = phi_beta_cur.t();

    // Update beta 
    for (int j = 0; j < p; j++) {
      beta_knots_cur.col(j) = update_beta_r(Y, X, sigmasq_beta_cur(j), phi_beta_cur(j), w_knots_cur, tausq_cur, knots);
    }
    
    
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