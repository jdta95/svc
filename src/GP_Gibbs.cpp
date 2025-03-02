#include <RcppArmadillo.h>
#include "common_fxns.h"


// [[Rcpp::export]]
Rcpp::List GP_Gibbs(
    const arma::vec& Y,               // Response vector
    const arma::mat& X,               // Design matrix
    const arma::mat& s,               // Spatial locations (n x 2 matrix)
    const arma::mat& knots,           // Knot locations (m x 2 matrix)
    arma::mat beta_knots_start,       // Initial value for coefficients: beta 1, ..., beta p
    arma::vec w_knots_start,          // Initial value for w
    arma::vec phi_beta_start,         // Initial value for phi: beta 1, ..., beta p
    double phi_w_start,               // Initial value for phi: w
    arma::vec sigmasq_beta_start,     // Initial value for sigma^2: beta 1, ..., beta p
    double sigmasq_w_start,           // Initial value for sigma^2: w
    double tausq_start,               // Initial value for tau^2
    const arma::vec& phi_beta_proposal_sd, // Proposal standard deviation for phi: beta 1, ..., beta p
    double phi_w_proposal_sd,         // Proposal standard deviation for phi: w
    const arma::vec& a_beta,          // Prior shape parameter for sigma^2: beta 1, ..., beta p
    const arma::vec& b_beta,          // Prior scale parameter for sigma^2: beta 1, ..., beta p
    double a_w,                       // Prior shape parameter for sigma^2: w
    double b_w,                       // Prior scale parameter for sigma^2: w
    double a_t,                       // Prior shape parameter for tau^2
    double b_t,                       // Prior scale parameter for tau^2
    const arma::vec& lower_beta,      // Lower bound for phi: beta 1, ..., beta p
    const arma::vec& upper_beta,      // Upper bound for phi: beta 1, ..., beta p
    const arma::vec& lower_w,         // Lower bound for phi: w
    const arma::vec& upper_w,         // Upper bound for phi: w
    int mcmc = 1000                   // Number of MCMC iterations
) {
  int n = Y.n_elem;                   // Number of observations
  int p = X.n_cols;                   // Number of covariates
  int m = knots.n_rows;               // Number of knots
  
  // Initialize storage for posterior samples
  arma::cube beta_samples = arma::zeros(mcmc, n, p);  // Real beta at observed locations
  arma::mat w_samples = arma::zeros(mcmc, n);         // Real w at observed locations
  arma::mat beta_knots_samples = arma::zeros(mcmc, m, p); // beta* at knot locations
  arma::mat w_knots_samples = arma::zeros(mcmc, m);       // w* at knot locations
  arma::mat phi_beta_samples = arma::zeros(mcmc, p);  // phi_beta samples
  arma::vec phi_w_samples = arma::zeros(mcmc);        // phi_w samples
  arma::mat sigmasq_beta_samples = arma::zeros(mcmc, p); // sigma^2_beta samples
  arma::vec sigmasq_w_samples = arma::zeros(mcmc);    // sigma^2_w samples
  arma::vec tausq_samples = arma::zeros(mcmc);        // tau^2 samples
  
  // Initialize current values
  arma::mat beta_knots_cur = beta_knots_start;        // Current beta* at knots
  arma::vec w_knots_cur = w_knots_start;              // Current w* at knots
  arma::vec phi_beta_cur = phi_beta_start;            // Current phi_beta
  double phi_w_cur = phi_w_start;                     // Current phi_w
  arma::vec sigmasq_beta_cur = sigmasq_beta_start;    // Current sigma^2_beta
  double sigmasq_w_cur = sigmasq_w_start;             // Current sigma^2_w
  double tausq_cur = tausq_start;                     // Current tau^2
  
  // Create lower and upper bounds matrices
  arma::mat phi_beta_bounds = arma::join_horiz(lower_beta, upper_beta);
  arma::mat phi_w_bounds = arma::join_horiz(lower_w, upper_w);
  
  // Gibbs sampling loop
  for (int i = 0; i < mcmc; i++) {
    // Update phi_w
    phi_w_cur = phi_RW(knots, w_knots_cur, sigmasq_w_cur, phi_w_cur, phi_w_proposal_sd, phi_w_bounds);
    phi_w_samples(i) = phi_w_cur;
    
    // Update phi_beta for each covariate
    for (int j = 0; j < p; j++) {
      phi_beta_cur(j) = phi_RW(knots, beta_knots_cur.col(j), sigmasq_beta_cur(j), phi_beta_cur(j), phi_beta_proposal_sd(j), phi_beta_bounds.row(j));
    }
    phi_beta_samples.row(i) = phi_beta_cur.t();
    
    // Update beta* at knot locations
    for (int j = 0; j < p; j++) {
      beta_knots_cur.col(j) = update_beta_r(tausq_cur, sigmasq_beta_cur(j), phi_beta_cur(j), knots, s);
    }
    
    // Update w* at knot locations
    w_knots_cur = update_w_s(Y, X, beta_knots_cur, sigmasq_w_cur, phi_w_cur, tausq_cur, knots);
    
    // Update sigma^2_beta for each covariate
    for (int j = 0; j < p; j++) {
      sigmasq_beta_cur(j) = update_sigma2_r(beta_knots_cur.col(j), a_beta(j), b_beta(j), phi_beta_cur(j), knots);
    }
    
    // Update sigma^2_w
    sigmasq_w_cur = update_sigma2_r(w_knots_cur, a_w, b_w, phi_w_cur, knots);
    
    // Calculate beta_tilde and w_tilde
    arma::vec beta_tilde = calc_beta_tilde(s, knots, phi_w_cur, beta_knots_cur);
    arma::vec w_tilde = calc_w_tilde(s, knots, phi_w_cur, w_knots_cur);
    
    // Store beta_tilde and w_tilde as the new beta and w
    beta_samples.slice(i) = beta_tilde;
    w_samples.row(i) = w_tilde.t();
    
    // Update sigmasq_beta using update_sigma2_r
    for (int j = 0; j < p; j++) {
      sigmasq_beta_cur(j) = update_sigma2_r(beta_knots_cur.col(j), a_beta(j), b_beta(j), phi_beta_cur(j), knots);
    }
    sigmasq_beta_samples.row(i) = sigmasq_beta_cur.t();
    
    // Update sigmasq_w using update_sigma2_r
    sigmasq_w_cur = update_sigma2_r(w_knots_cur, a_w, b_w, phi_w_cur, knots);
    sigmasq_w_samples(i) = sigmasq_w_cur;
    
    // Update tau^2 using update_tau2_r
    tausq_cur = update_tau2_r(Y, X, beta_tilde, w_tilde, a_t, b_t);
    tausq_samples(i) = tausq_cur;
  }
  
  // Return posterior samples
  Rcpp::List output;
  output["beta_samples"] = beta_samples;          // Real beta at observed locations
  output["w_samples"] = w_samples;                // Real w at observed locations
  output["phi_beta_samples"] = phi_beta_samples;  // phi_beta samples
  output["phi_w_samples"] = phi_w_samples;        // phi_w samples
  output["sigmasq_beta_samples"] = sigmasq_beta_samples; // sigma^2_beta samples
  output["sigmasq_w_samples"] = sigmasq_w_samples; // sigma^2_w samples
  output["tausq_samples"] = tausq_samples;         // tau^2 samples
  
  return output;
}  