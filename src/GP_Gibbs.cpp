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
  arma::cube beta_samples = arma::zeros(n, p, mcmc);  // Real beta at observed locations
  arma::mat w_samples = arma::zeros(mcmc, n);         // Real w at observed locations
  // arma::mat beta_knots_samples = arma::zeros(mcmc, m, p); // beta* at knot locations
  // arma::mat w_knots_samples = arma::zeros(mcmc, m);       // w* at knot locations
  phi_beta phi_beta(mcmc, p);  // phi_beta samples
  phi_w phi_w(mcmc);          // phi_w samples
  arma::mat sigmasq_beta_samples = arma::zeros(mcmc, p); // sigma^2_beta samples
  arma::vec sigmasq_w_samples = arma::zeros(mcmc);    // sigma^2_w samples
  arma::vec tausq_samples = arma::zeros(mcmc);        // tau^2 samples
  
  // Initialize current values
  arma::mat beta_knots_cur = beta_knots_start;        // Current beta* at knots
  arma::vec w_knots_cur = w_knots_start;              // Current w* at knots
  phi_beta.samples.row(0) = phi_beta_start.t();            // Current phi_beta
  phi_w.samples(0) = phi_w_start;                     // Current phi_w
  arma::vec sigmasq_beta_cur = sigmasq_beta_start;    // Current sigma^2_beta
  double sigmasq_w_cur = sigmasq_w_start;             // Current sigma^2_w
  double tausq_cur = tausq_start;                     // Current tau^2
  arma::mat X_knots = X.rows(arma::linspace<arma::uvec>(0, n - 1, m));  
  arma::vec Y_knots = Y.rows(arma::linspace<arma::uvec>(0, n - 1, m));
  
  // Create lower and upper bounds matrices
  arma::mat phi_beta_bounds = arma::join_horiz(lower_beta, upper_beta);
  arma::mat phi_w_bounds = arma::join_horiz(lower_w, upper_w);
  
  // Gibbs sampling loop
  for (int i = 0; i < mcmc; i++) {
    // Update phi_w
    phi_w.RWupdate(i + 1, knots, w_knots_cur, sigmasq_w_cur, phi_w.samples(i), phi_w_proposal_sd, phi_w_bounds);
    
    // Update phi_beta for each covariate
    for (int j = 0; j < p; j++) {
      phi_beta.RWupdate(i + 1, j, knots, beta_knots_cur.col(j), sigmasq_beta_cur(j), phi_beta.samples(i, j), phi_beta_proposal_sd(j), phi_beta_bounds.row(j));
    }
    
    // std::cout << "Check 10. " << std::endl;
    
    // Update beta* at knot locations
    for (int j = 0; j < p; j++) {
      beta_knots_cur.col(j) = update_beta_r(Y_knots, X_knots, beta_knots_cur, w_knots_cur, knots, sigmasq_beta_cur(j), phi_beta.samples(i + 1, j), tausq_cur, j);
    }
    
    // Update w* at knot locations
    w_knots_cur = update_w_s(Y_knots, X_knots, beta_knots_cur, sigmasq_w_cur, phi_w.samples(i + 1), tausq_cur, knots);
    
    // Update sigma^2_beta for each covariate
    for (int j = 0; j < p; j++) {
      sigmasq_beta_cur(j) = update_sigma2_r(beta_knots_cur.col(j), a_beta(j), b_beta(j), phi_beta.samples(i + 1, j), knots);
    }
    
    // std::cout << "Check 20. " << std::endl;
    
    // Update sigma^2_w
    sigmasq_w_cur = update_sigma2_r(w_knots_cur, a_w, b_w, phi_w.samples(i + 1), knots);
    
    // std::cout << "Check 30. " << std::endl;
    
    // Calculate beta_tilde and w_tilde
    arma::mat beta_tilde = arma::zeros(n, p); // Initialize beta_tilde
    
    for (int j = 0; j < p; j++) {
      // Calculate beta_tilde for each covariate
      beta_tilde.col(j) = calc_beta_tilde(s, knots, phi_beta.samples(i + 1, j), beta_knots_cur.col(j));
    }
    
    // arma::vec beta_tilde = calc_beta_tilde(s, knots, phi_w_cur, beta_knots_cur);
    
    // std::cout << "Check 35. " << std::endl;
    
    arma::vec w_tilde = calc_w_tilde(s, knots, phi_w.samples(i + 1), w_knots_cur);
    
    // std::cout << "Check 40. " << std::endl;
    
    // Store beta_tilde and w_tilde as the new beta and w
    beta_samples.slice(i) = beta_tilde;
    
    // std::cout << "Check 45. " << std::endl;
    
    w_samples.row(i) = w_tilde.t();
    
    // std::cout << "Check 50. " << std::endl;
    
    // Update sigmasq_beta using update_sigma2_r
    for (int j = 0; j < p; j++) {
      sigmasq_beta_cur(j) = update_sigma2_r(beta_knots_cur.col(j), a_beta(j), b_beta(j), phi_beta.samples(i + 1, j), knots);
    }
    sigmasq_beta_samples.row(i) = sigmasq_beta_cur.t();
    
    // std::cout << "Check 60. " << std::endl;
    
    // Update sigmasq_w using update_sigma2_r
    sigmasq_w_cur = update_sigma2_r(w_knots_cur, a_w, b_w, phi_w.samples(i + 1), knots);
    sigmasq_w_samples(i) = sigmasq_w_cur;
    
    // std::cout << "Check 70. " << std::endl;
    
    // Update tau^2 using update_tau2_r
    tausq_cur = update_tau2_r(Y, X, beta_tilde, w_tilde, a_t, b_t);
    tausq_samples(i) = tausq_cur;
    
    // std::cout << "Check 80. " << std::endl;
  }
  
  phi_beta.samples.shed_row(0);
  phi_beta.acceptance.shed_row(0);
  phi_w.samples.shed_row(0);
  phi_w.acceptance.shed_row(0);
  
  // Return posterior samples
  Rcpp::List output;
  output["beta_samples"] = beta_samples;          // Real beta at observed locations
  output["w_samples"] = w_samples;                // Real w at observed locations
  output["phi_beta_samples"] = phi_beta.samples;  // phi_beta samples
  output["phi_beta_acceptance"] = phi_beta.acceptance; // phi_beta acceptance matrix
  output["phi_w_samples"] = phi_w.samples;        // phi_w samples
  output["phi_w_acceptance"] = phi_w.acceptance;       // phi_w acceptance matrix
  output["sigmasq_beta_samples"] = sigmasq_beta_samples; // sigma^2_beta samples
  output["sigmasq_w_samples"] = sigmasq_w_samples; // sigma^2_w samples
  output["tausq_samples"] = tausq_samples;         // tau^2 samples
  
  return output;
}   
