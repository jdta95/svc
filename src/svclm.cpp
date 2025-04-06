// #include <RcppArmadillo.h>
// #include "svc_fxns.h"
// 
// // [[Rcpp::export]]
// Rcpp::List svclm(
//     const arma::vec& Y,               // Response vector
//     const arma::mat& X,               // Design matrix
//     const arma::mat& s,               // Spatial locations (n x 2 matrix)
//     const arma::vec& Y_knots,         // Response vector at knots
//     const arma::mat& X_knots,         // Design matrix at knots
//     const arma::mat& knots,           // Knot locations (m x 2 matrix)
//     const arma::mat& beta_knots_start,       // Initial value for coefficients: beta 1, ..., beta p
//     const arma::vec& phi_beta_start,        // Initial value for phi: beta 1, ..., beta p
//     const arma::vec& sigmasq_beta_start,     // Initial value for sigma^2: beta 1, ..., beta p
//     double tausq_start,               // Initial value for tau^2
//     const arma::vec& phi_beta_proposal_sd, // Proposal standard deviations for phi: beta 1, ..., beta p
//     const arma::vec& phi_beta_lower,      // Lower bounds for phi: beta 1, ..., beta p
//     const arma::vec& phi_beta_upper,      // Upper bounds for phi: beta 1, ..., beta p
//     const arma::vec& a_beta,            // Hyperparameter a for sigma^2_beta
//     const arma::vec& b_beta,            // Hyperparameter b for sigma^2_beta
//     double a_t,               // Hyperparameter a for tau^2
//     double b_t,               // Hyperparameter b for tau^2
//     unsigned int mcmc = 1000                   // Number of MCMC iterations
// ){
//   // Check input dimensions
//   arma::uword n = Y.n_elem; // number of observations
//   arma::uword p = X.n_cols; // number of predictors
//   arma::uword m = knots.n_rows; // number of knots
//   
//   // input errors
//   if (X_knots.n_rows != m || X_knots.n_cols != p) {
//     Rcpp::stop("Design matrix at knots 'X_knots' must be an m x p matrix.");
//   }
//   if (Y_knots.n_elem != m) {
//     Rcpp::stop("Response vector at knots 'Y_knots' must have length m.");
//   }
//   if (s.n_rows != n || s.n_cols != 2) {
//     Rcpp::stop("Spatial locations 's' must be an n x 2 matrix.");
//   }
//   if (knots.n_cols != 2) {
//     Rcpp::stop("Knot locations 'knots' must be an m x 2 matrix.");
//   }
//   if (beta_knots_start.n_cols != p) {
//     Rcpp::stop("Initial beta coefficients 'beta_knots_start' must have length p.");
//   }
//   if (phi_beta_start.n_elem != p) {
//     Rcpp::stop("Initial phi beta vector 'phi_beta_start' must have length p.");
//   }
//   if (sigmasq_beta_start.n_elem != p) {
//     Rcpp::stop("Initial sigmasq beta coefficients 'sigmasq_beta_start' must have length p.");
//   }
//   
//   // Rcpp::Rcout << "Check 10. " << std::endl;
//   
//   // initialize outputs
//   arma::cube beta_samples = arma::zeros(mcmc, n, p);  // posterior beta at observed locations
//   arma::mat phi_beta_samples = arma::zeros(mcmc, p); // mcmc x p matrix of samples for phi_beta's
//   arma::mat phi_beta_acceptance = arma::zeros(mcmc, p); // mcmc x p matrix to track acceptance for phi_beta
//   arma::mat sigmasq_beta_samples = arma::zeros(mcmc, p); // mcmc x p matrix of samples for sigmasq_beta
//   arma::vec tausq_samples = arma::zeros(mcmc); // mcmc vector of samples for tau^2
//   
//   // initialize parameters
//   arma::mat beta_knots_cur = beta_knots_start; // current beta coefficients
//   arma::vec sigmasq_beta_cur = sigmasq_beta_start;    // Current sigma^2_beta
//   double tausq_cur = tausq_start; // Current tau^2
// 
//   // Create lower and upper bounds matrices
//   arma::mat phi_beta_bounds = arma::join_horiz(phi_beta_lower, phi_beta_upper);
//   
//   // calculate constant parts of various calculations
//   // constant part of bigC
//   arma::mat const_bigC(m, m);
//   for (unsigned int i = 0; i < m; i++) {
//     for (unsigned int j = 0; j < m; j++) {
//       // calculate the distance^2
//       const_bigC(i, j) = (knots(i, 0) - knots(j, 0)) * (knots(i, 0) - knots(j, 0)) +
//         (knots(i, 1) - knots(j, 1)) * (knots(i, 1) - knots(j, 1));
//     }
//   }
//   const_bigC *= -0.5;
//   
//   // constant part of lilc
//   arma::mat const_lilc(m, n);
//   for (unsigned int i = 0; i < m; i++) {
//     for (unsigned int j = 0; j < n; j++) {
//       // calculate the distance between each s and each knot
//       const_lilc(i, j) = (s(j, 0) - knots(i, 0)) * (s(j, 0) - knots(i, 0)) +
//         (s(j, 1) - knots(i, 1)) * (s(j, 1) - knots(i, 1));
//     }
//   }
//   const_lilc *= -0.5;
//   
//   // constant X_knots_squared
//   arma::mat X_knots_squared = arma::pow(X_knots, 2);
//   
//   // // Initialize phi objects
//   std::vector<phi_beta> phi_beta_vec = initialize_phi_beta(
//     p,
//     mcmc,
//     phi_beta_start,
//     phi_beta_proposal_sd,
//     phi_beta_bounds,
//     const_bigC,
//     const_lilc
//   );
//   
//   // MCMC loop
//   for (unsigned int i = 0; i < mcmc; i++) {
//     
//     // Rcpp::Rcout << "Iteration: " << i + 1 << std::endl;
//     
//     for (unsigned int j = 0; j < p; j++) {
//       // Update phi_j via Random Walk Metropolis
//       phi_beta_vec.at(j).RWupdate(beta_knots_cur, sigmasq_beta_cur);
//       
//       // Update sigma^2_beta
//       sigmasq_beta_cur(j) = update_sigma2_r(beta_knots_cur.col(j), a_beta(j), b_beta(j), phi_beta_vec.at(j).C_phi_cur_inv);
// 
//       // update beta_knots_j
//       beta_knots_cur.col(j) = update_beta_r_knots(
//         Y_knots,
//         X_knots,
//         X_knots_squared,
//         beta_knots_cur,
//         phi_beta_vec.at(j).C_phi_cur_inv,
//         sigmasq_beta_cur(j),
//         tausq_cur,
//         j
//       );
//       
//       // predict beta_j at all locations
//       beta_samples.subcube(i, 0, j, i, n - 1, j) = calc_x_tilde(
//         phi_beta_vec.at(j).c_phi_cur,
//         phi_beta_vec.at(j).C_phi_cur_inv,
//         beta_knots_cur.col(j)
//       );
//     }
//     
//     sigmasq_beta_samples.row(i) = sigmasq_beta_cur.t(); // store the current sigma^2_beta samples
//     
//     // update tau^2
//     // Currently using only knots
//     // Should use all locations if we can get beta predictions
//     tausq_cur = update_tau2_r(Y_knots, X_knots, beta_knots_cur, a_t, b_t);
//     tausq_samples(i) = tausq_cur;
//   }
//   
//   for (unsigned int j = 0; j < p; j++) {
//     // Store the final phi_beta samples and acceptance rates
//     phi_beta_samples.col(j) = phi_beta_vec.at(j).samples;
//     phi_beta_acceptance.col(j) = phi_beta_vec.at(j).acceptance;
//   }
//   
//   Rcpp::List output;
//   output["beta_samples"] = beta_samples; // posterior beta samples at observed locations
//   output["phi_beta_samples"] = phi_beta_samples;
//   output["phi_beta_acceptance"] = phi_beta_acceptance;
//   output["sigmasq_beta_samples"] = sigmasq_beta_samples; // posterior sigma^2_beta samples
//   output["tausq_samples"] = tausq_samples; // posterior tau^2 samples
//   
//   return output;
// }
