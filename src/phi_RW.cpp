// #include <RcppArmadillo.h>
// #include "common_fxns.h"
// 
// // [[Rcpp::export]]
// double phi_RW(
//     arma::mat knots,
//     arma::vec w_knots, // Use w_knots
//     double sigmasq_cur,
//     double phi_cur,
//     double proposal_sd,
//     arma::mat phi_bounds
// ) {
// 
//   double phi_alt = par_huvtransf_back(par_huvtransf_fwd(
//     phi_cur, phi_bounds) + proposal_sd * arma::randn(), phi_bounds);
// 
//   arma::mat C_star_cur = calc_C(knots, phi_cur);
//   arma::mat C_star_alt = calc_C(knots, phi_alt);
//   
//   // Calculate the log density of (w | sigma^2, phi)
//   double curr_logdens = wGP_log_density(w_knots, sigmasq_cur, C_star_cur);
//   double prop_logdens = wGP_log_density(w_knots, sigmasq_cur, C_star_alt);
// 
//   // Calculate the Jacobian of the proposal from transformation
//   double jacobian  = calc_jacobian(phi_alt, phi_cur, phi_bounds);
// 
//   // Calculate the log acceptance probability
//   double logaccept = prop_logdens - curr_logdens + jacobian;
// 
//   bool accepted = do_I_accept(logaccept);
// 
//   if(accepted){
//     phi_cur = phi_alt;
//   }
// 
//   return phi_cur;
// }
