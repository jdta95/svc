#include <RcppArmadillo.h>
#include "common_fxns.h"

// arma::mat calc_C_tilde(
//     arma::mat s,
//     arma::mat knots,
//     double phi
// ) {
//   arma::mat C_star = calc_C(knots, phi);
//   arma::mat c = calc_c(s, knots, phi);
//   arma::mat C_tilde = c.t() * inv_Chol(C_star) * c;
//   
//   return C_tilde;
// }

// [[Rcpp::export]]
double phi_RW(
  arma::mat knots,
  arma::vec w_knots, // Use w_knots
  double sigmasq_cur,
  double phi_cur,
  double proposal_sd,
  arma::mat phi_bounds
) {
  
  // int accept_count = 0;
  
  // print a message here
  Rcpp::Rcout << "Check 1 good. " << std::endl;
  
  double phi_alt = par_huvtransf_back(par_huvtransf_fwd(
    phi_cur, phi_bounds) + proposal_sd * arma::randn(), phi_bounds);
  
  Rcpp::Rcout << "Check 2 good. " << std::endl;
  
  arma::mat C_star_cur = calc_C(knots, phi_cur);
  arma::mat C_star_alt = calc_C(knots, phi_alt);
  
  Rcpp::Rcout << "Check 3 good. " << std::endl;
  
  double curr_logdens = wGP_log_density(w_knots, sigmasq_cur, C_star_cur);
  double prop_logdens = wGP_log_density(w_knots, sigmasq_cur, C_star_alt);
  
  Rcpp::Rcout << "Check 4 good. " << std::endl;
  
  // make move
  double jacobian  = calc_jacobian(phi_alt, phi_cur, phi_bounds);
  
  Rcpp::Rcout << "Check 5 good. " << std::endl;
  
  double logaccept = prop_logdens - curr_logdens + jacobian;
  
  Rcpp::Rcout << "Check 6 good. " << std::endl;
  
  bool accepted = do_I_accept(logaccept);
  
  Rcpp::Rcout << "Check 7 good. " << std::endl;
  
  if(accepted){
    phi_cur = phi_alt;
    // accept_count++;
  }
  
  Rcpp::Rcout << "Check 8 good. " << std::endl;
  
  return phi_cur;
}
