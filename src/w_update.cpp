//#include <RcppArmadillo.h>
//#include "common_fxns.h"

//// [[Rcpp::depends(RcppArmadillo)]]

// Function to calculate w_tilde
//// [[Rcpp::export]]
//arma::vec calc_w_tilde(const arma::mat& s, const arma::mat& knots, double phi_w, const arma::vec& w_star) {
  //  arma::mat c_transpose = calc_c(s, knots, phi_w).t();
//    arma::mat C_inv = inv_Chol(calc_C(knots, phi_w));
//    return c_transpose * C_inv * w_star;
//}