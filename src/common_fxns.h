#ifndef COMMON_FXNS_H
#define COMMON_FXNS_H

#include <RcppArmadillo.h>

arma::mat inv_Chol(
    arma::mat A
) {
  A = 0.5 * (A + A.t()); // Ensure symmetry
  A += arma::eye(A.n_rows, A.n_cols) * 1e-10; // Add a small value to the diagonal for numerical stability
  arma::mat R = arma::chol(A);
  arma::mat Rinv = arma::inv(R);
  arma::mat Ainv = Rinv * Rinv.t();
  return Ainv;
}

double logdet(
    arma::mat A
) {
  A = 0.5 * (A + A.t()); // Ensure symmetry
  A += arma::eye(A.n_rows, A.n_cols) * 1e-10; // Add a small value to the diagonal for numerical stability
  arma::mat R = arma::chol(A);
  double logdet = 2 * arma::accu(arma::log(arma::diagvec(R)));
  return logdet;
}

arma::mat calc_C(
    arma::mat s,
    double phi
) {
  int n = s.n_rows;
  arma::mat C(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      // calculate the distance
      C(i, j) = std::sqrt((s(i, 0) - s(j, 0)) * (s(i, 0) - s(j, 0)) +
        (s(i, 1) - s(j, 1)) * (s(i, 1) - s(j, 1)));
    }
  }
  C = exp(-0.5 * C / phi);
  return C;
}

arma::mat calc_c(
    arma::mat s,
    arma::mat knots,
    double phi
) {
  int n = s.n_rows;
  int m = knots.n_rows;
  arma::mat c(m, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      // calculate the distance between each s and each knot
      c(j, i) = std::sqrt((s(i, 0) - knots(j, 0)) * (s(i, 0) - knots(j, 0)) +
        (s(i, 1) - knots(j, 1)) * (s(i, 1) - knots(j, 1)));
    }
  }
  c = exp(-0.5 * c / phi);
  return c;
}

arma::mat calc_C_tilde(
    arma::mat s,
    arma::mat knots,
    double phi
) {
  arma::mat C_star = calc_C(knots, phi);
  arma::mat c = calc_c(s, knots, phi);
  arma::mat C_tilde = c.t() * inv_Chol(C_star) * c;
  
  return C_tilde;
}

double logit(double x, double l=0, double u=1){
  return -log( (u-l)/(x-l) -1.0 );
}

double logistic(double x, double l=0, double u=1){
  return l + (u-l)/(1.0+exp(-x));
}

double par_huvtransf_fwd(double par, const arma::mat& set_unif_bounds){
  par = logit(par, set_unif_bounds(0, 0), set_unif_bounds(0, 1));
  return par;
}

double par_huvtransf_back(double par, const arma::mat& set_unif_bounds){
  par = logistic(par, set_unif_bounds(0, 0), set_unif_bounds(0, 1));
  return par;
}

double normal_proposal_logitscale(const double& x, double l=0, double u=1){
  //return x - 2 * log(1 + exp(x));
  return -log(u-x) - log(x-l);
}

double calc_jacobian(double new_param, double param, const arma::mat& set_unif_bounds){
  // logit normal proposal
  double jac = normal_proposal_logitscale(param, set_unif_bounds(0, 0), set_unif_bounds(0, 1)) -
          normal_proposal_logitscale(new_param, set_unif_bounds(0, 0), set_unif_bounds(0, 1));
  
  return jac;
}

double wGP_log_density(
    arma::vec w,
    double sigmasq,
    arma::mat C_phi
) {
  arma::mat Sigma = sigmasq * C_phi;
  
  double log_likelihood =  -0.5 * logdet(Sigma) -
    0.5 * arma::as_scalar(w.t() * inv_Chol(Sigma) * w);
  
  return log_likelihood;
}

bool do_I_accept(double logaccept){
  double u = arma::randu();
  bool answer = exp(logaccept) > u;
  return answer;
}

// Update functions

double phi_RW(
    const arma::mat& knots,
    const arma::vec& w_knots, // Use w_knots
    double sigmasq_cur,
    double phi_cur,
    double proposal_sd,
    const arma::mat& phi_bounds
) {
  
  double phi_alt = par_huvtransf_back(par_huvtransf_fwd(
    phi_cur, phi_bounds) + proposal_sd * arma::randn(), phi_bounds);
  
  arma::mat C_star_cur = calc_C(knots, phi_cur);
  arma::mat C_star_alt = calc_C(knots, phi_alt);
  
  // Calculate the log density of (w | sigma^2, phi)
  double curr_logdens = wGP_log_density(w_knots, sigmasq_cur, C_star_cur);
  double prop_logdens = wGP_log_density(w_knots, sigmasq_cur, C_star_alt);
  
  // Calculate the Jacobian of the proposal from transformation
  double jacobian  = calc_jacobian(phi_alt, phi_cur, phi_bounds);
  
  // Calculate the log acceptance probability
  double logaccept = prop_logdens - curr_logdens + jacobian;
  
  bool accepted = do_I_accept(logaccept);
  
  if(accepted){
    phi_cur = phi_alt;
  }
  
  return phi_cur;
}

arma::vec update_beta_r(
    const arma::vec& Y_star,
    const arma::mat& X_star,
    const arma::mat& beta_star,
    const arma::mat& w_star,
    const arma::mat& Knots,
    double sigma_r_square,
    double phi_r,
    double tau_square
) {
  int m = Knots.n_rows;
  int p = X_star.n_cols;
  
  // Calculate the covariance matrix K based on the locations and phi_r
  arma::mat Sigma0 = sigma_r_square * calc_C(Knots, phi_r);
  arma::mat Sigma0inv = inv_Chol(Sigma0);
  
  arma::mat Sigma = tau_square * arma::eye(m, m);
  arma::mat Sigmainv = inv_Chol(Sigma);
  
  arma::mat Sigma1 = inv_Chol(Sigma0inv + Sigmainv);
  
  arma::vec sumXbeta = arma::sum(X_star % beta_star, 1);
  
  arma::vec mu1 = Sigma1 * (Sigmainv * (Y_star - sumXbeta - w_star));
  
  arma::vec beta_r_star = arma::mvnrnd(mu1, Sigma1);
}


double update_sigma2_r(const arma::vec& beta_r, double a_r, double b_r, double phi_r, const arma::mat& Knots) {
  // Calculate the covariance matrix K based on the locations and phi_r
  arma::mat C = calc_C(Knots, phi_r);
  
  // Calculate the inverse of the covariance matrix using Cholesky decomposition
  arma::mat C_inv = inv_Chol(C);
  
  // Calculate (beta_r' C^{-1} beta_r)
  double beta_r_sum_square = arma::as_scalar(beta_r.t() * C_inv * beta_r);
  
  // Calculate the shape and scale parameters for the inverse gamma distribution
  double a_r_post = a_r + beta_r.n_elem / 2.0;
  double b_r_post = b_r + beta_r_sum_square / 2.0;
  
  // Sample sigma_r^2 from the inverse gamma distribution
  return 1.0 / arma::randg(arma::distr_param(a_r_post, 1.0 / b_r_post));
}


double update_tau2_r(const arma::vec& Y, const arma::mat& X, const arma::vec& beta, const arma::vec& w, double a_t, double b_t) {
  arma::vec residual = Y - X * beta - w;
  double residual_sum_square = arma::dot(residual, residual);
  double a_t_post = a_t + Y.n_elem / 2.0;
  double b_t_post = b_t + residual_sum_square / 2.0;
  
  // Sample tau^2 from the inverse gamma distribution
  return 1.0 / arma::randg(arma::distr_param(a_t_post, 1.0 / b_t_post));
}


arma::vec update_w_s(const arma::vec& Y, const arma::mat& X, const arma::mat& beta_knots, double sigmasq_w, double phi_w, double tausq, const arma::mat& knots) {
  // Calculate the covariance matrix K based on the locations and phi_w
  arma::mat C = calc_C(knots, phi_w);
  
  // Calculate the inverse of the covariance matrix using Cholesky decomposition
  arma::mat C_inv = inv_Chol(C);
  
  // Calculate the variance parameter for the posterior distribution
  arma::mat post_var = inv_Chol(inv_Chol(tausq * arma::eye(C.n_rows, C.n_rows)) + inv_Chol(sigmasq_w * C));
  
  // Calculate the mean parameter for the posterior distribution
  arma::vec post_mean = post_var * (C_inv * (Y - X * beta_knots));
  
  // Sample w_s from the posterior distribution
  return post_mean + arma::chol(post_var) * arma::randn(C.n_rows);
}


arma::vec calc_w_tilde(const arma::mat& s, const arma::mat& knots, double phi_w, const arma::vec& w_star) {
  arma::mat c_transpose = calc_c(s, knots, phi_w).t();
  arma::mat C_inv = inv_Chol(calc_C(knots, phi_w));
  return c_transpose * C_inv * w_star;
}

arma::vec calc_beta_tilde(const arma::mat& s, const arma::mat& knots, double phi_beta, const arma::vec& beta_star) {
  arma::mat c_transpose = calc_c(s, knots, phi_beta).t();
  arma::mat C_inv = inv_Chol(calc_C(knots, phi_beta));
  return c_transpose * C_inv * beta_star;
}


#endif 
