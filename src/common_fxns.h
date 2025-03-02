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