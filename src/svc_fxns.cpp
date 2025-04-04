#include <RcppArmadillo.h>
#include "svc_fxns.h"

arma::mat calc_bigC(
    const arma::mat& const_bigC,
    double phi
) {
  return exp(const_bigC / phi);
}

arma::mat calc_lilc(
    const arma::mat& const_lilc,
    double phi
) {
  return exp(const_lilc / phi);
}

arma::mat get_R(
  const arma::mat& A
){
  arma::mat symA = 0.5 * (A + A.t()); // Ensure symmetry
  symA += arma::eye(A.n_rows, A.n_cols) * 1e-5; // Add a small value to the diagonal for numerical stability
  arma::mat R = arma::chol(symA);
  return R; // Return the Cholesky factor
}

arma::mat inv_Chol(
    const arma::mat& R
) {
  arma::mat Rinv = arma::inv(R);
  arma::mat Ainv = Rinv * Rinv.t();
  return Ainv;
}

double logdet(
    const arma::mat& R
) {
  double logdet = 2 * arma::accu(arma::log(arma::diagvec(R)));
  return logdet;
}

double logit(double x, double l, double u){
  return -log( (u-l)/(x-l) -1.0 );
}

double logistic(double x, double l, double u){
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

double normal_proposal_logitscale(const double& x, double l, double u){
  //return x - 2 * log(1 + exp(x));
  return -log(u-x) - log(x-l);
}

double calc_jacobian(double new_param, double param, const arma::mat& set_unif_bounds){
  // logit normal proposal
  double jac = normal_proposal_logitscale(param, set_unif_bounds(0, 0), set_unif_bounds(0, 1)) -
    normal_proposal_logitscale(new_param, set_unif_bounds(0, 0), set_unif_bounds(0, 1));
  
  return jac;
}

double GP_log_density(
    const arma::vec& x, // can be beta or w vector
    double sigmasq,
    const arma::mat& C_phi_inv,
    double C_phi_logdet
) {
  arma::mat Sigma_inv = C_phi_inv / sigmasq;
  double Sigma_logdet = C_phi_inv.n_cols * log(sigmasq) + C_phi_logdet;
  
  double log_likelihood =  -0.5 * Sigma_logdet -
    0.5 * arma::as_scalar(x.t() * Sigma_inv * x);
  
  return log_likelihood;
}

bool do_I_accept(double logaccept){
  double u = arma::randu();
  bool answer = exp(logaccept) > u;
  return answer;
}

phi_beta::phi_beta(
  unsigned int mcmc,
  unsigned int r_in,
  const arma::vec phi_start,
  const arma::vec proposal_sd_in,
  const arma::mat& phi_bounds_in,
  const arma::mat& const_bigC_in,
  const arma::mat& const_lilc_in
)
  : const_bigC(const_bigC_in),
    const_lilc(const_lilc_in)
{
  samples = arma::zeros(mcmc);
  acceptance = arma::zeros(mcmc);
  phi_cur = phi_start(r_in);
  C_phi_cur = calc_bigC(const_bigC_in, phi_cur);
  R_phi_cur = get_R(C_phi_cur);
  C_phi_cur_inv = inv_Chol(R_phi_cur);
  C_phi_cur_logdet = logdet(R_phi_cur);
  c_phi_cur = calc_lilc(const_lilc_in, phi_cur);
  proposal_sd = proposal_sd_in(r_in);
  phi_bounds = phi_bounds_in.row(r_in);
  iter = 0;
  r = r_in;
}

std::vector<phi_beta> initialize_phi_beta(
    unsigned int p,
    unsigned int mcmc,
    const arma::vec& phi_beta_start,
    const arma::vec& phi_beta_proposal_sd,
    const arma::mat& phi_beta_bounds,
    const arma::mat& const_bigC,
    const arma::mat& const_lilc
){
  std::vector<phi_beta> phi_beta_vec;
  
  for (unsigned int k = 0; k < p; k++) {
    phi_beta_vec.emplace_back(
      mcmc,
      k,
      phi_beta_start,
      phi_beta_proposal_sd,
      phi_beta_bounds,
      const_bigC,
      const_lilc
    );
  }
  return phi_beta_vec;
}

void phi_beta::RWupdate(
    const arma::mat& x_knots_mat, // can be beta or w vector
    const arma::vec& sigmasq_cur_vec
) {
  arma::vec x_knots = x_knots_mat.col(r);
  double sigmasq_cur = sigmasq_cur_vec(r);
  
  double phi_alt = par_huvtransf_back(par_huvtransf_fwd(
    phi_cur, phi_bounds) + proposal_sd * arma::randn(), phi_bounds);

  arma::mat C_phi_alt = calc_bigC(const_bigC, phi_alt);
  arma::mat R_phi_alt = get_R(C_phi_alt);
  arma::mat C_phi_alt_inv = inv_Chol(R_phi_alt);
  double C_phi_alt_logdet = logdet(R_phi_alt);

  // Calculate the log density of (w | sigma^2, phi)
  double logdens_cur = GP_log_density(x_knots, sigmasq_cur, C_phi_cur_inv, C_phi_cur_logdet);
  double logdens_alt = GP_log_density(x_knots, sigmasq_cur, C_phi_alt_inv, C_phi_alt_logdet);

  // Calculate the Jacobian of the proposal from transformation
  double jacobian  = calc_jacobian(phi_alt, phi_cur, phi_bounds);

  // Calculate the log acceptance probability
  double logaccept = logdens_alt - logdens_cur + jacobian;
  
  bool accepted = do_I_accept(logaccept);

  if(accepted){
    phi_cur = phi_alt;
    C_phi_cur = C_phi_alt;
    R_phi_cur = R_phi_alt;
    C_phi_cur_inv = C_phi_alt_inv;
    C_phi_cur_logdet = C_phi_alt_logdet;
    c_phi_cur = calc_lilc(const_lilc, phi_cur);
    acceptance(iter) = 1;
  }
  
  samples(iter) = phi_cur;
  iter++;
}

arma::vec update_beta_r_knots(
    const arma::vec& Y_knots,
    const arma::mat& X_knots,
    const arma::mat& X_knots_squared,
    const arma::mat& beta_knots,
    const arma::mat& C_phi_cur_inv,
    double sigma_r_square,
    double tau_square,
    int j
) {
  // Calculate the covariance matrix K based on the locations and phi_r
  arma::mat Sigma0inv = C_phi_cur_inv / sigma_r_square; // Inverse of the prior covariance matrix
  arma::mat Sigmainv = arma::diagmat(X_knots_squared.col(j)) / tau_square;
  arma::mat Sigma1 = inv_Chol(Sigma0inv + Sigmainv);
  
  arma::mat X_new = X_knots;
  X_new.shed_col(j); // Removes column j in-place
  
  arma::mat beta_new = beta_knots;
  beta_new.shed_col(j); // Removes column j in-place
  
  // Compute sum without column j
  arma::vec sumXbeta = arma::sum(X_new % beta_new, 1);
  arma::vec Y_tilde = Y_knots - sumXbeta;
  arma::vec result = Y_tilde / X_knots.col(j);
  arma::vec mu1 = Sigma1 * Sigmainv * result;
  
  arma::vec beta_r_knots = arma::mvnrnd(mu1, Sigma1);
  
  return beta_r_knots;
}

arma::vec calc_x_tilde(
    const arma::mat& lilc_cur,
    const arma::mat& bigC_inv_cur,
    const arma::vec& x_knots
) {
  
  return lilc_cur.t() * bigC_inv_cur * x_knots;
}

double update_sigma2_r(
    const arma::vec& beta_r,
    double a_r,
    double b_r,
    const arma::mat& C_inv
) {
  // Calculate (beta_r' C^{-1} beta_r)
  double beta_r_sum_square = arma::as_scalar(beta_r.t() * C_inv * beta_r);
  
  // Calculate the shape and scale parameters for the inverse gamma distribution
  double a_r_post = a_r + beta_r.n_elem / 2.0;
  double b_r_post = b_r + beta_r_sum_square / 2.0;
  
  //Sample sigma_r^2 from the inverse gamma distribution
  return 1.0 / arma::randg(arma::distr_param(a_r_post, 1.0 / b_r_post));
}

double update_tau2_r(const arma::vec& Y, const arma::mat& X, const arma::mat& beta, double a_t, double b_t) {
  arma::vec residual = Y - arma::sum(X % beta, 1);
  double residual_sum_square = arma::dot(residual, residual);
  double a_t_post = a_t + Y.n_elem / 2.0;
  double b_t_post = b_t + residual_sum_square / 2.0;
  
  // Sample tau^2 from the inverse gamma distribution
  return 1.0 / arma::randg(arma::distr_param(a_t_post, 1.0 / b_t_post));
}
