#ifndef SVC_FXNS_H
#define SVC_FXNS_H

#include <RcppArmadillo.h>

arma::mat calc_bigC(
    const arma::mat& const_bigC,
    double phi
);

arma::mat get_R(const arma::mat& A);

arma::mat inv_Chol(const arma::mat& R);

double logdet(const arma::mat& R);

double logit(double x, double l, double u);

double logistic(double x, double l, double u);

double par_huvtransf_fwd(double par, const arma::mat& set_unif_bounds);

double par_huvtransf_back(double par, const arma::mat& set_unif_bounds);

double normal_proposal_logitscale(const double& x, double l, double u);

double calc_jacobian(
    double new_param,
    double param,
    const arma::mat& set_unif_bounds
);

double GP_log_density(
    const arma::vec& x, // can be beta or w vector
    double sigmasq,
    const arma::mat& C_phi_inv,
    double C_phi_logdet
);

bool do_I_accept(double logaccept);

arma::mat calc_lilc(
    const arma::mat& const_lilc,
    double phi
);

class phi_beta {
public:
  arma::vec samples;
  arma::vec acceptance;
  double phi_cur;
  arma::mat C_phi_cur;
  arma::mat R_phi_cur;
  arma::mat C_phi_cur_inv;
  double C_phi_cur_logdet;
  arma::mat c_phi_cur;
  double proposal_sd;
  arma::mat phi_bounds;
  const arma::mat& const_bigC;
  const arma::mat& const_lilc;
  unsigned int iter;
  unsigned int r;
  
  phi_beta(
    unsigned int mcmc,
    unsigned int r_in,
    const arma::vec phi_start,
    const arma::vec proposal_sd_in,
    const arma::mat& phi_bounds_in,
    const arma::mat& const_bigC_in,
    const arma::mat& const_lilc_in
  );
  
  void RWupdate(
      const arma::mat& x_knots_mat, // can be beta or w vector
      const arma::vec& sigmasq_cur_vec
  );
};

std::vector<phi_beta> initialize_phi_beta(
    unsigned int p,
    unsigned int mcmc,
    const arma::vec& phi_start,
    const arma::vec& phi_proposal_sd,
    const arma::mat& phi_bounds,
    const arma::mat& const_bigC,
    const arma::mat& const_lilc
);

arma::vec update_beta_r_knots(
    const arma::vec& Y_knots,
    const arma::mat& X_knots,
    const arma::mat& X_knots_squared,
    const arma::mat& beta_knots,
    const arma::mat& C_phi_cur_inv,
    double sigma_r_square,
    double tau_square,
    int j
);

arma::vec calc_x_tilde(
    const arma::mat& lilc_cur,
    const arma::mat& bigC_inv_cur,
    const arma::vec& x_knots
);

double update_sigma2_r(
    const arma::vec& beta_r,
    double a_r,
    double b_r,
    const arma::mat& C_inv
);

double update_tau2_r(
    const arma::vec& Y,
    const arma::mat& X,
    const arma::mat& beta,
    double a_t,
    double b_t
);

#endif