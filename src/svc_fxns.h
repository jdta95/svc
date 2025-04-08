#ifndef SVC_FXNS_H
#define SVC_FXNS_H

#include <RcppArmadillo.h>

arma::mat calc_bigC(
    const arma::mat& const_bigC,
    double phi
);

arma::mat get_R(const arma::mat& A);

arma::mat inv_Chol(const arma::mat& A);

arma::mat inv_Chol_R(const arma::mat& R);

double logdet_R(const arma::mat& R);

arma::mat calc_K(
    const arma::vec& sigmasq,
    const arma::vec& phi,
    const arma::mat& const_K,
    arma::uword n,
    arma::uword p
);

double calc_logdens(
    const arma::mat& K,
    double tausq,
    const arma::mat& Z,
    const arma::mat& YX,
    arma::uword n,
    arma::uword p
);

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
    const arma::vec& x,
    double sigmasq,
    const arma::mat& C_phi_inv,
    double C_phi_logdet
);

bool do_I_accept(double logaccept);

arma::mat calc_lilc(
    const arma::mat& const_lilc,
    double phi
);

class phi {
public:
  arma::vec samples;
  arma::vec acceptance;
  double phi_cur;
  arma::mat C_phi_cur;
  arma::mat R_phi_cur;
  arma::mat C_phi_cur_inv;
  double C_phi_cur_logdet;
  arma::mat c_phi_cur;
  double paramsd;
  arma::mat phi_bounds;
  const arma::mat& const_bigC;
  const arma::mat& const_lilc;
  unsigned int accept_count;
  int iter;
  double accept_ratio;
  double alpha_star;
  double gamma;
  int g0;
  int i;
  double S;
  double Sigma;
  double prodparam;
  double eta;
  bool started;
  unsigned int r;
  
  phi(
    unsigned int mcmc,
    unsigned int r_in,
    const arma::vec& phi_start,
    const arma::vec& proposal_sd_in,
    const arma::mat& phi_bounds_in,
    const arma::mat& const_bigC_in,
    const arma::mat& const_lilc_in,
    double target_accept
  );
  
  void RWupdate(
      const arma::mat& x_knots_mat,
      const arma::vec& sigmasq_cur_vec
  );
  
  void adapt(
      double u,
      double alpha
  );
  
};

std::vector<phi> initialize_phi(
    unsigned int p,
    unsigned int mcmc,
    const arma::vec& phi_start,
    const arma::vec& phi_proposal_sd,
    const arma::mat& phi_bounds,
    const arma::mat& const_bigC,
    const arma::mat& const_lilc,
    double target_accept
);

arma::vec update_w_r_knots(
    const arma::vec& Y_knots,
    const arma::mat& X_knots,
    const arma::mat& X_knots_squared,
    const arma::mat& w_knots,
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
    const arma::vec& w_r,
    double a_r,
    double b_r,
    const arma::mat& C_inv
);

double update_tau2_r(
    const arma::vec& Y,
    const arma::mat& X,
    const arma::mat& w,
    double a_t,
    double b_t
);

#endif