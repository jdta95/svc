#ifndef COMMON_FXNS_H
#define COMMON_FXNS_H

#include <RcppArmadillo.h>

arma::mat inv_Chol(
    arma::mat A
);

double logdet(
    arma::mat A
);

arma::mat calc_C(
    arma::mat s,
    double phi
);

arma::mat calc_c(
    arma::mat s,
    arma::mat knots,
    double phi
);

arma::mat calc_C_tilde(
    arma::mat s,
    arma::mat knots,
    double phi
);

double logit(double x, double l=0, double u=1);

double logistic(double x, double l=0, double u=1);

double par_huvtransf_fwd(double par, const arma::mat& set_unif_bounds);

double par_huvtransf_back(double par, const arma::mat& set_unif_bounds);

double normal_proposal_logitscale(const double& x, double l=0, double u=1);

double calc_jacobian(double new_param, double param, const arma::mat& set_unif_bounds);

double wGP_log_density(
    arma::vec w,
    double sigmasq,
    arma::mat C_phi
);

bool do_I_accept(double logaccept);

// Update functions

double phi_RW(
    const arma::mat& knots,
    const arma::vec& w_knots, // Use w_knots
    double sigmasq_cur,
    double phi_cur,
    double proposal_sd,
    const arma::mat& phi_bounds
);

arma::vec update_beta_r(
    const arma::vec& Y_star,
    const arma::mat& X_star,
    const arma::mat& beta_star,
    const arma::mat& w_star,
    const arma::mat& Knots,
    double sigma_r_square,
    double phi_r,
    double tau_square
);

double update_sigma2_r(
    const arma::vec& beta_r,
    double a_r,
    double b_r,
    double phi_r,
    const arma::mat& Knots
);

double update_tau2_r(const arma::vec& Y, const arma::mat& X, const arma::mat& beta, const arma::vec& w, double a_t, double b_t);

arma::vec update_w_s(const arma::vec& Y, const arma::mat& X, const arma::mat& beta_knots, double sigmasq_w, double phi_w, double tausq, const arma::mat& knots);

arma::vec calc_w_tilde(const arma::mat& s, const arma::mat& knots, double phi_w, const arma::vec& w_star);

arma::mat calc_beta_tilde(const arma::mat& s, const arma::mat& knots, double phi_beta, const arma::mat& beta_star);

#endif 
