# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

GP_Gibbs <- function(Y, X, s, knots, beta_knots_start, w_knots_start, phi_beta_start, phi_w_start, sigmasq_beta_start, sigmasq_w_start, tausq_start, phi_beta_proposal_sd, phi_w_proposal_sd, a_beta, b_beta, a_w, b_w, a_t, b_t, lower_beta, upper_beta, lower_w, upper_w, mcmc = 1000L) {
    .Call(`_svc_GP_Gibbs`, Y, X, s, knots, beta_knots_start, w_knots_start, phi_beta_start, phi_w_start, sigmasq_beta_start, sigmasq_w_start, tausq_start, phi_beta_proposal_sd, phi_w_proposal_sd, a_beta, b_beta, a_w, b_w, a_t, b_t, lower_beta, upper_beta, lower_w, upper_w, mcmc)
}

GP_Gibbs_phi_test <- function(Y, X, s, knots, beta_knots_start, w_knots_start, phi_beta_start, phi_w_start, sigmasq_beta_start, sigmasq_w_start, tausq_start, phi_beta_proposal_sd, phi_w_proposal_sd, lower_beta, upper_beta, lower_w, upper_w, mcmc = 1000L) {
    .Call(`_svc_GP_Gibbs_phi_test`, Y, X, s, knots, beta_knots_start, w_knots_start, phi_beta_start, phi_w_start, sigmasq_beta_start, sigmasq_w_start, tausq_start, phi_beta_proposal_sd, phi_w_proposal_sd, lower_beta, upper_beta, lower_w, upper_w, mcmc)
}

