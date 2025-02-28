// #include <RcppArmadillo.h>
// #include "ramadapt.h"
// 
// using namespace std;
// 
// arma::mat inv_Chol(
//   arma::mat A
// ) {
//   A = 0.5 * (A + A.t()); // Ensure symmetry
//   A += arma::eye(A.n_rows, A.n_cols) * 1e-6; // Add a small value to the diagonal for numerical stability
//   arma::mat R = arma::chol(A);
//   arma::mat Rinv = arma::inv(R);
//   arma::mat Ainv = Rinv * Rinv.t();
//   return Ainv;
// }
// 
// double logdet(
//   arma::mat A
// ) {
//   A = 0.5 * (A + A.t()); // Ensure symmetry
//   A += arma::eye(A.n_rows, A.n_cols) * 1e-6; // Add a small value to the diagonal for numerical stability
//   arma::mat R = arma::chol(A);
//   double logdet = 2 * arma::accu(arma::log(arma::diagvec(R)));
//   return logdet;
// }
// 
// arma::mat calc_C_phi(
//     arma::mat X,
//     arma::vec phi,
//     int n
// ) {
//   arma::mat C_phi = arma::zeros(n, n);
//   
//   for(int i = 0; i < n; i++) {
//     for (int j = 0; j < n; j++) {
//       C_phi(i, j) = arma::as_scalar(
//         exp(
//           -0.5 * (X.row(i) - X.row(j)) * arma::diagmat(1 / phi) *
//             (X.row(i) - X.row(j)).t()
//         )
//       );
//     }
//   }
//   
//   return C_phi;
// }
// 
// double GP_log_density(
//     arma::vec w,
//     double sigmasq,
//     arma::mat C_phi
// ) {
//   arma::mat Sigma = sigmasq * C_phi;
//   
//   double log_likelihood =  -0.5 * logdet(Sigma) -
//     0.5 * arma::as_scalar(w.t() * inv_Chol(Sigma) * w);
//   
//   return log_likelihood;
// }
// 
// arma::vec beta_update(
//     arma::vec Y,
//     arma::mat X,
//     double tausq,
//     arma::vec w,
//     arma::mat V
// ) {
//   // Compute the posterior mean and covariance of beta
//   arma::mat Sigma = inv_Chol(X.t() * X / tausq + inv_Chol(V));
//   arma::vec mu = Sigma * (X.t() * (Y - w) / tausq);
//   // Draw beta from the multivariate normal distribution
//   arma::vec beta = arma::mvnrnd(mu, Sigma);
// 
//   return beta;
// }
// 
// arma::vec w_update(
//   arma::vec Y,
//   arma::mat X,
//   arma::vec beta,
//   double sigmasq,
//   double tausq,
//   arma::mat C_phi,
//   int n
// ) {
//   // Compute the posterior mean and covariance of w
//   arma::mat Sigma = inv_Chol(arma::eye(n, n) / tausq + inv_Chol(C_phi) / sigmasq);
//   arma::vec mu = Sigma * (Y - X * beta) / tausq;
//   // Draw w from the multivariate normal distribution
//   arma::vec w = arma::mvnrnd(mu, Sigma);
// 
//   return w;
// }
// 
// double sigmasq_update(
//   arma::vec w,
//   arma::mat C_phi,
//   double a_s,
//   double b_s,
//   int n
// ) {
//   // Compute the posterior parameters for sigmasq
//   double a_post = a_s + n / 2;
//   double b_post = b_s + 0.5 * arma::as_scalar(w.t() * inv_Chol(C_phi) * w);
// 
//   // Draw sigmasq from the inverse gamma distribution
//   double sigmasq = 1 / arma::randg(arma::distr_param(a_post, 1 / b_post));
// 
//   return sigmasq;
// }
// 
// double tausq_update(
//   arma::vec Y,
//   arma::mat X,
//   arma::vec beta,
//   arma::vec w,
//   double a_t,
//   double b_t,
//   int n
// ) {
//   // Compute the posterior parameters for tausq
//   double a_post = a_t + n / 2;
//   double b_post = b_t + 0.5 * arma::as_scalar((Y - X * beta - w).t() * (Y - X * beta - w));
// 
//   // Draw tausq from the inverse gamma distribution
//   double tausq = 1 / arma::randg(arma::distr_param(a_post, 1 / b_post));
// 
//   return tausq;
// }
// 
// //[[Rcpp::export]]
// Rcpp::List GP_Gibbs(
//     arma::vec Y,
//     arma::mat X,
//     double sigmasq_start,
//     double tausq_start,
//     arma::vec phi_start,
//     arma::vec w_start,
//     arma::vec beta_start,
//     arma::mat V,
//     double a_s,
//     double b_s,
//     double a_t,
//     double b_t,
//     arma::vec lower,
//     arma::vec upper,
//     int mcmc = 1000
// ){
//   int n = Y.n_elem;
//   int p = X.n_cols;
// 
//   arma::mat beta_samples = arma::zeros(mcmc, p);
//   arma::mat w_samples = arma::zeros(mcmc, n);
//   arma::vec sigmasq_samples = arma::zeros(mcmc);
//   arma::vec tausq_samples = arma::zeros(mcmc);
//   arma::mat phi_samples = arma::zeros(mcmc, p);
// 
//   arma::mat phi_bounds = arma::zeros(p, 2);
//   phi_bounds.col(0) = lower;
//   phi_bounds.col(1) = upper;
// 
//   arma::vec beta_cur = beta_start;
//   arma::vec w_cur = w_start;
//   double sigmasq_cur = sigmasq_start;
//   double tausq_cur = tausq_start;
//   arma::vec phi_cur = phi_start;
//   arma::mat C_phi = arma::zeros(n, n);
// 
//   arma::mat phi_metrop_sd = 0.05 * arma::eye(p, p);
//   RAMAdapt phi_adapt(p, phi_metrop_sd, 0.234);
// 
//   for(int i=0; i<mcmc; i++){
// 
//     phi_adapt.count_proposal();
// 
//     Rcpp::RNGScope scope;
// 
//     arma::vec U_update = arma::randn(p);
// 
//     arma::vec phi_alt = par_huvtransf_back(par_huvtransf_fwd(
//       phi_cur, phi_bounds) + phi_adapt.paramsd * U_update, phi_bounds);
// 
//     arma::mat C_phi_cur = calc_C_phi(X, phi_cur, n);
//     arma::mat C_phi_alt = calc_C_phi(X, phi_alt, n);
// 
//     double curr_logdens = GP_log_density(w_cur, sigmasq_cur, C_phi_cur);
//     double prop_logdens = GP_log_density(w_cur, sigmasq_cur, C_phi_alt);
//       
//     // make move
//     double jacobian  = calc_jacobian(phi_alt, phi_cur, phi_bounds);
//     double logaccept = prop_logdens - curr_logdens + jacobian;
// 
//     bool accepted = do_I_accept(logaccept);
// 
//     if(accepted){
//       phi_cur = phi_alt;
//       C_phi_cur = C_phi_alt;
//       phi_adapt.count_accepted();
//     }
// 
//     phi_adapt.update_ratios();
// 
//     if(true){
//       phi_adapt.adapt(U_update, exp(logaccept), i);
//     }
// 
//     phi_samples.row(i) = phi_cur.t();
//     beta_cur = beta_update(Y, X, tausq_cur, w_cur, V);
//     beta_samples.row(i) = beta_cur.t();
//     w_cur = w_update(Y, X, beta_cur, sigmasq_cur, tausq_cur, C_phi_cur, n);
//     w_samples.row(i) = w_cur.t();
//     sigmasq_cur = sigmasq_update(w_cur, C_phi_cur, a_s, b_s, n);
//     sigmasq_samples(i) = sigmasq_cur;
//     tausq_cur = tausq_update(Y, X, beta_cur, w_cur, a_t, b_t, n);
//     tausq_samples(i) = tausq_cur;
//   }
// 
//   Rcpp::List output;
//   output["beta"] = beta_samples;
//   output["w"] = w_samples;
//   output["sigmasq"] = sigmasq_samples;
//   output["tausq"] = tausq_samples;
//   output["phi"] = phi_samples;
//   phi_adapt.print_acceptance();
// 
//   return output;
//   
// }