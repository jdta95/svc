#include <RcppArmadillo.h>


// Transform (0, inf) to real line (log transform)
double log_transform(double x) {
  return log(x);  // f(x) = log(x)
}

// Transform real line to (0, inf) (exp transform)
double exp_transform(double x) {
  return exp(x);  // f⁻¹(x) = exp(x)
}


// Transform all parameters to real line
arma::vec transform_theta(const arma::vec& theta, const arma::vec& phi_bounds_lower, 
                          const arma::vec& phi_bounds_upper, unsigned int p) {
  arma::vec theta_trans(theta.n_elem);
  
  // Transform variance parameters (tau^2, sigma_1^2, ..., sigma_p^2)
  for(unsigned int i = 0; i <= p; i++) {
    theta_trans(i) = log_transform(theta(i));
  }
  
  // Transform phi parameters (phi_1, ..., phi_p)
  for(unsigned int i = 1; i <= p; i++) {
    theta_trans(p + i) = logit(theta(p + i), phi_bounds_lower(i-1), phi_bounds_upper(i-1));
  }
  
  return theta_trans;
}


// Transform back from real line to original scale
arma::vec inverse_transform_theta(const arma::vec& theta_trans, const arma::vec& phi_bounds_lower,
                                  const arma::vec& phi_bounds_upper, unsigned int p) {
  arma::vec theta(theta_trans.n_elem);
  
  // Transform back variance parameters
  for(unsigned int i = 0; i <= p; i++) {
    theta(i) = exp_transform(theta_trans(i));
  }
  
  // Transform back phi parameters
  for(unsigned int i = 1; i <= p; i++) {
    theta(p + i) = logistic(theta_trans(p + i), phi_bounds_lower(i-1), phi_bounds_upper(i-1));
  }
  
  return theta;
}


// Function to calculate Jacobian

double calc_jacobian_theta(const arma::vec& theta, const arma::vec& phi_bounds_lower,
                           const arma::vec& phi_bounds_upper, unsigned int p) {
  double jac = 0.0;
  
  // Jacobian for variance parameters (log transform)
  for(unsigned int i = 0; i <= p; i++) {
    jac += log(theta(i)); // Jacobian for log(x) is 1/x, so log-Jacobian is -log(x)
  }
  
  // Jacobian for phi parameters (logit transform)
  for(unsigned int i = 1; i <= p; i++) {
    double phi = theta(p + i);
    double l = phi_bounds_lower(i-1);
    double u = phi_bounds_upper(i-1);
    jac += normal_proposal_logitscale(phi, l, u);
  }
  
  return jac;
}


// Generate proposal functions

// Update variance parameters (tau^2, sigma_1^2, ..., sigma_p^2)
arma::vec update_variances(const arma::vec& variances_current, const arma::vec& proposal_sd) {
  arma::vec variances_prop(variances_current.n_elem);
  for(unsigned int i = 0; i < variances_current.n_elem; i++) {
    double current_trans = log_transform(variances_current(i));
    double prop_trans = current_trans + proposal_sd(i) * arma::randn();
    variances_prop(i) = exp_transform(prop_trans);
  }
  return variances_prop;
}

// Update phi parameters (phi_1, ..., phi_p)
arma::vec update_phis(const arma::vec& phis_current, const arma::vec& proposal_sd,
                      const arma::vec& phi_bounds_lower, const arma::vec& phi_bounds_upper) {
  arma::vec phis_prop(phis_current.n_elem);
  for(unsigned int i = 0; i < phis_current.n_elem; i++) {
    double current_trans = logit(phis_current(i), phi_bounds_lower(i), phi_bounds_upper(i));
    double prop_trans = current_trans + proposal_sd(i) * arma::randn();
    phis_prop(i) = logistic(prop_trans, phi_bounds_lower(i), phi_bounds_upper(i));
  }
  return phis_prop;
}



// Add these declarations in other utility functions
//double log_transform(double x);
//double exp_transform(double x);
//double logit(double x, double l, double u);
//double logistic(double x, double l, double u);
//arma::vec transform_theta(const arma::vec& theta, 
  //                        const arma::vec& phi_bounds_lower,
    //                      const arma::vec& phi_bounds_upper, 
      //                    unsigned int p);
//arma::vec inverse_transform_theta(const arma::vec& theta_trans,
  //                                const arma::vec& phi_bounds_lower,
    //                              const arma::vec& phi_bounds_upper,
      //                            unsigned int p);