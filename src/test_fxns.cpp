#include <RcppArmadillo.h>


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# compile package
Rcpp::compileAttributes()
# load in library for testing
devtools::load_all()

mcmc = 1000

phi_samples = numeric(mcmc)
phi_samples[1] = 2

# test the function
for (i in 2:mcmc) {
  phi_samples[i] = phi_RW(
    knots = coords,
    w_knots = w,
    sigmasq_cur = sigmasq_w,
    phi_cur = phi_samples[i - 1],
    proposal_sd = 0.001,
    phi_bounds = matrix(c(0.00001, 4), nrow = 1)
  )
}

plot(phi_samples)

*/
