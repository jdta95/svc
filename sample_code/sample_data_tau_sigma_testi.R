# Load necessary libraries
library(Rcpp)
library(ggplot2)
library(gridExtra)
library(MASS)
library(svc)

# Function to calculate covariance matrix
calc_C_phi <- function(coords, phi) {
  n <- nrow(coords)
  C_phi <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      C_phi[i, j] <- exp(-0.5 * (sum((coords[i,] - coords[j,])^2) / phi))
    }
  }
  return(C_phi)
}


# Simulate new data
set.seed(123)

lat <- seq(0, 10, by = 1)
lon <- seq(0, 10, by = 1)
coords <- as.matrix(expand.grid(lat, lon))
colnames(coords) <- c("lat", "lon")

n <- nrow(coords)
p <- 2

# Generating w (spatial random effects)
sigmasq_w <- 0.1  # Smaller variance
phi_w <- 2
C_w <- calc_C_phi(coords, phi_w)
w <- MASS::mvrnorm(1, rep(0, n), sigmasq_w * C_w)

# Generating beta's (spatially-varying coefficients)
sigmasq_1 <- 0.1  # Smaller variance
phi_1 <- 2
C_1 <- calc_C_phi(coords, phi_1)
beta_1 <- MASS::mvrnorm(1, rep(0, n), sigmasq_1 * C_1)

sigmasq_2 <- 0.1  # Smaller variance
phi_2 <- 2
C_2 <- calc_C_phi(coords, phi_2)
beta_2 <- MASS::mvrnorm(1, rep(0, n), sigmasq_2 * C_2)

# Generating X (covariates)
X_1 <- rnorm(n, 0, 1)  # Smaller variance
X_2 <- rnorm(n, 0, 1)  # Smaller variance

# Generating epsilon (measurement error)
tausq <- 0.01  # Smaller variance
epsilon <- rnorm(n, mean = 0, sd = sqrt(tausq))

# Response variable
Y <- X_1 * beta_1 + X_2 * beta_2 + w + epsilon

# Generate knots
k <- 3
lat_knots <- unique(lat)
lat_knots <- lat_knots[seq(1, length(lat_knots), by = k)]
lon_knots <- unique(lon)
lon_knots <- lon_knots[seq(1, length(lon_knots), by = k)]
knots <- as.matrix(expand.grid(lat_knots, lon_knots))

# Prior parameters
a_t <- 1; b_t <- 1  # Prior for tau^2
a_r <- 1; b_r <- 1  # Prior for sigma_r^2

# Define beta as a vector of coefficients for X_1 and X_2
beta <- c(1, 1)  # Example coefficients for X_1 and X_2

# Subset beta_1 and beta_2 to match the dimensions of the knots
beta_1_knots <- beta_1[1:nrow(knots)]  # Use the first 16 values (assuming 16 knots)
beta_2_knots <- beta_2[1:nrow(knots)]  # Use the first 16 values

# Number of MCMC iterations
mcmc <- 1000

# Storage for posterior samples
tau2_samples <- numeric(mcmc)
sigma2_1_samples <- numeric(mcmc)
sigma2_2_samples <- numeric(mcmc)

# Initial values
tau2 <- 1  # Initial value for tau^2
sigma2_1 <- 1  # Initial value for sigma^2_1
sigma2_2 <- 1  # Initial value for sigma^2_2

# Store initial values in the first iteration
tau2_samples[1] <- tau2
sigma2_1_samples[1] <- sigma2_1
sigma2_2_samples[1] <- sigma2_2

# Run the Gibbs sampler with SCALED data
output <- svc::GP_Gibbs_test_tau_sigma(
  Y = Y,
  X = cbind(X_1, X_2),
  s = coords,
  knots = knots,
  beta_knots_start = as.matrix(knots_df[, c("beta_1", "beta_2")]),
  w_knots_start = knots_df$w,
  phi_beta_start = true_phi[1:2],
  phi_w_start = true_phi[3],
  sigmasq_beta_start = rep(1, p),
  sigmasq_w_start = 1,
  tausq_start = 1,
  phi_beta_proposal_sd = rep(0.1, p),
  phi_w_proposal_sd = 0.1,
  a_beta = rep(3, p),  # More informative priors
  b_beta = rep(0.1, p),
  a_w = 3,
  b_w = 0.1,
  a_t = 3,             # More informative prior for τ²
  b_t = 0.1,
  lower_beta = rep(0.1, p),
  upper_beta = rep(10, p),
  lower_w = 0.1,
  upper_w = 10,
  mcmc = 5000
)


# Analysis of results
burnin <- 3000
samples_keep <- burnin:5000
# Analysis of results
burnin <- 3000
samples_keep <- burnin:5000

# Plot tau^2 convergence
plot(output$tausq_samples, type = "l", main = "tau^2 samples")
abline(h = true_tausq, col = "red", lwd = 2)
cat("True tau^2:", true_tausq, "\nEstimated:", mean(output$tausq_samples[samples_keep]), 
    "\n95% CI:", quantile(output$tausq_samples[samples_keep], c(0.025, 0.975)), "\n")

# Plot sigma^2_w convergence
plot(output$sigmasq_w_samples, type = "l", main = "sigma^2_w samples")
abline(h = true_sigmasq_w, col = "red", lwd = 2)
cat("True sigma^2_w:", true_sigmasq_w, "\nEstimated:", mean(output$sigmasq_w_samples[samples_keep]),
    "\n95% CI:", quantile(output$sigmasq_w_samples[samples_keep], c(0.025, 0.975)), "\n")

# Plot sigma^2_beta convergence
for (j in 1:p) {
  plot(output$sigmasq_beta_samples[,j], type = "l", 
       main = paste("sigma^2_beta", j, "samples"))
  abline(h = true_sigmasq_beta[j], col = "red", lwd = 2)
  cat(paste0("\nTrue sigma^2_beta[", j, "]:"), true_sigmasq_beta[j], 
      "\nEstimated:", mean(output$sigmasq_beta_samples[samples_keep,j]),
      "\n95% CI:", quantile(output$sigmasq_beta_samples[samples_keep,j], c(0.025, 0.975)), "\n")
}