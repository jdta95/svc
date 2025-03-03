# Load necessary libraries
library(Rcpp)
library(ggplot2)
library(gridExtra)

# Compile and load the C++ functions
Rcpp::sourceCpp("Update_tau_sq.cpp")
Rcpp::sourceCpp("Update_sigma_sq.cpp")
Rcpp::sourceCpp("phi_RW.cpp")

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

set.seed(456)  # Different seed for generating a new dataset

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
phi_1_samples <- numeric(mcmc)
phi_2_samples <- numeric(mcmc)

# Initial values
tau2 <- 1  # Initial value for tau^2
sigma2_1 <- 1  # Initial value for sigma^2_1
sigma2_2 <- 1  # Initial value for sigma^2_2
phi_1 <- phi_1  # Initial value for phi_1
phi_2 <- phi_2  # Initial value for phi_2

# Store initial values in the first iteration
tau2_samples[1] <- tau2
sigma2_1_samples[1] <- sigma2_1
sigma2_2_samples[1] <- sigma2_2
phi_1_samples[1] <- phi_1
phi_2_samples[1] <- phi_2

# Gibbs sampling loop
for (i in 2:mcmc) {
  # Update tau^2
  tau2 <- update_tau2_r(Y, cbind(X_1, X_2), beta, w, a_t, b_t)
  tau2_samples[i] <- tau2
  
  # Update sigma^2 for beta_1
  sigma2_1 <- update_sigma2_r(beta_1_knots, a_r, b_r, phi_1, knots)
  sigma2_1_samples[i] <- sigma2_1
  
  # Update sigma^2 for beta_2
  sigma2_2 <- update_sigma2_r(beta_2_knots, a_r, b_r, phi_2, knots)
  sigma2_2_samples[i] <- sigma2_2
  
  # Update phi_1
  phi_1 <- phi_RW(knots, beta_1_knots, sigma2_1, phi_1, 0.1, matrix(c(0.1, 10), ncol = 2))
  phi_1_samples[i] <- phi_1
  
  # Update phi_2
  phi_2 <- phi_RW(knots, beta_2_knots, sigma2_2, phi_2, 0.1, matrix(c(0.1, 10), ncol = 2))
  phi_2_samples[i] <- phi_2
}


# Print first few samples    
cat("First few tau^2 samples:", tau2_samples[1:5], "\n")
cat("First few sigma^2_1 samples:", sigma2_1_samples[1:5], "\n")
cat("First few sigma^2_2 samples:", sigma2_2_samples[1:5], "\n")
cat("First few phi_1 samples:", phi_1_samples[1:5], "\n")
cat("First few phi_2 samples:", phi_2_samples[1:5], "\n")


plot(tau2_samples, type = "l", main = "Trace of tau^2", xlab = "Iteration", ylab = "tau^2")
plot(sigma2_1_samples, type = "l", main = "Trace of sigma^2 for beta_1", xlab = "Iteration", ylab = "sigma^2_1")
plot(sigma2_2_samples, type = "l", main = "Trace of sigma^2 for beta_2", xlab = "Iteration", ylab = "sigma^2_2")
plot(phi_1_samples, type = "l", main = "Trace of phi for beta_1", xlab = "Iteration", ylab = "phi_1")
plot(phi_2_samples, type = "l", main = "Trace of phi for beta_2", xlab = "Iteration", ylab = "phi_2")
