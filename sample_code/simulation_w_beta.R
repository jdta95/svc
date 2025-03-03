# Load necessary libraries
library(Rcpp)
library(ggplot2)
library(gridExtra)

# Compile and load the C++ functions
Rcpp::sourceCpp("C:/Users/s-edw/Downloads/svc-main/svc-main/src/beta_update.cpp")

# Function to calculate covariance matrix
calc_C_phi <- function(coords, phi) {
  n <- nrow(coords)
  C_phi <- matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      C_phi[i, j] <- exp(-0.5 * (sum((
        coords[i, ] - coords[j, ]
      ) ^ 2) / phi))
    }
  }
  return(C_phi)
}

# Simulate new data
set.seed(123)

lat <- seq(0, 2, by = 1)
lon <- seq(0, 2, by = 1)
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
k <- 1
lat_knots <- unique(lat)
lat_knots <- lat_knots[seq(1, length(lat_knots), by = k)]
lon_knots <- unique(lon)
lon_knots <- lon_knots[seq(1, length(lon_knots), by = k)]
knots <- as.matrix(expand.grid(lat_knots, lon_knots))


# Subset beta_1 and beta_2 to match the dimensions of the knots
beta_1_knots <-
  beta_1[1:nrow(knots)]  # Use the first 16 values (assuming 16 knots)
beta_2_knots <- beta_2[1:nrow(knots)]  # Use the first 16 values

# Number of MCMC iterations
mcmc <- 10

# Storage for posterior samples
beta1_samples <- numeric(mcmc)
beta2_samples <- numeric(mcmc)
w_samples <- numeric(mcmc)

# Initial values
beta1 <- rep(1,n)  # Initial value for beta1
beta2 <- rep(1,n)  # Initial value for beta2
w_new <- rep(1,n)

# Store initial values in the first iteration
beta1_samples[1] <- beta1
beta2_samples[1] <- beta2
w_samples[1] <- w_new

# Gibbs sampling loop
for (i in 2:mcmc) {
  beta_1_knots <-
    beta1[1:nrow(knots)]
  
  #Update beta1
  beta1 <-
    update_beta(Y, X_1, beta_1_knots, w, knots, sigmasq_1, phi_1, tausq)
  beta1 <- calc_beta(coords, knots, phi_1, beta1)
  beta1_samples[i] <- beta1
  
  beta_1_knots <-
    beta2[1:nrow(knots)]
  
  #Update beta2
  beta2 <-
    update_beta(Y,  X_2, beta_2_knots, w, knots, sigmasq_2, phi_2, tausq)
  beta2 <- calc_beta(coords, knots, phi_2, beta2)
  beta2_samples[i] <- beta2
  
  #Update w
  w_new <-
    update_w(
      Y,
      cbind(X_1, X_2),
      cbind(beta_1_knots, beta_2_knots),
      sigmasq_w,
      phi_w,
      tausq,
      knots
    )
  w_new <- calc_w(coords, knots, phi_w, w_new)
  w_samples[i] <- w_new
}

# Print initial values
cat("Initial W:", w_samples[1], "\n")
cat("Initial beta1:", beta1_samples[1], "\n")
cat("Initial beta2:", beta22_samples[1], "\n")

# Print first few samples
cat("First few w:", w_samples[1:5], "\n")
cat("First few beta1 samples:", beta1_samples[1:5], "\n")
cat("First few beta2 samples:", beta2_samples[1:5], "\n")

# Plot posterior samples
par(mfrow = c(3, 1))
plot(
  w_samples,
  type = "l",
  main = "Trace of w",
  xlab = "Iteration",
  ylab = "w"
)
plot(
  beta1_samples,
  type = "l",
  main = "Trace of beta_1",
  xlab = "Iteration",
  ylab = "beta1"
)
plot(
  beta2_2_samples,
  type = "l",
  main = "Trace of beta_2",
  xlab = "Iteration",
  ylab = "beta2"
)