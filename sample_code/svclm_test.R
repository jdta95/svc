# Data prep
library(ggplot2)
library(gridExtra)

calc_C_phi = function(coords, phi) {
  n = nrow(coords)
  dists = as.matrix(dist(coords)) ^ 2
  C_phi = exp(-0.5 / phi * dists)
}

set.seed(123)

lat = seq(0, 30, by = 1)
lon = seq(0, 30, by = 1)

coords = as.matrix(expand.grid(lat, lon))
colnames(coords) = c("lat", "lon")

n = nrow(coords)

sigmasq_0 = 2
phi_0 = 8
C_0 = calc_C_phi(coords, phi_0)
mean_0 = 0
beta_0 = MASS::mvrnorm(1, rep(mean_0, n), sigmasq_0 * C_0)

sigmasq_1 = 6
phi_1 = 5
C_1 = calc_C_phi(coords, phi_1)
mean_1 = 0
beta_1 = MASS::mvrnorm(1, rep(mean_1, n), sigmasq_1 * C_1)

sigmasq_2 = 4
phi_2 = 3
C_2 = calc_C_phi(coords, phi_1)
mean_2 = 0
beta_2 = MASS::mvrnorm(1, rep(mean_2, n), sigmasq_2 * C_2)

# generating X
X_0 = rep(1, n)
X_1 = rnorm(n)
X_2 = rnorm(n)

# generating epsilon
tausq = 2
epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))

Y = X_0 * beta_0 + X_1 * beta_1 + X_2 * beta_2 + epsilon

# create plot of Y according to coordinates
Y_plot = ggplot(data = data.frame(coords, Y)) +
  geom_tile(aes(x = lon, y = lat, fill = Y)) +
  guides(fill = guide_legend(title = "Y")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

beta_0_plot = ggplot(data = data.frame(coords, beta_0)) +
  geom_tile(aes(x = lon, y = lat, fill = beta_0)) +
  guides(fill = guide_legend(title = "beta_0")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

beta_1_plot = ggplot(data = data.frame(coords, beta_1)) +
  geom_tile(aes(x = lon, y = lat, fill = beta_1)) +
  guides(fill = guide_legend(title = "beta_1")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

beta_2_plot = ggplot(data = data.frame(coords, beta_2)) +
  geom_tile(aes(x = lon, y = lat, fill = beta_2)) +
  guides(fill = guide_legend(title = "beta_2")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

epsilon_plot = ggplot(data = data.frame(coords, epsilon)) +
  geom_tile(aes(x = lon, y = lat, fill = epsilon)) +
  guides(fill = guide_legend(title = "epsilon")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

# show all plots in a grid
grid.arrange(Y_plot, beta_0_plot, beta_1_plot, beta_2_plot, epsilon_plot, ncol = 2)

# generate knots
## 1 in every k^2 point on a grid is a knot
k = 2

lat_knots = unique(lat)
lat_knots = lat_knots[seq(1, length(lat_knots), by = k)]

lon_knots = unique(lon)
lon_knots = lon_knots[seq(1, length(lon_knots), by = k)]

knots = as.matrix(expand.grid(lat_knots, lon_knots))

df = data.frame(coords, Y, X_0, X_1, X_2, beta_0, beta_1, beta_2, epsilon)

knots_df = data.frame(knots)
colnames(knots_df) = c("lat", "lon")

knots_df = merge(knots_df, df, by = c("lat", "lon"), all.x = TRUE, all.y = FALSE)

X = as.matrix(df[, c("X_0", "X_1", "X_2")])
Y_knots = knots_df$Y
X_knots = as.matrix(knots_df[, c("X_0", "X_1", "X_2")])
beta_knots = as.matrix(knots_df[, c("beta_0", "beta_1", "beta_2")])
knots = as.matrix(knots_df[, c("lat", "lon")])

Y_knots_plot = ggplot(data = knots_df) +
  geom_tile(aes(x = lon, y = lat, fill = Y)) +
  guides(fill = guide_legend(title = "Y knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

beta_0_knots_plot = ggplot(data = knots_df) +
  geom_tile(aes(x = lon, y = lat, fill = beta_0)) +
  guides(fill = guide_legend(title = "beta_0 knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

beta_1_knots_plot = ggplot(data = knots_df) +
  geom_tile(aes(x = lon, y = lat, fill = beta_1)) +
  guides(fill = guide_legend(title = "beta_1 knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

beta_2_knots_plot = ggplot(data = knots_df) +
  geom_tile(aes(x = lon, y = lat, fill = beta_2)) +
  guides(fill = guide_legend(title = "beta_2 knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

grid.arrange(Y_plot, Y_knots_plot, ncol = 2)
grid.arrange(beta_0_plot, beta_0_knots_plot, ncol = 2)
grid.arrange(beta_1_plot, beta_1_knots_plot, ncol = 2)
grid.arrange(beta_2_plot, beta_2_knots_plot, ncol = 2)

l = 0
u = 10

# Test function

# compile package
Rcpp::compileAttributes()
# load in library for testing
devtools::load_all()

mcmc = 3000

start = mcmc * 0.5
end = mcmc

p = ncol(X)

# test with out knowing true values
output = svc::svclm(
  Y = Y,
  X = X,
  s = coords,
  Y_knots = Y_knots,
  X_knots = X_knots,
  knots = knots,
  beta_knots_start = matrix(0, nrow = nrow(knots), ncol = p),
  phi_beta_start = rep(l + (u - l) / 2, p),
  sigmasq_beta_start = rep(2, p),
  tausq_start = 2,
  # phi_beta_proposal_sd = rep(0.5, p),
  phi_beta_proposal_sd = c(0.5, 0.5, 0.5),
  phi_beta_lower = rep(l, p),
  phi_beta_upper = rep(u, p),
  a_beta = rep(2, p),
  b_beta = rep(2, p),
  a_t = 2,
  b_t = 2,
  mcmc = mcmc
)

my_summary = function(x){
  return (c(
    quantile(x, 0.025),
    mean(x),
    quantile(x, 0.975)
  ))
}

# phi_0
plot(y = output$phi_beta_samples[1:end, 1], x = 1:end, type = "l", main = "Trace plot of phi_0", ylab = "phi_0", xlab = "Iteration")
mean(output$phi_beta_acceptance[, 1])
my_summary(output$phi_beta_samples[start:end, 1])
phi_0

# phi_1
plot(y = output$phi_beta_samples[1:end, 2], x = 1:end, type = "l", main = "Trace plot of phi_1", ylab = "phi_1", xlab = "Iteration")
mean(output$phi_beta_acceptance[, 2])
my_summary(output$phi_beta_samples[start:end, 2])
phi_1

# phi_2
plot(y = output$phi_beta_samples[1:end, 3], x = 1:end, type = "l", main = "Trace plot of phi_2", ylab = "phi_2", xlab = "Iteration")
mean(output$phi_beta_acceptance[, 3])
my_summary(output$phi_beta_samples[start:end, 3])
phi_2

# sigma_0
plot(y = output$sigmasq_beta_samples[1:end, 1], x = 1:end, type = "l", main = "Trace plot of sigmasq_0", ylab = "sigmasq_0", xlab = "Iteration")
my_summary(output$sigmasq_beta_samples[start:end, 1])
sigmasq_0

# sigma_1
plot(y = output$sigmasq_beta_samples[1:end, 2], x = 1:end, type = "l", main = "Trace plot of sigmasq_1", ylab = "sigmasq_1", xlab = "Iteration")
my_summary(output$sigmasq_beta_samples[start:end, 2])
sigmasq_1

# sigma_2
plot(y = output$sigmasq_beta_samples[1:end, 3], x = 1:end, type = "l", main = "Trace plot of sigmasq_2", ylab = "sigmasq_2", xlab = "Iteration")
my_summary(output$sigmasq_beta_samples[start:end, 3])
sigmasq_2

# tau^2
plot(y = output$tausq_samples[1:end], x = 1:end, type = "l", main = "Trace plot of tau^2", ylab = "tau^2", xlab = "Iteration")
my_summary(output$tausq_samples[start:end])
tausq

# beta_0
beta_0_pred_plot = ggplot(data = data.frame(coords, betahat0 = colMeans(output$beta_samples[start:end, , 1]))) +
  geom_tile(aes(x = lon, y = lat, fill = betahat0)) +
  guides(fill = guide_legend(title = "beta_0")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

# beta_1
beta_1_pred_plot = ggplot(data = data.frame(coords, betahat1 = colMeans(output$beta_samples[start:end, , 2]))) +
  geom_tile(aes(x = lon, y = lat, fill = betahat1)) +
  guides(fill = guide_legend(title = "beta_1")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

# beta_2
beta_2_pred_plot = ggplot(data = data.frame(coords, betahat2 = colMeans(output$beta_samples[start:end, , 3]))) +
  geom_tile(aes(x = lon, y = lat, fill = betahat2)) +
  guides(fill = guide_legend(title = "beta_2")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

grid.arrange(beta_0_plot, beta_0_pred_plot, ncol = 2)
grid.arrange(beta_1_plot, beta_1_pred_plot, ncol = 2)
grid.arrange(beta_2_plot, beta_2_pred_plot, ncol = 2)

# arbitrary beta_0(s)
s = n / 2
plot(y = output$beta_samples[1:end, s, 1], x = 1:end, type = "l", main = "Trace plot of beta_0", ylab = "beta_0", xlab = "Iteration")
my_summary(output$beta_samples[start:end, s, 1])
beta_0[s]

