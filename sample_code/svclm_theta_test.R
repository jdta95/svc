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

sigmasq_w = 2
phi_w = 5
C_w = calc_C_phi(coords, phi_w)
mean_w = 0
beta_w = MASS::mvrnorm(1, rep(mean_w, n), sigmasq_w * C_w)

sigmasq_1 = 2
phi_1 = 3
C_1 = calc_C_phi(coords, phi_1)
mean_1 = 0
beta_1 = MASS::mvrnorm(1, rep(mean_1, n), sigmasq_1 * C_1)

# generating X
X_w = rep(1, n)
X_1 = rnorm(n)

# generating epsilon
tausq = 1
epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))

Y = X_w * beta_w + X_1 * beta_1 + epsilon

# create plot of Y according to coordinates
Y_plot = ggplot(data = data.frame(coords, Y)) +
  geom_tile(aes(x = lon, y = lat, fill = Y)) +
  guides(fill = guide_legend(title = "Y")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

w_plot = ggplot(data = data.frame(coords, beta_w)) +
  geom_tile(aes(x = lon, y = lat, fill = beta_w)) +
  guides(fill = guide_legend(title = "w")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

beta_1_plot = ggplot(data = data.frame(coords, beta_1)) +
  geom_tile(aes(x = lon, y = lat, fill = beta_1)) +
  guides(fill = guide_legend(title = "beta_1")) +
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
grid.arrange(Y_plot, w_plot, beta_1_plot, epsilon_plot, ncol = 2)

# generate knots
## 1 in every k^2 point on a grid is a knot
k = 2

lat_knots = unique(lat)
lat_knots = lat_knots[seq(1, length(lat_knots), by = k)]

lon_knots = unique(lon)
lon_knots = lon_knots[seq(1, length(lon_knots), by = k)]

knots = as.matrix(expand.grid(lat_knots, lon_knots))

df = data.frame(coords, Y, X_w, X_1, beta_w, beta_1, epsilon)

knots_df = data.frame(knots)
colnames(knots_df) = c("lat", "lon")

knots_df = merge(knots_df, df, by = c("lat", "lon"), all.x = TRUE, all.y = FALSE)

X = as.matrix(df[, c("X_w", "X_1")])
Y_knots = knots_df$Y
X_knots = as.matrix(knots_df[, c("X_w", "X_1")])
beta_knots = as.matrix(knots_df[, c("beta_w", "beta_1")])
knots = as.matrix(knots_df[, c("lat", "lon")])

Y_knots_plot = ggplot(data = knots_df) +
  geom_tile(aes(x = lon, y = lat, fill = Y)) +
  guides(fill = guide_legend(title = "Y knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

w_knots_plot = ggplot(data = knots_df) +
  geom_tile(aes(x = lon, y = lat, fill = beta_w)) +
  guides(fill = guide_legend(title = "w knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

beta_1_knots_plot = ggplot(data = knots_df) +
  geom_tile(aes(x = lon, y = lat, fill = beta_1)) +
  guides(fill = guide_legend(title = "beta_1 knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

grid.arrange(Y_plot, Y_knots_plot, ncol = 2)
grid.arrange(w_plot, w_knots_plot, ncol = 2)
grid.arrange(beta_1_plot, beta_1_knots_plot, ncol = 2)

l = 0
u = 10

# Test function

# compile package
Rcpp::compileAttributes()
# load in library for testing
devtools::load_all()

mcmc = 1000

start = mcmc * 0.5
end = mcmc

p = ncol(X)

# tau^2, sigma^2, and phi check: GOOD
output = svc::svclm_theta(
  Y = Y,
  X = X,
  s = coords,
  Y_knots = Y_knots,
  X_knots = X_knots,
  knots = knots,
  beta_knots_start = beta_knots,
  phi_beta_start = c(phi_w, phi_1),
  sigmasq_beta_start = rep(2, p),
  tausq_start = 2,
  phi_beta_proposal_sd = rep(0.3, p),
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

# phi_w
plot(y = output$phi_beta_samples[1:end, 1], x = 1:end, type = "l", main = "Trace plot of phi_w", ylab = "phi_w", xlab = "Iteration")
mean(output$phi_beta_acceptance[start:end, 1])
my_summary(output$phi_beta_samples[start:end, 1])
phi_w

# phi_1
plot(y = output$phi_beta_samples[1:end, 2], x = 1:end, type = "l", main = "Trace plot of phi_1", ylab = "phi_1", xlab = "Iteration")
mean(output$phi_beta_acceptance[start:end, 2])
my_summary(output$phi_beta_samples[start:end, 2])
phi_1

# sigma_w
plot(y = output$sigmasq_beta_samples[1:end, 1], x = 1:end, type = "l", main = "Trace plot of sigmasq_w", ylab = "sigmasq_w", xlab = "Iteration")
mean(output$sigmasq_beta_samples[start:end, 1])
my_summary(output$sigmasq_beta_samples[start:end, 1])
sigmasq_w

# sigma_1
plot(y = output$sigmasq_beta_samples[1:end, 2], x = 1:end, type = "l", main = "Trace plot of sigmasq_1", ylab = "sigmasq_1", xlab = "Iteration")
mean(output$sigmasq_beta_samples[start:end, 2])
my_summary(output$sigmasq_beta_samples[start:end, 2])
sigmasq_1

# tau^2
plot(y = output$tausq_samples[1:end], x = 1:end, type = "l", main = "Trace plot of tau^2", ylab = "tau^2", xlab = "Iteration")
mean(output$tausq_samples[start:end])
my_summary(output$tausq_samples[start:end])
tausq

