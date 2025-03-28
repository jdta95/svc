# Data prep

library(ggplot2)
library(gridExtra)

set.seed(123)

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
p = 2

# generating w
sigmasq_w = 1
phi_w = 8
C_w = calc_C_phi(coords, phi_w)
w = MASS::mvrnorm(1, rep(0, n), sigmasq_w * C_w)

# generating beta's
sigmasq_1 = 1
phi_1 = 2
C_1 = calc_C_phi(coords, phi_1)
beta_1 = MASS::mvrnorm(1, rep(0, n), sigmasq_1 * C_1)

sigmasq_2 = 1
phi_2 = 2
C_2 = calc_C_phi(coords, phi_2)
beta_2 = MASS::mvrnorm(1, rep(0, n), sigmasq_2 * C_2)

# generating X
X_1 = rnorm(n, 0, 10)
X_2 = rnorm(n, 0, 10)

# generating epsilon
tausq = 0.0001
epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))

Y = X_1 * beta_1 + X_2 * beta_2 + w + epsilon

# create plot of Y according to coordinates
Y_plot = ggplot(data = data.frame(coords, Y)) +
  geom_tile(aes(x = lon, y = lat, fill = Y)) +
  guides(fill = guide_legend(title = "Y")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

w_plot = ggplot(data = data.frame(coords, w)) +
  geom_tile(aes(x = lon, y = lat, fill = w)) +
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
grid.arrange(Y_plot, w_plot, beta_1_plot, beta_2_plot, epsilon_plot, ncol = 3)

# generate knots
## 1 in every k^2 point on a grid is a knot
k = 5

lat_knots = unique(lat)
lat_knots = lat_knots[seq(1, length(lat_knots), by = k)]

lon_knots = unique(lon)
lon_knots = lon_knots[seq(1, length(lon_knots), by = k)]

knots = as.matrix(expand.grid(lat_knots, lon_knots))

df = data.frame(coords, Y, X_1, X_2, w, beta_1, beta_2, epsilon)

knots_df = data.frame(knots)
colnames(knots_df) = c("lat", "lon")

knots_df = merge(knots_df, df, by = c("lat", "lon"), all.x = TRUE, all.y = FALSE)

X = as.matrix(df[, c("X_1", "X_2")])
Y_star = knots_df$Y
X_star = as.matrix(knots_df[, c("X_1", "X_2")])

m = nrow(knots)
l = 0
u = 20

# Test function

# compile package
Rcpp::compileAttributes()
# load in library for testing
devtools::load_all()

output = svc::GP_Gibbs_phi_test(
  Y = Y,
  X = X,
  s = coords,
  knots = knots,
  beta_knots_start = as.matrix(knots_df[, c("beta_1", "beta_2")]),
  w_knots_start = knots_df$w,
  phi_beta_start = c(phi_1, phi_2),
  phi_w_start = l + (u - l) / 2,
  sigmasq_beta_start = c(sigmasq_1, sigmasq_2),
  sigmasq_w_start = sigmasq_w,
  tausq_start = tausq,
  phi_beta_proposal_sd = rep(0.1, p),
  phi_w_proposal_sd = 0.2,
  lower_beta = rep(l, p),
  upper_beta = rep(u, p),
  lower_w = l,
  upper_w = u,
  mcmc = 5000
)

start = 3000
end = 5000

# trace plot of phi_w
plot(y = output$phi_w_samples[start:end], x = start:end, type = "l", main = "Trace plot of phi_w", ylab = "phi_w", xlab = "Iteration")
mean(output$phi_w_samples[start:end])
quantile(output$phi_w_samples[start:end], c(0.025, 0.975))
phi_w

# # trace plot of phi_beta 1
# plot(y = output$phi_beta_samples[start:end, 1], x = start:end, type = "l", main = "Trace plot of phi_beta_1", ylab = "phi_beta_1", xlab = "Iteration")
# mean(output$phi_beta_samples[start:end, 1])
# phi_1
# 
# # trace plot of phi_beta 2
# plot(y = output$phi_beta_samples[start:end, 2], x = start:end, type = "l", main = "Trace plot of phi_beta_2", ylab = "phi_beta_2", xlab = "Iteration")
# mean(output$phi_beta_samples[start:end, 2])
# phi_2

