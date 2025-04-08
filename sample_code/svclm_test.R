# Data prep
library(ggplot2)
library(gridExtra)

calc_C_phi = function(coords, phi) {
  n = nrow(coords)
  dists = as.matrix(dist(coords)) ^ 2
  C_phi = exp(-0.5 / phi * dists)
}

set.seed(123)

lat = seq(0, 20, by = 1)
lon = seq(0, 20, by = 1)

coords = as.matrix(expand.grid(lat, lon))
colnames(coords) = c("lat", "lon")

n = nrow(coords)

sigmasq_0 = 2
phi_0 = 13
C_0 = calc_C_phi(coords, phi_0)
mean_0 = 3
w_0 = MASS::mvrnorm(1, rep(mean_0, n), sigmasq_0 * C_0)

sigmasq_1 = 6
phi_1 = 7
C_1 = calc_C_phi(coords, phi_1)
mean_1 = 0
w_1 = MASS::mvrnorm(1, rep(mean_1, n), sigmasq_1 * C_1)

sigmasq_2 = 4
phi_2 = 10
C_2 = calc_C_phi(coords, phi_2)
mean_2 = -2
w_2 = MASS::mvrnorm(1, rep(mean_2, n), sigmasq_2 * C_2)

# generating X
X_0 = rep(1, n)
X_1 = rnorm(n)
X_2 = rnorm(n)

# generating epsilon
tausq = 0.01
epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))

Y = X_0 * w_0 + X_1 * w_1 + X_2 * w_2 + epsilon

# create plot of Y according to coordinates
Y_plot = ggplot(data = data.frame(coords, Y)) +
  geom_tile(aes(x = lon, y = lat, fill = Y)) +
  guides(fill = guide_legend(title = "Y")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

w_0_plot = ggplot(data = data.frame(coords, w_0)) +
  geom_tile(aes(x = lon, y = lat, fill = w_0)) +
  guides(fill = guide_legend(title = "w_0")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

w_1_plot = ggplot(data = data.frame(coords, w_1)) +
  geom_tile(aes(x = lon, y = lat, fill = w_1)) +
  guides(fill = guide_legend(title = "w_1")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

w_2_plot = ggplot(data = data.frame(coords, w_2)) +
  geom_tile(aes(x = lon, y = lat, fill = w_2)) +
  guides(fill = guide_legend(title = "w_2")) +
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
grid.arrange(Y_plot, w_0_plot, w_1_plot, w_2_plot, epsilon_plot, ncol = 2)

# generate knots
## 1 in every k^2 point on a grid is a knot
k = 2

lat_knots = unique(lat)
lat_knots = lat_knots[seq(1, length(lat_knots), by = k)]

lon_knots = unique(lon)
lon_knots = lon_knots[seq(1, length(lon_knots), by = k)]

knots = as.matrix(expand.grid(lat_knots, lon_knots))

df = data.frame(coords, Y, X_0, X_1, X_2, w_0, w_1, w_2, epsilon)

knots_df = data.frame(knots)
colnames(knots_df) = c("lat", "lon")

knots_df = merge(knots_df, df, by = c("lat", "lon"), all.x = TRUE, all.y = FALSE)

X = as.matrix(df[, c("X_0", "X_1", "X_2")])
Y_knots = knots_df$Y
X_knots = as.matrix(knots_df[, c("X_0", "X_1", "X_2")])
w_knots = as.matrix(knots_df[, c("w_0", "w_1", "w_2")])
knots = as.matrix(knots_df[, c("lat", "lon")])

Y_knots_plot = ggplot(data = knots_df) +
  geom_tile(aes(x = lon, y = lat, fill = Y)) +
  guides(fill = guide_legend(title = "Y knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

w_0_knots_plot = ggplot(data = knots_df) +
  geom_tile(aes(x = lon, y = lat, fill = w_0)) +
  guides(fill = guide_legend(title = "w_0 knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

w_1_knots_plot = ggplot(data = knots_df) +
  geom_tile(aes(x = lon, y = lat, fill = w_1)) +
  guides(fill = guide_legend(title = "w_1 knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

w_2_knots_plot = ggplot(data = knots_df) +
  geom_tile(aes(x = lon, y = lat, fill = w_2)) +
  guides(fill = guide_legend(title = "w_2 knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

grid.arrange(Y_plot, Y_knots_plot, ncol = 2)
grid.arrange(w_0_plot, w_0_knots_plot, ncol = 2)
grid.arrange(w_1_plot, w_1_knots_plot, ncol = 2)
grid.arrange(w_2_plot, w_2_knots_plot, ncol = 2)

l = 0
u = 30

# Test function

# compile package
Rcpp::compileAttributes()
# load in library for testing
devtools::load_all()

mcmc = 2000

start = mcmc * 0.5
end = mcmc

p = ncol(X)

# test with uninformative priors
output = svc::svclm(
  Y = Y,
  X = X,
  s = coords,
  Y_knots = Y_knots,
  X_knots = X_knots,
  knots = knots,
  phi_lower = rep(l, p),
  phi_upper = rep(u, p),
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
plot(y = output$phi_samples[1:end, 1], x = 1:end, type = "l", main = "Trace plot of phi_0", ylab = "phi_0", xlab = "Iteration")
mean(output$phi_acceptance[, 1])
my_summary(output$phi_samples[start:end, 1])
phi_0

# phi_1
plot(y = output$phi_samples[1:end, 2], x = 1:end, type = "l", main = "Trace plot of phi_1", ylab = "phi_1", xlab = "Iteration")
mean(output$phi_acceptance[, 2])
my_summary(output$phi_samples[start:end, 2])
phi_1

# phi_2
plot(y = output$phi_samples[1:end, 3], x = 1:end, type = "l", main = "Trace plot of phi_2", ylab = "phi_2", xlab = "Iteration")
mean(output$phi_acceptance[, 3])
my_summary(output$phi_samples[start:end, 3])
phi_2

# sigma_0
plot(y = output$sigmasq_samples[1:end, 1], x = 1:end, type = "l", main = "Trace plot of sigmasq_0", ylab = "sigmasq_0", xlab = "Iteration")
my_summary(output$sigmasq_samples[start:end, 1])
sigmasq_0

# sigma_1
plot(y = output$sigmasq_samples[1:end, 2], x = 1:end, type = "l", main = "Trace plot of sigmasq_1", ylab = "sigmasq_1", xlab = "Iteration")
my_summary(output$sigmasq_samples[start:end, 2])
sigmasq_1

# sigma_2
plot(y = output$sigmasq_samples[1:end, 3], x = 1:end, type = "l", main = "Trace plot of sigmasq_2", ylab = "sigmasq_2", xlab = "Iteration")
my_summary(output$sigmasq_samples[start:end, 3])
sigmasq_2

# tau^2
plot(y = output$tausq_samples[1:end], x = 1:end, type = "l", main = "Trace plot of tau^2", ylab = "tau^2", xlab = "Iteration")
my_summary(output$tausq_samples[start:end])
tausq

# w_0
w_0_pred_plot = ggplot(data = data.frame(coords, what0 = colMeans(output$w_samples[start:end, , 1]))) +
  geom_tile(aes(x = lon, y = lat, fill = what0)) +
  guides(fill = guide_legend(title = "w_0")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

# w_1
w_1_pred_plot = ggplot(data = data.frame(coords, what1 = colMeans(output$w_samples[start:end, , 2]))) +
  geom_tile(aes(x = lon, y = lat, fill = what1)) +
  guides(fill = guide_legend(title = "w_1")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

# w_2
w_2_pred_plot = ggplot(data = data.frame(coords, what2 = colMeans(output$w_samples[start:end, , 3]))) +
  geom_tile(aes(x = lon, y = lat, fill = what2)) +
  guides(fill = guide_legend(title = "w_2")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

grid.arrange(w_0_plot, w_0_pred_plot, ncol = 2)
grid.arrange(w_1_plot, w_1_pred_plot, ncol = 2)
grid.arrange(w_2_plot, w_2_pred_plot, ncol = 2)

# arbitrary w_0(s)
s = n / 2
plot(y = output$w_samples[1:end, s, 1], x = 1:end, type = "l", main = "Trace plot of w_0", ylab = "w_0", xlab = "Iteration")
my_summary(output$w_samples[start:end, s, 1])
w_0[s]

