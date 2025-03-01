# # compile package
# Rcpp::compileAttributes()
# # load in library for testing
# devtools::load_all()

library(ggplot2)
library(gridExtra)

set.seed(123)

calc_C_phi = function(coords, phi) {
  n = nrow(coords)
  C_phi = matrix(0, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      # Calculate the distance between points i and j
      C_phi[i, j] = exp(-0.5 * t(coords[i,] - coords[j,]) %*% (coords[i,] - coords[j,]) / phi)
    }
  }
  return(C_phi)
}

set.seed(123)

lat = seq(0, 10, by = 1)
lon = seq(0, 10, by = 1)

coords = as.matrix(expand.grid(lat, lon))
colnames(coords) = c("lat", "lon")

n = nrow(coords)
p = 2

# generating w
sigmasq_w = 1
phi_w = 0.01
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

# generate knots
## 1 in every k^2 point on a grid is a knot
k = 3

lat_knots = unique(lat)
lat_knots = lat_knots[seq(1, length(lat_knots), by = k)]

lon_knots = unique(lon)
lon_knots = lon_knots[seq(1, length(lon_knots), by = k)]

knots = as.matrix(expand.grid(lat_knots, lon_knots))

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
  