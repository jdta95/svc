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
X = cbind(X_0, X_1, X_2)

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

save.image("data/svclm_test_data.RData")