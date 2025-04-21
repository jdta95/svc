# Data prep
library(ggplot2)
library(gridExtra)

load("data/svclm_test_data.RData")

# show all plots in a grid
grid.arrange(Y_plot, w_0_plot, w_1_plot, w_2_plot, epsilon_plot, ncol = 2)

# # compile package
# Rcpp::compileAttributes()
# # load in library for testing
# devtools::load_all()

# generate knots
## 1 in every k^2 point on a grid is a knot
## ensures knot data does not become mismatched
knots_list = svc::simpleknots(
  Y = Y,
  X = as.matrix(data.frame(X_0, X_1, X_2)),
  coords = coords,
  k = 2
  )

Y_knots = knots_list$Y_knots
X_knots = knots_list$X_knots
knots = knots_list$knots
knots_df = knots_list$knots_data_frame

# Y_knots_plot = ggplot(data = knots_df) +
#   geom_tile(aes(x = lon, y = lat, fill = Y)) +
#   guides(fill = guide_legend(title = "Y knots")) +
#   coord_fixed() + 
#   scico::scale_fill_scico(palette="bam") +
#   labs(x = "longitude", y = "latitude")
# 
# grid.arrange(Y_plot, Y_knots_plot, ncol = 2)

l = 0
u = 30

# Test function
mcmc = 3000
p = ncol(X)

# svclm with minimum required arguments
# uses full GP by default, not feasible for large data sets
full_time = system.time({
  output_full = svc::svclm(
    Y = Y,
    X = X,
    coords = coords,
    phi_lower = rep(l, p),
    phi_upper = rep(u, p),
    mcmc = mcmc
  )
})

# save(output_full, full_time, file = "data/fullrank_output.RData")
load("data/fullrank_output.RData")

# svclm with recommended arguments using knot data from simpleknots()
LR_time = system.time({
  output_LR = svc::svclm(
    Y = Y,
    X = X,
    coords = coords,
    Y_knots = Y_knots,
    X_knots = X_knots,
    knots = knots,
    phi_lower = rep(l, p),
    phi_upper = rep(u, p),
    mcmc = mcmc
  )
})

# save(output_LR, LR_time, file = "data/LR_output.RData")
load("data/LR_output.RData")

full_time
LR_time

my_summary = function(x){
  return (c(
    quantile(x, 0.025),
    mean(x),
    quantile(x, 0.975)
  ))
}

end = mcmc

# all plots
plot(y = output_LR$phi_samples[1:end, 1], x = 1:end, type = "l", main = "Trace plot of phi_0 LR", ylab = "phi_0", xlab = "Iteration")

plot(y = output_full$phi_samples[1:end, 1], x = 1:end, type = "l", main = "Trace plot of phi_0 full", ylab = "phi_0", xlab = "Iteration")

plot(y = output_LR$phi_samples[1:end, 2], x = 1:end, type = "l", main = "Trace plot of phi_1 LR", ylab = "phi_1", xlab = "Iteration")

plot(y = output_full$phi_samples[1:end, 2], x = 1:end, type = "l", main = "Trace plot of phi_1 full", ylab = "phi_1", xlab = "Iteration")

plot(y = output_LR$phi_samples[1:end, 3], x = 1:end, type = "l", main = "Trace plot of phi_2 LR", ylab = "phi_2", xlab = "Iteration")

plot(y = output_full$phi_samples[1:end, 3], x = 1:end, type = "l", main = "Trace plot of phi_2 full", ylab = "phi_2", xlab = "Iteration")

plot(y = output_LR$sigmasq_samples[1:end, 1], x = 1:end, type = "l", main = "Trace plot of sigmasq_0 LR", ylab = "sigmasq_0", xlab = "Iteration")

plot(y = output_full$sigmasq_samples[1:end, 1], x = 1:end, type = "l", main = "Trace plot of sigmasq_0 full", ylab = "sigmasq_0", xlab = "Iteration")

plot(y = output_LR$sigmasq_samples[1:end, 2], x = 1:end, type = "l", main = "Trace plot of sigmasq_1 LR", ylab = "sigmasq_1", xlab = "Iteration")

plot(y = output_full$sigmasq_samples[1:end, 2], x = 1:end, type = "l", main = "Trace plot of sigmasq_1 full", ylab = "sigmasq_1", xlab = "Iteration")

plot(y = output_LR$sigmasq_samples[1:end, 3], x = 1:end, type = "l", main = "Trace plot of sigmasq_2 LR", ylab = "sigmasq_2", xlab = "Iteration")

plot(y = output_full$sigmasq_samples[1:end, 3], x = 1:end, type = "l", main = "Trace plot of sigmasq_2 full", ylab = "sigmasq_2", xlab = "Iteration")

plot(y = output_LR$tausq_samples[1:end], x = 1:end, type = "l", main = "Trace plot of tau^2 LR", ylab = "tau^2", xlab = "Iteration")

plot(y = output_full$tausq_samples[1:end], x = 1:end, type = "l", main = "Trace plot of tau^2 full", ylab = "tau^2", xlab = "Iteration")

start = 2000

# phi_0 LR
mean(output_LR$phi_acceptance[, 1])
my_summary(output_LR$phi_samples[start:end, 1])
phi_0

# phi_0 full
mean(output_full$phi_acceptance[, 1])
my_summary(output_full$phi_samples[start:end, 1])
phi_0

# phi_1 LR
mean(output_LR$phi_acceptance[, 2])
my_summary(output_LR$phi_samples[start:end, 2])
phi_1

# phi_1 full
mean(output_full$phi_acceptance[, 2])
my_summary(output_full$phi_samples[start:end, 2])
phi_1

# phi_2 LR
mean(output_LR$phi_acceptance[, 3])
my_summary(output_LR$phi_samples[start:end, 3])
phi_2

# phi_2 full
mean(output_full$phi_acceptance[, 3])
my_summary(output_full$phi_samples[start:end, 3])
phi_2

# sigma_0 LR
my_summary(output_LR$sigmasq_samples[start:end, 1])
sigmasq_0

# sigma_0 full
my_summary(output_full$sigmasq_samples[start:end, 1])
sigmasq_0

# sigma_1 LR
my_summary(output_LR$sigmasq_samples[start:end, 2])
sigmasq_1

# sigma_1 full
my_summary(output_full$sigmasq_samples[start:end, 2])
sigmasq_1

# sigma_2 LR
my_summary(output_LR$sigmasq_samples[start:end, 3])
sigmasq_2

# sigma_2 full
my_summary(output_full$sigmasq_samples[start:end, 3])
sigmasq_2

# tau^2 LR
my_summary(output_LR$tausq_samples[start:end])
tausq

# tau^2 full
my_summary(output_full$tausq_samples[start:end])
tausq

# w results
unknowns_df = data.frame(coords, w_0, w_1, w_2, epsilon)
knot_unknowns_df = data.frame(knots)
colnames(knot_unknowns_df) = c("lat", "lon")
knot_unknowns_df = merge(knot_unknowns_df, unknowns_df, by = c("lat", "lon"), all.x = TRUE, all.y = FALSE)

w_0_knots_plot = ggplot(data = knot_unknowns_df) +
  geom_tile(aes(x = lon, y = lat, fill = w_0)) +
  guides(fill = guide_legend(title = "w_0 knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

w_1_knots_plot = ggplot(data = knot_unknowns_df) +
  geom_tile(aes(x = lon, y = lat, fill = w_1)) +
  guides(fill = guide_legend(title = "w_1 knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

w_2_knots_plot = ggplot(data = knot_unknowns_df) +
  geom_tile(aes(x = lon, y = lat, fill = w_2)) +
  guides(fill = guide_legend(title = "w_2 knots")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

# w_0 LR
w_0_LR_pred_plot = ggplot(data = data.frame(coords, what0 = colMeans(output_LR$w_samples[start:end, , 1]))) +
  geom_tile(aes(x = lon, y = lat, fill = what0)) +
  guides(fill = guide_legend(title = "w_0 hats LR")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

# w_1 LR
w_1_LR_pred_plot = ggplot(data = data.frame(coords, what1 = colMeans(output_LR$w_samples[start:end, , 2]))) +
  geom_tile(aes(x = lon, y = lat, fill = what1)) +
  guides(fill = guide_legend(title = "w_1 hats LR")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

# w_2 LR
w_2_LR_pred_plot = ggplot(data = data.frame(coords, what2 = colMeans(output_LR$w_samples[start:end, , 3]))) +
  geom_tile(aes(x = lon, y = lat, fill = what2)) +
  guides(fill = guide_legend(title = "w_2 hats LR")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

# w_0 full
w_0_full_pred_plot = ggplot(data = data.frame(coords, what0 = colMeans(output_full$w_samples[start:end, , 1]))) +
  geom_tile(aes(x = lon, y = lat, fill = what0)) +
  guides(fill = guide_legend(title = "w_0 hats full")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

# w_1 full
w_1_full_pred_plot = ggplot(data = data.frame(coords, what1 = colMeans(output_full$w_samples[start:end, , 2]))) +
  geom_tile(aes(x = lon, y = lat, fill = what1)) +
  guides(fill = guide_legend(title = "w_1 hats full")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

# w_2 full
w_2_full_pred_plot = ggplot(data = data.frame(coords, what2 = colMeans(output_full$w_samples[start:end, , 3]))) +
  geom_tile(aes(x = lon, y = lat, fill = what2)) +
  guides(fill = guide_legend(title = "w_2 hats full")) +
  coord_fixed() + 
  scico::scale_fill_scico(palette="bam") +
  labs(x = "longitude", y = "latitude")

grid.arrange(w_0_plot, w_0_full_pred_plot, w_0_knots_plot, w_0_LR_pred_plot, ncol = 2)
grid.arrange(w_1_plot, w_1_full_pred_plot, w_1_knots_plot, w_1_LR_pred_plot, ncol = 2)
grid.arrange(w_2_plot, w_2_full_pred_plot, w_2_knots_plot, w_2_LR_pred_plot, ncol = 2)

# # arbitrary w_0(s)
# s = n / 2
# plot(y = output_LR$w_samples[1:end, s, 1], x = 1:end, type = "l", main = "Trace plot of w_0", ylab = "w_0", xlab = "Iteration")
# my_summary(output_LR$w_samples[start:end, s, 1])
# w_0[s]

