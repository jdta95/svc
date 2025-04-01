library(data.table)
library(tidyverse)
library(meshed)
library(gridExtra)

# path where data is kept and saved
path_data = "data/"

# path for tables and figures
path_tabfig =
  "tables figures/"

load(paste0(path_data, "TEMP.meshed.123.RData"))

mcmc_summary = function(x) {
  out = c(quantile(x, 0.025), mean(x), quantile(x, 0.975))
  names(out)[2] = "mean"
  return(out)
}

thin = seq(mcmc_thin, mcmc_keep * mcmc_thin, by = mcmc_thin)

# v chain
v1_chain = sapply(meshout$v_mcmc, mean)

mcmc_summary(v1_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = v1_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("v chain")

# w chain
w1_chain = sapply(meshout$v_mcmc, mean) * meshout$lambda_mcmc[, , thin]

mcmc_summary(w1_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = w1_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("w chain")

# beta_0 (intercept) chain
b01_chain = meshout$beta_mcmc[1, 1, thin]

mcmc_summary(b01_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = b01_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("beta0 chain")

# w + beta_0 intercept chain
wb1_chain = w1_chain + b01_chain

mcmc_summary(wb1_chain)

a = ggplot(data = data.frame(x = 1:mcmc_keep, y = wb1_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("w + beta0 chain")

# beta_1 (elevation coefficient) chain
b11_chain = meshout$beta_mcmc[2, 1, thin]

mcmc_summary(b11_chain)

b = ggplot(data = data.frame(x = 1:mcmc_keep, y = b11_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("beta1 chain")

grid.arrange(a, b)

# tau squared
t1_chain = meshout$tausq_mcmc[1, thin]

mcmc_summary(t1_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = t1_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("tau squared chain")

# phi chain
phi1_chain = meshout$theta_mcmc[1, 1, thin]

mcmc_summary(phi1_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = phi1_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("phi chain")

# lambda chains
lambda11_chain = meshout$lambda_mcmc[1, 1, thin]

mcmc_summary(lambda11_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = lambda11_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("lambda chain")

# covariance matrix
## note distances range from 0 to 2
dists = seq(0, 2, by = 0.01)

cov_array = array(
  dim = c(mcmc_keep, length(dists)),
  dimnames = list(paste0("chain ", 1:mcmc_keep), paste0("h = ", dists))
)

for (i in 1:length(dists)) {
  for (j in 1:length(thin)) {
    cov_array[j, i] = 
      meshout$lambda_mcmc[, , thin[j]] *
      exp(-1 * meshout$theta_mcmc[1, , thin[j]] * dists[i]) *
      meshout$lambda_mcmc[, , thin[j]]
  }
}

cor_array = array(
  dim = c(mcmc_keep, length(dists)),
  dimnames = list(paste0("chain ", 1:mcmc_keep), paste0("h = ", dists))
)

for (j in 1:length(thin)) {
  cov2cor_j =
    1 / meshout$lambda_mcmc[, , thin[j]]
  for (i in 1:length(dists)) {
    cor_array[j, i] =
      cov2cor_j ^ 2 * meshout$lambda_mcmc[, , thin[j]] ^ 2 * 
      exp(-1 * meshout$theta_mcmc[1, , thin[j]] * dists[i])
  }
}

mcmc_summary(cor_array[ , "h = 0.01"])
mcmc_summary(cor_array[ , "h = 0.05"])

mcmc_summary(cov_array[ , "h = 0"])
mcmc_summary(cov_array[ , "h = 0.01"])
mcmc_summary(cov_array[ , "h = 0.05"])

## emis variance chain at h = 0
ggplot(data = data.frame(x = 1:mcmc_keep, y = cov_array[, "h = 0"]),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("var chain at h = 0")

# predictions
y_post_sample = meshout$yhat_mcmc %>%
  abind::abind(along = 3) %>%
  apply(1, mcmc_summary)

# indices for observations in testing set
pred.ind = dt[, which(is.na(temp) & !is.na(temp.obs))]

# mean square error for predictions
MSE = mean((dt[pred.ind, temp.obs] - y_post_sample[2, pred.ind]) ^ 2)

MSE

# potential colorscales to use
# scico::scale_fill_scico(palette="vik") + # bam, broc, cork, managua, vik

p1 = meshout$coordsdata %>%
  cbind(post_TEMP = y_post_sample[2, ]) %>%
  filter(forced_grid == 0) %>%
  ggplot(aes(lon, lat)) +
  geom_tile(aes(fill = post_TEMP)) +
  xlab("longitude") +
  ylab("latitude") +
  guides(fill = guide_legend(title = "predicted temperature")) +
  coord_fixed()

p2 = ggplot(dt, aes(lon, lat)) +
  geom_tile(aes(fill = temp.obs)) +
  xlab("longitude") +
  ylab("latitude") +
  guides(fill = guide_legend(title = "observed temperature")) +
  coord_fixed()

p3 = ggplot(dt, aes(lon, lat)) +
  geom_tile(aes(fill = temp)) +
  xlab("longitude") +
  ylab("latitude") +
  guides(fill = guide_legend(title = "training temperature")) +
  coord_fixed()

png(
  filename = paste0(path_tabfig, "TEMP.pred.map.png"),
  width = 1620,
  height = 540)

grid.arrange(p3, p2, p1, nrow = 1, ncol = 3)

dev.off()

