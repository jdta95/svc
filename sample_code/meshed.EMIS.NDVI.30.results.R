library(data.table)
library(tidyverse)
library(meshed)
library(gridExtra)

# path where data is kept and saved
path_data = "C:/Users/jonat/OneDrive/Documents/Michele Peruzzi/Meshed JSS/data/"

# path for tables and figures
path_tabfig =
  "C:/Users/jonat/OneDrive/Documents/Michele Peruzzi/Meshed JSS/tables figures/"

load(paste0(path_data, "EMIS.NDVI.meshed.123.RData"))

mcmc_summary = function(x) {
  out = c(quantile(x, 0.025), mean(x), quantile(x, 0.975))
  names(out)[2] = "mean"
  return(out)
}

thin = seq(mcmc_thin, mcmc_keep * mcmc_thin, by = mcmc_thin)

# v chains
v1_chain = sapply(meshout$v_mcmc, function(x) {colSums(x)[1] / nrow(x)})
v2_chain = sapply(meshout$v_mcmc, function(x) {colSums(x)[2] / nrow(x)})

mcmc_summary(v1_chain)
mcmc_summary(v2_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = v1_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("v1 chain")

ggplot(data = data.frame(x = 1:mcmc_keep, y = v2_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("v2 chain")

# w chains
w1_chain = sapply(meshout$w_mcmc, function(x) {colSums(x)[1] / nrow(x)})
w2_chain = sapply(meshout$w_mcmc, function(x) {colSums(x)[2] / nrow(x)})

mcmc_summary(w1_chain)
mcmc_summary(w2_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = w1_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("w1 chain")

ggplot(data = data.frame(x = 1:mcmc_keep, y = w2_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("w2 chain")

# beta_0 (intercept) chains
b01_chain = meshout$beta_mcmc[1, 1, thin]
b02_chain = meshout$beta_mcmc[1, 2, thin]

mcmc_summary(b01_chain)
mcmc_summary(b02_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = b01_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("beta1 chain")

ggplot(data = data.frame(x = 1:mcmc_keep, y = b02_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("beta2 chain")

# w + beta_0 intercept chains
wb1_chain = w1_chain + b01_chain
wb2_chain = w2_chain + b02_chain

mcmc_summary(wb1_chain)
mcmc_summary(wb2_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = wb1_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("w1 + beta01 chain")

ggplot(data = data.frame(x = 1:mcmc_keep, y = wb2_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("w2 + beta02 chain")

# beta_1 (elevation coefficient) chains
b11_chain = meshout$beta_mcmc[2, 1, thin]
b12_chain = meshout$beta_mcmc[2, 2, thin]

mcmc_summary(b11_chain)
mcmc_summary(b12_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = b11_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("beta11 chain")

ggplot(data = data.frame(x = 1:mcmc_keep, y = b12_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("beta12 chain")

# tau squared
t1_chain = meshout$tausq_mcmc[1, thin]
t2_chain = meshout$tausq_mcmc[2, thin]

mcmc_summary(t1_chain)
mcmc_summary(t2_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = t1_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("tau squared 1 chain")

ggplot(data = data.frame(x = 1:mcmc_keep, y = t2_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("tau squared 2 chain")

# phi chain
phi1_chain = meshout$theta_mcmc[1, 1, thin]
phi2_chain = meshout$theta_mcmc[1, 2, thin]

mcmc_summary(phi1_chain)
mcmc_summary(phi2_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = phi1_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("phi 1 chain")

ggplot(data = data.frame(x = 1:mcmc_keep, y = phi2_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("phi 2 chain")

# lambda chains
lambda11_chain = meshout$lambda_mcmc[1, 1, thin]
lambda12_chain = meshout$lambda_mcmc[1, 2, thin]
lambda21_chain = meshout$lambda_mcmc[2, 1, thin]
lambda22_chain = meshout$lambda_mcmc[2, 2, thin]

mcmc_summary(lambda11_chain)
# mcmc_summary(lambda12_chain)
mcmc_summary(lambda21_chain)
mcmc_summary(lambda22_chain)

ggplot(data = data.frame(x = 1:mcmc_keep, y = lambda11_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("lambda 11 chain")

# ggplot(data = data.frame(x = 1:mcmc_keep, y = lambda12_chain),
#        mapping = aes(x = x, y = y)) +
#   geom_line() +
#   ggtitle("lambda 12 chain")

ggplot(data = data.frame(x = 1:mcmc_keep, y = lambda21_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("lambda 21 chain")

ggplot(data = data.frame(x = 1:mcmc_keep, y = lambda22_chain),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("lambda 22 chain")

# covariance matrix
## note distances range from 1.4 to 36
dists = seq(0, 2, by = 0.01)

cov_array = array(
  dim = c(2, 2, mcmc_keep, length(dists)),
  dimnames = list(NULL, NULL, paste0("chain ", 1:mcmc_keep), paste0("h = ", dists))
)

for (i in 1:length(dists)) {
  for (j in 1:length(thin)) {
    cov_array[, , j, i] = 
      meshout$lambda_mcmc[, , thin[j]] %*%
      diag(exp(-1 * meshout$theta_mcmc[1, , thin[j]] * dists[i])) %*%
      t(meshout$lambda_mcmc[, , thin[j]])
  }
}

cor_array = array(
  dim = c(2, 2, mcmc_keep, length(dists)),
  dimnames = list(NULL, NULL, paste0("chain ", 1:mcmc_keep), paste0("h = ", dists))
)

for (j in 1:length(thin)) {
  cov2cor_j =
    1 / sqrt(diag(meshout$lambda_mcmc[, , thin[j]] %*%
                    t(meshout$lambda_mcmc[, , thin[j]])))
  for (i in 1:length(dists)) {
    cor_array[, , j, i] =
      diag(cov2cor_j) %*% meshout$lambda_mcmc[, , thin[j]] %*%
      diag(exp(-1 * meshout$theta_mcmc[1, , thin[j]] * dists[i])) %*%
      t(meshout$lambda_mcmc[, , thin[j]]) %*% diag(cov2cor_j)
  }
}

mcmc_summary(cor_array[1, 2, , "h = 0"])

mcmc_summary(cov_array[1, 1, , "h = 0.01"])
mcmc_summary(cov_array[1, 2, , "h = 0.01"])
mcmc_summary(cov_array[2, 2, , "h = 0.01"])

## emis variance chain at h = 0
ggplot(data = data.frame(x = 1:mcmc_keep, y = cov_array[1, 1, , "h = 0"]),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("var 1 chain at h = 0")

## emis-ndvi covariance chain at h = 0
ggplot(data = data.frame(x = 1:mcmc_keep, y = cov_array[1, 2, , "h = 0"]),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("covar chain at h = 0")

## ndvi variance chain at h = 0
ggplot(data = data.frame(x = 1:mcmc_keep, y = cov_array[2, 2, , "h = 0"]),
       mapping = aes(x = x, y = y)) +
  geom_line() +
  ggtitle("var 2 chain at h = 0")

# predictions
y_post_sample = meshout$yhat_mcmc %>%
  abind::abind(along = 3) %>%
  apply(1:2, mcmc_summary)

# indices for EMIS observations in testing set
pred.ind.emis = dt[, which(is.na(emis) & !is.na(emis.obs))]

# mean square error for EMIS predictions
MSE.emis = mean((dt[pred.ind.emis, emis.obs] - y_post_sample[2, pred.ind.emis, 1]) ^ 2)

MSE.emis

# indices for NDVI observations in testing set
pred.ind.ndvi = dt[, which(is.na(ndvi) & !is.na(ndvi.obs))]

# mean square error for NDVI predictions
MSE.ndvi = mean((dt[pred.ind.ndvi, ndvi.obs] - y_post_sample[2, pred.ind.ndvi, 2]) ^ 2)

MSE.ndvi

e1 = meshout$coordsdata %>%
  cbind(post_EMIS = y_post_sample[2, , 1]) %>%
  filter(forced_grid == 0) %>%
  ggplot(aes(lon, lat)) +
  geom_tile(aes(fill = post_EMIS)) +
  xlab("longitude") +
  ylab("latitude") +
  guides(fill = guide_legend(title = "predicted emissivity")) +
  coord_fixed()

e2 = ggplot(dt, aes(lon, lat)) +
  geom_tile(aes(fill = emis.obs)) +
  xlab("longitude") +
  ylab("latitude") +
  guides(fill = guide_legend(title = "observed emissivity")) +
  coord_fixed()

e3 = ggplot(dt, aes(lon, lat)) +
  geom_tile(aes(fill = emis)) +
  xlab("longitude") +
  ylab("latitude") +
  guides(fill = guide_legend(title = "training emissivity")) +
  coord_fixed()

n1 = meshout$coordsdata %>%
  cbind(post_NDVI = y_post_sample[2, , 2]) %>%
  filter(forced_grid == 0) %>%
  ggplot(aes(lon, lat)) +
  geom_tile(aes(fill = post_NDVI)) +
  xlab("longitude") +
  ylab("latitude") +
  guides(fill = guide_legend(title = "predicted NDVI")) +
  coord_fixed()

n2 = ggplot(dt, aes(lon, lat)) +
  geom_tile(aes(fill = ndvi.obs)) +
  xlab("longitude") +
  ylab("latitude") +
  guides(fill = guide_legend(title = "observed NDVI")) +
  coord_fixed()

n3 = ggplot(dt, aes(lon, lat)) +
  geom_tile(aes(fill = ndvi)) +
  xlab("longitude") +
  ylab("latitude") +
  guides(fill = guide_legend(title = "training NDVI")) +
  coord_fixed()

png(
  filename = paste0(path_tabfig, "EMIS.NDVI.pred.map.png"),
  width = 1620,
  height = 1080)

grid.arrange(e3, e2, e1, n3, n2, n1, nrow = 2, ncol = 3)

dev.off()

