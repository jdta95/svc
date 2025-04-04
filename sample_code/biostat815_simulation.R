library(varycoef)
library(ggplot2)
library(gridExtra)
library(Matrix)
set.seed(123)

calc_C_phi = function(coords, phi) {
  n = nrow(coords)
  dists = as.matrix(dist(coords)) ^ 2
  C_phi = exp(-0.5 / phi * dists)
}



### Purpose
# The purpose of this simulation study is compare the performance between our method
# and varycoef R package. Specfically, we will assess the bias, computational speed, average length of credible interval
# and coverage

## Setting
# 50 replications
# mcmc = 1000 for gibbs
# 100 or about 500 sample size
# Vary by number of betas (2 or 3) and priors (low or high)


## metric
### Bias ###
### 1/R (\sum_{r=1} (beta hat - true beta ))  


# How long to run in seconds

### Coverage ### - later
### Coverage rate of the nominal 95% crtedible interval 

### 95% Credible Interval ### - later



mcmc = 1000

### Scenario 1: (2 betas and "low" priors)
### 100 Sample size
gibbs1 = data.frame()
coef1 = data.frame()
for (i in 1:50) {
  
  lat = seq(0, 9, by = 1)
  lon = seq(0, 9, by = 1)
  
  coords = as.matrix(expand.grid(lat, lon))
  colnames(coords) = c("lat", "lon")
  
  n = nrow(coords)
  p = 2
  
  # generating w
  sigmasq_w = 1
  phi_w = 2
  
  C_w = calc_C_phi(coords, phi_w)
  w = MASS::mvrnorm(1, rep(0, n), nearPD(sigmasq_w * C_w,corr = FALSE,keepDiag = FALSE,eig.tol = 1e-6)$mat)
  
  # generating beta's
  sigmasq_1 = 1
  phi_1 = 2
  
  C_1 = calc_C_phi(coords, phi_1)
  beta_1_mean = 3
  beta_1 = MASS::mvrnorm(1, rep(beta_1_mean, n),nearPD(sigmasq_1 * C_1,corr = FALSE,keepDiag = FALSE,eig.tol = 1e-6)$mat )
  
  sigmasq_2 = 1
  phi_2 = 2
  
  C_2 = calc_C_phi(coords, phi_2)
  beta_2_mean = -5
  beta_2 = MASS::mvrnorm(1, rep(beta_2_mean, n),nearPD(sigmasq_2 * C_2,corr = FALSE,keepDiag = FALSE,eig.tol = 1e-6)$mat )
  
  # generating X
  X_1 = rnorm(n, 0, 10)
  X_2 = rnorm(n, 0, 10)
  
  # generating epsilon
  tausq = 0.0001
  epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))
  
  Y = X_1 * beta_1 + X_2 * beta_2 + w + epsilon
  
  k = 2
  
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
  u = 4
  start_time <- Sys.time()
  output = GP_Gibbs_beta(Y = Y,
                X = X,
                s = coords,
                knots = knots,
                beta_knots_start = matrix(0, nrow = m, ncol = p),
                w_knots_start = rep(0, m),
                phi_beta_start = rep(l + (u - l) / 2, p),
                phi_w_start = l + (u - l) / 2,
                sigmasq_beta_start = rep(1, p),
                sigmasq_w_start = 1,
                tausq_start = 1,
                phi_beta_proposal_sd = rep(0.1, p),
                phi_w_proposal_sd = 0.1,
                a_beta = rep(2, p),
                b_beta = rep(1, p),
                a_w = 2,
                b_w = 1,
                a_t = 2,
                b_t = 1e-4,
                lower_beta = rep(l, p),
                upper_beta = rep(u, p),
                lower_w = l,
                upper_w = u,
                mcmc = mcmc)
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  gibbs1 = rbind(gibbs1,c(mean(colMeans(output$beta_samples[, 1, 800:1000]) -3),mean(colMeans(output$beta_samples[, 2, 800:1000]) +5), time_taken ))
  start_time <- Sys.time()
  fit <- SVC_mle(Y, X, coords,as.matrix(w))
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  coef1 = rbind(coef1,c(fit$coefficients[1]-3,fit$coefficients[2]+5,time_taken))
  }


gibbs2 = data.frame()
coef2 = data.frame()
### 361 Sample size
for (i in 1:50) {
  lat = seq(0, 9, by = .4)
  lon = seq(0, 9, by = .4)
  
  coords = as.matrix(expand.grid(lat, lon))
  colnames(coords) = c("lat", "lon")
  
  n = nrow(coords)
  p = 2
  
  # generating w
  sigmasq_w = 1
  phi_w = 2
  C_w = calc_C_phi(coords, phi_w)
  w = MASS::mvrnorm(1, rep(0, n), sigmasq_w * C_w)
  
  # generating beta's
  sigmasq_1 = 1
  phi_1 = 2
  C_1 = calc_C_phi(coords, phi_1)
  beta_1_mean = 3
  beta_1 = MASS::mvrnorm(1, rep(beta_1_mean, n), sigmasq_1 * C_1)
  
  sigmasq_2 = 1
  phi_2 = 2
  C_2 = calc_C_phi(coords, phi_2)
  beta_2_mean = -5
  beta_2 = MASS::mvrnorm(1, rep(beta_2_mean, n), sigmasq_2 * C_2)
  
  # generating X
  X_1 = rnorm(n, 0, 10)
  X_2 = rnorm(n, 0, 10)
  
  # generating epsilon
  tausq = 0.0001
  epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))
  
  Y = X_1 * beta_1 + X_2 * beta_2 + w + epsilon
  
  k = 2
  
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
  u = 4
  
  start_time <- Sys.time()
  output = GP_Gibbs_beta(Y = Y,
                         X = X,
                         s = coords,
                         knots = knots,
                         beta_knots_start = matrix(0, nrow = m, ncol = p),
                         w_knots_start = rep(0, m),
                         phi_beta_start = rep(l + (u - l) / 2, p),
                         phi_w_start = l + (u - l) / 2,
                         sigmasq_beta_start = rep(1, p),
                         sigmasq_w_start = 1,
                         tausq_start = 1,
                         phi_beta_proposal_sd = rep(0.1, p),
                         phi_w_proposal_sd = 0.1,
                         a_beta = rep(2, p),
                         b_beta = rep(1, p),
                         a_w = 2,
                         b_w = 1,
                         a_t = 2,
                         b_t = 1e-4,
                         lower_beta = rep(l, p),
                         upper_beta = rep(u, p),
                         lower_w = l,
                         upper_w = u,
                         mcmc = mcmc)
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  gibbs2 = rbind(gibbs2,c(mean(colMeans(output$beta_samples[, 1, 800:1000]) -3),mean(colMeans(output$beta_samples[, 2, 800:1000]) +5), time_taken ))
  start_time <- Sys.time()
  fit <- SVC_mle(Y, X, coords,as.matrix(w))
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  coef2 = rbind(coef2,c(fit$coefficients[1]-3,fit$coefficients[2]+5,time_taken))
}


gibbs3 = data.frame()
coef3 = data.frame()
### Scenario 2: (2 betas and "high" priors)
### 100 Sample size
for (i in 1:50) {
  
  lat = seq(0, 9, by = 1)
  lon = seq(0, 9, by = 1)
  
  coords = as.matrix(expand.grid(lat, lon))
  colnames(coords) = c("lat", "lon")
  
  n = nrow(coords)
  p = 2
  
  # generating w
  sigmasq_w = 4
  phi_w = 8
   
  C_w = calc_C_phi(coords, phi_w)
  w = MASS::mvrnorm(1, rep(0, n), sigmasq_w * C_w)
  
  # generating beta's
  sigmasq_1 = 4
  phi_1 = 8
   
  C_1 = calc_C_phi(coords, phi_1)
  beta_1_mean = 3
  beta_1 = MASS::mvrnorm(1, rep(beta_1_mean, n), sigmasq_1 * C_1)
  
  sigmasq_2 = 4
  phi_2 = 8
   
  C_2 = calc_C_phi(coords, phi_2)
  beta_2_mean = -5
  beta_2 = MASS::mvrnorm(1, rep(beta_2_mean, n), sigmasq_2 * C_2)
  
  # generating X
  X_1 = rnorm(n, 0, 10)
  X_2 = rnorm(n, 0, 10)
  
  # generating epsilon
  tausq = 0.01
  epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))
  
  Y = X_1 * beta_1 + X_2 * beta_2 + w + epsilon
  
  k = 2
  
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
  u = 4
  
  start_time <- Sys.time()
  output = GP_Gibbs_beta(Y = Y,
                         X = X,
                         s = coords,
                         knots = knots,
                         beta_knots_start = matrix(0, nrow = m, ncol = p),
                         w_knots_start = rep(0, m),
                         phi_beta_start = rep(l + (u - l) / 2, p),
                         phi_w_start = l + (u - l) / 2,
                         sigmasq_beta_start = rep(1, p),
                         sigmasq_w_start = 1,
                         tausq_start = 1,
                         phi_beta_proposal_sd = rep(0.1, p),
                         phi_w_proposal_sd = 0.1,
                         a_beta = rep(2, p),
                         b_beta = rep(1, p),
                         a_w = 2,
                         b_w = 1,
                         a_t = 2,
                         b_t = 1e-4,
                         lower_beta = rep(l, p),
                         upper_beta = rep(u, p),
                         lower_w = l,
                         upper_w = u,
                         mcmc = mcmc)
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  gibbs3 = rbind(gibbs3,c(mean(colMeans(output$beta_samples[, 1, 800:1000]) -3),mean(colMeans(output$beta_samples[, 2, 800:1000]) +5), time_taken ))
  start_time <- Sys.time()
  fit <- SVC_mle(Y, X, coords,as.matrix(w))
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  coef3 = rbind(coef3,c(fit$coefficients[1]-3,fit$coefficients[2]+5,time_taken))
}

gibbs4 = data.frame()
coef4 = data.frame()
### 361 Sample size
for (i in 1:50) {
  lat = seq(0, 9, by = .4)
  lon = seq(0, 9, by = .4)
  
  coords = as.matrix(expand.grid(lat, lon))
  colnames(coords) = c("lat", "lon")
  
  n = nrow(coords)
  p = 2
  
  # generating w
  sigmasq_w = 4
  phi_w = 8
  C_w = calc_C_phi(coords, phi_w)
  w = MASS::mvrnorm(1, rep(0, n), sigmasq_w * C_w)
  
  # generating beta's
  sigmasq_1 = 4
  phi_1 = 8
  C_1 = calc_C_phi(coords, phi_1)
  beta_1_mean = 3
  beta_1 = MASS::mvrnorm(1, rep(beta_1_mean, n), sigmasq_1 * C_1)
  
  sigmasq_2 = 4
  phi_2 = 8
  C_2 = calc_C_phi(coords, phi_2)
  beta_2_mean = -5
  beta_2 = MASS::mvrnorm(1, rep(beta_2_mean, n), sigmasq_2 * C_2)
  
  # generating X
  X_1 = rnorm(n, 0, 10)
  X_2 = rnorm(n, 0, 10)
  
  # generating epsilon
  tausq = 0.01
  epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))
  
  Y = X_1 * beta_1 + X_2 * beta_2 + w + epsilon
  
  k = 2
  
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
  u = 4
  
  start_time <- Sys.time()
  output = GP_Gibbs_beta(Y = Y,
                         X = X,
                         s = coords,
                         knots = knots,
                         beta_knots_start = matrix(0, nrow = m, ncol = p),
                         w_knots_start = rep(0, m),
                         phi_beta_start = rep(l + (u - l) / 2, p),
                         phi_w_start = l + (u - l) / 2,
                         sigmasq_beta_start = rep(1, p),
                         sigmasq_w_start = 1,
                         tausq_start = 1,
                         phi_beta_proposal_sd = rep(0.1, p),
                         phi_w_proposal_sd = 0.1,
                         a_beta = rep(2, p),
                         b_beta = rep(1, p),
                         a_w = 2,
                         b_w = 1,
                         a_t = 2,
                         b_t = 1e-4,
                         lower_beta = rep(l, p),
                         upper_beta = rep(u, p),
                         lower_w = l,
                         upper_w = u,
                         mcmc = mcmc)
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  gibbs4 = rbind(gibbs4,c(mean(colMeans(output$beta_samples[, 1, 800:1000]) -3),mean(colMeans(output$beta_samples[, 2, 800:1000]) +5), time_taken ))
  start_time <- Sys.time()
  fit <- SVC_mle(Y, X, coords,as.matrix(w))
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  coef4 = rbind(coef4,c(fit$coefficients[1]-3,fit$coefficients[2]+5,time_taken))
}















gibbs5 = data.frame()
coef5 = data.frame()
### Scenario 3: (3 betas and "low" priors)
### 100 Sample size
for (i in 1:50) {
  
  lat = seq(0, 9, by = 1)
  lon = seq(0, 9, by = 1)
  
  coords = as.matrix(expand.grid(lat, lon))
  colnames(coords) = c("lat", "lon")
  
  n = nrow(coords)
  p = 2
  
  # generating w
  sigmasq_w = 1
  phi_w = 2
   
  C_w = calc_C_phi(coords, phi_w)
  w = MASS::mvrnorm(1, rep(0, n), sigmasq_w * C_w)
  
  # generating beta's
  sigmasq_1 = 1
  phi_1 = 2
   
  C_1 = calc_C_phi(coords, phi_1)
  beta_1_mean = 3
  beta_1 = MASS::mvrnorm(1, rep(beta_1_mean, n), sigmasq_1 * C_1)
  
  sigmasq_2 = 1
  phi_2 = 2
   
  C_2 = calc_C_phi(coords, phi_2)
  beta_2_mean = -5
  beta_2 = MASS::mvrnorm(1, rep(beta_2_mean, n), sigmasq_2 * C_2)
  
  sigmasq_3 = 1
  phi_3 = 2
   
  C_3 = calc_C_phi(coords, phi_3)
  beta_3_mean = -10
  beta_3 = MASS::mvrnorm(1, rep(beta_3_mean, n), sigmasq_3 * C_3)
  
  # generating X
  X_1 = rnorm(n, 0, 10)
  X_2 = rnorm(n, 0, 10)
  X_3 = rnorm(n, 0, 10)
  
  # generating epsilon
  tausq = 0.0001
  epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))
  
  Y = X_1 * beta_1 + X_2 * beta_2 + w + epsilon
  
  k = 2
  
  lat_knots = unique(lat)
  lat_knots = lat_knots[seq(1, length(lat_knots), by = k)]
  
  lon_knots = unique(lon)
  lon_knots = lon_knots[seq(1, length(lon_knots), by = k)]
  
  knots = as.matrix(expand.grid(lat_knots, lon_knots))
  
  df = data.frame(coords, Y, X_1, X_2,X_3, w, beta_1, beta_2, epsilon)
  
  knots_df = data.frame(knots)
  colnames(knots_df) = c("lat", "lon")
  
  knots_df = merge(knots_df, df, by = c("lat", "lon"), all.x = TRUE, all.y = FALSE)
  
  X = as.matrix(df[, c("X_1", "X_2","X_3")])
  Y_star = knots_df$Y
  X_star = as.matrix(knots_df[, c("X_1", "X_2","X_3")])
  
  m = nrow(knots)
  l = 0
  u = 4
  
  start_time <- Sys.time()
  output = GP_Gibbs_beta(Y = Y,
                         X = X,
                         s = coords,
                         knots = knots,
                         beta_knots_start = matrix(0, nrow = m, ncol = p),
                         w_knots_start = rep(0, m),
                         phi_beta_start = rep(l + (u - l) / 2, p),
                         phi_w_start = l + (u - l) / 2,
                         sigmasq_beta_start = rep(1, p),
                         sigmasq_w_start = 1,
                         tausq_start = 1,
                         phi_beta_proposal_sd = rep(0.1, p),
                         phi_w_proposal_sd = 0.1,
                         a_beta = rep(2, p),
                         b_beta = rep(1, p),
                         a_w = 2,
                         b_w = 1,
                         a_t = 2,
                         b_t = 1e-4,
                         lower_beta = rep(l, p),
                         upper_beta = rep(u, p),
                         lower_w = l,
                         upper_w = u,
                         mcmc = mcmc)
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  gibbs5 = rbind(gibbs5,c(mean(colMeans(output$beta_samples[, 1, 800:1000]) -3),mean(colMeans(output$beta_samples[, 2, 800:1000]) +5), time_taken ))
  start_time <- Sys.time()
  fit <- SVC_mle(Y, X, coords,as.matrix(w))
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  coef5 = rbind(coef5,c(fit$coefficients[1]-3,fit$coefficients[2]+5,time_taken))
}

gibbs6 = data.frame()
coef6 = data.frame()
### 361 Sample size
for (i in 1:50) {
  lat = seq(0, 9, by = .4)
  lon = seq(0, 9, by = .4)
  
  coords = as.matrix(expand.grid(lat, lon))
  colnames(coords) = c("lat", "lon")
  
  n = nrow(coords)
  p = 2
  
  # generating w
  sigmasq_w = 1
  phi_w = 2
   
  C_w = calc_C_phi(coords, phi_w)
  w = MASS::mvrnorm(1, rep(0, n), sigmasq_w * C_w)
  
  # generating beta's
  sigmasq_1 = 1
  phi_1 = 2
   
  C_1 = calc_C_phi(coords, phi_1)
  beta_1_mean = 3
  beta_1 = MASS::mvrnorm(1, rep(beta_1_mean, n), sigmasq_1 * C_1)
  
  sigmasq_2 = 1
  phi_2 = 2
   
  C_2 = calc_C_phi(coords, phi_2)
  beta_2_mean = -5
  beta_2 = MASS::mvrnorm(1, rep(beta_2_mean, n), sigmasq_2 * C_2)
  
  sigmasq_3 = 1
  phi_3 = 2
   
  C_3 = calc_C_phi(coords, phi_3)
  beta_3_mean = -10
  beta_3 = MASS::mvrnorm(1, rep(beta_3_mean, n), sigmasq_3 * C_3)
  
  # generating X
  X_1 = rnorm(n, 0, 10)
  X_2 = rnorm(n, 0, 10)
  X_3 = rnorm(n, 0, 10)
  
  # generating epsilon
  tausq = 0.0001
  epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))
  
  Y = X_1 * beta_1 + X_2 * beta_2 + w + epsilon
  
  k = 2
  
  lat_knots = unique(lat)
  lat_knots = lat_knots[seq(1, length(lat_knots), by = k)]
  
  lon_knots = unique(lon)
  lon_knots = lon_knots[seq(1, length(lon_knots), by = k)]
  
  knots = as.matrix(expand.grid(lat_knots, lon_knots))
  
  df = data.frame(coords, Y, X_1, X_2,X_3, w, beta_1, beta_2, epsilon)
  
  knots_df = data.frame(knots)
  colnames(knots_df) = c("lat", "lon")
  
  knots_df = merge(knots_df, df, by = c("lat", "lon"), all.x = TRUE, all.y = FALSE)
  
  X = as.matrix(df[, c("X_1", "X_2","X_3")])
  Y_star = knots_df$Y
  X_star = as.matrix(knots_df[, c("X_1", "X_2","X_3")])
  
  m = nrow(knots)
  l = 0
  u = 4
  
  start_time <- Sys.time()
  output = GP_Gibbs_beta(Y = Y,
                         X = X,
                         s = coords,
                         knots = knots,
                         beta_knots_start = matrix(0, nrow = m, ncol = p),
                         w_knots_start = rep(0, m),
                         phi_beta_start = rep(l + (u - l) / 2, p),
                         phi_w_start = l + (u - l) / 2,
                         sigmasq_beta_start = rep(1, p),
                         sigmasq_w_start = 1,
                         tausq_start = 1,
                         phi_beta_proposal_sd = rep(0.1, p),
                         phi_w_proposal_sd = 0.1,
                         a_beta = rep(2, p),
                         b_beta = rep(1, p),
                         a_w = 2,
                         b_w = 1,
                         a_t = 2,
                         b_t = 1e-4,
                         lower_beta = rep(l, p),
                         upper_beta = rep(u, p),
                         lower_w = l,
                         upper_w = u,
                         mcmc = mcmc)
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  gibbs6 = rbind(gibbs6,c(mean(colMeans(output$beta_samples[, 1, 800:1000]) -3),mean(colMeans(output$beta_samples[, 2, 800:1000]) +5), time_taken ))
  start_time <- Sys.time()
  fit <- SVC_mle(Y, X, coords,as.matrix(w))
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  coef6 = rbind(coef6,c(fit$coefficients[1]-3,fit$coefficients[2]+5,time_taken))
}


gibbs7 = data.frame()
coef7 = data.frame()
### Scenario 4: (3 betas and "high" priors)
### 100 Sample size
for (i in 1:50) {
  
  lat = seq(0, 9, by = 1)
  lon = seq(0, 9, by = 1)
  
  coords = as.matrix(expand.grid(lat, lon))
  colnames(coords) = c("lat", "lon")
  
  n = nrow(coords)
  p = 2
  
  # generating w
  sigmasq_w = 4
  phi_w = 8
   
  C_w = calc_C_phi(coords, phi_w)
  w = MASS::mvrnorm(1, rep(0, n), sigmasq_w * C_w)
  
  # generating beta's
  sigmasq_1 = 4
  phi_1 = 8
   
  C_1 = calc_C_phi(coords, phi_1)
  beta_1_mean = 3
  beta_1 = MASS::mvrnorm(1, rep(beta_1_mean, n), sigmasq_1 * C_1)
  
  sigmasq_2 = 4
  phi_2 = 8
   
  C_2 = calc_C_phi(coords, phi_2)
  beta_2_mean = -5
  beta_2 = MASS::mvrnorm(1, rep(beta_2_mean, n), sigmasq_2 * C_2)
  
  sigmasq_3 = 4
  phi_3 = 8
   
  C_3 = calc_C_phi(coords, phi_3)
  beta_3_mean = -10
  beta_3 = MASS::mvrnorm(1, rep(beta_3_mean, n), sigmasq_3 * C_3)
  
  # generating X
  X_1 = rnorm(n, 0, 10)
  X_2 = rnorm(n, 0, 10)
  X_3 = rnorm(n, 0, 10)
  
  # generating epsilon
  tausq = 0.01
  epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))
  
  Y = X_1 * beta_1 + X_2 * beta_2 + w + epsilon
  
  k = 2
  
  lat_knots = unique(lat)
  lat_knots = lat_knots[seq(1, length(lat_knots), by = k)]
  
  lon_knots = unique(lon)
  lon_knots = lon_knots[seq(1, length(lon_knots), by = k)]
  
  knots = as.matrix(expand.grid(lat_knots, lon_knots))
  
  df = data.frame(coords, Y, X_1, X_2,X_3, w, beta_1, beta_2, epsilon)
  
  knots_df = data.frame(knots)
  colnames(knots_df) = c("lat", "lon")
  
  knots_df = merge(knots_df, df, by = c("lat", "lon"), all.x = TRUE, all.y = FALSE)
  
  X = as.matrix(df[, c("X_1", "X_2","X_3")])
  Y_star = knots_df$Y
  X_star = as.matrix(knots_df[, c("X_1", "X_2","X_3")])
  
  m = nrow(knots)
  l = 0
  u = 4
  
  start_time <- Sys.time()
  output = GP_Gibbs_beta(Y = Y,
                         X = X,
                         s = coords,
                         knots = knots,
                         beta_knots_start = matrix(0, nrow = m, ncol = p),
                         w_knots_start = rep(0, m),
                         phi_beta_start = rep(l + (u - l) / 2, p),
                         phi_w_start = l + (u - l) / 2,
                         sigmasq_beta_start = rep(1, p),
                         sigmasq_w_start = 1,
                         tausq_start = 1,
                         phi_beta_proposal_sd = rep(0.1, p),
                         phi_w_proposal_sd = 0.1,
                         a_beta = rep(2, p),
                         b_beta = rep(1, p),
                         a_w = 2,
                         b_w = 1,
                         a_t = 2,
                         b_t = 1e-4,
                         lower_beta = rep(l, p),
                         upper_beta = rep(u, p),
                         lower_w = l,
                         upper_w = u,
                         mcmc = mcmc)
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  gibbs7 = rbind(gibbs7,c(mean(colMeans(output$beta_samples[, 1, 800:1000]) -3),mean(colMeans(output$beta_samples[, 2, 800:1000]) +5), time_taken ))
  start_time <- Sys.time()
  fit <- SVC_mle(Y, X, coords,as.matrix(w))
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  coef7 = rbind(coef7,c(fit$coefficients[1]-3,fit$coefficients[2]+5,time_taken))
}

gibbs8 = data.frame()
coef8 = data.frame()
### 361 Sample size
for (i in 1:50) {
  lat = seq(0, 9, by = .4)
  lon = seq(0, 9, by = .4)
  
  coords = as.matrix(expand.grid(lat, lon))
  colnames(coords) = c("lat", "lon")
  
  n = nrow(coords)
  p = 2
  
  # generating w
  sigmasq_w = 4
  phi_w = 8
   
  C_w = calc_C_phi(coords, phi_w)
  w = MASS::mvrnorm(1, rep(0, n), sigmasq_w * C_w)
  
  # generating beta's
  sigmasq_1 = 4
  phi_1 = 8
   
  C_1 = calc_C_phi(coords, phi_1)
  beta_1_mean = 3
  beta_1 = MASS::mvrnorm(1, rep(beta_1_mean, n), sigmasq_1 * C_1)
  
  sigmasq_2 = 4
  phi_2 = 8
   
  C_2 = calc_C_phi(coords, phi_2)
  beta_2_mean = -5
  beta_2 = MASS::mvrnorm(1, rep(beta_2_mean, n), sigmasq_2 * C_2)
  
  sigmasq_3 = 4
  phi_3 = 8
   
  C_3 = calc_C_phi(coords, phi_3)
  beta_3_mean = -10
  beta_3 = MASS::mvrnorm(1, rep(beta_3_mean, n), sigmasq_3 * C_3)
  
  # generating X
  X_1 = rnorm(n, 0, 10)
  X_2 = rnorm(n, 0, 10)
  X_3 = rnorm(n, 0, 10)
  
  # generating epsilon
  tausq = 0.01
  epsilon = rnorm(n, mean = 0, sd = sqrt(tausq))
  
  Y = X_1 * beta_1 + X_2 * beta_2 + w + epsilon
  
  k = 2
  
  lat_knots = unique(lat)
  lat_knots = lat_knots[seq(1, length(lat_knots), by = k)]
  
  lon_knots = unique(lon)
  lon_knots = lon_knots[seq(1, length(lon_knots), by = k)]
  
  knots = as.matrix(expand.grid(lat_knots, lon_knots))
  
  df = data.frame(coords, Y, X_1, X_2,X_3, w, beta_1, beta_2, epsilon)
  
  knots_df = data.frame(knots)
  colnames(knots_df) = c("lat", "lon")
  
  knots_df = merge(knots_df, df, by = c("lat", "lon"), all.x = TRUE, all.y = FALSE)
  
  X = as.matrix(df[, c("X_1", "X_2","X_3")])
  Y_star = knots_df$Y
  X_star = as.matrix(knots_df[, c("X_1", "X_2","X_3")])
  
  m = nrow(knots)
  l = 0
  u = 4
  
  start_time <- Sys.time()
  output = GP_Gibbs_beta(Y = Y,
                         X = X,
                         s = coords,
                         knots = knots,
                         beta_knots_start = matrix(0, nrow = m, ncol = p),
                         w_knots_start = rep(0, m),
                         phi_beta_start = rep(l + (u - l) / 2, p),
                         phi_w_start = l + (u - l) / 2,
                         sigmasq_beta_start = rep(1, p),
                         sigmasq_w_start = 1,
                         tausq_start = 1,
                         phi_beta_proposal_sd = rep(0.1, p),
                         phi_w_proposal_sd = 0.1,
                         a_beta = rep(2, p),
                         b_beta = rep(1, p),
                         a_w = 2,
                         b_w = 1,
                         a_t = 2,
                         b_t = 1e-4,
                         lower_beta = rep(l, p),
                         upper_beta = rep(u, p),
                         lower_w = l,
                         upper_w = u,
                         mcmc = mcmc)
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  gibbs8 = rbind(gibbs8,c(mean(colMeans(output$beta_samples[, 1, 800:1000]) -3),mean(colMeans(output$beta_samples[, 2, 800:1000]) +5), time_taken ))
  start_time <- Sys.time()
  fit <- SVC_mle(Y, X, coords,as.matrix(w))
  end_time <- Sys.time()
  time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  coef8 = rbind(coef8,c(fit$coefficients[1]-3,fit$coefficients[2]+5,time_taken))
}








