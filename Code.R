###################################
## GPD STRATIFICATO PER STAGIONE ##
###################################
library(extRemes)
library(GSODR)
library(tidyverse)
library(urca)
library(lubridate)
library(evd)

set.seed(123)

######### CARICAMENTO E PREPARAZIONE #########
station <- nearest_stations(LAT = 45.46, LON = 9.19, distance = 50)[2, 1]
milano <- get_GSOD(years = 1976:2018, station = station)
str(milano)
nrow(milano)
ncol(milano)
head(milano$YEARMODA)
tail(milano$YEARMODA)
n_years <- round(nrow(milano)/365)
n_years == 2018 - 1976 + 1

milano2 <- milano %>%
  filter(!is.na(MIN)) %>%
  mutate(
    NEG_MIN = -MIN,
    SEASON = factor(case_when(
      MONTH %in% c(12, 1, 2)  ~ "DJF",
      MONTH %in% c(3, 4, 5)   ~ "MAM",
      MONTH %in% c(6, 7, 8)   ~ "JJA",
      MONTH %in% c(9, 10, 11) ~ "SON"
    ), levels = c("DJF", "MAM", "JJA", "SON"))
  )

data_list <- split(milano2$NEG_MIN, milano2$SEASON)
seasons <- names(data_list)

######### ANALISI ESPLORATIVA #########

# Verifica stazionarietà
plot(milano$MIN, type = "l", col = "blue",
     xlab = "Giorni", ylab = "Tmin (°C)", main = "Dati originali - Milano")

par(mfrow = c(2, 2))
for (s in seasons) plot(data_list[[s]], type = "l", col = "blue",
                        xlab = "Giorni", main = s, ylab = "-Tmin (°C)")

adf_results <- lapply(data_list, function(x) ur.df(x, type = "trend", lags = 3))
for (s in names(adf_results)) {
  print(summary(adf_results[[s]]))
}

# Verifica della dipendenza
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
for (s in seasons) {
  acf(data_list[[s]], lag.max = 30)
  pacf(data_list[[s]], lag.max = 30)
}

######### STIMA DI GPD #########

# Scelta delle soglie
par(mfrow = c(2, 2))
for (s in seasons) mrlplot(data_list[[s]], main = s)

thresholds <- c(DJF = 5, MAM = 0.5, JJA = -12.5, SON = 0)

# stabilità al variare di r
r_values <- c(1, 2, 3, 4)
results <- expand.grid(r = r_values, Season = seasons, stringsAsFactors = FALSE)
results$u <- results$n_exc <- results$n_clust <- results$theta <- NA
results$sigma <- results$se_sigma <- results$xi <- results$se_xi <- NA

for (i in seq_len(nrow(results))) {
  r <- results$r[i]
  s <- results$Season[i]
  x <- data_list[[s]]
  u <- thresholds[s]
  
  # Declustering
  cluster_data <- evd::clusters(x, u = u, r = r, cmax = TRUE)
  
  n_exc <- sum(x > u)
  n_clust <- length(cluster_data)
  
  fit <- fevd(cluster_data, threshold = u, type = "GP")
  
  pars <- fit$results$par
  cov_mat <- solve(fit$results$hessian)
  ses <- sqrt(diag(cov_mat))
  
  results[i, c("u", "n_exc", "n_clust")] <- c(u, n_exc, n_clust)
  results$theta[i] <- round(n_clust / n_exc, 3)
  results$sigma[i] <- round(pars[1], 3)
  results$se_sigma[i] <- round(ses[1], 3)
  results$xi[i] <- round(pars[2], 3)
  results$se_xi[i] <- round(ses[2], 3)
}

print(results) 

# Fit GPD con dati declusterizzati son r selezionato
r_opt <- 2

gpd_fits <- list()
for (s in seasons) {
  x <- data_list[[s]]
  u <- thresholds[s]
  clusters <- evd::clusters(x, u = u, r = r_opt, cmax = TRUE)
  n_exc <- sum(x > u)
  n_clust <- length(clusters)
  fit <- fevd(clusters, threshold = u, type = "GP")
  
  gpd_fits[[s]] <- list(fit = fit, theta = n_clust / n_exc, 
                        n_exc = n_exc, n_clust = n_clust, u = u)
}

print(gpd_fits)

######### LRT: GPD vs ESPONENZIALE #########
lrt_results <- data.frame()

for (s in names(gpd_fits)) {
  fit_gpd <- gpd_fits[[s]]$fit
  fit_exp <- fevd(fit_gpd$x, threshold = fit_gpd$threshold, type = "Exponential")
  test <- lr.test(fit_exp, fit_gpd)
  
  lrt_results <- rbind(lrt_results, data.frame(
    Season = s,
    LR_Stat = round(test$statistic, 3),
    p_val_std = round(test$p.value, 4)
  ))
}

print(lrt_results)

######### MODELLI FINALI #########
# True a chi rifiuta H0 del LRT
use_gpd <- c(DJF = FALSE, MAM = FALSE, JJA = FALSE, SON = TRUE)

final_models <- list()
for (s in names(gpd_fits)) {
  gf <- gpd_fits[[s]]
  clusters <- gf$fit$x
  u <- gf$u
  
  if (use_gpd[s]) { 
    fit <- gf$fit
    sigma <- fit$results$par[["scale"]]
    xi <- fit$results$par[["shape"]]
    V <- solve(fit$results$hessian)
    type <- "GPD"
    
  } else {
    fit <- fevd(clusters, threshold = u, type = "Exponential")
    sigma <- fit$results$par[["scale"]]
    xi <- 0
    n <- length(clusters)
    v_scale <- sigma^2 / n
    V <- matrix(c(v_scale, 0, 0, 0), nrow = 2, ncol = 2)
    type <- "Exp"
  }
  
  final_models[[s]] <- list(
    type = type,
    sigma = sigma,
    xi = xi,
    V = V,
    u = u,
    theta = gf$theta,
    n_exc = gf$n_exc,
    n_clust = gf$n_clust,
    data = clusters
  )
}

print(final_models)

######### DIAGNOSTICA #########
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))
for (s in names(final_models)) {
  model <- final_models[[s]]
  y <- model$data - model$u
  sigma <- model$sigma
  xi <- model$xi
  
  if (abs(xi) > 1e-6) {
    z <- (1 / xi) * log(1 + xi * y / sigma)
  } else {
    z <- y / sigma
  }
  
  z <- sort(z)
  n <- length(z)
  p_emp <- (1:n) / (n + 1)
  q_theo <- -log(1 - p_emp)
  p_model <- 1 - exp(-z)
  
  plot(p_emp, p_model, pch = 19, cex = 0.8,
       main = paste("PP -", s), xlab = "Empirica", ylab = "Modello")
  abline(0, 1, col = "red", lwd = 2)
  
  plot(q_theo, z, pch = 19, cex = 0.8,
       main = paste("QQ -", s), xlab = "Exp(1)", ylab = "Standardizzati")
  abline(0, 1, col = "red", lwd = 2)
}
par(mfrow = c(1, 1))

######### RETURN LEVELS CON IC #########
return_level_ci <- function(model, T_years, n_per_year, zeta_u, alpha = 0.05) {
  sigma <- model$sigma
  xi <- model$xi
  u <- model$u
  V <- model$V
  
  m <- T_years * n_per_year * zeta_u * model$theta
  
  if (abs(xi) > 1e-6) {
    z_T <- u + (sigma / xi) * (m^xi - 1)
    
    # Calcolo delle derivate parziali per IC
    dz_dsigma <- (m^xi - 1) / xi
    dz_dxi <- (sigma / xi) * (m^xi * log(m) - (m^xi - 1) / xi)
  } else {
    z_T <- u + sigma * log(m)
    
    # Calcolo delle derivate parziali per IC
    dz_dsigma <- log(m)
    dz_dxi <- 0
  }
  
  # Metodo delta per IC
  grad <- c(dz_dsigma, dz_dxi)
  se_z <- sqrt(t(grad) %*% V %*% grad)
  z_crit <- qnorm(1 - alpha/2)
  
  data.frame(z_T = z_T, se = as.numeric(se_z),
             ci_lower = as.numeric(z_T - z_crit * se_z),
             ci_upper = as.numeric(z_T + z_crit * se_z))
}

# Selezione degli anni
T_ret <- c(10, 20, 50, 100)
rl_results <- data.frame()

for (s in names(final_models)) {
  mod <- final_models[[s]]
  n_obs <- length(data_list[[s]])
  n_per_year <- n_obs / n_years
  zeta_u <- mod$n_exc / n_obs
  
  for (T in T_ret) {
    rl <- return_level_ci(mod, T, n_per_year, zeta_u)
    rl_results <- rbind(rl_results, data.frame(
      Season = s, T_years = T,
      Tmin = round(-rl$z_T, 2),
      CI_lower = round(-rl$ci_upper, 2),
      CI_upper = round(-rl$ci_lower, 2)
    ))
  }
}

for(s in unique(rl_results$Season)) {
  print(rl_results %>% filter(Season == s))
}




