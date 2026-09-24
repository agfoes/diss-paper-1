# purpose: evaluate P1 and GELc performance with true multimodal X|Z in DGP

# libraries ----
library(nimble)
library(nimbleNoBounds)
library(coda)
library(tidyverse)
library(ICenCov)


datagen_gamma <- function(n = 500,
                          mu = 3,
                          gamma = 0.02,
                          alpha = 10,
                          phi = 0.02,
                          mode = 1,
                          delta = 0.5) {
  
  ## observed covariates for both the outcome and marginal regression models
  X1 <- rep(1, n)
  X2 <- rbinom(n = n, size = 1, prob = 0.8)
  X3 <- rnorm(n = n, mean = 40, sd = 1)
  
  ## outcome mean and variance
  if (mode == 1) {
    component_mean = 16
  } else if (mode == 2) {
    component_mean = c(10, 22)
  } else if (mode == 3) {
    component_mean = c(8, 16, 24)
  }
  
  comp <- sample(x = seq_len(mode), size = n,
                 replace = TRUE, prob = rep(1/mode, mode))
  
  mean_Z <- component_mean[comp] + delta * (X2 - 0.8)
  
  Z <- rnorm(n = n, mean = mean_Z, sd = 2)
  mean_Y <- exp(alpha + gamma * Z)
  var_Y <- phi*mean_Y^2
  
  Y <- rgamma(n, shape = (mean_Y^2 / var_Y), rate = (mean_Y / var_Y))
  
  ## GELc censoring mechanism
  CL <- rep(0, n)
  CR <- rep(0, n)
  
  if (mu == 0) {
    CL <- Z
    CR <- Z
  } else {
    censoring <- gelc_censoring(mu, n, Z)
    CL <- censoring$CL
    CR <- censoring$CR
  }
  
  return(data.frame(
    Y = Y,
    Z = Z,
    X1 = X1,
    X2 = X2,
    X3 = X3 - 40,
    CL = CL,
    CR = CR
  ))
}

gelc_censoring <- function(mu = 0,
                           n = 500,
                           Z = rep(0, 500)) {
  
  CL <- numeric(n)
  CR <- numeric(n)
  
  if (mu == 0) {
    CL <- Z
    CR <- Z
  } else {
    for (i in 1:n) {
      eps0 <- runif(1, min = 0, max = mu)
      visits <- c(0, eps0)
      
      # init first visit value for loop
      visit <- visits[2]
      
      while(max(visits) <= Z[i]) {
        gap <- abs(rnorm(1, mean = mu, sd = sqrt(0.75*mu)))
        visit <- visit + gap
        visits <- c(visits, visit)
      }
      
      CL[i] <- max(visits[visits <= Z[i]])
      CR[i] <- min(visits[visits > Z[i]])
    }
  }
  
  
  
  return(data.frame(CL = CL, CR = CR))
}




main_folder <- "/work/users/a/g/agfoes/P1"
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if (task_id <= 100) {
  rep <- task_id
  mode = 1
  scenario <- "unimodal"
} else if (task_id <= 200) {
  rep <- task_id - 100
  mode = 2
  scenario <- "bimodal"
} else if (task_id <= 300) {
  rep <- task_id - 200
  mode = 3
  scenario <- "trimodal"
}

set.seed(09142026 + 1000 * rep + as.integer(mode))
n <- 500
delta <- 0.5
dat <- datagen_gamma(n = n,
                     mu = 9,
                     gamma = 0.02,
                     alpha = 10,
                     phi = 0.02,
                     mode = mode,
                     delta = delta)


gelc_time <- system.time({
  
  fit_gelc <- icglm(
    Y ~ X2 + X3 + ic(CL, CR, 'Z'),
    family = Gamma(link = "log"),
    data = dat
  )
  
})['elapsed']

gelc_summary <- summary(fit_gelc)
gelc_coef <- gelc_summary$coefficients
res <- as.data.frame(gelc_coef) %>%
  mutate(rep = rep,
         n = n,
         mu = mu,
         scenario = scenario,
         mode = mode,
         delta = delta,
         corr_ZX3 = cor(dat$Z, dat$X3),
         corr_ZX2 = cor(dat$Z, dat$X2),
         mean_width = mean(dat$CR - dat$CL),
         median_width = median(dat$CR - dat$CL),
         sd_Z = sd(dat$Z))%>%
  rownames_to_column("parameter") %>%
  mutate(runtime = as.numeric(gelc_time))


write.csv(res, 
          file = file.path(main_folder, "results/multimodal_tests/mu_9_three_modes/gelc", paste0("results_", task_id, ".csv")),
          row.names = FALSE)
