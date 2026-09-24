# purpose of this code is to run simulations testing GELc performance with Z independent vs dependent on X with a gamma outcome - following other data generation specifications of GELc. 

# libraries ----
library(nimble)
library(nimbleNoBounds)
library(coda)
library(tidyverse)
library(ICenCov)


datagen_gelc_gamma <- function(n = 500,
                               mu = 3,
                               gamma = 0.02,
                               alpha = 10,
                               phi = 0.02,
                               dep = 0,
                               delta = 0.5) {
  
  ## observed covariates for both the outcome and marginal regression models
  X1 <- rep(1, n)
  X2 <- rbinom(n = n, size = 1, prob = 0.8)
  X3 <- rnorm(n = n, mean = 40, sd = 1)
  
  ## outcome mean and variance
  if (dep == 0) {
    rateZ <- 1/12
  } else if (dep == 1) {
    rateZ <- exp(log(1/12) + delta*(X2 - 0.8))
  } else if (dep == 2) {
    rateZ <- exp(log(1/12) + delta*(X3 - 40))
  } else if (dep == 3) {
    rateZ <- exp(log(1/12) + delta*(X2 - 0.8) + delta*(X3 - 40))
  }
  
  Z <- rexp(n, rate = rateZ)
  mean <- exp(alpha + gamma * Z)
  var <- phi*mean^2
  
  Y <- rgamma(n, shape = (mean^2 / var), rate = (mean / var))
  
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

if (task_id <= 500) {
  rep <- task_id
  dep <- 0
  scenario <- "independent"
} else if (task_id <= 1000 & task_id > 500) {
  rep <- task_id - 500
  dep <- 1
  scenario <- "binary dependent Z"
} else if (task_id > 1000 & task_id <= 1500) {
  rep <- task_id - 1000
  dep <- 2
  scenario <- "continuous dependent Z"
} else if (task_id > 1500 & task_id <= 2000) {
  rep <- task_id - 1500
  dep <- 3
  scenario <- "binary + continuous dependent Z"
}

set.seed(09142026 + 1000 * rep + as.integer(dep))
n <- 500
delta <- 1.5
dat <- datagen_gelc_gamma(n = n,
                          mu = 6,
                          gamma = 0.02,
                          alpha = 10,
                          phi = 0.02,
                          dep = dep,
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
         n = 500,
         mu = 6,
         scenario = scenario,
         dep = dep,
         delta = delta,
         corr_ZX3 = cor(dat$Z, dat$X3),
         corr_ZX2 = cor(dat$Z, dat$X2)) %>%
  rownames_to_column("parameter") %>%
  mutate(runtime = as.numeric(gelc_time))


write.csv(res, file = file.path(main_folder, "results/dependence_gelc/delta_1-5_GELc", paste0("results_", task_id, ".csv")))
