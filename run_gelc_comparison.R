library(nimble)
library(nimbleNoBounds)
library(coda)
library(tidyverse)
library(ICenCov)

dir <- '/work/users/a/g/agfoes/P1'
run_name <- "small_gelc_comparison_P1_gamma"

source(file.path(dir, 'R', 'GELc', 'config_gelc_comparison.R'))
source(file.path(dir, 'R', 'helpers', 'helpers_data.R'))
source(file.path(dir, 'R', 'helpers', 'helpers_nimble.R'))
source(file.path(dir, 'R', 'helpers', 'helpers_sampler.R'))
source(file.path(dir, "R", "model_code", "gamma_all_cens_int_model_code.R"))
source(file.path(dir, "R", "custom_block_sampler.R"))


result_dir <- file.path(
  dir,
  'results',
  run_name
)

summary_dir <- file.path(
  result_dir,
  'replicate_results'
)

data_dir <- file.path(
  result_dir,
  'data'
)

dir.create(
  summary_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  data_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

rep <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

if(is.na(rep)) {
  rep <- 1L
}

result_file <- file.path(
  summary_dir,
  sprintf('results_rep_%04d.csv', rep)
)

if(file.exists(result_file)) {
  message('Replication ', rep, ' already completed.')
  quit(save = 'no')
}

## ------------------------------------------------------------
## Generate complete data once - different censoring levels to be applied on top of this full dataset
## ------------------------------------------------------------

set.seed(
  09102026 +
    rep
)

data_full <- datagen_gelc_gamma(
  n = 500,
  mu = 0,
  gamma = 0.02,
  phi = 0.02
)

saveRDS(
  data_full,
  file.path(
    data_dir,
    sprintf('complete_data_rep_%04d.rds', rep)
  )
)

rep_results <- list()

row_id <- 0

## ------------------------------------------------------------
## Generate censoring for each mu
## ------------------------------------------------------------
mu_values <- c(3, 6, 9, 12)
n_values <- c(100, 300, 500)
for(mu in mu_values) {
  
  ## separate deterministic censoring seed
  set.seed(
    09102026 +
      1000000 +
      rep * 100 +
      mu
  )
  
  censoring <- gelc_censoring(
    n = 500,
    mu = mu,
    Z = data_full$Z
  )
  
  data_mu <- data_full %>%
    mutate(CL = censoring$CL,
           CR = censoring$CR)
  
  ## ----------------------------------------------------------
  ## Nested n = 100, 300, 500
  ## ----------------------------------------------------------
  
  for(n in n_values) {
    
    message(
      'rep = ', rep,
      ', n = ', n,
      ', mu = ', mu
    )
    
    data <- data_mu[seq_len(n), , drop = FALSE]
    
    ## ========================================================
    ## oracle
    ## ========================================================
    
    oracle_time <- system.time({
      fit_oracle <- glm(Y ~ -1 + X1 + Z, data = data, family = Gamma(link = "log"))
      })[['elapsed']]
    
    oracle_coef <- summary(fit_oracle)$coefficients
    
    row_id <- row_id + 1
    
    rep_results[[row_id]] <- data.frame(
      rep = rep,
      n = n,
      mu = mu,
      method = 'oracle',
      parameter = 'gamma',
      truth = 0.02,
      estimate = oracle_coef['Z', 'Estimate'],
      se = oracle_coef['Z', 'Std. Error'],
      lower = oracle_coef['Z', 'Estimate'] - 1.96 * oracle_coef['Z', 'Std. Error'],
      upper = oracle_coef['Z', 'Estimate'] + 1.96 * oracle_coef['Z', 'Std. Error'],
      runtime = as.numeric(oracle_time)
    )
    
    ## ========================================================
    ## GELc
    ## ========================================================
    
    gelc_time <- system.time({
      fit_gelc <- icglm(Y ~ -1 + X1 + ic(CL, CR, 'Z'), family = Gamma(link = "log"), data = data)
    })[['elapsed']]
      
    gelc_summary <- summary(fit_gelc)
    
    gelc_coef <- gelc_summary$coefficients
    
    row_id <- row_id + 1
    
    rep_results[[row_id]] <- data.frame(
      rep = rep,
      n = n,
      mu = mu,
      method = 'GELc',
      parameter = 'gamma',
      truth = 0.02,
      estimate = gelc_coef['Z', 'Estimate'],
      se = gelc_coef['Z', 'Std. Error'],
      lower = gelc_coef['Z', 'Estimate'] - 1.96 * gelc_coef['Z', 'Std. Error'],
      upper = gelc_coef['Z', 'Estimate'] + 1.96 * gelc_coef['Z', 'Std. Error'],
      runtime = as.numeric(gelc_time)
    )
    
    ## ========================================================
    ## P1
    ## ========================================================
    
    pz <- 1
    L <- 50
    
    Z <- rep(1, n)
    
    idx_obs <- which(data$Dobs == 1)
    idx_cen <- which(data$Dobs == 0)
    
    nobs <- length(idx_obs)
    ncen <- length(idx_cen)
    
    if(length(idx_obs) < 2) {
      idx_obs <- c(
        idx_obs,
        rep(1L, 2 - length(idx_obs))
      )
    }
    
    if(length(idx_cen) < 2) {
      idx_cen <- c(
        idx_cen,
        rep(1L, 2 - length(idx_cen))
      )
    }
    
    Nconstants <- list(
      L = L,
      p = pz,
      mu_gamma = 0,
      nobs = nobs,
      ncen = ncen,
      idx_obs = idx_obs,
      idx_cen = idx_cen,
      beta_mean = rep(0, pz + 1),
      beta_cov = diag(pz + 1)
    )
    
    Ndata <- list(
      y_obs = data$Y[data$Dobs == 1],
      y_cen = data$Y[data$Dobs == 0],
      z_obs = as.vector(data$X1[data$Dobs == 1]),
      z_cen = as.vector(data$X1[data$Dobs == 0]),
      x_obs = data$Z[data$Dobs == 1],
      CL = data$CL[data$Dobs == 0],
      CR = data$CR[data$Dobs == 0],
      constraint_data = rep(1, ncen)
    )

    Ninits <- nimble_inits(
      data = data,
      pz = pz,
      constants = Nconstants,
      Ndata = Ndata
    )
    Ninits[["gammaTilde"]] <- as.vector(Ninits[["gammaTilde"]])
    
    p1_time <- system.time({
      model <- nimbleModel(
        code = model_code,
        constants = Nconstants,
        data = Ndata,
        inits = Ninits
      )
      
      set.seed(
        config[['sim_seed']] +
          2000000 +
          rep * 10000 +
          n * 10 +
          mu
      )
      
      p1_result <- run_sampler(
        model = model,
        sampler_type = "AF_slice",
        pz = pz,
        L = L,
        niter = 100000,
        nburnin = 40000,
        thin = 1,
        vars_to_monitor = c("beta", "tau"),
        ncen = ncen,
        nobs = nobs
      )
    })[['elapsed']]
      
      beta_x_samples <-
        p1_result$samples[ , 'beta[1]']
      
      beta_x_hpd <- coda::HPDinterval(
        as.mcmc(beta_x_samples),
        prob = 0.95
      )
      
      row_id <- row_id + 1
      
      rep_results[[row_id]] <- data.frame(
        rep = rep,
        n = n,
        mu = mu,
        method = 'P1',
        parameter = 'gamma',
        truth = 0.02,
        estimate = mean(beta_x_samples),
        se = sd(beta_x_samples),
        lower = beta_x_hpd[1, 'lower'],
        upper = beta_x_hpd[1, 'upper'],
        runtime = as.numeric(p1_time)
      )
    
    rm(
      fit_oracle,
      fit_gelc,
      p1_result,
      model
    )
    
    gc()
  }
}

## ------------------------------------------------------------
## Save complete replication
## ------------------------------------------------------------

rep_results <- bind_rows(
  rep_results
) %>%
  mutate(
    bias = estimate - truth,
    squared_error =
      (estimate - truth)^2,
    covered =
      truth >= lower &
      truth <= upper,
    interval_width =
      upper - lower
  )

write.csv(
  rep_results,
  result_file,
  row.names = FALSE
)
