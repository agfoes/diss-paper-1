library(nimble)
library(nimbleNoBounds)
library(coda)
library(tidyverse)

dir <- '/work/users/a/g/agfoes/P1'

source(file.path(dir, 'R', 'GELc', 'config_gelc_comparison.R'))
source(file.path(dir, 'R', 'helpers', 'helpers_data.R'))
source(file.path(dir, 'R', 'helpers', 'helpers_nimble.R'))
source(file.path(dir, 'R', 'helpers', 'helpers_sampler.R'))
source(file.path(dir, 'R', 'helpers', 'helpers_summary.R'))
source(file.path(dir, 'R', 'custom_block_sampler.R'))

includeObs <- TRUE
includeCen <- TRUE
source(file.path(dir, 'R', 'model_code', 'model_normal_joint_betas.R'))

result_dir <- file.path(config[['project_dir']], 'results', config[['run_name']])
summary_dir <- file.path(result_dir, 'replicate_results')
data_dir <- file.path(result_dir, 'data')

dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

rep <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

if(is.na(rep)) {
  rep <- 1L
}

result_file <- file.path(summary_dir, sprintf('P1_results_rep_%04d.csv', rep))

if(file.exists(result_file)) {
  message('Replication ', rep, ' already completed.')
  quit(save = 'no')
}

## ------------------------------------------------------------
## Generate complete data ONCE
## ------------------------------------------------------------

set.seed(config[['sim_seed']] + rep)

data_full <- datagen_gelc_normal(
  n = max(config[['n_values']]),
  beta_x = config[['beta_x']],
  beta_z = config[['beta_z']],
  beta_0 = config[['beta_0']],
  tau = config[['tau']]
)

rep_results <- list()

row_id <- 0

## ------------------------------------------------------------
## Generate censoring for each mu
## ------------------------------------------------------------

for(mu in config[['mu_values']]) {
  
  ## separate deterministic censoring seed
  set.seed(
    config[['sim_seed']] + 1000000 + rep * 100 + mu
  )
  
  censoring <- gelc_censoring(
    x = data_full$X,
    mu = mu
  )
  
  data_mu <- cbind(
    data_full,
    censoring
  )
  
  ## temporary P1 workaround
  data_mu$CL[1] <- data_mu$X[1]
  data_mu$CR[1] <- data_mu$X[1]
  data_mu$Dobs[1] <- 1L
  
  ## ----------------------------------------------------------
  ## Nested n = 100, 300, 500
  ## ----------------------------------------------------------
  
  for(n in config[['n_values']]) {
    
    message(
      'rep = ', rep,
      ', n = ', n,
      ', mu = ', mu
    )
    
    data <- data_mu[seq_len(n), , drop = FALSE]
    
    ## ========================================================
    ## P1
    ## ========================================================
    
    Z <- as.matrix(data[, c('Z1', 'Z2')])
    
    idx_obs <- which(data$Dobs == 1)
    idx_cen <- which(data$Dobs == 0)
    
    nobs <- length(idx_obs)
    ncen <- length(idx_cen)
    
    idx_obs <- c(idx_obs, rep(1L, max(0, 2 - length(idx_obs))))
    idx_cen <- c(idx_cen, rep(1L, max(0, 2 - length(idx_cen))))
    
    constants <- list(
      L = config[['L']],
      p = config[['pz']],
      mu_gamma = config[['mu_gamma']],
      nobs = nobs,
      ncen = ncen,
      idx_obs = idx_obs,
      idx_cen = idx_cen,
      beta_mean = rep(0, config[['pz']] + 1),
      beta_cov = diag(config[['pz']] + 1)
    )
    
    Ndata <- list(
      y_obs = data$Y[data$Dobs == 1],
      y_cen = data$Y[data$Dobs == 0],
      z_obs = Z[data$Dobs == 1, , drop = FALSE],
      z_cen = Z[data$Dobs == 0, , drop = FALSE],
      x_obs = data$X[data$Dobs == 1],
      CL = data$CL[data$Dobs == 0],
      CR = data$CR[data$Dobs == 0],
      constraint_data = rep(1, ncen)
    )
    
    Ninits <- nimble_inits(
      data = data,
      pz = config[['pz']],
      constants = constants,
      Ndata = Ndata
    )
    
    model <- nimbleModel(
      code = model_code,
      constants = constants,
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
    
    p1_time <- system.time({
      
      p1_result <- run_sampler(
        model = model,
        sampler_type = config[['sampler']],
        pz = config[['pz']],
        L = config[['L']],
        niter = config[['niter']],
        nburnin = config[['nburnin']],
        thin = config[['thin']],
        vars_to_monitor = config[['key_vars']],
        ncen = ncen,
        nobs = nobs
      )
      
    })[['elapsed']]
    
    p1_summary <- posterior_summary(
      samples = p1_result$samples,
      sampler_name = config[['sampler']],
      truth = config[['truth']],
      key_vars = config[['key_vars']],
      runtime = p1_time
    )
    
    beta_x_summary <- p1_summary %>%
      filter(parameter == 'beta[1]')
    
    row_id <- row_id + 1
    
    rep_results[[row_id]] <- data.frame(
      rep = rep,
      n = n,
      mu = mu,
      method = 'P1',
      parameter = 'beta_x',
      truth = beta_x_summary$truth,
      estimate = beta_x_summary$mean,
      se = beta_x_summary$sd,
      lower = beta_x_summary$hpd_lower,
      upper = beta_x_summary$hpd_upper,
      runtime = beta_x_summary$runtime
    )
    
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

message('Completed replication ', rep)