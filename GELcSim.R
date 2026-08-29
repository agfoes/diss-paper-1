## ========================================================================
## GELc-style simulation comparison: local pilot
## ========================================================================

library(nimble)
library(nimbleNoBounds)
library(coda)
library(tidyverse)
library(ICenCov)

## project setup ------------------------------------------------------------------------

project_dir <- getwd()

source(file.path(project_dir, 'R', 'helpers', 'helpers_data.R'))
source(file.path(project_dir, 'R', 'helpers', 'helpers_nimble.R'))
source(file.path(project_dir, 'R', 'helpers', 'helpers_sampler.R'))

joint_model_code <- source(
  file.path(
    project_dir,
    'R',
    'model_code',
    'model_normal_joint_betas.R'
  )
)$value

## simulation settings ------------------------------------------------------------------------

sim_seed <- 20260828

nrep <- 2

n_values <- c(100, 300, 500)
mu_values <- c(3, 6, 9)

L <- 25

niter <- 100000
nburnin <- 40000
thin <- 1

beta_x_true <- 0.1
beta_0_true <- 2
beta_z_true <- 0
tau_true <- 1

## result storage ------------------------------------------------------------------------

result_dir <- file.path(
  project_dir,
  'results',
  'gelc_comparison_local'
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

result_file <- file.path(
  result_dir,
  'simulation_results.csv'
)

##  simulation loop ------------------------------------------------------------------------

sim_id <- 0

for(rep in seq_len(nrep)) {
  
  ## ---------------------------------------------------------------
  ## Generate one underlying n = 500 dataset per censoring setting.
  ##
  ## The GELc simulation generates n = 500 and then takes the first
  ## 100, 300, or 500 observations for the corresponding n scenario.
  ## ---------------------------------------------------------------
  
  for(mu in mu_values) {
    
    set.seed(
      sim_seed +
        rep * 1000 +
        mu * 10
    )
    
    data_full <- datagen_gelc_normal(
      n = 500,
      mu = mu,
      beta_x = beta_x_true,
      beta_z = beta_z_true,
      beta_0 = beta_0_true,
      tau = tau_true,
      force_observed = TRUE
    )
    
    for(n in n_values) {
      
      sim_id <- sim_id + 1
      
      message(
        'rep = ', rep,
        ', n = ', n,
        ', mu = ', mu,
        ', sim_id = ', sim_id
      )
      
      data <- data_full[seq_len(n), ]
      
      ## ============================================================
      ## Oracle model
      ## ============================================================
      
      oracle_time <- system.time({
        
        fit_oracle <- lm(
          Y ~ Z2 + X,
          data = data
        )
        
      })['elapsed']
      
      oracle_coef <- summary(fit_oracle)$coefficients
      
      oracle_row <- data.frame(
        sim_id = sim_id,
        rep = rep,
        n = n,
        mu = mu,
        method = 'oracle',
        parameter = 'beta_x',
        truth = beta_x_true,
        estimate = oracle_coef['X', 'Estimate'],
        se = oracle_coef['X', 'Std. Error'],
        lower = oracle_coef['X', 'Estimate'] -
          1.96 * oracle_coef['X', 'Std. Error'],
        upper = oracle_coef['X', 'Estimate'] +
          1.96 * oracle_coef['X', 'Std. Error'],
        runtime = as.numeric(oracle_time)
      )
      
      ## ============================================================
      ## GELc
      ## ============================================================
      
      gelc_time <- system.time({
        
        fit_gelc <- icglm(
          Y ~ Z2 + ic(CL, CR, 'X'),
          family = gaussian,
          data = data
        )
        
      })['elapsed']
      
      gelc_summary <- summary(fit_gelc)
      
      gelc_coef <- gelc_summary$coefficients
      
      gelc_row <- data.frame(
        sim_id = sim_id,
        rep = rep,
        n = n,
        mu = mu,
        method = 'GELc',
        parameter = 'beta_x',
        truth = beta_x_true,
        estimate = gelc_coef['X', 'Estimate'],
        se = gelc_coef['X', 'Std. Error'],
        lower = gelc_coef['X', 'Estimate'] -
          1.96 * gelc_coef['X', 'Std. Error'],
        upper = gelc_coef['X', 'Estimate'] +
          1.96 * gelc_coef['X', 'Std. Error'],
        runtime = as.numeric(gelc_time)
      )
      
      ## ============================================================
      ## P1
      ## ============================================================
      
      pz <- 2
      
      Z <- as.matrix(
        data[, c('Z1', 'Z2')]
      )
      
      idx_obs <- which(data$Dobs == 1)
      idx_cen <- which(data$Dobs == 0)
      
      nobs <- length(idx_obs)
      ncen <- length(idx_cen)
      
      ## NIMBLE singleton-index workaround
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
      
      constants <- list(
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
        
        z_obs = Z[
          data$Dobs == 1,
          ,
          drop = FALSE
        ],
        
        z_cen = Z[
          data$Dobs == 0,
          ,
          drop = FALSE
        ],
        
        x_obs = data$X[data$Dobs == 1],
        
        CL = data$CL[data$Dobs == 0],
        CR = data$CR[data$Dobs == 0],
        
        constraint_data = rep(1, ncen)
      )
      
      includeObs <- nobs > 0
      includeCen <- ncen > 0
      
      Ninits <- nimble_inits(
        data = data,
        pz = pz,
        constants = constants,
        Ndata = Ndata
      )
      
      model <- nimbleModel(
        code = joint_model_code,
        constants = constants,
        data = Ndata,
        inits = Ninits
      )
      
      if(!is.finite(model$calculate())) {
        stop(
          'Non-finite initial model calculation for sim_id = ',
          sim_id
        )
      }
      
      set.seed(
        sim_seed +
          1000000 +
          sim_id
      )
      
      fit_p1 <- run_sampler(
        model = model,
        sampler_type = 'custom_joint',
        pz = pz,
        L = L,
        niter = niter,
        nburnin = nburnin,
        thin = thin,
        vars_to_monitor = c(
          'beta',
          'tau'
        ),
        ncen = ncen,
        nobs = nobs
      )
      
      beta_x_samples <- fit_p1$samples[, 'beta[1]']
      
      beta_x_hpd <- HPDinterval(
        as.mcmc(beta_x_samples),
        prob = 0.95
      )
      
      p1_row <- data.frame(
        sim_id = sim_id,
        rep = rep,
        n = n,
        mu = mu,
        method = 'P1',
        parameter = 'beta_x',
        truth = beta_x_true,
        estimate = mean(beta_x_samples),
        se = sd(beta_x_samples),
        lower = beta_x_hpd[1, 'lower'],
        upper = beta_x_hpd[1, 'upper'],
        runtime = as.numeric(fit_p1$runtime)
      )
      
      ## ============================================================
      ## combine and save immediately
      ## ============================================================
      
      results <- bind_rows(
        oracle_row,
        gelc_row,
        p1_row
      ) %>%
        mutate(
          bias = estimate - truth,
          squared_error = (estimate - truth)^2,
          covered = truth >= lower &
            truth <= upper,
          interval_width = upper - lower
        )
      
      write.table(
        results,
        file = result_file,
        sep = ',',
        row.names = FALSE,
        col.names = !file.exists(result_file),
        append = file.exists(result_file)
      )
      
      rm(
        fit_oracle,
        fit_gelc,
        fit_p1,
        model,
        Ninits,
        Ndata
      )
      
      gc()
    }
  }
}

results <- read.csv(
  'results/gelc_comparison_local/simulation_results.csv'
)

dim(results)
head(results)
table(results$n, results$mu, results$method)

results[
  order(results$rep, results$mu, results$n, results$method),
  c(
    'rep', 'n', 'mu', 'method',
    'estimate', 'se',
    'lower', 'upper',
    'bias', 'covered',
    'runtime'
  )
]

library(dplyr)

results_summary <- results %>%
  group_by(n, mu, method) %>%
  summarise(
    nrep = n(),
    mean_estimate = mean(estimate),
    bias = mean(estimate - truth),
    empirical_sd = sd(estimate),
    mean_se = mean(se),
    rmse = sqrt(mean((estimate - truth)^2)),
    coverage = mean(covered),
    mean_interval_width = mean(interval_width),
    mean_runtime = mean(runtime),
    .groups = 'drop'
  )

results_summary
