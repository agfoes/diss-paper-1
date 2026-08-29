library(nimble)
library(nimbleNoBounds)
library(coda)
library(tidyverse)
library(ICenCov)

dir <- '/work/users/a/g/agfoes/P1'

source(
  file.path(
    dir,
    'R',
    'config_gelc_comparison.R'
  )
)

source(
  file.path(
    dir,
    'R',
    'helpers_data.R'
  )
)

source(
  file.path(
    dir,
    'R',
    'helpers_nimble.R'
  )
)

source(
  file.path(
    dir,
    'R',
    'helpers_sampler.R'
  )
)

includeObs <- TRUE
includeCen <- TRUE

joint_model_code <- source(
  file.path(
    dir,
    'R',
    'model_code',
    'model_normal_joint_betas.R'
  )
)$value

result_dir <- file.path(
  config[['project_dir']],
  'results',
  config[['run_name']]
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

rep <- as.integer(
  Sys.getenv('SLURM_ARRAY_TASK_ID')
)

if(is.na(rep)) {
  rep <- 1L
}

result_file <- file.path(
  summary_dir,
  sprintf(
    'results_rep_%04d.csv',
    rep
  )
)

if(file.exists(result_file)) {
  message(
    'Replication ',
    rep,
    ' already completed.'
  )
  quit(save = 'no')
}

## ------------------------------------------------------------
## Generate complete data ONCE
## ------------------------------------------------------------

set.seed(
  config[['sim_seed']] +
    rep
)

data_full <- datagen_gelc_normal(
  n = max(config[['n_values']]),
  beta_x = config[['beta_x']],
  beta_z = config[['beta_z']],
  beta_0 = config[['beta_0']],
  tau = config[['tau']]
)

saveRDS(
  data_full,
  file.path(
    data_dir,
    sprintf(
      'complete_data_rep_%04d.rds',
      rep
    )
  )
)

rep_results <- list()

row_id <- 0

## ------------------------------------------------------------
## Generate censoring for each mu
## ------------------------------------------------------------

for(mu in config[['mu_values']]) {
  
  ## separate deterministic censoring seed
  set.seed(
    config[['sim_seed']] +
      1000000 +
      rep * 100 +
      mu
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
    
    data <- data_mu[
      seq_len(n),
      ,
      drop = FALSE
    ]
    
    ## ========================================================
    ## oracle
    ## ========================================================
    
    oracle_time <- system.time({
      
      fit_oracle <- lm(
        Y ~ Z2 + X,
        data = data
      )
      
    })[['elapsed']]
    
    oracle_coef <-
      summary(fit_oracle)$coefficients
    
    row_id <- row_id + 1
    
    rep_results[[row_id]] <- data.frame(
      rep = rep,
      n = n,
      mu = mu,
      method = 'oracle',
      parameter = 'beta_x',
      truth = config[['beta_x']],
      estimate = oracle_coef[
        'X',
        'Estimate'
      ],
      se = oracle_coef[
        'X',
        'Std. Error'
      ],
      lower =
        oracle_coef['X', 'Estimate'] -
        1.96 *
        oracle_coef['X', 'Std. Error'],
      upper =
        oracle_coef['X', 'Estimate'] +
        1.96 *
        oracle_coef['X', 'Std. Error'],
      runtime = as.numeric(oracle_time)
    )
    
    ## ========================================================
    ## GELc
    ## ========================================================
    
    gelc_error <- FALSE
    
    gelc_time <- system.time({
      
      fit_gelc <- try(
        icglm(
          Y ~ Z2 + ic(CL, CR, 'X'),
          family = gaussian,
          data = data
        ),
        silent = TRUE
      )
      
    })[['elapsed']]
    
    if(inherits(fit_gelc, 'try-error')) {
      
      gelc_error <- TRUE
      
      row_id <- row_id + 1
      
      rep_results[[row_id]] <- data.frame(
        rep = rep,
        n = n,
        mu = mu,
        method = 'GELc',
        parameter = 'beta_x',
        truth = config[['beta_x']],
        estimate = NA_real_,
        se = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        runtime = as.numeric(gelc_time)
      )
      
    } else {
      
      gelc_summary <- summary(fit_gelc)
      
      gelc_coef <-
        gelc_summary$coefficients
      
      row_id <- row_id + 1
      
      rep_results[[row_id]] <- data.frame(
        rep = rep,
        n = n,
        mu = mu,
        method = 'GELc',
        parameter = 'beta_x',
        truth = config[['beta_x']],
        estimate =
          gelc_coef[
            'X',
            'Estimate'
          ],
        se =
          gelc_coef[
            'X',
            'Std. Error'
          ],
        lower =
          gelc_coef['X', 'Estimate'] -
          1.96 *
          gelc_coef['X', 'Std. Error'],
        upper =
          gelc_coef['X', 'Estimate'] +
          1.96 *
          gelc_coef['X', 'Std. Error'],
        runtime = as.numeric(gelc_time)
      )
    }
    
    ## ========================================================
    ## P1
    ## ========================================================
    
    pz <- config[['pz']]
    L <- config[['L']]
    
    Z <- as.matrix(
      data[, c('Z1', 'Z2')]
    )
    
    idx_obs <- which(
      data$Dobs == 1
    )
    
    idx_cen <- which(
      data$Dobs == 0
    )
    
    nobs <- length(idx_obs)
    ncen <- length(idx_cen)
    
    if(length(idx_obs) < 2) {
      idx_obs <- c(
        idx_obs,
        rep(
          1L,
          2 - length(idx_obs)
        )
      )
    }
    
    if(length(idx_cen) < 2) {
      idx_cen <- c(
        idx_cen,
        rep(
          1L,
          2 - length(idx_cen)
        )
      )
    }
    
    constants <- list(
      L = L,
      p = pz,
      mu_gamma =
        config[['mu_gamma']],
      nobs = nobs,
      ncen = ncen,
      idx_obs = idx_obs,
      idx_cen = idx_cen,
      beta_mean =
        rep(0, pz + 1),
      beta_cov =
        diag(pz + 1)
    )
    
    Ndata <- list(
      y_obs =
        data$Y[
          data$Dobs == 1
        ],
      
      y_cen =
        data$Y[
          data$Dobs == 0
        ],
      
      z_obs =
        Z[
          data$Dobs == 1,
          ,
          drop = FALSE
        ],
      
      z_cen =
        Z[
          data$Dobs == 0,
          ,
          drop = FALSE
        ],
      
      x_obs =
        data$X[
          data$Dobs == 1
        ],
      
      CL =
        data$CL[
          data$Dobs == 0
        ],
      
      CR =
        data$CR[
          data$Dobs == 0
        ],
      
      constraint_data =
        rep(1, ncen)
    )
    
    Ninits <- nimble_inits(
      data = data,
      pz = pz,
      constants = constants,
      Ndata = Ndata
    )
    
    p1_error <- FALSE
    
    p1_time <- system.time({
      
      p1_result <- try({
        
        model <- nimbleModel(
          code = joint_model_code,
          constants = constants,
          data = Ndata,
          inits = Ninits
        )
        
        if(!is.finite(model$calculate())) {
          stop(
            'Non-finite initial model calculation'
          )
        }
        
        set.seed(
          config[['sim_seed']] +
            2000000 +
            rep * 10000 +
            n * 10 +
            mu
        )
        
        run_sampler(
          model = model,
          sampler_type =
            config[['sampler']],
          pz = pz,
          L = L,
          niter =
            config[['niter']],
          nburnin =
            config[['nburnin']],
          thin =
            config[['thin']],
          vars_to_monitor =
            config[['key_vars']],
          ncen = ncen,
          nobs = nobs
        )
        
      }, silent = TRUE)
      
    })[['elapsed']]
    
    if(inherits(p1_result, 'try-error')) {
      
      p1_error <- TRUE
      
      row_id <- row_id + 1
      
      rep_results[[row_id]] <- data.frame(
        rep = rep,
        n = n,
        mu = mu,
        method = 'P1',
        parameter = 'beta_x',
        truth = config[['beta_x']],
        estimate = NA_real_,
        se = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        runtime = as.numeric(p1_time)
      )
      
    } else {
      
      beta_x_samples <-
        p1_result$samples[
          ,
          'beta[1]'
        ]
      
      beta_x_hpd <- HPDinterval(
        as.mcmc(beta_x_samples),
        prob = 0.95
      )
      
      row_id <- row_id + 1
      
      rep_results[[row_id]] <- data.frame(
        rep = rep,
        n = n,
        mu = mu,
        method = 'P1',
        parameter = 'beta_x',
        truth = config[['beta_x']],
        estimate =
          mean(beta_x_samples),
        se =
          sd(beta_x_samples),
        lower =
          beta_x_hpd[1, 'lower'],
        upper =
          beta_x_hpd[1, 'upper'],
        runtime =
          as.numeric(p1_time)
      )
    }
    
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

message(
  'Completed replication ',
  rep
)