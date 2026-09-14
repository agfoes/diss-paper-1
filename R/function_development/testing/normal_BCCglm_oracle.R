
library(dplyr)
library(nimble)
library(nimbleNoBounds)
library(nimbleMacros)

main_folder <- "/work/users/a/g/agfoes/P1"

source(file.path(main_folder, "R/function_development/BCCglm.R"))
source(file.path(main_folder, "R/function_development/buildModelCode.R"))
source(file.path(main_folder, "R/function_development/customSamplers.R"))
source(file.path(main_folder, "R/function_development/dg_glm.R"))
source(file.path(main_folder, "R/function_development/dynamicglm.R"))
source(file.path(main_folder, "R/function_development/nimbleInputs.R"))

sim_grid <- expand.grid(n = c(100, 500, 1000),
                        cens = c("0.2", "0.4", "0.8"))

slurm_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

set.seed(as.numeric(09022026 + as.numeric(slurm_id)))

for (i in 1:9) {
  n <- sim_grid[i, "n"]
  cens <- sim_grid[i, "cens"]
  
  set.seed(as.numeric(09022026 + as.numeric(slurm_id) + as.numeric(n)))
  
  data <- dg_glm(n = n,
                 cens = cens,
                 family = gaussian(link = 'identity'),
                 beta = c(1, 0.7, -1, 1) # truth in order of output: (0.7, 1, -1, 1)
                 )
  fit_data <- data$data
  oracle_data <- data$data %>%
    mutate(X = data$X)
  
  fit_time <- system.time({
    fit <- BCCglm(data = fit_data,
                  outcome_formula = Y ~ X + Z2 + Z3, # default adds intercept
                  family = gaussian(link = 'identity'),
                  covariate_formula = X ~ Z2 + Z3, # default adds intercept
                  censored_covariate = "X",
                  censoring_bounds = c("CL", "CR"))
  })[["elapsed"]]
  
  beta_summary <- summary(fit$samples)$statistics[1:4, ]
  
  oracle_time <- system.time({
    oracle_fit <- glm(Y ~ X + Z2 + Z3, data = oracle_data,
                      family = gaussian(link = 'identity'))
  })[["elapsed"]]
  
  oracle_summary <- as.data.frame(summary(oracle_fit)$coefficients)
  oracle_summary$term <- rownames(oracle_summary)
  rownames(oracle_summary) <- NULL
  
  oracle_summary <- oracle_summary[
    , c('term', 'Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
  ]
  
  # runtime results
  runtimes <- data.frame(
    method = c("oracle", "bccglm"),
    runtime = c(
      oracle_time,
      fit_time
    )
  )
  
  # save res
  write.csv(fit$samples,
            file = file.path(main_folder, "results/normal_BCCglm_oracle/samples", 
                             paste0("data_n_", n, 
                                    "_cens_", as.numeric(cens), 
                                    "_array_", as.numeric(slurm_id),
                                    ".csv")),
            row.names = FALSE)
  write.csv(beta_summary,
            file = file.path(main_folder, "results/normal_BCCglm_oracle/summaries", 
                             paste0("data_n_", n, 
                                    "_cens_", as.numeric(cens), 
                                    "_array_", as.numeric(slurm_id),
                                    ".csv")),
            row.names = FALSE)
  write.csv(oracle_summary,
            file = file.path(main_folder, "results/normal_BCCglm_oracle/summaries", 
                             paste0("data_n_", n, 
                                    "_cens_", as.numeric(cens), 
                                    "_array_", as.numeric(slurm_id),
                                    ".csv")),
            row.names = FALSE)
  write.csv(runtimes,
            file = file.path(main_folder, "results/normal_BCCglm_oracle/runtimes", 
                             paste0("data_n_", n, 
                                    "_cens_", as.numeric(cens), 
                                    "_array_", as.numeric(slurm_id),
                                    ".csv")),
            row.names = FALSE)
  
  # save data
  saveRDS(data, file = file.path(main_folder, "results/normal_BCCglm_oracle/data", 
                                 paste0("data_n_", n, 
                                        "_cens_", as.numeric(cens), 
                                        "_array_", as.numeric(slurm_id),
                                        ".rds")))
}
