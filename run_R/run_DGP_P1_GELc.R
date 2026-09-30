# libraries ----
library(dplyr)
library(nimble)
library(nimbleNoBounds)
library(nimbleMacros)
library(tibble)
library(ICenCov)

main_folder <- "/work/users/a/g/agfoes/P1"
results_folder <- file.path(main_folder, "results", "sims_09232026")

source(file.path(main_folder, "R/function_development/BCCglm.R"))
source(file.path(main_folder, "R/function_development/buildModelCode.R"))
source(file.path(main_folder, "R/function_development/customSamplers.R"))
source(file.path(main_folder, "R/function_development/dg_glm.R"))
source(file.path(main_folder, "R/function_development/dynamicglm.R"))
source(file.path(main_folder, "R/function_development/nimbleInputs.R"))


task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#if (task_id <= 500) {
#  rep <- task_id
#  model = "gelc"
#} else if (task_id <= 1000) {
#  rep <- task_id - 500
#  model = "p1"
#}

rep <- task_id
model = "p1"

set.seed(09232026 + 1000 * rep)

data <- dg_glm(n = 1000,
               outcome_formula = Y ~ X + Z2 + Z3,
               covariate_formula = X ~ Z2 + Z3,
               family = gaussian(link = "identity"),
               beta = c(1, 0.7, 1, -1),
               censoring = 0.80,
               z_spec = list(
                 Z2 = list(dist = 'bernoulli',
                           prob = 0.8),
                 Z3 = list(dist = 'normal',
                           mean = 2, sd = 1)
               ))

censor_type <- case_when(
  is.infinite(data$CL) ~ "left",
  is.infinite(data$CR) ~ "right",
  data$CL == data$CR ~ "observed",
  TRUE ~ "interval"
)

width <- data$CR[censor_type == "interval"] - 
         data$CL[censor_type == "interval"]

data_summary <- tibble(
  rep = rep,
  
  prop_observed = mean(censor_type == "observed"),
  prop_interval = mean(censor_type == "interval"),
  prop_left = mean(censor_type == "left"),
  prop_right = mean(censor_type == "right"),
  
  mean_interval_width = mean(width),
  median_interval_width = median(width),
  q90_interval_width = quantile(width, 0.90),
  
  x_mean = mean(data$X),
  x_sd = sd(data$X),
  x_min = min(data$X),
  x_q01 = quantile(data$X, 0.01),
  x_q99 = quantile(data$X, 0.99),
  x_max = max(data$X)
)

write.csv(data_summary,
          file = file.path(results_folder, "datasummary", paste0("data_summary_", rep, ".csv")),
          row.names = FALSE)


data$CL[is.infinite(data$CL) & data$CL < 0] <- -50
data$CR[is.infinite(data$CR) & data$CR > 0] <- 50

data$data$CL <- data$CL
data$data$CR <- data$CR

dat <- data$data


if (model == "gelc") {
  # gelc
  gelc_time <- system.time({
    
    fit_gelc <- icglm(
      Y ~ Z2 + Z3 + ic(CL, CR, 'X'),
      family = gaussian(link = "identity"),
      data = dat
    )
    
  })['elapsed']
  
  gelc_summary <- summary(fit_gelc)
  gelc_coef <- gelc_summary$coefficients
  res <- as.data.frame(gelc_coef) %>%
    rownames_to_column("parameter") %>%
    mutate(runtime = as.numeric(gelc_time),
           rep = rep)
  
} else {
  # p1
  p1_time <- system.time({
    
    fit <- BCCglm(data = dat,
                  outcome_formula = Y ~ Z2 + Z3 + X,
                  family = gaussian(link = 'identity'),
                  covariate_formula = X ~ Z2 + Z3,
                  censored_covariate = "X",
                  censoring_bounds = c("CL", "CR"))
    
  })['elapsed']
  
  p1_summary <- summary(fit$samples)$statistics[1:4, ]
  
  ess <- coda::effectiveSize(coda::as.mcmc(fit$samples[, 1:4]))
  
  res <- as.data.frame(p1_summary) %>%
    tibble::rownames_to_column("parameter") %>%
    mutate(runtime = as.numeric(p1_time),
           rep = rep,
           ESS = as.numeric(ess),
           ESS_per_sec = ESS / as.numeric(p1_time),
           rel_MCSE = `Time-series SE` / SD)
  
  beta_names <- paste0("beta[", 1:4, "]")
  saveRDS(
    fit$samples[, beta_names],
    file = file.path(results_folder, "samples", paste0("beta_samples_", rep, ".rds"))
  )
}


write.csv(res, 
          file = file.path(results_folder, "rawres2", paste0("results_", task_id, ".csv")),
          row.names = FALSE)

