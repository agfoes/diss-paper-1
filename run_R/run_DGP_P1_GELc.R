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

if (task_id <= 500) {
  rep <- task_id
  model = "gelc"
} else if (task_id <= 1000) {
  rep <- task_id - 500
  model = "p1"
}

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
  res <- as.data.frame(p1_summary) %>%
    tibble::rownames_to_column("parameter") %>%
    mutate(runtime = as.numeric(p1_time),
           rep = rep)
}


write.csv(res, 
          file = file.path(results_folder, paste0("results_", task_id, ".csv")),
          row.names = FALSE)

