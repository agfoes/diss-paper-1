# the purpose of this script is to test GELc model performance on increasing interval widths beyond what their published results include. P1 is currently built to have arbitrarily large censoring intervals whereas GELc appears to need finite (and beyond that, small-ish) censoring intervals. I want to quantify how "small-ish" these intervals need ot be and if there's a practical maximum interval width that their current implementation can handle.

# load packages and libraries ####
library(tidyverse)
library(ICenCov)

# load GELc functions and data generation functions ####
dir <- '/work/users/a/g/agfoes/P1'
source(file.path(dir, 'R', 'helpers', 'helpers_data.R'))

# slurm id and define results file path ####
rep <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

if(is.na(rep)) {
  rep <- 1L
}

results_file <- file.path(dir, "results", "gelc_stress_test", paste0("raw_results2/results_", rep, ".csv"))

# generate largest dataset once at max n value with zero censoring ####
set.seed(as.numeric(09042026 + rep*1000))

data_full <- datagen_gelc_gamma(
  n = 500,
  mu = 0,
  gamma = 0.02,
  phi = 0.02
)

rep_results <- list()

row_id <- 0

sim_sett <- expand.grid(n = c(100, 300, 500),
                        mu = c(0, 1, 3, 6, 9, 12, 20, 50, 100, 500),
                        array_id = c(1:1000)) %>%
  filter(array_id == rep)


# loop over sample sizes and censoring levels by modifying original dataset ####
for (i in 1:nrow(sim_sett)) {
  
  # get iteration values for n and mu
  n <- sim_sett$n[i]
  mu <- sim_sett$mu[i]
  
  # restrict data to first n rows
  subdata <- data_full[1:n, ]
  
  # apply censoring based on width mu
  cens_subdata <- gelc_censoring(mu = mu,
                                 n = n,
                                 Z = subdata$Z)
  data <- subdata %>%
    mutate(CL = cens_subdata$CL,
           CR = cens_subdata$CR)
  
  # fit GELc model
  fit <- tryCatch(
    gelc_time <- system.time({
      fit_gelc <- icglm(Y ~ -1 + X1 + ic(CL, CR, "Z"), 
                        family = Gamma(link = "log"),
                        data = data)
    })[["elapsed"]],
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    results_row <- data.frame(method = "GELc",
                              n = n,
                              mu = mu,
                              mean_int_width = mean(data$CR - data$CL),
                              success = FALSE, 
                              error = conditionMessage(fit),
                              runtime = NA,
                              
                              beta_0_est = NA,
                              beta_Z_est = NA,
                              
                              beta_0_sd = NA,
                              beta_Z_sd = NA)
  } else {
    coef_fit <- summary(fit_gelc)$coefficients
    
    results_row <- data.frame(method = "GELc",
                              n = n,
                              mu = mu,
                              mean_int_width = mean(data$CR - data$CL),
                              success = TRUE, 
                              error = NA,
                              runtime = gelc_time,
                              
                              beta_0_est = coef_fit["X1", "Estimate"],
                              beta_Z_est = coef_fit["Z", "Estimate"],
                              
                              beta_0_sd = coef_fit["X1", "Std. Error"],
                              beta_Z_sd = coef_fit["Z", "Std. Error"])
  }
  
  row_id <- row_id + 1
  
  rep_results[[row_id]] <- results_row
  
  ## fit oracle and save results
  oracle_time <- system.time({
    fit_oracle <- glm(Y ~ -1 + X1 + Z,
                    family = Gamma(link = "log"),
                    data = data)
  })[["elapsed"]]
  
  coef_oracle <- summary(fit_oracle)$coefficients
  
  results_row <- data.frame(method = "oracle",
                            n = n,
                            mu = mu,
                            mean_int_width = mean(data$CR - data$CL),
                            success = TRUE, 
                            error = NA,
                            runtime = oracle_time,
                            
                            beta_0_est = coef_oracle["X1", "Estimate"],
                            beta_Z_est = coef_fit["Z", "Estimate"],
                            
                            beta_0_sd = coef_oracle["X1", "Std. Error"],
                            beta_Z_sd = coef_oracle["Z", "Std. Error"])
  row_id <- row_id + 1
  rep_results[[row_id]] <- results_row
 
  
  save_results <- bind_rows(
    rep_results
  ) %>%
    mutate(rep = rep)
  
  write.table(
    save_results,
    file = results_file,
    sep = ",",
    append = file.exists(results_file),
    row.names = FALSE,
    col.names = !file.exists(results_file)
  )
   
}



