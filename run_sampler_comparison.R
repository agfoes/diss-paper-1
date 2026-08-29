
# This will be the only .R file called by the slurm_outline.sh script. This script 
# will call other static and dynamic .R files.

## load libraries
library(nimble)
library(truncnorm)
library(coda)
library(tidyverse)
library(nimbleNoBounds)

## source helper functions
dir <- "/work/users/a/g/agfoes/P1"
source(paste0(dir, "/R/sampler_comparison/config_sampler_comparison.R"))                # config
source(paste0(dir, "/R/helpers/helpers_data.R"))                                        # find_L, datagen
source(paste0(dir, "/R/helpers/helpers_nimble.R"))                                      # nimble_constants, nimble_data_years, nimble_inits
source(paste0(dir, "/R/helpers/helpers_sampler.R"))                                     # make_conf, run_sampler
source(paste0(dir, "/R/helpers/helpers_summary.R"))                                     # posterior_summary, data_censoring_summary
source(paste0(dir, "/R/model_code/model_normal_joint_betas.R"))                         # model_code
source(paste0(dir, "/R/custom_block_sampler.R"))
 
## read simulation setting file
result_dir <- file.path(config[["project_dir"]], "results", config[["run_name"]])
sim_grid <- read.csv(file = file.path(result_dir, "sim_grid.csv"))

## filter simulation settings to those matching this array's ID
idx <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(idx)) {idx = 1}
sim_grid_filt <- sim_grid %>% filter(array_id == idx)

## stop order if no simulation rows are assigned to this array
if (nrow(sim_grid_filt) == 0) {
  stop("No simulation rows found for array_id = ", idx)
}

## read/generate directories for saving results
data_dir <- file.path(result_dir, "data")
summary_dir <- file.path(result_dir, "summaries")
sample_dir <- file.path(result_dir, "posterior_samples")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

summary_file <- file.path(summary_dir, sprintf("posterior_summaries_array_%03d.csv", idx))
key_file <- file.path(sample_dir, sprintf("key_variables_array_%03d.csv", idx))
completion_file <- file.path(summary_dir, sprintf("completed_fits_array_%03d.csv", idx))
data_summary_file <- file.path(summary_dir, sprintf("data_summaries_array_%03d.csv", idx))

## check for if this set of simulations has already been started and needs completion
if (file.exists(completion_file)) {
  completed_fits <- read.csv(completion_file) %>%
    distinct(sim_id, sampler)
} else {
  completed_fits <- data.frame(
    sim_id = integer(),
    sampler = character()
  )
}

if (file.exists(data_summary_file)) {
  completed_data_summaries <- read.csv(data_summary_file) %>%
    distinct(sim_id)
} else {
  completed_data_summaries <- data.frame(
    sim_id = integer()
  )
}

## outer loop over simulation scenarios
for (i in seq_len(nrow(sim_grid_filt))) {
  
  ## get sim id and save level
  sim_id <- sim_grid_filt[["sim_id"]][i]
  save_level <- sim_grid_filt[["save_level"]][i]
  
  ## define file path for saving generated data
  data_file <- file.path(data_dir, sprintf("data_sim_%05d.rds", sim_id))
  
  ## get parameter values for planned censoring level
  alpha_start <- unlist(sim_grid_filt[i, c(
    "alpha_start_1",
    "alpha_start_2",
    "alpha_start_3"
  )])
  alpha_width <- unlist(sim_grid_filt[i, c(
    "alpha_width_1",
    "alpha_width_2",
    "alpha_width_3"
  )])
  
  ## if there are already some results saved, read in the associated dataset, otherwise set seed and generate data
  set.seed(as.numeric(config[["sim_seed"]]) + sim_id)
  if (file.exists(data_file)) {
    data <- readRDS(data_file)
  } else {
    data <- datagen(n = sim_grid_filt[["n"]][i],
                    p_z2 = 0.8, 
                    mu_z3 = 2, 
                    sd_z3 = 1,
                    comps = 3, 
                    w = c(0.3, 0.5, 0.2), 
                    gamma = matrix(c(180, 730, 550,
                                     255, 1095, 730,
                                     330, 1460, 1095),
                                   nrow = 3, byrow = TRUE), 
                    x_sd = c(365, 455, 550),
                    beta_x = 0.7/365.25, 
                    beta_z = c(0,1,-1), 
                    tau = 1,
                    beta_omega = matrix(c(0, 1.5, -0.5,
                                          0, 0.5, 0),
                                        nrow = 2, byrow = TRUE), 
                    alpha_start = alpha_start, 
                    alpha_width = alpha_width,
                    sigma_start = 30, 
                    sigma_width = 0.15,
                    rho = c(0.85, 0.9, 0.99), 
                    jitter_sd = 30, 
                    eps_obs = sim_grid_filt[["eps_obs"]][i])
    saveRDS(data, data_file)
  }
  
  ## summarize actual censoring rates and proportions, interval widths
  if (!(sim_id %in% completed_data_summaries$sim_id)) {
    data_summary <- data_censoring_summary(data = data, sim_id = sim_id)
    
    write.table(
      data_summary,
      file = data_summary_file,
      sep = ",",
      row.names = FALSE,
      col.names = !file.exists(data_summary_file),
      append = file.exists(data_summary_file)
    )
    
    completed_data_summaries <- bind_rows(completed_data_summaries, data.frame(sim_id = sim_id))
  }
  
  ## get nimble constants, data, inits
  pz <- config[["pz"]]
  sigma_bx <- config[["sigma_bx"]]
  sigma_bz <- config[["sigma_bz"]]
  mu_gamma <- config[["mu_gamma"]]
  
  Nconstants <- nimble_constants(data,
                                 L = config[["L"]], pz = pz,
                                 sigma_bx = sigma_bx,
                                 sigma_bz = sigma_bz,
                                 mu_gamma = mu_gamma)
  Ndata <- nimble_data_years(data,
                             constants = Nconstants)
  Ninits <- nimble_inits(data,
                         pz = pz,
                         constants = Nconstants,
                         Ndata = Ndata)

  ## inner loop over samplers with same dataset
  for (sampler in config[["samplers"]]) {
    
    ## if any samplers have already been fit, load and skip this value
    if (any(
      completed_fits$sim_id == sim_id & completed_fits$sampler == sampler)) {
      message("Skipping completed fit: sim_id = ",  sim_id, ", sampler = ", sampler)
      next
    }
    
    ## fit model and get samples
    model <- nimbleModel(
      code = model_code,
      data = Ndata,
      constants = Nconstants,
      inits = Ninits
    )
    
    vars_to_monitor <- if (save_level >= 2) {
      config[["all_vars"]]
    } else {
      config[["key_vars"]]
    }

    ## set seed again
    sampler_id <- match(sampler, config[["samplers"]])
    set.seed(as.numeric(config[["sim_seed"]]) + 1000000 + sim_id*1000 + sampler_id)
    
    fit <- run_sampler(model,
                       sampler_type = as.character(sampler),
                       pz = pz,
                       L = config[["L"]],
                       niter = config[["niter"]],
                       nburnin = config[["nburnin"]],
                       thin = config[["thin"]],
                       vars_to_monitor = vars_to_monitor,
                       nobs = Nconstants[["nobs"]], 
                       ncen = Nconstants[["ncen"]]
                       )
    
    ## get results and summaries
    posterior_summaries <- posterior_summary(
      samples = fit$samples,
      sampler_name = as.character(sampler),
      truth = config[["truth"]],
      key_vars = config[["key_vars"]],
      runtime = fit$runtime
    ) %>%
      mutate(
        sim_id = sim_id,
        .before = 1
      )
    
    ## always store posterior summaries
    write.table(
      posterior_summaries,
      file = summary_file,
      sep = ",",
      row.names = FALSE,
      col.names = !file.exists(summary_file),
      append = file.exists(summary_file)
    )
    
    ## save thinned key variable posterior samples
    if (save_level == 1) {
      
      saved_thin_samples <- seq(1, nrow(fit$samples), by = config[["save_key_thinning"]])
      
      key_samples <- as.data.frame(
        fit$samples[saved_thin_samples, 
                    c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "tau"), drop = FALSE]
      ) %>%
        mutate(
          sim_id = sim_id,
          sampler = as.character(sampler),
          sample_index = saved_thin_samples,
          .before = 1
        )
      
      write.table(
        key_samples,
        file = key_file,
        sep = ",",
        row.names = FALSE,
        col.names = !file.exists(key_file),
        append = file.exists(key_file)
      )
    }
    
    ## save thinned posterior samples for all variables
    if (save_level == 2) {
      
      saved_thin_samples <- seq(1, nrow(fit$samples), by = config[["save_full_thinning"]])
      
      full_samples <- as.data.frame(
        fit$samples[saved_thin_samples, , drop = FALSE]
      ) %>%
        mutate(
          sim_id = sim_id,
          sampler = as.character(sampler),
          iteration = saved_thin_samples,
          .before = 1
        )
      
      full_file <- file.path(
        sample_dir,
        sprintf("all_variables_array_%03d_sampler_%s.csv", idx, as.character(sampler))
      )
      
      write.table(
        full_samples,
        file = full_file,
        sep = ",",
        row.names = FALSE,
        col.names = !file.exists(full_file),
        append = file.exists(full_file)
      )
    }
    
    ## save full posterior samples for select sims when indicated
    if (save_level == 3) {
      
      trace_file <- file.path(sample_dir, sprintf("trace_samples_sim_%05d_sampler_%s.rds", sim_id, as.character(sampler)))
      saveRDS(fit$samples, file = trace_file)
    }
    
    ## mark fit completed
    completed_row <- data.frame(sim_id = sim_id, sampler = sampler)
    write.table(
      completed_row, file = completion_file,
      sep = ",", row.names = FALSE,
      col.names = !file.exists(completion_file),
      append = file.exists(completion_file)
    )
    completed_fits <- bind_rows(
      completed_fits,
      completed_row
    )
   
    ## data cleanup to free up memory after each sampler fit
    rm(fit, model, posterior_summaries)
    if (exists("key_samples")) {rm(key_samples)}
    if (exists("full_samples")) {rm(full_samples)}
    gc()
     
  }
  
  ## data cleanup to free up memory after each data (over all sampler values) fit
  rm(data, Nconstants, Ndata, Ninits)
  if (exists("data_summary")) {
    rm(data_summary)
  }
  gc()

}

