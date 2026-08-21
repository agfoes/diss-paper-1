# The purpose of this script is to create analytic datasets from completed simulation results. This script requires minimal
# changes over each simulation and automatically creates and saves useful datasets for further analysis.

library(tidyverse)

################################ user input required #####

dir <- "/work/users/a/g/agfoes/P1"
sim_name <- "sampler_comp_07292026"
config_file_name <- "config_sampler_comparison"

################################ user input required #####

source(file.path(dir, "R", "helpers_summary.R"))
source(file.path(dir, "results", sim_name, paste0(config_file_name, ".R")))


## define result directories
result_dir <- file.path(dir, "results", config[["run_name"]])                   # result folder
summaries_dir <- file.path(result_dir, "summaries")                             # sub-folder within result_dir with posterior summaries               ## always saved
data_dir <- file.path(result_dir, "data")                                       # sub-folder within result_dir with generated data (may be empty)     ## non-empty with save_level >= 1
posterior_samples_dir <- file.path(result_dir, "posterior_samples")             # sub-folder within result_dir with posterior samples (may be empty)  ## non-empty with save_level >= 2
analysis_dir <- file.path(result_dir, "analysis")                               # name of folder to save produced datasets and figures
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)                # create analysis result folder if it doesn't already exist



################################################################################ posterior_results

summary_files <- list.files(
  path = summaries_dir,
  pattern = "^posterior_summaries_array_[0-9]+\\.csv$",
  recursive = TRUE, full.names = TRUE
)

if (length(summary_files) == 0) {
  stop("No posterior summary files were found in: ", summaries_dir)
}

summary_list <- vector(mode = "list", length = length(summary_files))

for (i in seq_along(summary_files)) {
  summary_list[[i]] <- read.csv(summary_files[i], stringsAsFactors = FALSE)
}

posterior_results <- do.call(rbind, summary_list)

## read simulation grid
sim_grid <- read.csv(file.path(result_dir, "sim_grid.csv"), stringsAsFactors = FALSE)

## combine posterior summary results and simulation scenarios
posterior_results <- posterior_results %>%
  left_join(sim_grid, by = "sim_id", relationship = "many-to-one")

## check for duplicate simulation scenarios or results with no simulation ID in sim_grid
if (any(is.na(posterior_results$array_id))) {
  unmatched_ids <- posterior_results %>%
    filter(is.na(array_id)) %>%
    distinct(sim_id) %>%
    pull(sim_id)
  }

duplicate_summaries <- posterior_results %>%
  count(sim_id, L, parameter) %>%
  filter(n > 1)

posterior_results <- posterior_results %>%
  distinct(sim_id, L, parameter, .keep_all = TRUE)

## check observed vs. expected fit results
expected_fits <- tidyr::crossing(sim_id = sim_grid$sim_id, L = config[["L_values"]])
observed_fits <- posterior_results %>% distinct(sim_id, L)
missing_fits <- expected_fits %>% anti_join(observed_fits, by = c("sim_id", "L"))
extra_fits <- observed_fits %>% anti_join(expected_fits, by = c("sim_id", "L"))

message(
  'Expected fits: ', nrow(expected_fits),
  '\nCompleted fits: ', nrow(observed_fits),
  '\nMissing fits: ', nrow(missing_fits),
  '\nUnexpected fits: ', nrow(extra_fits)
)

write.csv(
  posterior_results,
  file.path(analysis_dir, 'posterior_results.csv'),
  row.names = FALSE
)

saveRDS(
  posterior_results, 
  file.path(analysis_dir, "posterior_results.rds")
)

################################################################################ master_summary
master_summary <- posterior_results %>%
  group_by(cens, n, L, sampler, parameter) %>%
  summarize(
    n_rep = n_distinct(sim_id),
    
    mean_estimate = safe_mean(mean),
    bias = safe_mean(mean - truth),
    empirical_sd = safe_sd(mean),
    rmse = sqrt(safe_mean(squared_error)),
    
    coverage = safe_mean(covered),
    
    avg_posterior_sd = safe_mean(sd),
    
    avg_ess = safe_mean(ess),
    median_ess = median(ess),
    q1_ess = quantile(ess, 0.25, na.rm = TRUE),
    q3_ess = quantile(ess, 0.75, na.rm = TRUE),
    min_ess = safe_min(ess),
    
    avg_ess_sec = safe_mean(ess_sec),
    median_ess_sec = median(ess_sec, na.rm = TRUE),
    q1_ess_sec = quantile(ess_sec, 0.25, na.rm = TRUE),
    q3_ess_sec = quantile(ess_sec, 0.75, na.rm = TRUE),
    
    avg_mcse = safe_mean(mcse),
    
    .groups = "drop"
  ) %>%
  mutate(
    cens = factor(
      cens,
      levels = c(20, 40, 80),
      labels = c(
        "20% Censored",
        "40% Censored",
        "80% Censored"
      )
    ),
    n = factor(n),
    L = factor(L),
    
    coverage_mcse = sqrt(
      coverage * (1 - coverage) / n_rep
    ),
    sd_ratio = avg_posterior_sd / empirical_sd,
    stand_bias = bias / avg_posterior_sd,
    abs_stand_bias = abs(stand_bias)
  )

write.csv(
  master_summary,
  file.path(analysis_dir, 'master_summary.csv'),
  row.names = FALSE
)

saveRDS(
  master_summary,
  file.path(analysis_dir, "master_summary.rds")
)

################################################################################ data_summaries
data_summary_files <- list.files(
  summaries_dir, pattern = "^data_summaries_array_[0-9]+\\.csv", 
  full.names = TRUE
)

if (length(data_summary_files) == 0) {
  stop("No data summary files were found in: ", summaries_dir)
}

data_summaries <- data_summary_files %>%
  map_dfr(read.csv) %>%
  distinct(sim_id, .keep_all = TRUE) %>%
  left_join(
    sim_grid %>%
      select(sim_id, cens, n, rep, array_id),
    by = "sim_id",
    relationship = "many-to-one",
    suffix = c("", "_grid")
  )

write.csv(
  data_summaries,
  file.path(analysis_dir, "data_summaries.csv"),
  row.names = FALSE
)

saveRDS(
  data_summaries,
  file.path(analysis_dir, "data_summaries.rds")
)

################################################################################ censoring_summary
censoring_summary <- data_summaries %>%
  group_by(cens, n) %>%
  summarize(
    n_rep = n_distinct(sim_id),
    
    avg_prop_censored = safe_mean(prop_censored),
    avg_prop_left = safe_mean(prop_left),
    avg_prop_right = safe_mean(prop_right),
    avg_prop_interval = safe_mean(prop_interval),
    
    avg_prop_left_given_censored = safe_mean(prop_left_given_censored),
    avg_prop_interval_given_censored = safe_mean(prop_interval_given_censored),
    avg_prop_right_given_censored = safe_mean(prop_right_given_censored),
    
    avg_mean_interval_width = safe_mean(mean_interval_width),
    avg_median_interval_width = safe_mean(median_interval_width),
    
    .groups = "drop"
  ) %>%
  mutate(
    cens = factor(
      cens,
      levels = c(20, 40, 80),
      labels = c(
        "20% Censored",
        "40% Censored",
        "80% Censored"
      )
    ),
    n = factor(n)
  )

write.csv(
  censoring_summary,
  file.path(analysis_dir, "censoring_summary.csv"),
  row.names = FALSE
)

saveRDS(
  censoring_summary,
  file.path(analysis_dir, "censoring_summary.rds")
)

################################################################################ runtime_summary

runtime_summary <- posterior_results %>%
  distinct(
    sim_id,
    cens,
    n,
    L,
    sampler,
    runtime
  ) %>%
  group_by(
    cens,
    n,
    L,
    sampler
  ) %>%
  summarize(
    n_fit = n(),
    mean_runtime = mean(runtime),
    sd_runtime = sd(runtime),
    median_runtime = median(runtime),
    q1_runtime = quantile(runtime, 0.25),
    q3_runtime = quantile(runtime, 0.75),
    min_runtime = min(runtime),
    max_runtime = max(runtime),
    .groups = 'drop'
  ) %>%
  mutate(
    mean_runtime_hr = mean_runtime / 3600,
    sd_runtime_hr = sd_runtime / 3600,
    median_runtime_hr = median_runtime / 3600,
    q1_runtime_hr = q1_runtime / 3600,
    q3_runtime_hr = q3_runtime / 3600,
    min_runtime_hr = min_runtime / 3600,
    max_runtime_hr = max_runtime / 3600
  )

write.csv(
  runtime_summary,
  file.path(analysis_dir, "runtime_summary.csv"),
  row.names = FALSE
)

saveRDS(
  runtime_summary,
  file.path(analysis_dir, "runtime_summary.rds")
)

################################################################################ trace plots (if samples were saved)

trace_files <- list.files(
  posterior_samples_dir,
  pattern = "trace_samples.*\\.rds",
  full.names = TRUE
)

if (length(trace_files) > 0) {
  trace_list <- vector(mode = "list", 
                       length = length(trace_files))
  
  for (i in seq_along(trace_files)) {
    trace_file <- trace_files[i]
    file_name <- basename(trace_file)
    
    sim_id <- as.integer(
      sub(
        '^.*sim_([0-9]+)_L_.*$',
        '\\1',
        file_name
      )
    )
    
    L_value <- as.integer(
      sub(
        '^.*_L_([0-9]+)\\.rds$',
        '\\1',
        file_name
      )
    )
    
    samples <- readRDS(trace_file)
    
    available_parameters <- intersect(
      config[["key_vars"]],
      colnames(samples)
    )
    
    sample_matrix <- as.matrix(
      samples[, available_parameters, drop = FALSE]
      )
    sample_data <- as.data.frame(
      sample_matrix, check.names = FALSE
      ) %>%
      mutate(
        iteration = seq_len(n()),
        sim_id = sim_id,
        L = L_value
      ) %>%
      pivot_longer(
        cols = all_of(available_parameters),
        names_to = "parameter",
        values_to = "value"
      )

    trace_list[[i]] <- sample_data
  }
  
  if (length(trace_list) > 0) {
    
    trace_samples <- bind_rows(trace_list) %>%
      left_join(
        sim_grid %>%
          select(
            sim_id,
            cens,
            n,
          ) %>%
          distinct(),
        by = 'sim_id'
      ) %>%
      mutate(sampler = config[["sampler"]])
    
    saveRDS(
      trace_samples,
      file.path(analysis_dir, "trace_samples.rds")
    )
    
    write.csv(
      trace_samples, 
      file.path(analysis_dir, "trace_samples.csv"),
      row.names = FALSE
    )
  }
}
