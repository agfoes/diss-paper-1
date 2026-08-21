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
  
  stop("The following sim_id values were not found in sim_grid: ", paste(unmatched_ids, collapse = ", "))
}

duplicate_summaries <- posterior_results %>%
  count(sim_id, sampler, parameter) %>%
  filter(n > 1)

if (nrow(duplicate_summaries) > 0) {
  warning(nrow(duplicate_summaries), " duplicated sim_id/L/parameter combinations were found") 
}

posterior_results <- posterior_results %>%
  distinct(sim_id, sampler, parameter, .keep_all = TRUE)

## check observed vs. expected fit results
expected_fits <- tidyr::crossing(sim_id = sim_grid$sim_id, sampler = config[["samplers"]])
observed_fits <- posterior_results %>% distinct(sim_id, sampler)
missing_fits <- expected_fits %>% anti_join(observed_fits, by = c("sim_id", "sampler"))
extra_fits <- observed_fits %>% anti_join(expected_fits, by = c("sim_id", "sampler"))

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
  group_by(cens, n, sampler, parameter) %>%
  summarize(
    n_rep = n_distinct(sim_id),
    
    mean_estimate = safe_mean(mean),
    bias = safe_mean(mean - truth),
    empirical_sd = safe_sd(mean),
    rmse = sqrt(safe_mean(squared_error)),
    
    coverage = safe_mean(covered),
    
    avg_posterior_sd = safe_mean(sd),
    avg_ess = safe_mean(ess),
    min_ess = safe_min(ess),
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
    sampler = as.character(sampler),
    coverage_mcse = sqrt(
      coverage * (1 - coverage) / n_rep
    )
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
