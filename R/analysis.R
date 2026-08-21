# This script will be called after simulations have been run to analyze results.

## load libraries
library(tidyverse)
library(ggplot2)

## source helper functions and simulation information
dir <- "/work/users/a/g/agfoes/P1"
source(paste0(dir, "/results/sims_07152026/config_07152026.R"))                 # config

## define result directories
result_dir <- file.path(dir, "results", config[["run_name"]])                   # result folder
summaries_dir <- file.path(result_dir, "summaries")                             # sub-folder within result_dir with posterior summaries               ## always saved
data_dir <- file.path(result_dir, "data")                                       # sub-folder within result_dir with generated data (may be empty)     ## non-empty with save_level >= 1
posterior_samples_dir <- file.path(result_dir, "posterior_samples")             # sub-folder within result_dir with posterior samples (may be empty)  ## non-empty with save_level >= 2
analysis_dir <- file.path(result_dir, "analysis")                               # name of folder to save produced datasets and figures
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)                # create analysis result folder if it doesn't already exist

## read in posterior summary results
summary_files <- list.files(
  path = summaries_dir,
  pattern = "^posterior_summaries_array_[0-9]+\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)
summary_list <- vector(mode = "list", length = length(summary_files))

for (i in seq_along(summary_files)) {
  summary_list[[i]] <- read.csv(summary_files[i], stringsAsFactors = FALSE)
}

summary_results <- do.call(rbind, summary_list)

## read simulation grid
sim_grid <- read.csv(file.path(result_dir, "sim_grid.csv"), stringsAsFactors = FALSE)

## combine posterior summary results and simulation scenarios
summary_results <- summary_results %>%
  left_join(sim_grid, by = "sim_id", relationship = "many-to-one")

## check for duplicate simulation scenarios or results with no simulation ID in sim_grid
if (any(is.na(summary_results$array_id))) {
  unmatched_ids <- summary_results %>%
    filter(is.na(array_id)) %>%
    distinct(sim_id) %>%
    pull(sim_id)
  
  stop("The following sim_id values were not found in sim_grid: ", paste(unmatched_ids, collapse = ", "))
}

duplicate_summaries <- summary_results %>%
  count(sim_id, L, parameter) %>%
  filter(n > 1)

if (nrow(duplicate_summaries > 0)) {
 warning(nrow(duplicate_summaries), " duplicated sim_id/L/parameter combinations were found") 
}

summary_results <- summary_results %>%
  distinct(sim_id, L, parameter, .keep_all = TRUE)

## check observed vs. expected fit results
expected_fits <- tidyr::crossing(sim_id = sim_grid$sim_id, L = config[["L_values"]])
observed_fits <- summary_results %>% distinct(sim_id, L)
missing_fits <- expected_fits %>% anti_join(observed_fits, by = c("sim_id", "L"))
extra_fits <- observed_fits %>% anti_join(expected_fits, by = c("sim_id", "L"))

message(
  'Expected fits: ', nrow(expected_fits),
  '\nCompleted fits: ', nrow(observed_fits),
  '\nMissing fits: ', nrow(missing_fits),
  '\nUnexpected fits: ', nrow(extra_fits)
)

## posterior summary result analysis
master_summary <- summary_results %>%
  group_by(cens, n, L, sampler, parameter) %>%
  summarize(
    n_rep = n_distinct(sim_id),
    
    mean_estimate = mean(mean),
    bias = mean(mean - truth),
    empirical_sd = sd(mean),
    rmse = sqrt(mean(squared_error)),
    
    coverage = mean(covered),
    
    avg_posterior_sd = mean(sd),
    avg_ess = mean(ess),
    min_ess = min(ess),
    avg_mcse = mean(mcse),
    
    .groups = "drop"
  )  %>%
  mutate(
    cens = factor(
      cens,
      levels = c(20, 40, 80),
      labels = c("20% Censored", "40% Censored", "80% Censored")
    ),
    n = factor(n),
    L = factor(L),
    coverage_mcse = sqrt(coverage*(1-coverage)/n_rep)
  )

# table for coverage probability
coverage_summary <- summary_results %>%
  group_by(cens, n, L, sampler, parameter) %>%
  summarize(
    n_rep = n(),
    coverage = mean(covered),
    .groups = "drop"
  )

# data censoring summaries and IC widths
data_summary_files <- list.files(
  summary_dir, pattern = "^data_summaries_array_[0-9]+//.csv", 
  full.names = TRUE
)

data_summaries <- data_summary_files %>%
  map_dfr(read.csv) %>%
  distinct(sim_id, .keep_all = TRUE) %>%
  left_join(
    sim_grid,
    by = c("sim_id", "n")
  )

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
  )
################################# save datasets
write.csv(
  summary_results,
  file.path(analysis_dir, 'summary_results.csv'),
  row.names = FALSE
)

write.csv(
  master_summary,
  file.path(analysis_dir, 'master_summary.csv'),
  row.names = FALSE
)

write.csv(
  coverage_summary,
  file.path(analysis_dir, 'coverage_summary.csv'),
  row.names = FALSE
)

write.csv(
  censoring_summary,
  file.path(analysis_dir, "censoring_summary.csv"),
  row.names = FALSE
)
#################################


################################# create and save figures

# figure for bias
bias_plot <- ggplot(summary_results, aes(x = factor(L), y = (mean - truth), color = factor(n))) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid(parameter ~ cens, scales = "free_y") +
  labs(x = "Truncation Level (L)", y = "Bias", title = "Bias across censoring levels and sample size") +
  theme_bw()
ggsave(filename = file.path(analysis_dir, "bias_plot.png"), plot = bias_plot,
       width = 10, height = 6, dpi = 300)


# figure for mean
posterior_mean_plot <- ggplot(summary_results, aes(x = factor(L), y = mean, color = factor(n))) +
  geom_boxplot() +
  geom_hline(aes(yintercept = truth), linetype = 2, color = "red")+
  facet_grid(parameter ~ cens,
             scales = "free_y") +
  labs(x = "Truncation Level (L)", y = "Mean") +
  theme_bw()
ggsave(file.path(analysis_dir, "posterior_mean_plot.png"), plot = posterior_mean_plot,
       width = 10, height = 6, dpi = 300)


# figure for empirical coverage - faceted line plot
coverage_plot <- ggplot(coverage_summary, aes(x = factor(L), y = coverage, color = factor(n), group = n)) +
  geom_hline(yintercept = 0.95, linetype = 2) +
  geom_line() +
  geom_point() +
  facet_grid(parameter ~ cens, scales = "free_y") +
  scale_y_continuous(limits = c(0.8, 1)) +
  labs(
    x = "Truncation Level (L)",
    y = "Empirical Coverage",
    color = "Sample Size",
    title = "Coverage across simulation settings"
  ) +
  theme_bw()
ggsave(file.path(analysis_dir, "coverage_plot.png"), plot = coverage_plot,
       width = 10, height = 6, dpi = 300)
#################################

################################# save RDS of combined data
saveRDS(summary_results, file.path(analysis_dir, "summary_results.rds"))
saveRDS(master_summary, file.path(analysis_dir, "master_summary.rds"))
