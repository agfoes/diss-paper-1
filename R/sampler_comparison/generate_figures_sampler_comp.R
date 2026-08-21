# The purpose of this script is to create and save standard figures from
# analytic datasets produced by compile_results.R.

library(tidyverse)

################################ user input required #####

dir <- "/work/users/a/g/agfoes/P1"
sim_name <- "sampler_comp_07292026"
config_file_name <- "config_sampler_comparison"

################################ user input required #####

source(file.path(dir, "results", sim_name, paste0(config_file_name, ".R")))

## define directories
result_dir <- file.path(dir, "results", config[["run_name"]])
analysis_dir <- file.path(result_dir, "analysis")
figures_dir <- file.path(analysis_dir, "figures")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

## read analytic datasets
posterior_results <- readRDS(file.path(analysis_dir, "posterior_results.rds"))
master_summary <- readRDS(file.path(analysis_dir, "master_summary.rds"))
data_summaries <- readRDS(file.path(analysis_dir, "data_summaries.rds"))
censoring_summary <- readRDS(file.path(analysis_dir, "censoring_summary.rds"))

## plot theme common to all plots
plot_theme <- theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank()
  )

## save figures as .pdf and .png
save_plot <- function(plot, name, width = 11, height = 8) {
  ggsave(file.path(figures_dir, paste0(name, ".png")),
         plot, width = width, height = height, dpi = 300)
  
  ggsave(file.path(figures_dir, paste0(name, ".pdf")),
         plot, width = width, height = height)
}

################################################################################ coverage plot

coverage_plot <- master_summary %>%
  ggplot(aes(sampler, coverage, color = n, group = n)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_errorbar(aes(ymin = pmax(coverage - 1.96*coverage_mcse, 0),
                    ymax = pmin(coverage + 1.96*coverage_mcse, 1)),
                width = 0.1) +
  geom_point() +
  geom_line() +
  facet_grid(parameter ~ cens, scales = "free_y") +
  #coord_cartesian(ylim = c(0,1)) +
  labs(
    x = "Sampler",
    y = "Coverage probability",
    color = "Sample size"
  ) +
  plot_theme

save_plot(coverage_plot, "coverage_plot")

################################################################################ bias plot

bias_plot <- master_summary %>%
  ggplot(aes(sampler, bias, color = n, group = n)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  geom_line() +
  facet_grid(parameter ~ cens, scales = "free_y") +
  labs(
    x = "Sampler",
    y = "Bias",
    color = "Sample size"
  ) +
  plot_theme

save_plot(bias_plot, "bias_plot")


################################################################################ bias plot with free y-axes (unstacked by sample size)

for (sample_size in levels(master_summary$n)) {
  
  bias_free_plot <- master_summary %>%
    filter(n == sample_size) %>%
    ggplot(aes(sampler, bias, group = 1)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point() +
    geom_line() +
    facet_wrap(
      vars(parameter, cens),
      scales = "free_y",
      ncol = 3,
      labeller = label_both
    ) +
    labs(
      title = paste0("Bias for n = ", sample_size),
      x = "Sampler",
      y = "Bias"
    ) +
    plot_theme
  
  save_plot(
    bias_free_plot,
    paste0("bias_free_y_n", sample_size)
  )
}

################################################################################ RMSE plot

rmse_plot <- master_summary %>%
  ggplot(
    aes(
      x = sampler,
      y = rmse,
      group = n,
      color = n
    )
  ) +
  geom_line() +
  geom_point() +
  facet_grid(
    rows = vars(parameter),
    cols = vars(cens),
    scales = 'free_y'
  ) +
  labs(
    x = 'Truncation level',
    y = 'RMSE',
    color = 'Sample size'
  ) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    strip.background = element_rect(fill = 'white'),
    panel.grid.minor = element_blank()
  )

save_plot(rmse_plot, "rmse_plot")

################################################################################ ESS
ess_plot <- master_summary %>%
  group_by(sampler, n, cens) %>%
  summarize(min_sampler_ess = min(min_ess),
            .groups = "drop") %>%
  ggplot(
    aes(
      x = cens,
      y = min_sampler_ess,
      fill = sampler
    )
  ) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    x = 'Censoring rate',
    y = 'Minimum ESS',
    fill = 'Sampler',
    title = 'MCMC Minimum Effective Sample Size'
  ) +
  facet_wrap(~n) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

save_plot(ess_plot, "ess_plot")
