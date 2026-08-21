# The purpose of this script is to create and save standard figures from
# analytic datasets produced by compile_results.R.

library(tidyverse)

################################ user input required #####

dir <- "/work/users/a/g/agfoes/P1"
sim_name <- "expanded_rep_sims_08022026"
config_file_name <- "config_08022026"

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
  ggplot(aes(L, coverage, color = n, group = n)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_errorbar(aes(ymin = pmax(coverage - 1.96*coverage_mcse, 0),
                    ymax = pmin(coverage + 1.96*coverage_mcse, 1)),
                width = 0.1) +
  geom_point() +
  geom_line() +
  facet_grid(parameter ~ cens, scales = "free_y") +
  #coord_cartesian(ylim = c(0,1)) +
  labs(
    x = "Truncation level",
    y = "Coverage probability",
    color = "Sample size"
  ) +
  plot_theme

save_plot(coverage_plot, "coverage_plot")

################################################################################ bias plot

bias_plot <- master_summary %>%
  ggplot(aes(L, bias, color = n, group = n)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point() +
  geom_line() +
  facet_grid(parameter ~ cens, scales = "free_y") +
  labs(
    x = "Truncation level",
    y = "Bias",
    color = "Sample size"
  ) +
  plot_theme

save_plot(bias_plot, "bias_plot")


################################################################################ bias plot with free y-axes (unstacked by sample size)

for (sample_size in levels(master_summary$n)) {
  
  bias_free_plot <- master_summary %>%
    filter(n == sample_size) %>%
    ggplot(aes(L, bias, group = 1)) +
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
      x = "Truncation level",
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
      x = L,
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

################################################################################ SD ratio heatmap
sd_ratio_heatmap <- master_summary %>%
  filter(L == '25') %>%
  ggplot(
    aes(
      x = cens,
      y = parameter,
      fill = sd_ratio
    )
  ) +
  geom_tile(
    color = 'white',
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = sprintf('%.2f', sd_ratio)),
    size = 3.5
  ) +
  facet_wrap(
    ~ n,
    nrow = 1,
    labeller = labeller(
      n = function(x) paste0('n = ', x)
    )
  ) +
  scale_fill_gradient2(
    low = '#3B6FB6',
    mid = 'white',
    high = '#B64B4B',
    midpoint = 1,
    limits = c(0.8, 1.4),
    oob = scales::squish,
    name = 'Posterior SD /\nEmpirical SD'
  ) +
  labs(
    x = 'Censoring scenario',
    y = NULL,
    title = 'Calibration of Posterior Uncertainty',
    subtitle = 'DPMM truncation level L = 25'
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = 'bold'),
    legend.position = 'right',
    plot.title = element_text(face = 'bold')
  )

save_plot(sd_ratio_heatmap, "sd_ratio_heatmap")

################################################################################ ESS
ess_plot <- master_summary %>%
  filter(L == '25') %>%
  ggplot(
    aes(
      x = cens,
      y = median_ess,
      group = n,
      color = n
    )
  ) +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = q1_ess,
      ymax = q3_ess
    ),
    width = 0.10,
    linewidth = 0.5
  ) +
  geom_point(size = 2.5) +
  facet_wrap(
    ~ parameter,
    scales = 'free_y',
    ncol = 3,
    axes = "all"
  ) +
  labs(
    x = 'Censoring rate',
    y = 'Median effective sample size',
    color = 'Sample size',
    title = 'MCMC Effective Sample Size',
    subtitle = 'DPMM truncation level L = 25'
  ) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

save_plot(ess_plot, "ess_plot")

################################################################################ ESS
ess_sec_plot <- master_summary %>%
  filter(L == '25') %>%
  ggplot(
    aes(
      x = cens,
      y = median_ess_sec,
      group = n,
      color = n
    )
  ) +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = q1_ess_sec,
      ymax = q3_ess_sec
    ),
    width = 0.10,
    linewidth = 0.5
  ) +
  geom_point(size = 2.5) +
  facet_wrap(
    ~ parameter,
    scales = 'free_y',
    ncol = 3,
    axes = "all"
  ) +
  labs(
    x = 'Censoring rate',
    y = 'Median ESS per second',
    color = 'Sample size',
    title = 'MCMC Computational Efficiency',
    subtitle = 'DPMM truncation level L = 25'
  ) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )

save_plot(ess_sec_plot, "ess_sec_plot")

################################################################################ Trace plots (if samples were saved)

trace_file <- file.path(analysis_dir, "trace_samples.rds")

if (file.exists(trace_file)) {
  trace_samples <- readRDS(trace_file)
  trace_plot_dir <- file.path(figures_dir, "trace_plots")
  dir.create(trace_plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  trace_fits <- trace_samples %>%
    distinct(
      sim_id, cens, n, L, sampler
    ) %>%
    arrange(sim_id, L)
  
  for (i in seq_len(nrow(trace_fits))) {
    fit_info <- trace_fits[i, ]
    
    trace_data <- trace_samples %>%
      filter(
        sim_id == fit_info[["sim_id"]],
        L == fit_info[["L"]],
        parameter %in% config[["key_vars"]]
      )
    
    trace_plot <- ggplot(trace_data, aes(x = iteration, y = value)) +
      geom_line(linewidth = 0.18, alpha = 0.75) +
      facet_wrap(
        ~ parameter,
        scales = "free_y",
        ncol = 1,
        labeller = label_parsed
      ) +
      labs(
        x = "MCMC iteration",
        y = "Posterior draw",
        title = "Posterior Trace Plots",
        subtitle = paste0(
          "Simulation ", fit_info[["sim_id"]],
          "; n = ", fit_info[["n"]],
          "; ", fit_info[["cens"]], "% censoring",
          "; L = ", fit_info[["L"]]
        )
      ) +
      theme_bw()
    
    save_plot(trace_plot, 
              sprintf("trace_plots/trace_sim%05d_L_%03d", 
                      fit_info[["sim_id"]], fit_info[["L"]]
                      )
              )
    
  }
  
}
