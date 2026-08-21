
# This script will generate and save a file containing simulation settings in a grid form to be 
# read by slurm arrays and other .R files. The number of arrays to be used needs
# to be reflected in this script as well as in the slurm file. The true number of clusters of 
# the censored covariate (L) needs to be specified in here.

library(dplyr)

dir <- "/work/users/a/g/agfoes/P1"
source(paste0(dir, "/R/config_sampler_comparison.R"))                                       # config

result_dir <- file.path(
  dir,
  "results",
  config[["run_name"]]
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

nrep = config[["nrep"]]
n = c(100, 1000)

## data frame of parameter values pre-determined for ideal censoring rates
cens_grid <- data.frame(
  cens = c('93', '80', '40', '20'),
  
  alpha_start_1 = c(1250, 1080, 750, 550),
  alpha_start_2 = c(500, 500, 350, 150),
  alpha_start_3 = c(390, 390, 390, 390),
  
  alpha_width_1 = c(log(400), log(900), log(800), log(1100)),
  alpha_width_2 = c(0.35, 0.50, 1.50, 2.00),
  alpha_width_3 = c(-0.2, -0.3, -0.2, -0.2),
  
  eps_obs = c(450, 400, 550, 1000)
) %>%
  filter(cens != "93")

## expand combinations of variables
sim_grid <- expand.grid(
  cens = cens_grid$cens,
  n = n,
  rep = seq_len(nrep),
  stringsAsFactors = FALSE,
  KEEP.OUT.ATTRS = FALSE
)

## combine sim and cens grids, add sim_id, order columns
sim_grid <- merge(sim_grid, cens_grid, by = "cens", sort = FALSE)  %>%
  mutate(
    array_id = rep(
      seq_len(config[["n_arrays"]]),
      length.out = n()
    ),
    sim_id = row_number(), .before = 1,
    save_level = config[["save_level"]],
    save_level = ifelse(
      sim_id %in% config[["trace_runs"]],
      3, save_level
      )
    ) %>%
  select(
    sim_id, array_id, cens, n, rep, save_level,
    alpha_start_1, alpha_start_2, alpha_start_3, 
    alpha_width_1, alpha_width_2, alpha_width_3, eps_obs
  )

## save sim_grid
write.csv(sim_grid,
          file = file.path(result_dir, "sim_grid.csv"),
          row.names = FALSE)
