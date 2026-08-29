library(dplyr)

dir <- '/work/users/a/g/agfoes/P1'

source(
  file.path(
    dir,
    'R',
    'GELc',
    'config_gelc_comparison.R'
  )
)

result_dir <- file.path(
  config[['project_dir']],
  'results',
  config[['run_name']]
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

sim_grid <- data.frame(
  rep = seq_len(config[['nrep']]),
  array_id = seq_len(config[['nrep']])
)

write.csv(
  sim_grid,
  file = file.path(
    result_dir,
    'sim_grid.csv'
  ),
  row.names = FALSE
)