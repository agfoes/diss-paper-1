
# This file answers the following questions:
# 1. Where are these results being saved?
# 2. What simulation setup file am I reading?
# 3. What posterior samples (if any) am I saving?

# This will be the only R script to change when simulations change (other than the script generating simulation setting grid).
# This script will call all other static R functions and source files found in P1/R/

config <- list(
  #project_dir = '/work/users/a/g/agfoes/P1',
  project_dir = "~/Library/CloudStorage/OneDrive-UniversityofNorthCarolinaatChapelHill/Dissertation/Diss-Paper1",
  run_name = "custom_sampler_test",
  
  pz = 3, # dimension of fully observed covariate vector (matching datagen specification)
  sigma_bx = 1, # prior sd
  sigma_bz = 1, # prior sd
  mu_gamma = 0, # prior mean
  
  L = 50,
  samplers = c("AF_slice", "custom_joint", "slice", "block"),
  
  niter = 100000,
  nburnin = 40000,
  thin = 1,
  
  nrep = 100,
  sim_seed = 20260729,
  n_arrays = 100,
  
  # save levels:
  ## 0: posterior summaries and generated data
  ## 1:   + key posterior samples (posterior samples for "key_vars" specified below with thinning)
  ## 2:    + full posterior samples (with thinning)
  ## 3:     + save ALL posterior samples for select runs (trace_runs) to produce trace plots
  
  save_level = 1,
  save_full_thinning = 10000,
  save_key_thinning = 1000,
  trace_runs = c(1:10), # these specific sim_ids will have save_level 4 assigned to them in sim_grid.R
  
  # this variable controls what parameters are summarized and reported on
  key_vars = c(
    "beta",
    "tau"
  ),
  
  # this variable controls what parameters are monitored and saved
  all_vars = c(
    "beta",
    "tau",
    "lalpha",
    "x_cen",
    "xi",
    "gammaTilde",
    "sigmasqTilde"
  ),
  
  truth = c(0.7, 0, 1, -1, 1)
)
