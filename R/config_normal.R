
# This file answers the following questions:
# 1. Where are these results being saved?
# 2. What simulation setup file am I reading?
# 3. What posterior samples (if any) am I saving?

# This will be the only R script to change when simulations change (other than the script generating simulation setting grid).
# This script will call all other static R functions and source files found in P1/R/

config <- list(
  project_dir = '/work/users/a/g/agfoes/P1',
  run_name = 'expanded_rep_sims_08022026',
  
  pz = 3, # dimension of fully observed covariate vector (matching datagen specification)
  sigma_bx = 1, # prior sd
  sigma_bz = 1, # prior sd
  mu_gamma = 0, # prior mean
  
  sampler = "AF_slice",
  L_values = c(3, 25),

  niter = 100000,
  nburnin = 40000,
  thin = 1,
  
  nrep = 500,
  sim_seed = 20260727,
  n_arrays = 500,
  
  # save levels:
  ## 0: posterior summaries and generated data
  ## 1:   + key posterior samples (posterior samples for "key_vars" specified below with thinning)
  ## 2:    + full posterior samples (with thinning)
  ## 3:     + save ALL posterior samples for select runs (trace_runs) to produce trace plots
  
  save_level = 1,
  save_full_thinning = 10000,
  save_key_thinning = 1000,
  trace_runs = c(1:5), # these specific sim_ids will have save_level 4 assigned to them in sim_grid.R
  
  # this variable controls what parameters are summarized and reported on
  key_vars = c(
    "betax", 
    "betaz[1]",
    "betaz[2]",
    "betaz[3]",
    "tau"
  ),
  
  # this variable controls what parameters are monitored and saved
  all_vars = c(
    "betax", 
    "betaz[1]",
    "betaz[2]",
    "betaz[3]",
    "tau",
    "lalpha",
    "x_cen",
    "xi",
    "gammaTilde",
    "sigmasqTilde"
  ),
  
  truth = c(
    betax = 0.7,
    'betaz[1]' = 0,
    'betaz[2]' = 1,
    'betaz[3]' = -1,
    tau = 1
  )
)
