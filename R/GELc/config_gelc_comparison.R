config <- list(
  project_dir = '/work/users/a/g/agfoes/P1',
  run_name = 'gelc_comparison_08292026',
  
  n_values = c(100, 300, 500),
  mu_values = c(3, 6, 9, 12),
  nrep = 1000,
  
  beta_0 = 2,
  beta_x = 0.1,
  beta_z = 0,
  tau = 1,
  
  pz = 2,
  L = 50,
  mu_gamma = 0,
  
  niter = 100000,
  nburnin = 40000,
  thin = 1,
  
  sampler = 'custom_joint',
  
  sim_seed = 20260829,
  
  key_vars = c(
    'beta[1]',
    'beta[2]',
    'beta[3]',
    'tau'
  ),
  
  truth = c(
    'beta[1]' = 0.1,
    'beta[2]' = 2,
    'beta[3]' = 0,
    tau = 1
  )
)