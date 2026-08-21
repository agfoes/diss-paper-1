#' make_conf
#' 
#' @description
#' Helper function to generate nimble conf object depending on sampler type
#' 
#' @param model 
#' @param sampler_type 
#'
#' @returns
#' @export
#'
#' @examples
make_conf <- function(model, sampler_type = 'default', pz, L, vars, nobs, ncen) {
  
  conf <- configureMCMC(model,
                        useConjugacy = TRUE)
  
  # add monitors
  conf$addMonitors(vars)
  conf$addMonitors("xi")
  
  # regression parameters and variance
  beta_targets <- paste0("beta[", seq_len(pz+1), "]")
  gamma_targets <- as.vector(outer(seq_len(L), seq_len(pz),
                                   function(h, j) paste0('gammaTilde[', h, ', ', j, ']')))
  sigma_targets <- paste0('sigmasqTilde[', seq_len(L), ']')
  
  # block sampler
  if (sampler_type == 'block') {
    
    # remove and add block for betas
    conf$removeSamplers("beta")
    conf$addSampler(
      target = "beta",
      type = 'RW_block',
      control = list(scale = 0.05, adaptInterval = 100)
    )
  }
  
  # slice sampler
  if (sampler_type == 'slice') {
    
    # remove and add slice for all betas
    conf$removeSampler("beta")
    for (target in beta_targets) {
      conf$addSampler(target = target, type = "slice")
    }
  }
  
  # AF slice sampler (block)
  if (sampler_type == 'AF_slice') {
    
    # remove and add block for betas
    conf$removeSamplers("beta")
    conf$addSampler(
      target = "beta",
      type = 'AF_slice'
    )
  }
  
  ## custom joint conjugate sampler
  if (sampler_type == "custom_joint") {
    
    conf$removeSamplers("beta")
    conf$addSampler(
      target = "beta",
      type = sampler_beta_conjugate_block,
      control = list(
        p = pz,
        nobs = as.integer(nobs),
        ncen = as.integer(ncen)
      )
    )
  }
  
  conf$printMonitors()
  conf$printSamplers("beta")
  return(conf)
}

#' make_conf
#' 
#' @description
#' Helper function to generate nimble conf object depending on sampler type
#' 
#' @param model 
#' @param sampler_type 
#'
#' @returns
#' @export
#'
#' @examples
make_conf_univariate <- function(model, sampler_type = 'default', pz, L, vars) {
  
  conf <- configureMCMC(model,
                        useConjugacy = TRUE)
  
  # add monitors
  conf$addMonitors(vars)
  conf$addMonitors("xi")
  
  # outcome regression parameters
  beta_targets <- c('betax', paste0('betaz[', seq_len(pz), ']'))
  
  # regression parameters and variance
  gamma_targets <- as.vector(outer(seq_len(L), seq_len(pz),
                                   function(h, j) paste0('gammaTilde[', h, ', ', j, ']')))
  sigma_targets <- paste0('sigmasqTilde[', seq_len(L), ']')
  
  # block sampler
  if (sampler_type == 'block') {
    
    # remove and add block for betas
    conf$removeSamplers(beta_targets)
    conf$addSampler(
      target = beta_targets,
      type = 'RW_block',
      control = list(scale = 0.05, adaptInterval = 100)
    )
  }
  
  # slice sampler
  if (sampler_type == 'slice') {
    
    # remove and add slice for all betas
    conf$removeSamplers(beta_targets)
    for (t in beta_targets) {
      conf$addSampler(t, type = 'slice')
    }
  }
  
  # AF slice sampler (block)
  if (sampler_type == 'AF_slice') {
    
    # remove and add block for betas
    conf$removeSamplers(beta_targets)
    conf$addSampler(
      target = beta_targets,
      type = 'AF_slice'
    )
  }
  
  conf$printMonitors()
  conf$printSamplers(beta_targets)
  return(conf)
}

#' run_sampler
#' 
#' @description
#' Helper function to call function to build and configure model based on
#' sampler type and run/obtain samples after model is built. 
#' 
#'
#' @param model 
#' @param sampler_type 
#' @param pz 
#' @param L 
#' @param niter 
#' @param nburnin 
#' @param thin 
#'
#' @returns
#' @export
#'
#' @examples
run_sampler_univariate <- function(model, 
                        sampler_type, 
                        pz, 
                        L, 
                        niter = 240000, 
                        nburnin = 72000, 
                        thin = 1,
                        vars_to_monitor) {
  conf <- make_conf(model, 
                    sampler_type, 
                    pz = pz, L = L, 
                    vars = vars_to_monitor)
  
  mcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(model)
  Cmcmc <- compileNimble(mcmc, project = Cmodel)
  
  runtime <- system.time({
    samples <- runMCMC(
      Cmcmc,
      niter = niter,
      nburnin = nburnin,
      thin = thin,
      nchains = 1,
      samplesAsCodaMCMC = TRUE
    )
  })

  list(
    samples = samples,
    runtime = runtime['elapsed']
  )
}

#' run_sampler
#' 
#' @description
#' Helper function to call function to build and configure model based on
#' sampler type and run/obtain samples after model is built. 
#' 
#'
#' @param model 
#' @param sampler_type 
#' @param pz 
#' @param L 
#' @param niter 
#' @param nburnin 
#' @param thin 
#'
#' @returns
#' @export
#'
#' @examples
run_sampler <- function(model, 
                        sampler_type, 
                        pz, 
                        L, 
                        niter = 240000, 
                        nburnin = 72000, 
                        thin = 1,
                        vars_to_monitor,
                        nobs,
                        ncen) {
  conf <- make_conf(model, 
                    sampler_type, 
                    pz = pz, L = L, 
                    vars = vars_to_monitor,
                    nobs, ncen)
  
  mcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(model)
  Cmcmc <- compileNimble(mcmc, project = Cmodel)
  
  runtime <- system.time({
    samples <- runMCMC(
      Cmcmc,
      niter = niter,
      nburnin = nburnin,
      thin = thin,
      nchains = 1,
      samplesAsCodaMCMC = TRUE
    )
  })
  
  list(
    samples = samples,
    runtime = runtime['elapsed']
  )
}

