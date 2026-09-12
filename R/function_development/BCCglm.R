BCCglm <- function(
    data,
    outcome_formula,
    covariate_formula,
    censored_covariate,
    censoring_bounds,
    family <- "Gamma(link = 'log')"
    ) {
  
  # remove censored variable from outcome formula
  observed_formula <- update(
    outcome_formula,
    paste(". ~ . -", paste(censored_covariate, collapse = " - "))
  )
  
  # create model frame, response, matrices
  outcome_mf <- model.frame(observed_formula, data = data)
  Y <- model.response(outcome_mf)
  Zy <- model.matrix(observed_formula, data = data) # automatically adds an intercept
  Zx <- model.matrix(delete.response(terms(covariate_formula)),
                     data = data) # automatically adds an intercept
  
  # calculate dimensions of design matrices
  pzy <- dim(Zy)[2] # includes intercept + dim of supplied uncensored covariates
  py <- pzy + 1 # add one for the censored covariate
  pzx <- dim(Zx)[2]
  
  # returns list of nobs, ncen, Ndata, Nconstants, Ninits
  nimble_inputs <- nimbleInputs(CL = data$CL,
                                CR = data$CR,
                                Zy = Zy,
                                Zx = Zx,
                                Y = Y)
  
  # outcome distribution family and link function for model code
  family_def <- as.character(family)
  model_code <- buildModelCode(family = as.character(family_def),
                               nobs = nimble_inputs$nobs,
                               ncen = nimble_inputs$ncen)
  
  nimbleOptions(enableMacros = TRUE)
  model <- nimbleModel(
    code = model_code,
    data = nimble_inputs$Ndata,
    constants = nimble_inputs$Nconstants,
    inits = nimble_inputs$Ninits
  )
  
  conf <- configureMCMC(model, useConjugacy = TRUE)
  conf$removeSamplers("beta")
  
  sampler_type = "AF_slice"
  if (family_def == "gaussian(link = 'identity')") {
    conf$addSampler(
      target = "beta",
      type = sampler_beta_normal_conjugate_block,
      control = list(
        pzy = pzy,
        ncen = as.integer(nimble_inputs$ncen),
        nobs = as.integer(nimble_inputs$nobs)
      )
    )
  } else if (sampler_type == 'block') {
    
    conf$addSampler(
      target = "beta",
      type = 'RW_block',
      control = list(scale = 0.05, 
                     adaptInterval = 100)
    )
  } else if (sampler_type == 'AF_slice') {
    
    conf$addSampler(
      target = "beta",
      type = 'AF_slice'
    )
  }
  
  mcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(model)
  Cmcmc <- compileNimble(mcmc, project = Cmodel)
  
  runtime <- system.time({
    samples <- runMCMC(
      Cmcmc,
      niter = 100,
      nburnin = 40,
      thin = 1,
      nchains = 1,
      samplesAsCodaMCMC = TRUE
    )
  })
  
  return(list(
    runtime = runtime,
    samples = samples
  ))
}