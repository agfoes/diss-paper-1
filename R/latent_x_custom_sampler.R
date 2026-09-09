sampler_latent_x_normal <- nimbleFunction(
  contains = sampler_BASE,
  
  ## runs when nimble constructs the sampler -- contains things that don't change on each iteration
  setup = function(model, mvSaved, target, control) {
    
    # model parameters
    p <- control$p
    ncen <- control$ncen

    # recalculate nodes depending on changes in x_cen
    calcNodes <- model$getDependencies(target)
  },
  
  ## runs at every MCMC iteration
  run = function() {
    
    # get cluster parameter values for current iteration
    gamma_wi <- model$gamma[model$idx_cen, 1:p]
    sigmasq_wi <- model$sigmasq[model$idx_cen]
    
    # calculate vector/matrix products
    for (k in 1:ncen) {
      gammaZ[k] <- inprod(model$gamma[model$idx_cen[k], 1:p],
                          model$z_cen[k, 1:p])
    }
    betaZ <- model$z_cen[1:ncen, 1:p] %*% model$beta[2:(p+1)]
    
    A <- model$tau * model$beta[1]^2 + 1 / sigmasq_wi
    B <- gammaZ / sigmasq_wi + model$tau * model$beta[1] * (model$y_cen[1:ncen] - betaZ)
 
    mean <- B / A
    sd <- sqrt(1 / A)
    
    x_cen_new <- rnorm(mean = mean, 
                       sd = sd)
    
    model[[target]] <<- x_cen_new
    
    ## update dependencies
    model$calculate(calcNodes)
    copy(
      from = model, 
      to = mvSaved, 
      row = 1,
      nodes = calcNodes, 
      logProb = TRUE
    )
  },
  
  methods = list(reset = function() {})
)