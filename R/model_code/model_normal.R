model_code <- nimbleCode({
  
  ### hyper-priors
  betax ~ dnorm(0, sd = sigma_bx)
  
  for (j in 1:p) {
    betaz[j] ~ dnorm(0, sd = sigma_bz)
  }
  
  tau ~ dgamma(1,1)
  #alpha ~ dgamma(1,1)
  lalpha ~ dLogGamma(1, 1)
  alpha <- exp(lalpha)
  
  ### stick-breaking
  for (h in 1:(L-1)) {
    v[h] ~ dbeta(1, alpha)
  }
  w[1:L] <- stick_breaking(v[1:(L-1)])
  
  ### component-level parameters
  for (h in 1:L) {
    sigmasqTilde[h] ~ dinvgamma(1,1)
    
    for (j in 1:p) {
      gammaTilde[h, j] ~ dnorm(mu_gamma, var = sigmasqTilde[h])
    }
  }
  
  ### subject-level model
  
  for (k in 1:nobs) {
    
    # cluster assignment
    xi[idx_obs[k]] ~ dcat(w[1:L])
    
    # assign component parameters based on cluster assignment
    for (j in 1:p) {
      gamma[idx_obs[k],j] <- gammaTilde[xi[idx_obs[k]],j]
    }
    sigmasq[idx_obs[k]] <- sigmasqTilde[xi[idx_obs[k]]]
    
    # latent X
    x_obs[k] ~ dnorm(
      mean = inprod(gamma[idx_obs[k], 1:p], z_obs[k, 1:p]),
      var = sigmasq[idx_obs[k]])
    
    # outcome model
    eta[idx_obs[k]] <- inprod(z_obs[k, 1:p], betaz[1:p]) + betax*x_obs[k]
    
    y_obs[k] ~ dnorm(
      mean = eta[idx_obs[k]],
      tau = tau
    )
  }
  
  
  for (k in 1:ncen) {
    
    # cluster assignment
    xi[idx_cen[k]] ~ dcat(w[1:L])
    
    # assign component parameters based on cluster assignment
    for (j in 1:p) {
      gamma[idx_cen[k],j] <- gammaTilde[xi[idx_cen[k]],j]
    }
    sigmasq[idx_cen[k]] <- sigmasqTilde[xi[idx_cen[k]]]
    
    # latent X
    x_cen[k] ~ dnorm(
      mean = inprod(gamma[idx_cen[k], 1:p], z_cen[k, 1:p]),
      var = sigmasq[idx_cen[k]])
    
    # censoring constraint
    constraint_data[k] ~ dconstraint(
      (x_cen[k] > CL[k] & x_cen[k] <= CR[k])
    )
    
    # outcome model
    eta[idx_cen[k]] <- inprod(z_cen[k, 1:p], betaz[1:p]) + betax*x_cen[k]
    
    y_cen[k] ~ dnorm(
      mean = eta[idx_cen[k]],
      tau = tau
    )
  }
})