
model_code <- nimbleCode({
  beta[1:(p+1)] ~ dmnorm(
    mean = beta_mean[1:(p+1)],
    cov = beta_cov[1:(p+1), 1:(p+1)]
  )
  
  tau ~ dgamma(1,1)
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
    gammaTilde[h] ~ dnorm(mu_gamma, var = sigmasqTilde[h])
  }
  
  for (k in 1:ncen) {
    
    # cluster assignment
    xi[idx_cen[k]] ~ dcat(w[1:L])
    
    # assign component parameters based on cluster assignment
    gamma[idx_cen[k]] <- gammaTilde[xi[idx_cen[k]]]
    sigmasq[idx_cen[k]] <- sigmasqTilde[xi[idx_cen[k]]]
    
    # latent X
    x_cen[k] ~ dnorm(
      mean = gamma[idx_cen[k]] * z_cen[k],
      var = sigmasq[idx_cen[k]])
    
    # outcome model
    eta[idx_cen[k]] <- z_cen[k] * beta[2] + beta[1] * x_cen[k]
    
    
    y_cen[k] ~ dnorm(
      mean = eta[idx_cen[k]],
      tau = tau
    )
  }
})
  