buildModelCode <- function(
    family,
    nobs,
    ncen
) {
  
  family_def <- as.character(family)
  
  base_block <- "
    beta[1:py] ~ dmnorm(
      mean = beta_mean[1:py],
      cov = beta_cov[1:py, 1:py]
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
      
      if (pzx > 1) {
        for (j in 1:pzx) {
          gammaTilde[h, j] ~ dnorm(mu_gamma, var = sigmasqTilde[h])
        }
      } else if (pzx == 1) {
        gammaTilde[h] ~ dnorm(mu_gamma, var = sigmasqTilde[h])
      }
      
      
    }
  "
  
  obs_block <- paste0("
    for (k in 1:nobs) {
      
      # cluster assignment
      xi[idx_obs[k]] ~ dcat(w[1:L])
      
      # assign component parameters based on cluster assignment
      if (pzx > 1) {
        for (j in 1:pzx) {
          gamma[idx_obs[k],j] <- gammaTilde[xi[idx_obs[k]],j]
        }
      } else if (pzx == 1) {
        gamma[idx_obs[k]] <- gammaTilde[xi[idx_obs[k]]]
      }
  
      sigmasq[idx_obs[k]] <- sigmasqTilde[xi[idx_obs[k]]]
      
      # latent X
      if (pzx > 1) {
        mean_x[idx_obs[k]] <- inprod(gamma[idx_obs[k], 1:pzx], zx_obs[k, 1:pzx])
      } else if (pzx == 1) {
        mean_x[idx_obs[k]] <- gamma[idx_obs[k]] * zx_obs[k]
      }
      
      x_obs[k] ~ dnorm(
        mean = mean_x[idx_obs[k]],
        var = sigmasq[idx_obs[k]])
      
      # outcome model
      if (pzy > 1) {
        eta[idx_obs[k]] <- beta[1]*x_obs[k] + inprod(zy_obs[k, 1:pzy], beta[2:py])
      } else if (pzy == 1) {
        eta[idx_obs[k]] <- beta[1]*x_obs[k] + zy_obs[k] * beta[2]
      } else if (pzy == 0) {
        eta[idx_obs[k]] <- beta[1]*x_obs[k]
      }
      
      
      y_obs[k] ~ DYNAMICGLM(
        eta = eta[idx_obs[k]],
        family = ", family_def, ",
        tau = tau
      )
    }
  "
    )
  
  cen_block <- paste0("
    for (k in 1:ncen) {
      
      # cluster assignment
      xi[idx_cen[k]] ~ dcat(w[1:L])
      
      # assign component parameters based on cluster assignment
      if (pzx > 1) {
        for (j in 1:pzx) {
          gamma[idx_cen[k],j] <- gammaTilde[xi[idx_cen[k]],j]
        }
      } else if (pzx == 1) {
        gamma[idx_cen[k]] <- gammaTilde[xi[idx_cen[k]]]
      }
      
      sigmasq[idx_cen[k]] <- sigmasqTilde[xi[idx_cen[k]]]
  
      # latent X
      if (pzx > 1) {
        mean_x[idx_cen[k]] <- inprod(gamma[idx_cen[k], 1:pzx], zx_cen[k, 1:pzx])
      } else if (pzx == 1) {
        mean_x[idx_cen[k]] <- gamma[idx_cen[k]] * zx_cen[k]
      }
      x_cen[k] ~ dnorm(
        mean = mean_x[idx_cen[k]],
        var = sigmasq[idx_cen[k]])
      
      # outcome model
      if (pzy > 1) {
        eta[idx_cen[k]] <- beta[1]*x_cen[k] + inprod(zy_cen[k, 1:pzy], beta[2:py])
      } else if (pzy == 1) {
        eta[idx_cen[k]] <- beta[1]*x_cen[k] + zy_cen[k] * beta[2]
      } else if (pzy == 0) {
        eta[idx_cen[k]] <- beta[1]*x_cen[k]
      }
      
      
      y_cen[k] ~ DYNAMICGLM(
        eta = eta[idx_cen[k]],
        family = ", family_def, ",
        tau = tau
      )
    }
  "
    )
  
  if (nobs > 0 & ncen > 0) {
    model_text <- paste("{", base_block, cen_block, obs_block, "}")
  } else if (nobs > 0) {
    model_text <- paste("{", base_block, obs_block, "}")
  } else if (ncen > 0) {
    model_text <- paste("{", base_block, cen_block, "}")
  }
  
  model_code <- eval(call("nimbleCode", parse(text = model_text)[[1]]))
  
  return(model_code)
  
}