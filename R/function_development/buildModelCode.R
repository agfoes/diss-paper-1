buildModelCode <- function(
    family,
    nobs,
    ncen
) {
  
  family_def <- as.character(family)
  
  # dimension changes for one-subject group sizes
  if (nobs == 1) {
    obs_idx <- 'idx_obs'
    obs_k <- '1'
    obs_loop_start <- ''
    obs_loop_end <- ''
  } else {
    obs_idx <- 'idx_obs[k]'
    obs_k <- 'k'
    obs_loop_start <- paste0('for (k in 1:', nobs, ') {')
    obs_loop_end <- '}'
  }
  if (ncen == 1) {
    cen_idx <- 'idx_cen'
    cen_k <- '1'
    cen_loop_start <- ''
    cen_loop_end <- ''
  } else {
    cen_idx <- 'idx_cen[k]'
    cen_k <- 'k'
    cen_loop_start <- paste0('for (k in 1:', ncen, ') {')
    cen_loop_end <- '}'
  }
  
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
  
  obs_block <- paste0(
      "
  ", obs_loop_start, "
  
    # cluster assignment
    xi[", obs_idx, "] ~ dcat(w[1:L])
  
    # assign component parameters based on cluster assignment
    if (pzx > 1) {
      for (j in 1:pzx) {
        gamma[", obs_idx, ", j] <- gammaTilde[xi[", obs_idx, "], j]
      }
    } else if (pzx == 1) {
      gamma[", obs_idx, "] <- gammaTilde[xi[", obs_idx, "]]
    }
  
    sigmasq[", obs_idx, "] <- sigmasqTilde[xi[", obs_idx, "]]
  
    # latent X
    if (pzx > 1) {
      mean_x[", obs_idx, "] <- inprod(
        gamma[", obs_idx, ", 1:pzx],
        zx_obs[", obs_k, ", 1:pzx]
      )
    } else if (pzx == 1) {
      mean_x[", obs_idx, "] <- gamma[", obs_idx, "] * zx_obs[", obs_k, "]
    }
  
    x_obs[", obs_k, "] ~ dnorm(
      mean = mean_x[", obs_idx, "],
      var = sigmasq[", obs_idx, "]
    )
  
    # outcome model
    if (pzy > 1) {
      eta[", obs_idx, "] <- beta[1] * x_obs[", obs_k, "] +
        inprod(zy_obs[", obs_k, ", 1:pzy], beta[2:py])
    } else if (pzy == 1) {
      eta[", obs_idx, "] <- beta[1] * x_obs[", obs_k, "] +
        zy_obs[", obs_k, "] * beta[2]
    } else if (pzy == 0) {
      eta[", obs_idx, "] <- beta[1] * x_obs[", obs_k, "]
    }
  
    y_obs[", obs_k, "] ~ DYNAMICGLM(
      eta = eta[", obs_idx, "],
      family = ", family_def, ",
      tau = tau
    )
  
  ", obs_loop_end, "
  "
  )
  
  cen_block <- paste0(
      "
  ", cen_loop_start, "
  
    # cluster assignment
    xi[", cen_idx, "] ~ dcat(w[1:L])
    
    # assign component parameters based on cluster assignment
    if (pzx > 1) {
      for (j in 1:pzx) {
        gamma[", cen_idx, ", j] <- gammaTilde[xi[", cen_idx, "], j]
      }
    } else if (pzx == 1) {
      gamma[", cen_idx, "] <- gammaTilde[xi[", cen_idx, "]]
    }
    
    sigmasq[", cen_idx, "] <- sigmasqTilde[xi[", cen_idx, "]]
  
    # latent X
    if (pzx > 1) {
      mean_x[", cen_idx, "] <- inprod(
        gamma[", cen_idx, ", 1:pzx],
        zx_cen[", cen_k, ", 1:pzx]
      )
    } else if (pzx == 1) {
      mean_x[", cen_idx, "] <- gamma[", cen_idx, "] * zx_cen[", cen_k, "]
    }
  
    x_cen[", cen_k, "] ~ dnorm(
      mean = mean_x[", cen_idx, "],
      var = sigmasq[", cen_idx, "]
    )
    
    # censoring constraint
    constraint_data[", cen_k, "] ~ dconstraint(
      (x_cen[", cen_k, "] > CL_cen[", cen_k, "] &
      x_cen[", cen_k, "] <= CR_cen[", cen_k, "])
    )
    
    # outcome model
    if (pzy > 1) {
      eta[", cen_idx, "] <- beta[1] * x_cen[", cen_k, "] +
        inprod(zy_cen[", cen_k, ", 1:pzy], beta[2:py])
    } else if (pzy == 1) {
      eta[", cen_idx, "] <- beta[1] * x_cen[", cen_k, "] +
        zy_cen[", cen_k, "] * beta[2]
    } else if (pzy == 0) {
      eta[", cen_idx, "] <- beta[1] * x_cen[", cen_k, "]
    }
    
    y_cen[", cen_k, "] ~ DYNAMICGLM(
      eta = eta[", cen_idx, "],
      family = ", family_def, ",
      tau = tau
    )
  
  ", cen_loop_end, "
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