sampler_beta_normal_conjugate_block <- nimbleFunction(
  contains = sampler_BASE,
  
  setup = function(model, mvSaved, target, control) {
    pzy <- control$pzy
    py <- pzy + 1
    
    nobs <- as.integer(control$nobs)
    ncen <- as.integer(control$ncen)
    n <- nobs + ncen
    
    hasObs <- nobs > 0
    hasCen <- ncen > 0
    
    if (hasObs) {
      xObsName <- "x_obs"
      yObsName <- "y_obs"
      zyObsName <- "zy_obs"
    } else {
      xObsName <- "x_cen"
      yObsName <- "y_cen"
      zyObsName <- "zy_cen"
    }
    
    if (hasCen) {
      xCenName <- "x_cen"
      yCenName <- "y_cen"
      zyCenName <- "zy_cen"
    } else {
      xCenName <- "x_obs"
      yCenName <- "y_obs"
      zyCenName <- "zy_obs"
    }
    
    s_beta_inv <- inverse(model$beta_cov[1:py, 1:py])
    beta_prior_num <- s_beta_inv %*% asCol(model$beta_mean[1:py])
    
    W_fixed <- matrix(0, nrow = n, ncol = py)
    Y <- numeric(n)
    
    if (hasObs) {
      for (k in 1:nobs) {
        
        W_fixed[k, 1] <- model[[xObsName]][k]
        
        if (pzy > 1) {
          for (j in 1:pzy) {
            W_fixed[k, j + 1] <- model[[zyObsName]][k, j]
          }
        } else if (pzy == 1) {
          W_fixed[k, 2] <- model[[zyObsName]][k]
        }
        
        Y[k] <- model[[yObsName]][k]
      }
    }
    
    if (hasCen) {
      for (k in 1:ncen) {
        
        i <- nobs + k
        
        ## leave W_fixed[i, 1] for x_cen
        
        if (pzy > 1) {
          for (j in 1:pzy) {
            W_fixed[i, j + 1] <- model[[zyCenName]][k, j]
          }
        } else if (pzy == 1) {
          W_fixed[i, 2] <- model[[zyCenName]][k]
        }
        
        Y[i] <- model[[yCenName]][k]
      }
    }
    
    calcNodes <- model$getDependencies(target)
  },
  
  run = function() {
    
    W <- W_fixed
    
    if (hasCen) {
      for (k in 1:ncen) {
        i <- nobs + k
        W[i, 1] <- model[[xCenName]][k]
      }
    }
    
    # V_beta = tau * W^T W + Sigma_beta^{-1}
    V_beta <- s_beta_inv + model$tau[1] * (t(W) %*% W)
    
    # W^T Y
    tWY <- (t(W)%*%Y)
    
    # tau * W^T Y + Sigma_beta^{-1} beta_0
    mu_beta_num <- model$tau[1] * tWY + beta_prior_num
    
    # V_beta^{-1}
    V_beta_inv <- inverse(V_beta)
    
    # mu_beta = V_beta^{-1} mu_beta_num
    mu_beta <- (V_beta_inv %*% mu_beta_num)
    
    ## joint gibbs draw
    V_beta_cholesky <- chol(V_beta)
    
    beta_new <- rmnorm_chol(1, 
                            mean = mu_beta[, 1],
                            cholesky = V_beta_cholesky, 
                            prec_param = TRUE)
    model[[target]] <<- beta_new
    
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