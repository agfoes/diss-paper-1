
## base sampler ####
sampler_beta_conjugate_block <- nimbleFunction(
  contains = sampler_BASE,
  
  ## runs when nimble constructs the sampler -- contains things that don't change on each iteration
  setup = function(model, mvSaved, target, control) {
    
    # model dimensions
    p <- control$p
    p_beta <- p + 1
    
    nobs <- as.integer(control$nobs)
    ncen <- as.integer(control$ncen)
    n <- nobs + ncen
    
    # prior distribution precision: Sigma_{beta}^{-1}
    s_beta_inv <- inverse(model$beta_cov[1:p_beta, 1:p_beta])
    
    # fixed contribution of prior hypterparameters on posterior mean: Sigma_{beta}^{-1} * beta_0
    beta_prior_num <- (s_beta_inv %*% asCol(model$beta_mean[1:p_beta]))
    
    # recalculate nodes depending on changes in beta
    calcNodes <- model$getDependencies(target)
  },
  
  ## runs at every MCMC iteration
  run = function() {
    
    ## build current iteration's design matrix (X,Z) and outcome Y
    W <- matrix(0, nrow = n, ncol = p_beta)
    Y <- numeric(n)
    
    ## observed subjects
    if (nobs > 0) {
      for (k in 1:nobs) {
        W[k, 1] <- model$x_obs[k]
        
        if (p > 1) {
          for (j in 1:p) {
          W[k, j+1] <- model$z_obs[k, j]
          }
        } else if(p == 1){
          W[k, p+1] <- model$z_obs[k, p]
        }
        
        Y[k] <- model$y_obs[k]
      }
    }
    
    ## censored subjects
    if (ncen > 0) {
      for (k in 1:ncen) {
        i <- nobs + ks
        
        W[i, 1] <- model$x_cen[k]
        
        if (p > 0) {
          for (j in 1:p) {
            W[i, j+1] <- model$z_cen[k, j]
          }
        } else if(p == 1) {
          W[i, p+1] <- model$z_cen[k, p]
        }
        
        
        Y[i] <- model$y_cen[k]
      }
    }
    
    ## full conditional for beta
    
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

## custom sampler with mixed observed and censored subjects ####
sampler_beta_conjugate_block_mixed <- nimbleFunction(
  contains = sampler_BASE,
  
  ## runs when nimble constructs the sampler -- contains things that don't change on each iteration
  setup = function(model, mvSaved, target, control) {
    
    # model dimensions
    p <- control$p
    p_beta <- p + 1
    
    nobs <- as.integer(control$nobs)
    ncen <- as.integer(control$ncen)
    n <- nobs + ncen
    
    # prior distribution precision: Sigma_{beta}^{-1}
    s_beta_inv <- inverse(model$beta_cov[1:p_beta, 1:p_beta])
    
    # fixed contribution of prior hypterparameters on posterior mean: Sigma_{beta}^{-1} * beta_0
    beta_prior_num <- (s_beta_inv %*% asCol(model$beta_mean[1:p_beta]))
    
    # recalculate nodes depending on changes in beta
    calcNodes <- model$getDependencies(target)
  },
  
  ## runs at every MCMC iteration
  run = function() {
    
    ## build current iteration's design matrix (X,Z) and outcome Y
    W <- matrix(0, nrow = n, ncol = p_beta)
    Y <- numeric(n)
    
    ## observed subjects
    if (nobs > 0) {
      for (k in 1:nobs) {
        W[k, 1] <- model$x_obs[k]
        
        if (p > 1) {
          for (j in 1:p) {
            W[k, j+1] <- model$z_obs[k, j]
          }
        } else if(p == 1){
          W[k, p+1] <- model$z_obs[k]
        }
        
        Y[k] <- model$y_obs[k]
      }
    }
    
    ## censored subjects
    if (ncen > 0) {
      for (k in 1:ncen) {
        i <- nobs + ks
        
        W[i, 1] <- model$x_cen[k]
        
        if (p > 1) {
          for (j in 1:p) {
            W[i, j+1] <- model$z_cen[k, j]
          }
        } else if(p == 1) {
          W[i, p+1] <- model$z_cen[k]
        }
        
        
        Y[i] <- model$y_cen[k]
      }
    }
    
    ## full conditional for beta
    
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

## custom sampler with only observed subjects ####
sampler_beta_conjugate_block_obs <- nimbleFunction(
  contains = sampler_BASE,
  
  ## runs when nimble constructs the sampler -- contains things that don't change on each iteration
  setup = function(model, mvSaved, target, control) {
    
    # model dimensions
    p <- control$p
    p_beta <- p + 1
    
    nobs <- as.integer(control$nobs)
    n <- nobs
    
    # prior distribution precision: Sigma_{beta}^{-1}
    s_beta_inv <- inverse(model$beta_cov[1:p_beta, 1:p_beta])
    
    # fixed contribution of prior hypterparameters on posterior mean: Sigma_{beta}^{-1} * beta_0
    beta_prior_num <- (s_beta_inv %*% asCol(model$beta_mean[1:p_beta]))
    
    # recalculate nodes depending on changes in beta
    calcNodes <- model$getDependencies(target)
  },
  
  ## runs at every MCMC iteration
  run = function() {
    
    ## build current iteration's design matrix (X,Z) and outcome Y
    W <- matrix(0, nrow = n, ncol = p_beta)
    Y <- numeric(n)
    
    ## observed subjects
    if (nobs > 0) {
      for (k in 1:nobs) {
        W[k, 1] <- model$x_obs[k]
        
        if (p > 1) {
          for (j in 1:p) {
            W[k, j+1] <- model$z_obs[k, j]
          }
        } else if(p == 1){
          W[k, p+1] <- model$z_obs[k, p]
        }
        
        Y[k] <- model$y_obs[k]
      }
    }
    
    ## full conditional for beta
    
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

## custom sampler for data in which ALL subjects are censored ####
sampler_beta_conjugate_block_cen <- nimbleFunction(
  contains = sampler_BASE,
  
  ## runs when nimble constructs the sampler -- contains things that don't change on each iteration
  setup = function(model, mvSaved, target, control) {
    
    # model dimensions
    p <- control$p
    p_beta <- p + 1
    
    ncen <- as.integer(control$ncen)
    n <- ncen
    
    # prior distribution precision: Sigma_{beta}^{-1}
    s_beta_inv <- inverse(model$beta_cov[1:p_beta, 1:p_beta])
    
    # fixed contribution of prior hypterparameters on posterior mean: Sigma_{beta}^{-1} * beta_0
    beta_prior_num <- (s_beta_inv %*% asCol(model$beta_mean[1:p_beta]))
    
    # recalculate nodes depending on changes in beta
    calcNodes <- model$getDependencies(target)
  },
  
  ## runs at every MCMC iteration
  run = function() {
    
    ## build current iteration's design matrix (X,Z) and outcome Y
    W <- matrix(0, nrow = n, ncol = p_beta)
    Y <- numeric(n)
    
    ## censored subjects
    if (ncen > 0) {
      for (k in 1:ncen) {
        i <- nobs + ks
        
        W[i, 1] <- model$x_cen[k]
        
        if (p > 0) {
          for (j in 1:p) {
            W[i, j+1] <- model$z_cen[k, j]
          }
        } else if(p == 1) {
          W[i, p+1] <- model$z_cen[k, p]
        }
        
        
        Y[i] <- model$y_cen[k]
      }
    }
    
    ## full conditional for beta
    
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

## mixed - intercept only ####
sampler_beta_conjugate_block_mixed_intercept <- nimbleFunction(
  contains = sampler_BASE,
  
  ## runs when nimble constructs the sampler -- contains things that don't change on each iteration
  setup = function(model, mvSaved, target, control) {
    
    # model dimensions
    p <- control$p
    p_beta <- p + 1
    
    nobs <- as.integer(control$nobs)
    ncen <- as.integer(control$ncen)
    n <- nobs + ncen
    
    # prior distribution precision: Sigma_{beta}^{-1}
    s_beta_inv <- inverse(model$beta_cov[1:p_beta, 1:p_beta])
    
    # fixed contribution of prior hypterparameters on posterior mean: Sigma_{beta}^{-1} * beta_0
    beta_prior_num <- (s_beta_inv %*% asCol(model$beta_mean[1:p_beta]))
    
    # recalculate nodes depending on changes in beta
    calcNodes <- model$getDependencies(target)
  },
  
  ## runs at every MCMC iteration
  run = function() {
    
    ## build current iteration's design matrix (X,Z) and outcome Y
    W <- matrix(0, nrow = n, ncol = p_beta)
    Y <- numeric(n)
    
    ## observed subjects
    if (nobs > 0) {
      for (k in 1:nobs) {
        W[k, 1] <- model$x_obs[k]
        W[k, p+1] <- model$z_obs[k]
        Y[k] <- model$y_obs[k]
      }
    }
    
    ## censored subjects
    if (ncen > 0) {
      for (k in 1:ncen) {
        i <- nobs + k
        
        W[i, 1] <- model$x_cen[k]
        W[i, p+1] <- model$z_cen[k]
        Y[i] <- model$y_cen[k]
      }
    }
    
    ## full conditional for beta
    
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


## censored only - intercept only ####
sampler_beta_conjugate_block_cens_intercept <- nimbleFunction(
  contains = sampler_BASE,
  
  ## runs when nimble constructs the sampler -- contains things that don't change on each iteration
  setup = function(model, mvSaved, target, control) {
    
    # model dimensions
    p <- control$p
    p_beta <- p + 1
    
    ncen <- as.integer(control$ncen)
    n <- ncen
    
    # prior distribution precision: Sigma_{beta}^{-1}
    s_beta_inv <- inverse(model$beta_cov[1:p_beta, 1:p_beta])
    
    # fixed contribution of prior hypterparameters on posterior mean: Sigma_{beta}^{-1} * beta_0
    beta_prior_num <- (s_beta_inv %*% asCol(model$beta_mean[1:p_beta]))
    
    # recalculate nodes depending on changes in beta
    calcNodes <- model$getDependencies(target)
  },
  
  ## runs at every MCMC iteration
  run = function() {
    
    ## build current iteration's design matrix (X,Z) and outcome Y
    W <- matrix(0, nrow = n, ncol = p_beta)
    Y <- numeric(n)
    
    ## censored subjects
    if (ncen > 0) {
      for (k in 1:ncen) {
        i <- k
        
        W[i, 1] <- model$x_cen[k]
        W[i, p+1] <- model$z_cen[k]
        Y[i] <- model$y_cen[k]
      }
    }
    
    ## full conditional for beta
    
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
