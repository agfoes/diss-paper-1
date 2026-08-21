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
    beta_prior_num <- (s_beta_inv %>% model$beta_mean[1:p_beta])[, 1]
    
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
        
        for (j in 1:p) {
          W[k, j+1] <- model$z_obs[k, j]
        }
        
        Y[k] <- model$y_obs[k]
      }
    }
    
    ## censored subjects
    if (ncen > 0) {
      for (k in 1:ncen) {
        i <- nobs + k         # all the uncensored subjects are first in order, then add censored subjects after (lose relative position in original vector but this shouldn't matter for this)
        
        W[i, 1] <- model$x_cen[k]
        
        for (j in 1:p) {
          W[i, j+1] <- model$z_cen[k, j]
        }
        
        Y[i] <- model$y_cen[k]
      }
    }
    
    ## full conditional for beta
    
    # V_beta = tau * W^T W + Sigma_beta^{-1}
    V_beta <- s_beta_inv + model$tau * (t(W) %*% W)
    
    # W^T Y
    tWY <- (t(W)%*%Y)[, 1]
    
    # tau * W^T Y + Sigma_beta^{-1} beta_0
    mu_beta_num <- model$tau * tWY + beta_prior_num
    
    # V_beta^{-1}
    V_beta_inv <- inverse(V_beta)
    
    # mu_beta = V_beta^{-1} mu_beta_num
    mu_beta <- (V_beta_inv %*% mu_beta_num)[, 1]
    
    ## joint gibbs draw
    V_beta_cholesky <- chol(V_beta)
    
    beta_new <- rmnorm_chol(1, mean = mu_beta,
                            cholesky = V_beta_cholesky, prec_param = TRUE)
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