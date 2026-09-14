dg_glm <- function(
    n = 1000,
    outcome_formula = Y ~ X + Z2 + Z3,
    covariate_formula = X ~ Z2 + Z3,
    family = gaussian(link = 'identity'),
    
    beta = c(1, 0.7, 1, -1),
    
    censoring = 0.80,
    z_spec = list(
      Z2 = list(dist = 'bernoulli',
                prob = 0.8),
      Z3 = list(dist = 'normal',
                mean = 2, sd = 1)
    )) {
  
  outcome_name <- all.vars(outcome_formula)[1]
  x_name <- all.vars(covariate_formula)[1]
  
  # create Z - uncensored covariates
  gen_z <- function(spec, n) {
    switch(
      spec$dist,
      "normal" = {rnorm(n, mean = spec$mean, sd = spec$sd)},
      "bernoulli" = {rbinom(n, size = 1, prob = spec$prob)},
      "uniform" = {runif(n, min = spec$min, max = spec$max)}
    )
  }
  Z <- as.data.frame(lapply(z_spec, gen_z, n = n))
  names(Z) <- names(z_spec)
  
  # covariate model design matrix
  x_terms <- delete.response(terms(covariate_formula))
  Zx <- model.matrix(x_terms, data = Z)
  pzx <- ncol(Zx)
  
  # DPMM parameters for X|Z
  comps <- 3
  w <- c(0.3, 0.5, 0.2)
  x_sd <- c(1, 1.2, 1.5)
  gamma <- matrix(c(0.5, 2, 1.5,
                    0.7, 3, 2,
                    0.9, 4, 3),
                  nrow = comps,
                  ncol = pzx,
                  byrow = TRUE,
                  dimnames = list(NULL, colnames(Zx)))
  
  # generate X
  cluster <- sample(seq_len(comps),
                    size = n,
                    replace = TRUE,
                    prob = w)
  
  x_mean <- rowSums(Zx * gamma[cluster, , drop = FALSE])
  X <- rnorm(n, mean = x_mean, sd = x_sd[cluster])
  
  model_data <- Z
  model_data[[x_name]] <- X
  
  # outcome design matrix
  y_terms <- delete.response(terms(outcome_formula))
  Zy <- model.matrix(y_terms, data = model_data)
  pzy <- ncol(Zy)
  
  eta <- as.numeric(Zy %*% beta)
  
  # generate outcome
  mu <- family$linkinv(eta)
  tau <- 1
  family_name <- family$family
  
  if (family_name == "gaussian") {
    Y <- rnorm(n, mean = mu, sd = 1/sqrt(tau))
  } else if (family_name == "binomial") {
    Y <- rbinom(n, size = 1, prob = mu)
  } else if (family_name == "Gamma") {
    Y <- rgamma(n, shape = tau, rate = tau / mu)
  } else if (family_name == "poisson") {
    Y <- rpois(n, lambda = mu)
  }
  
  # censoring mechanism
  beta_omega <- matrix(c(0, 1.5, -0.5,
                         0, 0.5, 0),
                       nrow = 2, byrow = TRUE)
  rho <- c(0.90, 0.95, 0.99)
  
  sigma_start <- 0.082
  sigma_width <- 0.15
  jitter_sd <- 0.082
  
  censoring_key <- as.character(censoring)
  censoring_pars <- switch(
    censoring_key,
    "0.2" = list(
      alpha_start = c(1.8, 0.8, 0.4),
      alpha_width = c(1.4, 2, -0.1),
      eps_obs = 1.4
    ),
    "0.4" = list(
      alpha_start = c(2, 1.0, 0.6),
      alpha_width = c(1.3, 1.5, -0.1),
      eps_obs = 1.1
    ),
    "0.8" = list(
      alpha_start = c(2.3, 1.4, 1.1),
      alpha_width = c(0.9, 1, -0.3),
      eps_obs = 1.0
    ),
  )
  alpha_start <- censoring_pars$alpha_start
  alpha_width <- censoring_pars$alpha_width
  eps_obs <- censoring_pars$eps_obs
  
  Zc <- model.matrix( ~ Z2 + Z3, data = Z)
  
  omega_high <- as.numeric(Zc %*% beta_omega[1, ])
  omega_med <- as.numeric(Zc %*% beta_omega[2, ])
  omega_low <- rep(0, n)
  
  exp_omega <- cbind(exp(omega_high),
                     exp(omega_med),
                     exp(omega_low))
  
  pi_obs <- exp_omega / rowSums(exp_omega)
  nu <- apply(pi_obs, 1, function(p) {
    sample(1:3, size = 1, prob = p)
  })
  
  # observation window
  T0 <- rnorm(n, mean = as.numeric(Zc %*% alpha_start), sd = sigma_start)
  W <- rlnorm(n, meanlog = as.numeric(Zc %*% alpha_width), sdlog = sigma_width)
  T1 <- T0 + W
  
  CL <- CR <- numeric(n)
  K_planned <- K_obs <- integer(n)
  V_star <- V_tilde <- V <- r <- vector("list", n)
  
  ## generate visit data
  for (i in seq_len(n)) {
    
    # initialize empty (zero-length) storage for visit vectors and observation indicator vector
    V_star_i <- V_tilde_i <- V_i <- numeric(0)       
    r_i <- integer(0)
    
    K_planned[i] <- floor(W[i])                                        # number of visits if none missed
    
    if (K_planned[i] > 0) {                         
      V_star_i <- T0[i] + seq_len(K_planned[i])                        # visits start as every year from observation window start
      
      V_tilde_i <- V_star_i + rnorm(K_planned[i], mean = 0, sd = jitter_sd)     # add noise to annual visit
      V_tilde_i <- sort(pmin(pmax(V_tilde_i, T0[i]), T1[i]))                    # sort visits in case ordering switched from added noise
      
      r_i <- rbinom(K_planned[i], size = 1, prob = rho[nu[i]])                  # sample observation indicator
      V_i <- V_tilde_i[r_i == 1]                                                # mask visits with observation indicator
    }
    
    # add in T0 as a forced observed visit (if a subject never comes to the enrollment visit they aren't in the study anyway so this makes sense)'
    V_star_i <- sort(c(T0[i], V_star_i))
    V_tilde_i <- sort(c(T0[i], V_tilde_i))
    r_i <- c(1, r_i)
    V_i <- sort(c(T0[i], V_i))
    
    K_obs[i] <- length(V_i)                                                     # number of observed visits
    V[[i]] <- V_i                                                               # observed visit vector
    V_star[[i]] <- V_star_i                                                     # planned visit vector
    V_tilde[[i]] <- V_tilde_i                                                   # jittered planned visit vector
    r[[i]] <- r_i                                                               # visit observation indicator vector
    
    ## covariate censoring indicators and endpoints
    
    if (X[i] < V_i[1]) {                                                        # if X is left censored
      CL[i] <- -Inf
      CR[i] <- V_i[1]
    } else if (X[i] > V_i[K_obs[i]]) {                                          # if X is right censored
      CL[i] <- V_i[K_obs[i]]
      CR[i] <- Inf
    } else {                                                                    # X either interval or uncensored
      right_idx <- which(V_i >= X[i])[1]                                        # smallest visit greater than X
      
      if (X[i] == V_i[right_idx]) {                                             # if X falls exactly on a visit
        CL[i] <- CR[i] <- X[i]
      } else {                                                                  # if X is interval censored
        left <- V_i[right_idx - 1]
        right <- V_i[right_idx]
        
        if ((right - left) <= eps_obs) {                                        # if X is interval censored with endpoints within tolerance, define as uncensored
          CL[i] <- CR[i] <- X[i]
        } else {
          CL[i] <- left
          CR[i] <- right
        }
      }
    }
  }
  
  dat <- Z
  dat[[outcome_name]] <- Y
  dat$CL <- CL
  dat$CR <- CR
  
  return(list(data = dat,
              
              X = X,
              Y = Y,
              Zx = Zx,
              Zy = Zy,
              
              cluster = cluster,
              
              CL = CL,
              CR = CR
  ))
}