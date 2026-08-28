
#' find_L
#' 
#' @description
#' This function finds the recommended truncation value for the DPMM clusters
#' based on sample size. See Ohlssen et al. (2007) for more details.
#' 
#'
#' @returns
#' @export
#'
#' @examples
find_L <- function(alpha_prob = 0.99,
                   alpha_shape = 1,
                   alpha_rate = 1,
                   eps = 0.01
                   ) {
  alpha_max <- qgamma(alpha_prob, shape = alpha_shape, rate = alpha_rate)
  
  L <- ceiling(1 + log(eps) / log(alpha_max / (1 + alpha_max)))
  
  return(L)
}



#' datagen
#'
#' @param n 
#' @param p_z2 
#' @param mu_z3 
#' @param sd_z3 
#' @param comps 
#' @param w 
#' @param gamma 
#' @param x_sd 
#' @param beta_x 
#' @param beta_z 
#' @param tau 
#' @param beta_omega 
#' @param alpha_start 
#' @param alpha_width 
#' @param sigma_start 
#' @param sigma_width 
#' @param jitter_sd 
#' @param rho 
#' @param eps_obs 
#'
#' @returns
#' @export
#'
#' @examples
datagen <- function(n = 100,                                                    # sample size
                    p_z2 = 0.8,                                                 # observed covariate Z2 ~ bern(p_z2)
                    mu_z3 = 5, sd_z3 = 1,                                       # mean and sd for observed covariate Z3 ~ N(mu_z3, sd_z3)
                    comps = 3, w = c(0.3, 0.5, 0.2),                            # number of components for mixture X|Z and probabilities for each component
                    gamma = matrix(c(0, 3, 1, 
                                     0, 1, 5, 
                                     0, 2, 2), 
                                   nrow = 3, byrow = TRUE)*365.25,              # linear coeff. for X|Z ~ N(Z*gamma, x_sd)
                    x_sd = c(0.5, 1, 1)*365.25,
                    beta_x = 0.4/365.25,                                        # coefficients for outcome linear predictor on X
                    beta_z = c(0.2, 0.1, 0.2),                                  # coefficients for outcome linear predictor on Z
                    tau = 1,                                                    # outcome precision parameter
                    beta_omega = matrix(c(0, 0.5, 0.25,
                                          0, 0.25, 0.10),
                                        nrow = 2, byrow = TRUE),                # coefficients for visit generation observation category
                    alpha_start = c(5, 2, -0.5)*365.25,                         # parameter for start of observation window T0 ~ N(alpha * Z, sigma_start)
                    alpha_width = c(log(6*365.25), 0.25, -0.05),                # parameter for width of observation window log(W) ~ N(alpha * Z, sigma_width)
                    sigma_start = 1, sigma_width = 0.15,
                    jitter_sd = 30,                                             # sd for random noise added to each visit
                    rho = c(0.85, 0.9, 0.99),                                   # vector of probabilities for observing a visit
                    eps_obs = 350                                               # tolerance window for defining uncensored observations in X
) {
  
  ## observed covariates
  Z1 <- rep(1, n)                                                               # intercept
  Z2 <- rbinom(n = n, size = 1, p_z2)                                           # Z2 ~ binomial(p_z2)
  Z3 <- rnorm(n, mean = mu_z3, sd = sd_z3)                                      # Z3 ~ normal(mu_z3, sd_z3)
  Z <- cbind(Z1, Z2, Z3)
  
  ## censored covariate
  cluster <- sample(1:comps, n, replace = TRUE, prob = w)                       # sample X|Z cluster assignment from weights w
  means <- rowSums(Z * gamma[cluster, 1:3])                                     # mean of X|Z determined by cluster assignment and gamma row
  X <- rnorm(n, mean = means, sd = x_sd[cluster])                               # X|Z ~ normal(mean[cluster], sd[cluster])
  
  ## outcome
  eta <- beta_x * X + Z %*% beta_z
  Y <- rnorm(n, eta, sd = 1/sqrt(tau))                                          # Y|X,Z ~ norm(eta, 1/tau)
  
  ## probability of observing a visit
  omega_high <- as.numeric(Z %*% beta_omega[1, ])                               # linear pred if subjects have high observation probability
  omega_med <- as.numeric(Z %*% beta_omega[2, ])                                # linear pred if subjects have medium observation probability
  omega_low <- rep(0, n)                                                        # linear pred if subjects have low observation probability
  
  exp_omega <- cbind(exp(omega_high), exp(omega_med), exp(omega_low))           # exponentiate linear predictors
  pi_obs <- exp_omega / rowSums(exp_omega)                                      # calculate probabilities for observation probability assignment
  
  nu <- apply(pi_obs, 1, function(p) {                                          # sample observation categories with ^ probabilities per subject
    sample(1:3, size = 1, prob = p)
  })
  
  ## observation window
  T0 <- rnorm(n, mean = as.numeric(Z %*% alpha_start), sd = sigma_start)        # observation window start visit
  W <- rlnorm(n, meanlog = as.numeric(Z %*% alpha_width), sd = sigma_width)     # width of observation window
  T1 <- T0 + W
  
  ## storage
  CL <- CR <- numeric(n)
  DL <- DR <- Dobs <- Dnone <- integer(n)
  K_planned <- K_obs <- integer(n)
  V_star <- V_tilde <- V <- r <- vector('list', n)
  
  ## generate visit data
  for (i in seq_len(n)) {
    
    # initialize empty (zero-length) storage for visit vectors and observation indicator vector
    V_star_i <- V_tilde_i <- V_i <- numeric(0)       
    r_i <- integer(0)
    
    K_planned[i] <- floor(W[i] / 365.25)                                        # number of visits if none missed
    
    if (K_planned[i] > 0) {                         
      V_star_i <- T0[i] + 365.25 * seq_len(K_planned[i])                        # visits start as every year from observation window start
      
      V_tilde_i <- V_star_i + rnorm(K_planned[i], mean = 0, sd = jitter_sd)     # add noise to annual visit
      V_tilde_i <- sort(pmin(pmax(V_tilde_i, T0[i]), T1[i]))                    # sort visits in case ordering switched from added noise
      
      r_i <- rbinom(K_planned[i], size = 1, prob = rho[nu[i]])                  # sample observation indicator
      V_i <- V_tilde_i[r_i == 1]                                                # mask visits with observation idicator
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
    
    if (length(V_i) == 0) {                                                     # if subject has NO observed visits
      CL[i] <- -Inf
      CR[i] <- Inf
      Dnone[i] <- 1
    } else if (X[i] < V_i[1]) {                                                 # if X is left censored
      CL[i] <- -Inf
      CR[i] <- V_i[1]
      DL[i] <- 1
    } else if (X[i] > V_i[K_obs[i]]) {                                          # if X is right censored
      CL[i] <- V_i[K_obs[i]]
      CR[i] <- Inf
      DR[i] <- 1
    } else {                                                                    # X either interval or uncensored
      right_idx <- which(V_i >= X[i])[1]                                        # smallest visit greater than X
      
      if (X[i] == V_i[right_idx]) {                                             # if X falls exactly on a visit
        CL[i] <- CR[i] <- X[i]
        Dobs[i] <- 1
      } else {                                                                  # if X is interval censored
        left <- V_i[right_idx - 1]
        right <- V_i[right_idx]
        
        if ((right - left) <= eps_obs) {                                        # if X is interval censored with endpoints within tolerance, define as uncensored
          CL[i] <- CR[i] <- X[i]
          Dobs[i] <- 1
        } else {
          CL[i] <- left
          CR[i] <- right
        }
      }
    }
  }
  
  return(list(
    Z1 = Z1,
    Z2 = Z2,
    Z3 = Z3,
    X = X,
    Y = Y,
    cluster = cluster,
    
    T0 = T0,
    T1 = T1,
    W = W,
    nu = nu,
    K_planned = K_planned,
    K_obs = K_obs,
    V_star = V_star,
    V_tilde = V_tilde,
    V = V,
    r = r,
    rho = rho,
    
    CL = CL,
    CR = CR,
    DL = DL,
    DR = DR,
    Dobs = Dobs,
    Dnone = Dnone
  ))
}

#' datagen_bernoulli
#'
#' @param n 
#' @param p_z2 
#' @param mu_z3 
#' @param sd_z3 
#' @param comps 
#' @param w 
#' @param gamma 
#' @param x_sd 
#' @param beta_x 
#' @param beta_z 
#' @param tau 
#' @param beta_omega 
#' @param alpha_start 
#' @param alpha_width 
#' @param sigma_start 
#' @param sigma_width 
#' @param jitter_sd 
#' @param rho 
#' @param eps_obs 
#'
#' @returns
#' @export
#'
#' @examples
datagen_bernoulli <- function(n = 100,                                                    # sample size
                    p_z2 = 0.8,                                                 # observed covariate Z2 ~ bern(p_z2)
                    mu_z3 = 5, sd_z3 = 1,                                       # mean and sd for observed covariate Z3 ~ N(mu_z3, sd_z3)
                    comps = 3, w = c(0.3, 0.5, 0.2),                            # number of components for mixture X|Z and probabilities for each component
                    gamma = matrix(c(0, 3, 1, 
                                     0, 1, 5, 
                                     0, 2, 2), 
                                   nrow = 3, byrow = TRUE)*365.25,              # linear coeff. for X|Z ~ N(Z*gamma, x_sd)
                    x_sd = c(0.5, 1, 1)*365.25,
                    beta_x = 0.4/365.25,                                        # coefficients for outcome linear predictor on X
                    beta_z = c(0.2, 0.1, 0.2),                                  # coefficients for outcome linear predictor on Z
                    tau = 1,                                                    # outcome precision parameter
                    beta_omega = matrix(c(0, 0.5, 0.25,
                                          0, 0.25, 0.10),
                                        nrow = 2, byrow = TRUE),                # coefficients for visit generation observation category
                    alpha_start = c(5, 2, -0.5)*365.25,                         # parameter for start of observation window T0 ~ N(alpha * Z, sigma_start)
                    alpha_width = c(log(6*365.25), 0.25, -0.05),                # parameter for width of observation window log(W) ~ N(alpha * Z, sigma_width)
                    sigma_start = 1, sigma_width = 0.15,
                    jitter_sd = 30,                                             # sd for random noise added to each visit
                    rho = c(0.85, 0.9, 0.99),                                   # vector of probabilities for observing a visit
                    eps_obs = 350                                               # tolerance window for defining uncensored observations in X
) {
  
  ## observed covariates
  Z1 <- rep(1, n)                                                               # intercept
  Z2 <- rbinom(n = n, size = 1, p_z2)                                           # Z2 ~ binomial(p_z2)
  Z3 <- rnorm(n, mean = mu_z3, sd = sd_z3)                                      # Z3 ~ normal(mu_z3, sd_z3)
  Z <- cbind(Z1, Z2, Z3)
  
  ## censored covariate
  cluster <- sample(1:comps, n, replace = TRUE, prob = w)                       # sample X|Z cluster assignment from weights w
  means <- rowSums(Z * gamma[cluster, 1:3])                                     # mean of X|Z determined by cluster assignment and gamma row
  X <- rnorm(n, mean = means, sd = x_sd[cluster])                               # X|Z ~ normal(mean[cluster], sd[cluster])
  
  ## outcome
  eta <- beta_x * X + Z %*% beta_z
  theta <- plogis(eta)
  Y <- rbinom(n, size = 1, prob = theta)                                        # Y|X,Z ~ binomial(n, logit(eta))
  
  ## probability of observing a visit
  omega_high <- as.numeric(Z %*% beta_omega[1, ])                               # linear pred if subjects have high observation probability
  omega_med <- as.numeric(Z %*% beta_omega[2, ])                                # linear pred if subjects have medium observation probability
  omega_low <- rep(0, n)                                                        # linear pred if subjects have low observation probability
  
  exp_omega <- cbind(exp(omega_high), exp(omega_med), exp(omega_low))           # exponentiate linear predictors
  pi_obs <- exp_omega / rowSums(exp_omega)                                      # calculate probabilities for observation probability assignment
  
  nu <- apply(pi_obs, 1, function(p) {                                          # sample observation categories with ^ probabilities per subject
    sample(1:3, size = 1, prob = p)
  })
  
  ## observation window
  T0 <- rnorm(n, mean = as.numeric(Z %*% alpha_start), sd = sigma_start)        # observation window start visit
  W <- rlnorm(n, meanlog = as.numeric(Z %*% alpha_width), sd = sigma_width)     # width of observation window
  T1 <- T0 + W
  
  ## storage
  CL <- CR <- numeric(n)
  DL <- DR <- Dobs <- Dnone <- integer(n)
  K_planned <- K_obs <- integer(n)
  V_star <- V_tilde <- V <- r <- vector('list', n)
  
  ## generate visit data
  for (i in seq_len(n)) {
    
    # initialize empty (zero-length) storage for visit vectors and observation indicator vector
    V_star_i <- V_tilde_i <- V_i <- numeric(0)       
    r_i <- integer(0)
    
    K_planned[i] <- floor(W[i] / 365.25)                                        # number of visits if none missed
    
    if (K_planned[i] > 0) {                         
      V_star_i <- T0[i] + 365.25 * seq_len(K_planned[i])                        # visits start as every year from observation window start
      
      V_tilde_i <- V_star_i + rnorm(K_planned[i], mean = 0, sd = jitter_sd)     # add noise to annual visit
      V_tilde_i <- sort(pmin(pmax(V_tilde_i, T0[i]), T1[i]))                    # sort visits in case ordering switched from added noise
      
      r_i <- rbinom(K_planned[i], size = 1, prob = rho[nu[i]])                  # sample observation indicator
      V_i <- V_tilde_i[r_i == 1]                                                # mask visits with observation idicator
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
    
    if (length(V_i) == 0) {                                                     # if subject has NO observed visits
      CL[i] <- -Inf
      CR[i] <- Inf
      Dnone[i] <- 1
    } else if (X[i] < V_i[1]) {                                                 # if X is left censored
      CL[i] <- -Inf
      CR[i] <- V_i[1]
      DL[i] <- 1
    } else if (X[i] > V_i[K_obs[i]]) {                                          # if X is right censored
      CL[i] <- V_i[K_obs[i]]
      CR[i] <- Inf
      DR[i] <- 1
    } else {                                                                    # X either interval or uncensored
      right_idx <- which(V_i >= X[i])[1]                                        # smallest visit greater than X
      
      if (X[i] == V_i[right_idx]) {                                             # if X falls exactly on a visit
        CL[i] <- CR[i] <- X[i]
        Dobs[i] <- 1
      } else {                                                                  # if X is interval censored
        left <- V_i[right_idx - 1]
        right <- V_i[right_idx]
        
        if ((right - left) <= eps_obs) {                                        # if X is interval censored with endpoints within tolerance, define as uncensored
          CL[i] <- CR[i] <- X[i]
          Dobs[i] <- 1
        } else {
          CL[i] <- left
          CR[i] <- right
        }
      }
    }
  }
  
  return(list(
    Z1 = Z1,
    Z2 = Z2,
    Z3 = Z3,
    X = X,
    Y = Y,
    cluster = cluster,
    
    T0 = T0,
    T1 = T1,
    W = W,
    nu = nu,
    K_planned = K_planned,
    K_obs = K_obs,
    V_star = V_star,
    V_tilde = V_tilde,
    V = V,
    r = r,
    rho = rho,
    
    CL = CL,
    CR = CR,
    DL = DL,
    DR = DR,
    Dobs = Dobs,
    Dnone = Dnone
  ))
}

#' datagen_poisson
#'
#' @param n 
#' @param p_z2 
#' @param mu_z3 
#' @param sd_z3 
#' @param comps 
#' @param w 
#' @param gamma 
#' @param x_sd 
#' @param beta_x 
#' @param beta_z 
#' @param tau 
#' @param beta_omega 
#' @param alpha_start 
#' @param alpha_width 
#' @param sigma_start 
#' @param sigma_width 
#' @param jitter_sd 
#' @param rho 
#' @param eps_obs 
#'
#' @returns
#' @export
#'
#' @examples
datagen_poisson <- function(n = 100,                                                    # sample size
                              p_z2 = 0.8,                                                 # observed covariate Z2 ~ bern(p_z2)
                              mu_z3 = 5, sd_z3 = 1,                                       # mean and sd for observed covariate Z3 ~ N(mu_z3, sd_z3)
                              comps = 3, w = c(0.3, 0.5, 0.2),                            # number of components for mixture X|Z and probabilities for each component
                              gamma = matrix(c(0, 3, 1, 
                                               0, 1, 5, 
                                               0, 2, 2), 
                                             nrow = 3, byrow = TRUE)*365.25,              # linear coeff. for X|Z ~ N(Z*gamma, x_sd)
                              x_sd = c(0.5, 1, 1)*365.25,
                              beta_x = 0.4/365.25,                                        # coefficients for outcome linear predictor on X
                              beta_z = c(0.2, 0.1, 0.2),                                  # coefficients for outcome linear predictor on Z
                              tau = 1,                                                    # outcome precision parameter
                              beta_omega = matrix(c(0, 0.5, 0.25,
                                                    0, 0.25, 0.10),
                                                  nrow = 2, byrow = TRUE),                # coefficients for visit generation observation category
                              alpha_start = c(5, 2, -0.5)*365.25,                         # parameter for start of observation window T0 ~ N(alpha * Z, sigma_start)
                              alpha_width = c(log(6*365.25), 0.25, -0.05),                # parameter for width of observation window log(W) ~ N(alpha * Z, sigma_width)
                              sigma_start = 1, sigma_width = 0.15,
                              jitter_sd = 30,                                             # sd for random noise added to each visit
                              rho = c(0.85, 0.9, 0.99),                                   # vector of probabilities for observing a visit
                              eps_obs = 350                                               # tolerance window for defining uncensored observations in X
) {
  
  ## observed covariates
  Z1 <- rep(1, n)                                                               # intercept
  Z2 <- rbinom(n = n, size = 1, p_z2)                                           # Z2 ~ binomial(p_z2)
  Z3 <- rnorm(n, mean = mu_z3, sd = sd_z3)                                      # Z3 ~ normal(mu_z3, sd_z3)
  Z <- cbind(Z1, Z2, Z3)
  
  ## censored covariate
  cluster <- sample(1:comps, n, replace = TRUE, prob = w)                       # sample X|Z cluster assignment from weights w
  means <- rowSums(Z * gamma[cluster, 1:3])                                     # mean of X|Z determined by cluster assignment and gamma row
  X <- rnorm(n, mean = means, sd = x_sd[cluster])                               # X|Z ~ normal(mean[cluster], sd[cluster])
  
  ## outcome
  eta <- beta_x * X + Z %*% beta_z
  theta <- exp(eta)
  Y <- rpois(n, lambda = theta)                                        # Y|X,Z ~ binomial(n, logit(eta))
  
  ## probability of observing a visit
  omega_high <- as.numeric(Z %*% beta_omega[1, ])                               # linear pred if subjects have high observation probability
  omega_med <- as.numeric(Z %*% beta_omega[2, ])                                # linear pred if subjects have medium observation probability
  omega_low <- rep(0, n)                                                        # linear pred if subjects have low observation probability
  
  exp_omega <- cbind(exp(omega_high), exp(omega_med), exp(omega_low))           # exponentiate linear predictors
  pi_obs <- exp_omega / rowSums(exp_omega)                                      # calculate probabilities for observation probability assignment
  
  nu <- apply(pi_obs, 1, function(p) {                                          # sample observation categories with ^ probabilities per subject
    sample(1:3, size = 1, prob = p)
  })
  
  ## observation window
  T0 <- rnorm(n, mean = as.numeric(Z %*% alpha_start), sd = sigma_start)        # observation window start visit
  W <- rlnorm(n, meanlog = as.numeric(Z %*% alpha_width), sd = sigma_width)     # width of observation window
  T1 <- T0 + W
  
  ## storage
  CL <- CR <- numeric(n)
  DL <- DR <- Dobs <- Dnone <- integer(n)
  K_planned <- K_obs <- integer(n)
  V_star <- V_tilde <- V <- r <- vector('list', n)
  
  ## generate visit data
  for (i in seq_len(n)) {
    
    # initialize empty (zero-length) storage for visit vectors and observation indicator vector
    V_star_i <- V_tilde_i <- V_i <- numeric(0)       
    r_i <- integer(0)
    
    K_planned[i] <- floor(W[i] / 365.25)                                        # number of visits if none missed
    
    if (K_planned[i] > 0) {                         
      V_star_i <- T0[i] + 365.25 * seq_len(K_planned[i])                        # visits start as every year from observation window start
      
      V_tilde_i <- V_star_i + rnorm(K_planned[i], mean = 0, sd = jitter_sd)     # add noise to annual visit
      V_tilde_i <- sort(pmin(pmax(V_tilde_i, T0[i]), T1[i]))                    # sort visits in case ordering switched from added noise
      
      r_i <- rbinom(K_planned[i], size = 1, prob = rho[nu[i]])                  # sample observation indicator
      V_i <- V_tilde_i[r_i == 1]                                                # mask visits with observation idicator
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
    
    if (length(V_i) == 0) {                                                     # if subject has NO observed visits
      CL[i] <- -Inf
      CR[i] <- Inf
      Dnone[i] <- 1
    } else if (X[i] < V_i[1]) {                                                 # if X is left censored
      CL[i] <- -Inf
      CR[i] <- V_i[1]
      DL[i] <- 1
    } else if (X[i] > V_i[K_obs[i]]) {                                          # if X is right censored
      CL[i] <- V_i[K_obs[i]]
      CR[i] <- Inf
      DR[i] <- 1
    } else {                                                                    # X either interval or uncensored
      right_idx <- which(V_i >= X[i])[1]                                        # smallest visit greater than X
      
      if (X[i] == V_i[right_idx]) {                                             # if X falls exactly on a visit
        CL[i] <- CR[i] <- X[i]
        Dobs[i] <- 1
      } else {                                                                  # if X is interval censored
        left <- V_i[right_idx - 1]
        right <- V_i[right_idx]
        
        if ((right - left) <= eps_obs) {                                        # if X is interval censored with endpoints within tolerance, define as uncensored
          CL[i] <- CR[i] <- X[i]
          Dobs[i] <- 1
        } else {
          CL[i] <- left
          CR[i] <- right
        }
      }
    }
  }
  
  return(list(
    Z1 = Z1,
    Z2 = Z2,
    Z3 = Z3,
    X = X,
    Y = Y,
    cluster = cluster,
    
    T0 = T0,
    T1 = T1,
    W = W,
    nu = nu,
    K_planned = K_planned,
    K_obs = K_obs,
    V_star = V_star,
    V_tilde = V_tilde,
    V = V,
    r = r,
    rho = rho,
    
    CL = CL,
    CR = CR,
    DL = DL,
    DR = DR,
    Dobs = Dobs,
    Dnone = Dnone
  ))
}

#' datagen_negbin
#'
#' @param n 
#' @param p_z2 
#' @param mu_z3 
#' @param sd_z3 
#' @param comps 
#' @param w 
#' @param gamma 
#' @param x_sd 
#' @param beta_x 
#' @param beta_z 
#' @param tau 
#' @param beta_omega 
#' @param alpha_start 
#' @param alpha_width 
#' @param sigma_start 
#' @param sigma_width 
#' @param jitter_sd 
#' @param rho 
#' @param eps_obs 
#'
#' @returns
#' @export
#'
#' @examples
datagen_negbin <- function(n = 100,                                                    # sample size
                            p_z2 = 0.8,                                                 # observed covariate Z2 ~ bern(p_z2)
                            mu_z3 = 5, sd_z3 = 1,                                       # mean and sd for observed covariate Z3 ~ N(mu_z3, sd_z3)
                            comps = 3, w = c(0.3, 0.5, 0.2),                            # number of components for mixture X|Z and probabilities for each component
                            gamma = matrix(c(0, 3, 1, 
                                             0, 1, 5, 
                                             0, 2, 2), 
                                           nrow = 3, byrow = TRUE)*365.25,              # linear coeff. for X|Z ~ N(Z*gamma, x_sd)
                            x_sd = c(0.5, 1, 1)*365.25,
                            beta_x = 0.4/365.25,                                        # coefficients for outcome linear predictor on X
                            beta_z = c(0.2, 0.1, 0.2),                                  # coefficients for outcome linear predictor on Z
                            tau = 1,                                                    # outcome precision parameter
                            beta_omega = matrix(c(0, 0.5, 0.25,
                                                  0, 0.25, 0.10),
                                                nrow = 2, byrow = TRUE),                # coefficients for visit generation observation category
                            alpha_start = c(5, 2, -0.5)*365.25,                         # parameter for start of observation window T0 ~ N(alpha * Z, sigma_start)
                            alpha_width = c(log(6*365.25), 0.25, -0.05),                # parameter for width of observation window log(W) ~ N(alpha * Z, sigma_width)
                            sigma_start = 1, sigma_width = 0.15,
                            jitter_sd = 30,                                             # sd for random noise added to each visit
                            rho = c(0.85, 0.9, 0.99),                                   # vector of probabilities for observing a visit
                            eps_obs = 350                                               # tolerance window for defining uncensored observations in X
) {
  
  ## observed covariates
  Z1 <- rep(1, n)                                                               # intercept
  Z2 <- rbinom(n = n, size = 1, p_z2)                                           # Z2 ~ binomial(p_z2)
  Z3 <- rnorm(n, mean = mu_z3, sd = sd_z3)                                      # Z3 ~ normal(mu_z3, sd_z3)
  Z <- cbind(Z1, Z2, Z3)
  
  ## censored covariate
  cluster <- sample(1:comps, n, replace = TRUE, prob = w)                       # sample X|Z cluster assignment from weights w
  means <- rowSums(Z * gamma[cluster, 1:3])                                     # mean of X|Z determined by cluster assignment and gamma row
  X <- rnorm(n, mean = means, sd = x_sd[cluster])                               # X|Z ~ normal(mean[cluster], sd[cluster])
  
  ## outcome
  eta <- beta_x * X + Z %*% beta_z
  theta <- exp(eta)
  Y <- rnegbin(n, lambda = theta)                                        # Y|X,Z ~ binomial(n, logit(eta))
  
  ## probability of observing a visit
  omega_high <- as.numeric(Z %*% beta_omega[1, ])                               # linear pred if subjects have high observation probability
  omega_med <- as.numeric(Z %*% beta_omega[2, ])                                # linear pred if subjects have medium observation probability
  omega_low <- rep(0, n)                                                        # linear pred if subjects have low observation probability
  
  exp_omega <- cbind(exp(omega_high), exp(omega_med), exp(omega_low))           # exponentiate linear predictors
  pi_obs <- exp_omega / rowSums(exp_omega)                                      # calculate probabilities for observation probability assignment
  
  nu <- apply(pi_obs, 1, function(p) {                                          # sample observation categories with ^ probabilities per subject
    sample(1:3, size = 1, prob = p)
  })
  
  ## observation window
  T0 <- rnorm(n, mean = as.numeric(Z %*% alpha_start), sd = sigma_start)        # observation window start visit
  W <- rlnorm(n, meanlog = as.numeric(Z %*% alpha_width), sd = sigma_width)     # width of observation window
  T1 <- T0 + W
  
  ## storage
  CL <- CR <- numeric(n)
  DL <- DR <- Dobs <- Dnone <- integer(n)
  K_planned <- K_obs <- integer(n)
  V_star <- V_tilde <- V <- r <- vector('list', n)
  
  ## generate visit data
  for (i in seq_len(n)) {
    
    # initialize empty (zero-length) storage for visit vectors and observation indicator vector
    V_star_i <- V_tilde_i <- V_i <- numeric(0)       
    r_i <- integer(0)
    
    K_planned[i] <- floor(W[i] / 365.25)                                        # number of visits if none missed
    
    if (K_planned[i] > 0) {                         
      V_star_i <- T0[i] + 365.25 * seq_len(K_planned[i])                        # visits start as every year from observation window start
      
      V_tilde_i <- V_star_i + rnorm(K_planned[i], mean = 0, sd = jitter_sd)     # add noise to annual visit
      V_tilde_i <- sort(pmin(pmax(V_tilde_i, T0[i]), T1[i]))                    # sort visits in case ordering switched from added noise
      
      r_i <- rbinom(K_planned[i], size = 1, prob = rho[nu[i]])                  # sample observation indicator
      V_i <- V_tilde_i[r_i == 1]                                                # mask visits with observation idicator
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
    
    if (length(V_i) == 0) {                                                     # if subject has NO observed visits
      CL[i] <- -Inf
      CR[i] <- Inf
      Dnone[i] <- 1
    } else if (X[i] < V_i[1]) {                                                 # if X is left censored
      CL[i] <- -Inf
      CR[i] <- V_i[1]
      DL[i] <- 1
    } else if (X[i] > V_i[K_obs[i]]) {                                          # if X is right censored
      CL[i] <- V_i[K_obs[i]]
      CR[i] <- Inf
      DR[i] <- 1
    } else {                                                                    # X either interval or uncensored
      right_idx <- which(V_i >= X[i])[1]                                        # smallest visit greater than X
      
      if (X[i] == V_i[right_idx]) {                                             # if X falls exactly on a visit
        CL[i] <- CR[i] <- X[i]
        Dobs[i] <- 1
      } else {                                                                  # if X is interval censored
        left <- V_i[right_idx - 1]
        right <- V_i[right_idx]
        
        if ((right - left) <= eps_obs) {                                        # if X is interval censored with endpoints within tolerance, define as uncensored
          CL[i] <- CR[i] <- X[i]
          Dobs[i] <- 1
        } else {
          CL[i] <- left
          CR[i] <- right
        }
      }
    }
  }
  
  return(list(
    Z1 = Z1,
    Z2 = Z2,
    Z3 = Z3,
    X = X,
    Y = Y,
    cluster = cluster,
    
    T0 = T0,
    T1 = T1,
    W = W,
    nu = nu,
    K_planned = K_planned,
    K_obs = K_obs,
    V_star = V_star,
    V_tilde = V_tilde,
    V = V,
    r = r,
    rho = rho,
    
    CL = CL,
    CR = CR,
    DL = DL,
    DR = DR,
    Dobs = Dobs,
    Dnone = Dnone
  ))
}

#' gelc_censoring
#'
#' @description
#' Generates finite interval-censoring bounds using the censoring mechanism
#' used in the GELc simulation study. The initial inspection gap is Uniform
#' on (0, mu), and subsequent gap widths are absolute Normal draws with
#' mean mu and standard deviation 0.75 * mu.
#'
#' @param x true censored covariate values
#' @param mu mean inspection gap
#'
#' @returns list containing CL and CR
#' @export
gelc_censoring <- function(x, mu) {
  
  n <- length(x)
  
  if (mu <= 0) {
    stop('mu must be greater than 0')
  }
  
  gap_initial <- runif(n, min = 0, max = mu)
  
  gap_generator <- function(nn) {
    abs(rnorm(nn, mean = mu, sd = 0.75 * mu))
  }
  
  gaps <- gap_generator(n * ceiling(max(x)))
  gap_matrix <- matrix(gaps, nrow = n)
  
  while (any(rowSums(cbind(0, gap_initial, gap_matrix)) < x)) {
    gaps <- c(gaps, gap_generator(n * 5))
    gap_matrix <- matrix(gaps, nrow = n)
  }
  
  inspection_times <- t(
    apply(
      cbind(0, gap_initial, gap_matrix),
      1,
      cumsum
    )
  )
  
  left_possible <- (x - inspection_times > -.Machine$double.eps^0.5)
  right_possible <- (inspection_times - x > -.Machine$double.eps^0.5)
  
  CL <- apply(left_possible * inspection_times, 1, max)
  
  CR <- apply(
    right_possible * inspection_times,
    1,
    function(z) min(z[z > 0])
  )
  
  if (any(CL > x) || any(CR < x)) {
    stop('Generated censoring interval does not contain X')
  }
  
  return(list(
    CL = CL,
    CR = CR
  ))
}


#' datagen_gelc_normal
#'
#' @description
#' Generates data for the Gaussian P1/GELc comparison using the censored
#' covariate and censoring mechanism from the GELc simulation design.
#'
#' @param n number of subjects
#' @param mu mean censoring interval gap
#' @param beta_x outcome coefficient for X
#' @param beta_z outcome coefficient for Z2
#' @param beta_0 outcome intercept
#' @param tau outcome precision
#' @param force_observed whether to force one exactly observed X value
#'
#' @returns data frame
#' @export
datagen_gelc_normal <- function(n = 500,
                                mu = 3,
                                beta_x = 0.1,
                                beta_z = 0,
                                beta_0 = 2,
                                tau = 1,
                                force_observed = TRUE) {
  
  ## censored covariate from GELc simulation
  X <- rexp(n, rate = 1 / 12)
  
  ## additional fully observed covariate required for current P1 implementation
  Z1 <- rep(1, n)
  Z2 <- rnorm(n)
  
  ## Gaussian outcome
  eta <- beta_0 + beta_z * Z2 + beta_x * X
  
  Y <- rnorm(
    n,
    mean = eta,
    sd = 1 / sqrt(tau)
  )
  
  ## GELc censoring mechanism
  censoring <- gelc_censoring(
    x = X,
    mu = mu
  )
  
  CL <- censoring$CL
  CR <- censoring$CR
  
  Dobs <- rep(0L, n)
  DL <- rep(0L, n)
  DR <- rep(0L, n)
  
  ## temporary workaround for current P1 implementation
  if (force_observed) {
    CL[1] <- X[1]
    CR[1] <- X[1]
    Dobs[1] <- 1L
  }
  
  return(data.frame(
    Y = Y,
    X = X,
    Z1 = Z1,
    Z2 = Z2,
    CL = CL,
    CR = CR,
    Dobs = Dobs,
    DL = DL,
    DR = DR
  ))
}