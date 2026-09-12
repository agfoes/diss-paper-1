#' BCCglm
#'
#' @param data 
#' @param outcome_formula 
#' @param covariate_formula 
#' @param censored_covariate 
#' @param censoring_bounds 
#' @param family 
#'
#' @returns
#' @export
#'
#' @examples
BCCglm <- function(
    data,
    outcome_formula,
    covariate_formula,
    censored_covariate,
    censoring_bounds,
    family <- "Gamma(link = 'log')"
    ) {
  
  # remove censored variable from outcome formula
  observed_formula <- update(
    outcome_formula,
    paste(". ~ . -", paste(censored_covariate, collapse = " - "))
  )
  
  # create model frame, response, matrices
  outcome_mf <- model.frame(observed_formula, data = data)
  Y <- model.response(outcome_mf)
  Zy <- model.matrix(observed_formula, data = data) # automatically adds an intercept
  Zx <- model.matrix(delete.response(terms(covariate_formula)),
                     data = data) # automatically adds an intercept
  
  # calculate dimensions of design matrices
  pzy <- dim(Zy)[2] # includes intercept + dim of supplied uncensored covariates
  py <- pzy + 1 # add one for the censored covariate
  pzx <- dim(Zx)[2]
  
  # returns list of nobs, ncen, Ndata, Nconstants, Ninits
  nimble_inputs <- nimbleInputs(CL = data$CL,
                                CR = data$CR,
                                Zy = Zy,
                                Zx = Zx,
                                Y = Y)
  
  # outcome distribution family and link function for model code
  family_def <- as.character(family)
  model_code <- buildModelCode(family = as.character(family_def),
                               nobs = nimble_inputs$nobs,
                               ncen = nimble_inputs$ncen)
  
  nimbleOptions(enableMacros = TRUE)
  model <- nimbleModel(
    code = model_code,
    data = nimble_inputs$Ndata,
    constants = nimble_inputs$Nconstants,
    inits = nimble_inputs$Ninits
  )
  
  conf <- configureMCMC(model, useConjugacy = TRUE)
  conf$removeSamplers("beta")
  
  sampler_type = "AF_slice"
  if (family_def == "gaussian(link = 'identity')") {
    conf$addSampler(
      target = "beta",
      type = sampler_beta_normal_conjugate_block,
      control = list(
        pzy = pzy,
        ncen = as.integer(nimble_inputs$ncen),
        nobs = as.integer(nimble_inputs$nobs)
      )
    )
  } else if (sampler_type == 'block') {
    
    conf$addSampler(
      target = "beta",
      type = 'RW_block',
      control = list(scale = 0.05, 
                     adaptInterval = 100)
    )
  } else if (sampler_type == 'AF_slice') {
    
    conf$addSampler(
      target = "beta",
      type = 'AF_slice'
    )
  }
  
  mcmc <- buildMCMC(conf)
  Cmodel <- compileNimble(model)
  Cmcmc <- compileNimble(mcmc, project = Cmodel)
  
  runtime <- system.time({
    samples <- runMCMC(
      Cmcmc,
      niter = 100,
      nburnin = 40,
      thin = 1,
      nchains = 1,
      samplesAsCodaMCMC = TRUE
    )
  })
  
  return(list(
    runtime = runtime,
    samples = samples
  ))
}

#' buildModelCode
#'
#' @param family 
#' @param nobs 
#' @param ncen 
#'
#' @returns
#' @export
#'
#' @examples
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

#' nimbleInputs
#'
#' @param obs_tolerance 
#' @param upper_bound 
#' @param lower_bound 
#' @param L 
#' @param CL 
#' @param CR 
#' @param Zy 
#' @param Zx 
#' @param Y 
#' @param sigma_bx 
#' @param sigma_bz 
#'
#' @returns
#' @export
#'
#' @examples
nimbleInputs <- function(
    obs_tolerance = 1e-8,
    upper_bound = 1e10,
    lower_bound = -1e10,
    L = 50,
    CL,
    CR,
    Zy,
    Zx,
    Y = Y,
    sigma_bx = 1,
    sigma_bz = 1
) {
  
  # create dimension values
  pzy <- dim(Zy)[2]
  py <- pzy + 1
  pzx <- dim(Zx)[2]
  
  # create censoring indicators, counts, and indices
  Dobs <- as.numeric(abs(CR - CL) <= obs_tolerance)
  DL <- as.numeric(CL == -Inf | CL == lower_bound)
  DR <- as.numeric(CR == Inf | CR == upper_bound)
  
  n <- dim(Zy)[1]
  nobs <- sum(Dobs)
  ncen <- n - nobs
  
  idx_obs <- array(which(Dobs == 1), dim = nobs)
  idx_cen <- array(which(Dobs == 0), dim = ncen)
  
  # adjust censoring bounds in case of infinite endpoints to provided lower/upper bounds
  CL[is.infinite(CL)] <- lower_bound
  CR[is.infinite(CR)] <- upper_bound
  CL_cen <- CL[idx_cen]
  CR_cen <- CR[idx_cen]
  
  # prior covariance for beta vector
  cov_beta <- diag(c(sigma_bx, sigma_bz), py)
  
  # initial values
  beta = rep(0, py)
  tau = 1
  lalpha = 0
  sigmasqTilde = rep(1, L)
  gammaTilde = matrix(rnorm(L*pzx, 0, 1), nrow = L, ncol = pzx)
  xi <- sample(1:L, n, replace = TRUE)
  v <- rbeta(L-1, 1, 1)
  
  if (ncen > 0) {
    x_cen <- rep(0, ncen)
    
    for (i in 1:ncen) {
      if (DL[idx_cen[i]] == 1) {
        x_cen[i] = CR_cen[i] - 0.1
      } else if (DR[idx_cen[i]] == 1) {
        x_cen[i] = CL_cen[i] + 0.1
      } else if (Dobs[idx_cen[i]] == 0) {
        x_cen[i] = 0.5 * (CL_cen[i] + CR_cen[i])
      }
    }
  }
  
  
  # separate observed and censored data
  x_obs <- CL[Dobs == 1]
  zy_obs <- Zy[Dobs == 1, , drop = FALSE]
  zy_cen <- Zy[Dobs == 0, , drop = FALSE]
  zx_obs <- Zx[Dobs == 1, , drop = FALSE]
  zx_cen <- Zx[Dobs == 0, , drop = FALSE]
  y_obs <- Y[Dobs == 1]
  y_cen <- Y[Dobs == 0]
  
  
  # nimble constants to be returned as list
  Nconstants <- list(L = L,
                     py = py,
                     pzy = pzy,
                     pzx = pzx,
                     beta_mean = rep(0, py),
                     beta_cov = cov_beta,
                     mu_gamma = 0,
                     idx_obs = idx_obs,
                     idx_cen = idx_cen
  )
  
  # nimble data to be returned as list
  Ndata <- list(y_obs = y_obs,
                y_cen = y_cen,
                zy_obs = zy_obs,
                zy_cen = zy_cen,
                zx_obs = zx_obs,
                zx_cen = zx_cen,
                CL_cen = CL_cen,
                CR_cen = CR_cen,
                x_obs = x_obs,
                constraint_data = rep(1, ncen)
  )
  
  # nimble initial values to be returned as list
  Ninits <- list(beta = beta,
                 tau = tau,
                 lalpha = lalpha,
                 sigmasqTilde = sigmasqTilde,
                 gammaTilde = gammaTilde,
                 xi = xi,
                 v = v,
                 x_cen = x_cen
  )
  
  # return list of data, constants, and initial values for nimble model
  return(list(ncen = ncen,
              nobs = nobs,
              Ndata = Ndata,
              Nconstants = Nconstants,
              Ninits = Ninits))
}

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