library(rjags)
library(coda)


JWglm <- function(
    data,
    outcome_formula,
    covariate_formula,
    censored_covariate,
    censoring_bounds = c('CL', 'CR'),
    lower_inf = -1e4,
    upper_inf = 1e4,
    niter = 50000,
    burnin = 5000,
    thin = 1
) {
  
  ## outcome
  y_name <- as.character(outcome_formula[[2]])
  y <- data[[y_name]]
  
  ## outcome model matrix
  outcome_terms <- attr(
    terms(outcome_formula),
    'term.labels'
  )
  
  outcome_terms <- outcome_terms[outcome_terms != censored_covariate]
  outcome_intercept <- attr(terms(outcome_formula), 'intercept')
  
  if (length(outcome_terms) == 0) {
    rhs <- if (outcome_intercept == 1) '1' else '0'
  } else {
    rhs <- paste(outcome_terms, collapse = ' + ')
    if (outcome_intercept == 0) {
      rhs <- paste('0 +', rhs)
    }
  }
  
  Zy_formula <- as.formula(paste('~', rhs))
  
  Zy <- model.matrix(Zy_formula, data = data)
  
  
  ## covariate model matrix
  Zx <- model.matrix(
    delete.response(terms(covariate_formula)),
    data = data
  )
  
  
  ## censoring bounds
  CL <- data[[censoring_bounds[1]]]
  CR <- data[[censoring_bounds[2]]]
  
  cuts <- cbind(CL, CR)
  cuts[is.infinite(cuts[, 1]), 1] <- lower_inf
  cuts[is.infinite(cuts[, 2]), 2] <- upper_inf
  
  x_init <- numeric(nrow(data))
  
  for (i in seq_len(nrow(data))) {
    lower <- cuts[i, 1]
    upper <- cust[i, 2]
    
    if (upper == upper_inf) {
      x_init[i] <- lower + 1
    } else if (lower == lower_inf) {
      x_init[i] <- upper - 1
    } else {
      x_init <- max(mean(c(lower, upper)), 1e-4)
    }
  }
  
  
  ## dimensions
  n <- nrow(data)
  py <- ncol(Zy)
  px <- ncol(Zx)
  
  
  ## JAGS data
  jags_data <- list(
    n = n,
    py = py,
    px = px,
    y = y,
    Zy = Zy,
    Zx = Zx,
    cuts = cuts
  )
  
  
  ## initial values
  jags_inits <- list(
    x_cen = x_init
  )
  
  
  ## model
  model <- jags.model(
    textConnection(jw_model),
    data = jags_data,
    inits = jags_inits,
    n.chains = 1
  )
  
  update(
    model,
    n.iter = burnin
  )
  
  samples <- coda.samples(
    model,
    variable.names = c(
      'beta',
      'beta_x',
      'gamma',
      'shape_y'
    ),
    n.iter = niter,
    thin = thin
  )
  
  return(samples)
}

jw_model <- "
data {
    for (i in 1:n) {
        d[i] <- 1
    }
}
model {
    for (i in 1:n) {
        ## Gamma outcome
        y[i] ~ dgamma(
            shape_y,
            shape_y / mu_y[i]
        )
        log(mu_y[i]) <- inprod(Zy[i, 1:py], beta[1:py]) + beta_x * x_cen[i]
        
        ## interval-censored covariate
        d[i] ~ dinterval(x_cen[i], cuts[i, 1:2])
        
        ## Exponential covariate model
        x_cen[i] ~ dexp(rate_x[i])
        log(rate_x[i]) <- inprod(Zx[i, 1:px], gamma[1:px])
    }
    ## outcome priors
    for (j in 1:py) {
        beta[j] ~ dnorm(0, 0.001)
    }
    beta_x ~ dnorm(0, 0.001)
    shape_y ~ dgamma(0.001, 0.001)
    ## covariate model priors
    for (j in 1:px) {
        gamma[j] ~ dnorm(0, 0.001)
    }
}
"