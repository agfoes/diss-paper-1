DYNAMICGLM <- buildMacro(
  function(stoch, 
           LHS, 
           eta, 
           family, 
           tau = NULL, 
           modelInfo, 
           .env
  ) {
    
    family_obj <- eval(family, envir = .env)
    family_name <- family_obj$family
    link_name <- family_obj$link
    
    eta_index <- eta[[3]] # grabs the index idx_cen[k] from eta[idx_cen[k]]
    theta_node <- substitute(
      theta[INDEX],
      list(INDEX = eta_index)
    )
    
    # use family$link to determine theta
    if (link_name == "identity") {
      link_code <- substitute(
        THETA <- ETA,
        list(
          THETA = theta_node,
          ETA = eta
        )
      )
    } else if (link_name == "log") {
      link_code <- substitute(
        log(THETA) <- ETA,
        list(
          THETA = theta_node,
          ETA = eta
        )
      )
    } else if (link_name == "logit") {
      link_code <- substitute(
        logit(THETA) <- ETA,
        list(
          THETA = theta_node,
          ETA = eta
        )
      )
    } else {
      stop("Unsupported link function.")
    }
    
    # use family$family to determine outcome model
    if (family_name == "gaussian") {
      family_code <- substitute(
        Y ~ dnorm(mean = THETA,
                  tau = TAU),
        list(
          Y = LHS,
          THETA = theta_node,
          TAU = tau
        )
      )
    } else if (family_name == "poisson") {
      family_code <- substitute(
        Y ~ dpois(lambda = THETA),
        list(
          Y = LHS,
          THETA = theta_node
        )
      )
    } else if (family_name == "Gamma") {
      family_code <- substitute(
        Y ~ dgamma(shape = TAU,
                   rate = TAU / THETA),
        list(
          Y = LHS,
          THETA = theta_node,
          TAU = tau
        )
      )
    } else if (family_name == "binomial") {
      family_code <- substitute(
        Y ~ dbern(prob = THETA),
        list(
          Y = LHS,
          THETA = theta_node
        )
      )
    } else {
      stop("Unsupported family distribution.")
    }
    
    # combine family and link into final model code
    new_code <- as.call(c(list(as.name("{")),
                          list(link_code),
                          list(family_code)))
    
    list(code = new_code,
         modelInfo = modelInfo)
  },
  use3pieces = TRUE,
  unpackArgs = TRUE
)