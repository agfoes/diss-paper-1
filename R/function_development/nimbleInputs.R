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