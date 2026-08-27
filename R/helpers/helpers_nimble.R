
#' nimble_constants
#'
#' @description
#' This function generates a list of constants (provided as function inputs) 
#' to be used in nimble model implementation.
#' 
#' @param data dataframe output from data generation function
#' @param L numeric truncation for DPMM components
#' @param pz numeric dimension of fully observed covariates
#' @param sigma_bx numeric variance of prior on beta_x ~ norm(0, sigma_bx)
#' @param sigma_bz numeric variance of prior on beta_z ~ norm(0, sigma_bz)
#' @param mu_gamma numeric mean of prior on gamma vector (X|Z ~ N(gamma*Z, s2tilde)) 
#'
#' @returns list of constants for nimble model implementation
#' @export
#'
#' @examples
nimble_constants <- function(data,
                             L = 5,
                             pz = pz,
                             sigma_bx = 1,
                             sigma_bz = 1,
                             mu_gamma = 0) {
  
  cov_beta = diag(c(sigma_bz, sigma_bz), (pz+1))
  
  idx_obs = which(data$Dobs == 1)
  idx_cen = which(data$Dobs == 0)
  
  # pad idx_cen/obs if they only have one entry - nimble doesn't like length-one vectors and it doesn't matter bc idx are always indexed, so the added entries aren't ever accessed
  if (length(idx_obs) < 2) {
    idx_obs <- c(idx_obs, rep(1, 2 - length(idx_obs)))
  }
  if (length(idx_cen) == 1) {
    idx_cen <- c(idx_cen, rep(1, 2 - length(idx_cen)))
  }
  
  return(list(
    L = L,
    p = pz,
    beta_mean = rep(0, (pz+1)), # added to allow for joint beta prior and non-zero centered beta prior
    #sigma_bx = sigma_bx,
    #sigma_bz = sigma_bz,
    beta_cov = cov_beta, # added to allow for joint beta prior
    mu_gamma = mu_gamma,
    nobs = sum(data$Dobs),
    ncen = length(data$Y) - sum(data$Dobs),
    idx_obs = idx_obs,
    idx_cen = idx_cen
  ))
}

#' nimble_data_years
#'
#' @description
#' This function generates data for use in nimble model implementation. 
#' In particular, this function generates a list with the censored covariate 
#' data transformed into the year scale instead of the day scale (as output 
#' from the data generation function).
#' 
#' @param data dataframe output from data generation function
#' @param constants list output from nimble_constants function
#'
#' @returns list of data with covariates in year-level for nimble model
#' implementation
#' @export
#'
#' @examples
nimble_data_years <- function(data, constants) {
  x_obs <- data$X[data$Dobs == 1]/365.25
  z = as.matrix(cbind(data$Z1, data$Z2, data$Z3))
  
  CL = data$CL[data$Dobs == 0]/365.25
  CL[is.infinite(CL)] <- min(data$X)/365.25 - 10/365.25
  CR = data$CR[data$Dobs == 0]/365.25
  CR[is.infinite(CR)] <- max(data$X)/365.25 + 10/365.25
  
  return(list(
    y_obs = data$Y[data$Dobs == 1],
    y_cen = data$Y[data$Dobs == 0],
    z_obs = z[data$Dobs == 1, , drop = FALSE],
    z_cen = z[data$Dobs == 0, , drop = FALSE],
    CL = CL,
    CR = CR,
    constraint_data = rep(1, constants$ncen),
    x_obs = x_obs
  ))
}

#' nimble_inits
#' 
#' @description
#' This function generates a list of initial values for nimble model
#' implementation. In particular, the censored covariate is initialized based
#' on whether the missing covariate is left, right, or interval censored.
#' 
#'
#' @param data dataframe output from data generation function
#' @param pz numeric dimension of fully observed covariates
#' @param constants list output from nimble_constants function
#' @param Ndata list output from nimble_data or nimble_data_years function
#'
#' @returns
#' @export
#'
#' @examples
nimble_inits <- function(data, pz, constants, Ndata) {
  
  L <- constants[["L"]]

  Ninits <- list(
    #betax = 0,
    #betaz = rep(0, pz),
    beta = rep(0, (pz+1)),
    tau = 1,
    
    lalpha = 0,
    
    sigmasqTilde = rep(1, L),
    gammaTilde = matrix(rnorm(L*pz, 0, 1), nrow = L, ncol = pz),
    
    xi = sample(1:L, length(data$Y), replace = TRUE),
    
    #x_cen = x_cen,
    
    v = rbeta(L - 1, 1, 1)
  )
  
  if (constants$ncen > 0) {
    idx_cen <- constants[["idx_cen"]]
    
    x_cen <- rep(0, constants$ncen)
    for (i in seq_len(constants$ncen)) {
      if (data$DL[idx_cen[i]] == 1) {
        x_cen[i] = Ndata$CR[i] - 0.1
      } else if (data$DR[idx_cen[i]] == 1) {
        x_cen[i] = Ndata$CL[i] + 0.1
      } else if (data$Dobs[idx_cen[i]] == 0) {
        x_cen[i] = 0.5*(Ndata$CL[i] + Ndata$CR[i])
      }
    }
    
    Ninits$x_cen <- x_cen
  }
  
  return(Ninits)
}

