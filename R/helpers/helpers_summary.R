
# This file contains functions to create posterior summaries and append results to result tables.

#' posterior_summary
#' 
#' @description
#' This function summarizes MCMC samples and produces a tibble with key metrics
#' and numerics for saving.
#' 
#'
#' @param samples 
#' @param sampler_name 
#' @param truth 
#' @param key_vars 
#'
#' @returns
#' @export
#'
#' @examples
posterior_summary <- function(samples, sampler_name, truth, key_vars, runtime) {
  
  samples_matrix <- as.matrix(samples)
  
  beta_cols <- grep(
    '^beta\\[[0-9]+\\]$',
    colnames(samples_matrix)
  )
  
#  colnames(samples_matrix)[beta_cols] <- c(
#    'betax',
#    paste0('betaz[', seq_len(pz), ']')
#  )
  
  ## occupied clusters
  xi_cols <- grep("^xi\\[", colnames(samples_matrix))
  occupied <- apply(samples_matrix[, xi_cols, drop = FALSE],
                    1,
                    function(x) length(unique(x)))
  
  ## parameter summaries
  samples <- samples[, c(colnames(samples_matrix)[beta_cols], "tau"), drop = FALSE]
  samples_mcmc <- coda::as.mcmc(samples)
  
  parameter <- colnames(samples)
  truth_values <- unname(truth[parameter])
  
  post_mean <- colMeans(samples)
  post_sd <- apply(samples, 2, sd)
  post_median <- apply(samples, 2, median)
  
  ess <- as.numeric(coda::effectiveSize(samples_mcmc))
  hpd <- coda::HPDinterval(samples_mcmc, prob = 0.95)
  
  tibble::tibble(
    sampler = sampler_name,
    parameter = parameter,
    truth = truth_values,
    mean = post_mean,
    sd = post_sd,
    median = post_median,
    hpd_lower = hpd[, 'lower'],
    hpd_upper = hpd[, 'upper'],
    ess = ess,
    mcse = post_sd / sqrt(ess),
    runtime = runtime,
    ess_sec = ess / runtime,
    squared_error = (post_mean - truth_values)^2,
    covered = (truth_values >= hpd[, 'lower'] &
      truth_values <= hpd[, 'upper']),
    mean_occupied = mean(occupied),
    sd_occupied = sd(occupied),
    median_occupied = median(occupied),
    min_occupied = min(occupied),
    max_occupied = max(occupied)
  )
}

safe_ratio <- function(numerator, denominator) {
  if (denominator == 0) {
    NA
  } else {
    numerator / denominator
  }
}

safe_mean <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    NA
  } else {
    mean(x, na.rm = TRUE)
  }
}

safe_median <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    NA
  } else {
    median(x, na.rm = TRUE)
  }
}

safe_sd <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) < 2) {
    NA
  } else {
    sd(x)
  }
}

safe_min <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    NA
  } else {
    min(x, na.rm = TRUE)
  }
}

safe_max <- function(x) {
  if (length(x) == 0 || all(is.na(x))) {
    NA
  } else {
    max(x, na.rm = TRUE)
  }
}

data_censoring_summary <- function(data, sim_id) {
  DI <- 1 - data$Dobs - data$DL - data$DR
  censored <- 1 - data$Dobs
  
  if (any(abs(data$Dobs + data$DL + data$DR + DI - 1) > 1e-8)) {
    stop("Censoring indicators do not sum to one for sim_id = ", sim_id)
  }
  
  interval_width <- data$CR[DI == 1] - data$CL[DI == 1]
  
  data.frame(
    sim_id = sim_id,
    n = length(data$X),
    
    prop_observed = mean(data$Dobs),
    prop_censored = mean(censored),
    
    prop_left = mean(data$DL),
    prop_right = mean(data$DR),
    prop_interval = mean(DI),
    
    prop_left_given_censored = safe_ratio(sum(data$DL), sum(censored)),
    prop_right_given_censored = safe_ratio(sum(data$DR), sum(censored)),
    prop_interval_given_censored = safe_ratio(sum(DI), sum(censored)),
    
    n_interval = sum(DI),
    
    mean_interval_width = safe_mean(interval_width),
    median_interval_width = safe_median(interval_width),
    sd_interval_width = safe_sd(interval_width),
    min_interval_width = safe_min(interval_width),
    max_interval_width = safe_max(interval_width)
  )
}


