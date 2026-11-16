# Gaussian kernel normalized on a line or on the plane, with fixed support.
lgcp_gaussian_kernel <- function(distance, width, cutoff, scaling) {
  relative_cutoff <- cutoff / width
  # Broad kernels need the central normal integral without CDF subtraction.
  # The integrated Taylor series has truncation error below machine precision here.
  # Log normalizers and scaled distances avoid squaring extreme widths.
  if (scaling == "line" && relative_cutoff < 0.01) {
    squared_cutoff <- relative_cutoff^2
    correction <- -squared_cutoff/6 + squared_cutoff^2/40 - squared_cutoff^3/336
    log_peak <- -log(2) - log(cutoff) - log1p(correction)
  } else if (scaling == "line") {
    central_mass <- stats::pnorm(relative_cutoff) - stats::pnorm(-relative_cutoff)
    log_peak <- -log(2*pi)/2 - log(width) - log(central_mass)
  } else if (relative_cutoff < sqrt(.Machine$double.eps)) {
    log_peak <- -log(pi) - 2*log(cutoff)
  } else {
    disk_mass <- -expm1(-0.5 * relative_cutoff^2)
    log_peak <- -log(2*pi) - 2*log(width) - log(disk_mass)
  }
  weight <- numeric(length(distance))
  inside <- distance < cutoff
  weight[inside] <- exp(log_peak - 0.5 * (distance[inside] / width)^2)
  weight
}

# The augmented village likelihood includes the activity state as well as
# the active-village survey. Masked activity signs contribute no information.
lgcp_exposure_loglik <- function(exposure, predictor, S) {
  value <- lgcp_survey(exposure, predictor, S$y, S$n, S$link)$loglik
  if (!is.null(S$activity_sign)) {
    observed <- S$activity_sign != 0L
    eta <- S$gate_base[observed] + S$gate_slope*exposure[observed]
    value <- value + sum(stats::plogis(S$activity_sign[observed]*eta, log.p=TRUE))
  }
  value
}

# A symmetric random walk in log width targets a uniform prior on log width.
# Fixed support preserves the kernel's sparsity pattern and its sum-tree layout.
lgcp_width_step <- function(width, counts, high, spread, exposure, predictor, S,
                            bounds, proposal_sd) {
  candidate <- exp(log(width) + stats::rnorm(1L, sd=proposal_sd))
  unchanged <- list(width=width, K=S$K, exposure=exposure, accepted=FALSE)
  if (!is.finite(candidate) || candidate < bounds[1] || candidate > bounds[2]) return(unchanged)
  K <- S$K
  K@x <- lgcp_gaussian_kernel(S$kernel_distance, candidate, S$kernel_cutoff, S$kernel_scaling)
  if (any(!is.finite(K@x))) return(unchanged)
  mass <- if (length(high)) (counts-high)*(1-spread)+high*(1+spread) else counts
  proposed_exposure <- as.numeric(K %*% mass)
  log_ratio <- lgcp_exposure_loglik(proposed_exposure, predictor, S) -
    lgcp_exposure_loglik(exposure, predictor, S)
  if (is.finite(log_ratio) && log(stats::runif(1L)) < log_ratio) {
    return(list(width=candidate, K=K, exposure=proposed_exposure, accepted=TRUE))
  }
  unchanged
}
