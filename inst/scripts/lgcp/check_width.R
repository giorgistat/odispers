# Validate the conditional Gaussian-width update against direct quadrature.
# The two cells carry shared, marked sources. The second village is inactive,
# so its activity likelihood still depends on exposure although n = y = 0.
# Run with the current package installed: Rscript inst/scripts/lgcp/check_width.R
library(odispers)
library(Matrix)
library(posterior)

bounds <- c(.3, 8)
cutoff <- 3
iterations <- 8000L
burnin <- 2000L
initial_widths <- c(.4, 1.5, 6)
seeds <- 711070L + seq_along(initial_widths)
distances <- c(.2, 1.5, 2.5, .4)
counts <- c(2L, 1L)
high <- c(1L, 1L)
spread <- .7
mass <- c(2, 1.7)
predictor <- c(-.3, .1)
gate_base <- c(-.2, .8)
gate_slope <- 3

kernel <- odispers:::lgcp_gaussian_kernel(distances, 1, cutoff, "line")
S <- list(K = Matrix(matrix(kernel, 2), sparse = TRUE),
          kernel_distance = distances, kernel_cutoff = cutoff, kernel_scaling = "line",
          y = c(1, 0), n = c(4, 0), link = "logit", gate_base = gate_base,
          gate_slope = gate_slope, activity_sign = c(1L, -1L))

# Independent density in log width: the prior is constant in this coordinate.
log_density <- function(log_width) vapply(log_width, function(value) {
  width <- exp(value)
  K <- matrix(dnorm(distances, sd = width) /
                (pnorm(cutoff/width) - pnorm(-cutoff/width)), 2)
  exposure <- as.numeric(K %*% mass)
  odds <- exposure * exp(predictor)
  probability <- odds / (1 + odds)
  activity <- gate_base + gate_slope * exposure
  dbinom(1, 4, probability[1], log = TRUE) + plogis(activity[1], log.p = TRUE) +
    plogis(activity[2], lower.tail = FALSE, log.p = TRUE)
}, numeric(1))
center <- max(log_density(seq(log(bounds[1]), log(bounds[2]), length.out = 101)))
density <- function(x) exp(log_density(x) - center)
posterior_integral <- function(f) {
  integrate(function(x) f(x)*density(x), log(bounds[1]), log(bounds[2]),
            rel.tol = 1e-10)$value
}
normalizer <- posterior_integral(function(x) rep(1, length(x)))
exact <- c(log_width = posterior_integral(identity)/normalizer,
           width = posterior_integral(exp)/normalizer)

# Independent chains check the same two moments, allowing for Monte Carlo error.
draws <- matrix(NA_real_, iterations - burnin, length(seeds))
for (chain in seq_along(seeds)) {
  set.seed(seeds[chain])
  width <- initial_widths[chain]
  S$K@x <- odispers:::lgcp_gaussian_kernel(distances, width, cutoff, "line")
  exposure <- as.numeric(S$K %*% mass)
  for (iteration in seq_len(iterations)) {
    step <- odispers:::lgcp_width_step(width, counts, high, spread, exposure,
                                      predictor, S, bounds, .5)
    width <- step$width
    S$K <- step$K
    exposure <- step$exposure
    if (iteration > burnin) draws[iteration - burnin, chain] <- width
  }
}
checks <- do.call(rbind, lapply(names(exact), function(statistic) {
  values <- if (statistic == "log_width") log(draws) else draws
  mcse <- mcse_mean(values)
  data.frame(statistic = statistic, quadrature = exact[[statistic]],
             mcmc = mean(values), mcse = mcse,
             standardized_difference = abs(mean(values) - exact[[statistic]])/mcse,
             rhat = rhat(values), ess = ess_bulk(values))
}))
print(checks)
stopifnot(all(checks$standardized_difference < 5), all(checks$rhat < 1.02))
