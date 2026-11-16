# Validate posterior site counts against direct integration on one river cell.
# Run with the current odispers installed. The two villages share every site.
library(odispers)
library(sf)

crs <- "+proj=aeqd +lat_0=0 +lon_0=30 +datum=WGS84 +units=km +no_defs"
river <- st_sf(HYRIV_ID = 1L, NEXT_DOWN = 0L, geometry = st_sfc(
  st_linestring(rbind(c(0, 0), c(2, 0))), crs = crs))
nodes <- river_nodes(river, spacing = 3, crs = crs)
graph <- river_graph(nodes)
survey <- st_as_sf(data.frame(x = c(1, 1.5), north = 0, y = c(2, 1), n = c(4, 3)),
                   coords = c("x", "north"), crs = crs)
phi <- 5
source_rate <- 0.8
sigma_u <- 0.55
range <- 2

# Only the midpoint field enters this likelihood. Its marginal is univariate
# Gaussian, so each count's prior mass needs just one numerical integral.
field_var <- as.numeric(graph$A %*% solve(as.matrix(wm_precision(graph, range, sigma_u))) %*%
                        Matrix::t(graph$A))
mu0 <- nodes$len * source_rate
kernel <- wendland(c(0, 0.5) / phi) / (2 / 3 * phi)
counts <- 0:80
normal_integral <- function(f) {
  integrate(function(z) dnorm(z) * f(mu0 * exp(sqrt(field_var) * z)),
            -10, 10, abs.tol = 1e-15, rel.tol = 1e-11)$value
}
count_mass <- vapply(counts, function(k) normal_integral(function(mu) dpois(k, mu)), 0)
prior_tail_bound <- normal_integral(function(mu) ppois(max(counts), mu, lower.tail = FALSE)) +
  2 * pnorm(-10)
village_likelihood <- vapply(counts, function(k) {
  prod(dbinom(survey$y, survey$n, -expm1(-k * kernel)))
}, 0)
posterior_mean <- function(mass) weighted.mean(counts, mass * village_likelihood)
exact <- posterior_mean(count_mass)

# Another exact check sums out the point process *before* integrating the
# Gaussian field. Expand the positive individuals' factors (1-exp(-E_i)).
# Each resulting joint survival term has a Poisson generating function.
# This example needs (y1+1)(y2+1) terms; this product grows across villages.
b <- as.matrix(expand.grid(lapply(survey$y, function(y) 0:y)))
mass_coeff <- (-1)^rowSums(b) * apply(b, 1, function(z) prod(choose(survey$y, z)))
survival_kernel <- as.numeric(sweep(b, 2, survey$n - survey$y, "+") %*% kernel)
evidence <- prod(choose(survey$n, survey$y)) * sum(mass_coeff * vapply(
  survival_kernel, function(h) normal_integral(function(mu) exp(mu * expm1(-h))), 0))
stopifnot(abs(evidence - sum(count_mass * village_likelihood)) < 1e-10)

# A second deterministic benchmark integrates the Gaussian conditional by
# ordinary Laplace at each count, then sums those masses against the same data.
laplace_mass <- vapply(counts, function(k) {
  u <- uniroot(function(u) k - mu0 * exp(u) - u / field_var, c(-30, 30), tol = 1e-12)$root
  mu <- mu0 * exp(u)
  dpois(k, mu) * exp(-u^2 / (2 * field_var)) / sqrt(1 + field_var * mu)
}, 0)

# Structured variational Bayes: fit a separate Gaussian q(u | k) to each
# possible count. Its optimal mean and variance satisfy m = V(k-t) and
# v = 1/(1/V+t), with t = mu0 exp(m+v/2) and V the prior field variance.
# The exponential of the optimised lower bound approximates P(N=k).
vb_mass <- vapply(counts, function(k) {
  residual <- function(log_t) {
    t <- exp(log_t)
    log_t - log(mu0) - field_var * (k - t) - 0.5 / (1 / field_var + t)
  }
  t <- exp(uniroot(residual, c(-30, 30), tol = 1e-12)$root)
  m <- field_var * (k - t)
  v <- 1 / (1 / field_var + t)
  kl <- ((m^2 + v) / field_var - 1 + log(field_var / v)) / 2
  exp(k * (log(mu0) + m) - t - lgamma(k + 1) - kl)
}, 0)
# Configurations with zero survey likelihood receive zero posterior mass.
count_summary <- function(mass) {
  w <- mass * village_likelihood
  avg <- posterior_mean(mass)
  c(mean = avg, sd = sqrt(weighted.mean((counts - avg)^2, w)),
    prob_one = w[counts == 1] / sum(w))
}
deterministic <- rbind(quadrature = count_summary(count_mass),
                       conditional_Laplace = count_summary(laplace_mass),
                       conditional_VB = count_summary(vb_mass))
print(deterministic, digits = 8)
print(c(prior_tail_bound = prior_tail_bound))
stopifnot(prior_tail_bound < 1e-9, all(vb_mass <= count_mass + 1e-12))

for (method in c("laplace", "slice")) {
  chains <- lapply(c(1106, 2106, 3106, 4106), function(seed) {
    set.seed(seed)
    fit <- odisp_fit_lgcp(y ~ 0, ~ 0, survey, nodes, graph, phi = phi, range = range,
                          sigma_u = sigma_u, source_rate = source_rate,
                          n_iter = 4000, burnin = 800, field_update = method)
    fit$draws[, "total_sites"]
  })
  # Use consecutive batches within each chain to account for autocorrelation.
  batches <- unlist(lapply(chains, function(x) tapply(x, rep(1:16, each = 200), mean)))
  draws <- unlist(chains)
  mcse <- sd(batches) / sqrt(length(batches))
  cat(method, "update\n")
  print(c(mean = mean(draws), sd = sd(draws), prob_one = mean(draws == 1), mcse = mcse))
  stopifnot(abs(mean(draws) - exact) < 4 * mcse)
}
