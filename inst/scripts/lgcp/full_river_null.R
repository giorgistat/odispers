# Independent null landscapes on the full river network. Conditional on a
# joint coefficient draw from the fit, redraw the graph field, every
# breeding site and all 900 survey outcomes. Replicate 1 is the example map
# and the preselected null dataset for a marked refit.
library(odispers)
library(Matrix)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 3:4) stop("supply inputs.rds, observed chain directory, output directory and optional marks")
marks <- if (length(args) == 4L) match.arg(args[4], c("unit", "two_point")) else "unit"
a <- readRDS(args[1])
out <- args[3]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
seeds <- c(8106, 8206, 8306, 8406)
parts <- lapply(seeds, function(seed) {
  fit <- readRDS(file.path(args[2], sprintf("%s_%d.rds", marks, seed)))
  hyper <- if (!is.null(fit$settings)) unlist(fit$settings[c("phi","range","sigma_u","source_rate")]) else
    c(phi=40,range=25,sigma_u=.5,source_rate=.04)
  # Fits saved before the logit option used cloglog.
  link <- if (is.null(fit$settings$link)) "cloglog" else fit$settings$link
  list(draws=fit$draws, hyper=hyper, link=link, occurrence=fit$settings$occurrence,
       occurrence_exposure=isTRUE(fit$settings$occurrence_exposure))
})
settings <- parts[[1]]$hyper
if (!all(vapply(parts,function(x) identical(settings,x$hyper),TRUE)))
  stop("chain hyperparameters differ")
link <- parts[[1]]$link
if (!link %in% c("logit", "cloglog") || !all(vapply(parts,function(x) identical(link,x$link),TRUE)))
  stop("chain survey links must agree")
occurrence <- parts[[1]]$occurrence
include_exposure <- parts[[1]]$occurrence_exposure
if (!all(vapply(parts,function(x) identical(x$occurrence_exposure,include_exposure),TRUE)))
  stop("chain occurrence exposure specifications must agree")
if (!all(vapply(parts,function(x) identical(deparse(x$occurrence),deparse(occurrence)),TRUE)))
  stop("chain occurrence formulas must agree")
params <- do.call(rbind,lapply(parts,"[[","draws"))

expected_theta <- if (length(a$source_covariates)) paste0("theta.",a$source_covariates) else character()
expected_beta <- c("beta.(Intercept)", if (length(a$survey_covariates)) paste0("beta.",a$survey_covariates))
if (!setequal(grep("^theta[.]",colnames(params),value=TRUE), expected_theta) ||
    !setequal(grep("^beta[.]",colnames(params),value=TRUE), expected_beta))
  stop("simulation design does not match fitted coefficients; supply the same inputs used to fit the model")

# The source count in a fitted draw describes the old landscape. Only its
# coefficients and mark spread generate the new landscape.
params <- params[, colnames(params) != "total_sites", drop = FALSE]

phi <- as.numeric(settings["phi"])
source_rate <- as.numeric(settings["source_rate"])
range <- as.numeric(settings["range"])
sigma_u <- as.numeric(settings["sigma_u"])
if (!is.null(a$phi) && a$phi != phi) stop("input radius differs from fitted radius")
p <- a$pairs
v <- a$villages
if (!is.null(occurrence)) {
  design <- model.frame(occurrence,sf::st_drop_geometry(v),na.action=na.fail)
  G <- model.matrix(occurrence,design)
  if (include_exposure) G <- cbind(G,exposure=numeric(nrow(v)))
  gate_names <- if (ncol(G)) paste0("occurrence.",colnames(G)) else character()
  gate_offset <- model.offset(design)
  if (is.null(gate_offset)) gate_offset <- numeric(nrow(v))
  stopifnot(setequal(gate_names,grep("^occurrence[.]",colnames(params),value=TRUE)))
}
if (p$max_range < phi) stop("survey pairs do not cover the full fitted radius")
inside <- p$r < phi
K <- sparseMatrix(i = p$i[inside], j = p$j[inside], x = wendland(p$r[inside] / phi) / (2 / 3 * phi),
                  dims = c(nrow(v), nrow(a$nodes)))
used <- sort(unique(p$j[inside]))
prior_ch <- Cholesky(wm_precision(a$graph, range, sigma_u), LDL = FALSE)

xy <- sf::st_coordinates(sf::st_transform(v, attr(a$nodes, "crs")))
dist <- as.matrix(stats::dist(xy))
near <- which(upper.tri(dist) & dist <= 40, arr.ind = TRUE)
summaries <- function(y) {
  prop <- y / v$n
  c(positive_fraction = sum(y) / sum(v$n), zero_surveys = sum(y == 0L),
    at_least_half_positive = sum(prop >= 0.5),
    nearby_absolute_difference = mean(abs(prop[near[, 1]] - prop[near[, 2]])))
}

simulation_seed <- if (marks == "unit") 621106 else 621107
set.seed(simulation_seed)
replicates <- 250L
sim_y <- matrix(0L, replicates, nrow(v))
sim_p <- matrix(0, replicates, nrow(v))
draw_ids <- integer(replicates)
if (length(a$source_covariates)) {
  Z <- as.matrix(as.data.frame(a$nodes)[, a$source_covariates, drop=FALSE])
  stopifnot(all(paste0("theta.", a$source_covariates) %in% colnames(params)))
}
if (length(a$survey_covariates)) {
  X <- as.matrix(sf::st_drop_geometry(v)[, a$survey_covariates, drop=FALSE])
  stopifnot(all(paste0("beta.", a$survey_covariates) %in% colnames(params)))
}
totals <- matrix(0, replicates, 3, dimnames = list(NULL, c("all_sites", "within_reach", "beta")))
for (b in seq_len(replicates)) {
  draw_ids[b] <- sample.int(nrow(params), 1L)
  parameters <- params[draw_ids[b], ]
  beta <- unname(parameters["beta.(Intercept)"])
  transmission <- rep(beta,nrow(v))
  if (length(a$survey_covariates))
    transmission <- transmission + as.numeric(X %*% parameters[paste0("beta.",a$survey_covariates)])
  d <- if (marks == "two_point") unname(parameters["mark_spread"]) else 0
  # The current precision is factored as P Q P' = L L'.
  z <- rnorm(a$graph$nv)
  u <- as.numeric(solve(prior_ch, solve(prior_ch, z, system = "Lt"), system = "Pt"))
  log_rate <- log(source_rate) + as.numeric(a$graph$A %*% u)
  if (length(a$source_covariates))
    log_rate <- log_rate + as.numeric(Z %*% parameters[paste0("theta.",a$source_covariates)])
  N <- rpois(nrow(a$nodes), a$nodes$len * exp(log_rate))
  H <- if (marks == "two_point") rbinom(length(N), N, 0.5) else integer(length(N))
  mass <- if (marks == "two_point") (N-H)*(1-d) + H*(1+d) else N
  E <- as.numeric(K %*% mass)
  log_exposure <- log(E) + transmission
  sim_p[b, ] <- if (link == "logit") plogis(log_exposure) else -expm1(-exp(log_exposure))
  active <- rep(1L,nrow(v))
  if (!is.null(occurrence)) {
    if (include_exposure) G[,"exposure"] <- E
    q <- plogis(gate_offset + as.numeric(G %*% parameters[gate_names]))
    active <- rbinom(nrow(v),1L,q)
    sim_p[b, ] <- sim_p[b, ] * active
  }
  sim_y[b, ] <- rbinom(nrow(v), v$n, sim_p[b, ])
  totals[b, ] <- c(sum(N), sum(N[used]), beta)
  if (b == 1L) {
    saveRDS(list(N = N, high = H, spread = d, u = u, rate_per_km = exp(log_rate), exposure = E,
                 probability = sim_p[b, ], beta = beta, transmission = transmission, parameters = parameters,
                 link = link, active = active), file.path(out, "first_landscape.rds"))
    saveRDS(sim_y[b, ], file.path(out, "first_y.rds"))
  }
}
checks <- t(apply(sim_y, 1L, summaries))
envelope <- data.frame(statistic = colnames(checks), observed = as.numeric(summaries(v$y)),
                       lo = apply(checks, 2, quantile, 0.025), median = apply(checks, 2, median),
                       hi = apply(checks, 2, quantile, 0.975), row.names = NULL)
saveRDS(list(y = sim_y, probability = sim_p, totals = totals, checks = checks,
             seeds = seeds, simulation_seed = simulation_seed, marks = marks, link = link,
             generating_draws = params[draw_ids,,drop=FALSE]),
        file.path(out, if (marks == "unit") "null_simulations.rds" else "marked_simulations.rds"))
write.csv(envelope, file.path(out, if (marks == "unit") "null_envelope.csv" else "marked_envelope.csv"), row.names = FALSE)
print(envelope)
cat(replicates, "replicates on", nrow(a$nodes), "cells and", nrow(v), "villages\n")

# Repeat surveys on the inferred landscape. Selecting one joint MCMC draw
# for all villages preserves shared sites and contrasts this check with the
# independent landscape draws above.
set.seed(simulation_seed + 100)
rep_y <- matrix(0L, replicates, nrow(v))
for (k in seq_along(seeds)) {
  fit <- readRDS(file.path(args[2], sprintf("%s_%d.rds", marks, seeds[k])))
  ids <- seq.int(k, replicates, by = length(seeds))
  rows <- sample.int(nrow(fit$prevalence_draws), length(ids), replace = TRUE)
  for (j in seq_along(ids))
    rep_y[ids[j], ] <- rbinom(nrow(v), v$n, fit$prevalence_draws[rows[j], ])
}
conditional <- t(apply(rep_y, 1, summaries))
conditional_envelope <- envelope
conditional_envelope[, c("lo","median","hi")] <- t(apply(conditional, 2, quantile, c(.025,.5,.975)))
write.csv(conditional_envelope, file.path(out,"conditional_envelope.csv"), row.names = FALSE)
saveRDS(rep_y, file.path(out,"repeat_survey_counts.rds"))
print(conditional_envelope)
