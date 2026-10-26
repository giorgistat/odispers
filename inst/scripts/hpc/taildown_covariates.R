###############################################################################
# Zero-inflated river-network model with a tail-down breeding field and
# RiverATLAS habitat covariates in the breeding intensity: all steps, from the
# data shipped with odispers to the fits, confidence intervals,
# cross-validation and prediction maps. Self-contained: it needs only the
# package (install the GitHub version, remotes::install_github("giorgistat/odispers"))
# and does not use the output of any other script.
#
# Model, for village i with n_i examined and y_i with nodules:
#   P(Y_i = 0) = 1 - pi_i + pi_i (1 - p_i)^n_i,   p_i = 1 - exp(-lambda_i)
#   logit pi_i     = c0 + c1 log AI_i + a_E  log E_i          (zone of transmission)
#   log lambda_i   = b0 + b1 log AI_i + gamma log E_i + Z_i   (force of infection)
#   E_i = sum_j K(||x_i - s_j|| / 40 km) len_j w_j            (exposure, Wendland kernel)
#   log w_j = theta log Q_j + beta' h_j + u(s_j)              (breeding intensity)
# u: tail-down Gaussian field on the river network (exponential covariance in
# river distance, range estimated), Z_i ~ N(0, sigma_z^2).
#
# Models (nested, breeding covariates h):
#   M0  log discharge only (the selected tail-down model without covariates)
#   M1  + log(1 + channel slope, dm/km) + min{log(1 + upstream inundation %), 1.5}
#   M2  + hard rock (reach on basement or volcanic rock, GLiM classes)
#
# Steps for each model:
#   1. fit: Laplace approximation of the marginal likelihood; sigma_u,
#      sigma_z, a_E, gamma and the range of u estimated jointly
#   2. confidence intervals: Wald intervals for the regression coefficients
#      (from the Laplace approximation, conditional on the other
#      parameters); for sigma_u, sigma_z, a_E, gamma and the range, Wald
#      intervals from the observed information of the approximate marginal
#      likelihood (central finite differences; log scale for standard
#      deviations and range)
#   3. 10-fold random cross-validation (folds as in the earlier analyses,
#      set.seed(2026)), with the parameters other than the latent effects at
#      the full-data estimates: predictive distribution of each held-out
#      count, log score, AUC and Brier score for any positive, expected
#      zeros, non-randomised PIT (mean distribution function and its largest
#      deviation from uniform, with a simulated p-value)
#   4. predictions on a 20 km grid over Sudan and South Sudan
# Then a comparison table of all fitted models.
#
# Usage:
#   Rscript taildown_covariates.R [out_dir] [models]
#     out_dir  folder for the results (default taildown_covariates_results)
#     models   comma-separated subset of M0,M1,M2 (default all, in order)
# Each model warm-starts from the largest nested model already fitted in
# out_dir (M2 from M1, else M0), so running M0 first is fastest. Fits are
# checkpointed: a job that is stopped resumes from its best point.
# Example SLURM job (one process; about 1-2 h per model on one core):
#   #!/bin/bash
#   #SBATCH --time=12:00:00 --mem=8G --cpus-per-task=1
#   module load R
#   Rscript taildown_covariates.R results M0,M1,M2
###############################################################################
suppressPackageStartupMessages({ library(odispers); library(sf); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "taildown_covariates_results"
run_models <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("M0", "M1", "M2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
phi <- 40           # km, dispersal range of the kernel
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

# ---- 0. data --------------------------------------------------------------------------
data(remo); data(sudan_rivers); data(sudan_aridity); data(sudan_borders); data(sudan_river_covariates)
v <- subset(remo, country == "SUDAN")
v$y <- round(v$y)                                   # a few counts in the file are not integers
v <- st_as_sf(v, coords = c("long", "lat"), crs = 4326, remove = FALSE)
v$AI <- extract_aridity(v, sudan_aridity)
stopifnot(all(is.finite(v$AI)), all(v$AI > 0))
rivers <- sudan_rivers[sudan_rivers$DIS_AV_CMS >= 0.1, ]
nodes <- river_nodes(rivers, spacing = 2)           # quadrature nodes every <= 2 km
rc <- as.data.frame(sudan_river_covariates)[match(nodes$HYRIV_ID, sudan_river_covariates$HYRIV_ID), ]
stopifnot(!anyNA(rc$HYRIV_ID))
hard_classes <- c(2, 4, 7, 8, 9, 10, 12, 13)        # GLiM: basic/acid/intermediate volcanic and plutonic, metamorphic, pyroclastic
nodes$log_slope  <- log1p(rc$sgr_dk_rav)            # channel slope, dm per km
nodes$inundation <- pmin(log1p(rc$inu_pc_ult), 1.5) # upstream long-term maximum inundation %, capped
nodes$hard_rock  <- as.numeric(rc$lit_cl_cmj %in% hard_classes)
pairs <- river_pairs(v, nodes, phi)
graph <- directed_graph(river_graph(nodes, vertex_spacing = 5), "taildown")
log_msg("data: %d villages (%d with positives), %d nodes, %d graph vertices", nrow(v), sum(v$y > 0),
        nrow(nodes), graph$nv)

breeding <- list(M0 = ~ log(DIS_AV_CMS),
                 M1 = ~ log(DIS_AV_CMS) + log_slope + inundation,
                 M2 = ~ log(DIS_AV_CMS) + log_slope + inundation + hard_rock)
parent <- list(M0 = character(0), M1 = "M0", M2 = c("M1", "M0"))
res_file <- function(m) file.path(out_dir, sprintf("fit_%s.rds", m))

# ---- helpers: cross-validation (predictive distributions) -----------------------------
subset_setup <- function(S, idx) {
  S$y <- S$y[idx]; S$n <- S$n[idx]; S$X <- S$X[idx, , drop = FALSE]; S$K <- S$K[idx, , drop = FALSE]
  S$xy <- S$xy[idx, , drop = FALSE]; S$has <- S$has[idx]; S$lchoose <- sum(lchoose(S$n, S$y)); S$lai <- S$lai[idx]
  S
}
gh <- local({ m <- 20; J <- matrix(0, m, m)
  for (k in 1:(m - 1)) J[k, k + 1] <- J[k + 1, k] <- sqrt(k / 2)
  e <- eigen(J, symmetric = TRUE); list(x = sqrt(2) * e$values, w = e$vectors[1, ]^2) })
g2 <- expand.grid(i = seq_along(gh$x), j = seq_along(gh$x))
g2x1 <- gh$x[g2$i]; g2x2 <- gh$x[g2$j]; g2w <- gh$w[g2$i] * gh$w[g2$j]
pmf_dE <- function(md, me, Sig, n) {      # predictive pmf, 20 x 20 Gauss-Hermite over (logit pi, log lambda)
  L <- t(chol(Sig + diag(1e-10, 2)))
  d <- md + L[1, 1] * g2x1; e <- me + L[2, 1] * g2x1 + L[2, 2] * g2x2
  pi_ <- stats::plogis(d); p <- -expm1(-exp(pmin(e, 30)))
  pm <- sapply(0:n, function(y) sum(g2w * pi_ * stats::dbinom(y, n, p)))
  pm[1] <- pm[1] + sum(g2w * (1 - pi_))
  pm
}
cv_zil <- function(fit, folds) {
  S_all <- fit$setup; fld <- fit$fld
  p <- ncol(S_all$X); q <- ncol(S_all$B)
  one <- function(train, test) {
    S_tr <- subset_setup(S_all, train); S_te <- subset_setup(S_all, test)
    x0 <- list(beta = fit$coefficients[seq_len(p)], c = fit$coefficients[p + seq_len(p)],
               theta = fit$coefficients[2 * p + seq_len(q)], u = fit$u, z = rep(0, length(train)))
    ft <- odispers:::odisp_zil_laplace(fld, S_tr, x0)
    ev <- odispers:::odisp_zil_eval(ft$x, S_tr, fld)
    lk <- odispers:::odisp_zil_lik(ev$d, ev$e, S_tr)
    ops <- odispers:::odisp_zil_ops(ft$x, ev, lk, S_tr, fld)
    xt <- ft$x; xt$z <- rep(0, length(test))
    et <- odispers:::odisp_zil_eval(xt, S_te, fld)
    Jt <- as.matrix(et$V %*% S_te$B); Ju <- as.matrix(et$V %*% S_te$A); nz <- length(train)
    lapply(seq_along(test), function(k) {
      jd <- c(rep(0, p), S_te$X[k, ], fld$aE * Jt[k, ], fld$aE * Ju[k, ], rep(0, nz))
      je <- c(S_te$X[k, ], rep(0, p), fld$gamma * Jt[k, ], fld$gamma * Ju[k, ], rep(0, nz))
      hd <- ops$HGNinv(jd); he <- ops$HGNinv(je)
      Sig <- matrix(c(sum(jd * hd), sum(jd * he), sum(je * hd), sum(je * he)), 2)
      Sig <- (Sig + t(Sig)) / 2; Sig[2, 2] <- Sig[2, 2] + fld$sig_z^2
      list(village = S_te$has[k], y = S_te$y[k], n = S_te$n[k], pmf = pmf_dE(et$d[k], et$e[k], Sig, S_te$n[k]))
    })
  }
  out <- unlist(lapply(sort(unique(folds)), function(k) one(which(folds != k), which(folds == k))), recursive = FALSE)
  out[order(sapply(out, `[[`, "village"))]
}
cv_scores <- function(pr, AI) {
  cond_cdf <- function(pmf) { p0 <- pmf[1]; c(0, cumsum(pmf[-1]) / (1 - p0)) }
  pit_fun <- function(Fm1, F0, u) ifelse(u <= Fm1, 0, ifelse(u >= F0, 1, (u - Fm1) / pmax(F0 - Fm1, 1e-12)))
  ug <- seq(0, 1, 0.01)
  mean_cdf <- function(ys, cdfs) rowMeans(sapply(seq_along(ys), function(i) pit_fun(cdfs[[i]][ys[i]], cdfs[[i]][ys[i] + 1], ug)))
  y <- sapply(pr, `[[`, "y"); p1 <- sapply(pr, function(z) 1 - z$pmf[1]); o1 <- y > 0
  r <- rank(p1); auc <- (sum(r[o1]) - sum(o1) * (sum(o1) + 1) / 2) / (sum(o1) * sum(!o1))
  pos <- Filter(function(z) z$y > 0, pr); cdfs <- lapply(pos, function(z) cond_cdf(z$pmf))
  Fbar <- mean_cdf(sapply(pos, `[[`, "y"), cdfs)
  set.seed(2026)
  sims <- replicate(1000, mean_cdf(sapply(pos, function(z) { q <- z$pmf[-1] / (1 - z$pmf[1]); sample.int(length(q), 1, prob = q) }), cdfs))
  dev <- max(abs(Fbar - ug)); dsim <- apply(abs(sims - ug), 2, max)
  pm <- sapply(pr, function(z) sum((seq_along(z$pmf) - 1) * z$pmf) / z$n); ob <- sapply(pr, function(z) z$y / z$n)
  humid <- AI >= 0.5
  list(table = data.frame(
         sum_log_score = sum(sapply(pr, function(z) log(max(z$pmf[z$y + 1], 1e-300)))),
         AUC = auc, Brier = mean((p1 - o1)^2), zeros_observed = sum(!o1), zeros_expected = sum(1 - p1),
         zeros_humid_observed = sum(!o1 & humid), zeros_humid_expected = sum((1 - p1)[humid]),
         cor_prevalence = stats::cor(pm, ob), pit_max_dev = dev, pit_p = (1 + sum(dsim >= dev)) / (1 + length(dsim))),
       pit = data.frame(u = ug, Fbar = Fbar, lower = apply(sims, 1, quantile, 0.025), upper = apply(sims, 1, quantile, 0.975)))
}

# ---- helper: intervals for the hyperparameters ----------------------------------------
hyper_intervals <- function(fit) {
  S <- fit$setup; cf <- fit$coefficients; p <- ncol(S$X); q <- ncol(S$B)
  x_mode <- list(beta = cf[seq_len(p)], c = cf[p + seq_len(p)], theta = cf[2 * p + seq_len(q)], u = fit$u, z = fit$z)
  p0 <- fit$par                                    # log sigma_u, log sigma_z, a_E, gamma, log range
  ll_at <- function(par) {
    R <- wm_precision(graph, exp(par[5]), 1)
    fld <- list(field = "u", R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE),
                sig_u = exp(par[1]), sig_z = exp(par[2]), aE = par[3], gamma = par[4])
    odispers:::odisp_zil_laplace(fld, S, x_mode)$logLik
  }
  k <- length(p0); h <- c(0.05, 0.05, 0.02, 0.02, 0.05)
  f0 <- ll_at(p0); H <- matrix(0, k, k); fp <- fm <- numeric(k)
  for (i in 1:k) { e <- replace(numeric(k), i, h[i]); fp[i] <- ll_at(p0 + e); fm[i] <- ll_at(p0 - e)
                   H[i, i] <- (fp[i] - 2 * f0 + fm[i]) / h[i]^2 }
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    ei <- replace(numeric(k), i, h[i]); ej <- replace(numeric(k), j, h[j])
    H[i, j] <- H[j, i] <- (ll_at(p0 + ei + ej) - ll_at(p0 + ei - ej) - ll_at(p0 - ei + ej) + ll_at(p0 - ei - ej)) /
      (4 * h[i] * h[j])
  }
  ev <- eigen(-H, symmetric = TRUE)$values
  V <- solve(-H); se <- sqrt(pmax(diag(V), 0)); z <- stats::qnorm(0.975)
  logsc <- c(TRUE, TRUE, FALSE, FALSE, TRUE)
  hyp <- data.frame(parameter = c("sigma_u", "sigma_z", "a_E", "gamma", "range_u (km)"),
                    estimate = ifelse(logsc, exp(p0), p0),
                    lower = ifelse(logsc, exp(p0 - z * se), p0 - z * se), upper = ifelse(logsc, exp(p0 + z * se), p0 + z * se))
  sec <- sqrt(diag(fit$vcov))
  reg <- data.frame(parameter = names(cf), estimate = cf, lower = cf - z * sec, upper = cf + z * sec)
  tab <- rbind(reg, hyp); rownames(tab) <- NULL
  list(table = tab, hessian = H, eigen = ev, gradient = (fp - fm) / (2 * h), logLik_check = f0)
}

# ---- prediction grid (20 km) ----------------------------------------------------------
crs_km <- attr(nodes, "crs")
grid <- st_make_grid(st_transform(sudan_borders, crs_km), cellsize = 20, what = "centers")
grid <- grid[lengths(st_intersects(grid, st_transform(st_union(sudan_borders), crs_km))) > 0]
grid <- st_sf(geometry = st_transform(grid, 4326))
grid$AI <- extract_aridity(grid, sudan_aridity, buffer_km = 10)
grid <- grid[is.finite(grid$AI) & grid$AI > 0, ]
grid_pairs <- river_pairs(grid, nodes, phi)

set.seed(2026)
folds <- sample(rep(1:10, length.out = nrow(v)))   # identical to the folds of the earlier analyses

# ---- run the models -------------------------------------------------------------------
for (m in run_models) {
  log_msg("==== %s: breeding %s", m, deparse(breeding[[m]]))
  prev <- NULL; prev_name <- NA
  for (pm_ in parent[[m]]) if (file.exists(res_file(pm_))) { prev <- readRDS(res_file(pm_))$fit; prev_name <- pm_; break }
  start <- c(log(7.5), log(0.5), 0.5, 0.1); r0 <- 1000; xs <- NULL
  if (!is.null(prev)) {
    start <- prev$par[1:4]; r0 <- prev$hyper[["range_u"]]; cf <- prev$coefficients
    q_new <- length(attr(terms(breeding[[m]]), "term.labels")); th <- cf[grep("^theta\\.", names(cf))]
    xs <- list(beta = cf[1:2], c = cf[3:4], theta = c(th, rep(0, q_new - length(th))), u = prev$u, z = prev$z)
    log_msg("warm start from the fit of %s", prev_name)
  }
  t0 <- Sys.time()
  fit <- odisp_fit_zil(y ~ log(AI), breeding[[m]], v, nodes, graph, pairs, phi = phi, range_u = r0,
                       start = start, estimate_range = TRUE, bfgs = TRUE, x_start = xs, multistart = FALSE,
                       checkpoint = file.path(out_dir, sprintf("checkpoint_%s.rds", m)))
  fit_time <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  print(fit)
  log_msg("%s fitted in %.1f min, logLik %.3f", m, fit_time, fit$logLik)
  slim <- function(f) { f$setup$graph <- NULL; f$fld$Rch <- NULL; f }
  saveRDS(list(fit = slim(fit), time = fit_time), res_file(m))

  t0 <- Sys.time(); ci <- hyper_intervals(fit)
  log_msg("%s intervals in %.1f min; eigenvalues of the observed information: %s", m,
          as.numeric(difftime(Sys.time(), t0, units = "mins")), paste(signif(ci$eigen, 3), collapse = " "))
  print(ci$table, digits = 3)
  write.csv(ci$table, file.path(out_dir, sprintf("intervals_%s.csv", m)), row.names = FALSE)

  t0 <- Sys.time(); pr <- cv_zil(fit, folds); sc <- cv_scores(pr, v$AI[fit$setup$has])
  log_msg("%s cross-validation in %.1f min", m, as.numeric(difftime(Sys.time(), t0, units = "mins")))
  print(sc$table, digits = 4)

  pg <- predict(fit, grid, nodes, graph, pairs = grid_pairs)
  xy <- st_coordinates(grid)
  pred <- data.frame(lon = xy[, 1], lat = xy[, 2], AI = grid$AI, pg)
  saveRDS(list(fit = slim(fit), time = fit_time, intervals = ci, cv = pr, cv_scores = sc, grid = pred), res_file(m))
}

# ---- comparison -----------------------------------------------------------------------
done <- Filter(file.exists, sapply(names(breeding), res_file))
if (length(done)) {
  cmp <- do.call(rbind, lapply(names(done), function(m) {
    r <- readRDS(done[[m]]); f <- r$fit; h <- f$hyper
    row <- data.frame(model = m, breeding = deparse(breeding[[m]]), logLik = f$logLik,
                      n_theta = sum(grepl("^theta\\.", names(f$coefficients))),
                      sigma_u = h[["sigma_u"]], sigma_z = h[["sigma_z"]], a_E = h[["a_E"]], gamma = h[["gamma"]],
                      range_km = h[["range_u"]])
    if (!is.null(r$cv_scores)) row <- cbind(row, r$cv_scores$table)
    row }))
  cmp$d_logLik <- cmp$logLik - cmp$logLik[1]
  print(cmp, digits = 4, row.names = FALSE)
  write.csv(cmp, file.path(out_dir, "comparison.csv"), row.names = FALSE)
}
log_msg("done")
