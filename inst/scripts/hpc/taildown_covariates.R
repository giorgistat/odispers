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
#   log w_j = beta' h_j + u(s_j)                              (breeding intensity)
# u: tail-down Gaussian field on the river network (exponential covariance in
# river distance), Z_i ~ N(0, sigma_z^2) independent village effects.
#
# Models (breeding covariates h):
#   M1  log discharge + log(1 + channel slope, dm/km) + min{log(1 + upstream inundation %), 1.5}
#   M3  log(1 + channel slope) + min{log(1 + upstream inundation %), 1.5}   (M1 without discharge)
#
# Estimation: maximum likelihood, with no prior on any parameter. The
# likelihood of the fixed parameters (regression coefficients of both parts
# and of the breeding intensity, sigma_u, sigma_z, a_E, gamma and the range of
# u) is the Laplace approximation of the integral over the random effects u
# and Z only (odisp_fit_zil_ml). All parameters are maximised jointly; Wald
# 95% intervals for every parameter come from the observed information
# (central finite differences), on the log scale for sigma_u, sigma_z and the
# range.
#
# Steps for each model: (1) maximum likelihood fit and intervals; (2) 10-fold
# random cross-validation (folds as in the earlier analyses, set.seed(2026)),
# with all parameters at the full-data estimates and the random effects
# refitted to the training villages: log score, AUC and Brier score for any
# positive, expected zeros, non-randomised PIT (mean distribution function and
# its largest deviation from uniform, with a simulated p-value); (3)
# predictions on a 20 km grid. Then a comparison table, with the likelihood
# ratio test of M3 against M1.
#
# Usage:
#   Rscript taildown_covariates.R [out_dir] [models]
#     out_dir  folder for the results (default taildown_covariates_results)
#     models   comma-separated subset of M1,M3 (default both, in order)
# M3 starts from the M1 fit when it is in out_dir. A file start_<model>.rds in
# out_dir (a list with coefficients, par, u, z) is used as starting point if
# present. The fit is checkpointed and every completed stage is saved: a job
# that is stopped resumes and skips the stages already done.
# Example SLURM job (one process per model; several hours per model):
#   #!/bin/bash
#   #SBATCH --time=24:00:00 --mem=8G --cpus-per-task=1
#   module load R
#   Rscript taildown_covariates.R results M1
###############################################################################
suppressPackageStartupMessages({ library(odispers); library(sf); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "taildown_covariates_results"
run_models <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("M1", "M3")
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

breeding <- list(M1 = ~ log(DIS_AV_CMS) + log_slope + inundation,
                 M3 = ~ log_slope + inundation)
parent <- list(M1 = character(0), M3 = "M1")
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
slim <- function(f) { f$setup$graph <- NULL; f$fld$Rch <- NULL; f }
for (m in run_models) {
  log_msg("==== %s: breeding %s", m, deparse(breeding[[m]]))
  res <- if (file.exists(res_file(m))) readRDS(res_file(m)) else list()
  if (is.null(res$fit)) {
    st <- list()
    sf <- file.path(out_dir, sprintf("start_%s.rds", m))
    if (file.exists(sf)) { st <- readRDS(sf); log_msg("starting values from %s", sf) } else
      for (pm_ in parent[[m]]) if (file.exists(res_file(pm_))) {
        f0 <- readRDS(res_file(pm_))$fit
        st <- list(coefficients = f0$coefficients, par = f0$par, u = f0$u, z = f0$z)
        log_msg("starting values from the fit of %s", pm_); break
      }
    fit <- odisp_fit_zil_ml(y ~ log(AI), breeding[[m]], v, nodes, graph, pairs, phi = phi, start = st,
                            checkpoint = file.path(out_dir, sprintf("checkpoint_%s.rds", m)))
    print(fit)
    log_msg("%s fitted in %.1f min (%d likelihood evaluations), logLik %.3f; eigenvalues of the observed information %s",
            m, fit$time, fit$n_eval, fit$logLik, paste(signif(range(fit$eigen), 3), collapse = " to "))
    res$fit <- slim(fit); saveRDS(res, res_file(m))
  } else {
    fit <- res$fit; fit$fld$Rch <- Matrix::Cholesky(fit$fld$R, LDL = FALSE, perm = TRUE)
    log_msg("%s: fit loaded (logLik %.3f)", m, fit$logLik)
  }
  write.csv(fit$intervals, file.path(out_dir, sprintf("intervals_%s.csv", m)), row.names = FALSE)

  if (is.null(res$cv_scores)) {
    t0 <- Sys.time(); res$cv <- cv_zil(fit, folds); res$cv_scores <- cv_scores(res$cv, v$AI[fit$setup$has])
    log_msg("%s cross-validation in %.1f min", m, as.numeric(difftime(Sys.time(), t0, units = "mins")))
    saveRDS(res, res_file(m))
  }
  print(res$cv_scores$table, digits = 4)

  if (is.null(res$grid)) {
    pg <- predict(fit, grid, nodes, graph, pairs = grid_pairs)
    xy <- st_coordinates(grid)
    res$grid <- data.frame(lon = xy[, 1], lat = xy[, 2], AI = grid$AI, pg)
    saveRDS(res, res_file(m))
  }
}

# ---- comparison -----------------------------------------------------------------------
done <- Filter(file.exists, sapply(names(breeding), res_file))
if (length(done)) {
  cmp <- lapply(names(done), function(m) {
    r <- readRDS(done[[m]]); f <- r$fit; h <- f$hyper
    row <- data.frame(model = m, breeding = deparse(breeding[[m]]), logLik = f$logLik,
                      n_theta = sum(grepl("^theta\\.", names(f$coefficients))),
                      sigma_u = h[["sigma_u"]], sigma_z = h[["sigma_z"]], a_E = h[["a_E"]], gamma = h[["gamma"]],
                      range_km = h[["range_u"]])
    if (!is.null(r$cv_scores)) row <- cbind(row, r$cv_scores$table)
    row })
  cols <- Reduce(union, lapply(cmp, names))
  cmp <- do.call(rbind, lapply(cmp, function(r) { r[setdiff(cols, names(r))] <- NA; r[cols] }))
  cmp$n_par <- cmp$n_theta + 4 + 5                     # coefficients of both parts, breeding, 5 others
  print(cmp, digits = 5, row.names = FALSE)
  if (all(c("M1", "M3") %in% cmp$model)) {
    lr <- 2 * (cmp$logLik[cmp$model == "M1"] - cmp$logLik[cmp$model == "M3"])
    log_msg("likelihood ratio test of M3 against M1 (log discharge): statistic %.3f on 1 df, p = %.3f",
            lr, stats::pchisq(max(lr, 0), 1, lower.tail = FALSE))
  }
  write.csv(cmp, file.path(out_dir, "comparison.csv"), row.names = FALSE)
}
log_msg("done")
