###############################################################################
# M3 with a migration kernel: exposure from local dispersal (Wendland kernel,
# range 40 km) and long-range migration of the vector (Wendland kernel, range
# 300 km),
#   E_i(w) = (1 - w) E_i^loc + w E_i^mig,
# with the migration share w estimated by maximum likelihood together with
# the 11 parameters of M3 (no priors; Laplace approximation over u and Z).
# Everything else is as in M3 (taildown_covariates.R).
#
# Stages (run in this order; each is checkpointed and skipped when done):
#   m3       refit M3 from the published estimates, so that M3 and the
#            migration model are maximised by the same code (likelihood
#            ratio test on the same footing)
#   mig      fit the migration model, w estimated, started from M3
#   profile  fit with w held at the value given as third argument (profile
#            likelihood interval for w); run one job per value, e.g. as a
#            SLURM array over w in 0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.5
#   cv       10-fold cross-validation of M3 (refit) and of the migration model
#            (same folds as the published M3)
#   predict  predictions on a 20 km grid from both models, with maps
#   summary  parameter estimates with 95% intervals, likelihood ratio test,
#            profile interval for w, comparison table (log-likelihood, AIC and
#            CV scores of both models), PIT figure, and the share of each
#            village's exposure due to migration,
#            s_i = w E_mig / {(1 - w) E_loc + w E_mig}
#   final    after the profile jobs: full fit with w estimated (and its intervals), started
#            from the profile fit with the highest log-likelihood; it replaces the
#            migration fit and is followed by cv, predict and summary in the same job
#   all      m3, mig, cv, predict and summary
#
# Outputs in <out_dir>:
#   estimates_M3_refit.csv                               M3 refit estimates (intervals: published report)
#   intervals_M3_migration.csv                           migration model estimates and 95% intervals
#   profile_w.csv                                        profile log-likelihood for w
#   comparison_migration.csv                             logLik, AIC, CV scores
#   pit_migration.png                                    mean PIT distribution functions
#   grid_predictions.csv, maps_migration.png             20 km grid predictions and maps
#   migration_share_villages.csv                         s_i at the surveyed villages
#
# Usage:
#   Rscript m3_migration.R <out_dir> <stage> [w]
# The likelihood evaluations of each fit run in parallel on SLURM_CPUS_PER_TASK cores
# (or ODISPERS_CORES); each core needs ~2 GB for M3 and ~4 GB for the migration model.
# On SLURM use m3_migration_job.sh / submit_m3_migration.sh (12 cores, 64 GB per fit).
# On a laptop: ODISPERS_CORES=2 Rscript m3_migration.R res m3
###############################################################################
suppressPackageStartupMessages({ library(odispers); library(sf); library(Matrix) })
# PNG files without an X11 display (compute nodes): use the cairo device
if (capabilities("cairo")) options(bitmapType = "cairo")
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "m3_migration_results"
stage   <- if (length(args) >= 2) args[2] else "all"
w_fixed <- if (length(args) >= 3) as.numeric(args[3]) else NA
mig_range <- 300    # km, range of the migration kernel (long-range dispersal of S. damnosum s.l.)
# cores for the parallel finite-difference evaluations: from SLURM, else the
# environment variable ODISPERS_CORES, else 1 (e.g. on a laptop)
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("ODISPERS_CORES", "1")))
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
pairs <- river_pairs(v, nodes, mig_range)        # pairs up to the migration range
graph <- directed_graph(river_graph(nodes, vertex_spacing = 5), "taildown")
log_msg("data: %d villages (%d with positives), %d nodes, %d graph vertices", nrow(v), sum(v$y > 0),
        nrow(nodes), graph$nv)

breeding <- list(M3 = ~ log_slope + inundation)
parent <- list(M3 = character(0))
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
  if (length(folds) != S_all$m_all || anyNA(folds))
    stop("supply one non-missing fold label per original village")
  # The setup omits unreachable villages with zero counts. Keep the original
  # fold assignment for each retained village, in the setup's row order.
  folds <- folds[S_all$has]
  if (length(unique(folds)) < 2L) stop("cross-validation needs at least two folds among reachable villages")
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


set.seed(2026)
folds <- sample(rep(1:10, length.out = nrow(v)))   # identical to the folds of M3

# prediction grid (20 km), as for M3; pairs to the local range only (the migratory
# exposure is computed by predict() from blocks of river nodes)
crs_km <- attr(nodes, "crs")
grid <- st_make_grid(st_transform(sudan_borders, crs_km), cellsize = 20, what = "centers")
grid <- grid[lengths(st_intersects(grid, st_transform(st_union(sudan_borders), crs_km))) > 0]
grid <- st_sf(geometry = st_transform(grid, 4326))
grid$AI <- extract_aridity(grid, sudan_aridity, buffer_km = 10)
grid <- grid[is.finite(grid$AI) & grid$AI > 0, ]

slim <- function(f) { f$setup$graph <- NULL; f$fld$Rch <- NULL; f }
load_fit <- function(f) { x <- readRDS(f); x$fit$fld$Rch <- Matrix::Cholesky(x$fit$fld$R, LDL = FALSE, perm = TRUE); x }
m3_file <- file.path(out_dir, "fit_M3_refit.rds"); mig_file <- file.path(out_dir, "fit_M3_migration.rds")
breeding <- ~ log_slope + inundation

# published M3 estimates (inst/reports/river_model_covariates/results/intervals_M3.csv),
# as starting values; written out here because inst/reports is not part of the installed package
cf_pub <- c("lambda.(Intercept)" = -2.01176703209239, "lambda.log(AI)" = 0.535328302137372,
            "pi.(Intercept)" = 0.489747479889456, "pi.log(AI)" = 2.1661244154256,
            "theta.log_slope" = 1.83234009776781, "theta.inundation" = -2.29993043276964,
            "sigma_u" = 9.43722035886309, "sigma_z" = 0.489024217504059, "a_E" = 0.371065721321979,
            "gamma" = 0.0693967006890684, "range_u (km)" = 897.526926743231)
start_m3 <- list(coefficients = cf_pub[grep("^(lambda|pi|theta)\\.", names(cf_pub))],
                 par = c(log(cf_pub[["sigma_u"]]), log(cf_pub[["sigma_z"]]), cf_pub[["a_E"]], cf_pub[["gamma"]],
                         log(cf_pub[["range_u (km)"]])))

if (stage %in% c("m3", "all") && !file.exists(m3_file)) {
  log_msg("==== M3 refit")
  # only the maximised log-likelihood is needed (the intervals of M3 are in the published report),
  # so no Hessian: intervals = FALSE
  fit <- odisp_fit_zil_ml(y ~ log(AI), breeding, v, nodes, graph, pairs, phi = phi, start = start_m3,
                          intervals = FALSE, cores = n_cores, checkpoint = file.path(out_dir, "checkpoint_M3_refit.rds"))
  print(fit); saveRDS(list(fit = slim(fit)), m3_file)
}

if (stage %in% c("mig", "all") && !file.exists(mig_file)) {
  log_msg("==== M3 + migration kernel (%g km), w estimated", mig_range)
  st <- start_m3
  if (file.exists(m3_file)) { f0 <- readRDS(m3_file)$fit; st <- list(coefficients = f0$coefficients, par = f0$par, u = f0$u, z = f0$z) }
  fit <- odisp_fit_zil_ml(y ~ log(AI), breeding, v, nodes, graph, pairs, phi = phi, start = st,
                          migration_range = mig_range, w_start = 0.05,
                          cores = n_cores, checkpoint = file.path(out_dir, "checkpoint_M3_migration.rds"))
  print(fit); saveRDS(list(fit = slim(fit)), mig_file)
}

if (stage == "profile") {
  stopifnot(is.finite(w_fixed), w_fixed > 0, w_fixed < 1)
  pf <- file.path(out_dir, sprintf("profile_w_%g.rds", w_fixed))
  if (!file.exists(pf)) {
    log_msg("==== profile: w fixed at %g", w_fixed)
    # started from the M3 refit (w = 0), which nests every fit with w fixed; the
    # migration fit is not used as a start, since it may not be at its maximum
    st <- start_m3
    if (file.exists(m3_file)) { f0 <- readRDS(m3_file)$fit; st <- list(coefficients = f0$coefficients, par = f0$par[1:5], u = f0$u, z = f0$z) }
    fit <- odisp_fit_zil_ml(y ~ log(AI), breeding, v, nodes, graph, pairs, phi = phi, start = st,
                            migration_range = mig_range, fixed_w = w_fixed, intervals = FALSE,
                            cores = n_cores, checkpoint = file.path(out_dir, sprintf("checkpoint_profile_%g.rds", w_fixed)))
    saveRDS(data.frame(w = w_fixed, logLik = fit$logLik), pf)
    # estimates at this w, used by the final stage as the start of the full fit
    saveRDS(list(w = w_fixed, logLik = fit$logLik, coefficients = fit$coefficients, par = fit$par[1:5],
                 u = fit$u, z = fit$z), file.path(out_dir, sprintf("profile_start_w_%g.rds", w_fixed)))
    log_msg("w = %g: logLik %.3f", w_fixed, fit$logLik)
  }
}

if (stage == "final") {
  # full fit of the migration model (w estimated, with intervals) started from the best
  # profile fit; it replaces fit_M3_migration.rds, and the cv, predict and summary stages
  # follow in the same job
  sf <- list.files(out_dir, "^profile_start_w_.*\\.rds$", full.names = TRUE)
  if (!length(sf)) stop("no profile_start_w_*.rds files: run the profile stage first")
  starts <- lapply(sf, readRDS)
  ll <- sapply(starts, `[[`, "logLik"); b <- starts[[which.max(ll)]]
  ll0 <- readRDS(m3_file)$fit$logLik
  print(data.frame(w = sapply(starts, `[[`, "w"), logLik = ll)[order(sapply(starts, `[[`, "w")), ], digits = 7, row.names = FALSE)
  log_msg("best profile fit: w = %g, logLik %.3f (M3: %.3f, difference %.3f)", b$w, b$logLik, ll0, b$logLik - ll0)
  if (b$logLik <= ll0 + 0.1) log_msg("NOTE: no value of w improves on M3; the migration component is not supported by the data")
  if (file.exists(mig_file)) file.rename(mig_file, file.path(out_dir, "fit_M3_migration_first.rds"))
  gp_old <- file.path(out_dir, "grid_predictions.csv"); if (file.exists(gp_old)) file.remove(gp_old)
  log_msg("==== M3 + migration kernel (%g km), w estimated, started from the profile fit at w = %g", mig_range, b$w)
  fit <- odisp_fit_zil_ml(y ~ log(AI), breeding, v, nodes, graph, pairs, phi = phi,
                          start = list(coefficients = b$coefficients, par = c(b$par, stats::qlogis(b$w)), u = b$u, z = b$z),
                          migration_range = mig_range, cores = n_cores,
                          checkpoint = file.path(out_dir, "checkpoint_M3_migration_final.rds"))
  print(fit); saveRDS(list(fit = slim(fit)), mig_file)
  if (fit$logLik < b$logLik - 0.1) log_msg("WARNING: the full fit (%.3f) is below the best profile fit (%.3f)", fit$logLik, b$logLik)
  stage <- "final_outputs"
}

if (stage %in% c("cv", "all", "final_outputs")) {
  for (ff in c(m3_file, mig_file)) if (file.exists(ff)) {
    res <- load_fit(ff)
    if (is.null(res$cv_scores)) {
      t0 <- Sys.time(); res$cv <- cv_zil(res$fit, folds); res$cv_scores <- cv_scores(res$cv, v$AI[res$fit$setup$has])
      log_msg("%s: cross-validation in %.1f min", basename(ff), as.numeric(difftime(Sys.time(), t0, units = "mins")))
      res$fit <- slim(res$fit); saveRDS(res, ff)
    }
    print(res$cv_scores$table, digits = 4)
  }
}

if (stage %in% c("predict", "all", "final_outputs") && file.exists(m3_file) && file.exists(mig_file)) {
  gp_file <- file.path(out_dir, "grid_predictions.csv")
  if (!file.exists(gp_file)) {
    grid_pairs <- river_pairs(grid, nodes, phi)
    xy <- st_coordinates(grid)
    p0 <- predict(load_fit(m3_file)$fit, grid, nodes, graph, pairs = grid_pairs)
    p1 <- predict(load_fit(mig_file)$fit, grid, nodes, graph, pairs = grid_pairs)
    gp <- data.frame(lon = xy[, 1], lat = xy[, 2], AI = grid$AI,
                     M3_E = p0$E, M3_pi = p0$pi, M3_p_zone = p0$p_zone, M3_prevalence = p0$prevalence,
                     mig_E = p1$E, mig_E_loc = p1$E_loc, mig_E_mig = p1$E_mig, mig_s_mig = p1$s_mig,
                     mig_pi = p1$pi, mig_p_zone = p1$p_zone, mig_prevalence = p1$prevalence)
    write.csv(gp, gp_file, row.names = FALSE)
    log_msg("grid predictions written (%d cells)", nrow(gp))
  }
  gp <- read.csv(gp_file)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    long <- rbind(data.frame(lon = gp$lon, lat = gp$lat, panel = "M3: expected prevalence", value = gp$M3_prevalence),
                  data.frame(lon = gp$lon, lat = gp$lat, panel = "M3 + migration: expected prevalence", value = gp$mig_prevalence),
                  data.frame(lon = gp$lon, lat = gp$lat, panel = "Difference (migration - M3)", value = gp$mig_prevalence - gp$M3_prevalence),
                  data.frame(lon = gp$lon, lat = gp$lat, panel = "Share of exposure due to migration", value = gp$mig_s_mig))
    border <- st_transform(sudan_borders, 4326)
    mk <- function(d, lab, fill) ggplot2::ggplot() +
      ggplot2::geom_tile(data = d, ggplot2::aes(lon, lat, fill = value), width = 0.19, height = 0.19) +
      ggplot2::geom_sf(data = border, fill = NA, colour = "grey30", linewidth = 0.2) + fill +
      ggplot2::labs(title = lab, x = NULL, y = NULL, fill = NULL) + ggplot2::theme_minimal(base_size = 9)
    seq_fill <- ggplot2::scale_fill_viridis_c(option = "magma", direction = -1, na.value = "grey85")
    div_fill <- ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0)
    pl <- list(mk(long[long$panel == unique(long$panel)[1], ], unique(long$panel)[1], seq_fill),
               mk(long[long$panel == unique(long$panel)[2], ], unique(long$panel)[2], seq_fill),
               mk(long[long$panel == unique(long$panel)[3], ], unique(long$panel)[3], div_fill),
               mk(long[long$panel == unique(long$panel)[4], ], unique(long$panel)[4],
                  ggplot2::scale_fill_viridis_c(limits = c(0, 1), na.value = "grey85")))
    grDevices::png(file.path(out_dir, "maps_migration.png"), width = 10, height = 9, units = "in", res = 200)
    grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
    for (k in 1:4) print(pl[[k]], vp = grid::viewport(layout.pos.row = (k - 1) %/% 2 + 1, layout.pos.col = (k - 1) %% 2 + 1))
    grDevices::dev.off()
    log_msg("maps written to maps_migration.png")
  }
}

if (stage %in% c("summary", "all", "final_outputs") && file.exists(m3_file) && file.exists(mig_file)) {
  f0 <- readRDS(m3_file)$fit; f1 <- readRDS(mig_file)$fit
  lr <- 2 * (f1$logLik - f0$logLik)
  p <- 0.5 * stats::pchisq(max(lr, 0), 1, lower.tail = FALSE)        # w = 0 is on the boundary
  log_msg("M3 logLik %.3f; migration model logLik %.3f; w = %.4f", f0$logLik, f1$logLik, f1$hyper[["w"]])
  log_msg("likelihood ratio %.2f, p = %.3g (50:50 mixture of chi2_0 and chi2_1)", lr, p)
  prof <- do.call(rbind, lapply(list.files(out_dir, "^profile_w_.*\\.rds$", full.names = TRUE), readRDS))
  if (!is.null(prof)) {
    prof <- rbind(data.frame(w = 0, logLik = f0$logLik), prof, data.frame(w = f1$hyper[["w"]], logLik = f1$logLik))
    prof <- prof[order(prof$w), ]; prof$drop <- max(prof$logLik) - prof$logLik
    print(prof, digits = 6)
    inside <- prof$w[prof$drop <= stats::qchisq(0.95, 1) / 2]
    log_msg("approximate 95%% profile interval for w (grid): %.3g to %.3g", min(inside), max(inside))
    write.csv(prof, file.path(out_dir, "profile_w.csv"), row.names = FALSE)
  }
  write.csv(f0$intervals[, c("parameter", "estimate")], file.path(out_dir, "estimates_M3_refit.csv"), row.names = FALSE)
  write.csv(f1$intervals, file.path(out_dir, "intervals_M3_migration.csv"), row.names = FALSE)
  print(f1$intervals, digits = 4)
  # share of each village's exposure due to migration, s_i = w E_mig / {(1 - w) E_loc + w E_mig},
  # at the estimates and the mode of u
  w_hat <- f1$hyper[["w"]]; Su <- f1$setup
  b <- exp(as.vector(Su$B %*% f1$coefficients[grep("^theta\\.", names(f1$coefficients))]) + as.vector(Su$A %*% f1$u))
  E_loc <- as.vector(Su$K1 %*% b); E_mig <- as.vector(Su$K2 %*% b)
  s_mig <- w_hat * E_mig / ((1 - w_hat) * E_loc + w_hat * E_mig)
  log_msg("migration share of exposure s_i: median %.3f, quartiles %.3f-%.3f; villages with s_i > 0.5: %d of %d",
          stats::median(s_mig), stats::quantile(s_mig, 0.25), stats::quantile(s_mig, 0.75), sum(s_mig > 0.5), length(s_mig))
  vv <- sf::st_drop_geometry(v)[Su$has, c("long", "lat", "y", "n")]
  write.csv(data.frame(vv, E_loc = E_loc, E_mig = E_mig, s_mig = s_mig),
            file.path(out_dir, "migration_share_villages.csv"), row.names = FALSE)
  tab <- data.frame(model = c("M3", "M3 + migration"), logLik = c(f0$logLik, f1$logLik), n_par = c(11, 12),
                    w = c(0, f1$hyper[["w"]]))
  tab$AIC <- -2 * tab$logLik + 2 * tab$n_par
  r0 <- readRDS(m3_file); r1 <- readRDS(mig_file)
  if (!is.null(r0$cv_scores) && !is.null(r1$cv_scores)) {
    tab <- cbind(tab, rbind(r0$cv_scores$table, r1$cv_scores$table))
    # paired difference in log score (migration - M3), with a 1-degree block bootstrap SE
    ls_v <- function(cv) sapply(cv, function(z) log(max(z$pmf[z$y + 1], 1e-300)))
    d_ls <- ls_v(r1$cv) - ls_v(r0$cv)
    cell <- with(sf::st_drop_geometry(v)[r1$fit$setup$has, ], paste(floor(long), floor(lat)))
    ids <- split(seq_along(d_ls), cell); set.seed(1)
    bs <- replicate(2000, sum(d_ls[unlist(ids[sample(length(ids), replace = TRUE)])]))
    log_msg("CV log score, migration - M3: %.2f (block-bootstrap SE %.2f)", sum(d_ls), stats::sd(bs))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      pit <- rbind(cbind(model = "M3", r0$cv_scores$pit), cbind(model = "M3 + migration", r1$cv_scores$pit))
      tst <- data.frame(model = c("M3", "M3 + migration"),
                        lab = sprintf("max deviation %.3f\np = %.3f", c(r0$cv_scores$table$pit_max_dev, r1$cv_scores$table$pit_max_dev),
                                      c(r0$cv_scores$table$pit_p, r1$cv_scores$table$pit_p)))
      g <- ggplot2::ggplot(pit, ggplot2::aes(u)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "grey85") +
        ggplot2::geom_abline(linetype = 2) + ggplot2::geom_line(ggplot2::aes(y = Fbar), colour = "#1b7837", linewidth = 0.8) +
        ggplot2::geom_text(data = tst, ggplot2::aes(x = 0.03, y = 0.95, label = lab), hjust = 0, vjust = 1, size = 3) +
        ggplot2::facet_wrap(~ model) + ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        ggplot2::labs(x = "u", y = "Mean PIT distribution function") + ggplot2::theme_minimal(base_size = 10)
      ggplot2::ggsave(file.path(out_dir, "pit_migration.png"), g, width = 7.5, height = 4, dpi = 200, bg = "white")
      log_msg("PIT figure written to pit_migration.png")
    }
  }
  print(tab, digits = 5); write.csv(tab, file.path(out_dir, "comparison_migration.csv"), row.names = FALSE)
}
