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
#   cv       10-fold cross-validation of the migration model (same folds as M3)
#   summary  likelihood ratio test, profile interval, comparison table and the
#            share of each village's exposure due to migration,
#            s_i = w E_mig / {(1 - w) E_loc + w E_mig}
#   all      m3, mig, cv and summary
#
# Usage:
#   Rscript m3_migration.R <out_dir> <stage> [w]
# Example SLURM jobs:
#   sbatch --time=24:00:00 --mem=12G --wrap "Rscript m3_migration.R res all"
#   sbatch --array=1-8 --time=24:00:00 --mem=12G --wrap \
#     "W=(0.005 0.01 0.02 0.05 0.1 0.2 0.3 0.5); Rscript m3_migration.R res profile \${W[\$((SLURM_ARRAY_TASK_ID-1))]}"
###############################################################################
suppressPackageStartupMessages({ library(odispers); library(sf); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "m3_migration_results"
stage   <- if (length(args) >= 2) args[2] else "all"
w_fixed <- if (length(args) >= 3) as.numeric(args[3]) else NA
mig_range <- 300    # km, range of the migration kernel (long-range dispersal of S. damnosum s.l.)
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

slim <- function(f) { f$setup$graph <- NULL; f$fld$Rch <- NULL; f }
load_fit <- function(f) { x <- readRDS(f); x$fit$fld$Rch <- Matrix::Cholesky(x$fit$fld$R, LDL = FALSE, perm = TRUE); x }
m3_file <- file.path(out_dir, "fit_M3_refit.rds"); mig_file <- file.path(out_dir, "fit_M3_migration.rds")
breeding <- ~ log_slope + inundation

# published M3 estimates, as starting values
m3_pub <- read.csv(system.file("reports/river_model_covariates/results/intervals_M3.csv", package = "odispers"))
cf_pub <- setNames(m3_pub$estimate, m3_pub$parameter)
start_m3 <- list(coefficients = cf_pub[grep("^(lambda|pi|theta)\\.", names(cf_pub))],
                 par = c(log(cf_pub[["sigma_u"]]), log(cf_pub[["sigma_z"]]), cf_pub[["a_E"]], cf_pub[["gamma"]],
                         log(cf_pub[["range_u (km)"]])))

if (stage %in% c("m3", "all") && !file.exists(m3_file)) {
  log_msg("==== M3 refit")
  fit <- odisp_fit_zil_ml(y ~ log(AI), breeding, v, nodes, graph, pairs, phi = phi, start = start_m3,
                          checkpoint = file.path(out_dir, "checkpoint_M3_refit.rds"))
  print(fit); saveRDS(list(fit = slim(fit)), m3_file)
}

if (stage %in% c("mig", "all") && !file.exists(mig_file)) {
  log_msg("==== M3 + migration kernel (%g km), w estimated", mig_range)
  st <- start_m3
  if (file.exists(m3_file)) { f0 <- readRDS(m3_file)$fit; st <- list(coefficients = f0$coefficients, par = f0$par, u = f0$u, z = f0$z) }
  fit <- odisp_fit_zil_ml(y ~ log(AI), breeding, v, nodes, graph, pairs, phi = phi, start = st,
                          migration_range = mig_range, w_start = 0.05,
                          checkpoint = file.path(out_dir, "checkpoint_M3_migration.rds"))
  print(fit); saveRDS(list(fit = slim(fit)), mig_file)
}

if (stage == "profile") {
  stopifnot(is.finite(w_fixed), w_fixed > 0, w_fixed < 1)
  pf <- file.path(out_dir, sprintf("profile_w_%g.rds", w_fixed))
  if (!file.exists(pf)) {
    log_msg("==== profile: w fixed at %g", w_fixed)
    st <- start_m3
    if (file.exists(mig_file)) { f0 <- readRDS(mig_file)$fit; st <- list(coefficients = f0$coefficients, par = f0$par[1:5], u = f0$u, z = f0$z) }
    fit <- odisp_fit_zil_ml(y ~ log(AI), breeding, v, nodes, graph, pairs, phi = phi, start = st,
                            migration_range = mig_range, fixed_w = w_fixed, intervals = FALSE,
                            checkpoint = file.path(out_dir, sprintf("checkpoint_profile_%g.rds", w_fixed)))
    saveRDS(data.frame(w = w_fixed, logLik = fit$logLik), pf)
    log_msg("w = %g: logLik %.3f", w_fixed, fit$logLik)
  }
}

if (stage %in% c("cv", "all") && file.exists(mig_file)) {
  res <- load_fit(mig_file)
  if (is.null(res$cv_scores)) {
    t0 <- Sys.time(); res$cv <- cv_zil(res$fit, folds); res$cv_scores <- cv_scores(res$cv, v$AI[res$fit$setup$has])
    log_msg("cross-validation in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
    res$fit <- slim(res$fit); saveRDS(res, mig_file)
  }
  print(res$cv_scores$table, digits = 4)
}

if (stage %in% c("summary", "all") && file.exists(m3_file) && file.exists(mig_file)) {
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
  write.csv(f1$intervals, file.path(out_dir, "intervals_M3_migration.csv"), row.names = FALSE)
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
  tab <- data.frame(model = c("M3", "M3 + migration"), logLik = c(f0$logLik, f1$logLik), n_par = c(11, 12))
  tab$AIC <- -2 * tab$logLik + 2 * tab$n_par
  res1 <- readRDS(mig_file)
  if (!is.null(res1$cv_scores)) tab <- cbind(tab, rbind(NA, res1$cv_scores$table))
  print(tab, digits = 5); write.csv(tab, file.path(out_dir, "comparison_migration.csv"), row.names = FALSE)
}
