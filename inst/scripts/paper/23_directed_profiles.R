###############################################################################
# Fair comparison of the directed breeding fields in the zero-inflated model
# (tail-down, tail-up K1 with discharge weights, and their sum); the
# symmetric Whittle-Matern field, whose variance is inflated at the river
# ends by the Kirchhoff conditions, is no longer considered.
#
#   profile <type>  profile likelihood over a dense grid of ranges of u (each
#                   point maximised over the other parameters); adds to the
#                   grid fits of 21_fit_directed.R
#   joint <type>    range estimated jointly with the other parameters, from
#                   the maximum of the profile (for the sum of the two fields,
#                   a separate range for each part and the share of variance
#                   phi, from the tail-down and tail-up estimates)
#   cv              random-fold cross-validation of the joint fits
#   compare         profile curves (natural spline in log range), estimates,
#                   likelihood intervals for the range, scores and figures
#
# Every fit writes a checkpoint after each improvement, so an interrupted run
# resumes where it stopped. Usage: Rscript 23_directed_profiles.R <step> [type]
###############################################################################

library(odispers)
library(sf)
library(Matrix)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

args <- commandArgs(trailingOnly = TRUE)
step <- args[1]; type <- if (length(args) > 1) args[2] else NA
grids <- list(taildown = c(100, 200, 400, 800, 1600, 3200, 6400),
              tailup = c(25, 50, 100, 200, 400, 800, 3200))
types <- c("taildown", "tailup", "tailup+taildown")
labels <- c(taildown = "Tail-down", tailup = "Tail-up (K1, discharge weights)",
            "tailup+taildown" = "Tail-up + tail-down")
slug <- function(t) gsub("\\+", "_", t)
grid_file <- function(t) file.path(out_dir, sprintf("directed_%s.rds", slug(t)))
joint_file <- function(t) file.path(out_dir, sprintf("directed_%s_joint.rds", slug(t)))
ck_file <- function(t, what) file.path(out_dir, sprintf("ck_%s_%s.rds", slug(t), what))

inp <- load_result("inputs.rds"); v <- inp$villages
net <- load_result("inputs_full_network.rds")
nodes <- net$nodes
graph <- river_graph(nodes, vertex_spacing = 5)
fit_one <- function(t, r, start, ...) {
  f <- odisp_fit_zil(y ~ log(AI), ~ log(DIS_AV_CMS), v, nodes, directed_graph(graph, t), net$pairs,
                     phi = phi_main, range_u = r, start = start, bfgs = TRUE, ...)
  f$setup$graph <- NULL; f$fld$Rch <- NULL
  f
}
grid_fits <- function(t) if (file.exists(grid_file(t))) readRDS(grid_file(t)) else list()

# ---- profile over the range --------------------------------------------------------
if (step == "profile") {
  fits <- grid_fits(type)
  for (r in grids[[type]]) {
    key <- as.character(r)
    if (!is.null(fits[[key]])) next
    near <- as.numeric(names(fits))[which.min(abs(log(as.numeric(names(fits))) - log(r)))]
    st <- fits[[as.character(near)]]$par
    cat(sprintf("\n==== %s, range %g km (start from the fit at %g km) ====\n", type, r, near))
    t0 <- Sys.time()
    f <- fit_one(type, r, st, checkpoint = ck_file(type, key))
    f$time <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    print(f)
    fits <- grid_fits(type); fits[[key]] <- f; saveRDS(fits, grid_file(type))
  }
}

# ---- joint estimation of the range ----------------------------------------------------
if (step == "joint") {
  if (type == "tailup+taildown") {
    td <- readRDS(joint_file("taildown")); tu <- readRDS(joint_file("tailup"))
    st <- c(td$par[c("sigma_u", "sigma_z", "a_E", "gamma")], mix = stats::qlogis(0.05))
    r0 <- td$hyper[["range_u"]]; ru0 <- tu$hyper[["range_u"]]
  } else {
    fits <- grid_fits(type); ll <- sapply(fits, `[[`, "logLik")
    b <- fits[[which.max(ll)]]; st <- b$par; r0 <- b$hyper[["range_u"]]; ru0 <- r0
  }
  cat(sprintf("\n==== %s, joint estimation from range %g km ====\n", type, r0))
  t0 <- Sys.time()
  f <- fit_one(type, r0, st, estimate_range = TRUE, range_up = ru0, checkpoint = ck_file(type, "joint"))
  f$time <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  print(f); print(f$hyper)
  saveRDS(f, joint_file(type))
}

# ---- cross-validation (random folds of 12_cv_calibration.R) ---------------------------
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
pmf_dE <- function(md, me, Sig, n) {
  L <- t(chol(Sig + diag(1e-10, 2)))
  d <- md + L[1, 1] * g2x1; e <- me + L[2, 1] * g2x1 + L[2, 2] * g2x2
  pi_ <- stats::plogis(d); p <- -expm1(-exp(pmin(e, 30)))
  pm <- sapply(0:n, function(y) sum(g2w * pi_ * stats::dbinom(y, n, p)))
  pm[1] <- pm[1] + sum(g2w * (1 - pi_))
  pm
}
cv_zil <- function(fit, folds) {
  S_all <- fit$setup; fld <- fit$fld
  fld$Rch <- Matrix::Cholesky(fld$R, LDL = FALSE, perm = TRUE)
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
      list(village = test[k], y = S_te$y[k], n = S_te$n[k], pmf = pmf_dE(et$d[k], et$e[k], Sig, S_te$n[k]))
    })
  }
  out <- unlist(lapply(sort(unique(folds)), function(k) { cat(k, ""); one(which(folds != k), which(folds == k)) }),
                recursive = FALSE)
  cat("\n")
  out[order(sapply(out, `[[`, "village"))]
}
cv_file <- file.path(out_dir, "cv_directed_joint.rds")
if (step == "cv") {
  folds <- load_result("cv_calibration.rds")$folds[["random folds"]]
  preds <- if (file.exists(cv_file)) readRDS(cv_file) else list()
  for (t in types) {
    if (!file.exists(joint_file(t)) || !is.null(preds[[t]])) next
    cat("cross-validation,", t, ": ")
    preds[[t]] <- cv_zil(readRDS(joint_file(t)), folds)
    saveRDS(preds, cv_file)
  }
}

# ---- comparison -------------------------------------------------------------------------
if (step == "compare") {
  joint <- list(); for (t in types) if (file.exists(joint_file(t))) joint[[t]] <- readRDS(joint_file(t))
  preds <- readRDS(cv_file)
  # profile curves and likelihood intervals for the range
  prof <- NULL; rtab <- NULL
  for (t in c("taildown", "tailup")) {
    fits <- grid_fits(t)
    d <- data.frame(model = labels[[t]], range = as.numeric(names(fits)), logLik = sapply(fits, `[[`, "logLik"))
    d <- d[order(d$range), ]; prof <- rbind(prof, d)
    sp <- stats::splinefun(log(d$range), d$logLik, method = "natural")
    o <- stats::optimize(sp, range(log(d$range)), maximum = TRUE)
    thr <- o$objective - stats::qchisq(0.95, 1) / 2
    lo <- tryCatch(exp(stats::uniroot(function(x) sp(x) - thr, c(log(min(d$range)), o$maximum))$root), error = function(e) NA)
    hi <- tryCatch(exp(stats::uniroot(function(x) sp(x) - thr, c(o$maximum, log(max(d$range))))$root), error = function(e) NA)
    rtab <- rbind(rtab, data.frame(model = labels[[t]], range_profile = round(exp(o$maximum)),
                                   max_profile = round(o$objective, 2), lower95 = round(lo), upper95 = round(hi),
                                   range_joint = if (!is.null(joint[[t]])) round(joint[[t]]$hyper[["range_u"]]) else NA,
                                   logLik_joint = if (!is.null(joint[[t]])) round(joint[[t]]$logLik, 2) else NA))
  }
  cat("\nRange of u: maximum of the spline-interpolated profile, 95% likelihood interval, and joint estimate:\n")
  print(rtab, row.names = FALSE)
  pc <- do.call(rbind, lapply(split(prof, prof$model), function(d) {
    sp <- stats::splinefun(log(d$range), d$logLik, method = "natural")
    xx <- seq(log(min(d$range)), log(max(d$range)), length.out = 200)
    data.frame(model = d$model[1], range = exp(xx), logLik = sp(xx)) }))
  p_prof <- ggplot(prof, aes(range, logLik, colour = model)) + geom_line(data = pc) + geom_point(size = 2) +
    scale_x_log10() + labs(x = "Range of the breeding field (km, log scale)", y = "Approximate log-likelihood",
                           colour = NULL) + theme_minimal(base_size = 11) + theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "directed_profiles.png"), p_prof, width = 7, height = 4.5, dpi = 200, bg = "white")

  cond_cdf <- function(pmf) { p0 <- pmf[1]; c(0, cumsum(pmf[-1]) / (1 - p0)) }
  pit_fun <- function(Fm1, F0, u) ifelse(u <= Fm1, 0, ifelse(u >= F0, 1, (u - Fm1) / pmax(F0 - Fm1, 1e-12)))
  mean_pit_hist <- function(ys, cdfs, J = 10) {
    ub <- (0:J) / J
    Fbar <- rowMeans(sapply(seq_along(ys), function(i) pit_fun(cdfs[[i]][ys[i]], cdfs[[i]][ys[i] + 1], ub)))
    J * diff(Fbar)
  }
  auc <- function(pp, o) { r <- rank(pp); n1 <- sum(o); n0 <- sum(!o); (sum(r[o]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
  set.seed(2026)
  cls <- aridity_classes(v$AI[joint[[1]]$setup$has])
  scores <- list(); pit_tab <- list()
  for (t in names(joint)) {
    b <- joint[[t]]; pr <- preds[[t]]; h_ <- b$hyper
    row <- data.frame(model = labels[[t]], logLik = round(b$logLik, 2),
                      range_km = round(h_[["range_u"]]),
                      range_up_km = if ("range_up" %in% names(h_)) round(h_[["range_up"]]) else NA,
                      phi = if ("phi" %in% names(h_)) round(h_[["phi"]], 3) else NA,
                      sigma_u = round(h_[["sigma_u"]], 2), sigma_z = round(h_[["sigma_z"]], 3),
                      a_E = round(h_[["a_E"]], 3), gamma = round(h_[["gamma"]], 3),
                      theta = round(b$coefficients[["theta.log(DIS_AV_CMS)"]], 3))
    if (!is.null(pr)) {
      p1 <- sapply(pr, function(z) 1 - z$pmf[1]); o1 <- sapply(pr, `[[`, "y") > 0
      pos <- Filter(function(z) z$y > 0, pr); cdfs <- lapply(pos, function(z) cond_cdf(z$pmf))
      hh <- mean_pit_hist(sapply(pos, `[[`, "y"), cdfs)
      sims <- replicate(500, { ys <- sapply(pos, function(z) { qq <- z$pmf[-1] / (1 - z$pmf[1]); sample.int(length(qq), 1, prob = qq) })
                               mean_pit_hist(ys, cdfs) })
      pit_tab[[t]] <- data.frame(model = labels[[t]], mid = (seq_along(hh) - 0.5) / 10, density = hh,
                                 lower = apply(sims, 1, stats::quantile, 0.025), upper = apply(sims, 1, stats::quantile, 0.975))
      pm <- sapply(pr, function(z) sum((seq_along(z$pmf) - 1) * z$pmf) / z$n); ob <- sapply(pr, function(z) z$y / z$n)
      row <- cbind(row, data.frame(AUC = round(auc(p1, o1), 3), Brier = round(mean((p1 - o1)^2), 4),
                                   log_score = round(mean(sapply(pr, function(z) log(max(z$pmf[z$y + 1], 1e-300)))), 4),
                                   sum_log_score = round(sum(sapply(pr, function(z) log(max(z$pmf[z$y + 1], 1e-300)))), 1),
                                   zeros_expected = round(sum(1 - p1), 1),
                                   zeros_humid = round(sum((1 - p1)[cls == "Sub-humid/humid"]), 1),
                                   cor_cv = round(stats::cor(pm, ob), 3),
                                   pit_bins_out = sum(hh < pit_tab[[t]]$lower | hh > pit_tab[[t]]$upper)))
    }
    scores[[t]] <- row
  }
  scores <- do.call(rbind, lapply(scores, function(r) { r[setdiff(Reduce(union, lapply(scores, names)), names(r))] <- NA; r }))
  rownames(scores) <- NULL
  cat("\nJoint fits and random-fold cross-validation (550 villages without positives, 14 humid):\n")
  print(scores, row.names = FALSE)
  if (length(pit_tab)) {
    pit_tab <- do.call(rbind, pit_tab); pit_tab$model <- factor(pit_tab$model, unique(pit_tab$model))
    p_pit <- ggplot(pit_tab, aes(mid, density)) + geom_col(width = 0.1, fill = "#7fbf7b", colour = "white") +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.04, colour = "grey30") +
      geom_hline(yintercept = 1, linetype = 2) + facet_wrap(~ model, nrow = 1) +
      labs(x = "Non-randomised PIT (villages with at least one positive)", y = "Density") + theme_minimal(base_size = 10)
    ggsave(file.path(fig_dir, "directed_joint_pit.png"), p_pit, width = 11, height = 3.4, dpi = 200, bg = "white")
  }
  saveRDS(list(profile = prof, range = rtab, scores = scores, pit = pit_tab), file.path(out_dir, "directed_comparison.rds"))
}
