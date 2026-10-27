###############################################################################
# Geostatistical comparators of the river-network model M3: the same
# zero-inflated binomial model, with the exposure log E replaced by Euclidean
# Gaussian processes that use no river information. Self-contained: it needs
# only the package (remotes::install_github("giorgistat/odispers")) and does
# not use the output of any other script.
#
# Model, for village i with n_i examined and y_i with nodules:
#   P(Y_i = 0) = 1 - pi_i + pi_i (1 - p_i)^n_i,   p_i = 1 - exp(-lambda_i)
#   G1 (one process):
#     logit pi_i   = c0 + c1 log AI_i + alpha S(x_i)
#     log lambda_i = b0 + b1 log AI_i + gamma S(x_i) + Z_i
#   G2 (two independent processes):
#     logit pi_i   = c0 + c1 log AI_i + alpha S1(x_i)
#     log lambda_i = b0 + b1 log AI_i + gamma S2(x_i) + Z_i
# S, S1, S2: zero mean, unit variance, exponential correlation
# exp(-2 d / r) in straight-line distance d (km), the same parametrisation of
# the range r as the tail-down field of M3 (correlation 0.14 at d = r); in G2
# S1 and S2 share the range r. Z_i ~ N(0, sigma_z^2) independent village
# effects, as in M3. Constraints: alpha > 0 (G1, sign of S); alpha, gamma > 0
# (G2). Parameters: G1 c0, c1, b0, b1, alpha, gamma, r, sigma_z (8);
# G2 the same with gamma > 0 (8).
#
# Estimation: maximum likelihood, with no prior on any parameter. The
# likelihood of the fixed parameters is the Laplace approximation of the
# integral over the random effects only (the processes at the 900 villages and
# Z), in whitened form: S = L w with C = L L', w ~ N(0, I), Z = sigma_z zeta.
# The inner problem can have several local modes: each evaluation starts the
# inner Newton from the mode at the previous evaluation and from the mode at
# the best point so far (and from zero at the first evaluation and at each
# point where the Hessian is computed), and takes the mode with the highest
# log integrand; the finite differences for the Hessian all start from the mode at
# the centre point (one smooth branch). The log-likelihood is maximised by
# BFGS, then by Newton steps with the full finite-difference Hessian; Wald 95%
# intervals for every parameter come from the observed information, on the
# log scale for the positive parameters. The log-likelihood includes the
# binomial coefficients and is on the same scale as that of
# odisp_fit_zil_ml (M3).
#
# Steps for each model: (1) maximum likelihood fit and intervals; (2) 10-fold
# random cross-validation with the folds of M3 (set.seed(2026)): parameters at
# the full-data estimates, random effects refitted to the training villages,
# the processes at a held-out village predicted by kriging from the training
# villages, predictive distribution of the count from the Gaussian
# approximation of (logit pi, log lambda), a new village effect and 20 x 20
# Gauss-Hermite quadrature: log score, AUC, Brier score, expected zeros,
# non-randomised PIT; (3) plug-in predictions on the 20 km grid of M3.
#
# Outputs in out_dir, all small (a few MB at most):
#   gp_<model>.rds        list: estimates, intervals, vcov, logLik, AIC, modes,
#                         cross-validated predictive pmfs and scores, grid
#   intervals_<model>.csv estimates and 95% intervals
#   cv_<model>.csv        cross-validation scores
#   pit_<model>.csv       mean PIT distribution function with 95% band
#   gp_checkpoint_<model>.rds  resume file (can be deleted at the end)
# Send gp_<model>.rds (or the csv files) for the summaries.
#
# Usage:
#   Rscript geostat_gp_models.R [out_dir] [models] [test]
#     out_dir  folder for the results (default geostat_gp_results)
#     models   comma-separated subset of G1,G2 (default both, in order)
#     test     optional: a short run to check the set-up (few iterations,
#              two folds); its output is not for use
# A stopped job resumes from the checkpoint and skips the stages already done.
# Example SLURM job (one process per model; a few hours each):
#   #!/bin/bash
#   #SBATCH --time=12:00:00 --mem=8G --cpus-per-task=1
#   module load R      # an optimised BLAS (OpenBLAS, MKL) makes this much faster
#   Rscript geostat_gp_models.R results G1
###############################################################################
suppressPackageStartupMessages({ library(odispers); library(sf) })
args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "geostat_gp_results"
run_models <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else c("G1", "G2")
test <- length(args) >= 3 && args[3] == "test"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

# ---- 0. data (as for M3) --------------------------------------------------------------
data(remo); data(sudan_rivers); data(sudan_aridity); data(sudan_borders)
v <- subset(remo, country == "SUDAN")
v$y <- round(v$y)                                   # a few counts in the file are not integers
v <- st_as_sf(v, coords = c("long", "lat"), crs = 4326, remove = FALSE)
v$AI <- extract_aridity(v, sudan_aridity)
stopifnot(all(is.finite(v$AI)), all(v$AI > 0))
rivers <- sudan_rivers[sudan_rivers$DIS_AV_CMS >= 0.1, ]
crs_km <- odispers:::local_aeqd(rivers)             # the projection (km) of the M3 analysis
xy <- st_coordinates(st_transform(st_geometry(v), crs_km))[, 1:2]
Dm <- as.matrix(stats::dist(xy))
X <- stats::model.matrix(~ log(AI), st_drop_geometry(v))
S_all <- list(y = v$y, n = v$n, X = X)
nv_ <- nrow(v)
log_msg("data: %d villages (%d with positives)", nv_, sum(v$y > 0))

spec <- list(G1 = list(K = 1, ke = 1, gamma_pos = FALSE),
             G2 = list(K = 2, ke = 2, gamma_pos = TRUE))
nm_par <- c("lambda.(Intercept)", "lambda.log(AI)", "pi.(Intercept)", "pi.log(AI)",
            "log_alpha", "gamma", "log_range", "log_sigma_z")
unpack <- function(th, sp) list(b = th[1:2], c = th[3:4], alpha = exp(th[5]),
                                 gamma = if (sp$gamma_pos) exp(th[6]) else th[6],
                                 r = exp(th[7]), sz = exp(th[8]))
corr <- function(D, r) exp(-2 * D / r)

# ---- Laplace approximation --------------------------------------------------------------
# observed information of (d, e) per village, floored to be positive definite (as in M3)
wfloor <- function(lk) {
  Hdd <- -lk$hdd; Hde <- -lk$hde; Hee <- -lk$hee
  tr <- Hdd + Hee; disc <- sqrt(pmax((Hdd - Hee)^2 / 4 + Hde^2, 0))
  e1 <- tr / 2 + disc; e2 <- tr / 2 - disc
  v1a <- ifelse(abs(Hde) > 1e-12, e1 - Hee, ifelse(Hdd >= Hee, 1, 0))
  v1b <- ifelse(abs(Hde) > 1e-12, Hde, ifelse(Hdd >= Hee, 0, 1))
  nr <- sqrt(v1a^2 + v1b^2); v1a <- v1a / nr; v1b <- v1b / nr
  f1 <- pmax(e1, 1e-6); f2 <- pmax(e2, 1e-6)
  list(dd = f1 * v1a^2 + f2 * v1b^2, ee = f1 * v1b^2 + f2 * v1a^2, de = (f1 - f2) * v1a * v1b)
}
# latent vector (w_1, [w_2], zeta), all N(0, I); L lower Cholesky factor of the correlation
lap_gp <- function(q, sp, S, L, v0, maxit = 100, tol = 1e-9) {
  n <- length(S$y); K <- sp$K
  bw <- function(k) (k - 1) * n + seq_len(n); iz <- K * n + seq_len(n)
  md <- as.vector(S$X %*% q$c); me <- as.vector(S$X %*% q$b)
  fobj <- function(v) {
    s1 <- as.vector(L %*% v[bw(1)]); se <- if (sp$ke == 1) s1 else as.vector(L %*% v[bw(sp$ke)])
    d <- md + q$alpha * s1; e <- me + q$gamma * se + q$sz * v[iz]
    lk <- odispers:::odisp_zil_lik(d, e, S)
    list(f = if (is.finite(lk$ll)) lk$ll - 0.5 * sum(v^2) else -Inf, lk = lk, d = d, e = e)
  }
  grad <- function(v, lk) {
    g <- -v
    g[bw(1)] <- g[bw(1)] + q$alpha * as.vector(crossprod(L, lk$gd))
    g[bw(sp$ke)] <- g[bw(sp$ke)] + q$gamma * as.vector(crossprod(L, lk$ge))
    g[iz] <- g[iz] + q$sz * lk$ge
    g
  }
  hess <- function(W) {                              # I + J' W J
    H <- diag((K + 1) * n); a <- q$alpha; gm <- q$gamma; sz <- q$sz
    tL <- t(L)
    if (sp$ke == 1) {
      H[bw(1), bw(1)] <- H[bw(1), bw(1)] + crossprod(L, (a^2 * W$dd + 2 * a * gm * W$de + gm^2 * W$ee) * L)
      B <- sweep(tL, 2, sz * (a * W$de + gm * W$ee), "*")
      H[bw(1), iz] <- B; H[iz, bw(1)] <- t(B)
    } else {
      H[bw(1), bw(1)] <- H[bw(1), bw(1)] + a^2 * crossprod(L, W$dd * L)
      H[bw(2), bw(2)] <- H[bw(2), bw(2)] + gm^2 * crossprod(L, W$ee * L)
      B12 <- a * gm * crossprod(L, W$de * L); H[bw(1), bw(2)] <- B12; H[bw(2), bw(1)] <- t(B12)
      B1 <- sweep(tL, 2, a * sz * W$de, "*"); H[bw(1), iz] <- B1; H[iz, bw(1)] <- t(B1)
      B2 <- sweep(tL, 2, gm * sz * W$ee, "*"); H[bw(2), iz] <- B2; H[iz, bw(2)] <- t(B2)
    }
    H[cbind(iz, iz)] <- H[cbind(iz, iz)] + sz^2 * W$ee
    H
  }
  v <- v0; fx <- fobj(v)
  if (!is.finite(fx$f)) { v <- 0 * v0; fx <- fobj(v) }
  if (!is.finite(fx$f)) return(list(logLik = -Inf, logjoint = -Inf))
  for (it in seq_len(maxit)) {
    R <- chol(hess(wfloor(fx$lk))); g <- grad(v, fx$lk)
    dl <- backsolve(R, backsolve(R, g, transpose = TRUE)); slope <- sum(g * dl)
    if (!is.finite(slope) || slope < tol) break
    step <- 1
    repeat { fs <- fobj(v + step * dl)
      if (is.finite(fs$f) && fs$f >= fx$f + 1e-4 * step * slope) break
      step <- step / 2; if (step < 1e-10) break }
    if (!(is.finite(fs$f) && fs$f > fx$f)) break
    v <- v + step * dl; fx <- fs
  }
  R <- chol(hess(wfloor(fx$lk)))
  list(logLik = fx$f - sum(log(diag(R))), logjoint = fx$f, v = v, R = R, d = fx$d, e = fx$e, iterations = it)
}

# ---- cross-validation helpers (as in the M3 script) ---------------------------------------
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

# ---- maximum likelihood fit ---------------------------------------------------------------
fit_gp <- function(m, th0, ck_file, maxit = 200, refine = 6) {
  sp <- spec[[m]]; n <- nv_; nlat <- (sp$K + 1) * n
  st <- new.env(); st$x <- rep(0, nlat); st$best <- -Inf; st$best_th <- th0; st$best_x <- NULL
  st$stage <- "optimise"; st$refine_th <- NULL; st$cache <- list(); st$n_eval <- 0
  if (file.exists(ck_file)) {
    ck <- readRDS(ck_file); th0 <- ck$theta; st$x <- ck$x; st$best_x <- ck$x; st$stage <- ck$stage
    st$cache <- ck$cache; if (!is.null(ck$refine_th)) st$refine_th <- ck$refine_th
    if (st$stage == "refine" && !is.null(ck$refine_th)) th0 <- ck$refine_th
    log_msg("resuming %s from checkpoint (stage %s), logLik %.3f, %d values cached", m, st$stage, ck$logLik, length(st$cache))
  }
  save_ck <- function() if (is.finite(st$best))
    saveRDS(list(theta = st$best_th, x = st$best_x, logLik = st$best, stage = st$stage,
                 refine_th = st$refine_th, cache = st$cache), ck_file)
  Lcache <- new.env()
  getL <- function(r) { k <- sprintf("%.10f", r)
    if (is.null(Lcache[[k]])) { if (length(ls(Lcache)) > 20) rm(list = ls(Lcache), envir = Lcache)
      Lcache[[k]] <- t(chol(corr(Dm, r) + diag(1e-8, n))) }
    Lcache[[k]] }
  lap <- function(th, branch = NULL, zero = st$n_eval == 0) {
    q <- unpack(th, sp); L <- getL(q$r)
    starts <- if (!is.null(branch)) list(branch) else {
      s <- list(st$x); if (!is.null(st$best_x)) s <- c(s, list(st$best_x)); if (zero) s <- c(s, list(rep(0, nlat))); s }
    starts <- starts[!duplicated(lapply(starts, function(z) round(z, 10)))]
    r <- list(logLik = -Inf, logjoint = -Inf)
    for (x0 in starts) { rk <- tryCatch(lap_gp(q, sp, S_all, L, x0), error = function(e) list(logLik = -Inf, logjoint = -Inf))
      if (is.finite(rk$logjoint) && rk$logjoint > r$logjoint) r <- rk }
    st$n_eval <- st$n_eval + 1
    if (is.finite(r$logLik)) {
      if (is.null(branch)) st$x <- r$v
      if (r$logLik > st$best) { st$best <- r$logLik; st$best_th <- th; st$best_x <- r$v; save_ck() }
    }
    if (st$n_eval %% 25 == 0) log_msg("  %d evaluations, best logLik %.3f", st$n_eval, st$best)
    r
  }
  key_of <- function(th) paste(sprintf("%.8f", th), collapse = ",")
  ll <- function(th, branch = NULL, tag = "") {
    key <- paste0(tag, key_of(th))
    if (!is.null(st$cache[[key]])) return(st$cache[[key]])
    val <- lap(th, branch)$logLik; st$cache[[key]] <- val
    if (length(st$cache) %% 20 == 0) save_ck()
    val
  }
  t0 <- Sys.time()
  if (st$stage == "optimise") {
    o <- stats::optim(th0, function(t) { f <- ll(t); if (is.finite(f)) -f else 1e10 }, method = "BFGS",
                      control = list(maxit = maxit, reltol = 1e-10, ndeps = rep(1e-3, length(th0))))
    th <- if (st$best > -o$value) st$best_th else o$par
    st$stage <- "refine"; save_ck()
    log_msg("%s optimiser: logLik %.3f after %d evaluations (convergence code %d)", m, st$best, st$n_eval, o$convergence)
  } else th <- th0
  hfd <- function(th, h = 0.01) {
    r0 <- lap(th, zero = TRUE); st$x <- r0$v; xc <- r0$v; f0 <- r0$logLik
    tag <- paste0("branch[", key_of(th), "]:"); lb <- function(t) ll(t, branch = xc, tag = tag)
    d <- length(th); H <- matrix(0, d, d); fp <- fm <- numeric(d)
    for (i in 1:d) { e <- replace(numeric(d), i, h); fp[i] <- lb(th + e); fm[i] <- lb(th - e)
                     H[i, i] <- (fp[i] - 2 * f0 + fm[i]) / h^2 }
    for (i in 1:(d - 1)) for (j in (i + 1):d) {
      ei <- replace(numeric(d), i, h); ej <- replace(numeric(d), j, h)
      H[i, j] <- H[j, i] <- (lb(th + ei + ej) - lb(th + ei - ej) - lb(th - ei + ej) + lb(th - ei - ej)) / (4 * h^2)
    }
    list(th = th, f0 = f0, g = (fp - fm) / (2 * h), H = H, r0 = r0)
  }
  hs <- NULL
  for (it in seq_len(refine)) {
    st$refine_th <- th; save_ck()
    hs <- hfd(th)
    e <- eigen(-hs$H, symmetric = TRUE)
    Bm <- e$vectors %*% diag(pmax(e$values, 1e-3 * max(e$values)), length(th)) %*% t(e$vectors)
    step <- as.vector(solve(Bm, hs$g)); if (sqrt(sum(step^2)) > 1) step <- step / sqrt(sum(step^2))
    gain <- sum(hs$g * step) - 0.5 * sum(step * (Bm %*% step))
    log_msg("%s refinement %d: logLik %.4f, predicted gain %.4f", m, it, hs$f0, gain)
    if (gain < 0.01) break
    tt <- 1; ok <- FALSE
    while (tt >= 1 / 16) { if (ll(th + tt * step) > hs$f0 + 1e-3) { ok <- TRUE; break }; tt <- tt / 2 }
    if (!ok) break
    th <- th + tt * step; hs <- NULL
  }
  if (is.null(hs)) { st$refine_th <- th; save_ck(); hs <- hfd(th) }
  fin <- hs$r0; names(th) <- nm_par
  V <- tryCatch(solve(-hs$H), error = function(e) matrix(NA, length(th), length(th)))
  dimnames(V) <- list(nm_par, nm_par)
  se <- sqrt(pmax(diag(V), 0)); z <- stats::qnorm(0.975); lo <- th - z * se; hi <- th + z * se
  logsc <- c(FALSE, FALSE, FALSE, FALSE, TRUE, sp$gamma_pos, TRUE, TRUE)
  tab <- data.frame(parameter = c(nm_par[1:4], "alpha", "gamma", "range_r (km)", "sigma_z"),
                    estimate = ifelse(logsc, exp(th), th), se = ifelse(logsc, NA, se),
                    lower = ifelse(logsc, exp(lo), lo), upper = ifelse(logsc, exp(hi), hi))
  k <- length(th)
  list(model = m, spec = sp, par = th, vcov = V, intervals = tab, logLik = fin$logLik, n_par = k,
       AIC = -2 * fin$logLik + 2 * k, gradient = hs$g, eigen = eigen(-hs$H, symmetric = TRUE, only.values = TRUE)$values,
       mode = fin$v, d = fin$d, e = fin$e, n_eval = st$n_eval,
       time = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}

# ---- cross-validation ---------------------------------------------------------------------
cv_gp <- function(fit, folds) {
  sp <- fit$spec; q <- unpack(fit$par, sp); n <- nv_
  Lf <- t(chol(corr(Dm, q$r) + diag(1e-8, n)))
  s_full <- lapply(seq_len(sp$K), function(k) as.vector(Lf %*% fit$mode[(k - 1) * n + seq_len(n)]))
  zeta_full <- fit$mode[sp$K * n + seq_len(n)]
  one <- function(tr, te) {
    ntr <- length(tr); C <- corr(Dm, q$r)
    Ltr <- t(chol(C[tr, tr] + diag(1e-8, ntr)))
    S_tr <- list(y = S_all$y[tr], n = S_all$n[tr], X = S_all$X[tr, , drop = FALSE])
    warm <- c(unlist(lapply(s_full, function(s) forwardsolve(Ltr, s[tr]))), zeta_full[tr])
    r <- NULL
    for (x0 in list(warm, 0 * warm)) { rk <- lap_gp(q, sp, S_tr, Ltr, x0)
      if (is.null(r) || rk$logjoint > r$logjoint) r <- rk }
    A <- forwardsolve(Ltr, C[tr, te, drop = FALSE])                 # n_tr x n_te
    cv_ <- pmax(1 - colSums(A^2), 0)                                # kriging variance of the unit process
    bw <- function(k) (k - 1) * ntr + seq_len(ntr); nl <- (sp$K + 1) * ntr
    Jd <- matrix(0, nl, length(te)); Je <- matrix(0, nl, length(te))
    Jd[bw(1), ] <- q$alpha * A; Je[bw(sp$ke), ] <- q$gamma * A
    Bd <- backsolve(r$R, Jd, transpose = TRUE); Be <- backsolve(r$R, Je, transpose = TRUE)
    w1 <- r$v[bw(1)]; we <- r$v[bw(sp$ke)]
    md <- as.vector(S_all$X[te, , drop = FALSE] %*% q$c) + q$alpha * as.vector(crossprod(A, w1))
    me <- as.vector(S_all$X[te, , drop = FALSE] %*% q$b) + q$gamma * as.vector(crossprod(A, we))
    Vd <- colSums(Bd^2) + q$alpha^2 * cv_
    Ve <- colSums(Be^2) + q$gamma^2 * cv_ + q$sz^2
    Cde <- colSums(Bd * Be) + (if (sp$ke == 1) q$alpha * q$gamma * cv_ else 0)
    lapply(seq_along(te), function(k) {
      Sig <- matrix(c(Vd[k], Cde[k], Cde[k], Ve[k]), 2)
      list(village = te[k], y = S_all$y[te[k]], n = S_all$n[te[k]], pmf = pmf_dE(md[k], me[k], Sig, S_all$n[te[k]]))
    })
  }
  ks <- sort(unique(folds)); if (test) ks <- ks[1:2]
  out <- unlist(lapply(ks, function(k) one(which(folds != k), which(folds == k))), recursive = FALSE)
  out[order(sapply(out, `[[`, "village"))]
}

# ---- prediction grid (20 km, as for M3) ----------------------------------------------------
grid <- st_make_grid(st_transform(sudan_borders, crs_km), cellsize = 20, what = "centers")
grid <- grid[lengths(st_intersects(grid, st_transform(st_union(sudan_borders), crs_km))) > 0]
grid <- st_sf(geometry = st_transform(grid, 4326))
grid$AI <- extract_aridity(grid, sudan_aridity, buffer_km = 10)
grid <- grid[is.finite(grid$AI) & grid$AI > 0, ]
xy_g <- st_coordinates(st_transform(st_geometry(grid), crs_km))[, 1:2]
X_g <- stats::model.matrix(~ log(AI), st_drop_geometry(grid))
predict_gp <- function(fit) {                       # plug-in: kriged mean of the processes, Z = 0
  sp <- fit$spec; q <- unpack(fit$par, sp); n <- nv_
  Rc <- chol(corr(Dm, q$r) + diag(1e-8, n))
  Dg <- sqrt(outer(xy_g[, 1], xy[, 1], "-")^2 + outer(xy_g[, 2], xy[, 2], "-")^2)
  Cg <- corr(Dg, q$r)
  sk <- function(k) as.vector(Cg %*% backsolve(Rc, fit$mode[(k - 1) * n + seq_len(n)]))
  s1 <- sk(1); se <- if (sp$ke == 1) s1 else sk(sp$ke)
  pi_ <- stats::plogis(as.vector(X_g %*% q$c) + q$alpha * s1)
  pz <- -expm1(-exp(as.vector(X_g %*% q$b) + q$gamma * se))
  ll_ <- st_coordinates(grid)
  data.frame(lon = ll_[, 1], lat = ll_[, 2], AI = grid$AI, S_pi = s1, S_lambda = se,
             pi = pi_, p_zone = pz, prevalence = pi_ * pz)
}

set.seed(2026)
folds <- sample(rep(1:10, length.out = nrow(v)))   # identical to the folds of M3

# ---- run --------------------------------------------------------------------------------
glm_start <- function() {
  c_ <- coef(stats::glm(I(y > 0) ~ log(AI), binomial, data = st_drop_geometry(v)))
  pos <- st_drop_geometry(v)[v$y > 0, ]
  b_ <- coef(stats::glm(cbind(y, n - y) ~ log(AI), binomial("cloglog"), data = pos))
  unname(c(b_, c_))
}
res_file <- function(m) file.path(out_dir, sprintf("gp_%s.rds", m))
for (m in run_models) {
  log_msg("==== %s", m)
  res <- if (file.exists(res_file(m))) readRDS(res_file(m)) else list()
  if (is.null(res$fit)) {
    th0 <- c(glm_start(), log(1), if (spec[[m]]$gamma_pos) log(0.5) else 0.5, log(300), log(0.5))
    if (m == "G2" && file.exists(res_file("G1"))) {             # start from G1 when available
      p1 <- readRDS(res_file("G1"))$fit$par; th0 <- p1; th0[6] <- log(max(abs(p1[6]), 0.1))
      log_msg("starting values from the fit of G1")
    }
    fit <- fit_gp(m, th0, file.path(out_dir, sprintf("gp_checkpoint_%s.rds", m)),
                  maxit = if (test) 2 else 200, refine = if (test) 0 else 6)
    print(fit$intervals, digits = 3)
    log_msg("%s fitted in %.1f min (%d evaluations), logLik %.3f, AIC %.2f; eigenvalues of the observed information %s",
            m, fit$time, fit$n_eval, fit$logLik, fit$AIC, paste(signif(range(fit$eigen), 3), collapse = " to "))
    res$fit <- fit; saveRDS(res, res_file(m))
  } else log_msg("%s: fit loaded (logLik %.3f)", m, res$fit$logLik)
  write.csv(res$fit$intervals, file.path(out_dir, sprintf("intervals_%s.csv", m)), row.names = FALSE)

  if (is.null(res$cv_scores)) {
    t0 <- Sys.time(); res$cv <- cv_gp(res$fit, folds)
    vil <- sapply(res$cv, `[[`, "village")
    res$cv_scores <- cv_scores(res$cv, v$AI[vil])
    log_msg("%s cross-validation in %.1f min", m, as.numeric(difftime(Sys.time(), t0, units = "mins")))
    saveRDS(res, res_file(m))
  }
  print(res$cv_scores$table, digits = 4)
  write.csv(res$cv_scores$table, file.path(out_dir, sprintf("cv_%s.csv", m)), row.names = FALSE)
  write.csv(res$cv_scores$pit, file.path(out_dir, sprintf("pit_%s.csv", m)), row.names = FALSE)

  if (is.null(res$grid)) { res$grid <- predict_gp(res$fit); saveRDS(res, res_file(m)) }
  log_msg("%s: results in %s (%.1f MB)", m, res_file(m), file.size(res_file(m)) / 2^20)
}
log_msg("done")
