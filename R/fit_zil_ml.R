#' Maximum likelihood for the zero-inflated river-network model
#'
#' Fits the zero-inflated model of \code{\link{odisp_fit_zil}} by maximum
#' likelihood. The likelihood of all the fixed parameters, the regression
#' coefficients of both parts and of the breeding intensity, \eqn{\sigma_u},
#' \eqn{\sigma_Z}, \eqn{a_E}, \eqn{\gamma} and the range of \eqn{u}, is the
#' Laplace approximation of the integral over the random effects only, the
#' breeding field \eqn{u} and the village effects \eqn{Z}:
#' \deqn{L(\beta, \psi) \approx \int p(y \mid \beta, u, Z)\, p(u \mid \psi)\, p(Z \mid \psi)\, du\, dZ.}
#' No prior is placed on any parameter. For given parameters the mode of
#' \eqn{(u, Z)} is found by the Newton method of \code{odisp_zil_laplace} with
#' the regression coefficients held fixed, and the likelihood is maximised over
#' all the parameters by a trust-region Newton method with finite-difference
#' derivatives and BFGS updates, followed by Newton refinement with the full
#' finite-difference Hessian. Wald intervals for every parameter come from the
#' observed information, the negative Hessian of the log-likelihood at the
#' maximum.
#'
#' The inner problem can have several local modes of \eqn{(u, Z)}. Each
#' evaluation starts the inner Newton from the mode at the previous evaluation
#' and from the mode at the best point so far, and the Laplace approximation is
#' taken at the mode with the higher log-posterior. The finite differences for
#' the Hessian all start from the mode at the point where the Hessian is taken,
#' so that they follow a single smooth branch.
#'
#' @inheritParams odisp_fit_zil
#' @param start Named list with optional starting values: \code{coefficients}
#'   (named as in an \code{odisp_zilfit}, matched by name, 0 when absent),
#'   \code{par} (log sigma_u, log sigma_z, a_E, gamma, log range), \code{u}
#'   and \code{z}.
#' @param refine Maximum number of Newton refinement steps.
#' @param intervals Compute the observed information and Wald intervals.
#' @return An object of class \code{c("odisp_zilml", "odisp_zilfit")}: the
#'   estimates \code{coefficients} and \code{hyper}, \code{logLik}, the
#'   parameter vector \code{par} on the optimisation scale, \code{vcov} (on
#'   that scale), the table \code{intervals}, the modes \code{u} and \code{z},
#'   and \code{setup} and \code{fld} for prediction and cross-validation (in
#'   \code{setup} the coefficients are held fixed at their estimates).
#' @export
odisp_fit_zil_ml <- function(formula, breeding, data, nodes, graph, pairs = NULL, phi = 40,
                             trials = "n", aridity = "AI", start = list(), trace = TRUE,
                             checkpoint = NULL, refine = 6, intervals = TRUE) {
  S <- odisp_u_setup(formula, breeding, data, nodes, graph, pairs, phi, trials, aridity, prior_sd = 1)
  S$m_c <- c(stats::qlogis(min(max(mean(S$y > 0), 0.05), 0.95)), rep(0, ncol(S$X) - 1))
  p <- ncol(S$X); q <- ncol(S$B); n <- length(S$y)
  nm_cf <- c(paste0("lambda.", colnames(S$X)), paste0("pi.", colnames(S$X)), paste0("theta.", S$nodes_names))
  nm_ps <- c("log_sigma_u", "log_sigma_z", "a_E", "gamma", "log_range")
  # starting values: coefficients by name (else the fit without u and Z), parameters
  cf0 <- c(S$m_beta, S$m_c, S$m_theta); names(cf0) <- nm_cf
  if (!is.null(start$coefficients)) {
    hit <- intersect(names(start$coefficients), nm_cf); cf0[hit] <- start$coefficients[hit]
    miss <- setdiff(grep("^theta\\.", nm_cf, value = TRUE), names(start$coefficients)); cf0[miss] <- 0
  }
  ps0 <- if (!is.null(start$par)) unname(start$par[1:5]) else c(log(5), log(0.6), 0.8, 0.2, log(800))
  names(ps0) <- nm_ps
  sc <- c(rep(0.4, length(cf0)), 1, 1, 0.4, 0.4, 1)          # optimisation scale: steps of 0.02 / 0.05
  th0 <- c(cf0, ps0) / sc
  idx_b <- seq_len(p); idx_c <- p + seq_len(p); idx_t <- 2 * p + seq_len(q); idx_p <- 2 * p + q + 1:5
  state <- new.env()
  state$x <- list(beta = cf0[idx_b], c = cf0[idx_c], theta = cf0[idx_t],
                  u = if (!is.null(start$u)) start$u else rep(0, graph$nv),
                  z = if (!is.null(start$z)) start$z else rep(0, n))
  state$best <- -Inf; state$n_eval <- 0; state$stage <- "optimise"; state$cache <- list(); state$refine_th <- NULL
  cache_file <- if (!is.null(checkpoint)) sub("\\.rds$", "_evaluations.rds", checkpoint) else NULL
  if (!is.null(checkpoint) && file.exists(checkpoint)) {
    ck <- readRDS(checkpoint)
    if (length(ck$theta) == length(th0)) { th0 <- ck$theta; state$x <- ck$x
      if (!is.null(ck$stage)) state$stage <- ck$stage
      # in the refinement, resume at the point where the last Hessian was started (its finite
      # differences are in the cache), not at the best perturbed point
      if (state$stage == "refine" && !is.null(ck$refine_th)) { th0 <- ck$refine_th; state$refine_th <- th0 }
      if (trace) cat(sprintf("  resuming from checkpoint (stage %s), logLik %.3f\n", state$stage, ck$logLik)) }
  }
  # likelihood values already computed (e.g. by a job that was stopped), keyed by parameter value
  if (!is.null(cache_file) && file.exists(cache_file)) {
    state$cache <- readRDS(cache_file)
    if (trace) cat(sprintf("  %d likelihood values loaded from %s\n", length(state$cache), cache_file))
  }
  save_ck <- function() if (!is.null(checkpoint) && is.finite(state$best))
    saveRDS(list(theta = state$best_th, x = state$best_x, logLik = state$best, stage = state$stage,
                 refine_th = state$refine_th), checkpoint)
  make_fld <- function(ps) {
    R <- wm_precision(graph, exp(ps[5]), 1)
    list(field = "u", R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE),
         sig_u = exp(ps[1]), sig_z = exp(ps[2]), aE = ps[3], gamma = ps[4], range_u = exp(ps[5]))
  }
  # The objective of the inner problem (u, Z given the parameters) can have several local
  # modes. The Laplace approximation is taken at the mode with the highest log-posterior
  # among those reached from the starting points: the mode at the last evaluation and the
  # mode at the best point so far. With `branch` given, the inner Newton starts from that
  # mode only, so that nearby evaluations (finite differences) follow one smooth branch.
  lap <- function(th, branch = NULL) {
    v <- th * sc; cf <- v[seq_along(cf0)]; ps <- v[idx_p]
    S2 <- S; S2$m_beta <- cf[idx_b]; S2$m_c <- cf[idx_c]; S2$m_theta <- cf[idx_t]
    S2$s2 <- 1e-10                                           # coefficients held fixed
    fld <- make_fld(ps)
    starts <- if (!is.null(branch)) list(branch) else {
      s <- list(state$x); if (!is.null(state$best_x) && !identical(state$best_x$u, state$x$u)) s <- c(s, list(state$best_x)); s }
    r <- list(logLik = -Inf, logpost = -Inf)
    for (x in starts) {
      x$beta <- cf[idx_b]; x$c <- cf[idx_c]; x$theta <- cf[idx_t]
      rk <- tryCatch(odisp_zil_laplace(fld, S2, x), error = function(e) list(logLik = -Inf, logpost = -Inf))
      if (is.finite(rk$logpost) && rk$logpost > r$logpost) r <- rk
    }
    state$n_eval <- state$n_eval + 1
    if (is.finite(r$logLik)) {
      if (is.null(branch)) state$x <- r$x
      if (r$logLik > state$best) {
        state$best <- r$logLik; state$best_th <- th; state$best_x <- r$x
        save_ck()
      }
    }
    if (trace && state$n_eval %% 10 == 0) cat(sprintf("  %d evaluations, best logLik %.3f\n", state$n_eval, state$best))
    r$fld <- fld; r$S2 <- S2
    r
  }
  key_of <- function(th) paste(sprintf("%.8f", th), collapse = ",")
  ll <- function(th, branch = NULL, tag = "") {
    key <- paste0(tag, key_of(th))
    if (!is.null(state$cache[[key]])) return(state$cache[[key]])
    val <- lap(th, branch)$logLik
    state$cache[[key]] <- val
    if (!is.null(cache_file) && length(state$cache) %% 5 == 0) saveRDS(state$cache, cache_file)
    val
  }
  t0 <- Sys.time()
  if (state$stage == "optimise") {
    th <- odisp_newton_max(ll, th0, bfgs = TRUE, maxit = 50)
    if (state$best > ll(th)) th <- state$best_th
    state$stage <- "refine"; save_ck()
  } else th <- th0
  if (trace) cat(sprintf("  optimiser: logLik %.3f after %d evaluations\n", state$best, state$n_eval))
  # Newton refinement and observed information (central differences, step 0.05 on the optimisation scale)
  # mode of (u, Z) at th: the best of the starting points; all the finite differences around
  # th start from it (one branch of the inner mode)
  centre <- function(th) {
    r0 <- lap(th); state$x <- r0$x
    list(x = r0$x, f0 = r0$logLik, tag = paste0("branch[", key_of(th), "]:"))
  }
  hess <- function(th, h = 0.05) {
    ce <- centre(th); f0 <- ce$f0
    lb <- function(t) ll(t, branch = ce$x, tag = ce$tag)
    d <- length(th); H <- matrix(0, d, d); fp <- fm <- numeric(d)
    for (i in 1:d) { e <- replace(numeric(d), i, h); fp[i] <- lb(th + e); fm[i] <- lb(th - e)
                     H[i, i] <- (fp[i] - 2 * f0 + fm[i]) / h^2 }
    for (i in 1:(d - 1)) for (j in (i + 1):d) {
      ei <- replace(numeric(d), i, h); ej <- replace(numeric(d), j, h)
      H[i, j] <- H[j, i] <- (lb(th + ei + ej) - lb(th + ei - ej) - lb(th - ei + ej) + lb(th - ei - ej)) / (4 * h^2)
    }
    list(th = th, f0 = f0, g = (fp - fm) / (2 * h), H = H, x = ce$x)
  }
  hs <- NULL
  for (it in seq_len(refine)) {
    if (!intervals && it > 1) break
    state$refine_th <- th; save_ck()
    hs <- hess(th)
    e <- eigen(-hs$H, symmetric = TRUE)
    Bm <- e$vectors %*% diag(pmax(e$values, 1e-3 * max(e$values)), length(th)) %*% t(e$vectors)
    step <- as.vector(solve(Bm, hs$g))
    if (sqrt(sum(step^2)) > 1) step <- step / sqrt(sum(step^2))
    gain <- sum(hs$g * step) - 0.5 * sum(step * (Bm %*% step))
    if (trace) cat(sprintf("  refinement %d: logLik %.3f, predicted gain %.3f\n", it, hs$f0, gain))
    if (gain < 0.01) break
    t <- 1; ok <- FALSE
    while (t >= 1 / 16) { if (ll(th + t * step) > hs$f0 + 1e-3) { ok <- TRUE; break }; t <- t / 2 }
    if (!ok) break
    th <- th + t * step
  }
  if (!is.null(cache_file)) saveRDS(state$cache, cache_file)
  if (is.null(hs) || !identical(hs$th, th)) hs <- hess(th)
  fin <- lap(th, branch = hs$x)                             # the mode at which the Hessian was taken
  v <- th * sc; cf <- v[seq_along(cf0)]; names(cf) <- nm_cf; ps <- v[idx_p]; names(ps) <- nm_ps
  # observed information on the original scale of (coefficients, log sigma_u, log sigma_z, a_E, gamma, log range)
  Hn <- hs$H / outer(sc, sc)
  ev <- eigen(-Hn, symmetric = TRUE, only.values = TRUE)$values
  V <- tryCatch(solve(-Hn), error = function(e) matrix(NA, length(v), length(v)))
  dimnames(V) <- list(c(nm_cf, nm_ps), c(nm_cf, nm_ps))
  se <- sqrt(pmax(diag(V), 0)); z <- stats::qnorm(0.975)
  est <- c(cf, ps); lo <- est - z * se; hi <- est + z * se
  logsc <- c(rep(FALSE, length(cf)), TRUE, TRUE, FALSE, FALSE, TRUE)
  tab <- data.frame(parameter = c(nm_cf, "sigma_u", "sigma_z", "a_E", "gamma", "range_u (km)"),
                    estimate = ifelse(logsc, exp(est), est), se = ifelse(logsc, NA, se),
                    lower = ifelse(logsc, exp(lo), lo), upper = ifelse(logsc, exp(hi), hi))
  rownames(tab) <- NULL
  fld <- fin$fld
  structure(list(coefficients = cf, vcov = V, intervals = tab, eigen = ev, gradient = hs$g / sc,
                 hyper = c(sigma_u = exp(ps[[1]]), sigma_z = exp(ps[[2]]), a_E = ps[[3]], gamma = ps[[4]],
                           range_u = exp(ps[[5]])),
                 field_type = if (is.null(graph$type)) "wm" else graph$type, par = ps, logLik = fin$logLik,
                 u = fin$x$u, z = fin$x$z, E = exp(fin$a), pi = stats::plogis(fin$d), p = -expm1(-exp(fin$e)),
                 setup = fin$S2, fld = fld, n_eval = state$n_eval,
                 time = as.numeric(difftime(Sys.time(), t0, units = "mins"))),
            class = c("odisp_zilml", "odisp_zilfit"))
}

#' @export
print.odisp_zilml <- function(x, ...) {
  cat("Zero-inflated river-network model, maximum likelihood (Laplace approximation over u and Z)\n")
  print(x$intervals, digits = 3)
  cat("log-likelihood:", round(x$logLik, 3), "\n")
  invisible(x)
}
