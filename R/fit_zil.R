#' Zero-inflated river-network model with regressions in both parts
#'
#' A village lies in a zone of local transmission with probability
#' \eqn{\pi(x)}; outside, there are no positives, and inside, the number of
#' positives is binomial with prevalence \eqn{p(x) = 1 - e^{-\lambda(x)}}. By default both
#' parts are regressions on the logarithm of the exposure from the river
#' network and on the covariates of \code{formula} (for example log aridity),
#' each with its own coefficients:
#' \deqn{\mathrm{logit}\,\pi(x) = d(x)^\top a + a_E \log E(x), \qquad
#' \log\lambda(x) = d(x)^\top b + \gamma \log E(x) + Z(x),}
#' where \eqn{E(x) = \int_\Gamma K_\phi(\|x-s\|) \exp\{z(s)^\top\theta + u(s)\}
#' ds}, \eqn{u} is a Whittle-Matern field on the river network and \eqn{Z} an
#' independent village effect. Either regression can instead use raw exposure.
#' At zero exposure the active-village risk is zero by definition. A raw-exposure
#' activity logit then equals its village covariate predictor.
#'
#' The coefficients \eqn{a}, \eqn{b} and \eqn{\theta}, the field and the
#' village effects are integrated out by a Laplace approximation, the
#' coefficients with diffuse Gaussian priors. For fixed \eqn{a_E} and
#' \eqn{\gamma} the two predictors of each village are linear in the latent
#' vector apart from \eqn{\log E}, so the algorithm of \code{\link{odisp_fit_u}}
#' applies with two predictors per village; the Hessian used in the
#' approximation is the observed information of the two predictors with
#' eigenvalues bounded away from zero. \eqn{\sigma_u}, \eqn{\sigma_Z},
#' \eqn{a_E} and \eqn{\gamma} maximise the approximate marginal likelihood.
#'
#' @inheritParams odisp_fit_zi
#' @param start Optional starting values \code{c(log sigma_u, log sigma_z,
#'   a_E, gamma)}.
#' @param fixed Optional named list fixing some of \code{sigma_u},
#'   \code{sigma_z}, \code{a_E}, \code{gamma}.
#' @param range_up Range of the tail-up component when the graph combines
#'   tail-up and tail-down fields. Used as a starting value when
#'   \code{estimate_range = TRUE}, otherwise held fixed.
#' @param occurrence_scale Use log exposure (default) or raw, linear exposure in
#'   the activity logit.
#' @param risk_scale Use log exposure (default) or raw, linear exposure in the
#'   log hazard within an active village. At exactly zero exposure the risk
#'   is defined to be zero. With linear exposure its right-hand limit can differ.
#' @return An object of class \code{odisp_zilfit}.
#' @export
odisp_fit_zil <- function(formula, breeding, data, nodes, graph, pairs = NULL, phi = 40,
                          range_u = 800, trials = "n", aridity = "AI", prior_sd = 10,
                          start = NULL, fixed = list(), trace = TRUE, estimate_range = FALSE,
                          checkpoint = NULL, bfgs = FALSE, range_up = range_u, x_start = NULL,
                          multistart = TRUE, occurrence_scale = c("log","linear"),
                          risk_scale = c("log","linear")) {
  occurrence_scale <- match.arg(occurrence_scale)
  risk_scale <- match.arg(risk_scale)
  S <- odisp_u_setup(formula, breeding, data, nodes, graph, pairs, phi, trials, aridity, prior_sd)
  n <- length(S$y)
  S$m_c <- c(stats::qlogis(min(max(mean(S$y > 0), 0.05), 0.95)), rep(0, ncol(S$X) - 1))
  # with a sum of tail-up and tail-down fields, the share of variance of the
  # tail-up field, phi = plogis(mix), is a further parameter
  has_mix <- identical(graph$type, "tailup+taildown")
  R <- wm_precision(graph, range_u, 1)
  base <- list(field = "u", R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE))
  if (is.null(start)) start <- c(log(3), log(0.6), 1, 0.3)
  if (has_mix && length(start) == 4) start <- c(start, 0)
  nm_start <- c("sigma_u", "sigma_z", "a_E", "gamma", if (has_mix) "mix")
  start <- start[seq_along(nm_start)]; names(start) <- nm_start
  if (estimate_range) start <- c(start, range_u = log(range_u))
  # tail-up + tail-down with estimated ranges: a separate range for the
  # tail-up part, started from range_up
  if (estimate_range && has_mix) start <- c(start, range_up = log(range_up))
  tr <- list(sigma_u = log, sigma_z = log, a_E = identity, gamma = identity, mix = stats::qlogis,
             range_u = log, range_up = log)
  free <- setdiff(names(start), names(fixed))
  x0 <- list(beta = S$m_beta, c = S$m_c, theta = S$m_theta, u = rep(0, graph$nv), z = rep(0, n))
  if (!is.null(x_start)) x0 <- x_start                      # latent starting point, e.g. a neighbouring fit
  # resume from a checkpoint: best parameters and latent mode found so far
  if (!is.null(checkpoint) && file.exists(checkpoint)) {
    ck <- readRDS(checkpoint)
    ck_scale <- if (is.null(ck$occurrence_scale)) "log" else ck$occurrence_scale
    if (!identical(ck_scale,occurrence_scale)) stop("checkpoint has a different occurrence scale")
    ck_risk <- if (is.null(ck$risk_scale)) "log" else ck$risk_scale
    if (!identical(ck_risk,risk_scale)) stop("checkpoint has a different risk scale")
    if (all(names(ck$par) %in% names(start)) && length(ck$x$u) == graph$nv) {
      start[names(ck$par)] <- ck$par; x0 <- ck$x
      if (trace) cat(sprintf("  resuming from checkpoint, logLik %.3f\n", ck$logLik))
    }
  }
  make_fld <- function(par) {
    full <- start; full[free] <- par
    for (nm in names(fixed)) full[nm] <- tr[[nm]](fixed[[nm]])
    fb <- base
    r_u <- if (estimate_range) exp(full[["range_u"]]) else range_u
    if (has_mix || estimate_range) {
      mix <- if (has_mix) stats::plogis(full[["mix"]]) else 0.5
      r_up <- if (has_mix) {
        if (estimate_range) exp(full[["range_up"]]) else range_up
      } else r_u
      fb$R <- wm_precision(graph, r_u, 1, mix = mix, range_up = r_up)
      fb$range_up <- r_up
      fb$Rch <- Matrix::update(base$Rch, fb$R)
      if (has_mix) fb$mix <- mix
    }
    c(fb, list(sig_u = exp(full[["sigma_u"]]), sig_z = exp(full[["sigma_z"]]),
               aE = full[["a_E"]], gamma = full[["gamma"]], range_u = r_u,
               occurrence_scale=occurrence_scale, risk_scale=risk_scale))
  }
  state <- new.env(); state$x <- x0; state$cache <- list(); state$best <- NULL
  eval_par <- function(par) {
    key <- paste(sprintf("%.17g", par), collapse = ",")
    if (!is.null(state$cache[[key]])) {
      if (!is.null(state$best) && identical(state$best$key, key)) return(state$best$res)
      return(state$cache[[key]])
    }
    fld <- make_fld(par)
    lap <- function(x) tryCatch(odisp_zil_laplace(fld, S, x), error = function(e) list(logLik = -Inf, logpost = -Inf))
    res <- lap(state$x)
    # the posterior of the latent effects can have several modes: also start
    # from the mode of the best point so far and keep the higher mode
    if (multistart && !is.null(state$best) && !identical(state$best$res$x, state$x)) {
      r2 <- lap(state$best$res$x)
      if (isTRUE(r2$logpost > res$logpost + 1e-6)) res <- r2
    }
    if (is.finite(res$logLik)) state$x <- res$x
    if (is.finite(res$logLik) && (is.null(state$best) || res$logLik > state$best$res$logLik)) {
      state$best <- list(key = key, par = par, res = res)
      if (!is.null(checkpoint)) saveRDS(list(par = par, x = res$x, logLik = res$logLik,
                                            occurrence_scale=occurrence_scale,risk_scale=risk_scale), checkpoint)
    }
    state$cache[[key]] <- list(logLik = res$logLik)          # latent modes are kept only for the best point
    if (trace) cat(sprintf("  sigma_u %.3f  sigma_z %.3f  a_E %.3f  gamma %.3f%s%s%s  logLik = %.3f\n",
                           fld$sig_u, fld$sig_z, fld$aE, fld$gamma,
                           if (has_mix) sprintf("  phi %.3f", fld$mix) else "",
                           if (estimate_range) sprintf("  range %.0f", fld$range_u) else "",
                           if (estimate_range && has_mix) sprintf("  range_up %.0f", fld$range_up) else "", res$logLik))
    if (is.null(res$x)) res$x <- state$x
    res
  }
  par <- start[free]
  if (length(par)) {
    # A unit change in a raw-exposure slope can be enormous on some fitted rivers.
    # Rescale only the optimizer coordinates; returned slopes retain their original units.
    units <- stats::setNames(rep(1,length(par)),names(par))
    raw_slopes <- intersect(c(if (occurrence_scale=="linear") "a_E",
                              if (risk_scale=="linear") "gamma"),free)
    if (length(raw_slopes)) {
      log_E <- odisp_zil_eval(state$x,S,make_fld(par))$a
      units[raw_slopes] <- pmin(1,pmax(abs(par[raw_slopes]),exp(-max(log_E)),1e-12))
    }
    par <- units * odisp_newton_max(function(p) eval_par(p*units)$logLik,par/units,bfgs=bfgs)
  }
  res <- eval_par(par)
  if (is.null(res$Vfix)) {                                  # make sure the returned point has its mode
    par <- state$best$par
    res <- state$best$res
  }
  fld <- make_fld(par)
  range_u <- fld$range_u
  if (trace) cat(sprintf("  converged: logLik = %.3f, %d evaluations\n", res$logLik, length(state$cache)))
  p <- ncol(S$X)
  nm <- c(paste0("lambda.", colnames(S$X)), paste0("pi.", colnames(S$X)),
          if (length(S$nodes_names)) paste0("theta.", S$nodes_names))
  est <- c(res$x$beta, res$x$c, res$x$theta); names(est) <- nm
  dimnames(res$Vfix) <- list(nm, nm)
  structure(list(coefficients = est, vcov = res$Vfix,
                 hyper = c(sigma_u = fld$sig_u, sigma_z = fld$sig_z, a_E = fld$aE, gamma = fld$gamma,
                           range_u = range_u, if (has_mix) c(phi = fld$mix, range_up = fld$range_up)),
                 field_type = if (is.null(graph$type)) "wm" else graph$type,
                 par = par, fixed = fixed, logLik = res$logLik, u = res$x$u, z = res$x$z,
                 E = exp(res$a), pi = stats::plogis(res$d), p = -expm1(-exp(res$e)),
                 setup = S, fld = fld),
            class = "odisp_zilfit")
}

# log-likelihood and derivatives with respect to d = logit(pi) and e = log(lambda)
odisp_zil_lik <- function(d, e, S) {
  n <- S$n; y <- S$y
  if (any(!is.finite(d)) || any(!is.finite(e))) return(list(ll=-Inf))
  lp <- stats::plogis(d, log.p=TRUE)
  l1p <- stats::plogis(d, lower.tail=FALSE, log.p=TRUE)
  pi_ <- exp(lp)
  lam <- exp(e)
  ld <- le <- ldd <- lde <- lee <- ll <- numeric(length(y))
  pos <- y > 0
  if (any(pos)) {
    l <- lam[pos]; nn <- n[pos]; yy <- y[pos]
    logp <- log1mexp(l)
    tiny <- e[pos] < log(.Machine$double.eps)
    logp[tiny] <- e[pos][tiny]
    failure <- numeric(length(l))
    has_failure <- nn>yy
    failure[has_failure] <- (nn-yy)[has_failure]*l[has_failure]
    ll[pos] <- lp[pos]+yy*logp-failure+lchoose(nn,yy)
    if (any(!is.finite(ll[pos]))) return(list(ll=-Inf))
    score <- lgcp_cloglog_score(l,yy,nn)
    ld[pos] <- exp(l1p[pos]); ldd[pos] <- -pi_[pos]*ld[pos]
    le[pos] <- score$score; lee[pos] <- -score$curvature
  }
  zr <- !pos & n>0
  if (any(zr)) {
    # Given a zero, compute both posterior state probabilities on the log scale.
    log_k <- log(n[zr])+e[zr]
    k <- exp(log_k)
    inactive <- l1p[zr]
    active <- lp[zr]-k
    lf <- pmax(inactive,active)+log1p(exp(-abs(inactive-active)))
    lr <- active-lf; ls <- inactive-lf
    r <- exp(lr); s <- exp(ls)
    rk <- exp(lr+log_k)
    ll[zr] <- lf
    ld[zr] <- r*exp(l1p[zr])-s*pi_[zr]
    le[zr] <- -rk
    ldd[zr] <- r*s-pi_[zr]*exp(l1p[zr])
    lde[zr] <- -rk*s
    lee[zr] <- -rk+exp(lr+ls+2*log_k)
  }
  tot <- sum(ll) + 0
  list(ll = if (is.finite(tot)) tot else -Inf, gd = ld, ge = le, hdd = ldd, hde = lde, hee = lee)
}

odisp_zil_vec <- function(x) c(x$beta, x$c, x$theta, x$u, x$z)
odisp_zil_unvec <- function(v, x) {
  p <- length(x$beta); pc <- length(x$c); q <- length(x$theta); nu <- length(x$u)
  list(beta = v[seq_len(p)], c = v[p + seq_len(pc)], theta = v[p + pc + seq_len(q)],
       u = v[p + pc + q + seq_len(nu)], z = v[p + pc + q + nu + seq_along(x$z)])
}

odisp_zil_eval <- function(x, S, fld) {
  ev <- odisp_zi_eval(list(beta = x$beta, theta = x$theta, u = x$u, z = x$z), S)
  lin_c <- as.vector(S$X %*% x$c)
  gate_exposure <- if (identical(fld$occurrence_scale,"linear")) exp(ev$a) else ev$a
  risk_exposure <- if (identical(fld$risk_scale,"linear")) exp(ev$a) else ev$a
  list(a = ev$a, V = ev$V, d = lin_c + fld$aE * gate_exposure, e = ev$b + fld$gamma * risk_exposure)
}

odisp_zil_prior <- function(x, S, fld) {
  (sum((x$beta - S$m_beta)^2) + sum((x$c - S$m_c)^2) + sum((x$theta - S$m_theta)^2)) / S$s2 +
    sum(x$z^2) / fld$sig_z^2 + as.numeric(Matrix::crossprod(x$u, fld$R %*% x$u)) / fld$sig_u^2
}

odisp_zil_ops <- function(x, ev, lk, S, fld) {
  n <- length(S$y); p <- ncol(S$X); q <- ncol(S$B); nu <- length(x$u)
  ib <- seq_len(p); ic <- p + seq_len(p); it <- 2 * p + seq_len(q)
  iu <- 2 * p + q + seq_len(nu); iz <- 2 * p + q + nu + seq_len(n)
  linear_gate <- identical(fld$occurrence_scale,"linear")
  linear_risk <- identical(fld$risk_scale,"linear")
  # Derivatives of the two predictors with respect to log E.
  aE <- if (linear_gate) fld$aE * exp(ev$a) else fld$aE
  gm <- if (linear_risk) fld$gamma * exp(ev$a) else fld$gamma
  Jt <- as.matrix(ev$V %*% S$B)
  Ju <- ev$V %*% S$A
  LW <- Matrix::solve(fld$Rch, Matrix::solve(fld$Rch, Matrix::t(Ju), system = "P"), system = "L")
  Maa <- S$s2 * tcrossprod(Jt) + fld$sig_u^2 * as.matrix(Matrix::crossprod(LW))
  XX <- S$s2 * tcrossprod(S$X)
  Mdd <- XX + (if(linear_gate) tcrossprod(aE) else aE^2) * Maa
  Mee <- XX + (if(linear_risk) tcrossprod(gm) else gm^2) * Maa + diag(fld$sig_z^2, n)
  Mde <- if (linear_risk) outer(rep_len(aE,n),gm)*Maa else aE*gm*Maa
  # observed information of (d, e) per village, floored to be positive definite
  Hdd <- -lk$hdd; Hde <- -lk$hde; Hee <- -lk$hee
  tr <- Hdd + Hee; disc <- sqrt(pmax((Hdd - Hee)^2 / 4 + Hde^2, 0))
  e1 <- tr / 2 + disc; e2 <- tr / 2 - disc
  fl <- 1e-6
  angle <- atan2(2*Hde,Hdd-Hee)/2
  v1a <- cos(angle); v1b <- sin(angle)
  f1 <- pmax(e1, fl); f2 <- pmax(e2, fl)
  Wdd <- f1 * v1a^2 + f2 * v1b^2; Wee <- f1 * v1b^2 + f2 * v1a^2; Wde <- (f1 - f2) * v1a * v1b
  l11 <- sqrt(Wdd); l21 <- Wde / l11; l22 <- sqrt(pmax(Wee - l21^2, fl))
  # C = I + L' M L, M = [[Mdd, Mde], [Mde', Mee]], L = per village [[l11, 0], [l21, l22]]
  C11 <- outer(l11, l11) * Mdd + outer(l11, l21) * Mde + outer(l21, l11) * t(Mde) + outer(l21, l21) * Mee
  C12 <- outer(l11, l22) * Mde + outer(l21, l22) * Mee
  C22 <- outer(l22, l22) * Mee
  C <- rbind(cbind(C11, C12), cbind(t(C12), C22)) + diag(2 * n)
  cC <- chol((C + t(C)) / 2)
  Ja <- function(v) as.vector(Jt %*% v[it] + Ju %*% v[iu])
  Jv <- function(v) { av <- Ja(v)
    list(d = as.vector(S$X %*% v[ic]) + aE * av, e = as.vector(S$X %*% v[ib]) + gm * av + v[iz]) }
  Jtw <- function(wd, we) { wa <- aE * wd + gm * we
    c(as.vector(crossprod(S$X, we)), as.vector(crossprod(S$X, wd)), as.vector(crossprod(Jt, wa)),
      as.vector(Matrix::crossprod(Ju, wa)), we) }
  Pinv <- function(r) c(S$s2 * r[ib], S$s2 * r[ic], S$s2 * r[it],
                        fld$sig_u^2 * as.vector(Matrix::solve(fld$Rch, r[iu])), fld$sig_z^2 * r[iz])
  Pv <- function(v) c(v[ib] / S$s2, v[ic] / S$s2, v[it] / S$s2, as.vector(fld$R %*% v[iu]) / fld$sig_u^2,
                      v[iz] / fld$sig_z^2)
  HGNinv <- function(r) {
    s0 <- Pinv(r); jv <- Jv(s0)
    t1 <- l11 * jv$d + l21 * jv$e; t2 <- l22 * jv$e
    w <- backsolve(cC, backsolve(cC, c(t1, t2), transpose = TRUE))
    w1 <- w[seq_len(n)]; w2 <- w[n + seq_len(n)]
    s0 - Pinv(Jtw(l11 * w1, l21 * w1 + l22 * w2))
  }
  ga <- aE * lk$gd + gm * lk$ge                                     # d loglik / d log E
  cc <- as.vector(Matrix::crossprod(ev$V, ga))
  Cv <- function(v) {
    nuv <- as.vector(S$B %*% v[it]) + as.vector(S$A %*% v[iu])
    cx <- cc * nuv - as.vector(Matrix::crossprod(ev$V, ga * as.vector(ev$V %*% nuv)))
    if (linear_gate || linear_risk) {
      # For either raw term exp(log E) also contributes curvature along log E.
      extra <- (if (linear_gate) lk$gd*aE else 0) + (if (linear_risk) lk$ge*gm else 0)
      cx <- cx + as.vector(Matrix::crossprod(ev$V, extra * as.vector(ev$V %*% nuv)))
    }
    out <- numeric(length(v))
    out[it] <- as.vector(crossprod(S$B, cx)); out[iu] <- as.vector(Matrix::crossprod(S$A, cx))
    out
  }
  Hv <- function(v) { jv <- Jv(v)
    Pv(v) + Jtw(Hdd * jv$d + Hde * jv$e, Hde * jv$d + Hee * jv$e) - Cv(v) }
  list(HGNinv = HGNinv, Hv = Hv, Jtw = Jtw, Pv = Pv, Jt = Jt, Ju = Ju, cC = cC,
       ib = ib, ic = ic, it = it)
}

odisp_zil_laplace <- function(fld, S, x, maxit = 50, tol = 1e-7) {
  m <- odisp_zil_vec(list(beta = S$m_beta, c = S$m_c, theta = S$m_theta, u = 0 * x$u, z = 0 * x$z))
  f_of <- function(x) {
    ev <- odisp_zil_eval(x, S, fld); lk <- odisp_zil_lik(ev$d, ev$e, S)
    v <- lk$ll - 0.5 * odisp_zil_prior(x, S, fld)
    list(f = if (is.finite(v)) v else -Inf, ev = ev, lk = lk)
  }
  fx <- f_of(x)
  if (!is.finite(fx$f)) {
    x <- list(beta = S$m_beta, c = S$m_c, theta = S$m_theta, u = 0 * x$u, z = 0 * x$z); fx <- f_of(x)
  }
  if (!is.finite(fx$f)) stop("non-finite objective at the starting point")
  for (it in seq_len(maxit)) {
    ops <- odisp_zil_ops(x, fx$ev, fx$lk, S, fld)
    xv <- odisp_zil_vec(x)
    gv <- ops$Jtw(fx$lk$gd, fx$lk$ge) - ops$Pv(xv - m)
    pc <- odisp_u_pcg(gv, ops)
    slope <- sum(gv * pc$d)
    if (!is.finite(slope) || slope <= 0) { pc$d <- ops$HGNinv(gv); slope <- sum(gv * pc$d) }
    if (slope < tol) break
    step <- 1
    repeat {
      xs <- odisp_zil_unvec(xv + step * pc$d, x); fs <- f_of(xs)
      if (fs$f >= fx$f + 1e-4 * step * slope || step < 1e-8) break
      step <- step / 2
    }
    if (fs$f < fx$f) break
    gain <- fs$f - fx$f
    x <- xs; fx <- fs
    if (gain < tol) break
  }
  ops <- odisp_zil_ops(x, fx$ev, fx$lk, S, fld)
  logdet <- 2 * sum(log(diag(ops$cC)))
  idx <- c(ops$ib, ops$ic, ops$it)
  Vf <- sapply(idx, function(k) { e <- numeric(length(odisp_zil_vec(x))); e[k] <- 1; ops$HGNinv(e)[idx] })
  list(logLik = fx$f - 0.5 * logdet, logpost = fx$f, x = x, a = fx$ev$a, d = fx$ev$d, e = fx$ev$e, iterations = it,
       Vfix = (Vf + t(Vf)) / 2)
}

#' @export
print.odisp_zilfit <- function(x, ...) {
  h <- x$hyper
  cat("Zero-inflated river-network model (logit for the zone of transmission, log for the force of infection)\n")
  term <- if (identical(x$fld$occurrence_scale,"linear")) "E" else "log E"
  risk_term <- if (identical(x$fld$risk_scale,"linear")) "E" else "log E"
  cat(sprintf("logit(pi): %.3f times %s; log(lambda): %.3f times %s\n",h[["a_E"]],term,h[["gamma"]],risk_term))
  cat(sprintf("Field u: range %g km, sd %.3f; village effect sd %.3f\n", h[["range_u"]], h[["sigma_u"]], h[["sigma_z"]]))
  print(round(cbind(estimate = x$coefficients, se = sqrt(diag(x$vcov))), 3))
  cat("Laplace log-likelihood:", round(x$logLik, 2), "\n")
  invisible(x)
}

#' Predict from the zero-inflated model with regressions in both parts
#'
#' Plug-in prediction at the posterior mode, with the village effect set to
#' zero: probability of a zone of transmission \code{pi}, prevalence within a
#' zone \code{p_zone} and expected prevalence \code{prevalence}. For a model
#' fitted with a migration kernel (\code{migration_range} in
#' \code{\link{odisp_fit_zil_ml}}), the exposure is the mixture
#' \eqn{(1 - w) E^{loc} + w E^{mig}}, and the local and migratory exposures and
#' the share of exposure due to migration, \code{s_mig}, are also returned.
#' \code{pairs} need only reach the local range \code{phi}: the migratory
#' exposure is computed from the breeding mass of the river nodes summed in
#' square blocks of side \code{migration_block} km, which is accurate because
#' the migration kernel varies little over a block.
#'
#' @param object An \code{odisp_zilfit} object.
#' @param migration_block Side (km) of the blocks used for the migratory exposure.
#' @inheritParams predict.odisp_zifit
#' @export
predict.odisp_zilfit <- function(object, newdata, nodes, graph, pairs = NULL, migration_block = 5, ...) {
  S <- object$setup; phi <- S$phi; h <- object$hyper
  mr <- object$migration_range
  wm <- if ("w" %in% names(h)) h[["w"]] else 0
  mig <- !is.null(mr) && wm > 0
  pairs <- odisp_prediction_pairs(newdata, nodes, pairs, phi)
  keep <- pairs$r < phi
  i <- pairs$i[keep]; j <- pairs$j[keep]; r <- pairs$r[keep]
  wend <- function(t) ifelse(t < 1, (1 - t)^4 * (1 + 4 * t), 0)
  Bn <- odisp_prediction_matrix(S$breeding, as.data.frame(nodes)[j, , drop = FALSE],
                                 S$breeding_spec, colnames(S$B), drop_intercept = TRUE)
  p <- ncol(S$X); q <- length(S$m_theta); cf <- object$coefficients
  w <- exp(as.vector(Bn %*% cf[2 * p + seq_len(q)]) + as.vector(graph$A[j, , drop = FALSE] %*% object$u))
  sum_by <- function(x) as.vector(tapply(c(x, rep(0, pairs$m)), c(i, seq_len(pairs$m)), sum))
  E_loc <- sum_by(nodes$len[j] * wend(r / phi) / (2 / 3 * phi) * w)
  if (mig) {
    # breeding mass of every node, summed in blocks (mass-weighted centroids)
    Ball <- odisp_prediction_matrix(S$breeding, as.data.frame(nodes), S$breeding_spec, colnames(S$B),
                                    drop_intercept = TRUE)
    mass <- nodes$len * exp(as.vector(Ball %*% cf[2 * p + seq_len(q)]) + as.vector(graph$A %*% object$u))
    blk <- paste(floor(nodes$x / migration_block), floor(nodes$y / migration_block))
    M <- as.vector(tapply(mass, blk, sum))
    cx <- as.vector(tapply(mass * nodes$x, blk, sum)) / M; cy <- as.vector(tapply(mass * nodes$y, blk, sum)) / M
    ok <- M > 0; M <- M[ok]; cx <- cx[ok]; cy <- cy[ok]
    xy <- sf::st_coordinates(sf::st_transform(sf::st_geometry(newdata), attr(nodes, "crs")))[, 1:2, drop = FALSE]
    E_mig <- numeric(nrow(xy))
    for (ch in split(seq_len(nrow(xy)), ceiling(seq_len(nrow(xy)) / 200))) {
      D <- sqrt(outer(xy[ch, 1], cx, "-")^2 + outer(xy[ch, 2], cy, "-")^2)
      E_mig[ch] <- as.vector((wend(D / mr) / (2 / 3 * mr)) %*% M)
    }
    E <- (1 - wm) * E_loc + wm * E_mig
  } else E <- E_loc
  X <- odisp_prediction_matrix(S$formula, sf::st_drop_geometry(newdata),
                                S$prevalence_spec, colnames(S$X))
  la <- log(pmax(E, 1e-300))
  gate_base <- as.vector(X %*% cf[p + seq_len(p)])
  pi_ <- if (identical(object$fld$occurrence_scale,"linear")) {
    stats::plogis(gate_base + h[["a_E"]]*E)
  } else ifelse(E > 0, stats::plogis(gate_base + h[["a_E"]] * la), 0)
  risk_exposure <- if (identical(object$fld$risk_scale,"linear")) E else la
  pz <- ifelse(E > 0, -expm1(-exp(as.vector(X %*% cf[seq_len(p)]) + h[["gamma"]] * risk_exposure)), 0)
  out <- data.frame(E = E, pi = pi_, p_zone = pz, prevalence = pi_ * pz, zero_exposure = E == 0)
  if (mig) { out$E_loc <- E_loc; out$E_mig <- E_mig; out$s_mig <- ifelse(E > 0, wm * E_mig / E, NA) }
  out
}
