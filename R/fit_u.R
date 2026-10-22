#' Fit the dispersal model with a spatial random effect in the exposure or in the intensity
#'
#' Extends the model of \code{\link{odisp_fit}} (with the range fixed) by a
#' Gaussian random effect, either a Whittle-Matern field \eqn{u} on the river
#' network in the breeding intensity (the exposure) or a Gaussian process
#' \eqn{S} on the plane in the intensity of transmission, optionally with an
#' independent village effect \eqn{Z}:
#' \deqn{w(s) = \exp\{z(s)^\top\theta + u(s)\}, \quad
#' E(x) = \int_\Gamma K_\phi(\|x-s\|)\, w(s)\, ds, \quad
#' \lambda(x) = E(x) \exp\{d(x)^\top\beta + S(x) + Z(x)\},}
#' with \eqn{Y \sim \mathrm{Bin}(n, 1 - e^{-\lambda})}. The field \eqn{u} has
#' \eqn{\alpha = 1} (exponential covariance along the network; see
#' \code{\link{wm_precision}}); \eqn{S} has exponential covariance
#' \eqn{\sigma_S^2 \exp(-2d/\rho_S)} in straight-line distance, so that both
#' ranges are the distance at which the correlation is \eqn{e^{-2}}; and
#' \eqn{Z \sim N(0, \sigma_Z^2)}.
#'
#' The latent vector, made of the regression coefficients, the field and the
#' village effects, is integrated out by a Laplace approximation, with the
#' regression coefficients given diffuse Gaussian priors (standard deviation
#' \code{prior_sd}) centred on the estimates of the model without random
#' effects. The mode is found by a truncated Newton method: the exact Hessian,
#' which includes the curvature of \eqn{\log E} in \eqn{(\theta, u)}, is
#' inverted by conjugate gradients preconditioned with the Gauss-Newton
#' (Fisher) Hessian. That Hessian is the prior precision plus a term of rank
#' equal to the number of villages, so it is inverted by the Woodbury
#' identity, using a sparse Cholesky factor of the precision of \eqn{u} or a
#' dense Cholesky factor of the correlation matrix of \eqn{S}, computed once
#' for each range. The log-determinant in the approximation uses the
#' Gauss-Newton Hessian. The standard deviations maximise the approximate
#' marginal likelihood for each value of the range, so the range can be
#' profiled by giving several values.
#'
#' @inheritParams odisp_design
#' @param graph An \code{odisp_graph} object built from the same \code{nodes}.
#' @param phi Dispersal range (km), fixed.
#' @param field The spatial random effect: \code{"u"} (on the river network,
#'   in the exposure), \code{"S"} (on the plane, in the intensity of
#'   transmission), \code{"uS"} (both) or \code{"none"}.
#' @param range Range(s) of the field (km) to fit; one fit per value. For
#'   \code{field = "uS"}, a vector \code{c(range_u, range_S)} or a two-column
#'   matrix with one pair of ranges per row.
#' @param z Logical; include the independent village effect \eqn{Z}.
#' @param prior_sd Prior standard deviation of \eqn{\beta} and \eqn{\theta}.
#' @param start Optional starting values of the log standard deviations of
#'   \eqn{u}, \eqn{S} and \eqn{Z}, in that order, for those in the model.
#' @param trace Logical; print progress.
#' @param range_u,u Deprecated: \code{range_u} is \code{range} with
#'   \code{field = "u"}, and \code{u = FALSE} is \code{field = "none"}.
#' @return An object of class \code{odisp_ufit} (or a list of them, one per
#'   value of \code{range}, with the profile in attribute \code{profile}).
#' @export
odisp_fit_u <- function(formula, breeding, data, nodes, graph, pairs = NULL, phi = 40,
                        field = c("u", "S", "uS", "none"), range = 100, z = TRUE, trials = "n",
                        aridity = "AI", prior_sd = 10, start = NULL, trace = TRUE,
                        range_u = NULL, u = TRUE) {
  field <- match.arg(field)
  if (!is.null(range_u)) range <- range_u
  if (!u) field <- "none"
  if (field == "none" && !z) stop("the model needs a field or the village effect")
  S <- odisp_u_setup(formula, breeding, data, nodes, graph, pairs, phi, trials, aridity, prior_sd)
  if (field == "none") return(odisp_u_optimise(S, "none", NULL, TRUE, start, trace))
  nf <- if (field == "uS") 2 else 1
  if (is.null(start)) {
    sf <- rep(log(if (nf == 2) 1.5 else 2), nf)
    if (z) {
      if (trace) cat("village effect only (to initialise sigma_z)\n")
      f_z <- odisp_u_optimise(S, "none", NULL, TRUE, NULL, trace)
      # the field typically takes over half of the village effect
      start <- c(sf, log(f_z$hyper[["sigma_z"]] / 2))
    } else start <- sf
  }
  ranges <- if (field == "uS") matrix(range, ncol = 2, dimnames = list(NULL, c("u", "S"))) else matrix(range, ncol = 1)
  fits <- lapply(seq_len(nrow(ranges)), function(k) {
    r <- ranges[k, ]
    if (trace) cat(sprintf("field %s, range = %s km%s\n", field, paste(r, collapse = " (u), "),
                           if (z) " (with village effect)" else ""))
    if (trace && field == "uS") cat("  [ranges: u along the network, S straight-line]\n")
    f <- odisp_u_optimise(S, field, r, z, start, trace)
    start <<- f$log_hyper
    f
  })
  if (length(fits) == 1) return(fits[[1]])
  names(fits) <- apply(ranges, 1, paste, collapse = "_")
  prof <- data.frame(ranges, logLik = sapply(fits, `[[`, "logLik"))
  names(prof)[seq_len(ncol(ranges))] <- if (field == "uS") c("range_u", "range_s") else "range"
  if (field != "uS") prof$sigma_field <- sapply(fits, function(f) f$hyper[["sigma_field"]])
  if (field %in% c("u", "uS")) prof$sigma_u <- sapply(fits, function(f) f$hyper[["sigma_u"]])
  if (field %in% c("S", "uS")) prof$sigma_s <- sapply(fits, function(f) f$hyper[["sigma_s"]])
  prof$sigma_z <- sapply(fits, function(f) f$hyper[["sigma_z"]])
  rownames(prof) <- NULL
  attr(fits, "profile") <- prof
  fits
}

odisp_u_setup <- function(formula, breeding, data, nodes, graph, pairs, phi, trials, aridity, prior_sd) {
  if (!inherits(graph, "odisp_graph")) stop("'graph' must come from river_graph()")
  if (nrow(graph$A) != nrow(nodes)) stop("'graph' was built from different nodes")
  if (is.null(pairs)) pairs <- river_pairs(data, nodes, phi)
  # starting point and prior centre: the model without u and Z
  f0 <- suppressWarnings(odisp_fit(formula, breeding, data = data, nodes = nodes, pairs = pairs,
                                   max_range = pairs$max_range, trials = trials, aridity = aridity,
                                   fixed = list(gamma = 0, eta0 = log(phi))))
  d <- f0$design
  keep <- d$r < phi
  i <- d$i[keep]; jj <- d$j[keep]; r <- d$r[keep]
  used_all <- sort(unique(pairs$j))              # d$j indexes these
  jnode <- used_all[jj]
  used <- sort(unique(jnode))
  t <- r / phi
  kx <- nodes$len[jnode] * (1 - t)^4 * (1 + 4 * t) / (2 / 3 * phi)
  has <- sort(unique(i))
  if (any(d$y[-has] > 0)) stop("villages with positives are beyond the reach of every river")
  K <- Matrix::sparseMatrix(i = match(i, has), j = match(jnode, used), x = kx,
                            dims = c(length(has), length(used)))
  B <- d$Z[match(used, used_all), , drop = FALSE]
  cf <- coef(f0)
  p <- ncol(d$D); q <- ncol(d$Z)
  xy <- sf::st_coordinates(sf::st_transform(sf::st_geometry(data), attr(nodes, "crs")))[has, 1:2, drop = FALSE]
  list(y = d$y[has], n = d$n[has], X = d$D[has, , drop = FALSE], K = K, B = B, xy = xy, lai = d$logAI[has],
       A = graph$A[used, , drop = FALSE], graph = graph, used = used, has = has, m_all = d$m,
       m_beta = cf[seq_len(p)], m_theta = cf[p + seq_len(q)], s2 = prior_sd^2,
       lchoose = sum(lchoose(d$n[has], d$y[has])), fit0 = f0, phi = phi,
       formula = formula, breeding = breeding, nodes_names = colnames(d$Z))
}

# log-likelihood, score and Fisher information with respect to eta = log(lambda)
odisp_u_lik <- function(eta, S) {
  if (any(!is.finite(eta)) || max(eta) > 30) return(list(ll = -Inf, g = NA, D = NA))
  lam <- exp(eta)
  logp <- log1mexp(lam)
  t1 <- ifelse(S$y > 0, S$y * logp, 0)
  t2 <- ifelse(S$n > S$y, (S$n - S$y) * lam, 0)
  ll <- S$lchoose + sum(t1 - t2)
  if (is.na(ll)) ll <- -Inf
  g <- lam * (ifelse(S$y > 0, S$y / exp(logp), 0) - S$n)
  D <- S$n * lam^2 / expm1(lam)
  D[!is.finite(D)] <- 0
  list(ll = ll, g = g, D = D)
}

odisp_u_eta <- function(x, S) {
  w <- exp(as.vector(S$B %*% x$theta) + if (!is.null(x$u)) as.vector(S$A %*% x$u) else 0)
  KW <- S$K %*% Matrix::Diagonal(x = w)
  E <- Matrix::rowSums(KW)
  eta <- log(E) + as.vector(S$X %*% x$beta)
  if (!is.null(x$s)) eta <- eta + x$s
  if (!is.null(x$z)) eta <- eta + x$z
  list(eta = eta, V = Matrix::Diagonal(x = 1 / E) %*% KW, E = E)
}

# Random-effect structure for given standard deviations: `fld` is a list
# with the field type, its unit-variance precision R (u, sparse) and Cholesky
# factor Rch, or the Cholesky factor LC of the correlation matrix (S), and
# the standard deviations sig_u, sig_s (fields) and sig_z (village effect).
odisp_u_prior <- function(x, S, fld) {
  q <- (sum((x$beta - S$m_beta)^2) + sum((x$theta - S$m_theta)^2)) / S$s2
  if (!is.null(x$z)) q <- q + sum(x$z^2) / fld$sig_z^2
  if (!is.null(x$u)) q <- q + as.numeric(Matrix::crossprod(x$u, fld$R %*% x$u)) / fld$sig_u^2
  if (!is.null(x$s)) q <- q + sum(backsolve(fld$LC, x$s, transpose = TRUE)^2) / fld$sig_s^2
  q
}

odisp_u_vec <- function(x) c(x$beta, x$theta, x$u, x$s, x$z)
odisp_u_unvec <- function(v, x) {
  p <- length(x$beta); q <- length(x$theta); nu <- length(x$u); ns <- length(x$s); nz <- length(x$z)
  list(beta = v[seq_len(p)], theta = v[p + seq_len(q)],
       u = if (nu) v[p + q + seq_len(nu)] else NULL,
       s = if (ns) v[p + q + nu + seq_len(ns)] else NULL,
       z = if (nz) v[p + q + nu + ns + seq_len(nz)] else NULL)
}

# Operators at the current x, on flat latent vectors (beta, theta, u, s, z):
#   J v      Jacobian of eta;  J'w its transpose;
#   Pinv     inverse prior precision;  Pv prior precision times v;
#   HGNinv   inverse of the Gauss-Newton (Fisher) Hessian P + J'DJ, by Woodbury;
#   H        exact negative Hessian of the log posterior, P + J' Dobs J - C,
#            where C is the curvature of log E in (theta, u).
odisp_u_ops <- function(x, ev, lk, S, fld) {
  n <- length(S$y); p <- ncol(S$X); q <- ncol(S$B)
  nu <- length(x$u); ns <- length(x$s); nz <- length(x$z)
  ib <- seq_len(p); it <- p + seq_len(q); iu <- p + q + seq_len(nu)
  is <- p + q + nu + seq_len(ns); iz <- p + q + nu + ns + seq_len(nz)
  Jt <- as.matrix(ev$V %*% S$B)
  M <- S$s2 * (tcrossprod(S$X) + tcrossprod(Jt))
  if (nz) M <- M + diag(fld$sig_z^2, n)
  if (ns) M <- M + fld$sig_s^2 * crossprod(fld$LC)
  Ju <- NULL
  if (nu) {
    Ju <- ev$V %*% S$A
    LW <- Matrix::solve(fld$Rch, Matrix::solve(fld$Rch, Matrix::t(Ju), system = "P"), system = "L")
    M <- M + fld$sig_u^2 * as.matrix(Matrix::crossprod(LW))
  }
  sD <- sqrt(lk$D)
  cB <- chol(diag(n) + outer(sD, sD) * M)
  lam <- exp(ev$eta)
  Dobs <- ifelse(S$y > 0, lam^2 * S$y * exp(-lam) / (-expm1(-lam))^2, 0) - lk$g
  cc <- as.vector(Matrix::crossprod(ev$V, lk$g))
  Jv <- function(v) {
    out <- as.vector(S$X %*% v[ib] + Jt %*% v[it])
    if (nu) out <- out + as.vector(Ju %*% v[iu])
    if (ns) out <- out + v[is]
    if (nz) out <- out + v[iz]
    out
  }
  Jtw <- function(w) c(as.vector(crossprod(S$X, w)), as.vector(crossprod(Jt, w)),
                       if (nu) as.vector(Matrix::crossprod(Ju, w)), if (ns) w, if (nz) w)
  Pinv <- function(r) c(S$s2 * r[ib], S$s2 * r[it],
                        if (nu) fld$sig_u^2 * as.vector(Matrix::solve(fld$Rch, r[iu])),
                        if (ns) fld$sig_s^2 * as.vector(crossprod(fld$LC, fld$LC %*% r[is])),
                        if (nz) fld$sig_z^2 * r[iz])
  Pv <- function(v) c(v[ib] / S$s2, v[it] / S$s2,
                      if (nu) as.vector(fld$R %*% v[iu]) / fld$sig_u^2,
                      if (ns) backsolve(fld$LC, backsolve(fld$LC, v[is], transpose = TRUE)) / fld$sig_s^2,
                      if (nz) v[iz] / fld$sig_z^2)
  HGNinv <- function(r) {
    s0 <- Pinv(r)
    w <- sD * backsolve(cB, backsolve(cB, sD * Jv(s0), transpose = TRUE))
    s0 - Pinv(Jtw(w))
  }
  Cv <- function(v) {
    nuv <- as.vector(S$B %*% v[it]) + if (nu) as.vector(S$A %*% v[iu]) else 0
    cx <- cc * nuv - as.vector(Matrix::crossprod(ev$V, lk$g * as.vector(ev$V %*% nuv)))
    out <- numeric(length(v))
    out[it] <- as.vector(crossprod(S$B, cx))
    if (nu) out[iu] <- as.vector(Matrix::crossprod(S$A, cx))
    out
  }
  Hv <- function(v) Pv(v) + Jtw(Dobs * Jv(v)) - Cv(v)
  list(HGNinv = HGNinv, Hv = Hv, Jtw = Jtw, Pv = Pv, M = M, cB = cB, sD = sD, Jt = Jt)
}

# Truncated Newton direction: preconditioned conjugate gradients on H d = g,
# stopped at negative curvature (Steihaug); the first iterate is the
# Gauss-Newton step.
odisp_u_pcg <- function(g, ops, maxit = 50, rtol = 1e-4) {
  d <- numeric(length(g)); r <- g; z <- ops$HGNinv(r); pdir <- z
  rz <- sum(r * z); r0 <- sqrt(rz)
  for (k in seq_len(maxit)) {
    Hp <- ops$Hv(pdir); curv <- sum(pdir * Hp)
    if (!is.finite(curv) || curv <= 0) { if (k == 1) d <- z; break }
    al <- rz / curv
    d <- d + al * pdir; r <- r - al * Hp
    z <- ops$HGNinv(r); rz_new <- sum(r * z)
    if (sqrt(abs(rz_new)) < rtol * r0) break
    pdir <- z + (rz_new / rz) * pdir; rz <- rz_new
  }
  list(d = d, iterations = k)
}

# Laplace approximation for the random-effect structure `fld` (see
# odisp_u_prior). The mode is found by truncated Newton with the exact
# Hessian; the log-determinant uses the Gauss-Newton (Fisher) Hessian, which
# is positive definite.
odisp_u_laplace <- function(fld, S, x, maxit = 50, tol = 1e-7) {
  m <- odisp_u_vec(list(beta = S$m_beta, theta = S$m_theta, u = 0 * x$u, s = 0 * x$s, z = 0 * x$z))
  f_of <- function(x) {
    ev <- odisp_u_eta(x, S); lk <- odisp_u_lik(ev$eta, S)
    v <- lk$ll - 0.5 * odisp_u_prior(x, S, fld)
    list(f = if (is.nan(v)) -Inf else v, ev = ev, lk = lk)
  }
  fx <- f_of(x)
  if (!is.finite(fx$f)) stop("non-finite objective at the starting point")
  cg <- 0
  for (it in seq_len(maxit)) {
    ops <- odisp_u_ops(x, fx$ev, fx$lk, S, fld)
    xv <- odisp_u_vec(x)
    gv <- ops$Jtw(fx$lk$g) - ops$Pv(xv - m)                    # gradient of the log posterior
    pc <- odisp_u_pcg(gv, ops); cg <- cg + pc$iterations
    slope <- sum(gv * pc$d)
    if (!is.finite(slope) || slope <= 0) { pc$d <- ops$HGNinv(gv); slope <- sum(gv * pc$d) }
    if (slope < tol) break                                     # Newton decrement
    step <- 1
    repeat {
      xs <- odisp_u_unvec(xv + step * pc$d, x); fs <- f_of(xs)
      if (fs$f >= fx$f + 1e-4 * step * slope || step < 1e-8) break
      step <- step / 2
    }
    if (fs$f < fx$f) break
    gain <- fs$f - fx$f
    x <- xs; fx <- fs
    if (gain < tol) break
  }
  ops <- odisp_u_ops(x, fx$ev, fx$lk, S, fld)
  logdet <- 2 * sum(log(diag(ops$cB)))
  # covariance of (beta, theta) given the hyperparameters (Gauss-Newton)
  Jf <- cbind(S$X, ops$Jt)
  G <- backsolve(ops$cB, ops$sD * Jf, transpose = TRUE)
  Vf <- S$s2 * diag(ncol(Jf)) - S$s2^2 * crossprod(G)
  list(logLik = fx$f - 0.5 * logdet, x = x, eta = fx$ev$eta, E = fx$ev$E, iterations = it,
       cg_iterations = cg, Vfix = Vf, converged = it < maxit)
}

# Maximise a smooth function of a few parameters by Newton's method with
# central finite-difference derivatives (step h) and a trust region (steps of
# at most `radius`). Stops when the predicted gain is below `tol` or after
# `maxit` iterations. Suited to the approximate marginal likelihood, which is
# smooth and close to quadratic in the log standard deviations but computed
# with a small numerical error, and whose evaluations are cheap near the
# current point (warm start) and expensive far from it.
odisp_newton_max <- function(f, x, h = 0.05, radius = 0.5, tol = 1e-3, maxit = 20) {
  d <- length(x); fx <- f(x)
  for (it in seq_len(maxit)) {
    E <- diag(h, d)
    fp <- sapply(seq_len(d), function(k) f(x + E[, k]))
    fm <- sapply(seq_len(d), function(k) f(x - E[, k]))
    g <- (fp - fm) / (2 * h)
    H <- diag((fp - 2 * fx + fm) / h^2, d)
    if (d > 1) for (k in 1:(d - 1)) for (l in (k + 1):d) {
      H[k, l] <- H[l, k] <- (f(x + E[, k] + E[, l]) - fp[k] - fp[l] + fx) / h^2
    }
    # Newton step with the Hessian made negative definite (eigenvalues
    # replaced by -max(|lambda|, 1)), so that directions in which the function
    # is flat or convex, such as a standard deviation close to zero, give a
    # step along the gradient, limited by the trust region
    ev <- eigen(H, symmetric = TRUE)
    Hm <- ev$vectors %*% diag(-pmax(abs(ev$values), 1), d) %*% t(ev$vectors)
    step <- -solve(Hm, g)
    if (sqrt(sum(step^2)) > radius) step <- step * radius / sqrt(sum(step^2))
    gain_pred <- sum(g * step) + 0.5 * sum(step * (Hm %*% step))
    if (gain_pred < tol) break
    repeat {
      fn <- f(x + step)
      if (fn > fx || sqrt(sum(step^2)) < 1e-3) break
      step <- step / 2
    }
    if (fn <= fx) break
    x <- x + step; done <- fn - fx < tol; fx <- fn
    if (done) break
  }
  x
}

# exponential correlation matrix for S, range rho (correlation exp(-2) at rho)
odisp_s_corr <- function(xy, rho) exp(-2 * as.matrix(stats::dist(xy)) / rho)

odisp_u_optimise <- function(S, field, range, z, start, trace) {
  n <- length(S$y)
  has_u <- field %in% c("u", "uS"); has_s <- field %in% c("S", "uS")
  r_u <- if (has_u) unname(range[1]) else NA
  r_s <- if (has_s) unname(range[length(range)]) else NA
  x0 <- list(beta = S$m_beta, theta = S$m_theta,
             u = if (has_u) rep(0, S$graph$nv) else NULL,
             s = if (has_s) rep(0, n) else NULL,
             z = if (z) rep(0, n) else NULL)
  base <- list(field = field)
  if (has_u) {
    base$R <- wm_precision(S$graph, r_u, 1)
    base$Rch <- Matrix::Cholesky(base$R, LDL = FALSE, perm = TRUE)
  }
  if (has_s) base$LC <- chol(odisp_s_corr(S$xy, r_s) + diag(1e-8, n))
  make_fld <- function(lh) {
    fld <- base; k <- 0
    if (has_u) { k <- k + 1; fld$sig_u <- exp(lh[k]) }
    if (has_s) { k <- k + 1; fld$sig_s <- exp(lh[k]) }
    if (z) { k <- k + 1; fld$sig_z <- exp(lh[k]) }
    fld
  }
  state <- new.env(); state$x <- x0; state$cache <- list()
  # approximate log-likelihood at log standard deviations lh; results are
  # cached, so points revisited by the optimiser are not recomputed
  eval_lh <- function(lh) {
    key <- paste(sprintf("%.6f", lh), collapse = ",")
    if (!is.null(state$cache[[key]])) return(state$cache[[key]])
    res <- odisp_u_laplace(make_fld(lh), S, state$x)
    state$x <- res$x; state$cache[[key]] <- res
    if (trace) cat(sprintf("  sigma = %s  logLik = %.3f  (%d iterations)\n",
                           paste(sprintf("%.3f", exp(lh)), collapse = ", "), res$logLik, res$iterations))
    res
  }
  # The approximate likelihood is computed to about 1e-6, so the optimisers
  # stop when the log-likelihood changes by less than about 1e-3: greater
  # precision in the standard deviations has no practical meaning.
  npar <- has_u + has_s + z
  if (is.null(start)) start <- rep(log(0.5), npar)
  if (npar == 1) {
    lo <- start - 3; hi <- start + 3
    lh <- stats::optimize(function(l) -eval_lh(l)$logLik, c(lo, hi), tol = 1e-3)$minimum
  } else {
    lh <- odisp_newton_max(function(lh) eval_lh(lh)$logLik, start)
  }
  res <- eval_lh(lh)
  fld <- make_fld(lh)
  if (trace) cat(sprintf("  converged: sigma = %s, logLik = %.3f, %d evaluations\n",
                         paste(sprintf("%.3f", exp(lh)), collapse = ", "), res$logLik,
                         length(state$cache)))
  est <- c(res$x$beta, res$x$theta)
  names(est) <- c(paste0("beta.", colnames(S$X)), paste0("theta.", S$nodes_names))
  dimnames(res$Vfix) <- list(names(est), names(est))
  sig_field <- if (field == "u") fld$sig_u else if (field == "S") fld$sig_s else NA
  structure(list(coefficients = est, vcov = res$Vfix, field = field,
                 hyper = c(sigma_u = if (has_u) fld$sig_u else 0, range_u = r_u,
                           sigma_s = if (has_s) fld$sig_s else 0, range_s = r_s,
                           sigma_z = if (z) fld$sig_z else 0,
                           # single field: its sd and range
                           sigma_field = if (field %in% c("u", "S")) sig_field else NA,
                           range = if (field == "u") r_u else if (field == "S") r_s else NA),
                 log_hyper = lh, logLik = res$logLik, u = res$x$u, s = res$x$s,
                 z = if (z) res$x$z else rep(0, n), eta = res$eta, E = res$E,
                 fitted = 1 - exp(-exp(res$eta)), setup = S, converged = res$converged),
            class = "odisp_ufit")
}

#' @export
print.odisp_ufit <- function(x, ...) {
  cat("River-network dispersal model, range fixed at", x$setup$phi, "km\n")
  h <- x$hyper; fld <- if (is.null(x$field)) "u" else x$field
  if (fld %in% c("u", "uS"))
    cat(sprintf("Field u on the river network (in the exposure): range %g km, sd %.3f\n", h[["range_u"]], h[["sigma_u"]]))
  if (fld %in% c("S", "uS"))
    cat(sprintf("Gaussian process S on the plane (in the intensity): range %g km, sd %.3f\n", h[["range_s"]], h[["sigma_s"]]))
  if (fld == "none") cat("No spatial field\n")
  cat(if (h[["sigma_z"]] > 0) sprintf("Village effect sd %.3f\n", h[["sigma_z"]]) else "No village effect\n")
  se <- sqrt(diag(x$vcov))
  print(round(cbind(estimate = x$coefficients, se = se), 3))
  cat("Laplace log-likelihood:", round(x$logLik, 2), "\n")
  invisible(x)
}

#' Predict from a model with a spatial random effect
#'
#' Plug-in prediction at the posterior mode of the regression coefficients
#' and of the field, with the village effect set to zero. For the field
#' \eqn{S} on the plane, its value at the new locations is the conditional
#' mean given its mode at the villages (simple kriging).
#'
#' @param object An \code{odisp_ufit} object.
#' @param newdata An \code{sf} object of prediction locations with the
#'   covariates of the prevalence formula.
#' @param nodes,graph The network used in the fit.
#' @param pairs Optional \code{odisp_pairs} for \code{newdata}.
#' @param ... Not used.
#' @return A data frame with exposure \code{E}, force of infection
#'   \code{lambda}, \code{prevalence} and \code{zero_exposure}; attribute
#'   \code{u_nodes} gives the mode of \eqn{u} at every node, and column
#'   \code{S} the predicted \eqn{S} when that field is in the model.
#' @export
predict.odisp_ufit <- function(object, newdata, nodes, graph, pairs = NULL, ...) {
  S <- object$setup; phi <- S$phi
  if (is.null(pairs)) pairs <- river_pairs(newdata, nodes, phi)
  keep <- pairs$r < phi
  i <- pairs$i[keep]; j <- pairs$j[keep]; t <- pairs$r[keep] / phi
  Bn <- stats::model.matrix(S$breeding, as.data.frame(nodes)[j, , drop = FALSE])
  Bn <- Bn[, colnames(Bn) != "(Intercept)", drop = FALSE]
  uj <- if (!is.null(object$u)) as.vector(graph$A[j, , drop = FALSE] %*% object$u) else 0
  q <- length(S$m_theta)
  w <- exp(as.vector(Bn %*% object$coefficients[length(S$m_beta) + seq_len(q)]) + uj)
  E <- as.vector(tapply(c(nodes$len[j] * (1 - t)^4 * (1 + 4 * t) / (2 / 3 * phi) * w, rep(0, pairs$m)),
                        c(i, seq_len(pairs$m)), sum))
  X <- stats::model.matrix(stats::delete.response(stats::terms(S$formula)), sf::st_drop_geometry(newdata))
  lin <- as.vector(X %*% object$coefficients[seq_len(length(S$m_beta))])
  s_new <- NULL
  if (!is.null(object$s)) {
    xy_new <- sf::st_coordinates(sf::st_transform(sf::st_geometry(newdata), attr(nodes, "crs")))[, 1:2, drop = FALSE]
    rho <- if ("range_s" %in% names(object$hyper)) object$hyper[["range_s"]] else object$hyper[["range"]]
    d_nv <- sqrt(outer(xy_new[, 1], S$xy[, 1], "-")^2 + outer(xy_new[, 2], S$xy[, 2], "-")^2)
    C <- odisp_s_corr(S$xy, rho) + diag(1e-8, nrow(S$xy))
    s_new <- as.vector(exp(-2 * d_nv / rho) %*% solve(C, object$s))
    lin <- lin + s_new
  }
  lam <- E * exp(lin)
  out <- data.frame(E = E, lambda = lam, prevalence = 1 - exp(-lam), zero_exposure = E == 0)
  if (!is.null(s_new)) out$S <- s_new
  if (!is.null(object$u)) attr(out, "u_nodes") <- as.vector(graph$A %*% object$u)
  out
}
