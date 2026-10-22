#' Fit the dispersal model with a Whittle-Matern breeding field on the network
#'
#' Extends the model of \code{\link{odisp_fit}} (with the range fixed) by a
#' Gaussian Whittle-Matern field \eqn{u} on the river network in the breeding
#' intensity and an independent village effect \eqn{Z}:
#' \deqn{w(s) = \exp\{z(s)^\top\theta + u(s)\}, \quad
#' E(x) = \int_\Gamma K_\phi(\|x-s\|)\, w(s)\, ds, \quad
#' \lambda(x) = E(x) \exp\{d(x)^\top\beta + Z(x)\},}
#' with \eqn{Y \sim \mathrm{Bin}(n, 1 - e^{-\lambda})}. The field \eqn{u} has
#' \eqn{\alpha = 1} (exponential covariance along the network; see
#' \code{\link{wm_precision}}), standard deviation \code{sigma_u} and range
#' \code{range_u}; \eqn{Z \sim N(0, \sigma_Z^2)}.
#'
#' The latent vector \eqn{(\beta, \theta, u, Z)} is integrated out by a
#' Laplace approximation, with \eqn{\beta} and \eqn{\theta} given diffuse
#' Gaussian priors (standard deviation \code{prior_sd}) centred on the
#' estimates of the model without \eqn{u} and \eqn{Z}. The mode is found by
#' a truncated Newton method: the exact Hessian, which includes the curvature
#' of \eqn{\log E} in \eqn{(\theta, u)}, is inverted by conjugate gradients
#' preconditioned with the Gauss-Newton (Fisher) Hessian. That Hessian is the
#' prior precision plus a term of rank equal to the number of villages, so it
#' is inverted by the Woodbury identity with one sparse Cholesky factor of
#' the precision of \eqn{u}, computed once for each range. The
#' log-determinant in the approximation uses the Gauss-Newton Hessian.
#' \code{sigma_u} and \code{sigma_z} maximise the approximate marginal
#' likelihood for each value of \code{range_u}, so the range can be profiled
#' by giving several values.
#'
#' @inheritParams odisp_design
#' @param graph An \code{odisp_graph} object built from the same \code{nodes}.
#' @param phi Dispersal range (km), fixed.
#' @param range_u Range(s) of \eqn{u} (km) to fit; one fit per value.
#' @param u Logical; \code{FALSE} fits the model with the village effect only.
#' @param prior_sd Prior standard deviation of \eqn{\beta} and \eqn{\theta}.
#' @param start Optional starting values \code{c(log_sigma_u, log_sigma_z)}.
#' @param trace Logical; print progress.
#' @return An object of class \code{odisp_ufit} (or a list of them, one per
#'   value of \code{range_u}, with the profile in attribute \code{profile}).
#' @export
odisp_fit_u <- function(formula, breeding, data, nodes, graph, pairs = NULL, phi = 40,
                        range_u = 100, u = TRUE, trials = "n", aridity = "AI",
                        prior_sd = 10, start = NULL, trace = TRUE) {
  S <- odisp_u_setup(formula, breeding, data, nodes, graph, pairs, phi, trials, aridity, prior_sd)
  if (!u) return(odisp_u_optimise(S, NULL, start, trace))
  if (is.null(start)) {
    if (trace) cat("village effect only (to initialise sigma_z)\n")
    f_z <- odisp_u_optimise(S, NULL, NULL, trace)
    start <- c(log(0.5), log(f_z$hyper[["sigma_z"]]))
  }
  fits <- lapply(range_u, function(r) {
    if (trace) cat("range_u =", r, "km\n")
    f <- odisp_u_optimise(S, r, start, trace)
    start <<- f$log_hyper
    f
  })
  if (length(fits) == 1) return(fits[[1]])
  names(fits) <- range_u
  attr(fits, "profile") <- data.frame(range_u = range_u,
                                      logLik = sapply(fits, `[[`, "logLik"),
                                      sigma_u = sapply(fits, function(f) f$hyper["sigma_u"]),
                                      sigma_z = sapply(fits, function(f) f$hyper["sigma_z"]))
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
  list(y = d$y[has], n = d$n[has], X = d$D[has, , drop = FALSE], K = K, B = B,
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
  list(eta = log(E) + as.vector(S$X %*% x$beta) + x$z, V = Matrix::Diagonal(x = 1 / E) %*% KW, E = E)
}

odisp_u_prior <- function(x, S, R, sig) {
  q <- (sum((x$beta - S$m_beta)^2) + sum((x$theta - S$m_theta)^2)) / S$s2 + sum(x$z^2) / sig[2]^2
  if (!is.null(R)) q <- q + as.numeric(Matrix::crossprod(x$u, R %*% x$u)) / sig[1]^2
  q
}

odisp_u_vec <- function(x) c(x$beta, x$theta, x$u, x$z)
odisp_u_unvec <- function(v, x) {
  p <- length(x$beta); q <- length(x$theta); nu <- length(x$u)
  list(beta = v[seq_len(p)], theta = v[p + seq_len(q)],
       u = if (nu) v[p + q + seq_len(nu)] else NULL, z = v[p + q + nu + seq_along(x$z)])
}

# Operators at the current x, on flat latent vectors (beta, theta, u, z):
#   J v      Jacobian of eta;  J'w its transpose;
#   Pinv     inverse prior precision;
#   HGNinv   inverse of the Gauss-Newton (Fisher) Hessian P + J'DJ, by Woodbury;
#   H        exact negative Hessian of the log posterior, P + J' Dobs J - C,
#            where C is the curvature of log E in (theta, u).
odisp_u_ops <- function(x, ev, lk, S, Rch, R, sig) {
  n <- length(S$y); p <- ncol(S$X); q <- ncol(S$B); nu <- length(x$u)
  ib <- seq_len(p); it <- p + seq_len(q); iu <- p + q + seq_len(nu); iz <- p + q + nu + seq_len(n)
  Jt <- as.matrix(ev$V %*% S$B)
  M <- S$s2 * (tcrossprod(S$X) + tcrossprod(Jt)) + diag(sig[2]^2, n)
  Ju <- NULL
  if (nu) {
    Ju <- ev$V %*% S$A
    LW <- Matrix::solve(Rch, Matrix::solve(Rch, Matrix::t(Ju), system = "P"), system = "L")
    M <- M + sig[1]^2 * as.matrix(Matrix::crossprod(LW))
  }
  sD <- sqrt(lk$D)
  cB <- chol(diag(n) + outer(sD, sD) * M)
  lam <- exp(ev$eta)
  Dobs <- ifelse(S$y > 0, lam^2 * S$y * exp(-lam) / (-expm1(-lam))^2, 0) - lk$g
  cc <- as.vector(Matrix::crossprod(ev$V, lk$g))
  Jv <- function(v) as.vector(S$X %*% v[ib] + Jt %*% v[it]) + v[iz] +
    if (nu) as.vector(Ju %*% v[iu]) else 0
  Jtw <- function(w) c(as.vector(crossprod(S$X, w)), as.vector(crossprod(Jt, w)),
                       if (nu) as.vector(Matrix::crossprod(Ju, w)), w)
  Pinv <- function(r) c(S$s2 * r[ib], S$s2 * r[it],
                        if (nu) sig[1]^2 * as.vector(Matrix::solve(Rch, r[iu])), sig[2]^2 * r[iz])
  Pv <- function(v) c(v[ib] / S$s2, v[it] / S$s2,
                      if (nu) as.vector(R %*% v[iu]) / sig[1]^2, v[iz] / sig[2]^2)
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
  list(HGNinv = HGNinv, Hv = Hv, M = M, cB = cB, sD = sD, Jt = Jt)
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

# Laplace approximation for given range (prior precision R of u with unit sd
# and its Cholesky factor Rch) and standard deviations sig = c(sigma_u,
# sigma_z). The mode is found by truncated Newton with the exact Hessian; the
# log-determinant uses the Gauss-Newton (Fisher) Hessian, which is positive
# definite.
odisp_u_laplace <- function(sig, S, R, Rch, x, maxit = 50, tol = 1e-7) {
  f_of <- function(x) {
    ev <- odisp_u_eta(x, S); lk <- odisp_u_lik(ev$eta, S)
    v <- lk$ll - 0.5 * odisp_u_prior(x, S, R, sig)
    list(f = if (is.nan(v)) -Inf else v, ev = ev, lk = lk)
  }
  fx <- f_of(x)
  if (!is.finite(fx$f)) stop("non-finite objective at the starting point")
  cg <- 0
  for (it in seq_len(maxit)) {
    ops <- odisp_u_ops(x, fx$ev, fx$lk, S, Rch, R, sig)
    gl <- list(beta = as.vector(crossprod(S$X, fx$lk$g)) - (x$beta - S$m_beta) / S$s2,
               theta = as.vector(crossprod(ops$Jt, fx$lk$g)) - (x$theta - S$m_theta) / S$s2,
               u = if (!is.null(R)) as.vector(Matrix::crossprod(fx$ev$V %*% S$A, fx$lk$g) - R %*% x$u / sig[1]^2),
               z = fx$lk$g - x$z / sig[2]^2)
    gv <- odisp_u_vec(gl)
    pc <- odisp_u_pcg(gv, ops); cg <- cg + pc$iterations
    slope <- sum(gv * pc$d)
    if (!is.finite(slope) || slope <= 0) { pc$d <- ops$HGNinv(gv); slope <- sum(gv * pc$d) }
    if (slope < tol) break                                # Newton decrement
    xv <- odisp_u_vec(x); step <- 1
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
  ops <- odisp_u_ops(x, fx$ev, fx$lk, S, Rch, R, sig)
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
    ev <- eigen(H, symmetric = TRUE)
    if (all(ev$values < 0)) step <- -solve(H, g) else step <- g / max(abs(diag(H)), 1)
    if (sqrt(sum(step^2)) > radius) step <- step * radius / sqrt(sum(step^2))
    gain_pred <- sum(g * step) + 0.5 * sum(step * (H %*% step))
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

odisp_u_optimise <- function(S, range_u, start, trace) {
  x0 <- list(beta = S$m_beta, theta = S$m_theta,
             u = if (!is.null(range_u)) rep(0, S$graph$nv) else NULL, z = rep(0, length(S$y)))
  if (!is.null(range_u)) {
    R <- wm_precision(S$graph, range_u, 1)
    Rch <- Matrix::Cholesky(R, LDL = FALSE, perm = TRUE)
  } else { R <- NULL; Rch <- NULL }
  state <- new.env(); state$x <- x0; state$cache <- list()
  # approximate log-likelihood at log standard deviations lh; results are
  # cached, so points revisited by the optimiser are not recomputed
  eval_lh <- function(lh) {
    key <- paste(sprintf("%.6f", lh), collapse = ",")
    if (!is.null(state$cache[[key]])) return(state$cache[[key]])
    sig <- if (is.null(R)) c(0, exp(lh)) else exp(lh)
    res <- odisp_u_laplace(sig, S, R, Rch, state$x)
    state$x <- res$x; state$cache[[key]] <- res
    if (trace) cat(sprintf("  sigma = %s  logLik = %.3f  (%d iterations)\n",
                           paste(sprintf("%.3f", sig[sig > 0]), collapse = ", "), res$logLik, res$iterations))
    res
  }
  obj <- function(lh) -eval_lh(lh)$logLik
  # The approximate likelihood is computed to about 1e-6, so the optimisers
  # stop when the log-likelihood changes by less than about 1e-3: greater
  # precision in the standard deviations has no practical meaning.
  if (is.null(R)) {
    o <- stats::optimize(obj, c(log(0.01), log(5)), tol = 1e-3)
    lh <- o$minimum
  } else {
    if (is.null(start)) start <- c(log(0.5), log(0.5))
    lh <- odisp_newton_max(function(lh) eval_lh(lh)$logLik, start)
  }
  sig <- if (is.null(R)) c(0, exp(lh)) else exp(lh)
  res <- eval_lh(lh)
  if (trace) cat(sprintf("  converged: sigma = %s, logLik = %.3f, %d evaluations\n",
                         paste(sprintf("%.3f", sig[sig > 0]), collapse = ", "), res$logLik,
                         length(state$cache)))
  p <- length(S$m_beta); q <- length(S$m_theta)
  est <- c(res$x$beta, res$x$theta)
  names(est) <- c(paste0("beta.", colnames(S$X)), paste0("theta.", S$nodes_names))
  dimnames(res$Vfix) <- list(names(est), names(est))
  structure(list(coefficients = est, vcov = res$Vfix,
                 hyper = c(sigma_u = sig[1], sigma_z = sig[2], range_u = if (is.null(range_u)) NA else range_u),
                 log_hyper = lh, logLik = res$logLik, u = res$x$u, z = res$x$z, eta = res$eta, E = res$E,
                 fitted = 1 - exp(-exp(res$eta)), setup = S, converged = res$converged),
            class = "odisp_ufit")
}

#' @export
print.odisp_ufit <- function(x, ...) {
  cat("River-network dispersal model with a Whittle-Matern breeding field\n")
  cat("Dispersal range fixed at", x$setup$phi, "km\n")
  h <- x$hyper
  if (is.na(h["range_u"])) cat("No field u; village effect sd", round(h["sigma_z"], 3), "\n")
  else cat(sprintf("Field u: range %g km, sd %.3f; village effect sd %.3f\n", h["range_u"], h["sigma_u"], h["sigma_z"]))
  se <- sqrt(diag(x$vcov))
  print(round(cbind(estimate = x$coefficients, se = se), 3))
  cat("Laplace log-likelihood:", round(x$logLik, 2), "\n")
  invisible(x)
}

#' Predict from a model with a breeding field on the network
#'
#' Plug-in prediction at the posterior mode of \eqn{(\beta, \theta, u)},
#' with the village effect set to zero.
#'
#' @param object An \code{odisp_ufit} object.
#' @param newdata An \code{sf} object of prediction locations with the
#'   covariates of the prevalence formula.
#' @param nodes,graph The network used in the fit.
#' @param pairs Optional \code{odisp_pairs} for \code{newdata}.
#' @param ... Not used.
#' @return A data frame with exposure \code{E}, force of infection
#'   \code{lambda}, \code{prevalence} and \code{zero_exposure}; attribute
#'   \code{u_nodes} gives the mode of \eqn{u} at every node.
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
  lam <- E * exp(as.vector(X %*% object$coefficients[seq_len(length(S$m_beta))]))
  out <- data.frame(E = E, lambda = lam, prevalence = 1 - exp(-lam), zero_exposure = E == 0)
  if (!is.null(object$u)) attr(out, "u_nodes") <- as.vector(graph$A %*% object$u)
  out
}
