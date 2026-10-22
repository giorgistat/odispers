#' Zero-inflated river-network model with a breeding field
#'
#' A village is either outside any zone of local transmission, with no
#' positives, or inside one, with binomial counts. Both the probability of
#' being inside a zone and the force of infection depend on the exposure from
#' the river network:
#' \deqn{\pi(x) = 1 - \exp\{-\alpha E(x) \mathrm{AI}(x)^\delta\}, \qquad
#' \lambda(x) = E(x)^\gamma \exp\{d(x)^\top\beta + Z(x)\},}
#' \deqn{P(Y = 0) = 1 - \pi + \pi (1 - p)^n, \quad
#' P(Y = y) = \pi\, \mathrm{Bin}(y; n, p) \ (y > 0), \quad p = 1 - e^{-\lambda},}
#' with \eqn{E(x) = \int_\Gamma K_\phi(\|x-s\|) \exp\{z(s)^\top\theta + u(s)\} ds}
#' as in \code{\link{odisp_fit_u}}, \eqn{u} a Whittle-Matern field on the river
#' network and \eqn{Z} an independent village effect. The aridity index
#' \eqn{\mathrm{AI}} in \eqn{\pi} lets climate determine, with the rivers,
#' whether a village lies in a zone of transmission (\eqn{\delta = 0} by
#' default, i.e. rivers only). With \eqn{\gamma = 0}
#' the rivers determine only whether a village is in a transmission zone and
#' the intensity of transmission alone determines prevalence within it; with
#' \eqn{\gamma = 1} exposure also scales the force of infection; and as
#' \eqn{\alpha \to \infty} the model reduces to that of
#' \code{\link{odisp_fit_u}}.
#'
#' The latent vector \eqn{(\beta, \theta, u, Z)} is integrated out by a
#' Laplace approximation with the same algorithm as \code{\link{odisp_fit_u}},
#' applied to the two predictors \eqn{\log E} and \eqn{d^\top\beta + Z} of
#' each village; the Hessian used for the preconditioner and the
#' log-determinant is the observed information of the two predictors, with
#' eigenvalues floored at a small positive value. \eqn{\sigma_u},
#' \eqn{\sigma_Z}, \eqn{\alpha}, \eqn{\gamma} and \eqn{\delta} maximise
#' the approximate marginal likelihood for the given range of \eqn{u}.
#'
#' @inheritParams odisp_fit_u
#' @param range_u Range of \eqn{u} (km).
#' @param start Optional starting values \code{c(log sigma_u, log sigma_z,
#'   log alpha, logit gamma, delta)}.
#' @param fixed Optional named list fixing some of \code{sigma_u},
#'   \code{sigma_z}, \code{alpha}, \code{gamma}, \code{delta} (on their
#'   natural scale). By default \code{delta = 0}: aridity is not in \eqn{\pi}.
#' @return An object of class \code{odisp_zifit}.
#' @export
odisp_fit_zi <- function(formula, breeding, data, nodes, graph, pairs = NULL, phi = 40,
                         range_u = 800, trials = "n", aridity = "AI", prior_sd = 10,
                         start = NULL, fixed = list(delta = 0), trace = TRUE) {
  S <- odisp_u_setup(formula, breeding, data, nodes, graph, pairs, phi, trials, aridity, prior_sd)
  n <- length(S$y)
  R <- wm_precision(graph, range_u, 1)
  base <- list(field = "u", R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE))
  if (is.null(start)) start <- c(log(2), log(1), log(1), 0, 0)
  if (length(start) == 4) start <- c(start, 0)
  names(start) <- c("sigma_u", "sigma_z", "alpha", "gamma", "delta")
  tr <- list(sigma_u = log, sigma_z = log, alpha = log, gamma = stats::qlogis, delta = identity)
  free <- setdiff(names(start), names(fixed))
  make_fld <- function(par) {
    full <- start
    full[free] <- par
    for (nm in names(fixed)) full[nm] <- tr[[nm]](fixed[[nm]])
    c(base, list(sig_u = exp(full[["sigma_u"]]), sig_z = exp(full[["sigma_z"]]),
                 alpha = exp(full[["alpha"]]), gamma = stats::plogis(full[["gamma"]]),
                 delta = full[["delta"]]))
  }
  x0 <- list(beta = S$m_beta, theta = S$m_theta, u = rep(0, graph$nv), s = NULL, z = rep(0, n))
  state <- new.env(); state$x <- x0; state$cache <- list()
  eval_par <- function(par) {
    key <- paste(sprintf("%.6f", par), collapse = ",")
    if (!is.null(state$cache[[key]])) return(state$cache[[key]])
    fld <- make_fld(par)
    res <- tryCatch(odisp_zi_laplace(fld, S, state$x), error = function(e) list(logLik = -Inf))
    if (is.finite(res$logLik)) state$x <- res$x
    state$cache[[key]] <- res
    if (trace) cat(sprintf("  sigma_u %.3f  sigma_z %.3f  alpha %.4g  gamma %.3f  delta %.3f  logLik = %.3f\n",
                           fld$sig_u, fld$sig_z, fld$alpha, fld$gamma, fld$delta, res$logLik))
    res
  }
  par <- start[free]
  if (length(par)) par <- odisp_newton_max(function(p) eval_par(p)$logLik, par)
  res <- eval_par(par)
  fld <- make_fld(par)
  if (trace) cat(sprintf("  converged: logLik = %.3f, %d evaluations\n", res$logLik, length(state$cache)))
  est <- c(res$x$beta, res$x$theta)
  names(est) <- c(paste0("beta.", colnames(S$X)), paste0("theta.", S$nodes_names))
  dimnames(res$Vfix) <- list(names(est), names(est))
  a <- res$a; b <- res$b
  structure(list(coefficients = est, vcov = res$Vfix,
                 hyper = c(sigma_u = fld$sig_u, sigma_z = fld$sig_z, alpha = fld$alpha,
                           gamma = fld$gamma, delta = fld$delta, range_u = range_u),
                 par = par, fixed = fixed, logLik = res$logLik, u = res$x$u, z = res$x$z,
                 E = exp(a), pi = -expm1(-fld$alpha * exp(a + fld$delta * S$lai)),
                 p = -expm1(-exp(fld$gamma * a + b)), setup = S, fld = fld),
            class = "odisp_zifit")
}

# Log-likelihood of the zero-inflated model and its first and second
# derivatives with respect to a = log E and b = d'beta + Z, village by village.
odisp_zi_lik <- function(a, b, S, gamma, alpha, delta = 0) {
  n <- S$n; y <- S$y
  if (any(!is.finite(a)) || any(!is.finite(b))) return(list(ll = -Inf))
  eta <- gamma * a + b
  off <- if (delta != 0) delta * S$lai else 0                       # aridity in pi
  if (max(eta) > 30 || max(a + log(alpha) + off) > 30) return(list(ll = -Inf))
  lam <- exp(eta); kap <- exp(log(alpha) + a + off)
  pi_ <- -expm1(-kap); logpi <- log1mexp(kap)
  la <- lb <- laa <- lab <- lbb <- ll <- numeric(length(y))
  pos <- y > 0
  if (any(pos)) {
    k <- kap[pos]; l <- lam[pos]; nn <- n[pos]; yy <- y[pos]
    logp <- log1mexp(l)
    g <- l * (yy / exp(logp) - nn)                                   # d binomial / d eta
    hneg <- l^2 * yy * exp(-l) / (-expm1(-l))^2 - g                  # - d2 binomial / d eta2
    r <- k / expm1(k)                                                # d log pi / d a
    rp <- k * (expm1(k) - k * exp(k)) / expm1(k)^2                   # d2 log pi / d a2
    rp[!is.finite(rp)] <- 0; r[!is.finite(r)] <- 0
    ll[pos] <- logpi[pos] + yy * logp - (nn - yy) * l + lchoose(nn, yy)
    la[pos] <- r + gamma * g; lb[pos] <- g
    laa[pos] <- rp - gamma^2 * hneg; lab[pos] <- -gamma * hneg; lbb[pos] <- -hneg
  }
  zr <- !pos
  if (any(zr)) {
    k <- kap[zr]; l <- lam[zr]; nn <- n[zr]; p_ <- pi_[zr]
    lA <- -k; lBq <- -nn * l
    lf <- pmax(lA, log(p_) + lBq) + log1p(exp(-abs(lA - (log(p_) + lBq))))   # log f
    wA <- exp(lA - lf); wB <- exp(lBq - lf)                          # A / f, Bq / f
    A <- exp(lA); Bq <- exp(lBq); c2 <- -expm1(lBq)
    Fa <- -wA * c2 * k - p_ * nn * wB * gamma * l
    Fb <- -p_ * nn * wB * l
    Faa <- wA * c2 * k^2 - wA * c2 * k - 2 * nn * gamma * l * k * A * wB +
      nn^2 * gamma^2 * l^2 * p_ * wB - nn * p_ * gamma^2 * l * wB
    Fab <- -nn * l * k * A * wB + nn^2 * p_ * gamma * l^2 * wB - nn * p_ * gamma * l * wB
    Fbb <- nn^2 * p_ * l^2 * wB - nn * p_ * l * wB
    ll[zr] <- lf
    la[zr] <- Fa; lb[zr] <- Fb
    laa[zr] <- Faa - Fa^2; lab[zr] <- Fab - Fa * Fb; lbb[zr] <- Fbb - Fb^2
  }
  tot <- sum(ll)
  list(ll = if (is.finite(tot)) tot else -Inf, ga = la, gb = lb, haa = laa, hab = lab, hbb = lbb)
}

odisp_zi_eval <- function(x, S) {
  w <- exp(as.vector(S$B %*% x$theta) + as.vector(S$A %*% x$u))
  KW <- S$K %*% Matrix::Diagonal(x = w)
  E <- Matrix::rowSums(KW)
  list(a = log(E), b = as.vector(S$X %*% x$beta) + x$z, V = Matrix::Diagonal(x = 1 / E) %*% KW)
}

# Operators at x for the two predictors (a, b): latent vector (beta, theta, u, z).
odisp_zi_ops <- function(x, ev, lk, S, fld) {
  n <- length(S$y); p <- ncol(S$X); q <- ncol(S$B); nu <- length(x$u)
  ib <- seq_len(p); it <- p + seq_len(q); iu <- p + q + seq_len(nu); iz <- p + q + nu + seq_len(n)
  Jt <- as.matrix(ev$V %*% S$B)
  Ju <- ev$V %*% S$A
  LW <- Matrix::solve(fld$Rch, Matrix::solve(fld$Rch, Matrix::t(Ju), system = "P"), system = "L")
  Maa <- S$s2 * tcrossprod(Jt) + fld$sig_u^2 * as.matrix(Matrix::crossprod(LW))
  Mbb <- S$s2 * tcrossprod(S$X) + diag(fld$sig_z^2, n)
  # observed information of (a, b) per village, floored to be positive definite
  Haa <- -lk$haa; Hab <- -lk$hab; Hbb <- -lk$hbb
  tr <- Haa + Hbb; dt <- Haa * Hbb - Hab^2
  disc <- sqrt(pmax((Haa - Hbb)^2 / 4 + Hab^2, 0))
  e1 <- tr / 2 + disc; e2 <- tr / 2 - disc
  fl <- 1e-6
  # eigenvectors of each 2 x 2 block, to rebuild it with floored eigenvalues
  v1a <- ifelse(abs(Hab) > 1e-12, e1 - Hbb, ifelse(Haa >= Hbb, 1, 0))
  v1b <- ifelse(abs(Hab) > 1e-12, Hab, ifelse(Haa >= Hbb, 0, 1))
  nr <- sqrt(v1a^2 + v1b^2); v1a <- v1a / nr; v1b <- v1b / nr
  f1 <- pmax(e1, fl); f2 <- pmax(e2, fl)
  Waa <- f1 * v1a^2 + f2 * v1b^2; Wbb <- f1 * v1b^2 + f2 * v1a^2; Wab <- (f1 - f2) * v1a * v1b
  # Cholesky factor of each block: W = L L', L = [[l11, 0], [l21, l22]]
  l11 <- sqrt(Waa); l21 <- Wab / l11; l22 <- sqrt(pmax(Wbb - l21^2, fl))
  # C = I + L' M L with M = blockdiag(Maa, Mbb) in the order (a, b)
  # L' M L: [ l11 Maa l11 + l21 Mbb l21 , l21 Mbb l22 ; l22 Mbb l21 , l22 Mbb l22 ]
  C11 <- outer(l11, l11) * Maa + outer(l21, l21) * Mbb
  C12 <- outer(l21, l22) * Mbb
  C22 <- outer(l22, l22) * Mbb
  C <- rbind(cbind(C11, C12), cbind(t(C12), C22)) + diag(2 * n)
  cC <- chol(C)
  Jv <- function(v) list(a = as.vector(Jt %*% v[it] + Ju %*% v[iu]),
                         b = as.vector(S$X %*% v[ib]) + v[iz])
  Jtw <- function(wa, wb) c(as.vector(crossprod(S$X, wb)), as.vector(crossprod(Jt, wa)),
                            as.vector(Matrix::crossprod(Ju, wa)), wb)
  Pinv <- function(r) c(S$s2 * r[ib], S$s2 * r[it],
                        fld$sig_u^2 * as.vector(Matrix::solve(fld$Rch, r[iu])), fld$sig_z^2 * r[iz])
  Pv <- function(v) c(v[ib] / S$s2, v[it] / S$s2, as.vector(fld$R %*% v[iu]) / fld$sig_u^2,
                      v[iz] / fld$sig_z^2)
  HGNinv <- function(r) {
    s0 <- Pinv(r); jv <- Jv(s0)
    t1 <- l11 * jv$a + l21 * jv$b; t2 <- l22 * jv$b                        # L' J s0
    w <- backsolve(cC, backsolve(cC, c(t1, t2), transpose = TRUE))
    w1 <- w[seq_len(n)]; w2 <- w[n + seq_len(n)]
    s0 - Pinv(Jtw(l11 * w1, l21 * w1 + l22 * w2))                          # J' L w
  }
  cc <- as.vector(Matrix::crossprod(ev$V, lk$ga))
  Cv <- function(v) {
    nuv <- as.vector(S$B %*% v[it]) + as.vector(S$A %*% v[iu])
    cx <- cc * nuv - as.vector(Matrix::crossprod(ev$V, lk$ga * as.vector(ev$V %*% nuv)))
    out <- numeric(length(v))
    out[it] <- as.vector(crossprod(S$B, cx)); out[iu] <- as.vector(Matrix::crossprod(S$A, cx))
    out
  }
  Hv <- function(v) {
    jv <- Jv(v)
    Pv(v) + Jtw(Haa * jv$a + Hab * jv$b, Hab * jv$a + Hbb * jv$b) - Cv(v)
  }
  list(HGNinv = HGNinv, Hv = Hv, Jtw = Jtw, Pv = Pv, Jt = Jt, Ju = Ju, cC = cC,
       l11 = l11, l21 = l21, l22 = l22)
}

odisp_zi_laplace <- function(fld, S, x, maxit = 50, tol = 1e-7) {
  m <- odisp_u_vec(list(beta = S$m_beta, theta = S$m_theta, u = 0 * x$u, s = NULL, z = 0 * x$z))
  prior <- function(x) odisp_u_prior(x, S, fld)
  f_of <- function(x) {
    ev <- odisp_zi_eval(x, S)
    lk <- odisp_zi_lik(ev$a, ev$b, S, fld$gamma, fld$alpha, if (is.null(fld$delta)) 0 else fld$delta)
    v <- lk$ll - 0.5 * prior(x)
    list(f = if (is.finite(v)) v else -Inf, ev = ev, lk = lk)
  }
  fx <- f_of(x)
  if (!is.finite(fx$f)) {                           # restart from the prior centre
    x <- list(beta = S$m_beta, theta = S$m_theta, u = 0 * x$u, s = NULL, z = 0 * x$z); fx <- f_of(x)
  }
  if (!is.finite(fx$f)) stop("non-finite objective at the starting point")
  for (it in seq_len(maxit)) {
    ops <- odisp_zi_ops(x, fx$ev, fx$lk, S, fld)
    xv <- odisp_u_vec(x)
    gv <- ops$Jtw(fx$lk$ga, fx$lk$gb) - ops$Pv(xv - m)
    pc <- odisp_u_pcg(gv, ops)
    slope <- sum(gv * pc$d)
    if (!is.finite(slope) || slope <= 0) { pc$d <- ops$HGNinv(gv); slope <- sum(gv * pc$d) }
    if (slope < tol) break
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
  ops <- odisp_zi_ops(x, fx$ev, fx$lk, S, fld)
  logdet <- 2 * sum(log(diag(ops$cC)))
  # covariance of (beta, theta) from the floored Hessian, by Woodbury
  p <- ncol(S$X); q <- ncol(S$B)
  Vf <- sapply(seq_len(p + q), function(k) {
    e <- numeric(length(odisp_u_vec(x))); e[k] <- 1
    ops$HGNinv(e)[seq_len(p + q)]
  })
  list(logLik = fx$f - 0.5 * logdet, x = x, a = fx$ev$a, b = fx$ev$b, iterations = it,
       Vfix = (Vf + t(Vf)) / 2)
}

#' @export
print.odisp_zifit <- function(x, ...) {
  h <- x$hyper
  cat("Zero-inflated river-network model, range of dispersal", x$setup$phi, "km\n")
  dl <- if ("delta" %in% names(h)) h[["delta"]] else 0
  cat(sprintf("pi(x) = 1 - exp(-alpha E(x) AI(x)^delta), alpha = %.4g, delta = %.3f;  lambda(x) = E(x)^gamma T(x), gamma = %.3f\n",
              h[["alpha"]], dl, h[["gamma"]]))
  cat(sprintf("Field u: range %g km, sd %.3f; village effect sd %.3f\n", h[["range_u"]], h[["sigma_u"]], h[["sigma_z"]]))
  print(round(cbind(estimate = x$coefficients, se = sqrt(diag(x$vcov))), 3))
  cat("Laplace log-likelihood:", round(x$logLik, 2), "\n")
  invisible(x)
}
