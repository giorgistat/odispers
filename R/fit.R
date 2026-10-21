#' Build the design for a river-network dispersal model
#'
#' Collects everything the likelihood needs: the outcome, the covariates for
#' prevalence and for breeding, the aridity at each location, and the
#' location-node pairs with their distances. Usually called through
#' \code{\link{odisp_fit}}; it is exported so that the same design can be
#' used for simulation with \code{\link{odisp_simulate}}.
#'
#' @param formula Formula for the prevalence part, with the number of positive
#'   individuals on the left and the covariates \eqn{d(x)} on the right, for
#'   example \code{y ~ 1}. The intercept of the model belongs here.
#' @param breeding One-sided formula for the breeding intensity on the network,
#'   evaluated on the columns of \code{nodes}, for example
#'   \code{~ log(DIS_AV_CMS)}. Any intercept is dropped, since it is not
#'   identifiable separately from the intercept in \code{formula}. Use
#'   \code{~ 0} or \code{~ 1} for constant breeding along the network.
#' @param data An \code{sf} object of survey locations containing the outcome,
#'   the number examined, the aridity index and the prevalence covariates.
#' @param nodes An \code{odisp_network} object from \code{\link{river_nodes}}.
#' @param pairs An \code{odisp_pairs} object from \code{\link{river_pairs}}
#'   computed for \code{data}; computed if \code{NULL}.
#' @param trials Name of the column in \code{data} with the number examined.
#' @param aridity Name of the column in \code{data} with the aridity index
#'   (precipitation / evapotranspiration, positive values).
#' @param max_range Maximum dispersal range (km) used when \code{pairs} is
#'   computed here.
#' @return An object of class \code{odisp_design}.
#' @export
odisp_design <- function(formula, breeding = ~ 0, data, nodes, pairs = NULL,
                         trials = "n", aridity = "AI", max_range = 50) {
  if (!inherits(data, "sf")) stop("'data' must be an sf object of survey locations")
  if (is.null(pairs)) pairs <- river_pairs(data, nodes, max_range)
  if (pairs$m != nrow(data)) stop("'pairs' was computed for a different set of locations")
  df <- sf::st_drop_geometry(data)

  mf <- stats::model.frame(formula, df, na.action = stats::na.fail)
  y <- stats::model.response(mf)
  D <- stats::model.matrix(formula, mf)
  n <- df[[trials]]
  if (is.null(n)) stop("column '", trials, "' not found in data")
  if (any(y > n) || any(y < 0)) stop("counts must satisfy 0 <= y <= n")
  ai <- df[[aridity]]
  if (is.null(ai)) stop("column '", aridity, "' not found in data")
  if (any(!is.finite(ai)) || any(ai <= 0)) stop("aridity must be positive and finite at every location")

  # breeding covariates, only for nodes that are within range of some location
  used <- sort(unique(pairs$j))
  Z <- stats::model.matrix(breeding, as.data.frame(nodes)[used, , drop = FALSE])
  Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]
  jmap <- match(pairs$j, used)

  structure(list(y = as.numeric(y), n = as.numeric(n), D = D, Z = Z,
                 logAI = log(ai), i = pairs$i, j = jmap, r = pairs$r,
                 len = nodes$len[used], m = nrow(df), max_range = pairs$max_range,
                 formula = formula, breeding = breeding, trials = trials,
                 aridity = aridity),
            class = "odisp_design")
}

# names of all parameters in the order used internally
odisp_par_names <- function(design, background) {
  c(paste0("beta.", colnames(design$D)),
    if (ncol(design$Z)) paste0("theta.", colnames(design$Z)),
    "eta0", "gamma",
    if (background) "log_eps")
}

# Exposure, its derivatives and the force of infection for one parameter vector
odisp_eval <- function(par, design, background, deriv = FALSE) {
  p <- ncol(design$D); q <- ncol(design$Z)
  beta <- par[seq_len(p)]
  theta <- par[p + seq_len(q)]
  eta0 <- par[p + q + 1]; gamma <- par[p + q + 2]
  eps <- if (background) exp(par[p + q + 3]) else 0

  phi <- exp(eta0 + gamma * design$logAI)
  i <- design$i; j <- design$j
  phi_i <- phi[i]
  t <- design$r / phi_i
  w <- if (q) as.vector(exp(design$Z %*% theta)) else rep(1, length(design$len))
  base <- design$len[j] * w[j] / phi_i          # Delta_j w_j / phi_i
  inside <- t < 1
  kt <- numeric(length(t)); kt[inside] <- (1 - t[inside])^4 * (1 + 4 * t[inside])
  contrib <- base * kt

  sum_by <- function(v) { out <- numeric(design$m); s <- rowsum(v, i, reorder = FALSE)
                          out[as.integer(rownames(s))] <- s[, 1]; out }
  E <- sum_by(contrib)
  h <- as.vector(exp(design$D %*% beta))
  lambda <- (E + eps) * h
  res <- list(phi = phi, E = E, h = h, lambda = lambda, eps = eps)
  if (!deriv) return(res)

  # dE/deta0 = - sum Delta w phi^-1 [k(t) + t k'(t)],  k + t k' = (1-t)^3 (1 + 3t - 24 t^2)
  g <- numeric(length(t))
  ti <- t[inside]
  g[inside] <- (1 - ti)^3 * (1 + 3 * ti - 24 * ti^2)
  res$dE_deta0 <- -sum_by(base * g)
  res$dE_dgamma <- res$dE_deta0 * design$logAI
  if (q) {
    res$dE_dtheta <- vapply(seq_len(q), function(k) sum_by(contrib * design$Z[j, k]),
                            numeric(design$m))
    dim(res$dE_dtheta) <- c(design$m, q)
  }
  res
}

# negative log-likelihood and gradient, as functions of the free parameters
odisp_objective <- function(design, background, fixed) {
  nm <- odisp_par_names(design, background)
  free <- stats::setNames(!(nm %in% names(fixed)), nm)
  full <- function(x) { v <- numeric(length(nm)); names(v) <- nm
                        v[free] <- x; v[names(fixed)] <- unlist(fixed); v }
  y <- design$y; n <- design$n; pos <- y > 0
  big <- 1e12
  penalty <- 1e4   # keeps phi below max_range, where the precomputed pairs are complete

  fn <- function(x) {
    ev <- odisp_eval(full(x), design, background)
    lam <- ev$lambda
    if (any(!is.finite(lam)) || any(lam[pos] <= 0)) return(big)
    ll <- -(n - y) * lam
    ll[pos] <- ll[pos] + y[pos] * log1mexp(lam[pos])
    excess <- pmax(0, log(ev$phi / design$max_range))
    out <- -sum(ll) + penalty * sum(excess^2)
    if (!is.finite(out)) big else out
  }
  gr <- function(x) {
    par <- full(x)
    ev <- odisp_eval(par, design, background, deriv = TRUE)
    lam <- ev$lambda
    if (any(!is.finite(lam)) || any(lam[pos] <= 0)) return(rep(0, sum(free)))
    dl <- -(n - y)                                 # d loglik / d lambda
    dl[pos] <- dl[pos] + y[pos] / expm1(lam[pos])
    p <- ncol(design$D); q <- ncol(design$Z)
    g <- c(colSums(dl * lam * design$D),
           if (q) colSums(dl * ev$h * ev$dE_dtheta),
           sum(dl * ev$h * ev$dE_deta0),
           sum(dl * ev$h * ev$dE_dgamma),
           if (background) sum(dl * ev$h * ev$eps))
    excess <- pmax(0, log(ev$phi / design$max_range))
    pen_eta <- 2 * penalty * sum(excess)
    pen_gam <- 2 * penalty * sum(excess * design$logAI)
    k_eta <- ncol(design$D) + ncol(design$Z) + 1
    g[k_eta] <- g[k_eta] - pen_eta
    g[k_eta + 1] <- g[k_eta + 1] - pen_gam
    -unname(g[free])
  }
  list(fn = fn, gr = gr, full = full, free = free, names = nm)
}

#' Fit a river-network dispersal model by maximum likelihood
#'
#' Fits the model without latent random fields:
#' \deqn{w(s) = \exp\{z(s)^\top\theta\}, \qquad
#'       E(x) = \int_\Gamma \phi(x)^{-1} k(\|x - s\|/\phi(x))\, w(s)\, ds,}
#' \deqn{\log\phi(x) = \eta_0 + \gamma \log AI(x), \qquad
#'       Y_i \sim Bin(n_i, 1 - \exp[-\{E(x_i) + \epsilon\} \exp\{d(x_i)^\top\beta\}]),}
#' where \eqn{k} is the Wendland kernel and \eqn{\epsilon \ge 0} an optional
#' background exposure. Without a background term, prevalence is exactly zero
#' at locations further than \eqn{\phi(x)} from every river node.
#'
#' @inheritParams odisp_design
#' @param background Logical; if \code{TRUE}, a background exposure
#'   \eqn{\epsilon > 0} is estimated (as \code{log_eps}). This prevents a
#'   single positive survey far from any river from forcing a large range.
#' @param fixed Named list of parameters held fixed, for example
#'   \code{list(gamma = 0)} for a range that does not depend on aridity.
#'   Parameter names are those returned by \code{coef()}.
#' @param start Optional named vector of starting values.
#' @param gamma_starts Starting values for \code{gamma} when neither
#'   \code{start} nor \code{fixed} sets it; the optimisation is run from each
#'   and the best solution is kept, which guards against local optima.
#' @param control List passed to \code{\link[stats]{optim}}.
#' @param design Optional pre-built \code{odisp_design}; if given, the
#'   arguments used to build it are ignored.
#' @return An object of class \code{odisp_fit}.
#' @export
odisp_fit <- function(formula, breeding = ~ 0, data, nodes, pairs = NULL,
                      trials = "n", aridity = "AI", max_range = 50,
                      background = FALSE, fixed = list(), start = NULL,
                      gamma_starts = c(0, 0.5, 1), control = list(maxit = 500),
                      design = NULL) {
  cl <- match.call()
  if (is.null(design))
    design <- odisp_design(formula, breeding, data, nodes, pairs, trials, aridity, max_range)
  obj <- odisp_objective(design, background, fixed)
  nm <- obj$names
  if (length(setdiff(names(fixed), nm)))
    stop("unknown fixed parameter(s): ", paste(setdiff(names(fixed), nm), collapse = ", "))

  # starting values; the range starts large enough for every location with
  # positives to reach a river node
  dmin <- rep(Inf, design$m)
  dm <- tapply(design$r, design$i, min)
  dmin[as.integer(names(dm))] <- dm
  pos_reach <- design$y > 0 & is.finite(dmin)
  pbar <- sum(design$y) / sum(design$n)
  set_intercept <- function(st) {
    if ("beta.(Intercept)" %in% nm && !("beta.(Intercept)" %in% names(fixed))) {
      ev <- odisp_eval(st, design, background); reach <- ev$E + ev$eps > 0
      if (any(reach)) st["beta.(Intercept)"] <- log(-log(1 - pbar) / mean(ev$E[reach] + ev$eps))
    }
    st
  }
  make_start <- function(g0) {
    st <- stats::setNames(numeric(length(nm)), nm)
    st["gamma"] <- g0
    st["eta0"] <- log(0.5 * design$max_range) - g0 * mean(design$logAI)
    if (any(pos_reach))
      st["eta0"] <- max(st["eta0"], max(log(1.05 * dmin[pos_reach]) - g0 * design$logAI[pos_reach]))
    if (background) st["log_eps"] <- log(1e-3)
    st[names(fixed)] <- unlist(fixed)
    set_intercept(st)
  }
  if (!is.null(start)) {
    st <- make_start(if ("gamma" %in% names(fixed)) fixed$gamma
                     else if ("gamma" %in% names(start)) start[["gamma"]] else 0)
    keep <- intersect(names(start), setdiff(nm, names(fixed)))
    st[keep] <- start[keep]
    if (obj$fn(st[obj$free]) >= 1e12) st <- make_start(st[["gamma"]])
    starts <- list(st)
  } else if ("gamma" %in% names(fixed)) {
    starts <- list(make_start(fixed$gamma))
  } else {
    starts <- lapply(gamma_starts, make_start)
  }
  ok <- vapply(starts, function(st) obj$fn(st[obj$free]) < 1e12, logical(1))
  if (!any(ok)) {
    st <- starts[[1]]
    bad <- which(design$y > 0 & odisp_eval(st, design, background)$lambda <= 0)
    stop(length(bad), " location(s) with positive counts have no river node within max_range. ",
         "Increase max_range or set background = TRUE. Locations: ",
         paste(utils::head(bad, 10), collapse = ", "))
  }
  opts <- lapply(starts[ok], function(st)
    stats::optim(st[obj$free], obj$fn, obj$gr, method = "BFGS",
                 control = utils::modifyList(list(maxit = 500), control)))
  opt <- opts[[which.min(vapply(opts, `[[`, numeric(1), "value"))]]
  H <- stats::optimHess(opt$par, obj$fn, obj$gr)
  V <- tryCatch(solve(H), error = function(e) matrix(NA_real_, nrow(H), ncol(H)))
  dimnames(V) <- list(names(opt$par), names(opt$par))
  est <- obj$full(opt$par)
  ev <- odisp_eval(est, design, background)

  cap_binds <- max(ev$phi) > design$max_range
  if (cap_binds)
    warning(sprintf("fitted dispersal range reaches %.1f km at some locations, beyond max_range = %g km; ",
                    max(ev$phi), design$max_range),
            "exposure there is truncated. Refit with a larger max_range.")
  if (opt$convergence != 0) warning("optim did not report convergence (code ", opt$convergence, ")")

  ev_ll <- -expm1(-ev$lambda)
  ll <- sum(lgamma(design$n + 1) - lgamma(design$y + 1) - lgamma(design$n - design$y + 1)) +
    sum(ifelse(design$y > 0, design$y * log(ev_ll), 0)) - sum((design$n - design$y) * ev$lambda)
  structure(list(coefficients = est, vcov = V, free = obj$free, fixed = fixed,
                 logLik = ll, df = sum(obj$free), nobs = design$m,
                 phi = ev$phi, exposure = ev$E, lambda = ev$lambda,
                 fitted = -expm1(-ev$lambda), background = background,
                 range_capped = cap_binds,
                 design = design, optim = opt, call = cl),
            class = "odisp_fit")
}
