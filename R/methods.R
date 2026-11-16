#' @export
coef.odisp_fit <- function(object, ...) object$coefficients

#' @export
vcov.odisp_fit <- function(object, ...) object$vcov

#' @export
logLik.odisp_fit <- function(object, ...) {
  structure(object$logLik, df = object$df, nobs = object$nobs, class = "logLik")
}

#' @export
fitted.odisp_fit <- function(object, ...) object$fitted

#' @export
nobs.odisp_fit <- function(object, ...) object$nobs

#' @export
print.odisp_fit <- function(x, ...) {
  cat("River-network dispersal model (no latent fields)\n")
  cat("Log-likelihood:", format(x$logLik, digits = 6), " AIC:", format(stats::AIC(x), digits = 6), "\n\n")
  print(round(x$coefficients, 4))
  cat("\nDispersal range (km) at the survey locations:",
      sprintf("min %.1f, median %.1f, max %.1f", min(x$phi), stats::median(x$phi), max(x$phi)), "\n")
  invisible(x)
}

#' @export
summary.odisp_fit <- function(object, ...) {
  est <- object$coefficients
  se <- rep(NA_real_, length(est)); names(se) <- names(est)
  se[names(object$free)[object$free]] <- sqrt(pmax(diag(object$vcov), 0))
  tab <- cbind(Estimate = est, `Std. Error` = se, `z value` = est / se,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(est / se)))
  tab[!object$free, 2:4] <- NA
  structure(list(table = tab, fixed = object$fixed, logLik = object$logLik,
                 AIC = stats::AIC(object), phi = object$phi,
                 n_zero_exposure = sum(object$exposure == 0),
                 m = object$nobs, background = object$background,
                 phi_max = object$design$phi_max),
            class = "summary.odisp_fit")
}

#' @export
print.summary.odisp_fit <- function(x, ...) {
  cat("River-network dispersal model (no latent fields)\n\n")
  stats::printCoefmat(x$table, na.print = "fixed", P.values = TRUE, has.Pvalue = TRUE)
  cat("\nLog-likelihood:", format(x$logLik, digits = 6), "  AIC:", format(x$AIC, digits = 6), "\n")
  cat(sprintf("Dispersal range at survey locations: min %.1f km, median %.1f km, max %.1f km\n",
              min(x$phi), stats::median(x$phi), max(x$phi)))
  if (!is.null(x$phi_max))
    cat(sprintf("Range bounded above by %g km; %d of %d locations within 10%% of the bound\n",
                x$phi_max, sum(x$phi > 0.9 * x$phi_max), x$m))
  cat(sprintf("Locations with zero exposure: %d of %d%s\n", x$n_zero_exposure, x$m,
              if (x$background) " (before the background term)" else ""))
  invisible(x)
}

#' Predict exposure and prevalence at new locations
#'
#' @param object An \code{odisp_fit} object.
#' @param newdata An \code{sf} object of points with the aridity column and
#'   the prevalence covariates used in the fit.
#' @param nodes The \code{odisp_network} object used in the fit.
#' @param pairs Optional \code{odisp_pairs} for \code{newdata}; computed if
#'   \code{NULL}, using the largest dispersal range at the new locations.
#'   Supplied pairs must cover that range.
#' @param ... Not used.
#' @return A data frame with the dispersal range \code{phi} (km), the
#'   exposure \code{E}, the force of infection \code{lambda}, the prevalence
#'   \code{prevalence}, and \code{zero_exposure}, which is \code{TRUE} where no
#'   river node is within reach.
#' @export
predict.odisp_fit <- function(object, newdata, nodes, pairs = NULL, ...) {
  d0 <- object$design
  df <- sf::st_drop_geometry(newdata)
  ai <- df[[d0$aridity]]
  if (is.null(ai) || any(!is.finite(ai)) || any(ai <= 0))
    stop("aridity must be positive and finite at every location")
  phi <- odisp_range(object$coefficients[["eta0"]], object$coefficients[["gamma"]],
                      log(ai), d0$phi_max)$phi
  pairs <- odisp_prediction_pairs(newdata, nodes, pairs, phi)
  used <- sort(unique(pairs$j))
  # Prediction has no response or trial counts; reuse only the fitted design
  # conventions and replace its covariates and location-node pairs.
  des <- d0
  des$D <- odisp_prediction_matrix(d0$formula, df, d0$prevalence_spec, colnames(d0$D))
  des$Z <- odisp_prediction_matrix(d0$breeding, as.data.frame(nodes)[used, , drop = FALSE],
                                    d0$breeding_spec, colnames(d0$Z), drop_intercept = TRUE)
  ord <- order(pairs$i, pairs$r)
  des$m <- nrow(df); des$logAI <- log(ai)
  des$i <- pairs$i[ord]; des$j <- match(pairs$j[ord], used); des$r <- pairs$r[ord]
  des$start <- as.integer(c(0L, cumsum(tabulate(des$i, nbins = des$m))))
  des$len <- nodes$len[used]; des$max_range <- pairs$max_range
  des$y <- des$n <- NULL
  ev <- odisp_eval(object$coefficients, des, object$background)
  data.frame(phi = ev$phi, E = ev$E, lambda = ev$lambda,
             prevalence = -expm1(-ev$lambda), zero_exposure = ev$E == 0)
}

#' Simulate survey outcomes from a dispersal model
#'
#' @param design An \code{odisp_design} object.
#' @param par Named parameter vector, with the names produced by
#'   \code{coef()} on a fit using the same design.
#' @param background Logical; whether \code{par} contains \code{log_eps}.
#' @param nsim Number of simulated data sets.
#' @return A matrix with one column of simulated positive counts per data set.
#' @export
odisp_simulate <- function(design, par, background = FALSE, nsim = 1) {
  nm <- odisp_par_names(design, background)
  miss <- setdiff(nm, names(par))
  if (length(miss)) stop("missing parameter(s): ", paste(miss, collapse = ", "))
  p <- -expm1(-odisp_eval(par[nm], design, background)$lambda)
  matrix(stats::rbinom(design$m * nsim, design$n, p), design$m, nsim)
}

#' @export
simulate.odisp_fit <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  odisp_simulate(object$design, object$coefficients, object$background, nsim)
}

#' Profile log-likelihood for one parameter
#'
#' Refits the model with \code{param} fixed at each of \code{values}, keeping
#' any parameters fixed in the original fit, and returns the maximised
#' log-likelihood. Useful for checking whether the dispersal range is
#' identified, e.g. \code{odisp_profile(fit, "eta0", log(c(5, 10, 20, 40)))}.
#'
#' @param fit An \code{odisp_fit} object.
#' @param param Name of the parameter to profile.
#' @param values Values at which to fix it.
#' @return A data frame with the values, the profile log-likelihood and
#'   \code{range_capped}, which is \code{TRUE} where the fitted range exceeds
#'   \code{max_range}; log-likelihoods for those rows are not reliable.
#' @export
odisp_profile <- function(fit, param, values) {
  if (!param %in% names(coef(fit))) stop("unknown parameter: ", param)
  st <- coef(fit)
  res <- lapply(values, function(val) {
    fx <- utils::modifyList(fit$fixed, stats::setNames(list(val), param))
    f <- tryCatch(suppressWarnings(odisp_fit(design = fit$design, background = fit$background,
                                             fixed = fx, start = st)),
                  error = function(e) NULL)
    if (is.null(f)) c(NA_real_, NA) else c(f$logLik, f$range_capped)
  })
  res <- do.call(rbind, res)
  data.frame(value = values, logLik = res[, 1], range_capped = as.logical(res[, 2]))
}
