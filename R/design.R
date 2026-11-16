# Build covariates without dropping rows, and retain their training encoding.
# Saved terms also retain prediction information for bases such as poly().
odisp_model_matrix <- function(formula, data, spec = NULL, drop_intercept = FALSE) {
  if (!is.null(spec) && nrow(data) == 0L) {
    X <- matrix(numeric(0), 0L, length(spec$columns), dimnames = list(NULL, spec$columns))
    return(list(matrix = X, spec = spec))
  }
  tt <- if (is.null(spec)) formula else spec$terms
  if (!is.null(spec)) {
    # model.frame resets levels; the saved contrasts are restored below.
    # Remove input contrast attributes to avoid a misleading "dropped" warning.
    for (nm in intersect(names(spec$xlevels), names(data))) {
      if (is.factor(data[[nm]])) attr(data[[nm]], "contrasts") <- NULL
    }
  }
  mf <- stats::model.frame(tt, data, na.action = stats::na.fail,
                           xlev = if (!is.null(spec)) spec$xlevels else NULL)
  tt <- attr(mf, "terms")
  X <- stats::model.matrix(tt, mf, contrasts.arg = if (!is.null(spec)) spec$contrasts else NULL)
  contrasts <- attr(X, "contrasts")
  if (drop_intercept) X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (nrow(X) != nrow(data) || any(!is.finite(X)))
    stop("covariates must be finite at every row")
  if (is.null(spec)) {
    spec <- list(terms = stats::delete.response(tt), xlevels = stats::.getXlevels(tt, mf),
                 contrasts = contrasts, columns = colnames(X))
  }
  list(matrix = X, frame = mf, spec = spec)
}

odisp_prediction_matrix <- function(formula, data, spec, columns, drop_intercept = FALSE) {
  out <- odisp_model_matrix(stats::delete.response(stats::terms(formula)), data,
                            spec, drop_intercept)
  # Older saved fits lack factor levels and contrasts. Numeric models remain
  # usable, but reconstructing a categorical encoding would risk wrong answers.
  if (is.null(spec) && length(out$spec$xlevels))
    stop("this saved fit lacks factor encoding; refit before predicting categorical covariates")
  if (!identical(colnames(out$matrix), columns))
    stop("prediction covariates do not match those of the fit")
  out$matrix
}

odisp_prediction_pairs <- function(newdata, nodes, pairs, range) {
  if (any(!is.finite(range)) || any(range <= 0))
    stop("prediction ranges must be positive and finite")
  required <- max(range)
  if (is.null(pairs)) return(river_pairs(newdata, nodes, required))
  if (pairs$m != nrow(newdata)) stop("'pairs' was computed for a different set of locations")
  # exp(log(phi)) can exceed phi by a few floating-point rounding units.
  if (required > pairs$max_range * (1 + 8 * .Machine$double.eps))
    stop("the max_range of 'pairs' is too small for prediction; recompute pairs with max_range >= ", required)
  pairs
}

# The fit and prediction must use the same smooth upper bound on range.
odisp_range <- function(eta0, gamma, logAI, phi_max = NULL) {
  g_lin <- eta0 + gamma * logAI
  if (is.null(phi_max)) return(list(phi = exp(g_lin), dlogphi = rep(1, length(g_lin))))
  zz <- log(phi_max) - g_lin
  sp <- ifelse(zz > 30, zz, log1p(exp(pmin(zz, 30))))
  list(phi = phi_max * exp(-sp), dlogphi = stats::plogis(zz))
}
