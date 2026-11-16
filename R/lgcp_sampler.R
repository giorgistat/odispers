# Seed a feasible joint point pattern, sharing a site whenever it already
# reaches a positive village. Transposing once lets us visit only the cells
# in reach, rather than materializing a dense row for every positive village.
lgcp_initial_counts <- function(K, y) {
  village_cells <- methods::as(Matrix::t(K), "CsparseMatrix")
  counts <- integer(ncol(K))
  for (village in which(y > 0)) {
    first <- village_cells@p[village] + 1L
    last <- village_cells@p[village + 1L]
    entries <- if (first <= last) seq.int(first, last) else integer()
    entries <- entries[village_cells@x[entries] > 0]
    if (!length(entries)) {
      stop("a village with positives has no positive kernel weight; ",
           "check the kernel width and river-cell resolution")
    }
    cells <- village_cells@i[entries] + 1L
    if (!any(counts[cells] > 0L)) {
      strongest <- which.max(village_cells@x[entries])
      counts[cells[strongest]] <- 1L
    }
  }
  counts
}

# A single vector contains the breeding coefficients followed by the graph
# field. Given cell counts, prevalence drops out of this Gaussian block.
lgcp_poisson_loglik <- function(x, counts, S) {
  eta <- S$log_mass + as.numeric(S$B %*% x)
  sum(counts * eta - exp(eta))
}

lgcp_field_log_density <- function(x, counts, S) {
  lgcp_poisson_loglik(x, counts, S) - sum(x * as.numeric(S$Q %*% x)) / 2
}

lgcp_field_laplace <- function(counts, S) {
  # Starting at the prior mean makes this proposal a function of counts alone.
  x <- numeric(ncol(S$B))
  for (iter in seq_len(80L)) {
    mu <- exp(S$log_mass + as.numeric(S$B %*% x))
    H <- Matrix::forceSymmetric(S$Q + Matrix::crossprod(
      S$B, Matrix::Diagonal(x = mu) %*% S$B))
    ch <- Matrix::Cholesky(H, LDL = FALSE)
    score <- Matrix::crossprod(S$B, counts - mu) - S$Q %*% x
    step <- as.numeric(Matrix::solve(ch, score))
    gain <- sum(as.numeric(score) * step)
    if (gain <= 1e-12) return(list(mode = x, H = H, ch = ch))

    # Evaluate the log density *difference* along this Newton direction.
    # Subtracting two full log densities can stall on large graphs.
    projected <- as.numeric(S$B %*% step)
    prior_curve <- sum(step * as.numeric(S$Q %*% step))
    scale <- 1
    repeat {
      d <- scale * projected
      improvement <- scale * gain - scale^2 * prior_curve / 2 - sum(mu * (expm1(d) - d))
      if (is.finite(improvement) && improvement >= 1e-4 * scale * gain) break
      scale <- scale / 2
      if (scale < 1e-12) stop("LGCP conditional mode line search failed")
    }
    x <- x + scale * step
  }
  stop("LGCP conditional mode did not converge")
}

lgcp_gaussian_draw <- function(ch, n) {
  # Matrix factors the permuted precision as L L'. Apply P' after L'^(-1).
  innovation <- Matrix::solve(ch, stats::rnorm(n), system = "Lt")
  as.numeric(Matrix::solve(ch, innovation, system = "Pt"))
}

lgcp_field_step <- function(x, counts, S, laplace) {
  candidate <- laplace$mode + lgcp_gaussian_draw(laplace$ch, length(x))
  log_weight <- function(v) {
    delta <- v - laplace$mode
    lgcp_field_log_density(v, counts, S) + sum(delta * as.numeric(laplace$H %*% delta)) / 2
  }
  diff <- log_weight(candidate) - log_weight(x)
  accept <- is.finite(diff) && log(stats::runif(1L)) < diff
  list(x = if (accept) candidate else x, accepted = accept)
}

# A force of exactly zero is valid for villages with y = 0. Use the log
# survival directly for failures even when infection probability rounds to 1.
lgcp_survey <- function(exposure, predictor, y, n, link = "cloglog") {
  log_force <- rep(-Inf, length(y))
  reached <- exposure > 0
  log_force[reached] <- log(exposure[reached]) + predictor[reached]
  cases <- y > 0
  negatives <- n > y
  if (identical(link, "logit")) {
    log_success <- stats::plogis(log_force[cases], log.p = TRUE)
    log_failure <- stats::plogis(log_force[negatives], lower.tail = FALSE, log.p = TRUE)
    return(list(probability = stats::plogis(log_force),
                loglik = sum(y[cases] * log_success) + sum((n-y)[negatives] * log_failure)))
  }
  force <- exp(log_force)
  log_prob <- log1mexp(force[cases])
  # Below machine precision log(1-exp(-force)) equals log(force) to
  # working precision, including when the positive force underflows.
  tiny <- log_force[cases] < log(.Machine$double.eps)
  log_prob[tiny] <- log_force[cases][tiny]
  list(probability = -expm1(-force),
       loglik = sum(y[cases] * log_prob) - sum((n - y)[negatives] * force[negatives]))
}

# Below a positive village must be active. For an observed zero, Bayes' rule
# adds the active-state log probability of n negative tests to the prior log odds.
lgcp_activity_probability <- function(eta, exposure, predictor, y, n, link) {
  zero <- y == 0 & n > 0
  log_force <- log(exposure[zero]) + predictor[zero]
  log_failure <- if (link == "logit") stats::plogis(log_force, lower.tail=FALSE, log.p=TRUE) else -exp(log_force)
  eta[zero] <- eta[zero] + n[zero] * log_failure
  probability <- stats::plogis(eta)
  probability[y > 0] <- 1
  probability
}

lgcp_count_sweep <- function(counts, mu, predictor, S, block_size, high = integer(), spread = 0,
                             relocation_group = integer()) {
  lgcp_local_sweep(S$K, S$kernel_tree, counts, high, mu, predictor, S$y, S$n,
                   sample.int(length(counts)), min(block_size, length(counts)), spread, relocation_group,
                   identical(S$link, "logit"), as.numeric(S$gate_base),
                   as.integer(S$activity_sign), if(is.null(S$gate_slope)) 0 else S$gate_slope)
}

# Uniform-prior slice on [0,1]. Conditional on the two type counts the survey
# exposure is affine in the mark spread. Positive summands preserve tiny risks.
lgcp_mark_step <- function(counts, high, spread, predictor, S) {
  low_E <- as.numeric(S$K %*% (counts - high))
  high_E <- as.numeric(S$K %*% high)
  evaluate <- function(d) (1 - d) * low_E + (1 + d) * high_E
  loglik <- function(E) lgcp_exposure_loglik(E, predictor, S)
  threshold <- loglik(evaluate(spread)) + log(stats::runif(1))
  lo <- 0
  hi <- 1
  repeat {
    d <- stats::runif(1, lo, hi)
    E <- evaluate(d)
    if (loglik(E) >= threshold) return(list(spread = d, exposure = E))
    if (d < spread) lo <- d else hi <- d
  }
}

# Elliptical slice sampling for a centred Gaussian prior. The independent
# prior draw defines the ellipse and only the conditional likelihood is used.
lgcp_elliptical_slice <- function(x, direction, loglik) {
  threshold <- loglik(x) + log(stats::runif(1))
  force(direction)
  angle <- stats::runif(1, 0, 2 * pi)
  lo <- angle - 2 * pi
  hi <- angle
  repeat {
    candidate <- x * cos(angle) + direction * sin(angle)
    if (loglik(candidate) >= threshold) return(candidate)
    if (angle < 0) lo <- angle else hi <- angle
    angle <- stats::runif(1, lo, hi)
  }
}

lgcp_beta_step <- function(beta, exposure, S) {
  if (!length(beta)) return(beta)
  loglik <- function(b) {
    predictor <- S$survey_offset + as.numeric(S$X %*% b)
    lgcp_survey(exposure, predictor, S$y, S$n, S$link)$loglik
  }
  lgcp_elliptical_slice(beta, stats::rnorm(length(beta), sd = S$prior_sd), loglik)
}

# Score and negative curvature with respect to the log infection hazard.
# The limiting success score is one when a positive hazard underflows.
lgcp_cloglog_score <- function(h, y, n) {
  score <- curvature <- numeric(length(h))
  fail <- n > y
  score[fail] <- -(n-y)[fail] * h[fail]
  curvature[fail] <- -score[fail]
  pos <- y > 0 & h < 700
  t <- h[pos]
  prob <- -expm1(-t)
  b <- ratio <- numeric(length(t))
  small <- t < 1e-5
  b[small] <- 1-t[small]/2+t[small]^2/12
  ratio[small] <- t[small]/2+t[small]^2/12
  b[!small] <- t[!small]*exp(-t[!small])/prob[!small]
  ratio[!small] <- t[!small]/prob[!small]-1
  score[pos] <- score[pos] + y[pos]*b
  curvature[pos] <- curvature[pos] + y[pos]*b*ratio
  list(score=score, curvature=curvature)
}

lgcp_beta_newton <- function(beta, exposure, S) {
  local <- function(x) {
    eta <- S$survey_offset + as.numeric(S$X %*% x)
    survey <- lgcp_survey(exposure, eta, S$y, S$n, S$link)
    value <- survey$loglik - sum(x^2)/(2*S$prior_sd^2)
    if (!is.finite(value)) return(NULL)
    if (identical(S$link, "logit")) {
      # Evaluate the failure probability separately when p rounds to one.
      q <- stats::plogis(log(exposure) + eta, lower.tail = FALSE)
      p <- survey$probability
      deriv <- list(score = S$y*q - (S$n-S$y)*p, curvature = S$n*p*q)
    } else {
      reached <- exposure > 0
      h <- numeric(length(exposure))
      h[reached] <- exp(log(exposure[reached]) + eta[reached])
      deriv <- lgcp_cloglog_score(h,S$y,S$n)
    }
    score <- as.numeric(crossprod(S$X,deriv$score)) - x/S$prior_sd^2
    H <- crossprod(S$X,deriv$curvature*S$X) + diag(1/S$prior_sd^2,length(x))
    ch <- chol(H)
    center <- x + as.numeric(backsolve(ch, forwardsolve(t(ch), score)))
    list(ch=ch, center=center, value=value)
  }
  lgcp_regression_proposal(beta,local)
}
