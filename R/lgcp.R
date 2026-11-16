#' Fit a latent log-Gaussian Cox breeding process to village prevalence
#'
#' An experimental sampler for the integrated river model with discrete
#' breeding sites. Conditional on the metric-graph Gaussian field, the site
#' count in river cell j is Poisson with mean length times intensity. The
#' log intensity is log(source_rate) plus the breeding predictor and the
#' projected field. Each site contributes a dispersal kernel centred at
#' its cell midpoint, multiplied by its strength. The same counts generate exposures at all villages.
#' The binomial survey predictor is \eqn{\log(E_i)+d_i^T\beta}. The default
#' complementary log-log link gives \eqn{p_i=1-\exp[-E_i\exp(d_i^T\beta)]}.
#' With a logit link the infection odds are \eqn{E_i\exp(d_i^T\beta)}.
#' Under either link zero exposure gives zero infection probability.
#'
#' A conditional Laplace approximation proposes the graph field and breeding
#' coefficients given counts. A Metropolis correction accounts for the error
#' in that approximation. Alternatively, elliptical slice sampling updates the
#' same block directly. Count updates use the joint survey likelihood and
#' prevalence coefficients use elliptical slice sampling, combined with Metropolis-
#' corrected local Gaussian proposals for multivariate regression in \code{block_slice}. All coefficients have
#' independent centred Gaussian priors. Inference is conditional on the stated
#' kernel cutoff and field hyperparameters. An optional Gaussian kernel width
#' is sampled jointly using a uniform prior on its logarithm within stated bounds.
#' The midpoint rule approximates
#' continuous source locations and the field is interpolated along graph edges.
#'
#' @param formula Prevalence formula with the number positive on the left.
#' @param breeding One-sided formula for log source density per km. Use
#'   \code{~ 0} to fix the baseline, or \code{~ 1} to estimate an intercept
#'   relative to \code{source_rate}. A free source intercept and a free
#'   transmission intercept can be weakly identified from prevalence alone.
#' @inheritParams odisp_design
#' @param data An \code{sf} object of survey locations containing the response,
#'   number examined and covariates used in \code{formula}.
#' @param scaling Kernel normalisation. For \code{"line"} divide the Wendland
#'   kernel by \eqn{2\phi/3}, for \code{"area"} by \eqn{\pi\phi^2/7}.
#'   A Gaussian kernel is normalized over its truncated line or disk support.
#' @param graph River graph built from these nodes in precisely the same order.
#' @param phi Fixed kernel radius in km.
#' @param range Practical range (km) passed to \code{\link{wm_precision}}.
#' @param sigma_u Positive graph field scale passed to \code{\link{wm_precision}}.
#' @param source_rate Baseline source intensity per km at zero field and zero
#'   breeding predictor. Mean site strength is one.
#' @param prior_sd Prior standard deviation for regression coefficients.
#' @param n_iter Total number of sweeps for one chain. Use \code{set.seed()}
#'   before fitting for reproducibility and run independent chains to assess mixing.
#' @param burnin Number of initial sweeps to discard.
#' @param block_size Maximum number of cells refreshed jointly in a count proposal.
#' @param keep_counts Save the retained joint count draws as an integer matrix.
#'   Columns correspond to the rows of \code{cells}.
#' @param field_update Use a Metropolis-corrected conditional Laplace update or
#'   elliptical slice sampling with one factorisation of the Gaussian prior.
#'   On a Whittle-Matern or tail-down graph, \code{"block_slice"} updates
#'   conditionally independent river blocks with separate slice angles.
#'   The slice methods report \code{NA} for field Metropolis acceptance.
#' @param marks With \code{"unit"} each site has strength one. With
#'   \code{"two_point"} each site independently has strength 1-d or 1+d
#'   with equal probability. The spread d has a Uniform(0,1) prior.
#' @param keep_prevalence Save the retained village prevalence draws for
#'   posterior summaries and spatial mixing diagnostics.
#' @param relocate Also propose moving sites between cells in overlapping
#'   spatial tiles, keeping site number and strength fixed during each move.
#' @param link Village binomial link, either \code{"cloglog"} or \code{"logit"}.
#'   The Poisson site intensity always uses a log link.
#' @param occurrence Optional one-sided formula for the log odds that a village
#'   is in a local transmission zone. Outside the zone all results are zero.
#'   Inside the zone the binomial survey uses the specified exposure and link.
#'   One village state is shared by everyone examined. With \code{NULL} all
#'   villages use the binomial survey model. Occurrence coefficients have the
#'   same centred Gaussian prior as the other regression coefficients.
#' @param occurrence_exposure Include the current marked kernel exposure as a
#'   linear covariate in the village activity logit. Requires an occurrence
#'   formula. The coefficient uses the same Gaussian prior as other coefficients.
#' @param kernel Dispersal kernel: the default Wendland or a Gaussian truncated
#'   at the fixed cutoff \code{phi}.
#' @param width Standard deviation (km) of the untruncated Gaussian kernel.
#'   Used as its initial value when \code{estimate_width = TRUE}.
#' @param estimate_width Sample the Gaussian width jointly with the latent
#'   sources and coefficients. Requires \code{kernel = "gaussian"}. Width
#'   proposals include both the survey and exposure-dependent activity likelihoods.
#' @param width_bounds Positive lower and upper bounds for a uniform prior on
#'   log width. Width proposals occur every four sweeps; their scale adapts only
#'   during burn-in. Posterior width draws are stored as \code{kernel_width}.
#' @return A list with retained draws of regression coefficients and total sites
#'   in cells within reach, cell posterior summaries, posterior mean exposure
#'   and prevalence in survey row order, optional joint count draws and sampler
#'   acceptance rates. For two-point marks the fit also reports cell mean
#'   marked mass and high-site counts, and saves matching high-site draws
#'   when \code{keep_counts} is true. An occurrence fit includes the posterior
#'   village activity probabilities in \code{fitted$pr_active}. Its prevalence
#'   draws include the sampled village state, shared across individuals.
#'   Cells outside kernel support integrate out conditional
#'   on the field. The count acceptance rate uses proposals that changed a count.
#'   Acceptance rates describe the updates and do not establish convergence.
#' @export
odisp_fit_lgcp <- function(formula, breeding = ~ 0, data, nodes, graph, pairs = NULL,
                           phi = 40, range = 100, sigma_u = 1, source_rate = 0.1,
                           trials = "n", scaling = c("line", "area"), prior_sd = 2,
                           n_iter = 2000L, burnin = floor(n_iter / 2),
                           block_size = 8L, keep_counts = FALSE,
                           field_update = c("laplace", "slice", "block_slice"), marks = c("unit", "two_point"),
                           keep_prevalence = FALSE, relocate = FALSE,
                           link = c("cloglog", "logit"), occurrence = NULL,
                           occurrence_exposure = FALSE,
                           kernel = c("wendland", "gaussian"), width = phi/3,
                           estimate_width = FALSE, width_bounds = c(phi/60, 4*phi)) {
  kernel <- match.arg(kernel)
  link <- match.arg(link)
  marks <- match.arg(marks)
  marked <- marks == "two_point"
  scaling <- match.arg(scaling)
  field_update <- match.arg(field_update)
  # Validate user controls before building the graph and kernel matrices.
  positive <- list(phi = phi, range = range, sigma_u = sigma_u,
                   source_rate = source_rate, prior_sd = prior_sd)
  for (name in names(positive)) {
    v <- positive[[name]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v <= 0)
      stop(name, " must be a finite positive scalar")
  }
  for (v in list(n_iter, burnin, block_size)) {
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v != floor(v) ||
        v < 0 || v > .Machine$integer.max)
      stop("n_iter, burnin and block_size must be non-negative integer scalars")
  }
  if (n_iter <= burnin || block_size < 1) stop("require n_iter > burnin and block_size >= 1")
  if (!is.logical(keep_counts) || length(keep_counts) != 1L || is.na(keep_counts))
    stop("keep_counts must be TRUE or FALSE")
  if (!is.logical(keep_prevalence) || length(keep_prevalence) != 1L || is.na(keep_prevalence))
    stop("keep_prevalence must be TRUE or FALSE")
  if (!is.logical(relocate) || length(relocate) != 1L || is.na(relocate))
    stop("relocate must be TRUE or FALSE")
  if (!is.logical(occurrence_exposure) || length(occurrence_exposure)!=1L || is.na(occurrence_exposure))
    stop("occurrence_exposure must be TRUE or FALSE")
  if (occurrence_exposure && is.null(occurrence)) stop("occurrence_exposure requires an occurrence formula")
  if (!is.logical(estimate_width) || length(estimate_width)!=1L || is.na(estimate_width))
    stop("estimate_width must be TRUE or FALSE")
  if (estimate_width && kernel != "gaussian") stop("estimate_width requires the Gaussian kernel")
  if (kernel == "gaussian" && (!is.numeric(width) || length(width)!=1L || !is.finite(width) || width<=0))
    stop("width must be a finite positive scalar")
  if (estimate_width && (!is.numeric(width_bounds) || length(width_bounds)!=2L ||
      any(!is.finite(width_bounds)) || width_bounds[1]<=0 || width_bounds[2]<=width_bounds[1] ||
      width<width_bounds[1] || width>width_bounds[2]))
    stop("width_bounds must be positive increasing bounds containing the initial width")
  S <- lgcp_setup(formula, breeding, data, nodes, graph, pairs, phi, range,
                  sigma_u, source_rate, trials, scaling, prior_sd, kernel, width)
  S$link <- link
  # The chain state consists of counts, optional strong-site counts, the
  # Gaussian field/breeding coefficients, and village regression coefficients.
  initial_width <- width
  width_proposal_sd <- .15
  width_every <- 4L
  width_accepts <- width_proposals <- width_attempts <- 0L
  high <- if (marked) integer(length(S$used)) else integer()
  spread <- if (marked) stats::runif(1L) else 0
  keep <- n_iter - burnin
  n_cell <- length(S$used)
  counts <- lgcp_initial_counts(S$K, S$y)
  if (marked) high <- stats::rbinom(n_cell, counts, 0.5)
  state <- numeric(ncol(S$B))
  beta <- numeric(ncol(S$X))
  has_occurrence <- !is.null(occurrence)
  gamma <- numeric()
  if (has_occurrence) {
    if (!inherits(occurrence, "formula") || length(occurrence) != 2L)
      stop("occurrence must be a one-sided formula")
    gate <- odisp_model_matrix(occurrence, sf::st_drop_geometry(data))
    if (occurrence_exposure) {
      if ("exposure" %in% colnames(gate$matrix)) stop("exposure is reserved for the latent exposure coefficient")
      gate$matrix <- cbind(gate$matrix, exposure=numeric(length(S$y)))
    }
    offset <- stats::model.offset(gate$frame)
    if (is.null(offset)) offset <- numeric(length(S$y))
    if (any(!is.finite(offset))) stop("occurrence offset must be finite")
    trials_original <- S$n
    G <- list(X=gate$matrix, survey_offset=offset, n=as.integer(trials_original > 0),
              link="logit", prior_sd=prior_sd)
    gamma <- numeric(ncol(G$X))
    mark_mass <- if (marked) (counts-high)*(1-spread) + high*(1+spread) else counts
    previous_exposure <- as.numeric(S$K %*% mark_mass)
    active_sum <- numeric(length(S$y))
  }
  n_theta <- length(S$theta_names)
  # Allocate retained draws and running posterior summaries once.
  coef_names <- c(if (length(beta)) paste0("beta.", colnames(S$X)),
                  if (n_theta) paste0("theta.", S$theta_names),
                  if (length(gamma)) paste0("occurrence.", colnames(G$X)))
  draw_names <- c(coef_names, "total_sites", if (marked) "mark_spread",
                  if (kernel == "gaussian") "kernel_width")
  draws <- matrix(0, keep, length(draw_names), dimnames = list(NULL, draw_names))
  count_draws <- if (keep_counts) matrix(0L, keep, n_cell, dimnames = list(NULL, S$used))
  high_draws <- if (keep_counts && marked) count_draws
  prevalence_draws <- if (keep_prevalence) matrix(0, keep, length(S$y))
  mean_high <- mean_mass <- numeric(n_cell)
  mean_count <- occupied <- intensity <- numeric(n_cell)
  mean_exposure <- prevalence <- numeric(length(S$y))
  field_accepts <- count_accepts <- count_proposals <- 0
  move_accepts <- move_proposals <- 0
  tiles <- list(integer())
  if (relocate) {
    grid <- expand.grid(width = phi * c(0.35, 1), dx = c(0, 0.5), dy = c(0, 0.5))
    tiles <- lapply(seq_len(nrow(grid)), function(k) {
      gx <- floor(nodes$x[S$used] / grid$width[k] + grid$dx[k])
      gy <- floor(nodes$y[S$used] / grid$width[k] + grid$dy[k])
      id <- paste(gx, gy)
      match(id, unique(id))
    })
  }
  laplace <- NULL
  if (field_update == "block_slice")
    partitions <- list(lgcp_block_setup(S, graph, 0L), lgcp_block_setup(S, graph, 6L))
  if (field_update == "slice") prior_ch <- Matrix::Cholesky(S$Q, LDL = FALSE)

  # Only the field update changes these Poisson means. Reuse them in the
  # next count sweep and in the retained intensity summary.
  mu <- exp(S$log_mass + as.numeric(S$B %*% state))
  for (iter in seq_len(n_iter)) {
    predictor <- S$survey_offset + as.numeric(S$X %*% beta)
    # Condition village activity on its survey, then update its regression.
    if (has_occurrence) {
      if (occurrence_exposure) G$X[,"exposure"] <- previous_exposure
      gate_eta <- G$survey_offset + as.numeric(G$X %*% gamma)
      prob <- lgcp_activity_probability(gate_eta, previous_exposure, predictor,
                                        S$y, trials_original, link)
      active <- stats::rbinom(length(prob), 1L, prob)
      # A masked village contributes neither a binomial outcome nor an observed activity label.
      G$y <- active * G$n
      if (length(gamma)) gamma <- if (iter %% 4L == 0L) {
        lgcp_beta_step(gamma, rep(1,length(prob)), G)
      } else {
        lgcp_beta_newton(gamma, rep(1,length(prob)), G)
      }
      masked <- trials_original == 0
      if (occurrence_exposure) {
        S$gate_slope <- unname(gamma[ncol(G$X)])
        # Calculate the fixed part directly to avoid subtracting a large latent term.
        fixed_cols <- seq_len(ncol(G$X)-1L)
        S$gate_base <- G$survey_offset + as.numeric(G$X[,fixed_cols,drop=FALSE] %*% gamma[fixed_cols])
        S$activity_sign <- as.integer((2L*active-1L)*G$n)
      } else if (any(masked)) {
        q <- stats::plogis(G$survey_offset[masked] + as.numeric(G$X[masked,,drop=FALSE] %*% gamma))
        active[masked] <- stats::rbinom(length(q), 1L, q)
      }
      S$n <- trials_original * active
    }
    # Counts, marks and kernel width share the same village exposure.
    sweep <- lgcp_count_sweep(counts, mu, predictor, S, block_size, high, spread,
                              tiles[[1L + (iter-1L) %% length(tiles)]])
    if (marked) {
      high <- sweep$high
      update <- lgcp_mark_step(sweep$counts, high, spread, predictor, S)
      spread <- update$spread
      sweep$exposure <- update$exposure
    }
    if (estimate_width && iter %% width_every == 0L) {
      width_update <- lgcp_width_step(width, sweep$counts, high, spread,
        sweep$exposure, predictor, S, width_bounds, width_proposal_sd)
      width <- width_update$width
      S$K <- width_update$K
      sweep$exposure <- width_update$exposure
      width_attempts <- width_attempts + 1L
      if (iter <= burnin) {
        # Freeze adaptation before retained sampling to preserve a fixed kernel.
        rate <- min(.05, 1/sqrt(width_attempts))
        width_proposal_sd <- exp(log(width_proposal_sd)+rate*(width_update$accepted-.35))
        width_proposal_sd <- min(max(width_proposal_sd,.005),1)
      } else {
        width_proposals <- width_proposals + 1L
        width_accepts <- width_accepts + width_update$accepted
      }
    }
    # Given counts, the field sees a Poisson likelihood; village outcomes
    # enter through those counts, rather than a second binomial factor.
    if (field_update == "laplace") {
      if (is.null(laplace) || !identical(sweep$counts, counts))
        laplace <- lgcp_field_laplace(sweep$counts, S)
      field <- lgcp_field_step(state, sweep$counts, S, laplace)
      state <- field$x
      if (iter > burnin) field_accepts <- field_accepts + field$accepted
    } else if (field_update == "block_slice") {
      state <- lgcp_block_step(state, sweep$counts, S, partitions[[1L + (iter - 1L) %% 2L]])
    } else {
      poisson_loglik <- function(x) lgcp_poisson_loglik(x, sweep$counts, S)
      state <- lgcp_elliptical_slice(state, lgcp_gaussian_draw(prior_ch, length(state)),
                                    poisson_loglik)
    }
    mu <- exp(S$log_mass + as.numeric(S$B %*% state))
    counts <- sweep$counts
    if (has_occurrence) previous_exposure <- sweep$exposure
    if (field_update == "block_slice" && length(beta) > 1L && iter %% 4L != 0L) {
      beta <- lgcp_beta_newton(beta, sweep$exposure, S)
    } else {
      beta <- lgcp_beta_step(beta, sweep$exposure, S)
    }
    if (iter <= burnin) next

    # Save one coherent joint state after all updates in this sweep.
    row <- iter - burnin
    draws[row, ] <- c(beta, state[seq_len(n_theta)], gamma, sum(counts), if (marked) spread,
                      if (kernel == "gaussian") width)
    if (marked) {
      if (keep_counts) high_draws[row, ] <- high
      mean_high <- mean_high + high
      mean_mass <- mean_mass + ((1 - spread) * (counts - high) + (1 + spread) * high)
    }
    if (keep_counts) count_draws[row, ] <- counts
    mean_count <- mean_count + counts
    occupied <- occupied + (counts > 0)
    intensity <- intensity + mu / S$length
    mean_exposure <- mean_exposure + sweep$exposure
    predictor <- S$survey_offset + as.numeric(S$X %*% beta)
    probability <- lgcp_survey(sweep$exposure, predictor, S$y, S$n, link)$probability
    if (has_occurrence) {
      if (occurrence_exposure && any(masked)) {
        q <- stats::plogis(S$gate_base[masked] + S$gate_slope*sweep$exposure[masked])
        active[masked] <- stats::rbinom(length(q),1L,q)
      }
      probability <- probability * active
      active_sum <- active_sum + active
    }
    prevalence <- prevalence + probability
    if (keep_prevalence) prevalence_draws[row, ] <- probability
    count_accepts <- count_accepts + sweep$accepts
    count_proposals <- count_proposals + sweep$proposals
    move_accepts <- move_accepts + sweep$moves
    move_proposals <- move_proposals + sweep$move_proposals
  }

  cells <- data.frame(node = S$used, mean_count = mean_count / keep,
                      pr_occupied = occupied / keep, mean_rate_per_km = intensity / keep)
  if (marked) {
    cells$mean_high_count <- mean_high / keep
    cells$mean_marked_mass <- mean_mass / keep
  }
  fit <- list(call = match.call(), marks = marks, high_draws = high_draws, field_update = field_update,
       settings = list(phi = phi, range = range, sigma_u = sigma_u, source_rate = source_rate,
                       prior_sd = prior_sd, scaling = scaling, n_iter = n_iter, burnin = burnin,
                       block_size = block_size, relocate = relocate, link = link,
                       occurrence = occurrence, occurrence_exposure=occurrence_exposure,
                       kernel=kernel, width=if(kernel=="gaussian" && !estimate_width) width else NULL,
                       estimate_width=estimate_width,
                       width_bounds=if(estimate_width) width_bounds else NULL), draws = draws,
       count_draws = count_draws, prevalence_draws = prevalence_draws, cells = cells,
       fitted = data.frame(exposure = mean_exposure / keep, prevalence = prevalence / keep),
       acceptance = c(field = if (field_update == "laplace") field_accepts / keep else NA_real_,
                      count = if (count_proposals) count_accepts / count_proposals else NA_real_))
  if (relocate) fit$acceptance["relocate"] <- if (move_proposals) move_accepts / move_proposals else NA_real_
  if (kernel == "gaussian") fit$initial_width <- initial_width
  if (estimate_width) {
    fit$acceptance["width"] <- if (width_proposals) width_accepts/width_proposals else NA_real_
    fit$width_proposal_sd <- width_proposal_sd
  }
  if (has_occurrence) fit$fitted$pr_active <- active_sum / keep
  structure(fit, class = "odisp_lgcpfit")
}

#' @export
print.odisp_lgcpfit <- function(x, ...) {
  cat("Latent river LGCP with", x$field_update, "field updates\n")
  cat("Village binomial link:", if (is.null(x$settings$link)) "cloglog" else x$settings$link, "\n")
  cat(nrow(x$draws), "retained draws;", nrow(x$cells), "cells within reach\n")
  cat("Posterior mean sites within reach:", format(mean(x$draws[, "total_sites"])), "\n")
  invisible(x)
}

#' @export
fitted.odisp_lgcpfit <- function(object, ...) object$fitted$prevalence

