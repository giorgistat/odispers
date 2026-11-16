# Conditional elliptical-slice updates on two colours of the river graph. Both
# colours move on every scan. A second partition shifts all band boundaries.
lgcp_block_setup <- function(S, graph, shift, width = 12L) {
  if (!is.null(graph$type) && !graph$type %in% c("wm", "taildown"))
    stop("block_slice requires a Whittle-Matern or tail-down river graph")
  labels <- lgcp_graph_blocks(graph$from, graph$to, graph$nv, width, shift)
  n_theta <- length(S$theta_names)
  field <- n_theta + seq_len(graph$nv)
  Q <- S$Q[field, field, drop = FALSE]
  A <- graph$A[S$used, , drop = FALSE]
  lapply(unique(labels$colour), function(colour) {
    idx <- which(labels$colour == colour)
    other <- which(labels$colour != colour)
    group <- match(labels$group[idx], unique(labels$group[idx]))
    self <- Matrix::drop0(A[, idx, drop = FALSE])
    nz <- Matrix::summary(self)
    row_group <- integer(nrow(self))
    row_group[nz$i] <- group[nz$j]
    if (any(row_group[nz$i] != group[nz$j]))
      stop("a cell couples two proposed independent blocks")
    rows <- which(row_group > 0L)
    list(idx = idx, other = other, ch = Matrix::Cholesky(Q[idx, idx, drop = FALSE], LDL = FALSE),
         cross = Q[idx, other, drop = FALSE], vertex_group = group,
         rows = rows, row_group = row_group[rows], groups = max(group),
         self = self[rows, , drop = FALSE], rest = A[rows, other, drop = FALSE])
  })
}

# With site counts fixed the breeding coefficients form a low-dimensional
# Poisson regression. Use its local curvature to propose a joint Newton step.
# Include the reverse Gaussian proposal density in the Metropolis ratio.
lgcp_breeding_step <- function(theta, Z, offset, counts, prior_sd) {
  cross_counts <- as.numeric(crossprod(Z, counts))
  local <- function(x) {
    mu <- exp(offset + as.numeric(Z %*% x))
    if (any(!is.finite(mu))) return(NULL)
    score <- cross_counts - as.numeric(crossprod(Z, mu)) - x / prior_sd^2
    H <- crossprod(Z, mu * Z) + diag(1 / prior_sd^2, length(x))
    ch <- chol(H)
    center <- x + as.numeric(backsolve(ch, forwardsolve(t(ch), score)))
    value <- sum(x * cross_counts) - sum(mu) - sum(x^2) / (2 * prior_sd^2)
    list(ch = ch, center = center, value = value)
  }
  lgcp_regression_proposal(theta, local)
}

# Local Gaussian proposals are generally asymmetric. Evaluate the proposal
# again at the candidate to correct both its centre and its normalisation.
lgcp_regression_proposal <- function(theta, local) {
  propose_logdens <- function(x, from) {
    sum(log(diag(from$ch))) - sum((from$ch %*% (x-from$center))^2) / 2
  }
  forward <- local(theta)
  if (is.null(forward)) stop("non-finite likelihood at current regression coefficients")
  candidate <- forward$center + as.numeric(backsolve(forward$ch, stats::rnorm(length(theta))))
  reverse <- local(candidate)
  if (is.null(reverse)) return(theta)
  ratio <- reverse$value - forward$value + propose_logdens(theta, reverse) -
    propose_logdens(candidate, forward)
  if (log(stats::runif(1)) < ratio) candidate else theta
}

lgcp_block_step <- function(state, counts, S, blocks) {
  theta_idx <- seq_along(S$theta_names)
  field_idx <- length(theta_idx) + seq_len(ncol(S$B) - length(theta_idx))
  x <- state[field_idx]
  offset <- S$log_mass
  if (length(theta_idx)) {
    Z <- S$breeding_matrix
    river <- as.numeric(S$field_projection %*% x)
    state[theta_idx] <- lgcp_breeding_step(state[theta_idx], Z, offset + river, counts, S$prior_sd)
    offset <- offset + as.numeric(Z %*% state[theta_idx])
  }
  for (b in blocks) {
    # x[idx] | x[other] has precision Q[idx,idx] and mean -Q[idx,idx]^{-1} Q[idx,other] x[other].
    center <- as.numeric(Matrix::solve(b$ch, -b$cross %*% x[b$other]))
    residual <- x[b$idx] - center
    noise <- lgcp_gaussian_draw(b$ch, length(b$idx))
    base <- offset[b$rows] + as.numeric(b$self %*% center + b$rest %*% x[b$other])
    angle <- lgcp_block_angles(b$row_group, b$groups, base, as.numeric(b$self %*% residual),
                               as.numeric(b$self %*% noise), counts[b$rows])
    x[b$idx] <- center + residual * angle[b$vertex_group, 1] + noise * angle[b$vertex_group, 2]
  }
  state[field_idx] <- x
  state
}
