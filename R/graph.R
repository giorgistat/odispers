#' Metric graph of a river network for a Whittle-Matern field
#'
#' Builds the graph on which the breeding field \eqn{u} is defined. Every
#' reach is an edge between its upstream and downstream vertices; the
#' downstream vertex of a reach is the upstream vertex of the reach it flows
#' into (\code{NEXT_DOWN} in HydroRIVERS), so confluences join the reaches.
#' Reaches longer than \code{vertex_spacing} are split into equal edges by
#' extra vertices. The quadrature nodes from \code{\link{river_nodes}} lie
#' along the edges, and \eqn{u} at a node is the linear interpolation of
#' \eqn{u} at the two ends of its edge (matrix \code{A}).
#'
#' @param nodes An \code{odisp_network} object from \code{\link{river_nodes}},
#'   not aggregated, with a column \code{HYRIV_ID}.
#' @param next_down Downstream reach identifier for each node, or \code{NULL}
#'   to use the column \code{NEXT_DOWN} of \code{nodes}. Identifiers not
#'   present in \code{nodes} (including 0) mark outlets.
#' @param vertex_spacing Maximum edge length (km).
#' @return An object of class \code{odisp_graph}: a list with the edge list
#'   \code{from}, \code{to}, \code{length} (km), the number of vertices
#'   \code{nv}, and the sparse projection matrix \code{A} (nodes x vertices).
#' @export
river_graph <- function(nodes, next_down = NULL, vertex_spacing = 5) {
  if (!inherits(nodes, "odisp_network")) stop("'nodes' must come from river_nodes()")
  if (anyNA(nodes$line)) stop("aggregated nodes cannot be used to build the graph")
  if (is.null(next_down)) next_down <- nodes$NEXT_DOWN
  if (is.null(next_down)) stop("supply 'next_down' or a NEXT_DOWN column in 'nodes'")
  rl <- rle(nodes$line)
  if (anyDuplicated(rl$values)) stop("nodes of each line must be contiguous")
  nl <- length(rl$values)
  first <- cumsum(c(1L, rl$lengths))[seq_len(nl)]
  id <- nodes$HYRIV_ID[first]
  if (anyDuplicated(id)) stop("each reach must be a single line")
  L <- nodes$len[first] * rl$lengths                     # reach length
  m <- pmax(1L, as.integer(ceiling(L / vertex_spacing))) # edges per reach
  nxt <- match(next_down[first], id)                     # NA = outlet

  # vertex numbering: upstream vertex of each reach, then interior vertices,
  # then one vertex for each outlet
  up <- seq_len(nl)
  n_int <- m - 1L
  int_start <- nl + c(0L, cumsum(n_int))[seq_len(nl)]
  out_id <- rep(NA_integer_, nl)
  is_out <- is.na(nxt)
  out_id[is_out] <- nl + sum(n_int) + seq_len(sum(is_out))
  down <- ifelse(is_out, out_id, up[pmax(nxt, 1L)])
  nv <- nl + sum(n_int) + sum(is_out)

  # vertex sequence along reach l: up, interior..., down
  vert_at <- function(l, s) {  # s = 0..m[l], vectorised over l and s
    ifelse(s == 0L, up[l], ifelse(s == m[l], down[l], int_start[l] + s))
  }
  le <- rep(seq_len(nl), m)
  se <- sequence(m) - 1L
  from <- vert_at(le, se); to <- vert_at(le, se + 1L)
  elen <- (L / m)[le]

  # projection of the nodes onto the vertices
  ln <- rep(seq_len(nl), rl$lengths)
  k <- sequence(rl$lengths)
  t <- (k - 0.5) * nodes$len / (L[ln] / m[ln])           # position in edge units
  s <- pmin(floor(t), m[ln] - 1L)
  fr <- t - s
  A <- Matrix::sparseMatrix(i = rep(seq_len(nrow(nodes)), 2),
                            j = c(vert_at(ln, s), vert_at(ln, s + 1L)),
                            x = c(1 - fr, fr), dims = c(nrow(nodes), nv))
  structure(list(from = as.integer(from), to = as.integer(to), length = elen, nv = nv,
                 A = A, vertex_spacing = vertex_spacing),
            class = "odisp_graph")
}

#' @export
print.odisp_graph <- function(x, ...) {
  cat("River graph:", format(x$nv, big.mark = ","), "vertices,",
      format(length(x$from), big.mark = ","), "edges,",
      format(round(sum(x$length)), big.mark = ","), "km\n")
  invisible(x)
}

#' Precision matrix of a Whittle-Matern field on a river graph
#'
#' Exact precision at the vertices of the Gaussian Whittle-Matern field with
#' \eqn{\alpha = 1} (exponential covariance along each edge, Kirchhoff
#' conditions at the vertices) of Bolin, Simas and Wallin (2024). Each edge
#' of length \eqn{\ell} contributes
#' \deqn{\frac{1}{\sigma^2(1-\rho^2)}\begin{pmatrix}(1+\rho^2)/2 & -\rho\\
#' -\rho & (1+\rho^2)/2\end{pmatrix}, \qquad \rho = e^{-\kappa\ell},}
#' so on a long chain the field has variance \eqn{\sigma^2} and correlation
#' \eqn{e^{-\kappa d}} at distance \eqn{d} along the river. The variance is
#' larger at river ends and smaller at confluences, as implied by the
#' Kirchhoff conditions.
#'
#' @param graph An \code{odisp_graph} object.
#' @param range Practical range (km), \eqn{2/\kappa}: correlation \eqn{e^{-2}
#'   \approx 0.14} at this distance along the network.
#' @param sigma Standard deviation \eqn{\sigma}.
#' @return A sparse symmetric matrix.
#' @export
wm_precision <- function(graph, range, sigma = 1) {
  kappa <- 2 / range
  rho <- exp(-kappa * graph$length)
  c1 <- 1 / (-expm1(-2 * kappa * graph$length))           # 1 / (1 - rho^2)
  d <- c1 * (1 + rho^2) / 2
  o <- -rho * c1
  Q <- Matrix::sparseMatrix(i = c(graph$from, graph$to, pmin(graph$from, graph$to)),
                            j = c(graph$from, graph$to, pmax(graph$from, graph$to)),
                            x = c(d, d, o), dims = c(graph$nv, graph$nv),
                            symmetric = TRUE)
  Q / sigma^2
}
