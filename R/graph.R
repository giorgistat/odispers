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
  estart <- c(0L, cumsum(m))[seq_len(nl)]
  q <- if (!is.null(nodes$DIS_AV_CMS)) nodes$DIS_AV_CMS[first][le] else rep(1, length(le))
  structure(list(from = as.integer(from), to = as.integer(to), length = elen, nv = nv,
                 A = A, vertex_spacing = vertex_spacing, type = "wm",
                 edge_reach = le, edge_q = q,
                 node_edge = as.integer(estart[ln] + s + 1L), node_fr = fr),
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
wm_precision <- function(graph, range, sigma = 1, mix = 0.5) {
  type <- if (is.null(graph$type)) "wm" else graph$type
  if (type != "wm") return(directed_precision(graph, range, mix) / sigma^2)
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

#' Directed (tail-up and tail-down) fields on a river graph
#'
#' Turns a river graph into the graph of a directed Gaussian field of Bolin,
#' Simas, Karlsson Strandh and Wallin (2026), the solution of the first-order
#' equation \eqn{(\kappa + \partial_\Gamma)(\tau u) = \dot W} with a forward
#' condition at the vertices and stationary anchoring at the sources. Along
#' each edge \eqn{u} is an Ornstein-Uhlenbeck process with variance
#' \eqn{\sigma^2} and correlation \eqn{e^{-\kappa d}}, as in the
#' Whittle-Matern field of \code{\link{wm_precision}}; the types differ at the
#' vertices.
#' \describe{
#'   \item{\code{"taildown"}}{Edges reversed, so that dependence propagates
#'     upstream from the outlets. The field is continuous, with stationary
#'     variance \eqn{\sigma^2} and covariance \eqn{\sigma^2 e^{-\kappa d}} at
#'     river distance \eqn{d}, also between points that are not
#'     flow-connected.}
#'   \item{\code{"tailup"}}{Edges oriented along the flow, with the
#'     flow-weighted condition K1: below a confluence the field starts from
#'     the average of the incoming values weighted by their discharge
#'     (\code{DIS_AV_CMS} of the incoming reaches). Points are correlated only
#'     if one is downstream of the other; tributaries are independent above
#'     their confluence.}
#'   \item{\code{"tailup+taildown"}}{Sum of independent tail-up and tail-down
#'     fields with variances \eqn{\phi\sigma^2} and
#'     \eqn{(1-\phi)\sigma^2}, with the same range.}
#' }
#' The precision is \eqn{Q = L^\top W L}, where \eqn{L} contains the
#' autoregressive recursions along the edges and \eqn{W} their conditional
#' precisions. For the tail-up field the latent vector holds the value at
#' each source and at the downstream end of each edge, so that incoming
#' values at a confluence are distinct; the matrix \code{A} interpolates
#' along each edge between its upstream and downstream values.
#'
#' @param graph An \code{odisp_graph} object from \code{\link{river_graph}}.
#' @param type One of \code{"taildown"}, \code{"tailup"} or
#'   \code{"tailup+taildown"}.
#' @return An \code{odisp_graph} object with the projection \code{A} and the
#'   number of latent values \code{nv} of the directed field, for use in
#'   \code{\link{odisp_fit_zil}}.
#' @export
directed_graph <- function(graph, type = c("taildown", "tailup", "tailup+taildown")) {
  type <- match.arg(type)
  if (is.null(graph$node_edge)) stop("rebuild the graph with the current river_graph()")
  ne <- length(graph$from); nn <- length(graph$node_edge)
  out <- graph; out$type <- type
  # tail-down: one value per vertex; each non-outlet vertex is the upstream
  # end of exactly one edge
  if (any(duplicated(graph$from))) stop("a vertex has more than one downstream edge")
  out$outlets <- setdiff(seq_len(graph$nv), graph$from)
  # tail-up: values at the sources and at the downstream end of each edge
  indeg <- tabulate(graph$to, graph$nv)
  src <- which(indeg == 0L)
  src_id <- integer(graph$nv); src_id[src] <- seq_along(src)
  ns <- length(src)
  q <- pmax(graph$edge_q, 1e-6)
  qin <- as.vector(tapply(c(q, rep(0, graph$nv)), c(graph$to, seq_len(graph$nv)), sum))
  p_in <- q / qin[graph$to]                            # share of each edge at its downstream vertex
  v <- graph$from
  is_src <- src_id[v] > 0L
  # incoming edges at the upstream vertex of each edge
  inc <- split(seq_len(ne), factor(graph$to, levels = seq_len(graph$nv)))
  ti <- c(seq_len(ne)[is_src], rep(seq_len(ne)[!is_src], lengths(inc[v[!is_src]])))
  tj <- c(src_id[v[is_src]], ns + unlist(inc[v[!is_src]], use.names = FALSE))
  tx <- c(rep(1, sum(is_src)), p_in[unlist(inc[v[!is_src]], use.names = FALSE)])
  out$tail_map <- Matrix::sparseMatrix(i = ti, j = tj, x = tx, dims = c(ne, ns + ne))
  out$n_src <- ns
  fr <- graph$node_fr; ee <- graph$node_edge
  A_up <- Matrix::Diagonal(nn, 1 - fr) %*% out$tail_map[ee, , drop = FALSE] +
    Matrix::sparseMatrix(i = seq_len(nn), j = ns + ee, x = fr, dims = c(nn, ns + ne))
  A_down <- graph$A
  out$nv_up <- ns + ne; out$nv_down <- graph$nv
  if (type == "taildown") { out$A <- A_down; out$nv <- graph$nv }
  if (type == "tailup") { out$A <- methods::as(A_up, "CsparseMatrix"); out$nv <- ns + ne }
  if (type == "tailup+taildown") { out$A <- cbind(A_up, A_down); out$nv <- ns + ne + graph$nv }
  out
}

# precision (sigma = 1) of the directed fields of directed_graph()
directed_precision <- function(graph, range, mix = 0.5) {
  kappa <- 2 / range
  rho <- exp(-kappa * graph$length)
  w <- 1 / (-expm1(-2 * kappa * graph$length))           # 1 / (1 - rho^2)
  ne <- length(graph$from)
  Q_down <- function() {
    nv <- graph$nv_down; no <- length(graph$outlets)
    L <- Matrix::sparseMatrix(i = c(graph$from, graph$from, graph$outlets),
                              j = c(graph$from, graph$to, graph$outlets),
                              x = c(rep(1, ne), -rho, rep(1, no)), dims = c(nv, nv))
    Wd <- numeric(nv); Wd[graph$from] <- w; Wd[graph$outlets] <- 1
    Matrix::forceSymmetric(Matrix::crossprod(L, Matrix::Diagonal(nv, Wd) %*% L))
  }
  Q_up <- function() {
    ns <- graph$n_src; n <- ns + ne
    L <- rbind(cbind(Matrix::Diagonal(ns), Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(ns, ne))),
               cbind(Matrix::sparseMatrix(i = integer(0), j = integer(0), dims = c(ne, ns)), Matrix::Diagonal(ne)) -
                 Matrix::Diagonal(ne, rho) %*% graph$tail_map)
    Matrix::forceSymmetric(Matrix::crossprod(L, Matrix::Diagonal(n, c(rep(1, ns), w)) %*% L))
  }
  Q <- switch(graph$type,
              taildown = Q_down(),
              tailup = Q_up(),
              "tailup+taildown" = Matrix::bdiag(Q_up() / mix, Q_down() / (1 - mix)))
  methods::as(Matrix::forceSymmetric(Q), "CsparseMatrix")
}
