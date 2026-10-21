#' Quadrature nodes along a river network
#'
#' Splits every river line into pieces of length at most \code{spacing} km
#' and returns one node at the midpoint of each piece, together with the
#' piece length. The nodes and lengths are the quadrature used to compute the
#' exposure integral \eqn{E(x) = \int_\Gamma K_\phi(\|x - s\|) w(s) ds}.
#'
#' @param rivers An \code{sf} object with LINESTRING or MULTILINESTRING
#'   geometries, for example a clip of HydroRIVERS. Attributes of each line
#'   are copied to the nodes it generates.
#' @param spacing Maximum length (km) of each piece; should be well below the
#'   smallest dispersal range of interest.
#' @param crs Projected coordinate reference system in kilometres. Defaults to
#'   an azimuthal equidistant projection centred on the network.
#' @return An object of class \code{odisp_network}: a data frame with node
#'   coordinates \code{x}, \code{y} (km), piece length \code{len} (km), the
#'   attributes of the parent line, and attribute \code{crs}.
#' @export
river_nodes <- function(rivers, spacing = 2, crs = NULL) {
  if (!inherits(rivers, "sf")) stop("'rivers' must be an sf object")
  if (is.null(crs)) crs <- local_aeqd(rivers)
  r <- sf::st_transform(rivers, crs)
  r <- suppressWarnings(sf::st_cast(sf::st_cast(r, "MULTILINESTRING"), "LINESTRING"))
  attrs <- sf::st_drop_geometry(r)
  xy <- sf::st_coordinates(r)
  L1 <- xy[, "L1"]
  # cumulative distance along each line, computed for all lines at once
  first <- c(TRUE, L1[-1] != L1[-length(L1)])
  seg <- c(0, sqrt(diff(xy[, "X"])^2 + diff(xy[, "Y"])^2))
  seg[first] <- 0
  cum <- stats::ave(seg, L1, FUN = cumsum)
  L <- as.vector(tapply(seg, L1, sum))
  ids <- as.integer(names(tapply(seg, L1, sum)))
  # lay the lines end to end on one axis, with a gap of 1 km between lines,
  # so that a single interpolation places the nodes on every line
  offset <- c(0, cumsum(L + 1))[seq_along(L)]
  pos <- cum + offset[match(L1, ids)]
  keep <- L > 0
  n <- ifelse(keep, pmax(1L, ceiling(L / spacing)), 0L)
  line <- rep(ids, n)
  k <- sequence(n)
  target <- (k - 0.5) * rep(L / pmax(n, 1), n) + rep(offset, n)
  nodes <- data.frame(line = line,
                      x = stats::approx(pos, xy[, "X"], xout = target, ties = "ordered")$y,
                      y = stats::approx(pos, xy[, "Y"], xout = target, ties = "ordered")$y,
                      len = rep(L / pmax(n, 1), n))
  nodes <- cbind(nodes, attrs[nodes$line, , drop = FALSE])
  rownames(nodes) <- NULL
  attr(nodes, "crs") <- crs
  class(nodes) <- c("odisp_network", "data.frame")
  nodes
}

#' @export
print.odisp_network <- function(x, ...) {
  cat("River network quadrature:", nrow(x), "nodes,",
      format(round(sum(x$len)), big.mark = ","), "km of river\n")
  cat("Projection:", attr(x, "crs"), "\n")
  invisible(x)
}

#' Aggregate quadrature nodes onto a coarser grid
#'
#' Merges nodes that fall in the same square cell of side \code{cell} km and,
#' optionally, in the same class of one or more numeric attributes (binned on
#' the log scale). Each merged node is placed at the length-weighted centroid
#' of the nodes it replaces, carries their total length, and the
#' length-weighted mean of the attributes in \code{vars}. This reduces the
#' number of nodes, and therefore of location-node pairs, at the cost of a
#' quadrature error of the order of \code{cell}, which is small when
#' \code{cell} is well below the dispersal range. Only the attributes in
#' \code{vars} are kept.
#'
#' @param nodes An \code{odisp_network} object.
#' @param cell Side (km) of the aggregation cells.
#' @param vars Names of numeric node attributes to carry over (for example
#'   \code{"DIS_AV_CMS"}); nodes are merged only within the same class.
#' @param log_width Width, on the natural log scale, of the classes used for
#'   \code{vars}; \code{0.5} merges nodes whose values differ by less than a
#'   factor of about 1.6.
#' @return An \code{odisp_network} object with fewer nodes.
#' @export
aggregate_nodes <- function(nodes, cell = 4, vars = NULL, log_width = 0.5) {
  if (!inherits(nodes, "odisp_network")) stop("'nodes' must come from river_nodes()")
  key <- paste(floor(nodes$x / cell), floor(nodes$y / cell))
  for (v in vars) {
    val <- nodes[[v]]
    if (!is.numeric(val) || any(val <= 0, na.rm = TRUE)) stop("'", v, "' must be numeric and positive")
    key <- paste(key, floor(log(val) / log_width))
  }
  g <- match(key, unique(key))
  L <- as.vector(rowsum(nodes$len, g))
  wmean <- function(z) as.vector(rowsum(z * nodes$len, g)) / L
  out <- data.frame(line = NA_integer_, x = wmean(nodes$x), y = wmean(nodes$y), len = L)
  for (v in vars) out[[v]] <- exp(wmean(log(nodes[[v]])))
  attr(out, "crs") <- attr(nodes, "crs")
  class(out) <- c("odisp_network", "data.frame")
  out
}
