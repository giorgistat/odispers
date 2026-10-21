#' Pairs of locations and river nodes within a maximum range
#'
#' Finds, for every location, the network nodes within \code{max_range} km.
#' The likelihood is evaluated only over these pairs, so \code{max_range}
#' must exceed the largest dispersal range that the fitted model can reach;
#' \code{\link{odisp_fit}} warns if it does not.
#'
#' @param points An \code{sf} object of point locations (villages or a
#'   prediction grid).
#' @param nodes An \code{odisp_network} object from \code{\link{river_nodes}}.
#' @param max_range Maximum distance (km) at which a node can contribute.
#' @return An object of class \code{odisp_pairs}: a list with integer vectors
#'   \code{i} (location) and \code{j} (node), distances \code{r} (km), the
#'   number of locations \code{m}, and \code{max_range}.
#' @export
river_pairs <- function(points, nodes, max_range = 50) {
  if (!inherits(nodes, "odisp_network")) stop("'nodes' must come from river_nodes()")
  P <- sf::st_coordinates(sf::st_transform(sf::st_geometry(points), attr(nodes, "crs")))[, 1:2, drop = FALSE]
  nx <- nodes$x; ny <- nodes$y
  ord <- order(nx); sx <- nx[ord]
  ii <- jj <- rr <- vector("list", nrow(P))
  for (p in seq_len(nrow(P))) {
    lo <- findInterval(P[p, 1] - max_range, sx) + 1L
    hi <- findInterval(P[p, 1] + max_range, sx)
    if (hi < lo) next
    cand <- ord[lo:hi]
    cand <- cand[abs(ny[cand] - P[p, 2]) < max_range]
    if (!length(cand)) next
    d <- sqrt((nx[cand] - P[p, 1])^2 + (ny[cand] - P[p, 2])^2)
    keep <- d < max_range
    ii[[p]] <- rep.int(p, sum(keep)); jj[[p]] <- cand[keep]; rr[[p]] <- d[keep]
  }
  structure(list(i = as.integer(unlist(ii)), j = as.integer(unlist(jj)),
                 r = unlist(rr), m = nrow(P), max_range = max_range),
            class = "odisp_pairs")
}

#' @export
print.odisp_pairs <- function(x, ...) {
  n_with <- length(unique(x$i))
  cat(format(length(x$r), big.mark = ","), "location-node pairs within",
      x$max_range, "km;", n_with, "of", x$m, "locations have at least one node in range\n")
  invisible(x)
}

#' Distance from each location to the nearest network node
#'
#' @inheritParams river_pairs
#' @return Numeric vector of distances in km.
#' @export
nearest_node_distance <- function(points, nodes) {
  P <- sf::st_coordinates(sf::st_transform(sf::st_geometry(points), attr(nodes, "crs")))[, 1:2, drop = FALSE]
  vapply(seq_len(nrow(P)), function(p)
    sqrt(min((nodes$x - P[p, 1])^2 + (nodes$y - P[p, 2])^2)), numeric(1))
}
