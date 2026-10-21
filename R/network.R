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
  line_id <- xy[, "L1"]
  split_xy <- split(seq_len(nrow(xy)), line_id)

  out <- vector("list", length(split_xy))
  for (k in seq_along(split_xy)) {
    idx <- split_xy[[k]]
    px <- xy[idx, "X"]; py <- xy[idx, "Y"]
    seglen <- sqrt(diff(px)^2 + diff(py)^2)
    cum <- c(0, cumsum(seglen))
    L <- cum[length(cum)]
    if (L <= 0) next
    n <- max(1L, ceiling(L / spacing))
    target <- (seq_len(n) - 0.5) * L / n
    out[[k]] <- data.frame(line = as.integer(names(split_xy)[k]),
                           x = stats::approx(cum, px, xout = target, ties = "ordered")$y,
                           y = stats::approx(cum, py, xout = target, ties = "ordered")$y,
                           len = rep(L / n, n))
  }
  nodes <- do.call(rbind, out)
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
