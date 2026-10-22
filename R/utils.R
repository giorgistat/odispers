#' Wendland kernel with compact support on [0, 1]
#'
#' Evaluates \eqn{k(t) = (1 - t)_+^4 (1 + 4t)}, the Wendland function of
#' smoothness class \eqn{C^2}, which is zero for \eqn{t \ge 1}.
#'
#' @param t Non-negative numeric vector of scaled distances.
#' @return Numeric vector of kernel values.
#' @export
wendland <- function(t) {
  out <- numeric(length(t))
  inside <- t < 1
  ti <- t[inside]
  out[inside] <- (1 - ti)^4 * (1 + 4 * ti)
  out
}

#' Azimuthal equidistant projection centred on a set of features
#'
#' Returns a PROJ string for an azimuthal equidistant projection in
#' kilometres, centred on the bounding box of \code{x}. Distances from the
#' centre are exact and distortion is small over regions of a few thousand
#' kilometres, which makes it suitable for computing dispersal distances.
#'
#' @param x An \code{sf} or \code{sfc} object.
#' @return A character string giving the projection.
#' @export
local_aeqd <- function(x) {
  bb <- sf::st_bbox(sf::st_transform(x, 4326))
  lon0 <- round(unname((bb["xmin"] + bb["xmax"]) / 2), 3)
  lat0 <- round(unname((bb["ymin"] + bb["ymax"]) / 2), 3)
  sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +datum=WGS84 +units=km +no_defs", lat0, lon0)
}

# log(1 - exp(-x)) computed stably for x > 0
log1mexp <- function(x) {
  ifelse(x < log(2), log(-expm1(-x)), log1p(-exp(-x)))
}

#' Aridity index at point locations
#'
#' Extracts the Global Aridity Index (Zomer et al., 2022) at each point, as
#' the mean over a small buffer so that points on a river bank or near the
#' edge of the raster do not get missing values. Values stored as integers
#' multiplied by 10,000 are rescaled automatically.
#'
#' @param points An \code{sf} object of points.
#' @param raster A \code{terra::SpatRaster}, a packed raster such as \code{\link{sudan_aridity}}, or the path to the GeoTIFF.
#' @param buffer_km Radius (km) of the buffer over which the index is averaged.
#' @return Numeric vector of aridity index values.
#' @export
extract_aridity <- function(points, raster, buffer_km = 2) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("package 'terra' is required")
  if (is.character(raster)) raster <- terra::rast(raster)
  if (inherits(raster, "PackedSpatRaster")) raster <- terra::unwrap(raster)
  pts <- sf::st_transform(sf::st_geometry(points), 4326)
  crs_km <- local_aeqd(pts)
  buf <- sf::st_transform(sf::st_buffer(sf::st_transform(pts, crs_km), buffer_km), 4326)
  bb <- sf::st_bbox(buf)
  r <- terra::crop(raster, terra::ext(bb["xmin"] - 0.1, bb["xmax"] + 0.1, bb["ymin"] - 0.1, bb["ymax"] + 0.1))
  if (terra::global(r, "max", na.rm = TRUE)[1, 1] > 100) r <- r / 10000
  terra::extract(r, terra::vect(buf), fun = mean, na.rm = TRUE)[, 2]
}
