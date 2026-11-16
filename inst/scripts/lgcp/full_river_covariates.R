# Prepare intensity and transmission predictors for both Sudan LGCP models.
# Input is the full 40-km survey/river design. If it lacks the two HydroATLAS
# attributes, recreate them by matching the reach ID in the shipped dataset.
library(odispers)
args <- commandArgs(trailingOnly = TRUE)
coordinates <- "--coordinates" %in% args
args <- args[args != "--coordinates"]
rad_arg <- grep("^--radius=",args,value=TRUE)
if (length(rad_arg)>1L) stop("supply a single --radius in km")
args <- args[!args %in% rad_arg]
if (length(args) != 2L) stop("supply full inputs.rds and output directory, optionally --coordinates")
a <- readRDS(args[1])
out <- args[2]
radius <- if (length(rad_arg)) as.numeric(sub("^--radius=","",rad_arg)) else if (!is.null(a$phi)) a$phi else 40
if (length(radius)!=1L || !is.finite(radius) || radius<=0) stop("invalid dispersal radius")
if (is.null(a$pairs) || is.null(a$pairs$max_range) || a$pairs$max_range < radius)
  a$pairs <- river_pairs(a$villages, a$nodes, radius)
a$phi <- radius
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!all(c("log_slope", "inundation") %in% names(a$nodes))) {
  data(sudan_river_covariates, package = "odispers")
  hydro <- sudan_river_covariates[match(a$nodes$HYRIV_ID, sudan_river_covariates$HYRIV_ID), ]
  stopifnot(!anyNA(hydro$HYRIV_ID))
  a$nodes$log_slope <- log1p(hydro$sgr_dk_rav)
  a$nodes$inundation <- pmin(log1p(hydro$inu_pc_ult), 1.5)
}
if (!"AI" %in% names(a$villages)) {
  data(sudan_aridity, package = "odispers")
  a$villages$AI <- extract_aridity(a$villages, sudan_aridity)
}
standardise <- function(x, w, name, support) {
  stopifnot(all(is.finite(x)), all(is.finite(w)), all(w > 0))
  center <- weighted.mean(x, w)
  scale <- sqrt(weighted.mean((x-center)^2, w))
  stopifnot(scale>0)
  list(z=(x-center)/scale, row=data.frame(predictor=name, support=support, center=center, scale=scale))
}
Z <- list(z_flow=log(a$nodes$DIS_AV_CMS), z_slope=a$nodes$log_slope, z_flood=a$nodes$inundation)
if (coordinates) {
  xy <- data.frame(x=a$nodes$x, y=a$nodes$y)
  points <- sf::st_as_sf(xy, coords=c("x","y"), crs=attr(a$nodes,"crs"))
  lonlat <- sf::st_coordinates(sf::st_transform(points,4326))
  Z$z_lon <- lonlat[,1]
  Z$z_lat <- lonlat[,2]
}
scales <- list()
for (name in names(Z)) {
  z <- standardise(Z[[name]], a$nodes$len, name, "all river cells, weighted by length")
  a$nodes[[name]] <- z$z
  scales[[name]] <- z$row
}
z <- standardise(log(a$villages$AI), rep(1,nrow(a$villages)), "z_logAI", "survey villages")
a$villages$z_logAI <- z$z
scales$z_logAI <- z$row
a$source_covariates <- names(Z)
a$survey_covariates <- "z_logAI"
write.csv(do.call(rbind,scales), file.path(out,"covariate_scales.csv"), row.names=FALSE)
saveRDS(a, file.path(out,"inputs.rds"), compress=FALSE)
