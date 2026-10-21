# A small synthetic study area: three rivers and a set of villages, in
# geographic coordinates near the equator (about 111 km per degree).
toy_setup <- function(m = 300, seed = 1) {
  set.seed(seed)
  rivers <- sf::st_sf(
    flow = c(50, 5, 0.5),
    geometry = sf::st_sfc(
      sf::st_linestring(rbind(c(30.0, 0.0), c(30.3, 0.6), c(30.2, 1.2), c(30.5, 2.0))),
      sf::st_linestring(rbind(c(30.3, 0.6), c(30.9, 0.8), c(31.4, 0.5))),
      sf::st_linestring(rbind(c(30.2, 1.2), c(29.6, 1.5), c(29.2, 1.4))),
      crs = 4326))
  lon <- runif(m, 29.0, 31.6); lat <- runif(m, -0.1, 2.1)
  villages <- sf::st_as_sf(data.frame(lon = lon, lat = lat, n = 40L,
                                      AI = exp(-2 + 1.6 * (lat - 0) / 2.1)),
                           coords = c("lon", "lat"), crs = 4326)
  nodes <- river_nodes(rivers, spacing = 1)
  pairs <- river_pairs(villages, nodes, max_range = 80)
  list(rivers = rivers, villages = villages, nodes = nodes, pairs = pairs)
}
