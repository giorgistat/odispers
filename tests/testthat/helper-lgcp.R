# One 2-km reach, surveyed twice within range and once outside.
lgcp_test_data <- function() {
  crs <- "+proj=aeqd +lat_0=0 +lon_0=30 +datum=WGS84 +units=km +no_defs"
  line <- sf::st_linestring(rbind(c(0, 0), c(2, 0)))
  rivers <- sf::st_sf(HYRIV_ID = 1L, NEXT_DOWN = 0L,
                      geometry = sf::st_sfc(line, crs = crs))
  nodes <- river_nodes(rivers, spacing = 3, crs = crs)
  v <- sf::st_as_sf(data.frame(x = c(1, 1.5, 10), north = 0, n = c(4, 3, 4),
                               y = c(2, 1, 0)), coords = c("x", "north"), crs = crs)
  list(nodes = nodes, graph = river_graph(nodes), data = v)
}

