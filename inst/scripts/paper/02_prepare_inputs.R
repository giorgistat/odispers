###############################################################################
# Survey locations with aridity, river-network quadrature nodes, village-node
# pairs and a prediction grid. Writes paper_results/inputs.rds.
###############################################################################

library(odispers)
library(sf)
source("00_settings.R")

data(remo)            # REMO nodule survey data
data(sudan_rivers)    # HydroRIVERS reaches with discharge >= 0.1 m3/s
data(sudan_aridity)   # Global Aridity Index v3.1 (packed terra raster, AI x 10,000)
data(sudan_borders)   # borders of Sudan and South Sudan

v <- subset(remo, country == "SUDAN")
v$y <- round(v$y)                       # a few counts in the file are not integers
v <- st_as_sf(v, coords = c("long", "lat"), crs = 4326, remove = FALSE)
v$AI <- extract_aridity(v, sudan_aridity)
stopifnot(all(is.finite(v$AI)), all(v$AI > 0))

bb <- st_bbox(st_buffer(st_as_sfc(st_bbox(v)), 1))           # 1 degree margin
s2 <- sf_use_s2(FALSE)                                         # planar test, as GDAL's filter
rivers <- suppressMessages(st_filter(sudan_rivers, st_as_sfc(bb), .predicate = st_intersects))
invisible(sf_use_s2(s2))
rivers <- rivers[rivers$DIS_AV_CMS >= min_flow, c("HYRIV_ID", "NEXT_DOWN", "DIS_AV_CMS")]

nodes <- river_nodes(rivers, spacing = spacing)
pairs <- river_pairs(v, nodes, max_range = max_range)
print(nodes); print(pairs)
v$dist_river_km <- nearest_node_distance(v, nodes)

# prediction grid over the two countries
sudans <- sudan_borders
crs_km <- attr(nodes, "crs")
grid <- st_make_grid(st_transform(sudans, crs_km), cellsize = grid_km, what = "centers")
grid <- grid[lengths(st_intersects(grid, st_transform(st_union(sudans), crs_km))) > 0]
grid <- st_sf(geometry = st_transform(grid, 4326))
grid$AI <- extract_aridity(grid, sudan_aridity, buffer_km = grid_km / 2)
grid <- grid[is.finite(grid$AI) & grid$AI > 0, ]
grid_pairs <- river_pairs(grid, nodes, max_range = max_range)

saveRDS(list(villages = v, rivers = st_drop_geometry(rivers), nodes = nodes, pairs = pairs,
             grid = grid, grid_pairs = grid_pairs),
        file.path(out_dir, "inputs.rds"))
