###############################################################################
# Survey locations with aridity, river-network quadrature nodes, village-node
# pairs and a prediction grid. Writes paper_results/inputs.rds.
###############################################################################

library(odispers)
library(sf)
source("00_settings.R")

remo <- read.csv(remo_file)
v <- subset(remo, country == "SUDAN")
v$y <- round(v$y)                       # a few counts in the file are not integers
v <- st_as_sf(v, coords = c("long", "lat"), crs = 4326, remove = FALSE)
v$AI <- extract_aridity(v, ai_file)
stopifnot(all(is.finite(v$AI)), all(v$AI > 0))

bb <- st_bbox(st_buffer(st_as_sfc(st_bbox(v)), 1))           # 1 degree margin
rivers <- st_read(rivers_file, wkt_filter = st_as_text(st_as_sfc(bb)), quiet = TRUE)
rivers <- rivers[rivers$DIS_AV_CMS >= min_flow, c("HYRIV_ID", "NEXT_DOWN", "DIS_AV_CMS")]

nodes <- river_nodes(rivers, spacing = spacing)
pairs <- river_pairs(v, nodes, max_range = max_range)
print(nodes); print(pairs)
v$dist_river_km <- nearest_node_distance(v, nodes)

# prediction grid over the two countries
sudans <- sudan_boundaries()
crs_km <- attr(nodes, "crs")
grid <- st_make_grid(st_transform(sudans, crs_km), cellsize = grid_km, what = "centers")
grid <- grid[lengths(st_intersects(grid, st_transform(st_union(sudans), crs_km))) > 0]
grid <- st_sf(geometry = st_transform(grid, 4326))
grid$AI <- extract_aridity(grid, ai_file, buffer_km = grid_km / 2)
grid <- grid[is.finite(grid$AI) & grid$AI > 0, ]
grid_pairs <- river_pairs(grid, nodes, max_range = max_range)

saveRDS(list(villages = v, rivers = st_drop_geometry(rivers), nodes = nodes, pairs = pairs,
             grid = grid, grid_pairs = grid_pairs),
        file.path(out_dir, "inputs.rds"))
