###############################################################################
# Maps for the comparison of directed breeding fields (23_directed_profiles.R
# and 24_tailup_joint_warm.R): posterior mode of the breeding field u along
# the rivers and expected prevalence pi x p on the 20 km grid, for the joint
# fits of the tail-down, tail-up and tail-up + tail-down models.
###############################################################################
library(odispers)
library(sf)
library(Matrix)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
inp <- load_result("inputs.rds"); net <- load_result("inputs_full_network.rds")
nodes <- net$nodes
graph <- river_graph(nodes, vertex_spacing = 5)
data(sudan_borders)
types <- c("taildown", "tailup", "tailup+taildown")
labels <- c(taildown = "Tail-down", tailup = "Tail-up (K1, discharge weights)",
            "tailup+taildown" = "Tail-up + tail-down")
fits <- lapply(types, function(t) readRDS(file.path(out_dir, sprintf("directed_%s_joint.rds", gsub("\\+", "_", t)))))
names(fits) <- types
ll <- st_coordinates(st_transform(st_as_sf(as.data.frame(nodes)[, c("x", "y")], coords = c("x", "y"),
                                           crs = attr(nodes, "crs")), 4326))
sel <- seq(1, nrow(nodes), by = 3)
grid <- inp$grid; xy <- st_coordinates(grid)
um <- NULL; pm <- NULL
for (t in types) {
  g <- directed_graph(graph, t); b <- fits[[t]]
  um <- rbind(um, data.frame(model = labels[[t]], lon = ll[sel, 1], lat = ll[sel, 2],
                             u = as.vector(g$A[sel, , drop = FALSE] %*% b$u)))
  pg <- predict(b, grid, nodes, g, pairs = net$grid_pairs)
  pm <- rbind(pm, data.frame(model = labels[[t]], lon = xy[, 1], lat = xy[, 2], prevalence = pg$prevalence,
                             pi = pg$pi, zero = pg$zero_exposure))
  cat(t, ": expected prevalence on the grid, median", round(median(pg$prevalence), 4),
      " 99th percentile", round(quantile(pg$prevalence, 0.99), 3), "\n")
}
um$model <- factor(um$model, labels); pm$model <- factor(pm$model, labels)
um <- um[order(abs(um$u)), ]
lim <- max(abs(quantile(um$u, c(0.005, 0.995))))
bb <- sf::st_bbox(sudan_borders)
m_u <- ggplot() + geom_sf(data = sudan_borders, fill = "grey97", colour = "grey40", linewidth = 0.3) +
  geom_point(data = um, aes(lon, lat, colour = pmax(pmin(u, lim), -lim)), size = 0.08) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b", midpoint = 0, limits = c(-lim, lim),
                         name = "u (mode)") + facet_wrap(~ model, nrow = 1) +
  coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], expand = FALSE) + theme_minimal(base_size = 9) +
  labs(x = NULL, y = NULL)
ggsave(file.path(fig_dir, "directed_joint_u_maps.png"), m_u, width = 13, height = 5, dpi = 200, bg = "white")
m_p <- ggplot() + geom_tile(data = subset(pm, zero), aes(lon, lat), fill = "grey80", width = 0.2, height = 0.2) +
  geom_tile(data = subset(pm, !zero), aes(lon, lat, fill = prevalence), width = 0.2, height = 0.2) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Expected\nprevalence", limits = c(0, NA)) +
  geom_sf(data = sudan_borders, fill = NA, colour = "grey25", linewidth = 0.3) + facet_wrap(~ model, nrow = 1) +
  coord_sf(expand = FALSE) + theme_minimal(base_size = 9) + labs(x = NULL, y = NULL)
ggsave(file.path(fig_dir, "directed_joint_prevalence_maps.png"), m_p, width = 13, height = 5, dpi = 200, bg = "white")
# differences between the tail-down and the other fields on the grid
w <- split(pm$prevalence, pm$model)
cat("correlation of expected prevalence on the grid: tail-down vs tail-up", round(cor(w[[1]], w[[2]]), 3),
    "; tail-down vs combined", round(cor(w[[1]], w[[3]]), 4), "\n")
cat("max |difference| tail-down vs combined:", round(max(abs(w[[1]] - w[[3]])), 4), "\n")
saveRDS(list(u = um, prevalence = pm), file.path(out_dir, "directed_maps.rds"))
