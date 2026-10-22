###############################################################################
# Figures for the zero-inflated model of 13_fit_zero_inflated.R (gamma free,
# no aridity in pi): maps of the probability of a transmission zone and of
# the expected prevalence on the prediction grid, and the random-fold
# cross-validated calibration of 14_cv_zero_inflated.R against the model
# without zero inflation.
###############################################################################

library(odispers)
library(sf)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

inp <- load_result("inputs.rds")
zi <- load_result("zi_fits.rds")$fits$free
cv <- readRDS(file.path(out_dir, "cv_zero_inflated.rds"))
nodes <- inp$nodes; v <- inp$villages
graph <- river_graph(nodes, vertex_spacing = 5)
data(sudan_borders)

# ---- maps on the prediction grid ------------------------------------------------
grid <- inp$grid
pg <- predict(zi, grid, nodes, graph, pairs = inp$grid_pairs)
grid <- cbind(grid, pg)
xy <- st_coordinates(grid); grid$lon <- xy[, 1]; grid$lat <- xy[, 2]
base <- list(geom_sf(data = sudan_borders, fill = NA, colour = "grey25", linewidth = 0.3),
             coord_sf(expand = FALSE), theme_minimal(base_size = 10), labs(x = NULL, y = NULL))
m_pi <- ggplot(grid) + geom_tile(aes(lon, lat, fill = pi), width = 0.2, height = 0.2) +
  scale_fill_viridis_c(name = "P(transmission\nzone)", limits = c(0, 1)) + base +
  labs(title = "Probability of lying in a zone of local transmission")
ggsave(file.path(fig_dir, "zi_pi_map.png"), m_pi, width = 6.5, height = 7, dpi = 200, bg = "white")
m_prev <- ggplot(grid) + geom_tile(aes(lon, lat, fill = prevalence), width = 0.2, height = 0.2) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Expected\nprevalence", limits = c(0, NA)) +
  base + labs(title = "Expected nodule prevalence, zero-inflated model")
ggsave(file.path(fig_dir, "zi_prevalence_map.png"), m_prev, width = 6.5, height = 7, dpi = 200, bg = "white")

# ---- random-fold calibration ----------------------------------------------------
pit <- subset(cv$pit, scheme == "random folds")
pit$model <- factor(ifelse(grepl("zero-inflated", pit$analysis), "zero-inflated", "no zero inflation"),
                    c("no zero inflation", "zero-inflated"))
p_pit <- ggplot(pit, aes(mid, density)) +
  geom_col(width = 0.1, fill = "#7fbf7b", colour = "white") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.04, colour = "grey30") +
  geom_hline(yintercept = 1, linetype = 2) + facet_wrap(~ model) +
  labs(x = "Probability integral transform", y = "Density") + theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "zi_cv_pit_random.png"), p_pit, width = 8, height = 3.6, dpi = 200, bg = "white")
rel <- subset(cv$reliability, scheme == "random folds")
rel$model <- ifelse(grepl("zero-inflated", rel$analysis), "zero-inflated", "no zero inflation")
p_rel <- ggplot(rel, aes(pred, obs, colour = model)) +
  geom_abline(linetype = 2) + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, alpha = 0.6) +
  geom_point(size = 2) + geom_line() + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = c("zero-inflated" = "#1b7837", "no zero inflation" = "grey50"), name = NULL) +
  labs(x = "Predicted probability of at least one positive", y = "Observed proportion") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "zi_cv_reliability_random.png"), p_rel, width = 5.5, height = 5.5, dpi = 200, bg = "white")

cat("Random-fold scores:\n"); print(subset(cv$scores, grepl("^random", analysis)), row.names = FALSE)
cat("\nZero villages by aridity class (random folds):\n"); print(subset(cv$zeros, grepl("^random", analysis)), row.names = FALSE)
cat("\nPredicted prevalence on the grid: share with pi < 0.05:", round(mean(grid$pi < 0.05), 3),
    "; maximum expected prevalence:", round(max(grid$prevalence), 3), "\n")
