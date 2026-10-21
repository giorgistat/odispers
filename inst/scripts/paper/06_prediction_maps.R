###############################################################################
# Maps of predicted prevalence, exposure and intensity of transmission.
###############################################################################

library(odispers); library(sf); library(ggplot2)
source("00_settings.R")
out <- readRDS(file.path(out_dir, "spline_fits.rds")); f <- out[[as.character(phi_main)]]$fit; cp <- out[[as.character(phi_main)]]$cp
inp <- readRDS(file.path(out_dir, "inputs.rds")); v <- inp$villages; nodes <- inp$nodes
grid <- inp$grid[, c("AI")]; grid$h <- pmax(log(grid$AI) - cp, 0)
pg <- predict(f, grid, nodes, pairs = inp$grid_pairs)
g <- cbind(grid, pg); xy <- st_coordinates(g); g$lon <- xy[, 1]; g$lat <- xy[, 2]
g$T <- ifelse(g$E > 0, g$lambda / g$E, exp(coef(f)["beta.(Intercept)"] + coef(f)["beta.log(AI)"] * log(g$AI) + coef(f)["beta.h"] * g$h))
cat("zero-exposure share:", round(mean(g$zero_exposure), 3), " max prevalence:", round(max(g$prevalence), 3), "\n")
sud <- sudan_boundaries()
base <- list(geom_sf(data = sud, fill = NA, colour = "grey25", linewidth = 0.3), coord_sf(expand = FALSE),
             theme_minimal(base_size = 10), labs(x = NULL, y = NULL))
tile <- function(var) geom_tile(aes(lon, lat, fill = .data[[var]]), width = 0.2, height = 0.2)
vd <- st_drop_geometry(v)
m1 <- ggplot(g) + tile("prevalence") +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Predicted\nprevalence", limits = c(0, NA)) +
  geom_point(data = vd, aes(long, lat), shape = 1, size = 0.5, colour = "grey30", alpha = 0.5) +
  base + labs(title = "Predicted nodule prevalence")
ggsave(file.path(fig_dir, "pred_prevalence.png"), m1, width = 6.5, height = 7, dpi = 200, bg = "white")
g$Erel <- ifelse(g$E > 0, g$E / median(g$E[g$E > 0]), NA)
m2 <- ggplot(g) + tile("Erel") +
  scale_fill_viridis_c(option = "mako", direction = -1, trans = "log10", na.value = "grey85",
                       name = "Exposure\n(relative to\nmedian)", oob = scales::squish, limits = c(0.01, 100)) +
  base + labs(title = "Exposure from the river network", subtitle = "Grey: no river within 40 km (zero exposure)")
ggsave(file.path(fig_dir, "pred_exposure.png"), m2, width = 6.5, height = 7, dpi = 200, bg = "white")
g$Trel <- g$T / exp(coef(f)["beta.(Intercept)"] + coef(f)["beta.log(AI)"] * log(0.5) + coef(f)["beta.h"] * max(log(0.5) - cp, 0))
m3 <- ggplot(g) + tile("Trel") +
  scale_fill_viridis_c(option = "viridis", trans = "log10", name = "Transmission\nintensity\n(relative to\nAI = 0.5)",
                       oob = scales::squish, limits = c(1e-3, 3)) +
  base + labs(title = "Intensity of transmission (aridity)")
ggsave(file.path(fig_dir, "pred_intensity.png"), m3, width = 6.5, height = 7, dpi = 200, bg = "white")
saveRDS(g, file.path(out_dir, "prediction_grid.rds"))
