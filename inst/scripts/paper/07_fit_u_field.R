###############################################################################
# Whittle-Matern breeding field u on the river network (alpha = 1), plus an
# independent village effect, added to the main model of 03_fit_spline.R
# (range 40 km, breeding quadratic in log discharge, aridity spline with the
# change point estimated there). The range of u is profiled over a grid.
###############################################################################

library(odispers)
library(sf)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

inp <- load_result("inputs.rds")
spl <- load_result("spline_fits.rds")[[as.character(phi_main)]]
v <- inp$villages; nodes <- inp$nodes
v$h <- pmax(log(v$AI) - spl$cp, 0)
pairs <- river_pairs(v, nodes, phi_main)

graph <- river_graph(nodes, vertex_spacing = 5)
print(graph)

form <- y ~ log(AI) + h
fit_z <- odisp_fit_u(form, breeding, v, nodes, graph, pairs, phi = phi_main, u = FALSE)
range_grid <- c(25, 50, 100, 200, 400, 800, 1600)   # the profile peaks at 400-800 km
fits_u <- odisp_fit_u(form, breeding, v, nodes, graph, pairs, phi = phi_main, range_u = range_grid)
prof <- attr(fits_u, "profile")
best <- fits_u[[which.max(prof$logLik)]]
# the model without u and without the village effect, in the same Laplace
# framework (village effect sd fixed at a negligible value), for comparison
S <- fit_z$setup
fit_none <- odispers:::odisp_u_laplace(list(field = "none", sig_z = 1e-4), S,
  list(beta = S$m_beta, theta = S$m_theta, u = NULL, s = NULL, z = rep(0, length(S$y))))

comparison <- data.frame(
  model = c("no u, no village effect", "village effect only",
            sprintf("u (range %g km) + village effect", range_grid)),
  logLik = round(c(fit_none$logLik, fit_z$logLik, prof$logLik), 2),
  sigma_u = round(c(0, 0, prof$sigma_field), 3),
  sigma_z = round(c(0, fit_z$hyper["sigma_z"], prof$sigma_z), 3))
cat("Laplace log-likelihoods (beta and theta integrated with diffuse priors);",
    "maximised log-likelihood of the model without u and Z:", round(spl$fit$logLik, 2), "\n")
print(comparison, row.names = FALSE)
print(best)

slopes <- function(f) {
  cf <- f$coefficients; V <- f$vcov
  b <- cf["beta.log(AI)"]; a <- b + cf["beta.h"]
  se_b <- sqrt(V["beta.log(AI)", "beta.log(AI)"])
  se_a <- sqrt(V["beta.log(AI)", "beta.log(AI)"] + V["beta.h", "beta.h"] + 2 * V["beta.log(AI)", "beta.h"])
  c(below = unname(b), se_below = unname(se_b), above = unname(a), se_above = unname(se_a))
}
cf0 <- coef(spl$fit); V0 <- vcov(spl$fit)
sl <- rbind("no u (ML)" = c(cf0["beta.log(AI)"], sqrt(V0["beta.log(AI)", "beta.log(AI)"]),
                             cf0["beta.log(AI)"] + cf0["beta.h"],
                             sqrt(V0["beta.log(AI)", "beta.log(AI)"] + V0["beta.h", "beta.h"] + 2 * V0["beta.log(AI)", "beta.h"])),
            "village effect only" = slopes(fit_z), "u + village effect" = slopes(best))
colnames(sl) <- c("slope below", "se", "slope above", "se")
cat("\nAridity slopes either side of the change point:\n"); print(round(sl, 2))

# mode of u along the rivers
un <- predict(best, v[1, ], nodes, graph)
nodes$u <- attr(un, "u_nodes")
ll <- st_coordinates(st_transform(st_as_sf(as.data.frame(nodes)[, c("x", "y")], coords = c("x", "y"),
                                           crs = attr(nodes, "crs")), 4326))
nodes$lon <- ll[, 1]; nodes$lat <- ll[, 2]
near <- unique(pairs$j)                         # nodes within reach of a village
data(sudan_borders)
sud <- sudan_borders
vd <- st_drop_geometry(v)
m_u <- ggplot() +
  geom_sf(data = sud, fill = "grey97", colour = "grey40", linewidth = 0.3) +
  geom_point(data = as.data.frame(nodes)[near, ], aes(lon, lat, colour = u), size = 0.25) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b", midpoint = 0,
                         name = "u (mode)\nlog breeding\nmultiplier") +
  geom_point(data = vd, aes(long, lat), shape = 1, size = 0.6, colour = "grey20") +
  coord_sf(expand = FALSE) + theme_minimal(base_size = 10) + labs(x = NULL, y = NULL,
  title = "Breeding field on the river network",
  subtitle = sprintf("Whittle-Matern, range %g km; rivers within 40 km of a village", best$hyper["range_u"]))
ggsave(file.path(fig_dir, "u_field_map.png"), m_u, width = 6.5, height = 7, dpi = 200, bg = "white")

# predictions (plug-in at the mode, village effect zero)
v$fit_u <- 1 - exp(-best$E * exp(as.vector(best$setup$X %*% best$coefficients[1:3])))[match(seq_len(nrow(v)), best$setup$has)]
v$fit_0 <- spl$fit$fitted
cat("correlation observed vs predicted: no u", round(cor(v$y / v$n, v$fit_0), 3),
    "; with u", round(cor(v$y / v$n, v$fit_u, use = "complete.obs"), 3), "\n")
cls <- aridity_classes(v$AI)
print(round(cbind(observed = tapply(v$y, cls, sum) / tapply(v$n, cls, sum),
                  no_u = tapply(v$fit_0 * v$n, cls, sum) / tapply(v$n, cls, sum),
                  with_u = tapply(v$fit_u * v$n, cls, sum) / tapply(v$n, cls, sum)), 3))
grid <- inp$grid; grid$h <- pmax(log(grid$AI) - spl$cp, 0)
pg <- predict(best, grid, nodes, graph, pairs = inp$grid_pairs)
grid$prevalence <- pg$prevalence
xy <- st_coordinates(grid); grid$lon <- xy[, 1]; grid$lat <- xy[, 2]
m_p <- ggplot(grid) +
  geom_tile(aes(lon, lat, fill = prevalence), width = 0.2, height = 0.2) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Predicted\nprevalence", limits = c(0, NA)) +
  geom_sf(data = sud, fill = NA, colour = "grey25", linewidth = 0.3) +
  coord_sf(expand = FALSE) + theme_minimal(base_size = 10) + labs(x = NULL, y = NULL,
  title = "Predicted nodule prevalence with the breeding field u")
ggsave(file.path(fig_dir, "pred_prevalence_u.png"), m_p, width = 6.5, height = 7, dpi = 200, bg = "white")

saveRDS(list(fit_z = fit_z, fits_u = fits_u, best = best, comparison = comparison, slopes = sl,
             villages = v, grid = grid), file.path(out_dir, "u_fits.rds"))
