###############################################################################
# Simplified model with the breeding field u: aridity linear in log AI in the
# intensity of transmission, and breeding intensity linear in log discharge,
#   w(s) = exp{theta log Q(s) + u(s)},  T(x) = exp{beta0 + beta1 log AI(x) + Z(x)},
# with the dispersal range fixed at 40 km. The range of u is profiled. The
# results are compared with the spline/quadratic model of 07_fit_u_field.R.
###############################################################################

library(odispers)
library(sf)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

inp <- readRDS(file.path(out_dir, "inputs.rds"))
v <- inp$villages; nodes <- inp$nodes
pairs <- river_pairs(v, nodes, phi_main)
graph <- river_graph(nodes, vertex_spacing = 5)

form_lin <- y ~ log(AI)
breeding_lin <- ~ log(DIS_AV_CMS)

# without latent effects (maximum likelihood), for reference
fit_ml <- suppressWarnings(odisp_fit(form_lin, breeding_lin, data = v, nodes = nodes, pairs = pairs,
                                     max_range = phi_main, fixed = list(gamma = 0, eta0 = log(phi_main))))
print(summary(fit_ml)$table[, 1:2])

fit_z <- odisp_fit_u(form_lin, breeding_lin, v, nodes, graph, pairs, phi = phi_main, u = FALSE)
range_grid <- c(100, 200, 400, 800, 1600)
fits_u <- odisp_fit_u(form_lin, breeding_lin, v, nodes, graph, pairs, phi = phi_main, range_u = range_grid)
prof <- attr(fits_u, "profile")
best <- fits_u[[which.max(prof$logLik)]]

comparison <- data.frame(
  model = c("village effect only", sprintf("u (range %g km) + village effect", range_grid)),
  logLik = round(c(fit_z$logLik, prof$logLik), 2),
  sigma_u = round(c(0, prof$sigma_u), 3),
  sigma_z = round(c(fit_z$hyper["sigma_z"], prof$sigma_z), 3))
cat("\nLinear model, Laplace log-likelihoods (beta and theta integrated with diffuse priors):\n")
print(comparison, row.names = FALSE)
cat("Maximised log-likelihood without u and Z:", round(fit_ml$logLik, 2), "\n\n")
print(fit_z)
print(best)

# comparison with the spline/quadratic model at the same range of u
uq <- file.path(out_dir, "u_fits.rds")
if (file.exists(uq)) {
  fq <- readRDS(uq)$fits_u
  r <- as.character(best$hyper["range_u"])
  if (r %in% names(fq))
    cat(sprintf("\nSame range of u (%s km): Laplace log-likelihood linear %.2f, spline/quadratic %.2f\n",
                r, best$logLik, fq[[r]]$logLik))
}

# mode of u along the rivers
un <- predict(best, v[1, ], nodes, graph)
nodes$u <- attr(un, "u_nodes")
ll <- st_coordinates(st_transform(st_as_sf(as.data.frame(nodes)[, c("x", "y")], coords = c("x", "y"),
                                           crs = attr(nodes, "crs")), 4326))
nodes$lon <- ll[, 1]; nodes$lat <- ll[, 2]
near <- unique(pairs$j)
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
  title = "Breeding field on the river network (linear model)",
  subtitle = sprintf("Whittle-Matern, range %g km; rivers within 40 km of a village", best$hyper["range_u"]))
ggsave(file.path(fig_dir, "u_field_map_linear.png"), m_u, width = 6.5, height = 7, dpi = 200, bg = "white")

# predictions (plug-in at the mode, village effect zero)
p <- length(best$setup$m_beta)
lam <- best$E * exp(as.vector(best$setup$X %*% best$coefficients[seq_len(p)]))
v$fit_u <- (1 - exp(-lam))[match(seq_len(nrow(v)), best$setup$has)]
v$fit_ml <- fit_ml$fitted
v$p_obs <- v$y / v$n
cat("correlation observed vs predicted: without u", round(cor(v$p_obs, v$fit_ml), 3),
    "; with u", round(cor(v$p_obs, v$fit_u, use = "complete.obs"), 3), "\n")
cls <- aridity_classes(v$AI)
print(round(cbind(observed = tapply(v$y, cls, sum) / tapply(v$n, cls, sum),
                  without_u = tapply(v$fit_ml * v$n, cls, sum) / tapply(v$n, cls, sum),
                  with_u = tapply(v$fit_u * v$n, cls, sum) / tapply(v$n, cls, sum)), 3))

p_ov <- ggplot(st_drop_geometry(v), aes(fit_u, p_obs, colour = cls)) +
  geom_abline(linetype = 2) + geom_point(alpha = 0.55, size = 1.6) +
  scale_colour_manual(values = c("#a6611a", "#dfc27d", "#7fbf7b", "#1b7837"), name = "Aridity class") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted prevalence", y = "Observed prevalence",
       title = "Observed against predicted prevalence (linear model with u)",
       subtitle = "Plug-in at the posterior mode, village effect set to zero; in-sample") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "obs_vs_pred_u_linear.png"), p_ov, width = 8, height = 6, dpi = 180, bg = "white")

grid <- inp$grid
pg <- predict(best, grid, nodes, graph, pairs = inp$grid_pairs)
grid$prevalence <- pg$prevalence
xy <- st_coordinates(grid); grid$lon <- xy[, 1]; grid$lat <- xy[, 2]
m_p <- ggplot(grid) +
  geom_tile(aes(lon, lat, fill = prevalence), width = 0.2, height = 0.2) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Predicted\nprevalence", limits = c(0, NA)) +
  geom_sf(data = sud, fill = NA, colour = "grey25", linewidth = 0.3) +
  coord_sf(expand = FALSE) + theme_minimal(base_size = 10) + labs(x = NULL, y = NULL,
  title = "Predicted nodule prevalence (linear model with u)")
ggsave(file.path(fig_dir, "pred_prevalence_u_linear.png"), m_p, width = 6.5, height = 7, dpi = 200, bg = "white")

saveRDS(list(fit_ml = fit_ml, fit_z = fit_z, fits_u = fits_u, best = best, comparison = comparison,
             villages = v, grid = grid), file.path(out_dir, "u_fits_linear.rds"))
