###############################################################################
# Zero-inflated model with aridity in the probability of a transmission zone:
#   pi(x)     = 1 - exp(-alpha E(x) AI(x)^delta)
#   lambda(x) = E(x)^gamma exp(beta0 + beta1 log AI + Z)
# so that climate as well as the rivers determines whether a village lies in
# a zone of transmission, and u is left with the effects of the rivers. All of
# sigma_u, sigma_z, alpha, gamma and delta are estimated (range of u 800 km).
# The fit is then assessed by cross-validation as in 14_cv_zero_inflated.R.
###############################################################################

library(odispers)
library(sf)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

inp <- load_result("inputs.rds")
v <- inp$villages; nodes <- inp$nodes
pairs <- river_pairs(v, nodes, phi_main)
graph <- river_graph(nodes, vertex_spacing = 5)

# start from the zero-inflated fit without aridity in pi, if available
f13 <- file.path(out_dir, "zi_fits.rds")
start <- if (file.exists(f13)) c(readRDS(f13)$fits$free$par, delta = 1) else c(log(3), log(0.6), log(0.3), 0, 1)
fit <- odisp_fit_zi(y ~ log(AI), ~ log(DIS_AV_CMS), v, nodes, graph, pairs, phi = phi_main,
                    range_u = 800, start = start, fixed = list())
print(fit)
cls <- aridity_classes(v$AI)
cat("\nMean probability of being in a transmission zone, by aridity class:\n")
print(round(tapply(fit$pi, cls[fit$setup$has], mean), 3))

vd <- st_drop_geometry(v)[fit$setup$has, ]; vd$pi <- fit$pi
data(sudan_borders)
m_pi <- ggplot() +
  geom_sf(data = sudan_borders, fill = "grey97", colour = "grey40", linewidth = 0.3) +
  geom_point(data = vd, aes(long, lat, colour = pi, shape = y > 0), size = 1.4) +
  scale_colour_viridis_c(name = "P(transmission\nzone)", limits = c(0, 1)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "positives", `FALSE` = "no positives"),
                     name = NULL) +
  coord_sf(expand = FALSE) + theme_minimal(base_size = 10) + labs(x = NULL, y = NULL,
  title = "Probability that a village lies in a zone of local transmission",
  subtitle = sprintf("pi = 1 - exp(-alpha E AI^delta), alpha = %.3g, delta = %.2f; gamma = %.2f",
                     fit$hyper[["alpha"]], fit$hyper[["delta"]], fit$hyper[["gamma"]]))
ggsave(file.path(fig_dir, "zi_aridity_pi_villages.png"), m_pi, width = 6.5, height = 7, dpi = 200, bg = "white")
saveRDS(list(fits = list(free = fit)), file.path(out_dir, "zi_aridity_fits.rds"))

# cross-validated calibration, with the same folds
zi_file <- "zi_aridity_fits.rds"; zi_label <- "zero-inflated, aridity in pi"; zi_suffix <- "_zi_aridity"
source(system.file("scripts/paper/14_cv_zero_inflated.R", package = "odispers"))
