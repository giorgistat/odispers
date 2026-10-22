###############################################################################
# Zero-inflated river-network model:
#   pi(x)     = 1 - exp(-alpha E(x))                    (in a transmission zone)
#   lambda(x) = E(x)^gamma exp(beta0 + beta1 log AI + Z)  (force of infection)
#   P(Y = 0) = 1 - pi + pi (1 - p)^n,  P(Y = y > 0) = pi Bin(y; n, p)
# with breeding linear in log discharge, the field u on the river network
# (range 800 km, the estimate of 08_fit_u_linear.R) and the village effect Z.
# Fitted with gamma free, gamma = 0 (exposure only in pi) and gamma = 1, and
# compared with the model without zero inflation (08).
###############################################################################

library(odispers)
library(sf)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

inp <- load_result("inputs.rds")
v <- inp$villages; nodes <- inp$nodes
pairs <- river_pairs(v, nodes, phi_main)
graph <- river_graph(nodes, vertex_spacing = 5)
range_u <- 800

fit <- function(...) odisp_fit_zi(y ~ log(AI), ~ log(DIS_AV_CMS), v, nodes, graph, pairs,
                                  phi = phi_main, range_u = range_u, ...)
start <- c(log(2.5), log(0.8), log(1), 0)
fits <- list()
cat("gamma free\n");  fits$free   <- fit(start = start)
st <- fits$free$par
cat("gamma = 0\n");   fits$gamma0 <- fit(start = c(st[1:3], 0), fixed = list(gamma = 1e-6))
cat("gamma = 1\n");   fits$gamma1 <- fit(start = c(st[1:3], 0), fixed = list(gamma = 1 - 1e-6))
for (f in fits) print(f)

f08 <- file.path(out_dir, "u_fits_linear.rds")
ll08 <- if (file.exists(f08)) { b <- readRDS(f08)$best; b$logLik } else NA
tab <- data.frame(model = c("no zero inflation (08)", "zero-inflated, gamma free",
                            "zero-inflated, gamma = 0", "zero-inflated, gamma = 1"),
                  logLik = round(c(ll08, sapply(fits, `[[`, "logLik")), 2),
                  alpha = c(NA, signif(sapply(fits, function(f) f$hyper[["alpha"]]), 3)),
                  gamma = c(1, round(sapply(fits, function(f) f$hyper[["gamma"]]), 3)),
                  sigma_u = c(if (file.exists(f08)) round(b$hyper[["sigma_u"]], 3) else NA,
                              round(sapply(fits, function(f) f$hyper[["sigma_u"]]), 3)),
                  sigma_z = c(if (file.exists(f08)) round(b$hyper[["sigma_z"]], 3) else NA,
                              round(sapply(fits, function(f) f$hyper[["sigma_z"]]), 3)))
cat("\nComparison (Laplace log-likelihoods):\n"); print(tab, row.names = FALSE)

# where are villages predicted to lie outside a transmission zone?
best <- fits$free
cls <- aridity_classes(v$AI)
cat("\nMean probability of being in a transmission zone, by aridity class (gamma free):\n")
print(round(tapply(best$pi, cls[best$setup$has], mean), 3))
vd <- st_drop_geometry(v)[best$setup$has, ]; vd$pi <- best$pi
data(sudan_borders)
m_pi <- ggplot() +
  geom_sf(data = sudan_borders, fill = "grey97", colour = "grey40", linewidth = 0.3) +
  geom_point(data = vd, aes(long, lat, colour = pi, shape = y > 0), size = 1.4) +
  scale_colour_viridis_c(name = "P(transmission\nzone)", limits = c(0, 1)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "positives", `FALSE` = "no positives"),
                     name = NULL) +
  coord_sf(expand = FALSE) + theme_minimal(base_size = 10) + labs(x = NULL, y = NULL,
  title = "Probability that a village lies in a zone of local transmission",
  subtitle = sprintf("pi = 1 - exp(-alpha E), alpha = %.3g; gamma = %.2f", best$hyper[["alpha"]], best$hyper[["gamma"]]))
ggsave(file.path(fig_dir, "zi_pi_villages.png"), m_pi, width = 6.5, height = 7, dpi = 200, bg = "white")

saveRDS(list(fits = fits, comparison = tab), file.path(out_dir, "zi_fits.rds"))
