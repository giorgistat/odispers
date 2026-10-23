###############################################################################
# Figures for the stand-alone report on the zero-inflated model of 17_fit_zil.R
# (log E and log AI in both parts): prediction maps on the 20 km grid (river
# exposure, probability of a zone of transmission, prevalence within a zone,
# expected prevalence, and the area beyond the reach of any river), the mode
# of the breeding field u on the rivers, and the random-fold cross-validated
# non-randomised PIT (Czado, Gneiting and Held, 2009) among villages with at
# least one positive and reliability of P(Y >= 1), for this model alone.
###############################################################################

library(odispers)
library(sf)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

r <- load_result("zil_fit.rds"); fit <- r$fit; grid <- r$grid
net <- load_result("inputs_full_network.rds")
inp <- load_result("inputs.rds")
data(sudan_borders)
nm <- "zero-inflated, log E and log AI in both parts"
set.seed(2026)

base <- list(geom_sf(data = sudan_borders, fill = NA, colour = "grey25", linewidth = 0.3),
             coord_sf(expand = FALSE), theme_minimal(base_size = 10), labs(x = NULL, y = NULL))
g <- st_drop_geometry(grid)
gz <- subset(g, zero_exposure)                      # beyond 40 km of every river reach
cat(sprintf("Grid cells beyond the reach of any river (E = 0): %d of %d (%.1f%%)\n",
            nrow(gz), nrow(g), 100 * mean(g$zero_exposure)))
cat("Villages with E = 0:", sum(fit$E == 0), "\n")
cls <- aridity_classes(g$AI)
cat("Share of cells with E = 0 by aridity class:\n"); print(round(tapply(g$zero_exposure, cls, mean), 3))
grey0 <- geom_tile(data = gz, aes(lon, lat), fill = "grey80", width = 0.2, height = 0.2)
tile <- function(var) geom_tile(data = subset(g, !zero_exposure), aes(lon, lat, fill = .data[[var]]),
                                width = 0.2, height = 0.2)
note <- "Grey: more than 40 km from every river reach (E = 0, prevalence 0 by construction)"

m_E <- ggplot() + grey0 + tile("E") +
  scale_fill_viridis_c(option = "mako", direction = -1, trans = "log10", name = "Exposure E\n(log scale)",
                       labels = function(x) format(x, scientific = TRUE, digits = 1)) +
  base + labs(title = "River exposure E(x)", subtitle = note)
ggsave(file.path(fig_dir, "zil_exposure_map.png"), m_E, width = 6.5, height = 7, dpi = 200, bg = "white")

m_pi <- ggplot() + grey0 + tile("pi") + scale_fill_viridis_c(name = "P(transmission\nzone)", limits = c(0, 1)) +
  base + labs(title = "Probability of lying in a zone of local transmission", subtitle = note)
ggsave(file.path(fig_dir, "zil_pi_map.png"), m_pi, width = 6.5, height = 7, dpi = 200, bg = "white")

m_pz <- ggplot() + grey0 + tile("p_zone") +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Prevalence\nwithin a zone", limits = c(0, NA)) +
  base + labs(title = "Prevalence within a zone of transmission, 1 - exp(-lambda)", subtitle = note)
ggsave(file.path(fig_dir, "zil_pzone_map.png"), m_pz, width = 6.5, height = 7, dpi = 200, bg = "white")

m_prev <- ggplot() + grey0 + tile("prevalence") +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Expected\nprevalence", limits = c(0, NA)) +
  base + labs(title = "Expected nodule prevalence, pi x (1 - exp(-lambda))", subtitle = note)
ggsave(file.path(fig_dir, "zil_prevalence_map.png"), m_prev, width = 6.5, height = 7, dpi = 200, bg = "white")

# mode of u along the rivers (full network), interpolated from the graph vertices to the nodes
nodes <- net$nodes
graph <- river_graph(nodes, vertex_spacing = 5)
nodes$u <- as.vector(graph$A %*% fit$u)
ll <- st_coordinates(st_transform(st_as_sf(as.data.frame(nodes)[, c("x", "y")], coords = c("x", "y"),
                                           crs = attr(nodes, "crs")), 4326))
nd <- data.frame(lon = ll[, 1], lat = ll[, 2], u = nodes$u)
nd <- nd[order(abs(nd$u)), ]
vd <- st_drop_geometry(inp$villages)[fit$setup$has, ]
lim <- max(abs(quantile(nd$u, c(0.005, 0.995))))
m_u <- ggplot() + geom_sf(data = sudan_borders, fill = "grey97", colour = "grey40", linewidth = 0.3) +
  geom_point(data = nd, aes(lon, lat, colour = pmax(pmin(u, lim), -lim)), size = 0.12) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b", midpoint = 0, limits = c(-lim, lim),
                         name = "u (mode)\nlog breeding\nmultiplier") +
  geom_point(data = vd, aes(long, lat), shape = 1, size = 0.6, colour = "grey20") +
  coord_sf(xlim = sf::st_bbox(sudan_borders)[c(1, 3)], ylim = sf::st_bbox(sudan_borders)[c(2, 4)], expand = FALSE) +
  theme_minimal(base_size = 10) + labs(x = NULL, y = NULL,
  title = "Breeding field u on the river network", subtitle = "Whittle-Matern, range 800 km; circles: surveyed villages")
ggsave(file.path(fig_dir, "zil_u_map.png"), m_u, width = 6.5, height = 7, dpi = 200, bg = "white")

# ---- cross-validated calibration of this model only ---------------------------
pr <- load_result("cv_zil.rds")$preds[[nm]]
cond_cdf <- function(pmf) { p0 <- pmf[1]; c(0, cumsum(pmf[-1]) / (1 - p0)) }
pit_fun <- function(Fm1, F0, u) ifelse(u <= Fm1, 0, ifelse(u >= F0, 1, (u - Fm1) / pmax(F0 - Fm1, 1e-12)))
mean_pit_hist <- function(ys, cdfs, J = 10) {
  ub <- (0:J) / J
  Fbar <- rowMeans(sapply(seq_along(ys), function(i) pit_fun(cdfs[[i]][ys[i]], cdfs[[i]][ys[i] + 1], ub)))
  J * diff(Fbar)
}
pos <- Filter(function(z) z$y > 0, pr); cdfs <- lapply(pos, function(z) cond_cdf(z$pmf))
h <- mean_pit_hist(sapply(pos, `[[`, "y"), cdfs)
sims <- replicate(1000, { ys <- sapply(pos, function(z) { qq <- z$pmf[-1] / (1 - z$pmf[1]); sample.int(length(qq), 1, prob = qq) })
                          mean_pit_hist(ys, cdfs) })
pt <- data.frame(mid = (seq_along(h) - 0.5) / length(h), density = h,
                 lower = apply(sims, 1, stats::quantile, 0.025), upper = apply(sims, 1, stats::quantile, 0.975))
cat("\nNon-randomised PIT (villages with positives, n =", length(pos), "):\n"); print(round(pt, 3))
p_pit <- ggplot(pt, aes(mid, density)) + geom_col(width = 0.1, fill = "#7fbf7b", colour = "white") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.04, colour = "grey30") +
  geom_hline(yintercept = 1, linetype = 2) +
  labs(x = "Non-randomised probability integral transform", y = "Density",
       title = "PIT among villages with at least one positive",
       subtitle = "10 random folds; bars: 95% range under calibration (1000 simulations)") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "zil_only_pit.png"), p_pit, width = 5.8, height = 4.2, dpi = 200, bg = "white")

p1 <- sapply(pr, function(z) 1 - z$pmf[1]); o1 <- sapply(pr, `[[`, "y") > 0
br <- unique(stats::quantile(p1, seq(0, 1, 0.1))); b <- cut(p1, br, include.lowest = TRUE)
rel <- data.frame(pred = tapply(p1, b, mean), obs = tapply(o1, b, mean), n = as.vector(table(b)))
rel$lower <- stats::qbeta(0.025, rel$obs * rel$n + 0.5, (1 - rel$obs) * rel$n + 0.5)
rel$upper <- stats::qbeta(0.975, rel$obs * rel$n + 0.5, (1 - rel$obs) * rel$n + 0.5)
p_rel <- ggplot(rel, aes(pred, obs)) + geom_abline(linetype = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, colour = "#1b7837", alpha = 0.6) +
  geom_point(size = 2, colour = "#1b7837") + geom_line(colour = "#1b7837") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted probability of at least one positive", y = "Observed proportion",
       title = "Reliability of P(Y >= 1)", subtitle = "10 random folds; deciles of the prediction") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "zil_only_reliability.png"), p_rel, width = 4.8, height = 5, dpi = 200, bg = "white")
saveRDS(list(pit = pt, reliability = rel), file.path(out_dir, "zil_report.rds"))

# ---- exposure E along the river network -------------------------------------------
# E evaluated at points on the rivers (every 4th quadrature node, about 4-8 km
# apart, within Sudan and South Sudan), in chunks to limit memory
sel <- seq(1, nrow(nodes), by = 4)
pts <- st_as_sf(data.frame(x = nodes$x[sel], y = nodes$y[sel]), coords = c("x", "y"), crs = attr(nodes, "crs"))
inside <- lengths(st_intersects(st_transform(pts, st_crs(sudan_borders)), sudan_borders)) > 0
sel <- sel[inside]; pts <- pts[inside, ]
chunks <- split(seq_along(sel), ceiling(seq_along(sel) / 5000))
Er <- unlist(lapply(chunks, function(k) {
  pk <- pts[k, ]; pk$AI <- 1          # AI is not used by E; needed only for the design matrix
  predict(fit, pk, nodes, graph, pairs = river_pairs(pk, nodes, fit$setup$phi))$E
}))
er <- data.frame(lon = ll[sel, 1], lat = ll[sel, 2], E = Er)
cat(sprintf("\nExposure on the rivers (%d points): quantiles 1%%, 50%%, 99%%: %s\n", nrow(er),
            paste(signif(quantile(er$E, c(0.01, 0.5, 0.99)), 3), collapse = ", ")))
er <- er[order(er$E), ]
m_Er <- ggplot() + geom_sf(data = sudan_borders, fill = "grey97", colour = "grey40", linewidth = 0.3) +
  geom_point(data = er, aes(lon, lat, colour = E), size = 0.15) +
  scale_colour_viridis_c(option = "magma", direction = -1, trans = "log10", name = "Exposure E\n(log scale)",
                         labels = function(x) format(x, scientific = TRUE, digits = 1)) +
  coord_sf(xlim = sf::st_bbox(sudan_borders)[c(1, 3)], ylim = sf::st_bbox(sudan_borders)[c(2, 4)], expand = FALSE) +
  theme_minimal(base_size = 10) + labs(x = NULL, y = NULL, title = "River exposure E along the river network",
  subtitle = "Posterior mode, at points on the rivers 4-8 km apart")
ggsave(file.path(fig_dir, "zil_exposure_rivers.png"), m_Er, width = 6.5, height = 7, dpi = 200, bg = "white")
saveRDS(er, file.path(out_dir, "zil_exposure_rivers.rds"))
