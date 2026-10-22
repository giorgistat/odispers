###############################################################################
# Spatial correlation in the exposure or in the intensity of transmission?
# Same fixed effects as 08_fit_u_linear.R (aridity linear in log AI, breeding
# linear in log discharge, range 40 km), with the spatial random effect either
#   E: a Whittle-Matern field u on the river network, in the breeding intensity
#   T: a Gaussian process S on the plane, in the intensity of transmission
# each with and without the independent village effect Z. The ranges of u and
# S are profiled. All models use the same Laplace approximation and priors, so
# their approximate log-likelihoods are comparable.
###############################################################################

library(odispers)
library(sf)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

inp <- load_result("inputs.rds")
v <- inp$villages; nodes <- inp$nodes
pairs <- river_pairs(v, nodes, phi_main)
graph <- river_graph(nodes, vertex_spacing = 5)

form_lin <- y ~ log(AI)
breeding_lin <- ~ log(DIS_AV_CMS)
fit <- function(...) odisp_fit_u(form_lin, breeding_lin, v, nodes, graph, pairs, phi = phi_main, ...)

range_s <- c(10, 25, 50, 100, 200, 400, 800, 1600)   # km, straight-line distance
range_u <- c(200, 400, 800, 1600)           # km, along the river network

fits <- list()
fits$T1 <- fit(field = "S", range = range_s, z = FALSE)      # S only
fits$T2 <- fit(field = "S", range = range_s, z = TRUE)       # S + Z
fits$E1 <- fit(field = "u", range = range_u, z = FALSE)      # u only
uq <- file.path(out_dir, "u_fits_linear.rds")                # u + Z, from 08
fits$E2 <- if (file.exists(uq)) readRDS(uq)$fits_u else fit(field = "u", range = c(100, range_u), z = TRUE)

labels <- c(T1 = "S on the plane, in T", T2 = "S on the plane, in T, + village effect Z",
            E1 = "u on the river network, in E", E2 = "u on the river network, in E, + village effect Z")
profiles <- do.call(rbind, lapply(names(fits), function(m) {
  pr <- attr(fits[[m]], "profile")
  if (is.null(pr$range)) pr$range <- pr$range_u
  if (is.null(pr$sigma_field)) pr$sigma_field <- pr$sigma_u
  data.frame(model = m, description = labels[[m]], range = pr$range, logLik = round(pr$logLik, 2),
             sigma_field = round(pr$sigma_field, 3), sigma_z = round(pr$sigma_z, 3))
}))
cat("\nApproximate log-likelihoods by model and range of the field:\n")
print(profiles, row.names = FALSE)

best <- lapply(fits, function(f) f[[which.max(sapply(f, `[[`, "logLik"))]])
# hyperparameter of a fit, also for fits made by earlier versions of odisp_fit_u()
hyp <- function(b, name, old) if (name %in% names(b$hyper)) b$hyper[[name]] else b$hyper[[old]]
v$p_obs <- v$y / v$n
summary_tab <- do.call(rbind, lapply(names(best), function(m) {
  b <- best[[m]]; cf <- b$coefficients; se <- sqrt(diag(b$vcov))
  p <- length(b$setup$m_beta)
  lam <- b$E * exp(as.vector(b$setup$X %*% cf[seq_len(p)]) + if (!is.null(b$s)) b$s else 0)
  pred <- rep(NA, nrow(v)); pred[b$setup$has] <- 1 - exp(-lam)
  data.frame(model = m, description = labels[[m]], range = hyp(b, "range", "range_u"),
             logLik = round(b$logLik, 2),
             sigma_field = round(hyp(b, "sigma_field", "sigma_u"), 3),
             sigma_z = round(b$hyper[["sigma_z"]], 3),
             beta_logAI = sprintf("%.2f (%.2f)", cf["beta.log(AI)"], se["beta.log(AI)"]),
             theta_logQ = sprintf("%.2f (%.2f)", cf["theta.log(DIS_AV_CMS)"], se["theta.log(DIS_AV_CMS)"]),
             cor_obs_pred = round(cor(v$p_obs, pred, use = "complete.obs"), 3))
}))
cat("\nBest range for each model:\n")
print(summary_tab, row.names = FALSE)

# Gaussian process S (mode at the villages, kriged to the grid) and predicted
# prevalence for the best model with S
bs <- if (best$T2$logLik > best$T1$logLik) best$T2 else best$T1
grid <- inp$grid
pg <- predict(bs, grid, nodes, graph, pairs = inp$grid_pairs)
grid$S <- pg$S; grid$prevalence <- pg$prevalence
xy <- st_coordinates(grid); grid$lon <- xy[, 1]; grid$lat <- xy[, 2]
data(sudan_borders)
sud <- sudan_borders
base <- list(geom_sf(data = sud, fill = NA, colour = "grey25", linewidth = 0.3), coord_sf(expand = FALSE),
             theme_minimal(base_size = 10), labs(x = NULL, y = NULL))
m_s <- ggplot(grid) +
  geom_tile(aes(lon, lat, fill = S), width = 0.2, height = 0.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "grey92", high = "#b2182b", midpoint = 0,
                       name = "S (mode)\nlog multiplier\nof intensity") +
  geom_point(data = st_drop_geometry(v), aes(long, lat), shape = 1, size = 0.5, colour = "grey30") +
  base + labs(title = "Gaussian process S in the intensity of transmission",
              subtitle = sprintf("Range %g km%s", bs$hyper[["range"]],
                                 if (bs$hyper[["sigma_z"]] > 0) ", with village effect" else ""))
ggsave(file.path(fig_dir, "s_field_map.png"), m_s, width = 6.5, height = 7, dpi = 200, bg = "white")
m_p <- ggplot(grid) +
  geom_tile(aes(lon, lat, fill = prevalence), width = 0.2, height = 0.2) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Predicted\nprevalence", limits = c(0, NA)) +
  base + labs(title = "Predicted nodule prevalence with S in the intensity")
ggsave(file.path(fig_dir, "pred_prevalence_s.png"), m_p, width = 6.5, height = 7, dpi = 200, bg = "white")

# profiles of the four models against the range of the field
pp <- ggplot(profiles, aes(range, logLik, colour = description)) + geom_line() + geom_point() +
  scale_x_log10() +
  scale_colour_manual(values = c("#b2182b", "#ef8a62", "#2166ac", "#67a9cf"), name = NULL) +
  labs(x = "Range of the field (km, log scale; along the network for u, straight-line for S)",
       y = "Laplace log-likelihood", title = "Spatial random effect in the exposure or in the intensity") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(file.path(fig_dir, "profile_s_vs_u.png"), pp, width = 7, height = 5.5, dpi = 200, bg = "white")

saveRDS(list(fits = fits, best = best, profiles = profiles, summary = summary_tab, grid = grid),
        file.path(out_dir, "s_vs_u_fits.rds"))
