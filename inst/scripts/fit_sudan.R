###############################################################################
# Fit the river-network dispersal model (no latent fields) to REMO nodule
# prevalence data for Sudan and South Sudan.
#
# Inputs (set the paths below):
#   - REMO data (Remo_AllCountries.csv), columns country, village, lat, long, n, y
#   - HydroRIVERS Africa shapefile (HydroRIVERS_v10_af.shp)
#   - Global Aridity Index v3.1, annual GeoTIFF (from Global-AI_ET0__annual_v3_1.zip)
#
# Outputs: model comparison table, estimates, profile likelihood for gamma,
# and maps of predicted prevalence, dispersal range and zero-exposure area,
# all written to out_dir.
###############################################################################

library(odispers)
library(sf)
library(terra)
library(ggplot2)

remo_file   <- "Remo_AllCountries.csv"
rivers_file <- "HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp"
ai_file     <- list.files(".", pattern = "^ai.*yr.*\\.tif$", recursive = TRUE,
                          full.names = TRUE, ignore.case = TRUE)[1]
out_dir     <- "odispers_results"
dir.create(out_dir, showWarnings = FALSE)

min_flow  <- 0.1   # m3/s; reaches with less flow are dropped
spacing   <- 2     # km between quadrature nodes
max_range <- 100   # km; largest dispersal range the model can reach (increase if the fit warns)
grid_km   <- 20    # resolution of the prediction grid

# ---- 1. Data -------------------------------------------------------------------
remo <- read.csv(remo_file)
v <- subset(remo, country == "SUDAN")
v$y <- round(v$y)                       # a few counts in the file are not integers
v <- st_as_sf(v, coords = c("long", "lat"), crs = 4326, remove = FALSE)

bb <- st_bbox(st_buffer(st_as_sfc(st_bbox(v)), 1))           # 1 degree margin
rivers <- st_read(rivers_file, wkt_filter = st_as_text(st_as_sfc(bb)), quiet = TRUE)
rivers <- rivers[rivers$DIS_AV_CMS >= min_flow, c("HYRIV_ID", "DIS_AV_CMS")]

v$AI <- extract_aridity(v, ai_file)
stopifnot(all(is.finite(v$AI)), all(v$AI > 0))

# ---- 2. Network quadrature and village-node pairs -------------------------------
nodes <- river_nodes(rivers, spacing = spacing)
print(nodes)
pairs <- river_pairs(v, nodes, max_range = max_range)
print(pairs)
v$dist_river_km <- nearest_node_distance(v, nodes)

# ---- 3. Fit a sequence of models ---------------------------------------------------
fit_one <- function(breeding, fixed = list(), background = FALSE)
  odisp_fit(y ~ 1, breeding, data = v, nodes = nodes, pairs = pairs,
            max_range = max_range, fixed = fixed, background = background)

fits <- list(
  "constant range"                     = fit_one(~ 0, fixed = list(gamma = 0)),
  "range depends on aridity"           = fit_one(~ 0),
  "+ log flow"                         = fit_one(~ log(DIS_AV_CMS)),
  "+ log flow + (log flow)^2"          = fit_one(~ log(DIS_AV_CMS) + I(log(DIS_AV_CMS)^2))
)
comparison <- data.frame(model = names(fits),
                         parameters = sapply(fits, function(f) attr(logLik(f), "df")),
                         logLik = round(sapply(fits, logLik), 2),
                         AIC = round(sapply(fits, AIC), 2))
print(comparison, row.names = FALSE)
write.csv(comparison, file.path(out_dir, "model_comparison.csv"), row.names = FALSE)

best <- fits[[which.min(comparison$AIC)]]
print(summary(best))
capture.output(summary(best), file = file.path(out_dir, "best_model_summary.txt"))

# ---- 4. How well is the aridity effect identified? --------------------------------
g_hat <- coef(best)["gamma"]
prof <- odisp_profile(best, "gamma", g_hat + seq(-0.6, 0.6, by = 0.15))
print(prof)
p_prof <- ggplot(prof, aes(value, logLik - max(logLik, na.rm = TRUE))) +
  geom_line() + geom_point() +
  geom_hline(yintercept = -qchisq(0.95, 1) / 2, linetype = 2) +
  labs(x = expression(gamma ~ "(effect of log aridity on log dispersal range)"),
       y = "Profile log-likelihood (relative to maximum)",
       title = "Profile likelihood for the aridity effect",
       subtitle = "Dashed line: approximate 95% confidence limit") +
  theme_minimal()
ggsave(file.path(out_dir, "profile_gamma.png"), p_prof, width = 6, height = 4, dpi = 200)

# Fitted dispersal range as a function of aridity
ai_seq <- exp(seq(log(min(v$AI)), log(max(v$AI)), length.out = 100))
p_phi <- ggplot(data.frame(AI = ai_seq,
                           phi = exp(coef(best)["eta0"] + coef(best)["gamma"] * log(ai_seq))),
                aes(AI, phi)) +
  geom_line() + scale_x_log10() +
  labs(x = "Aridity index (log scale; smaller = drier)", y = "Dispersal range (km)",
       title = "Fitted dispersal range against aridity") +
  theme_minimal()
ggsave(file.path(out_dir, "range_vs_aridity.png"), p_phi, width = 6, height = 4, dpi = 200)

# ---- 5. Observed against fitted ----------------------------------------------------
v$fitted <- fitted(best)
v$p_obs <- v$y / v$n
p_obs <- ggplot(v, aes(fitted, p_obs)) +
  geom_abline(linetype = 2) + geom_point(alpha = 0.4) +
  labs(x = "Fitted prevalence", y = "Observed prevalence", title = "Observed against fitted") +
  theme_minimal()
ggsave(file.path(out_dir, "observed_vs_fitted.png"), p_obs, width = 5, height = 5, dpi = 200)

# ---- 6. Prediction on a grid ---------------------------------------------------------
sudans <- st_as_sf(rnaturalearth::ne_countries(scale = "medium", returnclass = "sf"))
sudans <- sudans[sudans$admin %in% c("Sudan", "South Sudan"), ]
crs_km <- attr(nodes, "crs")
grid <- st_make_grid(st_transform(sudans, crs_km), cellsize = grid_km, what = "centers")
grid <- grid[lengths(st_intersects(grid, st_transform(st_union(sudans), crs_km))) > 0]
grid <- st_sf(geometry = st_transform(grid, 4326))
grid$AI <- extract_aridity(grid, ai_file, buffer_km = grid_km / 2)
grid <- grid[is.finite(grid$AI) & grid$AI > 0, ]
pred <- predict(best, grid, nodes)
grid <- cbind(grid, pred)
xy <- st_coordinates(grid)
grid$lon <- xy[, 1]; grid$lat <- xy[, 2]
saveRDS(list(fits = fits, best = best, grid = grid, villages = v), file.path(out_dir, "results.rds"))

map_base <- list(
  geom_sf(data = sudans, fill = NA, colour = "grey30", linewidth = 0.3),
  coord_sf(expand = FALSE), theme_minimal(base_size = 10), labs(x = NULL, y = NULL))

m_prev <- ggplot(grid) +
  geom_tile(aes(lon, lat, fill = prevalence), width = grid_km / 100, height = grid_km / 100) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Predicted\nprevalence") +
  map_base + labs(title = "Predicted nodule prevalence")
ggsave(file.path(out_dir, "map_predicted_prevalence.png"), m_prev, width = 7, height = 7.5, dpi = 200)

m_phi <- ggplot(grid) +
  geom_tile(aes(lon, lat, fill = phi), width = grid_km / 100, height = grid_km / 100) +
  scale_fill_viridis_c(name = "Dispersal\nrange (km)") +
  map_base + labs(title = "Fitted dispersal range")
ggsave(file.path(out_dir, "map_dispersal_range.png"), m_phi, width = 7, height = 7.5, dpi = 200)

m_zero <- ggplot(grid) +
  geom_tile(aes(lon, lat, fill = zero_exposure), width = grid_km / 100, height = grid_km / 100) +
  scale_fill_manual(values = c(`TRUE` = "#bdbdbd", `FALSE` = "#3182bd"),
                    labels = c(`TRUE` = "No river within reach (prevalence 0)",
                               `FALSE` = "Within reach of a river"), name = NULL) +
  geom_point(data = st_drop_geometry(v), aes(long, lat, colour = y > 0), size = 0.7) +
  scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "black"),
                      labels = c(`TRUE` = "Infected village", `FALSE` = "No infection found"), name = NULL) +
  map_base + labs(title = "Areas beyond the reach of any river") +
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "map_zero_exposure.png"), m_zero, width = 7, height = 8, dpi = 200)

cat("Results written to", out_dir, "\n")
