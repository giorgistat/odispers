###############################################################################
# Exploratory analysis of REMO nodule prevalence in Sudan / South Sudan
# in relation to the river network (HydroRIVERS) and aridity (Global-AI v3.1)
#
# Inputs come with the package: data(remo), data(sudan_rivers),
# data(sudan_aridity) and data(sudan_borders).
#
# All figures are written to the folder "figures/", and a village table with
# every derived variable to "sudan_villages_explore.csv".
#
# install.packages(c("sf", "terra", "ggplot2", "dplyr", "tidyr",
#                    "rnaturalearth", "rnaturalearthdata", "scales"))
###############################################################################

library(sf)
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)


# ---- 0. Settings ---------------------------------------------------------------
data(remo, sudan_rivers, sudan_aridity, sudan_borders, package = "odispers")
out_dir     <- "figures"
dir.create(out_dir, showWarnings = FALSE)

xlim <- c(21, 40); ylim <- c(2.5, 23.5)                 # study area (lon, lat)
aeqd <- "+proj=aeqd +lat_0=13 +lon_0=31 +datum=WGS84 +units=m +no_defs"  # metres
sf_use_s2(FALSE)

save_fig <- function(p, name, w, h) {
  print(p)
  ggsave(file.path(out_dir, name), p, width = w, height = h, dpi = 200)
}

# ---- 1. Survey data --------------------------------------------------------------
v <- remo %>%
  filter(country == "SUDAN") %>%
  mutate(vid = row_number(), p = y / n, infected = y > 0)
v_sf   <- st_as_sf(v, coords = c("long", "lat"), crs = 4326, remove = FALSE)
v_sf_m <- st_transform(v_sf, aeqd)
cat(nrow(v), "villages\n")

# ---- 2. River network --------------------------------------------------------------
# DIS_AV_CMS = long-term average discharge (m3/s), modelled by HydroSHEDS.
# Reaches with flow below 0.1 m3/s (many with flow 0) are dropped.
bbox <- st_bbox(c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),
                crs = st_crs(4326))
rivers <- suppressMessages(st_filter(sudan_rivers, st_as_sfc(bbox), .predicate = st_intersects))
rivers <- rivers[rivers$DIS_AV_CMS >= 0.1, c("HYRIV_ID", "DIS_AV_CMS")]
rivers_m <- st_transform(rivers, aeqd)
cat(nrow(rivers), "river reaches with flow >= 0.1 m3/s\n")

# ---- 3. Aridity raster -------------------------------------------------------------
ai <- crop(unwrap(sudan_aridity), ext(xlim[1], xlim[2], ylim[1], ylim[2]))
if (global(ai, "max", na.rm = TRUE)[1, 1] > 100) ai <- ai / 10000   # stored as AI x 10,000
names(ai) <- "AI"

ai_breaks <- c(0, 0.05, 0.20, 0.50, Inf)
ai_labels <- c("Hyper-arid\n(AI < 0.05)", "Arid\n(0.05-0.20)",
               "Semi-arid\n(0.20-0.50)", "Sub-humid / humid\n(AI > 0.50)")
ai_cols   <- c("#a6611a", "#dfc27d", "#c7e9c0", "#41ab5d")

# ---- 4. Derived variables for each village -------------------------------------------
# (a) aridity: mean AI within 2 km (avoids missing values for villages on the river bank)
v_buf2 <- st_transform(st_buffer(v_sf_m, 2000), 4326)
v$AI <- terra::extract(ai, vect(v_buf2), fun = mean, na.rm = TRUE)[, "AI"]
v$ai_class <- cut(v$AI, ai_breaks, labels = ai_labels, include.lowest = TRUE)

# (b) distance (km) to the nearest river with flow >= threshold
dist_to_river <- function(th) {
  r <- rivers_m[rivers_m$DIS_AV_CMS >= th, ]
  as.numeric(st_distance(v_sf_m, r[st_nearest_feature(v_sf_m, r), ], by_element = TRUE)) / 1000
}
v$dist_small_km  <- dist_to_river(0.1)   # all streams with some flow
v$dist_medium_km <- dist_to_river(1)     # medium and major rivers
v$dist_major_km  <- dist_to_river(10)    # major rivers

# (c) rivers within 10 km: flow of the largest river, and length of river
buf10 <- st_buffer(v_sf_m[, "vid"], 10000)
near  <- rivers_m[lengths(st_intersects(rivers_m, st_union(buf10))) > 0, ]
pieces <- st_intersection(near, buf10)
pieces$len_km <- as.numeric(st_length(pieces)) / 1000
feat <- st_drop_geometry(pieces) %>%
  group_by(vid) %>%
  summarise(max_flow_10km = max(DIS_AV_CMS),
            len_flow1_10km = sum(len_km[DIS_AV_CMS >= 1]), .groups = "drop")
v <- v %>% left_join(feat, by = "vid") %>%
  mutate(max_flow_10km = coalesce(max_flow_10km, 0),
         len_flow1_10km = coalesce(len_flow1_10km, 0))

write.csv(v, "sudan_villages_explore.csv", row.names = FALSE)
print(summary(v[, c("AI", "dist_small_km", "dist_medium_km", "dist_major_km", "max_flow_10km")]))
print(table(v$ai_class, useNA = "ifany"))

# ---- 5. Layers for maps ----------------------------------------------------------------
sudans <- sudan_borders
ai_df <- as.data.frame(aggregate(ai, fact = 4, fun = "mean", na.rm = TRUE), xy = TRUE, na.rm = TRUE)
ai_df$ai_class <- cut(ai_df$AI, ai_breaks, labels = ai_labels, include.lowest = TRUE)
rivers_draw <- rivers[order(rivers$DIS_AV_CMS), ]                  # big rivers drawn last
major <- rivers[rivers$DIS_AV_CMS >= 10, ]

map_theme <- theme_minimal(base_size = 10) +
  theme(panel.background = element_rect(fill = "grey97", colour = NA))
village_layers <- list(
  geom_point(data = subset(v, !infected), aes(long, lat), shape = 4, size = 0.7, colour = "grey15"),
  geom_point(data = subset(v, infected), aes(long, lat), size = 1.1, colour = "red")
)

# Map 1: prevalence
m1 <- ggplot() +
  geom_sf(data = sudans, fill = "white", colour = "grey40", linewidth = 0.3) +
  geom_sf(data = major, colour = "#6baed6", linewidth = 0.4) +
  geom_point(data = subset(v, !infected), aes(long, lat), shape = 4, size = 0.8, colour = "grey30") +
  geom_point(data = subset(v, infected), aes(long, lat, colour = p), size = 1.5) +
  scale_colour_viridis_c(option = "inferno", direction = -1, limits = c(0, 0.5),
                         oob = scales::squish, name = "Nodule\nprevalence") +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "REMO nodule prevalence, Sudan and South Sudan",
       subtitle = "Crosses = no nodules found; blue = major rivers (flow >= 10 m³/s)",
       x = NULL, y = NULL) + map_theme
save_fig(m1, "map1_prevalence.png", 7.5, 8)

# Map 2: river network coloured by flow
m2 <- ggplot() +
  geom_sf(data = rivers_draw, aes(colour = DIS_AV_CMS, linewidth = DIS_AV_CMS)) +
  scale_colour_distiller(palette = "Blues", direction = 1, trans = "log10",
                         limits = c(0.1, 3000), oob = scales::squish,
                         name = "Average flow\n(m³/s, log scale)") +
  scale_linewidth(range = c(0.05, 0.9), trans = "log10", guide = "none") +
  geom_sf(data = sudans, fill = NA, colour = "grey30", linewidth = 0.3) +
  village_layers +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "River network coloured by average flow (HydroRIVERS)",
       subtitle = "Red = at least one infected person; crosses = no infection found",
       x = NULL, y = NULL) + map_theme
save_fig(m2, "map2_river_flow.png", 7.5, 8)

# Map 3: aridity with major rivers and villages
m3 <- ggplot() +
  geom_raster(data = ai_df, aes(x, y, fill = ai_class)) +
  scale_fill_manual(values = setNames(ai_cols, ai_labels), name = "Aridity index", drop = FALSE) +
  geom_sf(data = major, colour = "#2166ac", linewidth = 0.35) +
  geom_sf(data = sudans, fill = NA, colour = "grey20", linewidth = 0.3) +
  village_layers +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  labs(title = "Aridity index (precipitation / evapotranspiration), major rivers and villages",
       subtitle = "Red = at least one infected person; crosses = no infection found",
       x = NULL, y = NULL) + map_theme
save_fig(m3, "map3_aridity.png", 8, 8)

# ---- 6. Exploratory plots ------------------------------------------------------------------
# Helper: % of villages infected, by aridity class and a grouping variable
bar_by_aridity <- function(data, group, group_title, title, fills) {
  tab <- data %>%
    filter(!is.na(ai_class), !is.na(.data[[group]])) %>%
    group_by(ai_class, g = .data[[group]]) %>%
    summarise(villages = n(), infected = sum(infected), .groups = "drop") %>%
    filter(villages >= 5) %>%                                  # hide groups with < 5 villages
    mutate(pct = 100 * infected / villages,
           label = paste0(round(pct), "%\nn=", villages))
  print(tab, n = Inf)
  ggplot(tab, aes(x = ai_class, y = pct, fill = g)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    geom_text(aes(label = label, group = g), position = position_dodge(width = 0.85),
              vjust = -0.3, size = 2.6, lineheight = 0.85) +
    scale_fill_manual(values = fills, name = group_title, drop = FALSE) +
    scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 20)) +
    labs(x = "Aridity class of the village (driest on the left)",
         y = "Villages with at least one infected person (%)", title = title,
         caption = "Labels: % of villages infected and number of villages. Groups with fewer than 5 villages not shown.") +
    theme_minimal(base_size = 11) +
    theme(plot.caption = element_text(hjust = 0, size = 8), panel.grid.major.x = element_blank())
}
dist_fills <- c("Under 2 km" = "#1f4e79", "2-5 km" = "#3f7fb5", "5-15 km" = "#7fb2dc", "Over 15 km" = "#c6def1")
dist_cut <- function(x) cut(x, c(0, 2, 5, 15, Inf), include.lowest = TRUE,
                            labels = c("Under 2 km", "2-5 km", "5-15 km", "Over 15 km"))

# Plot 1: infected vs uninfected villages, distance to the nearest major river
v$status <- factor(ifelse(v$infected, "Infected\n(at least one nodule)", "Not infected\n(no nodules found)"),
                   levels = c("Not infected\n(no nodules found)", "Infected\n(at least one nodule)"))
meds <- v %>% group_by(status) %>% summarise(med = median(dist_major_km), .groups = "drop")
set.seed(1)
p1 <- ggplot(v, aes(status, dist_major_km)) +
  geom_jitter(aes(colour = ai_class), width = 0.28, height = 0, size = 1.2, alpha = 0.7) +
  geom_errorbar(data = meds, aes(x = status, ymin = med, ymax = med), inherit.aes = FALSE,
                width = 0.7, linewidth = 1) +
  geom_text(data = meds, aes(x = status, y = med, label = paste0("median\n", round(med), " km")),
            inherit.aes = FALSE, nudge_x = 0.48, size = 3.2, lineheight = 0.9) +
  scale_colour_manual(values = setNames(ai_cols, ai_labels), name = "Aridity class") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(x = NULL, y = "Distance from village to the nearest major river (km)",
       title = "Infected villages are, on average, further from major rivers than uninfected villages",
       caption = "Each dot is a village. Major river = average flow at least 10 m³/s.") +
  theme_minimal(base_size = 11) + theme(plot.caption = element_text(hjust = 0, size = 8))
save_fig(p1, "plot1_distance_by_infection.png", 8.5, 5.5)

# Plot 2: infection by aridity class and distance to the nearest major river
v$dist_major_group <- dist_cut(v$dist_major_km)
p2 <- bar_by_aridity(v, "dist_major_group", "Distance to nearest\nmajor river",
                     "Infection by aridity class and distance to the nearest major river",
                     dist_fills)
save_fig(p2, "plot2_infection_by_aridity_and_distance.png", 10, 5.5)

# Plot 3: the same for three definitions of "river"
long_d <- v %>%
  select(vid, infected, ai_class, dist_major_km, dist_medium_km, dist_small_km) %>%
  pivot_longer(starts_with("dist_"), names_to = "network", values_to = "dist") %>%
  mutate(network = factor(network,
                          levels = c("dist_major_km", "dist_medium_km", "dist_small_km"),
                          labels = c("Major rivers only (flow >= 10 m³/s)",
                                     "Medium and major rivers (flow >= 1 m³/s)",
                                     "All streams with some flow (flow >= 0.1 m³/s)")),
         dist_group = dist_cut(dist))
tab3 <- long_d %>%
  filter(!is.na(ai_class)) %>%
  group_by(network, ai_class, dist_group) %>%
  summarise(villages = n(), infected = sum(infected), .groups = "drop") %>%
  filter(villages >= 5) %>%
  mutate(pct = 100 * infected / villages, label = paste0(round(pct), "%\nn=", villages))
p3 <- ggplot(tab3, aes(ai_class, pct, fill = dist_group)) +
  geom_col(position = position_dodge(width = 0.85), width = 0.8) +
  geom_text(aes(label = label, group = dist_group), position = position_dodge(width = 0.85),
            vjust = -0.3, size = 2, lineheight = 0.85) +
  facet_wrap(~network, ncol = 1) +
  scale_fill_manual(values = dist_fills, name = "Distance to nearest\nriver of this type") +
  scale_y_continuous(limits = c(0, 115), breaks = seq(0, 100, 25)) +
  labs(x = "Aridity class of the village", y = "Villages with at least one infected person (%)",
       title = "Distance to the nearest river explains little, however rivers are defined") +
  theme_minimal(base_size = 10) + theme(panel.grid.major.x = element_blank())
save_fig(p3, "plot3_three_river_definitions.png", 10, 11)

# Plot 4: infection by aridity class and size of the largest river within 10 km
v$river_size <- cut(v$max_flow_10km, c(-Inf, 1, 10, 100, Inf),
                    labels = c("Under 1", "1-10", "10-100", "Over 100"))
p4 <- bar_by_aridity(v, "river_size", "Flow of the largest river\nwithin 10 km (m³/s)",
                     "Large rivers are not where infection is highest",
                     c("Under 1" = "#c6def1", "1-10" = "#7fb2dc", "10-100" = "#3f7fb5", "Over 100" = "#1f4e79"))
save_fig(p4, "plot4_infection_by_river_size.png", 10, 5.5)

# Plot 5: prevalence against aridity
p5 <- ggplot(v, aes(AI, p)) +
  geom_point(aes(size = n, colour = dist_major_km < 5), alpha = 0.5) +
  geom_smooth(method = "loess", span = 0.5, se = FALSE, colour = "black") +
  scale_x_log10(breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)) +
  scale_colour_manual(values = c(`TRUE` = "#2166ac", `FALSE` = "#d6604d"),
                      labels = c(`TRUE` = "Within 5 km of a major river", `FALSE` = "Further than 5 km"),
                      name = NULL) +
  scale_size_area(max_size = 3, guide = "none") +
  labs(x = "Aridity index at the village (log scale; smaller = drier)",
       y = "Observed nodule prevalence",
       title = "Observed prevalence against aridity") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
save_fig(p5, "plot5_prevalence_vs_aridity.png", 8.5, 5.5)

# Plot 6: how far rivers reach at each level of aridity
# Share of villages infected against aridity (in bins), separately for villages
# near to and far from a major river. The gap between the lines shows how much
# being close to a river matters at that level of aridity.
v$ai_bin <- cut(v$AI, unique(quantile(v$AI, seq(0, 1, 0.125), na.rm = TRUE)), include.lowest = TRUE)
tab6 <- v %>%
  filter(!is.na(ai_bin)) %>%
  mutate(near = ifelse(dist_major_km < 5, "Within 5 km of a major river", "Further than 5 km")) %>%
  group_by(ai_bin, near) %>%
  summarise(AI_mid = median(AI), villages = n(), share = mean(infected), .groups = "drop") %>%
  filter(villages >= 5)
p6 <- ggplot(tab6, aes(AI_mid, 100 * share, colour = near)) +
  geom_line() + geom_point(aes(size = villages)) +
  scale_x_log10(breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)) +
  scale_colour_manual(values = c("Within 5 km of a major river" = "#2166ac", "Further than 5 km" = "#d6604d"),
                      name = NULL) +
  scale_size_area(max_size = 4, name = "Villages") +
  labs(x = "Aridity index (median of each group of villages; log scale)",
       y = "Villages with at least one infected person (%)",
       title = "Infection rises steeply with wetness, whether or not a village is near a major river",
       caption = "Villages grouped into eight equal-sized groups by aridity; groups with fewer than 5 villages not shown.") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.caption = element_text(hjust = 0, size = 8))
save_fig(p6, "plot6_river_proximity_by_aridity.png", 8.5, 5.5)

cat("All figures written to", out_dir, "\n")
