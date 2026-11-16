# Map the downloaded locations and river context. Run from the package root:
#   python3 data-raw/external_blackflies/prepare_map_locations.py
#   Rscript data-raw/external_blackflies/map_locations.R
# If Pandoc is not on PATH, set RSTUDIO_PANDOC to its directory first.
# Inputs remain unchanged. Outputs are a self-contained HTML map, PNG maps,
# site GeoJSON and cropped river/country GeoJSON files in maps/ and map_data/.

library(sf)
library(ggplot2)
library(leaflet)
library(htmlwidgets)
library(htmltools)

root <- "data-raw/external_blackflies"
out <- file.path(root, "maps")
context <- file.path(root, "map_data")
sf_use_s2(FALSE)

locations <- read.csv(file.path(out, "locations.csv"), check.names = FALSE,
                      na.strings = "", encoding = "UTF-8")
locations$needs_review <- tolower(as.character(locations$needs_review)) == "true"
sites <- st_as_sf(locations, coords = c("longitude", "latitude"), crs = 4326,
                  remove = FALSE)

read_context <- function(name) {
  directory <- file.path(context, paste0("natural_earth_", name))
  file <- list.files(directory, pattern = "\\.shp$", full.names = TRUE)
  st_read(file, quiet = TRUE)
}

region <- st_bbox(c(xmin = 18, ymin = -10.5, xmax = 41, ymax = 21.5), crs = 4326)
countries <- suppressWarnings(st_crop(st_make_valid(read_context("countries")), region))
rivers <- suppressWarnings(st_crop(read_context("rivers"), region))
lakes <- suppressWarnings(st_crop(st_make_valid(read_context("lakes")), region))
drc_rivers <- st_read(file.path(context, "drc_rivers.geojson"), quiet = TRUE)

# Retain source coordinates but separate clear country/position conflicts.
# A 10 km margin allows for generalized borders and border-adjacent sites.
country_keys <- c(DRC = "Democratic Republic of the Congo", CAR = "Central African Republic",
                  Uganda = "Uganda", Kenya = "Kenya", Sudan = "Sudan",
                  "South Sudan" = "South Sudan", Ethiopia = "Ethiopia",
                  Tanzania = "United Republic of Tanzania")
for (country in unique(sites$country)) {
  index <- which(sites$country == country)
  expected <- countries[countries$ADMIN == country_keys[[country]], ]
  if (country == "Sudan") {
    # Museum records can use Sudan's historical boundaries.
    expected <- countries[countries$ADMIN %in% c("Sudan", "South Sudan"), ]
  }
  distance <- as.numeric(st_distance(st_transform(sites[index, ], 3857),
                                     st_union(st_transform(expected, 3857))))
  conflicting <- index[distance > 10000]
  if (length(conflicting)) {
    sites$needs_review[conflicting] <- TRUE
    sites$coordinate_note[conflicting] <- paste0(sites$coordinate_note[conflicting],
      "; reported coordinates lie >10 km outside the source country on Natural Earth borders")
  }
}
usable <- sites[!sites$needs_review, ]
review <- sites[sites$needs_review, ]
write.csv(st_drop_geometry(sites), file.path(out, "mapped_locations.csv"),
           row.names = FALSE, na = "")
write.csv(st_drop_geometry(review), file.path(out, "coordinates_needing_review.csv"),
           row.names = FALSE, na = "")

# HydroRIVERS already ships with the package. Show its finer reaches near the
# downloaded vector-survey locations, rather than embedding all 93,094 reaches.
river_data <- new.env()
load("data/sudan_rivers.rda", envir = river_data)
survey_sites <- usable[usable$data_type == "Fly / vector survey", ]
context_windows <- st_transform(st_buffer(st_transform(survey_sites, 3857), 45000), 4326)
near_survey <- lengths(st_intersects(river_data$sudan_rivers, context_windows)) > 0
local_rivers <- river_data$sudan_rivers[near_survey, ]
local_rivers <- st_transform(st_simplify(st_transform(local_rivers, 3857),
                                        dTolerance = 80), 4326)

# Geographic inputs are saved with the map so the offline result is inspectable.
for (name in c("countries", "rivers", "lakes", "local_rivers")) {
  st_write(get(name), file.path(context, paste0("regional_", name, ".geojson")),
           delete_dsn = TRUE, quiet = TRUE)
}
st_write(sites, file.path(out, "locations.geojson"), delete_dsn = TRUE, quiet = TRUE)

colours <- c("Fly / vector survey" = "#C15319",
             "Human infection survey" = "#714A99",
             "GBIF presence record" = "#267A60")
land <- "#F4F1E9"
water <- "#6B9FBF"

country_labels <- data.frame(
  country = c("SUDAN", "SOUTH SUDAN", "ETHIOPIA", "UGANDA", "KENYA", "DRC",
              "CENTRAL AFRICAN\nREPUBLIC", "CHAD", "TANZANIA", "ERITREA"),
  longitude = c(28, 29, 37.4, 32.8, 38, 25.8, 23, 20.5, 35.4, 38.5),
  latitude = c(16.1, 8.3, 9.0, 1.1, -1.4, -0.4, 5.2, 15.3, -6.8, 16.2))

base_plot <- function() {
  ggplot() +
    geom_sf(data = countries, fill = land, colour = "#B8B3A8", linewidth = 0.35) +
    geom_sf(data = lakes, fill = "#D8E8EF", colour = "#A9CBDC", linewidth = 0.15) +
    geom_sf(data = rivers, colour = water, linewidth = 0.32, alpha = 0.8) +
    scale_colour_manual(values = colours, name = NULL, drop = FALSE) +
    scale_shape_manual(values = c("Fly / vector survey" = 16,
                                  "Human infection survey" = 17,
                                  "GBIF presence record" = 15), name = NULL) +
    theme_minimal(base_family = "sans", base_size = 12) +
    theme(panel.grid.major = element_line(colour = "#E6E3DB", linewidth = 0.3),
          panel.background = element_rect(fill = "#F7FAFB", colour = NA),
          plot.background = element_rect(fill = "white", colour = NA),
          legend.position = "bottom", legend.box = "vertical",
          plot.title = element_text(face = "bold", size = 20),
          plot.subtitle = element_text(colour = "#555555", margin = margin(b = 12)),
          plot.caption = element_text(colour = "#666666", hjust = 0, size = 9),
          axis.title = element_blank(), plot.margin = margin(18, 20, 15, 20))
}

regional <- base_plot() +
  geom_text(data = country_labels, aes(longitude, latitude, label = country),
             size = 3.3, colour = "#7B776E", fontface = "bold", lineheight = 0.95) +
  geom_sf(data = usable, aes(colour = data_type, shape = data_type), size = 2,
          alpha = 0.8, stroke = 0.2) +
  annotate("rect", xmin = 30.45, xmax = 31.68, ymin = 1.58, ymax = 3.14,
           colour = "#3E4852", fill = NA, linewidth = 0.5, linetype = 2) +
  annotate("text", x = 30.25, y = 3.3, label = "DRC / Uganda close-up", hjust = 1,
           size = 3.0, colour = "#3E4852") +
  coord_sf(xlim = c(18, 41), ylim = c(-10.5, 21.5), expand = FALSE) +
  labs(title = "Rivers and downloaded survey locations",
       subtitle = "Sudan, South Sudan and neighboring countries",
       caption = paste0("One marker per reported coordinate within each study. Colours indicate data type, not prevalence.\n",
                        "Sources: downloaded studies / GBIF; river and boundary context: Natural Earth.\n",
                        "Missing coordinates and flagged outliers are omitted. Downloaded 31 October 2026."))
ggsave(file.path(out, "rivers_and_locations.png"), regional,
       width = 10, height = 12, dpi = 220, device = ragg::agg_png)

drc_box <- st_bbox(c(xmin = 30.45, ymin = 1.58, xmax = 31.68, ymax = 3.14), crs = 4326)
drc_sites <- usable[lengths(st_intersects(usable, st_as_sfc(drc_box))) > 0, ]
drc_detail <- suppressWarnings(st_crop(drc_rivers, drc_box))

# Label the longest mapped part of each named local river; all source lines
# remain visible, including unnamed streams.
river_names <- unique(na.omit(drc_detail$RiverName))
river_label_points <- lapply(river_names, function(name) {
  pieces <- st_cast(st_transform(drc_detail[drc_detail$RiverName == name &
                                            !is.na(drc_detail$RiverName), ], 32636),
                    "LINESTRING", warn = FALSE)
  longest <- pieces[which.max(st_length(pieces)), ]
  point <- st_line_sample(longest, sample = 0.5)
  st_sf(river = name, geometry = st_cast(st_transform(point, 4326), "POINT"))
})
river_label_points <- do.call(rbind, river_label_points)
river_label_xy <- cbind(st_drop_geometry(river_label_points), st_coordinates(river_label_points))

drc_map <- base_plot() +
  geom_sf(data = local_rivers, colour = "#8FB3C6", linewidth = 0.2, alpha = 0.75) +
  geom_sf(data = drc_detail, colour = "#4C94BA", linewidth = 0.5, alpha = 0.9) +
  geom_sf(data = drc_sites, aes(colour = data_type, shape = data_type), size = 2.2,
          alpha = 0.8) +
  ggrepel::geom_text_repel(data = river_label_xy, aes(X, Y, label = river),
                          colour = "#39708E", size = 3.0, fontface = "italic",
                          max.overlaps = 12, seed = 22, min.segment.length = 0) +
  annotate("text", x = 30.58, y = 2.42, label = "DRC", colour = "#8A857B",
           size = 5, fontface = "bold") +
  annotate("text", x = 31.52, y = 2.68, label = "UGANDA", colour = "#8A857B",
           size = 4, fontface = "bold") +
  annotate("text", x = 31.49, y = 1.77, label = "Lake Albert", colour = "#638DA5",
           size = 4, fontface = "italic", angle = 55) +
  coord_sf(xlim = c(30.45, 31.68), ylim = c(1.58, 3.14), expand = FALSE) +
  labs(title = "DRC / Uganda: rivers and survey locations",
       subtitle = "Kakoi–Koda focus and adjacent cross-border surveys",
       caption = paste("Local streams: Amaral et al. (2026) / HydroRIVERS; broader context: Natural Earth.",
                       "\nFly/vector surveys include river prospections and crab trapping. Human markers are survey locations.",
                       "\nReported coordinates are not snapped to rivers. Use the interactive map for site names and sources."))
ggsave(file.path(out, "drc_uganda_rivers_and_locations.png"), drc_map,
       width = 9, height = 11, dpi = 220, device = ragg::agg_png)

focus_map <- base_plot() +
  geom_sf(data = drc_detail, colour = "#4C94BA", linewidth = 0.6, alpha = 0.9) +
  geom_sf(data = drc_sites, aes(colour = data_type, shape = data_type),
          size = 2.5, alpha = 0.8) +
  ggrepel::geom_text_repel(data = river_label_xy, aes(X, Y, label = river),
                          colour = "#39708E", size = 3.4, fontface = "italic",
                          max.overlaps = 12, seed = 22, min.segment.length = 0) +
  coord_sf(xlim = c(30.84, 31.18), ylim = c(1.9, 2.24), expand = FALSE) +
  labs(title = "Kakoi–Koda: local rivers and survey sites",
       subtitle = "Ituri, Democratic Republic of the Congo",
       caption = paste0("Rivers: layer published with Amaral et al. (2026). Sites: Post et al. (2022) / Amaral et al. (2026).\n",
                        "Colours indicate data type, not prevalence. Reported coordinates are not snapped to rivers."))
ggsave(file.path(out, "drc_kakoi_koda_detail.png"), focus_map,
       width = 10, height = 10, dpi = 220, device = ragg::agg_png)

# Interactive layers retain study identity and show the source on click.
escape <- function(value) htmlEscape(ifelse(is.na(value), "", value))
popup <- function(data) paste0(
  "<b>", escape(data$site), "</b><br>", escape(data$study),
  "<br><b>Type:</b> ", escape(data$data_type),
  "<br><b>Country:</b> ", escape(data$country),
  "<br><b>Period:</b> ", escape(data$period),
  "<br><b>River:</b> ", escape(data$river),
  "<br><b>Source rows at this position:</b> ", data$source_rows,
  "<br><b>Coordinates:</b> ", sprintf("%.5f, %.5f", data$latitude, data$longitude),
  "<br>", escape(data$coordinate_note),
  ifelse(is.na(data$source_doi) | data$source_doi == "", "", paste0(
    '<br><a target="_blank" rel="noopener" href="https://doi.org/',
    escape(data$source_doi), '">Source publication</a>')))

map <- leaflet(width = "100%", height = "100vh",
                 options = leafletOptions(preferCanvas = TRUE, minZoom = 3, maxZoom = 17)) |>
  addPolygons(data = countries, fillColor = land, fillOpacity = 1,
               color = "#ADA99F", weight = 1, group = "Country boundaries",
               label = ~NAME, options = pathOptions(interactive = FALSE)) |>
  addPolygons(data = lakes, fillColor = "#D8E8EF", fillOpacity = 0.9,
               color = "#A2C5D6", weight = 0.5, group = "Major rivers and lakes") |>
  addPolylines(data = rivers, color = water, weight = 1.5, opacity = 0.9,
                group = "Major rivers and lakes", label = ~name) |>
  addPolylines(data = local_rivers, color = "#5999BA", weight = 1, opacity = 0.65,
                group = "Detailed rivers near survey sites") |>
  addPolylines(data = drc_rivers, color = "#3B8AB4", weight = 1.8, opacity = 0.9,
                group = "Detailed rivers near survey sites", label = ~RiverName) |>
  addScaleBar(position = "bottomleft", options = scaleBarOptions(imperial = FALSE))

studies <- unique(usable$study)
for (study in studies) {
  data <- usable[usable$study == study, ]
  map <- addCircleMarkers(map, data = data, radius = 5,
                            fillColor = unname(colours[data$data_type]),
                            fillOpacity = 0.8, color = "white", weight = 1,
                            popup = popup(data), label = data$site, group = study)
}
map <- addCircleMarkers(map, data = review, radius = 7, color = "#C62828",
                         fillColor = "#C62828", fillOpacity = 0.7, weight = 2,
                         popup = popup(review), group = "Reported coordinates needing review") |>
  hideGroup("Reported coordinates needing review") |>
  addLayersControl(overlayGroups = c(studies, "Major rivers and lakes",
                                    "Detailed rivers near survey sites",
                                    "Reported coordinates needing review"),
                     options = layersControlOptions(collapsed = FALSE)) |>
  addLegend(position = "bottomright", colors = colours, labels = names(colours),
              title = "Location type", opacity = 0.9) |>
  fitBounds(18.5, -9.5, 39.5, 20.2)

header <- HTML(paste0(
  '<div class="map-title"><h1>Rivers &amp; survey locations</h1>',
  '<p>Downloaded blackfly and onchocerciasis data</p>',
  '<div class="map-buttons"><button data-view="regional">Regional view</button>',
  '<button data-view="drc">DRC / Uganda</button>',
  '<button data-view="focus">Kakoi–Koda detail</button>',
  '<button data-view="atbara">Sudan / Ethiopia</button>',
  '<button data-view="south">South Sudan</button></div>',
  '<p class="map-note">Click a marker for its site, period and source. Colours show data type, not prevalence.<br>',
  'Missing coordinates are omitted. A reported outlier is available in the review layer.</p></div>'))
map <- addControl(map, header, position = "topleft")
style <- tags$style(HTML(paste0(
  'html,body{margin:0;height:100%;font-family:Arial,sans-serif}',
  '.map-title{background:rgba(255,255,255,.96);padding:14px 17px;max-width:405px;',
  'border-radius:9px;box-shadow:0 2px 14px #0002;color:#26343c}',
  '.map-title h1{font-size:20px;margin:0 0 5px;font-weight:700}',
  '.map-title p{margin:4px 0 10px;font-size:12px;color:#52616b}',
  '.map-buttons{display:flex;gap:5px;flex-wrap:wrap}',
  '.map-buttons button{background:#f4f1e9;border:1px solid #d5d2c9;',
  'border-radius:5px;padding:6px 8px;color:#26343c;cursor:pointer;font-size:12px}',
  '.map-buttons button:hover{background:#e7edf0}',
  '.map-title .map-note{font-size:11px;line-height:1.5;margin-top:10px;margin-bottom:0}',
  '.leaflet-control-layers{font-family:Arial,sans-serif;font-size:11px;line-height:1.5;',
  'border:0!important;border-radius:8px!important;box-shadow:0 2px 14px #0002!important;',
  'max-height:68vh;overflow:auto}',
  '.leaflet-popup-content{font-size:12px;line-height:1.6;max-width:330px}',
  '.leaflet-control-attribution{font-size:9px}',
  '@media(max-width:800px){.map-title{max-width:270px}.map-title h1{font-size:16px}',
  '.leaflet-control-layers{max-width:210px}}')))
map <- prependContent(map, style)
map <- onRender(map, "function(el,x){
  var map=this;
  var views={regional:[[ -9.5,18.5],[20.2,39.5]],
             drc:[[1.58,30.45],[3.14,31.68]],
             focus:[[1.9,30.84],[2.24,31.18]],
             atbara:[[11.3,35.4],[15.0,37.4]],
             south:[[4.6,29.6],[6.3,30.8]]};
  el.querySelectorAll('[data-view]').forEach(function(button){
    L.DomEvent.disableClickPropagation(button);
    button.addEventListener('click',function(){map.fitBounds(views[button.dataset.view],
      {paddingTopLeft:[25,175],paddingBottomRight:[260,70]});});
  });
}")
map <- addControl(map, HTML('<span>Rivers: Natural Earth / HydroRIVERS / Kakoi-Koda study. Locations: downloaded studies / GBIF.</span>'),
                    position = "bottomleft", className = "map-source")
saveWidget(map, file.path(out, "rivers_and_locations.html"), selfcontained = TRUE,
             title = "Rivers and survey locations")
cat("Saved", nrow(usable), "site entries,", nrow(local_rivers),
    "local HydroRIVERS reaches and", nrow(drc_rivers), "DRC river features.\n")
