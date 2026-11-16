# Colour DRC human-survey locations by observed positivity, separately for
# each diagnostic and source survey. Run prepare_drc_infection.R first.
# Outputs: self-contained HTML, three PNGs and site-summary GeoJSON in maps/.

library(sf)
library(ggplot2)
library(leaflet)
library(htmlwidgets)
library(htmltools)

root <- "data-raw/external_blackflies"
out <- file.path(root, "maps")
context <- file.path(root, "map_data")
rates <- read.csv(file.path(out, "drc_infection_proportions.csv"),
                  check.names = FALSE, encoding = "UTF-8")
locations <- read.csv(file.path(out, "mapped_locations.csv"),
                      check.names = FALSE, encoding = "UTF-8")
human <- locations[locations$data_type == "Human infection survey", ]
human <- st_as_sf(human, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
vectors <- locations[locations$data_type == "Fly / vector survey" &
                       locations$needs_review == FALSE, ]
vectors <- st_as_sf(vectors, coords = c("longitude", "latitude"), crs = 4326)
site_rates <- st_as_sf(rates, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
st_write(site_rates, file.path(out, "drc_infection_proportions.geojson"),
           delete_dsn = TRUE, quiet = TRUE)

countries <- st_read(file.path(context, "regional_countries.geojson"), quiet = TRUE)
lakes <- st_read(file.path(context, "regional_lakes.geojson"), quiet = TRUE)
rivers <- st_read(file.path(context, "drc_rivers.geojson"), quiet = TRUE)
warm_colours <- c("#F2C14E", "#E69B2E", "#D46632", "#B63B3A", "#7B1E4A")
positive_palette <- colorNumeric(warm_colours, domain = c(0, 1))
zero_colour <- "#151515"
missing_colour <- "#BDB9B0"
coordinate_key <- function(data) paste(data$latitude, data$longitude, sep = ":")
escape <- function(x) htmlEscape(ifelse(is.na(x), "", x))

survey_keys <- unique(paste(rates$diagnostic, rates$dataset_id, sep = "::"))
catalogue <- lapply(seq_along(survey_keys), function(number) {
  selected <- rates[paste(rates$diagnostic, rates$dataset_id, sep = "::") == survey_keys[number], ]
  list(group = paste0("human_", number), diagnostic = selected$diagnostic[1],
       dataset = selected$dataset_id[1], survey = selected$survey[1],
       sites = nrow(selected), tested = sum(selected$n_tested),
       positive = sum(selected$n_positive),
       community = selected$dataset_id[1] %in% c("7.1", "7.1.1", "8.1"))
})
default_group <- Filter(function(x) x$diagnostic == "Skin snip" && x$dataset == "2.1",
                         catalogue)[[1]]$group

popup <- function(data) paste0(
  "<b>", escape(data$site), "</b><br>",
  "<b>", data$n_positive, " / ", data$n_tested, " positive (",
  sprintf("%.1f", 100 * data$proportion_positive), "%)</b>",
  "<br><b>Test:</b> ", escape(data$diagnostic),
  "<br><b>Survey:</b> ", escape(data$survey),
  "<br><b>Design:</b> ", escape(data$survey_design),
  "<br><b>Coordinates:</b> ", sprintf("%.5f, %.5f", data$latitude, data$longitude),
  "<br>Proportion among tested source records; no missing results counted as negative.",
  '<br><a target="_blank" rel="noopener" href="https://doi.org/10.1371/journal.pntd.0014164">Source and dataset</a>')

map <- leaflet(width = "100%", height = "100vh",
                 options = leafletOptions(preferCanvas = TRUE, minZoom = 6, maxZoom = 17)) |>
  addPolygons(data = countries, fillColor = "#F4F1E9", fillOpacity = 1,
               color = "#B8B3A8", weight = 1, options = pathOptions(interactive = FALSE)) |>
  addPolygons(data = lakes, fillColor = "#D8E8EF", fillOpacity = 1,
               color = "#A9CBDC", weight = 0.6) |>
  addPolylines(data = rivers, color = "#559BC0", weight = 1.8,
                opacity = 0.85, label = ~RiverName, group = "Local rivers") |>
  addCircleMarkers(data = vectors, radius = 3, color = "#8A8177", weight = 1,
                     fillOpacity = 0, label = ~site, group = "Vector survey sites") |>
  hideGroup("Vector survey sites") |>
  addScaleBar(position = "bottomleft", options = scaleBarOptions(imperial = FALSE))

for (entry in catalogue) {
  selected <- site_rates[site_rates$diagnostic == entry$diagnostic &
                           site_rates$dataset_id == entry$dataset, ]
  untested <- human[!coordinate_key(human) %in% coordinate_key(selected), ]
  marker_colours <- ifelse(selected$proportion_positive == 0, zero_colour,
                            positive_palette(selected$proportion_positive))
  map <- addCircleMarkers(map, data = untested, radius = 4, color = "white", weight = 0.8,
                             fillColor = missing_colour, fillOpacity = 0.65,
                             popup = paste0("<b>", escape(untested$site),
                               "</b><br>No tested record at this coordinate for the selected survey/test."),
                             group = entry$group) |>
    addCircleMarkers(data = selected, radius = 7, color = "white", weight = 1.2,
                       fillColor = marker_colours, fillOpacity = 0.95,
                       label = paste0(selected$site, ": ", selected$n_positive, "/",
                         selected$n_tested, " (", round(100 * selected$proportion_positive, 1), "%)"),
                       popup = popup(selected), group = entry$group)
}
map <- hideGroup(map, vapply(catalogue, `[[`, "", "group")) |>
  showGroup(default_group) |>
  addLayersControl(overlayGroups = c("Local rivers", "Vector survey sites"),
                     options = layersControlOptions(collapsed = TRUE)) |>
  fitBounds(30.84, 1.90, 31.18, 2.24)

controls <- HTML(paste0(
  '<div class="infection-panel"><h1>DRC human infection surveys</h1>',
  '<p>Observed proportion positive at each reported location</p>',
  '<label>Diagnostic test<select data-diagnostic>',
  '<option>Skin snip</option><option>Ov16 antibodies</option><option>Nodules</option>',
  '</select></label><label>Survey<select data-survey></select></label>',
  '<p class="survey-summary"></p><p class="sampling-note"></p>',
  '<p class="map-note">Click a point for positives / tested. Different tests and surveys are kept separate.</p>',
  '<button data-reset>Reset to Kakoi–Koda</button></div>'))
legend <- HTML(paste0(
  '<div class="infection-legend"><b>Observed positive</b>',
  '<div class="gradient"></div><div class="ticks"><span>&gt;0%</span>',
  '<span>25%</span><span>50%</span><span>75%</span><span>100%</span></div>',
  '<div><i style="background:#151515"></i>0% positive</div>',
  '<div><i style="background:#BDB9B0"></i>No selected test result</div></div>'))
map <- addControl(map, controls, position = "topleft") |>
  addControl(legend, position = "bottomright")
css <- tags$style(HTML(paste0(
  'html,body{margin:0;height:100%;font-family:Arial,sans-serif;color:#28343b}',
  '.infection-panel,.infection-legend{background:rgba(255,255,255,.96);padding:15px 17px;',
  'border-radius:9px;box-shadow:0 2px 14px #0002}',
  '.infection-panel{width:330px}.infection-panel h1{font-size:20px;margin:0 0 7px}',
  '.infection-panel p{font-size:12px;line-height:1.5;margin:5px 0 10px}',
  '.infection-panel label{display:block;font-size:12px;font-weight:600;margin:10px 0}',
  '.infection-panel select{display:block;width:100%;margin-top:4px;padding:7px;border:1px solid #c9c6bf;',
  'border-radius:5px;background:white;color:#26343c;font-size:12px}',
  '.infection-panel .sampling-note{font-size:11px;background:#f5efe2;padding:8px;border-radius:4px}',
  '.infection-panel .map-note{font-size:11px;color:#666}',
  '.infection-panel button{background:#f4f1e9;border:1px solid #d5d2c9;padding:7px 10px;',
  'border-radius:5px;cursor:pointer}.infection-legend{font-size:12px;width:230px}',
  '.infection-legend .gradient{height:12px;margin-top:8px;background:linear-gradient(to right,',
  paste(warm_colours, collapse = ","), ')}',
  '.infection-legend .ticks{display:flex;justify-content:space-between;font-size:10px;margin:3px 0 8px}',
  '.infection-legend i{display:inline-block;width:10px;height:10px;border-radius:50%;margin:4px 8px 0 0}',
  '.leaflet-popup-content{font-size:12px;line-height:1.6;max-width:340px}',
  '@media(max-width:700px){.infection-panel{width:245px;padding:10px 12px}',
  '.infection-panel h1{font-size:16px}.infection-legend{width:180px}}')))
map <- prependContent(map, css)
map <- onRender(map, "function(el,x,data){
  var map=this, diagnostic=el.querySelector('[data-diagnostic]'),
      survey=el.querySelector('[data-survey]'), entries=data.catalogue,
      preferred={'Skin snip':'2.1','Ov16 antibodies':'7.1','Nodules':'2.1'};
  L.DomEvent.disableClickPropagation(el.querySelector('.infection-panel'));
  L.DomEvent.disableScrollPropagation(el.querySelector('.infection-panel'));
  function show(){
    entries.forEach(function(entry){map.layerManager.getLayerGroup(entry.group).removeFrom(map);});
    var entry=entries.find(function(e){return e.group===survey.value;});
    map.layerManager.getLayerGroup(entry.group).addTo(map);
    el.querySelector('.survey-summary').textContent=entry.sites+' locations · '+entry.positive+
      ' positive / '+entry.tested+' tested at mapped coordinates';
    el.querySelector('.sampling-note').textContent=entry.community?
      'Community survey. This selection describes the sampled ages and survey year. Ov16 measures antibodies, not current parasite positivity.':
      'Selected or case-control records. These observed proportions are not estimates of community prevalence.';
  }
  function chooseDiagnostic(){
    survey.replaceChildren();
    var available=entries.filter(function(e){return e.diagnostic===diagnostic.value;});
    available.forEach(function(e){var option=document.createElement('option');
      option.value=e.group;option.textContent=e.survey;survey.appendChild(option);});
    var initial=available.find(function(e){return e.dataset===preferred[diagnostic.value];});
    if(initial) survey.value=initial.group;
    show();
  }
  diagnostic.addEventListener('change',chooseDiagnostic);
  survey.addEventListener('change',show);
  el.querySelector('[data-reset]').addEventListener('click',function(){
    map.fitBounds([[1.9,30.84],[2.24,31.18]],{paddingTopLeft:[380,40],paddingBottomRight:[30,70]});});
  chooseDiagnostic();
}", data = list(catalogue = catalogue))
saveWidget(map, file.path(out, "drc_human_infection.html"), selfcontained = TRUE,
             title = "DRC human infection proportions")

# Static versions use the same point-level numerators and denominators.
plot_survey <- function(diagnostic, dataset_id, filename, title) {
  selected <- site_rates[site_rates$diagnostic == diagnostic & site_rates$dataset_id == dataset_id, ]
  other <- human[!coordinate_key(human) %in% coordinate_key(selected), ]
  positive <- selected[selected$proportion_positive > 0, ]
  zero <- selected[selected$proportion_positive == 0, ]
  labels <- st_drop_geometry(selected)
  labels$text <- sprintf("%.0f%%", 100 * labels$proportion_positive)
  if (nrow(labels) > 16) labels <- labels[FALSE, ]
  subtitle <- paste0(selected$survey[1], "\n", sum(selected$n_positive), " positive / ",
                      sum(selected$n_tested), " tested at ", nrow(selected), " mapped locations")
  plot <- ggplot() +
    geom_sf(data = countries, fill = "#F4F1E9", colour = "#B8B3A8", linewidth = 0.3) +
    geom_sf(data = lakes, fill = "#D8E8EF", colour = "#A9CBDC", linewidth = 0.2) +
    geom_sf(data = rivers, colour = "#6CA6C4", linewidth = 0.5) +
    geom_sf(data = other, fill = missing_colour, colour = "white", shape = 21,
             size = 2.5, stroke = 0.5, alpha = 0.55) +
    geom_sf(data = positive, aes(fill = proportion_positive), colour = "white",
             shape = 21, size = 4.5, stroke = 0.8) +
    geom_sf(data = zero, fill = zero_colour, colour = "white", shape = 21,
             size = 4.5, stroke = 0.8) +
    ggrepel::geom_text_repel(data = labels, aes(longitude, latitude, label = text),
                            size = 3.5, seed = 22, min.segment.length = 0,
                            box.padding = 0.5, point.padding = 0.6, max.overlaps = Inf) +
    scale_fill_gradientn(colours = warm_colours, limits = c(0, 1),
                          breaks = c(0, .25, .5, .75, 1),
                          labels = c(">0%", "25%", "50%", "75%", "100%"),
                          name = "Observed positive") +
    coord_sf(xlim = c(30.84, 31.18), ylim = c(1.90, 2.24), expand = FALSE) +
    labs(title = title, subtitle = subtitle,
         caption = paste0("Black: 0% positive. Grey: no selected-survey result. Values are positives / tested at each location.\n",
           if (dataset_id %in% c("7.1", "7.1.1", "8.1"))
             "Community survey; Ov16 detects antibodies. Source: Amaral et al. (2026), S1 Data."
           else "Selected/case-control records; not community prevalence. Source: Amaral et al. (2026), S1 Data.")) +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_blank(), panel.grid.major = element_line(colour = "#E6E3DB", linewidth = 0.25),
          plot.title = element_text(face = "bold", size = 20),
          plot.subtitle = element_text(colour = "#555555", lineheight = 1.3),
          plot.caption = element_text(hjust = 0, colour = "#555555", size = 9),
          legend.position = "bottom", legend.key.width = grid::unit(1.7, "cm"),
          plot.margin = margin(18, 20, 14, 20))
  ggsave(file.path(out, filename), plot, width = 10, height = 10, dpi = 220,
          device = ragg::agg_png)
}
plot_survey("Skin snip", "2.1", "drc_human_skin_snip.png", "DRC: skin-snip positivity")
plot_survey("Ov16 antibodies", "7.1", "drc_human_ov16_2016.png", "DRC: Ov16 antibody positivity")
plot_survey("Nodules", "2.1", "drc_human_nodules.png", "DRC: nodule positivity")
cat("Saved DRC positivity map with", length(catalogue), "separate test/survey selections.\n")
