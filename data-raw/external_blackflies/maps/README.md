# Rivers and locations

- [DRC human infection proportions](drc_human_infection.html): colour the human
  survey locations by positives / tested, with separate diagnostic and survey
  selectors. See [interpretation and denominators](#drc-human-infection-proportions).
- [Interactive map](rivers_and_locations.html): a self-contained HTML file with
  study layers, site/source popups and buttons for the regional view, DRC/Uganda,
  Kakoi–Koda, Sudan/Ethiopia and South Sudan. River and country geometry is
  embedded, so the map does not require an online basemap.
- [Regional PNG](rivers_and_locations.png).
- [DRC/Uganda PNG](drc_uganda_rivers_and_locations.png).
- [Kakoi–Koda detailed PNG](drc_kakoi_koda_detail.png).
- [Map locations, CSV](mapped_locations.csv) and [GeoJSON](locations.geojson).

## What is mapped

There are 234 site entries after combining repeated source observations at the
same coordinates within each study and data type. The default map shows 231
entries. Three reported-coordinate outliers are available in an initially hidden
review layer and in [coordinates_needing_review.csv](coordinates_needing_review.csv).
Site entries in different studies can overlap; these counts are not counts of
independent villages or participants.

## DRC human infection proportions

The [DRC positivity map](drc_human_infection.html) has 12 separate diagnostic /
survey selections. The diagnostic options are skin-snip positivity, Ov16
antibody positivity and palpable-nodule positivity. Each marker is
`number of source records with result 1 / number with result 0 or 1` at that
reported coordinate within the selected dataset. Missing results do not enter
the denominator. Black is zero positive; the yellow-to-burgundy scale represents
positive proportions; grey means no tested record for the selected combination.
Popups include the numerator, denominator, percentage, survey design and source.

The default skin-snip view uses the 2015 mapping dataset 2.1: 54 positives out
of 133 tested records at 10 georeferenced locations. The full source contains
173 skin-snip-tested records in that dataset; 40 lack usable coordinates and
are not placed on the map. This is selected/case-control material, so its
observed proportions should not be described as community prevalence.

The Ov16 selector defaults to the 2016 community survey 7.1: 280 positive out
of 915 tested records at 29 georeferenced coordinates. Ov16 positivity measures
antibodies, not current parasite detection. The nodule selector defaults to
the 2015 mapping data: 13 positives out of 55 tested records at six coordinates.
Different ages, diagnostic methods, years and study designs are not pooled.
Datasets 1.1 and 2.1 explicitly overlap in the source notes; they are alternative
selections and are never added together.

Static maps: [skin snip](drc_human_skin_snip.png),
[2016 Ov16](drc_human_ov16_2016.png), [nodules](drc_human_nodules.png).
The [site summaries](drc_infection_proportions.csv),
[GeoJSON](drc_infection_proportions.geojson) and
[complete denominator audit](drc_infection_denominators.csv) are saved alongside.
Only aggregated site results are added to the map. These are descriptive
summaries of the supplied records; no statistical model was fitted.

Rebuild from the repository root:

```sh
Rscript data-raw/external_blackflies/prepare_drc_infection.R
RSTUDIO_PANDOC=/Applications/quarto/bin/tools/aarch64 Rscript data-raw/external_blackflies/map_drc_infection.R
```

Orange markers show fly/vector survey locations, including larval searches and
crab trapping. Purple shows human infection survey locations. Green shows GBIF
presence records. The colours do not encode prevalence or abundance. A survey
marker does not imply a positive result.

Coordinates are taken from downloaded source tables and converted from
degrees/minutes where necessary. Human records are combined at their reported
coordinates and only site/survey metadata are included in the map. No sites were
geocoded, jittered or snapped to rivers. Source spelling and source-country
labels remain; older GBIF records labelled Sudan can lie in present-day South
Sudan.

## Coordinate checks and omissions

- Mograt's latitude is printed as 91 degrees in the 2016 Sudan paper. It is
  outside the geographic range and is omitted without correction.
- A DRC crab-trapping site's reported longitude, 30.05022, is an outlier relative
  to the other Kakoi–Koda sites and is held for review.
- Two cross-border-survey locations reported as DRC lie more than 10 km outside
  DRC on the Natural Earth borders. Their original positions are held in the
  review layer. This is a coordinate check, not a proposed correction.
- 659 rows of the DRC epidemiological CSV have no usable coordinates. They are
  excluded before combining repeated records.
- The downloaded tables do not provide coordinates for the two Chad sites,
  11 Merowe sites, five Maridi villages or 16 western-Uganda sentinel villages.
  They are not placed on the map from place names alone.
- The detailed list is in [unmapped_locations.csv](unmapped_locations.csv).

## River context and sources

Major rivers, lakes and borders are from [Natural Earth](https://www.naturalearthdata.com/downloads/).
The interactive map also embeds 3,258 HydroRIVERS reaches near survey locations,
from the package's existing `data/sudan_rivers.rda`. These cover the package's
Sudan/South Sudan study extent and are simplified by 80 m for display.

The DRC close-ups embed 190 local river/stream features from the river layer
published with [Amaral et al. (2026)](https://doi.org/10.1371/journal.pntd.0014164).
The [source web map](https://www.arcgis.com/apps/mapviewer/index.html?webmap=43d179502a234f0dbc38af712ad68190)
and river metadata are retained in `../map_data/`. The article describes the
river data as derived from OpenStreetMap and adjusted with field waypoints.
Credit OpenStreetMap contributors when reusing that layer.

The map layers provide geographic context; they are not a current-flow inventory
or a complete catalogue of vector breeding habitat. Exact download URLs and
checksums are in [map_data/sources.csv](../map_data/sources.csv). Study sources
are linked in [the dataset README](../README.md).

## Rebuild and verification

Run from the repository root:

```sh
python3 data-raw/external_blackflies/prepare_map_locations.py
RSTUDIO_PANDOC=/Applications/quarto/bin/tools/aarch64 Rscript data-raw/external_blackflies/map_locations.R
```

On another computer, set `RSTUDIO_PANDOC` to its Pandoc directory or put Pandoc
on `PATH`. R packages used: sf, ggplot2, leaflet, htmlwidgets, htmltools, ragg and
ggrepel. The source files in `map_data/` are retained; rebuilding needs no new
network downloads.

The PNGs were rendered and inspected. Embedded site counts, source-table row
counts, GeoJSON and JavaScript syntax were checked. Local browser access was
disabled in Work settings, so interactive button/layer behaviour could not be
tested in a browser.
