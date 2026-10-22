###############################################################################
# Builds the datasets shipped with odispers (data/*.rda) from the original
# sources. Run from a folder that contains:
#   Remo_AllCountries.csv                         REMO nodule survey data
#   HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp HydroRIVERS v1.0, Africa
#   an annual Global Aridity Index v3.1 GeoTIFF (e.g. ai_v31_yr.tif)
# and with Natural Earth (rnaturalearth) available.
###############################################################################

library(sf)
library(terra)

out <- "data"    # run from the package root, or change this
xlim <- c(21, 40); ylim <- c(2.5, 23.5)    # study area for Sudan and South Sudan
bb <- st_bbox(c(xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]), crs = st_crs(4326))

# REMO: all countries, as in the original file
remo <- read.csv("Remo_AllCountries.csv")
save(remo, file = file.path(out, "remo.rda"), compress = "xz")

# HydroRIVERS reaches in the study area with long-term average discharge >= 0.1 m3/s
sudan_rivers <- st_read("HydroRIVERS_v10_af_shp/HydroRIVERS_v10_af.shp",
                        wkt_filter = st_as_text(st_as_sfc(bb)), quiet = TRUE)
sudan_rivers <- sudan_rivers[sudan_rivers$DIS_AV_CMS >= 0.1, c("HYRIV_ID", "NEXT_DOWN", "DIS_AV_CMS")]
save(sudan_rivers, file = file.path(out, "sudan_rivers.rda"), compress = "xz")

# Global Aridity Index v3.1, annual, cropped to the study area; integer values
# equal to the aridity index x 10,000, as in the original GeoTIFF
ai_file <- list.files(".", pattern = "^ai.*yr.*\\.tif$", recursive = TRUE, full.names = TRUE,
                      ignore.case = TRUE)[1]
a <- crop(rast(ai_file), ext(xlim[1], xlim[2], ylim[1], ylim[2]))
names(a) <- "AI_x10000"
sudan_aridity <- wrap(as.int(a))
save(sudan_aridity, file = file.path(out, "sudan_aridity.rda"), compress = "xz")

# Natural Earth country borders (1:50m)
b <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
sudan_borders <- b[b$admin %in% c("Sudan", "South Sudan"), "admin"]
save(sudan_borders, file = file.path(out, "sudan_borders.rda"), compress = "xz")
