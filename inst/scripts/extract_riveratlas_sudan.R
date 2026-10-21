# Extract RiverATLAS attributes for Sudan and South Sudan
#
# 1. Download RiverATLAS (geodatabase, 2.3 GB) from
#      https://figshare.com/ndownloader/files/20087321
#    and unzip it.
# 2. Set `src` below to the unzipped .gdb folder (or to the Africa
#    shapefile, RiverATLAS_v10_af.shp, if you took the shapefile version).
# 3. Run this script. It writes riveratlas_sudan.csv (no geometry, roughly
#    20 MB), which links to our river network through HYRIV_ID. Upload that
#    file to the chat.

library(sf)

src <- "RiverATLAS_v10.gdb"          # <- edit this path

# Same extent as the HydroRIVERS network already used in the model
bbox <- st_as_sfc(st_bbox(c(xmin = 20.7, ymin = 2.1, xmax = 40.3, ymax = 23.7),
                          crs = st_crs(4326)))

layer <- st_layers(src)$name[1]
ra <- st_read(src, layer = layer, wkt_filter = st_as_text(bbox), quiet = TRUE)
ra <- st_drop_geometry(ra)

# Columns to keep. Matching is by pattern, so names absent from your version
# are skipped rather than causing an error.
keep <- c(
  "^HYRIV_ID$",
  "^dis_m3_",     # natural discharge: annual mean (pyr), monthly min (pmn), monthly max (pmx)
  "^sgr_dk_",     # stream gradient, decimetres per km
  "^ele_mt_",     # elevation
  "^slp_dg_",     # terrain slope
  "^run_mm_",     # land surface runoff
  "^inu_pc_",     # inundation extent
  "^riv_tc_",     # river area / volume
  "^ari_ix_",     # aridity index (catchment and upstream)
  "^for_pc_",     # forest cover
  "^lit_cl_"      # lithology
)
cols <- unique(unlist(lapply(keep, grep, x = names(ra), value = TRUE)))
ra <- ra[, cols]

write.csv(ra, "riveratlas_sudan.csv", row.names = FALSE)
cat("Wrote", nrow(ra), "reaches and", ncol(ra), "columns to riveratlas_sudan.csv\n")
print(cols)
