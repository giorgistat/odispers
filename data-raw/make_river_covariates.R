###############################################################################
# Builds data/sudan_river_covariates.rda: the RiverATLAS (HydroATLAS v1.0)
# attributes needed for the breeding-habitat part of the river models, for the
# same reaches as `sudan_rivers` (matched by HYRIV_ID). Run from the package
# root on a machine that has the RiverATLAS Africa shapefile; the result is
# small (about 1 MB) and is the only thing that has to be committed.
###############################################################################

library(sf)

# <- edit: the Africa shapefile of RiverATLAS v1.0 (unzipped)
src <- "/home/giorgi/Downloads/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp/RiverATLAS_v10_af.shp"

source_cols <- c(
  "sgr_dk_rav",                            # channel slope, dm per km (water-speed proxy)
  "dis_m3_pyr", "dis_m3_pmn", "dis_m3_pmx", # discharge: annual mean, min month, max month (m3/s)
  "run_mm_cyr",                            # runoff, mm per year
  "inu_pc_umn", "inu_pc_umx", "inu_pc_ult", # inundated area upstream, % (min, max, long-term max)
  "lka_pc_use",                            # lake area upstream, %
  "dor_pc_pva",                            # degree of regulation (dams), %
  "for_pc_use",                            # forest cover upstream, %
  "lit_cl_cmj",                            # lithology class (majority in the reach catchment)
  "ele_mt_cav",                            # elevation of the reach catchment, m
  "riv_tc_usu"                             # river volume upstream, thousand m3
)

# Read only the attribute table (no geometry) for the columns we need
layer <- tools::file_path_sans_ext(basename(src))
q <- sprintf("SELECT HYRIV_ID, %s FROM %s", paste(source_cols, collapse = ", "), layer)
ra <- suppressWarnings(st_read(src, query = q, quiet = TRUE))
ra <- as.data.frame(ra)
if ("geometry" %in% names(ra)) ra$geometry <- NULL

# HydroATLAS codes missing values as -9999
for (v in source_cols) ra[[v]][ra[[v]] == -9999] <- NA

# Keep exactly the reaches of the package river network, in the same order
load("data/sudan_rivers.rda")
sudan_river_covariates <- ra[match(sudan_rivers$HYRIV_ID, ra$HYRIV_ID), ]
stopifnot(!anyNA(sudan_river_covariates$HYRIV_ID),
          identical(sudan_river_covariates$HYRIV_ID, sudan_rivers$HYRIV_ID))
rownames(sudan_river_covariates) <- NULL

# Whole numbers are stored as integers to keep the file small
for (v in setdiff(names(sudan_river_covariates), c("dis_m3_pyr", "dis_m3_pmn", "dis_m3_pmx", "riv_tc_usu")))
  sudan_river_covariates[[v]] <- as.integer(sudan_river_covariates[[v]])

print(summary(sudan_river_covariates))
cat("missing values per column:\n"); print(colSums(is.na(sudan_river_covariates)))
save(sudan_river_covariates, file = "data/sudan_river_covariates.rda", compress = "xz")
cat("size:", round(file.size("data/sudan_river_covariates.rda") / 1e6, 2), "MB\n")
