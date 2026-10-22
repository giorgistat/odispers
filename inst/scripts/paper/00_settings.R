###############################################################################
# Settings shared by the scripts that reproduce the analysis in
# "A river-network exposure model for onchocerciasis prevalence"
# (range fixed at 40 km, breeding intensity quadratic in log discharge,
# aridity in the intensity of transmission as a linear spline in log AI).
#
# All inputs come with the package and are loaded with data() in the scripts
# that use them: data(remo), data(sudan_rivers), data(sudan_aridity) and
# data(sudan_borders).
###############################################################################

out_dir <- "paper_results"     # intermediate .rds files
fig_dir <- "paper_figures"     # figures used in the paper
dir.create(out_dir, showWarnings = FALSE)
dir.create(fig_dir, showWarnings = FALSE)

min_flow  <- 0.1   # m3/s; reaches with less flow are dropped
spacing   <- 2     # km between quadrature nodes
max_range <- 60    # km; pairs are computed up to this distance (>= every range used)
phi_main  <- 40    # km; dispersal range of the main analysis
phi_sens  <- 53    # km; sensitivity analysis (calibration to forest biting-rate decline)
grid_km   <- 20    # resolution of the prediction grid

breeding <- ~ log(DIS_AV_CMS) + I(log(DIS_AV_CMS)^2)

aridity_classes <- function(ai)
  cut(ai, c(0, 0.05, 0.2, 0.5, Inf),
      labels = c("Hyper-arid", "Arid", "Semi-arid", "Sub-humid/humid"))

# Results of earlier scripts, read from out_dir. The inputs are built by
# 02_prepare_inputs.R if they are missing; other results must be created by
# running the script that writes them.
load_result <- function(name) {
  f <- file.path(out_dir, name)
  if (!file.exists(f) && name == "inputs.rds") {
    message("paper_results/inputs.rds not found: running 02_prepare_inputs.R first")
    source(system.file("scripts/paper/02_prepare_inputs.R", package = "odispers"), local = new.env())
  }
  if (!file.exists(f)) {
    by <- c(spline_fits.rds = "03_fit_spline.R", u_fits.rds = "07_fit_u_field.R",
          u_fits_linear.rds = "08_fit_u_linear.R")[name]
    stop(sprintf("%s not found in %s: run %s first", name, normalizePath(out_dir),
                 if (is.na(by)) "the script that creates it" else by), call. = FALSE)
  }
  readRDS(f)
}
