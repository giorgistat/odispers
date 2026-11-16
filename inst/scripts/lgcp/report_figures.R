# Rebuild the LGCP report figures and numerical tables from its bundled data.
# Run from the repository root: Rscript inst/scripts/lgcp/report_figures.R
# Forecasts are archived from 7 November; field diagnostics use 10 November refits.
library(odispers)
library(sf)
report <- "inst/reports/lgcp"
cv <- read.csv(file.path(report, "data/heldout_comparison.csv"))
full <- read.csv(file.path(report, "data/full_data_comparison.csv"))
predictions <- read.csv(file.path(report, "data/heldout_predictions.csv"))
villages <- read.csv(file.path(report, "data/villages.csv"))
widths <- read.csv(file.path(report, "data/width_traces.csv"))
field_name <- c(wm = "Whittle--Mat\\'ern", taildown = "Tail-down")
field_colour <- c(wm = "#3C6F9D", taildown = "#B15C31")

# Check that each displayed score can be reconstructed from village forecasts.
for (row in seq_len(nrow(cv))) {
  d <- subset(predictions, field == cv$field[row] & range == cv$range[row])
  stopifnot(nrow(d) == nrow(villages), !anyDuplicated(d$village),
            isTRUE(all.equal(mean(d$brier), cv$brier[row])),
            isTRUE(all.equal(mean(d$crps), cv$crps[row])),
            isTRUE(all.equal(sum(d$log_score), cv$log_score[row])))
}
tex_rows <- function(rows, filename) {
  writeLines(c(paste0(rows, " \\\\"), "\\bottomrule"), file.path(report, filename))
}
tex_rows(vapply(seq_len(nrow(cv)), function(i) {
  d <- cv[i, ]
  sprintf("%s & %d & %.4f & %.4f & %.1f & %.1f", field_name[d$field],
          d$range, d$brier, d$crps, d$log_score, d$expected_nonzero)
}, ""), "heldout_rows.tex")
tex_rows(vapply(seq_len(nrow(cv)), function(i) {
  d <- cv[i, ]
  sprintf("%s & %d & %.2f & %.1f & %d/900 & %d/900", field_name[d$field],
          d$range, d$max_rhat, d$min_ess, d$prediction_rhat_flagged, d$prediction_ess_flagged)
}, ""), "diagnostic_rows.tex")
tex_rows(vapply(seq_len(nrow(full)), function(i) {
  d <- full[i, ]
  sprintf("%s & %d & %.1f & %.1f--%.1f & %.1f & %.2f", field_name[d$field],
          d$range, d$width_median, d$width_lo, d$width_hi, d$expected_nonzero,
          100*d$positive_fraction)
}, ""), "full_data_rows.tex")

# Compare the saved held-out forecasts on the same score scales.
png(file.path(report, "figures/heldout_scores.png"), 1800, 640, res = 180)
par(mfrow = c(1, 3), mar = c(4.5, 5.6, 2.5, 1), mgp = c(2.8, 1, 0), las = 1)
for (metric in c("brier", "crps", "mean_log_score")) {
  labels <- c(brier = "Brier score", crps = "Prevalence CRPS", mean_log_score = "Mean count log score")
  plot(range(cv$range), range(cv[[metric]]), type = "n", xlab = "River-field range (km)",
       ylab = "", xaxt = "n", main = if (metric == "mean_log_score") "Higher is better" else "Lower is better")
  mtext(labels[metric], side = 2, line = 4.2, las = 0)
  axis(1, c(200, 400, 800))
  for (field in names(field_colour)) {
    d <- cv[cv$field == field, ]
    lines(d$range, d[[metric]], col = field_colour[field], type = "b", pch = 19, lwd = 2)
  }
  if (metric == "brier") legend("topright", c("Whittle-Matern", "Tail-down"),
                               col = field_colour, pch = 19, lwd = 2, bty = "n", cex = .8)
}
dev.off()

# The separate field-monitoring experiment stores full-draw diagnostics and
# ACFs alongside display-thinned traces, so thinning never changes the estimates.
source("inst/scripts/lgcp/field_mixing_figure.R")
field_plot <- readRDS(file.path(report, "data/field_plot_data.rds"))
draw_lgcp_field_mixing(field_plot, file.path(report, "figures/field_mixing.png"))
field_diagnostics <- read.csv(file.path(report, "data/field_diagnostics.csv"))
field_variables <- c("field_mean_all", "field_mean_footprint", "field_mean_north",
                     "field_mean_south", "field_north_minus_south", field_plot$variables[3])
field_labels <- c("Whole-network mean", "Survey-footprint mean", "Northern mean",
                  "Southern mean", "North--south contrast", "Largest local $\\widehat R$")
tex_rows(vapply(seq_along(field_variables), function(i) {
  d <- field_diagnostics[field_diagnostics$variable == field_variables[i], ]
  sprintf("%s & %.3f & %.1f & %.1f & %.3f & %.3f", field_labels[i],
          d$rhat, d$ess_bulk, d$ess_tail, d$lag1_mean, d$lag100_mean)
}, ""), "field_mixing_rows.tex")

# Show how each full-data chain explores the kernel width; keep all regimes.
png(file.path(report, "figures/width_traces.png"), 1800, 1600, res = 180)
par(mfrow = c(3, 2), mar = c(4.5, 4, 2.4, 1), mgp = c(2.7, .8, 0), las = 1)
chain_colours <- c("#1B6CA8", "#D87929", "#2D8755", "#9A57A5")
for (range_km in c(200, 400, 800)) for (field in c("wm", "taildown")) {
  d <- widths[widths$field == field & widths$range == range_km, ]
  title <- paste(if (field == "wm") "Whittle-Matern" else "Tail-down", range_km, "km")
  plot(c(2000, 8000), c(1, 240), type = "n", log = "y", xlab = "Sweep",
       ylab = "Kernel width (km)", main = title)
  abline(h = 60, col = "grey70", lty = 2)
  for (chain in 1:4) {
    z <- d[d$chain == chain, ]
    lines(z$iteration, z$width_km, col = chain_colours[chain], lwd = .8)
  }
  if (range_km == 200 && field == "wm") legend("topleft", paste("Chain", 1:4),
    col = chain_colours, lwd = 2, bty = "n", cex = .7, ncol = 2)
}
dev.off()

# Survey-location maps: the right panel combines five separate held-out fits.
# It is not a single jointly simulated landscape or an unsurveyed-location map.
data(sudan_borders, package = "odispers")
borders <- st_transform(sudan_borders, 4326)
d <- subset(predictions, field == "taildown" & range == 800)
d <- d[match(seq_len(nrow(villages)), d$village), ]
prevalence_colours <- hcl.colors(100, "YlOrRd", rev = TRUE)
prevalence_colour <- function(p) prevalence_colours[1L + floor(pmin(p, 1)*99)]
png(file.path(report, "figures/heldout_map.png"), 1600, 1050, res = 180)
par(mfrow = c(1, 2), mar = c(3.5, 3.5, 3, 1))
for (panel in 1:2) {
  risk <- if (panel == 1) villages$y/villages$n else d$mean_prevalence
  plot(st_geometry(borders), col = "#F4F4F2", border = "grey60", axes = TRUE,
       xlim = range(villages$long) + c(-.5, .5), ylim = range(villages$lat) + c(-.5, .5),
       main = c("Observed fraction positive", "Held-out mean prevalence")[panel])
  points(villages$long, villages$lat, pch = 21, cex = .65, lwd = .25,
         bg = prevalence_colour(risk), col = "grey35")
  legend("topright", legend = c("0", "0.1", "0.3", "0.6", "1"),
         pt.bg = prevalence_colour(c(0, .1, .3, .6, 1)), pch = 21, bty = "n", cex = .8,
         title = "Prevalence")
}
dev.off()
