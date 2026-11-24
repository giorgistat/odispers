###############################################################################
# Inspect the outputs of m3_migration.R: prints every table with the checks
# that tell whether the fits make sense, and draws the maps (one PNG per map
# and a 2 x 3 panel, written with the cairo device so that it also works on a
# cluster without X11).
#
# Usage:
#   Rscript inspect_m3_migration.R [out_dir]        (default m3_migration_results)
# or, in R, set out_dir and source() it.
# Only the CSV files written by m3_migration.R are read (no refitting).
###############################################################################
suppressPackageStartupMessages({ library(odispers); library(sf); library(ggplot2) })
# PNG files without an X11 display (compute nodes): use the cairo device
if (capabilities("cairo")) options(bitmapType = "cairo")
if (!exists("out_dir")) {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args) >= 1) args[1] else "m3_migration_results"
}
stopifnot(dir.exists(out_dir))
fig_dir <- file.path(out_dir, "inspect"); dir.create(fig_dir, showWarnings = FALSE)
f <- function(x) file.path(out_dir, x)
rd <- function(x) if (file.exists(f(x))) read.csv(f(x)) else { cat("  !! missing:", x, "\n"); NULL }
hdr <- function(x) cat("\n", strrep("=", 78), "\n ", x, "\n", strrep("=", 78), "\n", sep = "")
flag <- function(ok, msg) cat(if (isTRUE(ok)) "  [ok]    " else "  [CHECK] ", msg, "\n", sep = "")

hdr("Files in the results folder")
print(data.frame(file = list.files(out_dir), size_kb = round(file.size(list.files(out_dir, full.names = TRUE)) / 1024, 1)),
      row.names = FALSE)

# ---- 1. parameter estimates ---------------------------------------------------------------
hdr("1. Parameter estimates: M3 refit and M3 + migration (95% Wald intervals)")
e0 <- rd("estimates_M3_refit.csv"); e1 <- rd("intervals_M3_migration.csv")
if (!is.null(e1)) {
  tab <- e1
  if (!is.null(e0)) tab$M3_refit <- e0$estimate[match(tab$parameter, e0$parameter)]
  num <- vapply(tab, is.numeric, TRUE); tab[num] <- lapply(tab[num], signif, 4)
  print(tab, row.names = FALSE)
  flag(!anyNA(e1$se), "standard errors available for every parameter of the migration model")
  if (!is.null(e0)) {
    common <- intersect(e0$parameter, e1$parameter)
    rel <- abs(e1$estimate[match(common, e1$parameter)] - e0$estimate[match(common, e0$parameter)]) /
      pmax(abs(e0$estimate[match(common, e0$parameter)]), 1e-3)
    cat("  largest relative changes from M3 (migration vs M3 refit):\n")
    print(head(sort(setNames(round(rel, 3), common), decreasing = TRUE), 5))
  }
}

# ---- 2. model comparison -----------------------------------------------------------------
hdr("2. Model comparison: log-likelihood, AIC, likelihood ratio test for w = 0, CV scores")
cmp <- rd("comparison_migration.csv")
if (!is.null(cmp)) {
  print(t(data.frame(cmp[, -1], row.names = cmp$model)), digits = 5)
  lr <- 2 * (cmp$logLik[2] - cmp$logLik[1])
  p <- 0.5 * pchisq(max(lr, 0), 1, lower.tail = FALSE)
  cat(sprintf("\n  w hat = %.4g; likelihood ratio = %.3f; p = %.3g (50:50 mixture of chi2_0 and chi2_1)\n",
              cmp$w[2], lr, p))
  flag(lr > -0.2, sprintf("migration model at least as good as M3 (difference in logLik %.3f; M3 is the case w = 0)",
                          lr / 2))
  if ("sum_log_score" %in% names(cmp)) {
    cat(sprintf("  CV log score, migration - M3: %.2f (block-bootstrap SE in the .out log)\n",
                diff(cmp$sum_log_score)))
    for (k in seq_len(nrow(cmp))) {
      flag(cmp$pit_p[k] > 0.05, sprintf("%s: PIT max deviation %.3f, p = %.3f", cmp$model[k], cmp$pit_max_dev[k], cmp$pit_p[k]))
      flag(abs(cmp$zeros_observed[k] - cmp$zeros_expected[k]) < 2 * sqrt(cmp$zeros_expected[k]),
           sprintf("%s: zeros observed %d, expected %.1f (humid: %d vs %.1f)", cmp$model[k], cmp$zeros_observed[k],
                   cmp$zeros_expected[k], cmp$zeros_humid_observed[k], cmp$zeros_humid_expected[k]))
    }
  } else cat("  (no CV scores in the table: the cv stage has not been run)\n")
}

# ---- 3. profile likelihood for w ---------------------------------------------------------
hdr("3. Profile log-likelihood for w")
prof <- rd("profile_w.csv")
if (!is.null(prof)) {
  prof <- prof[order(prof$w), ]
  print(prof, digits = 6, row.names = FALSE)
  inside <- prof$w[prof$drop <= qchisq(0.95, 1) / 2]
  cat(sprintf("  approximate 95%% profile interval for w (on the grid): %.3g to %.3g\n", min(inside), max(inside)))
  # a profile that rises again away from the maximum points to a fit that did not converge
  i_max <- which.max(prof$logLik); ll <- prof$logLik
  mono <- all(diff(ll[seq_len(i_max)]) >= -0.05) && all(diff(ll[i_max:length(ll)]) <= 0.05)
  flag(mono, "profile rises to its maximum and then falls (no bumps from fits that stopped early)")
  g <- ggplot(prof, aes(w, -drop)) + geom_line(colour = "grey50") + geom_point(size = 2) +
    geom_hline(yintercept = -qchisq(0.95, 1) / 2, linetype = 2, colour = "#b2182b") +
    scale_x_continuous(trans = "sqrt") +
    labs(x = "w (square-root scale)", y = "profile log-likelihood - maximum",
         title = "Profile likelihood for the migration share w", subtitle = "dashed line: 95% cut-off") +
    theme_minimal(base_size = 10)
  ggsave(file.path(fig_dir, "profile_w.png"), g, width = 6, height = 4, dpi = 200, bg = "white")
} else cat("  (run the profile stage for a grid of w and then the summary stage)\n")

# ---- 4. share of exposure due to migration at the villages -------------------------------
hdr("4. Share of each village's exposure due to migration, s_i")
sv <- rd("migration_share_villages.csv")
if (!is.null(sv)) {
  print(summary(sv$s_mig))
  cat(sprintf("  villages with s_i > 0.1: %d; > 0.5: %d (of %d)\n", sum(sv$s_mig > 0.1), sum(sv$s_mig > 0.5), nrow(sv)))
  sv$prev <- sv$y / sv$n
  cat("  mean observed prevalence by class of s_i:\n")
  cl <- cut(sv$s_mig, unique(c(-Inf, quantile(sv$s_mig, c(0.25, 0.5, 0.75)), Inf)))
  print(round(tapply(sv$prev, cl, mean), 3))
}

# ---- 5. grid predictions -----------------------------------------------------------------
hdr("5. Grid predictions (20 km)")
gp <- rd("grid_predictions.csv")
if (!is.null(gp)) {
  cat(sprintf("  %d grid cells\n", nrow(gp)))
  cols <- setdiff(names(gp), c("lon", "lat"))
  print(t(sapply(gp[cols], function(x) c(n_NA = sum(!is.finite(x)), min = min(x, na.rm = TRUE),
                                         median = median(x, na.rm = TRUE), max = max(x, na.rm = TRUE)))), digits = 3)
  flag(all(is.finite(gp$M3_prevalence)) && all(is.finite(gp$mig_prevalence)), "no missing predictions")
  d <- gp$mig_prevalence - gp$M3_prevalence
  cat(sprintf("  prevalence difference (migration - M3): mean %.4f, 1%%-99%% range %.4f to %.4f\n",
              mean(d, na.rm = TRUE), quantile(d, 0.01, na.rm = TRUE), quantile(d, 0.99, na.rm = TRUE)))
  cat(sprintf("  correlation of the two prevalence surfaces: %.4f\n", cor(gp$M3_prevalence, gp$mig_prevalence, use = "complete")))

  # ---- maps ----
  data(sudan_borders)
  border <- st_transform(sudan_borders, 4326)
  data(sudan_rivers)
  riv <- st_transform(sudan_rivers[sudan_rivers$DIS_AV_CMS >= 10, ], 4326)   # main rivers, for orientation
  bb <- st_bbox(border); zoom <- coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], expand = FALSE)
  # tile size from the grid itself (cells are 20 km in a projected system, so not exactly regular in degrees)
  tw <- median(diff(sort(unique(round(gp$lon, 3))))); th <- median(diff(sort(unique(round(gp$lat, 3)))))
  tw <- if (is.finite(tw) && tw > 0) min(tw * 1.05, 0.25) else 0.19; th <- if (is.finite(th) && th > 0) min(th * 1.05, 0.25) else 0.19
  base <- function(var, title, fill, pts = NULL) {
    g <- ggplot() + geom_tile(data = gp, aes(lon, lat, fill = .data[[var]]), width = tw, height = th) +
      geom_sf(data = riv, colour = "steelblue", linewidth = 0.2, alpha = 0.6) +
      geom_sf(data = border, fill = NA, colour = "grey20", linewidth = 0.3) + fill +
      labs(title = title, x = NULL, y = NULL, fill = NULL) + theme_minimal(base_size = 9) + zoom
    if (!is.null(pts)) g <- g + pts
    g
  }
  prev_lim <- c(0, max(c(gp$M3_prevalence, gp$mig_prevalence), na.rm = TRUE))
  seq_fill <- scale_fill_viridis_c(option = "magma", direction = -1, limits = prev_lim, na.value = "grey85")
  dmax <- max(abs(d), 1e-6, na.rm = TRUE)
  gp$diff <- d
  maps <- list(
    M3_prevalence  = base("M3_prevalence", "M3: predicted prevalence", seq_fill),
    mig_prevalence = base("mig_prevalence", "M3 + migration: predicted prevalence", seq_fill),
    difference     = base("diff", "Difference (migration - M3)",
                          scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                                               limits = c(-dmax, dmax))),
    share          = base("mig_s_mig", "Share of exposure due to migration",
                          scale_fill_viridis_c(limits = c(0, max(gp$mig_s_mig, 1e-6, na.rm = TRUE)), na.value = "grey85")),
    exposure       = base("M3_E", "M3: exposure E (log10)",
                          scale_fill_viridis_c(trans = "log10", na.value = "grey85")))
  if (!is.null(sv)) {
    maps$observed <- ggplot() + geom_sf(data = border, fill = "grey97", colour = "grey20", linewidth = 0.3) +
      geom_sf(data = riv, colour = "steelblue", linewidth = 0.2, alpha = 0.6) +
      geom_point(data = sv, aes(long, lat, colour = prev), size = 0.9) +
      scale_colour_viridis_c(option = "magma", direction = -1, limits = prev_lim, oob = scales::squish) +
      labs(title = "Observed prevalence (surveyed villages)", x = NULL, y = NULL, colour = NULL) + theme_minimal(base_size = 9) + zoom
  }
  for (nm in names(maps))
    ggsave(file.path(fig_dir, paste0("map_", nm, ".png")), maps[[nm]], width = 6, height = 5, dpi = 200, bg = "white")
  # all maps on one page
  png_type <- if (capabilities("cairo")) "cairo" else NULL
  fn <- file.path(fig_dir, "maps_all.png")
  if (is.null(png_type)) grDevices::png(fn, width = 15, height = 9, units = "in", res = 150) else
    grDevices::png(fn, width = 15, height = 9, units = "in", res = 150, type = png_type)
  grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 3)))
  for (k in seq_along(maps))
    print(maps[[k]], vp = grid::viewport(layout.pos.row = (k - 1) %/% 3 + 1, layout.pos.col = (k - 1) %% 3 + 1))
  grDevices::dev.off()
}

# PIT figure from the summary stage
if (file.exists(f("pit_migration.png"))) file.copy(f("pit_migration.png"), file.path(fig_dir, "pit_migration.png"), overwrite = TRUE)

hdr("Figures")
fl <- list.files(fig_dir, full.names = TRUE)
print(data.frame(file = fl, size_kb = round(file.size(fl) / 1024, 1)), row.names = FALSE)
