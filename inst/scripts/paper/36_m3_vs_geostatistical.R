###############################################################################
# Comparison of the river-network model M3 with the Euclidean geostatistical
# models G1 (one Gaussian process in both parts) and G2 (two independent
# processes), from the outputs of
#   inst/scripts/hpc/taildown_covariates.R   (fit_M3.rds)
#   inst/scripts/hpc/geostat_gp_models.R     (gp_G1.rds, gp_G2.rds)
# Produces: comparison table (csv), paired differences in cross-validated log
# score with standard errors (independent villages and 1 degree block
# bootstrap), mean PIT distribution functions and expected prevalence maps.
# Usage: Rscript 36_m3_vs_geostatistical.R <folder with fit_M3.rds> <folder with gp_G*.rds> [hessian_M3.rds]
###############################################################################
library(odispers)
library(sf)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
a <- commandArgs(trailingOnly = TRUE)
m3_dir <- a[1]; gp_dir <- a[2]; hess_file <- if (length(a) >= 3) a[3] else NA
data(remo); data(sudan_borders)
v <- subset(remo, country == "SUDAN")

R3 <- readRDS(file.path(m3_dir, "fit_M3.rds"))
G <- lapply(c(G1 = "G1", G2 = "G2"), function(m) readRDS(file.path(gp_dir, sprintf("gp_%s.rds", m))))
iv3 <- R3$fit$intervals
if (!is.na(hess_file) && file.exists(hess_file)) iv3 <- readRDS(hess_file)$intervals
labels <- c(M3 = "M3: river network (tail-down, slope, inundation)", G1 = "G1: one Euclidean GP", G2 = "G2: two Euclidean GPs")

# ---- comparison table -----------------------------------------------------------------------
row <- function(m, logLik, k, cv) cbind(data.frame(model = m, logLik = logLik, n_par = k, AIC = -2 * logLik + 2 * k), cv)
cmp <- rbind(row("M3", R3$fit$logLik, 11, R3$cv_scores$table),
             row("G1", G$G1$fit$logLik, G$G1$fit$n_par, G$G1$cv_scores$table),
             row("G2", G$G2$fit$logLik, G$G2$fit$n_par, G$G2$cv_scores$table))
print(cmp, digits = 5, row.names = FALSE)
write.csv(cmp, file.path(out_dir, "m3_vs_gp_comparison.csv"), row.names = FALSE)

# ---- paired differences in log score ----------------------------------------------------------
ls_of <- function(pr) sapply(pr, function(z) log(max(z$pmf[z$y + 1], 1e-300)))
vil <- sapply(R3$cv, `[[`, "village")
stopifnot(identical(vil, sapply(G$G1$cv, `[[`, "village")), identical(vil, sapply(G$G2$cv, `[[`, "village")))
ls <- cbind(M3 = ls_of(R3$cv), G1 = ls_of(G$G1$cv), G2 = ls_of(G$G2$cv))
cell <- paste(floor(v$long[vil]), floor(v$lat[vil]))
set.seed(2026)
cells <- unique(cell); idx <- split(seq_along(cell), cell)
boot <- function(d) replicate(2000, { s <- sample(cells, replace = TRUE); sum(d[unlist(idx[s], use.names = FALSE)]) })
pd <- do.call(rbind, lapply(c("G1", "G2"), function(m) {
  d <- ls[, "M3"] - ls[, m]; b <- boot(d)
  data.frame(comparison = sprintf("M3 - %s", m), difference = sum(d), se_iid = sd(d) * sqrt(length(d)),
             se_block = sd(b), lower_block = quantile(b, 0.025), upper_block = quantile(b, 0.975),
             share_villages_M3_better = mean(d > 0), row.names = NULL)
}))
print(pd, digits = 4)
write.csv(pd, file.path(out_dir, "m3_vs_gp_logscore_differences.csv"), row.names = FALSE)

# ---- mean PIT -----------------------------------------------------------------------------------
pit <- rbind(cbind(model = labels[["M3"]], R3$cv_scores$pit), cbind(model = labels[["G1"]], G$G1$cv_scores$pit),
             cbind(model = labels[["G2"]], G$G2$cv_scores$pit))
pit$model <- factor(pit$model, labels)
tests <- data.frame(model = factor(labels, labels), dev = cmp$pit_max_dev, p = cmp$pit_p)
p <- ggplot(pit, aes(u)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey85") +
  geom_abline(linetype = 2) + geom_line(aes(y = Fbar), colour = "#1b7837", linewidth = 0.8) +
  geom_text(data = tests, aes(x = 0.03, y = 0.95, label = sprintf("max deviation %.3f\np = %.3f", dev, p)),
            hjust = 0, vjust = 1, size = 3) +
  facet_wrap(~ model, nrow = 1, labeller = label_wrap_gen(30)) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "u", y = "Mean PIT distribution function") + theme_minimal(base_size = 10)
ggsave(file.path(fig_dir, "m3_vs_gp_pit.png"), p, width = 10, height = 3.9, dpi = 200, bg = "white")

# ---- prevalence maps ------------------------------------------------------------------------------
g3 <- R3$grid; g3$model <- labels[["M3"]]
gg <- lapply(c("G1", "G2"), function(m) { g <- G[[m]]$grid; g$zero_exposure <- FALSE; g$model <- labels[[m]]; g })
cols <- c("lon", "lat", "prevalence", "zero_exposure", "model")
pm <- rbind(g3[cols], gg[[1]][cols], gg[[2]][cols]); pm$model <- factor(pm$model, labels)
mp <- ggplot() + geom_tile(data = subset(pm, zero_exposure), aes(lon, lat), fill = "grey80", width = 0.2, height = 0.2) +
  geom_tile(data = subset(pm, !zero_exposure), aes(lon, lat, fill = prevalence), width = 0.2, height = 0.2) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Expected\nprevalence", limits = c(0, NA)) +
  geom_sf(data = sudan_borders, fill = NA, colour = "grey25", linewidth = 0.3) +
  facet_wrap(~ model, nrow = 1, labeller = label_wrap_gen(30)) +
  coord_sf(expand = FALSE) + theme_minimal(base_size = 9) + labs(x = NULL, y = NULL)
ggsave(file.path(fig_dir, "m3_vs_gp_prevalence_maps.png"), mp, width = 13, height = 5, dpi = 200, bg = "white")
key <- function(g) paste(round(g$lon, 4), round(g$lat, 4))
for (m in c("G1", "G2")) {
  j <- match(key(g3), key(G[[m]]$grid))
  cat(sprintf("prevalence on the grid, M3 vs %s: correlation %.3f; cells with zero exposure in M3 (%d): mean %s prevalence %.4f, max %.4f\n",
              m, cor(g3$prevalence, G[[m]]$grid$prevalence[j], use = "complete.obs"), sum(g3$zero_exposure), m,
              mean(G[[m]]$grid$prevalence[j][g3$zero_exposure], na.rm = TRUE), max(G[[m]]$grid$prevalence[j][g3$zero_exposure], na.rm = TRUE)))
}

# ---- estimates -----------------------------------------------------------------------------------
fmt <- function(e, l, u) sprintf("%.2f (%.2f, %.2f)", e, l, u)
for (m in c("M3", "G1", "G2")) {
  iv <- if (m == "M3") iv3 else G[[m]]$fit$intervals
  cat("\n", m, "\n"); cat(sprintf("%-20s %s\n", iv$parameter, fmt(iv$estimate, iv$lower, iv$upper)), sep = "")
}
