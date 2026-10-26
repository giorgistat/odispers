###############################################################################
# Figures and tables for the tail-down model with habitat covariates, from
# the output folder of inst/scripts/hpc/taildown_covariates.R:
#   - mean PIT distribution functions of the fitted models
#   - expected prevalence maps (M0 and the model with covariates)
#   - breeding intensity along the rivers implied by the covariates
#   - LaTeX rows of the estimates with 95% intervals
# Usage: Rscript 35_covariate_fit_figures.R <results folder>
###############################################################################
library(odispers)
library(sf)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
res_dir <- commandArgs(trailingOnly = TRUE)[1]
data(sudan_borders); data(sudan_rivers); data(sudan_river_covariates)
models <- c("M0", "M1", "M2", "M3")
labels <- c(M0 = "M0: discharge only", M1 = "M1: + slope, inundation", M2 = "M2: + slope, inundation, hard rock", M3 = "M3: slope, inundation (no discharge)")
R <- list(); for (m in models) { f <- file.path(res_dir, sprintf("fit_%s.rds", m))
  if (file.exists(f)) { r <- readRDS(f); if (!is.null(r$grid)) R[[m]] <- r } }   # completed models only
models <- names(R)

# PIT distribution functions
pit <- do.call(rbind, lapply(models, function(m) if (!is.null(R[[m]]$cv_scores))
  cbind(model = labels[[m]], R[[m]]$cv_scores$pit)))
tests <- do.call(rbind, lapply(models, function(m) if (!is.null(R[[m]]$cv_scores))
  data.frame(model = labels[[m]], dev = R[[m]]$cv_scores$table$pit_max_dev, p = R[[m]]$cv_scores$table$pit_p)))
if (!is.null(pit)) {
  pit$model <- factor(pit$model, labels[models]); tests$model <- factor(tests$model, labels[models])
  p <- ggplot(pit, aes(u)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey85") +
    geom_abline(linetype = 2) + geom_line(aes(y = Fbar), colour = "#1b7837", linewidth = 0.8) +
    geom_text(data = tests, aes(x = 0.03, y = 0.95, label = sprintf("max deviation %.3f\np = %.3f", dev, p)),
              hjust = 0, vjust = 1, size = 3) +
    facet_wrap(~ model, nrow = 1) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(x = "u", y = "Mean PIT distribution function") + theme_minimal(base_size = 10)
  ggsave(file.path(fig_dir, "cov_fit_pit.png"), p, width = 10, height = 3.9, dpi = 200, bg = "white")
}

# prevalence maps
pm <- do.call(rbind, lapply(models, function(m) cbind(model = labels[[m]], R[[m]]$grid)))
pm$model <- factor(pm$model, labels[models])
mp <- ggplot() + geom_tile(data = subset(pm, zero_exposure), aes(lon, lat), fill = "grey80", width = 0.2, height = 0.2) +
  geom_tile(data = subset(pm, !zero_exposure), aes(lon, lat, fill = prevalence), width = 0.2, height = 0.2) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Expected\nprevalence", limits = c(0, NA)) +
  geom_sf(data = sudan_borders, fill = NA, colour = "grey25", linewidth = 0.3) + facet_wrap(~ model, nrow = 1) +
  coord_sf(expand = FALSE) + theme_minimal(base_size = 9) + labs(x = NULL, y = NULL)
ggsave(file.path(fig_dir, "cov_fit_prevalence_maps.png"), mp, width = 13, height = 5, dpi = 200, bg = "white")
w <- split(pm$prevalence, pm$model)
if (length(w) > 1) for (k in 2:length(w))
  cat(sprintf("prevalence on the grid, %s vs %s: correlation %.3f, max |difference| %.3f\n", names(w)[1], names(w)[k],
              cor(w[[1]], w[[k]]), max(abs(w[[1]] - w[[k]]))))

# breeding multiplier from the covariates, exp(beta' h), along the rivers, for the last model
mlast <- tail(models, 1); cf <- R[[mlast]]$fit$coefficients
riv <- sudan_rivers[sudan_rivers$DIS_AV_CMS >= 0.1, ]
rc <- as.data.frame(sudan_river_covariates)[match(riv$HYRIV_ID, sudan_river_covariates$HYRIV_ID), ]
h <- cbind(log_slope = log1p(rc$sgr_dk_rav), inundation = pmin(log1p(rc$inu_pc_ult), 1.5),
           hard_rock = as.numeric(rc$lit_cl_cmj %in% c(2, 4, 7, 8, 9, 10, 12, 13)))
b <- cf[paste0("theta.", colnames(h))]; b <- b[!is.na(b)]
riv$eta <- as.vector(h[, sub("theta.", "", names(b)), drop = FALSE] %*% b)
riv$eta <- riv$eta - median(riv$eta)
lim <- quantile(abs(riv$eta), 0.99)
mb <- ggplot() + geom_sf(data = sudan_borders, fill = "grey97", colour = "grey40", linewidth = 0.3) +
  geom_sf(data = riv[order(abs(riv$eta)), ], aes(colour = pmax(pmin(eta, lim), -lim)), linewidth = 0.25) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b", midpoint = 0, limits = c(-lim, lim),
                         name = "log breeding\nmultiplier\n(vs median)") +
  coord_sf(expand = FALSE) + theme_minimal(base_size = 9) + labs(x = NULL, y = NULL)
ggsave(file.path(fig_dir, "cov_fit_breeding_multiplier.png"), mb, width = 7, height = 6.5, dpi = 200, bg = "white")

# table rows with intervals
fmt <- function(e, l, u) sprintf("%.2f (%.2f, %.2f)", e, l, u)
for (m in models) if (!is.null(R[[m]]$intervals)) { cat("\n", m, "\n"); print(R[[m]]$intervals$table, digits = 3) }
