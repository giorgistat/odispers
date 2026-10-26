###############################################################################
# Village-level scatter plots of the RiverATLAS covariates against the
# outcome, for all villages, with the kernel range fixed at 40 km: on the x
# axis the village average of each covariate over the rivers within reach
# (kernel-weighted, theta = 0; from 31_covariates_eda.R, cov_eda.rds), on the
# y axis (1) the observed prevalence y/n and (2) the empirical logit
# log{(y + 0.5)/(n - y + 0.5)}, with a loess smooth and its 95% band.
###############################################################################
library(odispers)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
e <- readRDS(file.path(out_dir, "cov_eda.rds"))
d <- e$villages; X <- e$xbar$main[e$ok, ]
labels <- c(log_slope = "log(1 + slope, dm/km)  [+]", steep_share = "share of river with slope >= 2 m/km  [+]",
            velocity = "log velocity proxy  [+]", permanence = "flow permanence (min/mean)  [+]",
            inundation = "log(1 + inundation %)  [-]", hard_rock = "share of river on hard rock  [+]",
            forest = "upstream forest cover %  [+]", lakes = "log(1 + lake %)  [?]",
            regulation = "log(1 + regulation %)  [?]", elevation = "elevation (100 m)  [none]",
            runoff = "log(1 + runoff, mm/yr)  [none]", log_discharge = "log discharge (in the model)")
long <- do.call(rbind, lapply(names(labels), function(k)
  data.frame(covariate = labels[[k]], x = X[, k], n = d$n, y = d$y, prevalence = d$y / d$n,
             elogit = log((d$y + 0.5) / (d$n - d$y + 0.5)), positive = d$y > 0)))
long$covariate <- factor(long$covariate, labels)
plot_y <- function(yvar, ylab) {
  ggplot(long, aes(x, .data[[yvar]])) +
    geom_point(aes(shape = positive), size = 0.7, alpha = 0.45, colour = "grey25") +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16), labels = c("no positives", "positives"), name = NULL) +
    geom_smooth(method = "loess", formula = y ~ x, span = 0.75, colour = "#1b7837", fill = "#1b7837", alpha = 0.2,
                linewidth = 0.9) +
    facet_wrap(~ covariate, scales = "free_x", ncol = 4) + theme_minimal(base_size = 9) +
    theme(legend.position = "bottom") +
    labs(x = "village average of the covariate over the rivers within 40 km", y = ylab)
}
ggsave(file.path(fig_dir, "cov_scatter_prevalence.png"), plot_y("prevalence", "observed prevalence (y / n)"),
       width = 11, height = 8.5, dpi = 180, bg = "white")
ggsave(file.path(fig_dir, "cov_scatter_elogit.png"),
       plot_y("elogit", "empirical logit, log{(y + 0.5) / (n - y + 0.5)}"),
       width = 11, height = 8.5, dpi = 180, bg = "white")
