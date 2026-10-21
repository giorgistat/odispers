###############################################################################
# Observed against predicted prevalence for the main model.
###############################################################################

library(odispers); library(sf); library(ggplot2)
source("00_settings.R")
out <- readRDS(file.path(out_dir, "spline_fits.rds")); f <- out[[as.character(phi_main)]]$fit
v <- st_drop_geometry(readRDS(file.path(out_dir, "inputs.rds"))$villages)
v$fitted <- f$fitted; v$obs <- v$y / v$n
v$cls <- cut(v$AI, c(0, 0.05, 0.2, 0.5, Inf), labels = c("Hyper-arid (AI < 0.05)", "Arid (0.05-0.20)", "Semi-arid (0.20-0.50)", "Sub-humid / humid (AI > 0.50)"))
cat("correlation obs vs fitted:", round(cor(v$obs, v$fitted), 3), "\n")
p <- ggplot(v, aes(fitted, obs, colour = cls)) +
  geom_abline(linetype = 2) + geom_point(alpha = 0.55, size = 1.6) +
  scale_colour_manual(values = c("#a6611a", "#dfc27d", "#7fbf7b", "#1b7837"), name = "Aridity class") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted prevalence", y = "Observed prevalence",
       title = "Observed against predicted prevalence",
       subtitle = sprintf("Spline in log aridity (change point at AI = %.3f), range fixed at %d km", exp(out[[as.character(phi_main)]]$cp), phi_main)) +
  theme_minimal(base_size = 11) + theme(legend.position = "right")
ggsave(file.path(fig_dir, "spline_obs_vs_pred.png"), p, width = 8, height = 6, dpi = 180)
