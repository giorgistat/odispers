###############################################################################
# Non-randomised PIT (Czado, Gneiting and Held, 2009) shown as a distribution
# function: the mean PIT distribution Fbar(u) = mean_i F_i(u | y_i) over the
# villages with at least one positive, against u, with the identity line
# (calibration) and a pointwise 95% band from counts simulated from the
# predictive distributions, for the random-fold cross-validation of the
# three directed fields (cv_directed_joint.rds from 23_directed_profiles.R).
# Also reports the largest deviation from the identity and its simulated
# p-value.
###############################################################################
library(odispers)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
preds <- readRDS(file.path(out_dir, "cv_directed_joint.rds"))
labels <- c(taildown = "Tail-down", tailup = "Tail-up (K1, discharge weights)",
            "tailup+taildown" = "Tail-up + tail-down")
ug <- seq(0, 1, by = 0.01)
cond_cdf <- function(pmf) { p0 <- pmf[1]; c(0, cumsum(pmf[-1]) / (1 - p0)) }
pit_fun <- function(Fm1, F0, u) ifelse(u <= Fm1, 0, ifelse(u >= F0, 1, (u - Fm1) / pmax(F0 - Fm1, 1e-12)))
mean_pit_cdf <- function(ys, cdfs) rowMeans(sapply(seq_along(ys), function(i) pit_fun(cdfs[[i]][ys[i]], cdfs[[i]][ys[i] + 1], ug)))
set.seed(2026)
out <- NULL; tests <- NULL
for (t in names(labels)) {
  pr <- preds[[t]]; if (is.null(pr)) next
  pos <- Filter(function(z) z$y > 0, pr); cdfs <- lapply(pos, function(z) cond_cdf(z$pmf))
  Fbar <- mean_pit_cdf(sapply(pos, `[[`, "y"), cdfs)
  sims <- replicate(1000, { ys <- sapply(pos, function(z) { qq <- z$pmf[-1] / (1 - z$pmf[1]); sample.int(length(qq), 1, prob = qq) })
                            mean_pit_cdf(ys, cdfs) })
  dev_obs <- max(abs(Fbar - ug)); dev_sim <- apply(abs(sims - ug), 2, max)
  tests <- rbind(tests, data.frame(model = labels[[t]], max_deviation = round(dev_obs, 4),
                                   p_value = round((1 + sum(dev_sim >= dev_obs)) / (1 + length(dev_sim)), 3)))
  out <- rbind(out, data.frame(model = labels[[t]], u = ug, Fbar = Fbar,
                               lower = apply(sims, 1, stats::quantile, 0.025), upper = apply(sims, 1, stats::quantile, 0.975)))
}
print(tests, row.names = FALSE)
out$model <- factor(out$model, labels[labels %in% out$model])
tests$model <- factor(tests$model, levels(out$model))
p <- ggplot(out, aes(u)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey85") +
  geom_abline(linetype = 2) + geom_line(aes(y = Fbar), colour = "#1b7837", linewidth = 0.8) +
  geom_text(data = tests, aes(x = 0.03, y = 0.95, label = sprintf("max deviation %.3f\np = %.3f", max_deviation, p_value)),
            hjust = 0, vjust = 1, size = 3) +
  facet_wrap(~ model, nrow = 1) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "u", y = "Mean PIT distribution function") + theme_minimal(base_size = 10)
ggsave(file.path(fig_dir, "directed_pit_cdf.png"), p, width = 10, height = 3.9, dpi = 200, bg = "white")
saveRDS(list(curves = out, tests = tests), file.path(out_dir, "directed_pit_cdf.rds"))
