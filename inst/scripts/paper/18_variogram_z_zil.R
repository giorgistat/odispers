###############################################################################
# Residual spatial correlation of the village effect Z in the zero-inflated
# model of 17_fit_zil.R (log E and log AI in both parts), compared with the
# model without zero inflation of 08_fit_u_linear.R. Empirical variogram of
# the predicted village effects (posterior modes), with a 95% envelope under
# spatial independence from 999 random permutations of the effects among the
# villages and a global permutation test (as in 09_variogram_z.R).
###############################################################################

library(odispers)
library(sf)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

zil <- load_result("zil_fit.rds")$fit
f08 <- file.path(out_dir, "u_fits_linear.rds")
n_perm <- 999
breaks <- c(0, 10, 25, 50, 75, 100, 150, 200, 250, 300)   # km
set.seed(2026)

variogram_perm <- function(z, xy, breaks, n_perm) {
  D <- as.matrix(stats::dist(xy))
  ut <- which(upper.tri(D))
  ii <- row(D)[ut]; jj <- col(D)[ut]; d <- D[ut]
  bin <- cut(d, breaks, include.lowest = TRUE)
  keep <- !is.na(bin); ii <- ii[keep]; jj <- jj[keep]; bin <- droplevels(bin[keep])
  vg <- function(zz) as.vector(tapply((zz[ii] - zz[jj])^2 / 2, bin, mean))
  obs <- vg(z)
  perm <- replicate(n_perm, vg(sample(z)))
  m <- rowMeans(perm); s <- apply(perm, 1, stats::sd)
  stat <- function(g) sum(((g - m) / s)^2)
  t_obs <- stat(obs); t_perm <- apply(perm, 2, stat)
  # one-sided test against positive correlation: semivariances below the
  # permutation mean at distances up to 100 km
  short <- breaks[-1][seq_along(obs)] <= 100
  stat1 <- function(g) sum(((m - g) / s)[short])
  t1_obs <- stat1(obs); t1_perm <- apply(perm, 2, stat1)
  mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
  list(table = data.frame(distance = mid[seq_len(nlevels(bin))], pairs = tabulate(bin, nlevels(bin)),
                          semivariance = obs, lower = apply(perm, 1, stats::quantile, 0.025),
                          upper = apply(perm, 1, stats::quantile, 0.975), perm_mean = m),
       p_value = (1 + sum(t_perm >= t_obs)) / (n_perm + 1),
       p_positive = (1 + sum(t1_perm >= t1_obs)) / (n_perm + 1))
}

# In the zero-inflated model Z acts only inside zones of transmission, so for
# villages without positives the data carry almost no information on it and
# its prediction is shrunk to about 0; as these villages cluster in space,
# the variogram over all villages is pulled down at short distances. The
# variogram is therefore also shown for the villages with positives only.
pos <- zil$setup$y > 0
models <- list("Zero-inflated, all villages" = list(z = zil$z, xy = zil$setup$xy),
               "Zero-inflated, villages with positives" = list(z = zil$z[pos], xy = zil$setup$xy[pos, ]))
if (file.exists(f08)) {
  b <- readRDS(f08)$best
  xy08 <- b$setup$xy
  if (is.null(xy08)) {            # fits made before the coordinates were stored
    inp <- load_result("inputs.rds")
    xy08 <- st_coordinates(st_transform(st_geometry(inp$villages), attr(inp$nodes, "crs")))[b$setup$has, 1:2]
  }
  models <- c(list("No zero inflation (u + village effect)" = list(z = b$z, xy = xy08)), models)
}
out <- lapply(names(models), function(nm) {
  vp <- variogram_perm(models[[nm]]$z, models[[nm]]$xy, breaks, n_perm)
  cat(sprintf("\n%s: variance of predicted Z %.4f; global permutation test p = %.3f; one-sided test of positive correlation up to 100 km p = %.3f\n",
              nm, stats::var(models[[nm]]$z), vp$p_value, vp$p_positive))
  print(round(vp$table, 3), row.names = FALSE)
  cbind(model = sprintf("%s\nglobal test p = %.3f; test of positive correlation p = %.3f",
                        nm, vp$p_value, vp$p_positive), vp$table)
})
tab <- do.call(rbind, out); tab$model <- factor(tab$model, levels = unique(tab$model))

p <- ggplot(tab, aes(distance)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey85") +
  geom_line(aes(y = perm_mean), colour = "grey50", linetype = 2) +
  geom_line(aes(y = semivariance), colour = "#1b7837") +
  geom_point(aes(y = semivariance), colour = "#1b7837", size = 1.8) +
  facet_wrap(~ model, scales = "free_y") +
  labs(x = "Distance between villages (km)", y = "Semivariance of the predicted village effect",
       title = "Empirical variogram of the village effect Z",
       subtitle = sprintf("Grey band: 95%% envelope under spatial independence (%d permutations); dashed: permutation mean",
                          n_perm)) +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "variogram_z_zil.png"), p, width = 13, height = 4.6, dpi = 200, bg = "white")
saveRDS(tab, file.path(out_dir, "variogram_z_zil.rds"))

# ---- observed against predicted prevalence ---------------------------------------
# in-sample: expected prevalence pi * p at the posterior mode with Z = 0;
# cross-validated: mean of the random-fold predictive distribution (17_fit_zil.R)
inp <- load_result("inputs.rds")
v <- inp$villages[zil$setup$has, ]
e0 <- log(-log1p(-pmin(zil$p, 1 - 1e-12))) - zil$z
pred_in <- zil$pi * (-expm1(-exp(e0)))
cvz <- load_result("cv_zil.rds")$preds[["zero-inflated, log E and log AI in both parts"]]
pred_cv <- rep(NA, nrow(inp$villages))
pred_cv[sapply(cvz, `[[`, "village")] <- sapply(cvz, function(z) sum((seq_along(z$pmf) - 1) * z$pmf) / z$n)
pred_cv <- pred_cv[zil$setup$has]
d <- rbind(data.frame(type = "In-sample (posterior mode, village effect 0)", pred = pred_in,
                      obs = zil$setup$y / zil$setup$n, AI = v$AI),
           data.frame(type = "Cross-validated (10 random folds, predictive mean)", pred = pred_cv,
                      obs = zil$setup$y / zil$setup$n, AI = v$AI))
d$type <- factor(d$type, unique(d$type))
d$class <- aridity_classes(d$AI)
r <- tapply(seq_len(nrow(d)), d$type, function(k) stats::cor(d$pred[k], d$obs[k]))
cat("\nCorrelation of observed and predicted prevalence:\n"); print(round(r, 3))
levels(d$type) <- sprintf("%s\ncorrelation %.2f", levels(d$type), r)
p2 <- ggplot(d, aes(pred, obs, colour = class)) +
  geom_abline(linetype = 2) + geom_point(alpha = 0.55, size = 1.5) +
  scale_colour_manual(values = c("#a6611a", "#dfc27d", "#7fbf7b", "#1b7837"), name = "Aridity class") +
  facet_wrap(~ type) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted prevalence", y = "Observed prevalence",
       title = "Observed against predicted prevalence, zero-inflated model") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "zil_obs_vs_pred.png"), p2, width = 10, height = 5.2, dpi = 200, bg = "white")
