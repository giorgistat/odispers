###############################################################################
# Variogram of the predicted village effect Z in the zero-inflated model with
# the tail-down breeding field (joint fit of 23_directed_profiles.R), among
# the villages with at least one positive, where Z is informed by the data,
# with a 95% envelope of spatial independence from 9999 permutations, a
# global test and a one-sided test of positive correlation up to 100 km (as
# in 18_variogram_z_zil.R).
###############################################################################
library(odispers)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
fit <- readRDS(file.path(out_dir, "directed_taildown_joint.rds"))
n_perm <- 9999
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
pos <- fit$setup$y > 0
z <- fit$z[pos]; xy <- fit$setup$xy[pos, ]
cat(sprintf("variance of predicted Z: villages with positives %.4f, without %.5f\n",
            stats::var(z), stats::var(fit$z[!pos])))
vp <- variogram_perm(z, xy, breaks, n_perm)
cat(sprintf("global permutation test p = %.4f; one-sided test of positive correlation up to 100 km p = %.4f\n",
            vp$p_value, vp$p_positive))
print(round(vp$table, 3), row.names = FALSE)
tab <- vp$table
p <- ggplot(tab, aes(distance)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey85") +
  geom_line(aes(y = perm_mean), colour = "grey50", linetype = 2) +
  geom_line(aes(y = semivariance), colour = "#1b7837") +
  geom_point(aes(y = semivariance), colour = "#1b7837", size = 1.8) +
  labs(x = "Distance between villages (km)", y = "Semivariance of the predicted village effect",
       title = "Village effect Z, tail-down model, villages with at least one positive (n = 350)",
       subtitle = sprintf("Global test p = %.3f; one-sided test of positive correlation p = %.3f; %d permutations",
                          vp$p_value, vp$p_positive, n_perm)) +
  theme_minimal(base_size = 10)
ggsave(file.path(fig_dir, "variogram_z_taildown.png"), p, width = 7, height = 4.6, dpi = 200, bg = "white")
saveRDS(vp, file.path(out_dir, "variogram_z_taildown.rds"))
