###############################################################################
# Is the village effect Z spatially correlated once the breeding field u is in
# the model? Empirical variogram of the predicted village effects (posterior
# modes) with a 95% envelope under spatial independence from random
# permutations of the effects among the village locations, and a global
# permutation test. Done for the linear model of 08_fit_u_linear.R with the
# village effect only and with u.
###############################################################################

library(odispers)
library(sf)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

res <- load_result("u_fits_linear.rds")
inp <- load_result("inputs.rds")
crs_km <- attr(inp$nodes, "crs")

n_perm <- 999
breaks <- c(0, 10, 25, 50, 75, 100, 150, 200, 250, 300)   # km
set.seed(2026)

# Empirical variogram of z at locations xy (km), for the distance bins in
# `breaks`, with a permutation envelope and a global test. The test statistic
# is the sum over bins of squared standardised differences between the
# variogram and the mean of the permutation variograms.
variogram_perm <- function(z, xy, breaks, n_perm) {
  D <- as.matrix(stats::dist(xy))
  ut <- which(upper.tri(D))
  ii <- row(D)[ut]; jj <- col(D)[ut]; d <- D[ut]
  bin <- cut(d, breaks, include.lowest = TRUE)
  keep <- !is.na(bin); ii <- ii[keep]; jj <- jj[keep]; bin <- droplevels(bin[keep])
  npairs <- tabulate(bin, nlevels(bin))
  vg <- function(zz) as.vector(tapply((zz[ii] - zz[jj])^2 / 2, bin, mean))
  obs <- vg(z)
  perm <- replicate(n_perm, vg(sample(z)))
  m <- rowMeans(perm); s <- apply(perm, 1, stats::sd)
  stat <- function(g) sum(((g - m) / s)^2)
  t_obs <- stat(obs); t_perm <- apply(perm, 2, stat)
  mid <- (breaks[-1] + breaks[-length(breaks)]) / 2
  list(table = data.frame(bin = levels(bin), distance = mid[seq_len(nlevels(bin))], pairs = npairs,
                          semivariance = obs, lower = apply(perm, 1, stats::quantile, 0.025),
                          upper = apply(perm, 1, stats::quantile, 0.975), perm_mean = m),
       p_value = (1 + sum(t_perm >= t_obs)) / (n_perm + 1))
}

models <- list("Village effect only" = res$fit_z, "With the breeding field u" = res$best)
out <- lapply(names(models), function(nm) {
  f <- models[[nm]]
  v <- res$villages[f$setup$has, ]
  xy <- st_coordinates(st_transform(st_geometry(v), crs_km))
  vp <- variogram_perm(f$z, xy, breaks, n_perm)
  cat(sprintf("\n%s: variance of predicted Z %.3f; global permutation test p = %.3f\n",
              nm, stats::var(f$z), vp$p_value))
  print(within(vp$table, { semivariance <- round(semivariance, 3); lower <- round(lower, 3)
                           upper <- round(upper, 3); perm_mean <- round(perm_mean, 3) }), row.names = FALSE)
  cbind(model = sprintf("%s (permutation test p = %.3f)", nm, vp$p_value), vp$table)
})
tab <- do.call(rbind, out)
tab$model <- factor(tab$model, levels = unique(tab$model))

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
ggsave(file.path(fig_dir, "variogram_z_envelope.png"), p, width = 10, height = 4.5, dpi = 200, bg = "white")
saveRDS(tab, file.path(out_dir, "variogram_z.rds"))
