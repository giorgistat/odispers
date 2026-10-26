###############################################################################
# Posterior mode of the breeding field u of the tail-down model fitted
# without habitat covariates (directed_taildown_joint.rds; breeding
# log w = theta log Q + u) against the covariates at the river nodes, to
# decide how discharge and the habitat covariates should enter E.
# Nodes: those within 40 km of a village (the only ones the data inform).
# Each panel: 2-d histogram of the nodes, weighted means of u in 20 bins of
# the covariate (weights: kernel weight summed over villages, i.e. how much
# the data see the node) and a weighted GAM smooth. For discharge, a second
# panel shows the total log breeding intensity theta log Q + u.
###############################################################################
library(odispers)
library(Matrix)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
net <- load_result("inputs_full_network.rds"); nodes <- net$nodes
data(sudan_river_covariates)
rc <- as.data.frame(sudan_river_covariates)[match(nodes$HYRIV_ID, sudan_river_covariates$HYRIV_ID), ]
fit <- readRDS(file.path(out_dir, "directed_taildown_joint.rds"))
S <- fit$setup; used <- S$used
g <- directed_graph(river_graph(nodes, vertex_spacing = 5), "taildown")
u <- as.vector(g$A[used, , drop = FALSE] %*% fit$u)
theta <- fit$coefficients[["theta.log(DIS_AV_CMS)"]]
wt <- colSums(S$K)                                  # kernel weight of each node summed over villages
lq <- log(nodes$DIS_AV_CMS[used])
covs <- list(
  "log discharge (m3/s): u" = list(x = lq, y = u),
  "log discharge (m3/s): theta log Q + u" = list(x = lq, y = theta * lq + u),
  "log(1 + slope, dm/km)" = list(x = log1p(rc$sgr_dk_rav[used]), y = u),
  "log(1 + upstream inundation %)" = list(x = log1p(rc$inu_pc_ult[used]), y = u),
  "hard rock (0 = sediments, 1 = basement/volcanic)" =
    list(x = as.numeric(rc$lit_cl_cmj[used] %in% c(2, 4, 7, 8, 9, 10, 12, 13)), y = u),
  "upstream forest cover % (not selected)" = list(x = rc$for_pc_use[used], y = u))
pts <- NULL; bins <- NULL; sm <- NULL; tab <- NULL
for (nm in names(covs)) {
  x <- covs[[nm]]$x; y <- covs[[nm]]$y
  pts <- rbind(pts, data.frame(panel = nm, x = x, y = y, w = wt))
  binary <- length(unique(x)) <= 2
  b <- if (binary) factor(x) else cut(x, unique(quantile(x, seq(0, 1, 0.05))), include.lowest = TRUE)
  for (lv in levels(b)) { s <- b == lv; if (!any(s)) next
    ww <- wt[s] / sum(wt[s]); mu <- sum(ww * y[s]); se <- sqrt(sum(ww * (y[s] - mu)^2) * sum(ww^2))
    bins <- rbind(bins, data.frame(panel = nm, x = sum(ww * x[s]), y = mu, lo = mu - 1.96 * se, hi = mu + 1.96 * se)) }
  if (!binary) {
    gm <- mgcv::gam(y ~ s(x, k = 10), weights = wt / mean(wt), data = data.frame(x = x, y = y), method = "REML")
    xs <- seq(quantile(x, 0.005), quantile(x, 0.995), length.out = 200)
    sm <- rbind(sm, data.frame(panel = nm, x = xs, y = predict(gm, data.frame(x = xs))))
  }
  tab <- rbind(tab, data.frame(panel = nm, weighted_cor = cov.wt(cbind(x, y), wt = wt / sum(wt), cor = TRUE)$cor[1, 2],
                               spearman = cor(x, y, method = "spearman")))
}
pts$panel <- factor(pts$panel, names(covs)); bins$panel <- factor(bins$panel, names(covs))
sm$panel <- factor(sm$panel, names(covs))
print(tab, digits = 2, row.names = FALSE)
p <- ggplot(pts, aes(x, y)) + geom_bin2d(bins = 60) +
  scale_fill_gradient(low = "grey90", high = "grey20", trans = "log10", name = "nodes") +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_line(data = sm, colour = "#b2182b", linewidth = 0.9) +
  geom_pointrange(data = bins, aes(ymin = lo, ymax = hi), colour = "#2166ac", size = 0.2) +
  facet_wrap(~ panel, scales = "free", ncol = 3) + theme_minimal(base_size = 9) +
  labs(x = "covariate at the river node", y = "posterior mode of u (or theta log Q + u)",
       caption = "Nodes within 40 km of a village. Blue: weighted means of u in 20 bins (weights = kernel weight summed over villages). Red: weighted GAM smooth.")
ggsave(file.path(fig_dir, "u_vs_covariates.png"), p, width = 12, height = 7.5, dpi = 180, bg = "white")
saveRDS(list(bins = bins, smooth = sm, cor = tab), file.path(out_dir, "u_vs_covariates.rds"))
