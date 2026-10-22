###############################################################################
# Zero-inflated model with regressions on river exposure and aridity in both
# parts:
#   logit pi(x)    = a0 + a_E log E(x) + a1 log AI(x)        (zone of transmission)
#   log lambda(x)  = b0 + gamma log E(x) + b1 log AI(x) + Z  (force of infection)
# with E from the river network (breeding linear in log discharge, field u,
# range 800 km), on the full packaged river network so that predictions cover
# the whole of Sudan and South Sudan. Fit, maps, and cross-validated
# calibration with the random folds of 12_cv_calibration.R, compared with the
# models without zero inflation and with the earlier zero-inflated models.
###############################################################################

library(odispers)
library(sf)
library(Matrix)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

inp <- load_result("inputs.rds")
v <- inp$villages
data(sudan_rivers); data(sudan_borders)

# ---- river network covering the whole prediction area --------------------------
f_full <- file.path(out_dir, "inputs_full_network.rds")
if (file.exists(f_full)) {
  net <- readRDS(f_full)
} else {
  rivers <- sudan_rivers[sudan_rivers$DIS_AV_CMS >= min_flow, ]
  nodes <- river_nodes(rivers, spacing = spacing)
  net <- list(nodes = nodes, pairs = river_pairs(v, nodes, phi_main),
              grid_pairs = river_pairs(inp$grid, nodes, phi_main))
  saveRDS(net, f_full)
}
nodes <- net$nodes; print(nodes)
graph <- river_graph(nodes, vertex_spacing = 5); print(graph)

# ---- fit ----------------------------------------------------------------------
fit <- odisp_fit_zil(y ~ log(AI), ~ log(DIS_AV_CMS), v, nodes, graph, net$pairs, phi = phi_main,
                     range_u = 800, start = c(log(4), log(0.5), 1, 0.2))
print(fit)
cls <- aridity_classes(v$AI)
cat("\nMean probability of a zone of transmission at the villages, by aridity class:\n")
print(round(tapply(fit$pi, cls[fit$setup$has], mean), 3))

# ---- maps -----------------------------------------------------------------------
grid <- inp$grid
pg <- predict(fit, grid, nodes, graph, pairs = net$grid_pairs)
grid <- cbind(grid, pg)
xy <- st_coordinates(grid); grid$lon <- xy[, 1]; grid$lat <- xy[, 2]
gcls <- aridity_classes(grid$AI)
d_v <- apply(st_distance(st_transform(grid, 32636), st_transform(v, 32636)), 1, min) / 1000
cat("\nExpected prevalence on the grid (median) by aridity class and distance from surveys:\n")
print(aggregate(cbind(pi, prevalence) ~ gcls + I(d_v > 100), data = data.frame(st_drop_geometry(grid), gcls, d_v),
                FUN = function(x) round(stats::median(x), 3)))
cat("Arid or hyper-arid cells with expected prevalence > 5%:",
    sum(grid$prevalence > 0.05 & gcls %in% c("Hyper-arid", "Arid")), "of", sum(gcls %in% c("Hyper-arid", "Arid")), "\n")
base <- list(geom_sf(data = sudan_borders, fill = NA, colour = "grey25", linewidth = 0.3),
             coord_sf(expand = FALSE), theme_minimal(base_size = 10), labs(x = NULL, y = NULL))
m_pi <- ggplot(grid) + geom_tile(aes(lon, lat, fill = pi), width = 0.2, height = 0.2) +
  scale_fill_viridis_c(name = "P(transmission\nzone)", limits = c(0, 1)) + base +
  labs(title = "Probability of lying in a zone of local transmission")
ggsave(file.path(fig_dir, "zil_pi_map.png"), m_pi, width = 6.5, height = 7, dpi = 200, bg = "white")
m_prev <- ggplot(grid) + geom_tile(aes(lon, lat, fill = prevalence), width = 0.2, height = 0.2) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Expected\nprevalence", limits = c(0, NA)) +
  base + labs(title = "Expected nodule prevalence")
ggsave(file.path(fig_dir, "zil_prevalence_map.png"), m_prev, width = 6.5, height = 7, dpi = 200, bg = "white")
vd <- st_drop_geometry(v)[fit$setup$has, ]; vd$pi <- fit$pi
m_v <- ggplot() + geom_sf(data = sudan_borders, fill = "grey97", colour = "grey40", linewidth = 0.3) +
  geom_point(data = vd, aes(long, lat, colour = pi, shape = y > 0), size = 1.4) +
  scale_colour_viridis_c(name = "P(transmission\nzone)", limits = c(0, 1)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "positives", `FALSE` = "no positives"),
                     name = NULL) +
  coord_sf(expand = FALSE) + theme_minimal(base_size = 10) + labs(x = NULL, y = NULL,
  title = "Probability that a village lies in a zone of local transmission")
ggsave(file.path(fig_dir, "zil_pi_villages.png"), m_v, width = 6.5, height = 7, dpi = 200, bg = "white")
saveRDS(list(fit = fit, grid = grid), file.path(out_dir, "zil_fit.rds"))

# ---- cross-validation (random folds of 12_cv_calibration.R) ----------------------
cv12 <- load_result("cv_calibration.rds")
folds <- cv12$folds[["random folds"]]
S_all <- fit$setup; fld <- fit$fld
subset_setup <- function(S, idx) {
  S$y <- S$y[idx]; S$n <- S$n[idx]; S$X <- S$X[idx, , drop = FALSE]; S$K <- S$K[idx, , drop = FALSE]
  S$xy <- S$xy[idx, , drop = FALSE]; S$has <- S$has[idx]; S$lchoose <- sum(lchoose(S$n, S$y)); S$lai <- S$lai[idx]
  S
}
gh <- local({ m <- 20; J <- matrix(0, m, m)
  for (k in 1:(m - 1)) J[k, k + 1] <- J[k + 1, k] <- sqrt(k / 2)
  e <- eigen(J, symmetric = TRUE); list(x = sqrt(2) * e$values, w = e$vectors[1, ]^2) })
g2 <- expand.grid(i = seq_along(gh$x), j = seq_along(gh$x))
g2x1 <- gh$x[g2$i]; g2x2 <- gh$x[g2$j]; g2w <- gh$w[g2$i] * gh$w[g2$j]
pmf_dE <- function(md, me, Sig, n) {
  L <- t(chol(Sig + diag(1e-10, 2)))
  d <- md + L[1, 1] * g2x1; e <- me + L[2, 1] * g2x1 + L[2, 2] * g2x2
  pi_ <- stats::plogis(d); p <- -expm1(-exp(pmin(e, 30)))
  pm <- sapply(0:n, function(y) sum(g2w * pi_ * stats::dbinom(y, n, p)))
  pm[1] <- pm[1] + sum(g2w * (1 - pi_))
  pm
}
p <- ncol(S_all$X); q <- ncol(S_all$B)
cv_fold <- function(train, test) {
  S_tr <- subset_setup(S_all, train); S_te <- subset_setup(S_all, test)
  x0 <- list(beta = fit$coefficients[seq_len(p)], c = fit$coefficients[p + seq_len(p)],
             theta = fit$coefficients[2 * p + seq_len(q)], u = fit$u, z = rep(0, length(train)))
  ft <- odispers:::odisp_zil_laplace(fld, S_tr, x0)
  ev <- odispers:::odisp_zil_eval(ft$x, S_tr, fld)
  lk <- odispers:::odisp_zil_lik(ev$d, ev$e, S_tr)
  ops <- odispers:::odisp_zil_ops(ft$x, ev, lk, S_tr, fld)
  xt <- ft$x; xt$z <- rep(0, length(test))
  et <- odispers:::odisp_zil_eval(xt, S_te, fld)
  Jt <- as.matrix(et$V %*% S_te$B); Ju <- as.matrix(et$V %*% S_te$A); nz <- length(train)
  lapply(seq_along(test), function(k) {
    jd <- c(rep(0, p), S_te$X[k, ], fld$aE * Jt[k, ], fld$aE * Ju[k, ], rep(0, nz))
    je <- c(S_te$X[k, ], rep(0, p), fld$gamma * Jt[k, ], fld$gamma * Ju[k, ], rep(0, nz))
    hd <- ops$HGNinv(jd); he <- ops$HGNinv(je)
    Sig <- matrix(c(sum(jd * hd), sum(jd * he), sum(je * hd), sum(je * he)), 2)
    Sig <- (Sig + t(Sig)) / 2; Sig[2, 2] <- Sig[2, 2] + fld$sig_z^2
    list(village = test[k], y = S_te$y[k], n = S_te$n[k], pmf = pmf_dE(et$d[k], et$e[k], Sig, S_te$n[k]))
  })
}
cat("cross-validation, random folds: ")
out <- unlist(lapply(sort(unique(folds)), function(k) { cat(k, ""); cv_fold(which(folds != k), which(folds == k)) }),
              recursive = FALSE)
cat("\n")
preds <- list("zero-inflated, log E and log AI in both parts" = out[order(sapply(out, `[[`, "village"))],
              "no zero inflation" = cv12$preds[["random folds | u + village effect"]])
fz <- file.path(out_dir, "cv_zero_inflated.rds")
if (file.exists(fz)) preds[["zero-inflated, log E in pi"]] <- readRDS(fz)$preds[["random folds | zero-inflated u + village effect"]]

cond_cdf <- function(pmf) { p0 <- pmf[1]; c(0, cumsum(pmf[-1]) / (1 - p0)) }
pit_fun <- function(Fm1, F0, u) ifelse(u <= Fm1, 0, ifelse(u >= F0, 1, (u - Fm1) / pmax(F0 - Fm1, 1e-12)))
mean_pit_hist <- function(ys, cdfs, J = 10) {
  ub <- (0:J) / J
  Fbar <- rowMeans(sapply(seq_along(ys), function(i) pit_fun(cdfs[[i]][ys[i]], cdfs[[i]][ys[i] + 1], ub)))
  J * diff(Fbar)
}
auc <- function(pp, o) { r <- rank(pp); n1 <- sum(o); n0 <- sum(!o); (sum(r[o]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
set.seed(2026)
pit_tab <- list(); rel_tab <- list(); scores <- list(); zero_tab <- list()
for (nm in names(preds)) {
  pr <- preds[[nm]]
  pos <- Filter(function(z) z$y > 0, pr); cdfs <- lapply(pos, function(z) cond_cdf(z$pmf))
  h <- mean_pit_hist(sapply(pos, `[[`, "y"), cdfs)
  sims <- replicate(200, { ys <- sapply(pos, function(z) { qq <- z$pmf[-1] / (1 - z$pmf[1]); sample.int(length(qq), 1, prob = qq) })
                           mean_pit_hist(ys, cdfs) })
  pit_tab[[nm]] <- data.frame(model = nm, mid = (seq_along(h) - 0.5) / length(h), density = h,
                              lower = apply(sims, 1, stats::quantile, 0.025), upper = apply(sims, 1, stats::quantile, 0.975))
  p1 <- sapply(pr, function(z) 1 - z$pmf[1]); o1 <- sapply(pr, `[[`, "y") > 0
  scores[[nm]] <- data.frame(model = nm, AUC = round(auc(p1, o1), 3), Brier = round(mean((p1 - o1)^2), 4),
                             zeros_expected = round(sum(1 - p1), 1),
                             mean_log_score = round(mean(sapply(pr, function(z) log(max(z$pmf[z$y + 1], 1e-300)))), 3))
  br <- unique(stats::quantile(p1, seq(0, 1, 0.1))); b <- cut(p1, br, include.lowest = TRUE)
  rel <- data.frame(model = nm, pred = tapply(p1, b, mean), obs = tapply(o1, b, mean), n = as.vector(table(b)))
  rel$lower <- stats::qbeta(0.025, rel$obs * rel$n + 0.5, (1 - rel$obs) * rel$n + 0.5)
  rel$upper <- stats::qbeta(0.975, rel$obs * rel$n + 0.5, (1 - rel$obs) * rel$n + 0.5)
  rel_tab[[nm]] <- rel
  zero_tab[[nm]] <- data.frame(model = nm, class = levels(cls), observed = as.vector(tapply(!o1, cls, sum)),
                               expected = round(as.vector(tapply(1 - p1, cls, sum)), 1))
}
scores <- do.call(rbind, scores); rownames(scores) <- NULL
cat("\nRandom-fold cross-validation (550 villages without positives observed):\n"); print(scores, row.names = FALSE)
cat("\nVillages without positives by aridity class:\n"); print(do.call(rbind, zero_tab), row.names = FALSE)
pit_tab <- do.call(rbind, pit_tab); rel_tab <- do.call(rbind, rel_tab)
show <- c("no zero inflation", "zero-inflated, log E and log AI in both parts")
pt <- subset(pit_tab, model %in% show); pt$model <- factor(pt$model, show)
p_pit <- ggplot(pt, aes(mid, density)) + geom_col(width = 0.1, fill = "#7fbf7b", colour = "white") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.04, colour = "grey30") +
  geom_hline(yintercept = 1, linetype = 2) + facet_wrap(~ model) +
  labs(x = "Probability integral transform", y = "Density") + theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "zil_cv_pit_random.png"), p_pit, width = 8, height = 3.6, dpi = 200, bg = "white")
rt <- subset(rel_tab, model %in% show)
p_rel <- ggplot(rt, aes(pred, obs, colour = model)) + geom_abline(linetype = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, alpha = 0.6) + geom_point(size = 2) + geom_line() +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = c("zero-inflated, log E and log AI in both parts" = "#1b7837", "no zero inflation" = "grey50"),
                      labels = c("zero-inflated, log E and log AI in both parts" = "zero-inflated", "no zero inflation" = "no zero inflation"),
                      name = NULL) +
  labs(x = "Predicted probability of at least one positive", y = "Observed proportion") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "zil_cv_reliability_random.png"), p_rel, width = 5.5, height = 5.5, dpi = 200, bg = "white")
saveRDS(list(preds = preds, scores = scores, zeros = do.call(rbind, zero_tab), pit = pit_tab, reliability = rel_tab),
        file.path(out_dir, "cv_zil.rds"))
