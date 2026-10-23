###############################################################################
# Directed breeding fields in the zero-inflated model of 17_fit_zil.R.
# The symmetric Whittle-Matern field u is replaced by the directed fields of
# Bolin, Simas, Karlsson Strandh and Wallin (2026), see ?directed_graph:
#   tail-down        dependence by river distance, stationary variance;
#   tail-up (K1)     dependence along the flow only, with discharge-weighted
#                    mixing at confluences (weights DIS_AV_CMS, omega = 1);
#   tail-up + tail-down, with the share of variance phi estimated.
# Each model is fitted at the range of u of the symmetric model (800 km) and
# at 200 and 3200 km (the combined model, with five parameters, at 800 km only), and the best fit of each is cross-validated with the
# random folds of 12_cv_calibration.R. Results are compared with the
# symmetric field (zil_fit.rds, cv_zil.rds).
#
# Usage: Rscript 21_fit_directed.R [model ...]  (models: symmetric, taildown,
# tailup, tailup+taildown; default: all four, the symmetric one at the ranges
# other than 800 km, which is taken from zil_fit.rds),
# then once more with argument "compare" to cross-validate and compare.
# Fits are saved after each range, so the script can be restarted.
###############################################################################

library(odispers)
library(sf)
library(Matrix)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

args <- commandArgs(trailingOnly = TRUE)
types <- c("taildown", "tailup", "tailup+taildown")
ranges <- c(200, 800, 3200)
todo <- if (length(args)) args else c("symmetric", types)

inp <- load_result("inputs.rds"); v <- inp$villages
net <- load_result("inputs_full_network.rds")
nodes <- net$nodes
graph <- river_graph(nodes, vertex_spacing = 5)
sym <- load_result("zil_fit.rds")$fit
fit_file <- function(type) file.path(out_dir, sprintf("directed_%s.rds", gsub("\\+", "_", type)))

# ---- fits -------------------------------------------------------------------------
for (type in intersect(todo, c("symmetric", types))) {
  dg <- if (type == "symmetric") graph else directed_graph(graph, type)
  cat(type, ":", format(dg$nv, big.mark = ","), "latent values of u\n")
  fits <- if (file.exists(fit_file(type))) readRDS(fit_file(type)) else list()
  if (type == "symmetric" && is.null(fits[["800"]])) fits[["800"]] <- sym
  # the combined model, with five parameters, is fitted at 800 km only
  for (r in if (type == "tailup+taildown") 800 else c(800, setdiff(ranges, 800))) {
    key <- as.character(r)
    if (!is.null(fits[[key]])) next
    cat(sprintf("\n==== %s, range %g km ====\n", type, r))
    t0 <- Sys.time()
    st <- if (length(fits)) fits[[1]]$par else c(sym$par, if (type == "tailup+taildown") 0)
    f <- odisp_fit_zil(y ~ log(AI), ~ log(DIS_AV_CMS), v, nodes, dg, net$pairs, phi = phi_main,
                       range_u = r, start = st)
    f$time <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    print(f); cat(sprintf("time %.1f min\n", f$time))
    f$setup$graph <- NULL                         # the graph is rebuilt when needed
    fits[[key]] <- f
    saveRDS(fits, fit_file(type))
  }
}
if (!"compare" %in% todo) quit(save = "no")

# ---- cross-validation of the best range of each model -------------------------------
cv12 <- load_result("cv_calibration.rds")
folds <- cv12$folds[["random folds"]]
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
cv_zil <- function(fit, folds) {
  S_all <- fit$setup; fld <- fit$fld
  fld$Rch <- Matrix::Cholesky(fld$R, LDL = FALSE, perm = TRUE)
  p <- ncol(S_all$X); q <- ncol(S_all$B)
  one <- function(train, test) {
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
  out <- unlist(lapply(sort(unique(folds)), function(k) { cat(k, ""); one(which(folds != k), which(folds == k)) }),
                recursive = FALSE)
  cat("\n")
  out[order(sapply(out, `[[`, "village"))]
}

labels <- c(symmetric = "Symmetric Whittle-Matern (current)", taildown = "Tail-down",
            tailup = "Tail-up (K1, discharge weights)", "tailup+taildown" = "Tail-up + tail-down")
best <- list(); prof <- NULL
for (type in c("symmetric", types)) {
  if (!file.exists(fit_file(type))) { if (type == "symmetric") { best$symmetric <- sym
    prof <- data.frame(model = "symmetric", range = 800, logLik = sym$logLik) }; next }
  fits <- readRDS(fit_file(type))
  ll <- sapply(fits, `[[`, "logLik")
  prof <- rbind(prof, data.frame(model = type, range = as.numeric(names(fits)), logLik = ll))
  best[[type]] <- fits[[which.max(ll)]]
}
cat("\nApproximate log-likelihood by range of the field:\n"); print(prof[order(prof$model, prof$range), ], row.names = FALSE)

cv_file <- file.path(out_dir, "cv_directed.rds")
preds <- if (file.exists(cv_file)) readRDS(cv_file)$preds else list()
if (best$symmetric$hyper[["range_u"]] == 800 && is.null(preds[["symmetric 800"]]))
  preds[["symmetric 800"]] <- load_result("cv_zil.rds")$preds[["zero-inflated, log E and log AI in both parts"]]
for (type in names(best)) {
  key <- paste(type, best[[type]]$hyper[["range_u"]])
  if (!is.null(preds[[key]])) next
  cat("cross-validation,", key, ": ")
  b <- best[[type]]
  preds[[key]] <- cv_zil(b, folds)
  saveRDS(list(preds = preds), cv_file)
}
names_cv <- sapply(names(best), function(t) paste(t, best[[t]]$hyper[["range_u"]]))

cond_cdf <- function(pmf) { p0 <- pmf[1]; c(0, cumsum(pmf[-1]) / (1 - p0)) }
pit_fun <- function(Fm1, F0, u) ifelse(u <= Fm1, 0, ifelse(u >= F0, 1, (u - Fm1) / pmax(F0 - Fm1, 1e-12)))
mean_pit_hist <- function(ys, cdfs, J = 10) {
  ub <- (0:J) / J
  Fbar <- rowMeans(sapply(seq_along(ys), function(i) pit_fun(cdfs[[i]][ys[i]], cdfs[[i]][ys[i] + 1], ub)))
  J * diff(Fbar)
}
auc <- function(pp, o) { r <- rank(pp); n1 <- sum(o); n0 <- sum(!o); (sum(r[o]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
set.seed(2026)
cls <- aridity_classes(v$AI[sym$setup$has])
scores <- list(); pit_tab <- list()
for (k in seq_along(names_cv)) {
  nm <- names(best)[k]; pr <- preds[[names_cv[k]]]; b <- best[[nm]]
  p1 <- sapply(pr, function(z) 1 - z$pmf[1]); o1 <- sapply(pr, `[[`, "y") > 0
  pos <- Filter(function(z) z$y > 0, pr); cdfs <- lapply(pos, function(z) cond_cdf(z$pmf))
  h <- mean_pit_hist(sapply(pos, `[[`, "y"), cdfs)
  sims <- replicate(500, { ys <- sapply(pos, function(z) { qq <- z$pmf[-1] / (1 - z$pmf[1]); sample.int(length(qq), 1, prob = qq) })
                           mean_pit_hist(ys, cdfs) })
  pit_tab[[nm]] <- data.frame(model = labels[[nm]], mid = (seq_along(h) - 0.5) / 10, density = h,
                              lower = apply(sims, 1, stats::quantile, 0.025), upper = apply(sims, 1, stats::quantile, 0.975))
  pred_mean <- sapply(pr, function(z) sum((seq_along(z$pmf) - 1) * z$pmf) / z$n)
  obs <- sapply(pr, function(z) z$y / z$n)
  scores[[nm]] <- data.frame(model = labels[[nm]], range_km = b$hyper[["range_u"]], logLik = round(b$logLik, 2),
                             sigma_u = round(b$hyper[["sigma_u"]], 2), sigma_z = round(b$hyper[["sigma_z"]], 3),
                             a_E = round(b$hyper[["a_E"]], 3), gamma = round(b$hyper[["gamma"]], 3),
                             phi = if ("phi" %in% names(b$hyper)) round(b$hyper[["phi"]], 2) else NA,
                             AUC = round(auc(p1, o1), 3), Brier = round(mean((p1 - o1)^2), 4),
                             log_score = round(mean(sapply(pr, function(z) log(max(z$pmf[z$y + 1], 1e-300)))), 3),
                             zeros_expected = round(sum(1 - p1), 1),
                             zeros_humid = round(sum((1 - p1)[cls == "Sub-humid/humid"]), 1),
                             cor_cv = round(stats::cor(pred_mean, obs), 3),
                             pit_bins_out = sum(h < pit_tab[[nm]]$lower | h > pit_tab[[nm]]$upper))
}
scores <- do.call(rbind, scores); rownames(scores) <- NULL
cat("\nComparison (random-fold cross-validation; 550 villages without positives, 14 of them humid):\n")
print(scores, row.names = FALSE)
pit_tab <- do.call(rbind, pit_tab); pit_tab$model <- factor(pit_tab$model, labels[names(best)])
p_pit <- ggplot(pit_tab, aes(mid, density)) + geom_col(width = 0.1, fill = "#7fbf7b", colour = "white") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.04, colour = "grey30") +
  geom_hline(yintercept = 1, linetype = 2) + facet_wrap(~ model, nrow = 1) +
  labs(x = "Non-randomised PIT (villages with at least one positive)", y = "Density") + theme_minimal(base_size = 10)
ggsave(file.path(fig_dir, "directed_pit.png"), p_pit, width = 12, height = 3.4, dpi = 200, bg = "white")
prof$label <- factor(labels[prof$model], labels)
p_prof <- ggplot(prof, aes(range, logLik, colour = label)) + geom_line() + geom_point(size = 2) +
  scale_x_log10() + labs(x = "Range of the breeding field (km, log scale)", y = "Approximate log-likelihood",
                         colour = NULL) + theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "directed_profile.png"), p_prof, width = 7, height = 4.5, dpi = 200, bg = "white")
saveRDS(list(preds = preds, scores = scores, profile = prof, pit = pit_tab), cv_file)

# ---- maps of the breeding field and of expected prevalence ---------------------------
data(sudan_borders)
ll <- st_coordinates(st_transform(st_as_sf(as.data.frame(nodes)[, c("x", "y")], coords = c("x", "y"),
                                           crs = attr(nodes, "crs")), 4326))
sel <- seq(1, nrow(nodes), by = 3)
grid <- inp$grid; xy <- st_coordinates(grid)
umaps <- list(); pmaps <- list()
for (nm in names(best)) {
  b <- best[[nm]]
  g <- if (nm == "symmetric") graph else directed_graph(graph, nm)
  umaps[[nm]] <- data.frame(model = labels[[nm]], lon = ll[sel, 1], lat = ll[sel, 2],
                            u = as.vector(g$A[sel, , drop = FALSE] %*% b$u))
  pg <- predict(b, grid, nodes, g, pairs = net$grid_pairs)
  pmaps[[nm]] <- data.frame(model = labels[[nm]], lon = xy[, 1], lat = xy[, 2], prevalence = pg$prevalence,
                            zero = pg$zero_exposure)
}
um <- do.call(rbind, umaps); um$model <- factor(um$model, labels[names(best)]); um <- um[order(abs(um$u)), ]
lim <- max(abs(quantile(um$u, c(0.005, 0.995))))
bb <- sf::st_bbox(sudan_borders)
m_u <- ggplot() + geom_sf(data = sudan_borders, fill = "grey97", colour = "grey40", linewidth = 0.3) +
  geom_point(data = um, aes(lon, lat, colour = pmax(pmin(u, lim), -lim)), size = 0.08) +
  scale_colour_gradient2(low = "#2166ac", mid = "grey85", high = "#b2182b", midpoint = 0, limits = c(-lim, lim),
                         name = "u (mode)") + facet_wrap(~ model, nrow = 1) +
  coord_sf(xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], expand = FALSE) + theme_minimal(base_size = 9) +
  labs(x = NULL, y = NULL)
ggsave(file.path(fig_dir, "directed_u_maps.png"), m_u, width = 14, height = 4.8, dpi = 200, bg = "white")
pm <- do.call(rbind, pmaps); pm$model <- factor(pm$model, labels[names(best)])
m_p <- ggplot() + geom_tile(data = subset(pm, zero), aes(lon, lat), fill = "grey80", width = 0.2, height = 0.2) +
  geom_tile(data = subset(pm, !zero), aes(lon, lat, fill = prevalence), width = 0.2, height = 0.2) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Expected\nprevalence", limits = c(0, NA)) +
  geom_sf(data = sudan_borders, fill = NA, colour = "grey25", linewidth = 0.3) + facet_wrap(~ model, nrow = 1) +
  coord_sf(expand = FALSE) + theme_minimal(base_size = 9) + labs(x = NULL, y = NULL)
ggsave(file.path(fig_dir, "directed_prevalence_maps.png"), m_p, width = 14, height = 4.8, dpi = 200, bg = "white")
