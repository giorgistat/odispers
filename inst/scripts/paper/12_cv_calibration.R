###############################################################################
# Cross-validated calibration of the model of 08_fit_u_linear.R (aridity
# linear in log AI, breeding linear in log discharge, field u on the river
# network, village effect Z, range of dispersal 40 km), compared with the
# model with the village effect only.
#
# For each held-out village the predictive distribution of the number of
# positives combines the Laplace (Gaussian) posterior of the latent effects
# fitted to the other villages, a new village effect and binomial sampling.
# The variance parameters and the range of u are fixed at their estimates
# from all the data. Ten random folds.
#
# Assessment, split as in a hurdle:
#   1. calibration among villages with at least one positive: non-randomised
#      PIT (Czado, Gneiting and Held, 2009, Biometrics 65, 1254-1261) of the
#      predictive distribution conditional on Y >= 1, with an envelope from
#      counts simulated from the predictive distributions;
#   2. zero against at least one positive: AUC, Brier score, reliability
#      diagram of P(Y >= 1), expected and observed numbers of zero villages;
#   3. overall: mean log predictive probability of the observed count.
###############################################################################

library(odispers)
library(sf)
library(Matrix)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

res08 <- load_result("u_fits_linear.rds")
inp <- load_result("inputs.rds")
v <- inp$villages; nodes <- inp$nodes
pairs <- river_pairs(v, nodes, phi_main)
graph <- river_graph(nodes, vertex_spacing = 5)
best <- res08$best; fz <- res08$fit_z
hyp <- function(b, name, old) if (name %in% names(b$hyper)) b$hyper[[name]] else b$hyper[[old]]

S_all <- odispers:::odisp_u_setup(y ~ log(AI), ~ log(DIS_AV_CMS), v, nodes, graph, pairs,
                                  phi_main, "n", "AI", 10)
stopifnot(length(S_all$has) == nrow(v))
range_u <- hyp(best, "range", "range_u")
R <- wm_precision(graph, range_u, 1)
models <- list(
  "u + village effect" = list(field = "u", R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE),
                              sig_u = hyp(best, "sigma_u", "sigma_u"), sig_z = best$hyper[["sigma_z"]]),
  "village effect only" = list(field = "none", sig_z = fz$hyper[["sigma_z"]]))
cat(sprintf("u: range %g km, sd %.3f, village effect sd %.3f; village effect only: sd %.3f\n",
            range_u, models[[1]]$sig_u, models[[1]]$sig_z, models[[2]]$sig_z))

subset_setup <- function(S, idx) {
  S$y <- S$y[idx]; S$n <- S$n[idx]; S$X <- S$X[idx, , drop = FALSE]; S$K <- S$K[idx, , drop = FALSE]
  S$xy <- S$xy[idx, , drop = FALSE]; S$has <- S$has[idx]; S$lchoose <- sum(lchoose(S$n, S$y))
  S
}

# Gauss-Hermite nodes and weights for integrals against N(0, 1)
gh <- local({ m <- 40; J <- matrix(0, m, m)
  for (k in 1:(m - 1)) J[k, k + 1] <- J[k + 1, k] <- sqrt(k / 2)
  e <- eigen(J, symmetric = TRUE); list(x = sqrt(2) * e$values, w = e$vectors[1, ]^2) })

# predictive probabilities P(Y = 0..n) for linear predictor ~ N(m, s2)
pred_pmf <- function(m, s2, n) {
  p <- 1 - exp(-exp(m + sqrt(s2) * gh$x))
  as.vector(sapply(0:n, function(y) sum(gh$w * stats::dbinom(y, n, p))))
}

cv_fold <- function(fld, train, test) {
  S_tr <- subset_setup(S_all, train); S_te <- subset_setup(S_all, test)
  x0 <- list(beta = best$coefficients[seq_along(S_all$m_beta)],
             theta = best$coefficients[length(S_all$m_beta) + seq_along(S_all$m_theta)],
             u = if (fld$field == "u") best$u else NULL, s = NULL, z = rep(0, length(train)))
  fit <- odispers:::odisp_u_laplace(fld, S_tr, x0)
  ev <- odispers:::odisp_u_eta(fit$x, S_tr); lk <- odispers:::odisp_u_lik(ev$eta, S_tr)
  ops <- odispers:::odisp_u_ops(fit$x, ev, lk, S_tr, fld)
  ev_te <- odispers:::odisp_u_eta(list(beta = fit$x$beta, theta = fit$x$theta, u = fit$x$u, s = NULL, z = NULL), S_te)
  Jt <- as.matrix(ev_te$V %*% S_te$B)
  Ju <- if (fld$field == "u") as.matrix(ev_te$V %*% S_te$A) else NULL
  lapply(seq_along(test), function(k) {
    j <- c(S_te$X[k, ], Jt[k, ], if (!is.null(Ju)) Ju[k, ], rep(0, length(train)))
    s2 <- sum(j * ops$HGNinv(j)) + fld$sig_z^2      # latent uncertainty + new village effect
    list(village = test[k], y = S_te$y[k], n = S_te$n[k], pmf = pred_pmf(ev_te$eta[k], s2, S_te$n[k]))
  })
}

set.seed(2026)
nv <- nrow(v)
folds <- list("random folds" = sample(rep(1:10, length.out = nv)))
preds <- list()
for (sch in names(folds)) for (m in names(models)) {
  f <- folds[[sch]]
  cat(sprintf("%s, %s: ", sch, m))
  out <- unlist(lapply(sort(unique(f)), function(k) { cat(k, ""); cv_fold(models[[m]], which(f != k), which(f == k)) }),
                recursive = FALSE)
  cat("\n")
  preds[[paste(sch, m, sep = " | ")]] <- out[order(sapply(out, `[[`, "village"))]
}

# ---- 1. non-randomised PIT among villages with at least one positive --------
# Conditional predictive distribution given Y >= 1, and the PIT function
# F(u | y) of Czado et al. (2009); the mean PIT is averaged over villages.
cond_cdf <- function(pmf) { p0 <- pmf[1]; c(0, cumsum(pmf[-1]) / (1 - p0)) }   # index y + 1, y = 0..n
pit_fun <- function(Fm1, F0, u) ifelse(u <= Fm1, 0, ifelse(u >= F0, 1, (u - Fm1) / pmax(F0 - Fm1, 1e-12)))
mean_pit_hist <- function(ys, cdfs, J = 10) {
  ub <- (0:J) / J
  Fbar <- rowMeans(sapply(seq_along(ys), function(i) pit_fun(cdfs[[i]][ys[i]], cdfs[[i]][ys[i] + 1], ub)))
  J * diff(Fbar)                                                 # density; 1 under calibration
}
pit_tab <- list(); env_tab <- list()
for (nm in names(preds)) {
  pr <- preds[[nm]]
  pos <- Filter(function(p) p$y > 0, pr)
  cdfs <- lapply(pos, function(p) cond_cdf(p$pmf))
  h <- mean_pit_hist(sapply(pos, `[[`, "y"), cdfs)
  # envelope: counts simulated from the conditional predictive distributions
  sims <- replicate(200, { ys <- sapply(pos, function(p) { q <- p$pmf[-1] / (1 - p$pmf[1]); sample.int(length(q), 1, prob = q) })
                           mean_pit_hist(ys, cdfs) })
  pit_tab[[nm]] <- data.frame(analysis = nm, bin = seq_along(h), mid = (seq_along(h) - 0.5) / length(h), density = h,
                              lower = apply(sims, 1, stats::quantile, 0.025), upper = apply(sims, 1, stats::quantile, 0.975))
}
pit_tab <- do.call(rbind, pit_tab)
pit_tab$scheme <- sub(" \\|.*", "", pit_tab$analysis); pit_tab$model <- sub(".*\\| ", "", pit_tab$analysis)

p_pit <- ggplot(pit_tab, aes(mid, density)) +
  geom_col(width = 0.1, fill = "#7fbf7b", colour = "white") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.04, colour = "grey30") +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_grid(model ~ scheme) +
  labs(x = "Probability integral transform", y = "Density",
       title = "Non-randomised PIT among villages with at least one positive",
       subtitle = "Predictive distribution conditional on Y >= 1; bars: 95% range under calibration (200 simulations)") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "cv_pit_positive.png"), p_pit, width = 8, height = 6, dpi = 200, bg = "white")

# ---- 2. zero against at least one positive -----------------------------------
auc <- function(p, o) { r <- rank(p); n1 <- sum(o); n0 <- sum(!o); (sum(r[o]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
cls <- cut(v$AI, c(0, 0.05, 0.2, 0.5, Inf), labels = c("Hyper-arid", "Arid", "Semi-arid", "Sub-humid/humid"))
rel_tab <- list(); scores <- list(); zero_tab <- list()
for (nm in names(preds)) {
  pr <- preds[[nm]]
  p1 <- sapply(pr, function(p) 1 - p$pmf[1]); o1 <- sapply(pr, `[[`, "y") > 0
  logs <- sapply(pr, function(p) log(max(p$pmf[p$y + 1], 1e-300)))
  scores[[nm]] <- data.frame(analysis = nm, AUC = auc(p1, o1), Brier = mean((p1 - o1)^2),
                             zeros_observed = sum(!o1), zeros_expected = sum(1 - p1),
                             mean_log_score = mean(logs))
  br <- unique(stats::quantile(p1, seq(0, 1, 0.1)))
  b <- cut(p1, br, include.lowest = TRUE)
  rel <- data.frame(pred = tapply(p1, b, mean), obs = tapply(o1, b, mean), n = as.vector(table(b)))
  rel$lower <- stats::qbeta(0.025, rel$obs * rel$n + 0.5, (1 - rel$obs) * rel$n + 0.5)
  rel$upper <- stats::qbeta(0.975, rel$obs * rel$n + 0.5, (1 - rel$obs) * rel$n + 0.5)
  rel_tab[[nm]] <- cbind(analysis = nm, rel)
  zero_tab[[nm]] <- data.frame(analysis = nm, class = levels(cls),
                               observed = as.vector(tapply(!o1, cls, sum)),
                               expected = round(as.vector(tapply(1 - p1, cls, sum)), 1))
}
scores <- do.call(rbind, scores); rownames(scores) <- NULL
cat("\nZero against at least one positive, and overall log score:\n")
print(within(scores, { AUC <- round(AUC, 3); Brier <- round(Brier, 4); zeros_expected <- round(zeros_expected, 1)
                       mean_log_score <- round(mean_log_score, 3) }), row.names = FALSE)
cat("\nZero villages by aridity class:\n"); print(do.call(rbind, zero_tab), row.names = FALSE)

rel_tab <- do.call(rbind, rel_tab)
rel_tab$scheme <- sub(" \\|.*", "", rel_tab$analysis); rel_tab$model <- sub(".*\\| ", "", rel_tab$analysis)
p_rel <- ggplot(rel_tab, aes(pred, obs, colour = model)) +
  geom_abline(linetype = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, alpha = 0.6) +
  geom_point(size = 2) + geom_line() +
  facet_wrap(~ scheme) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = c("u + village effect" = "#1b7837", "village effect only" = "grey50"), name = NULL) +
  labs(x = "Predicted probability of at least one positive (deciles)",
       y = "Observed proportion of villages with at least one positive",
       title = "Reliability of the predicted probability of at least one positive",
       subtitle = "Cross-validated; bars: 95% intervals") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, "cv_reliability_any_positive.png"), p_rel, width = 9, height = 5, dpi = 200, bg = "white")

saveRDS(list(preds = preds, folds = folds, pit = pit_tab, reliability = rel_tab, scores = scores,
             zeros = do.call(rbind, zero_tab)), file.path(out_dir, "cv_calibration.rds"))
