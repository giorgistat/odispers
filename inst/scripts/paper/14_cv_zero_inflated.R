###############################################################################
# Cross-validated calibration of the zero-inflated model of
# 13_fit_zero_inflated.R (gamma free), with the same folds and the same
# assessment as 12_cv_calibration.R, and comparison with the model without
# zero inflation (u + village effect).
#
# For a held-out village, (a, b) = (log E, beta'd) are jointly Gaussian under
# the Laplace posterior from the training villages; the new village effect is
# added to b; the predictive probabilities
#   P(Y = 0) = E[1 - pi + pi (1 - p)^n],  P(Y = y) = E[pi Bin(y; n, p)]
# are computed by two-dimensional Gauss-Hermite quadrature.
###############################################################################

library(odispers)
library(sf)
library(Matrix)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

# which zero-inflated fit to assess: set zi_file and zi_label before sourcing
# this script to assess another one (see 15_fit_zi_aridity.R)
if (!exists("zi_file")) zi_file <- "zi_fits.rds"
if (!exists("zi_label")) zi_label <- "zero-inflated"
if (!exists("zi_suffix")) zi_suffix <- "_zi"
zi <- load_result(zi_file)$fits$free
cv12 <- load_result("cv_calibration.rds")
inp <- load_result("inputs.rds")
v <- inp$villages; nodes <- inp$nodes
pairs <- river_pairs(v, nodes, phi_main)
graph <- river_graph(nodes, vertex_spacing = 5)
S_all <- zi$setup
fld <- zi$fld
if (is.null(fld$delta)) fld$delta <- 0
if (is.null(S_all$lai)) S_all$lai <- log(v$AI)[S_all$has]
cat(sprintf("%s model: alpha %.4g, gamma %.3f, delta %.3f, sigma_u %.3f, sigma_z %.3f\n",
            zi_label, fld$alpha, fld$gamma, fld$delta, fld$sig_u, fld$sig_z))

subset_setup <- function(S, idx) {
  S$y <- S$y[idx]; S$n <- S$n[idx]; S$X <- S$X[idx, , drop = FALSE]; S$K <- S$K[idx, , drop = FALSE]
  S$xy <- S$xy[idx, , drop = FALSE]; S$has <- S$has[idx]; S$lchoose <- sum(lchoose(S$n, S$y))
  S$lai <- S$lai[idx]
  S
}
gh <- local({ m <- 20; J <- matrix(0, m, m)
  for (k in 1:(m - 1)) J[k, k + 1] <- J[k + 1, k] <- sqrt(k / 2)
  e <- eigen(J, symmetric = TRUE); list(x = sqrt(2) * e$values, w = e$vectors[1, ]^2) })
g2 <- expand.grid(i = seq_along(gh$x), j = seq_along(gh$x))
g2x1 <- gh$x[g2$i]; g2x2 <- gh$x[g2$j]; g2w <- gh$w[g2$i] * gh$w[g2$j]

zi_pmf <- function(ma, mb, Sig, n, lai) {
  L <- t(chol(Sig + diag(1e-10, 2)))
  a <- ma + L[1, 1] * g2x1; b <- mb + L[2, 1] * g2x1 + L[2, 2] * g2x2
  pi_ <- -expm1(-fld$alpha * exp(a + fld$delta * lai)); p <- -expm1(-exp(fld$gamma * a + b))
  pm <- sapply(0:n, function(y) sum(g2w * pi_ * stats::dbinom(y, n, p)))
  pm[1] <- pm[1] + sum(g2w * (1 - pi_))
  pm
}

cv_fold <- function(train, test) {
  S_tr <- subset_setup(S_all, train); S_te <- subset_setup(S_all, test)
  x0 <- list(beta = zi$coefficients[seq_along(S_all$m_beta)],
             theta = zi$coefficients[length(S_all$m_beta) + seq_along(S_all$m_theta)],
             u = zi$u, s = NULL, z = rep(0, length(train)))
  fit <- odispers:::odisp_zi_laplace(fld, S_tr, x0)
  ev <- odispers:::odisp_zi_eval(fit$x, S_tr)
  lk <- odispers:::odisp_zi_lik(ev$a, ev$b, S_tr, fld$gamma, fld$alpha, fld$delta)
  ops <- odispers:::odisp_zi_ops(fit$x, ev, lk, S_tr, fld)
  ev_te <- odispers:::odisp_zi_eval(list(beta = fit$x$beta, theta = fit$x$theta, u = fit$x$u, s = NULL,
                                         z = rep(0, length(test))), S_te)
  Jt <- as.matrix(ev_te$V %*% S_te$B); Ju <- as.matrix(ev_te$V %*% S_te$A)
  p <- ncol(S_te$X); q <- ncol(S_te$B); nz <- length(train)
  lapply(seq_along(test), function(k) {
    ja <- c(rep(0, p), Jt[k, ], Ju[k, ], rep(0, nz))
    jb <- c(S_te$X[k, ], rep(0, q), rep(0, ncol(Ju)), rep(0, nz))
    ha <- ops$HGNinv(ja); hb <- ops$HGNinv(jb)
    Sig <- matrix(c(sum(ja * ha), sum(ja * hb), sum(jb * ha), sum(jb * hb)), 2)
    Sig <- (Sig + t(Sig)) / 2; Sig[2, 2] <- Sig[2, 2] + fld$sig_z^2
    list(village = test[k], y = S_te$y[k], n = S_te$n[k],
         pmf = zi_pmf(ev_te$a[k], ev_te$b[k], Sig, S_te$n[k], S_te$lai[k]))
  })
}

preds <- list()
for (sch in names(cv12$folds)) {
  f <- cv12$folds[[sch]]
  cat(sprintf("%s, %s: ", sch, zi_label))
  out <- unlist(lapply(sort(unique(f)), function(k) { cat(k, ""); cv_fold(which(f != k), which(f == k)) }),
                recursive = FALSE)
  cat("\n")
  preds[[paste(sch, paste(zi_label, "u + village effect"), sep = " | ")]] <- out[order(sapply(out, `[[`, "village"))]
}
# the model without zero inflation, from 12_cv_calibration.R
for (nm in grep("u \\+ village effect", names(cv12$preds), value = TRUE)) preds[[nm]] <- cv12$preds[[nm]]
# and the zero-inflated model without aridity in pi, if assessing another one
if (zi_file != "zi_fits.rds" && file.exists(file.path(out_dir, "cv_zero_inflated.rds"))) {
  prev <- readRDS(file.path(out_dir, "cv_zero_inflated.rds"))$preds
  for (nm in grep("^[^|]*\\| zero-inflated u", names(prev), value = TRUE)) preds[[nm]] <- prev[[nm]]
}

# ---- assessment, as in 12_cv_calibration.R ------------------------------------
cond_cdf <- function(pmf) { p0 <- pmf[1]; c(0, cumsum(pmf[-1]) / (1 - p0)) }
pit_fun <- function(Fm1, F0, u) ifelse(u <= Fm1, 0, ifelse(u >= F0, 1, (u - Fm1) / pmax(F0 - Fm1, 1e-12)))
mean_pit_hist <- function(ys, cdfs, J = 10) {
  ub <- (0:J) / J
  Fbar <- rowMeans(sapply(seq_along(ys), function(i) pit_fun(cdfs[[i]][ys[i]], cdfs[[i]][ys[i] + 1], ub)))
  J * diff(Fbar)
}
auc <- function(p, o) { r <- rank(p); n1 <- sum(o); n0 <- sum(!o); (sum(r[o]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
cls <- aridity_classes(v$AI)
set.seed(2026)
pit_tab <- list(); rel_tab <- list(); scores <- list(); zero_tab <- list()
for (nm in names(preds)) {
  pr <- preds[[nm]]
  pos <- Filter(function(p) p$y > 0, pr)
  cdfs <- lapply(pos, function(p) cond_cdf(p$pmf))
  h <- mean_pit_hist(sapply(pos, `[[`, "y"), cdfs)
  sims <- replicate(200, { ys <- sapply(pos, function(p) { q <- p$pmf[-1] / (1 - p$pmf[1]); sample.int(length(q), 1, prob = q) })
                           mean_pit_hist(ys, cdfs) })
  pit_tab[[nm]] <- data.frame(analysis = nm, mid = (seq_along(h) - 0.5) / length(h), density = h,
                              lower = apply(sims, 1, stats::quantile, 0.025), upper = apply(sims, 1, stats::quantile, 0.975))
  p1 <- sapply(pr, function(p) 1 - p$pmf[1]); o1 <- sapply(pr, `[[`, "y") > 0
  logs <- sapply(pr, function(p) log(max(p$pmf[p$y + 1], 1e-300)))
  scores[[nm]] <- data.frame(analysis = nm, AUC = round(auc(p1, o1), 3), Brier = round(mean((p1 - o1)^2), 4),
                             zeros_observed = sum(!o1), zeros_expected = round(sum(1 - p1), 1),
                             mean_log_score = round(mean(logs), 3))
  br <- unique(stats::quantile(p1, seq(0, 1, 0.1))); b <- cut(p1, br, include.lowest = TRUE)
  rel <- data.frame(pred = tapply(p1, b, mean), obs = tapply(o1, b, mean), n = as.vector(table(b)))
  rel$lower <- stats::qbeta(0.025, rel$obs * rel$n + 0.5, (1 - rel$obs) * rel$n + 0.5)
  rel$upper <- stats::qbeta(0.975, rel$obs * rel$n + 0.5, (1 - rel$obs) * rel$n + 0.5)
  rel_tab[[nm]] <- cbind(analysis = nm, rel)
  zero_tab[[nm]] <- data.frame(analysis = nm, class = levels(cls), observed = as.vector(tapply(!o1, cls, sum)),
                               expected = round(as.vector(tapply(1 - p1, cls, sum)), 1))
}
scores <- do.call(rbind, scores); rownames(scores) <- NULL
cat("\nZero against at least one positive, and overall log score:\n"); print(scores, row.names = FALSE)
cat("\nZero villages by aridity class:\n"); print(do.call(rbind, zero_tab), row.names = FALSE)

split_lab <- function(d) { d$scheme <- sub(" \\|.*", "", d$analysis)
  d$model <- sub(" u \\+ village effect$", "", sub(".*\\| ", "", d$analysis))
  d$model[d$model == "u + village effect"] <- "no zero inflation"; d }
pit_tab <- split_lab(do.call(rbind, pit_tab)); rel_tab <- split_lab(do.call(rbind, rel_tab))
pit_tab$model <- factor(pit_tab$model, unique(c("no zero inflation", "zero-inflated", pit_tab$model)))
p_pit <- ggplot(pit_tab, aes(mid, density)) +
  geom_col(width = 0.1, fill = "#7fbf7b", colour = "white") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.04, colour = "grey30") +
  geom_hline(yintercept = 1, linetype = 2) + facet_grid(model ~ scheme) +
  labs(x = "Probability integral transform", y = "Density",
       title = "Non-randomised PIT among villages with at least one positive (u + village effect)",
       subtitle = "Predictive distribution conditional on Y >= 1; bars: 95% range under calibration") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, paste0("cv_pit_positive", zi_suffix, ".png")), p_pit, width = 8,
       height = 2.6 * length(unique(pit_tab$model)) + 1, dpi = 200, bg = "white")
p_rel <- ggplot(rel_tab, aes(pred, obs, colour = model)) +
  geom_abline(linetype = 2) + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, alpha = 0.6) +
  geom_point(size = 2) + geom_line() + facet_wrap(~ scheme) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_colour_manual(values = c("zero-inflated" = "#7fbf7b", "no zero inflation" = "grey50",
                                 "zero-inflated, aridity in pi" = "#1b7837"), name = NULL) +
  labs(x = "Predicted probability of at least one positive (deciles)",
       y = "Observed proportion with at least one positive",
       title = "Reliability of the predicted probability of at least one positive",
       subtitle = "Cross-validated; bars: 95% intervals") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(fig_dir, paste0("cv_reliability_any_positive", zi_suffix, ".png")), p_rel, width = 9, height = 5,
       dpi = 200, bg = "white")
saveRDS(list(preds = preds, scores = scores, zeros = do.call(rbind, zero_tab), pit = pit_tab, reliability = rel_tab),
        file.path(out_dir, if (zi_suffix == "_zi") "cv_zero_inflated.rds" else paste0("cv", zi_suffix, ".rds")))
