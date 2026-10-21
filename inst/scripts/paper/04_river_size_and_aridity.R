###############################################################################
# Effect of river size on breeding intensity, share of exposure by river size,
# aridity effect on the intensity of transmission and change-point profile.
###############################################################################

library(odispers); library(sf); library(ggplot2)
source("00_settings.R")
inp <- readRDS(file.path(out_dir, "inputs.rds"))
out <- readRDS(file.path(out_dir, "spline_fits.rds")); f <- out[[as.character(phi_main)]]$fit; d <- f$design
v <- inp$villages; nodes <- inp$nodes; pairs <- inp$pairs
cf <- coef(f); V <- vcov(f)
# contributions of each pair to exposure
used <- sort(unique(pairs$j)); Q <- nodes$DIS_AV_CMS[used][d$j]
w <- exp(cf["theta.log(DIS_AV_CMS)"] * log(Q) + cf["theta.I(log(DIS_AV_CMS)^2)"] * log(Q)^2)
t <- d$r / phi_main; k <- ifelse(t < 1, (1 - t)^4 * (1 + 4 * t), 0)
contrib <- d$len[d$j] * w * k / (2/3 * phi_main)
qcls <- cut(Q, c(0, 1, 10, 100, Inf), labels = c("< 1", "1-10", "10-100", "> 100"), right = FALSE)
acls <- aridity_classes(v$AI)[d$i]
share <- prop.table(tapply(contrib, list(acls, qcls), sum), 1)
cat("Share of exposure by river size (rows: aridity class of village):\n"); print(round(share, 3))
cat("Overall:\n"); print(round(prop.table(tapply(contrib, qcls, sum)), 3))
# river length within reach by size class
len_share <- prop.table(tapply(d$len[d$j] * k, list(acls, qcls), sum), 1)
cat("Share of kernel-weighted river LENGTH by size class:\n"); print(round(len_share, 3))
# breeding intensity curve with delta-method CI, relative to 1 m3/s
lq <- seq(log(0.1), log(3000), length.out = 200)
X <- cbind(lq, lq^2); b <- cf[c("theta.log(DIS_AV_CMS)", "theta.I(log(DIS_AV_CMS)^2)")]
Vb <- V[c("theta.log(DIS_AV_CMS)", "theta.I(log(DIS_AV_CMS)^2)"), c("theta.log(DIS_AV_CMS)", "theta.I(log(DIS_AV_CMS)^2)")]
eta <- X %*% b; se <- sqrt(rowSums((X %*% Vb) * X))
cur <- data.frame(Q = exp(lq), w = exp(eta), lo = exp(eta - 1.96 * se), hi = exp(eta + 1.96 * se))
cat("relative breeding intensity at Q = 0.1, 1, 10, 28.7, 100, 1000:\n")
print(round(exp(cbind(log(c(0.1,1,10,28.7,100,1000)), log(c(0.1,1,10,28.7,100,1000))^2) %*% b), 3))
p1 <- ggplot(cur, aes(Q, w)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") + geom_line(linewidth = 0.9) +
  scale_x_log10(breaks = c(0.1, 1, 10, 100, 1000), labels = c("0.1", "1", "10", "100", "1000")) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 100), labels = c("0.01", "0.1", "1", "10", "100")) +
  labs(x = "Long-term average discharge of the reach (m³/s, log scale)",
       y = "Breeding intensity per km of river\n(relative to a reach of 1 m³/s, log scale)",
       title = "Estimated effect of river size on breeding intensity", subtitle = "Range fixed at 40 km; band: 95% confidence interval (binomial variation only)") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "river_size_effect.png"), p1, width = 7, height = 4.5, dpi = 200)
sh <- as.data.frame(as.table(share)); names(sh) <- c("aridity", "size", "share")
p2 <- ggplot(sh, aes(aridity, share, fill = size)) + geom_col() +
  scale_fill_manual(values = c("< 1" = "#c6def1", "1-10" = "#7fb2dc", "10-100" = "#3f7fb5", "> 100" = "#1f4e79"),
                    name = "Discharge of reach (m³/s)") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Aridity class of the village", y = "Share of exposure", title = "Where exposure comes from: share contributed by rivers of each size") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "exposure_share_by_size.png"), p2, width = 7.5, height = 4.5, dpi = 200)
# aridity effect curve and change-point profile
cp <- out[[as.character(phi_main)]]$cp; la <- seq(log(min(v$AI)), log(max(v$AI)), length.out = 300)
Xa <- cbind(la, pmax(la - cp, 0)); ba <- cf[c("beta.log(AI)", "beta.h")]; Va <- V[c("beta.log(AI)", "beta.h"), c("beta.log(AI)", "beta.h")]
ref <- c(log(0.5), max(log(0.5) - cp, 0)); Xr <- sweep(Xa, 2, ref)
ea <- Xr %*% ba; sea <- sqrt(rowSums((Xr %*% Va) * Xr))
ca <- data.frame(AI = exp(la), m = exp(ea), lo = exp(ea - 1.96 * sea), hi = exp(ea + 1.96 * sea))
p3 <- ggplot(ca, aes(AI, m)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") + geom_line(linewidth = 0.9) +
  geom_vline(xintercept = exp(cp), linetype = 2) +
  geom_rug(data = data.frame(AI = v$AI), aes(AI), inherit.aes = FALSE, alpha = 0.3) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.05, 0.2, 0.5, 1), labels = c("0.001", "0.01", "0.05", "0.2", "0.5", "1")) +
  scale_y_log10() +
  labs(x = "Aridity index (log scale; ticks: villages)", y = "Multiplier of the force of infection\n(relative to AI = 0.5, log scale)",
       title = "Effect of aridity on the intensity of transmission", subtitle = sprintf("Linear spline in log AI; dashed line: estimated change point (AI = %.3f)", exp(cp))) +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "aridity_spline.png"), p3, width = 7, height = 4.5, dpi = 200)
pr <- out[[as.character(phi_main)]]$prof
p4 <- ggplot(pr, aes(exp(cp), logLik - max(logLik))) + geom_line() + geom_point(size = 0.8) +
  geom_hline(yintercept = -qchisq(0.95, 1) / 2, linetype = 2) + scale_x_log10() + coord_cartesian(ylim = c(-60, 2)) +
  labs(x = "Change point (aridity index, log scale)", y = "Profile log-likelihood\n(relative to maximum)", title = "Profile likelihood for the change point") +
  theme_minimal(base_size = 11)
ggsave(file.path(fig_dir, "changepoint_profile.png"), p4, width = 6, height = 4, dpi = 200)
