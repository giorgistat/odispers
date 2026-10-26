###############################################################################
# The village averages of the river covariates (kernel range 40 km, theta = 0;
# from 31_covariates_eda.R) against the two parts of the zero-inflated model,
# for all villages:
#   (1) the probability of observing at least one positive: indicator of
#       y > 0 against the covariate, with decile proportions (Wilson 95%
#       intervals) and a binomial GAM smooth on the logit scale;
#   (2) prevalence among the villages with positives, on the complementary
#       log-log scale, log(-log(1 - p)), the scale of log lambda in the
#       model (p = 1 - exp(-lambda)), with p = (y + 0.5)/(n + 1), and a
#       quasi-binomial GAM smooth with the cloglog link (extra-binomial
#       variation between villages widens the band).
# A panel for log E0, the log exposure from river length alone (sum of kernel
# weights times river length within 40 km), shows the association of the
# exposure itself with the two parts.
# Under the model, a covariate acting only through E shifts log E by about
# beta * xbar, hence logit(pi) by a_E * beta * xbar and log(lambda) by
# gamma * beta * xbar; the GAM slopes are also summarised by straight-line
# (GLM) slopes per standard deviation of the village average.
###############################################################################
library(odispers)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
e <- readRDS(file.path(out_dir, "cov_eda.rds"))
d <- e$villages; X <- e$xbar$main[e$ok, ]
net <- load_result("inputs_full_network.rds"); nodes <- net$nodes; pr <- net$pairs
kern <- function(t) ifelse(t < 1, (1 - t)^4 * (1 + 4 * t), 0)
E0 <- tapply(nodes$len[pr$j] * kern(pr$r / 40), factor(pr$i, levels = seq_len(pr$m)), sum)
X <- cbind(log_E0 = log(as.vector(E0)[e$ok]), X)
labels <- c(log_E0 = "log E0: kernel-weighted river length",
            log_slope = "log(1 + slope, dm/km)  [+]", steep_share = "share of river with slope >= 2 m/km  [+]",
            velocity = "log velocity proxy  [+]", permanence = "flow permanence (min/mean)  [+]",
            inundation = "log(1 + inundation %)  [-]", hard_rock = "share of river on hard rock  [+]",
            forest = "upstream forest cover %  [+]", lakes = "log(1 + lake %)  [?]",
            regulation = "log(1 + regulation %)  [?]", elevation = "elevation (100 m)  [none]",
            runoff = "log(1 + runoff, mm/yr)  [none]", log_discharge = "log discharge (in the model)")
ci_wilson <- function(s, n) { p <- s / n; z <- 1.96; c0 <- p + z^2 / (2 * n); h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  cbind((c0 - h) / (1 + z^2 / n), (c0 + h) / (1 + z^2 / n)) }
pts <- NULL; bins <- NULL; sm1 <- NULL; sm2 <- NULL; slopes <- NULL
pos <- d$y > 0
for (k in names(labels)) {
  x <- X[, k]; lab <- labels[[k]]
  pts <- rbind(pts, data.frame(covariate = lab, x = x, pos = as.numeric(pos), n = d$n, y = d$y,
                               cll = log(-log(1 - (d$y + 0.5) / (d$n + 1)))))
  br <- unique(quantile(x, seq(0, 1, 0.1)))
  b <- if (length(br) > 2) cut(x, br, include.lowest = TRUE) else factor(x)
  for (lv in levels(b)) { s <- b == lv; if (!any(s)) next; w <- ci_wilson(sum(pos[s]), sum(s))
    bins <- rbind(bins, data.frame(covariate = lab, x = mean(x[s]), p = mean(pos[s]), lo = w[1], hi = w[2])) }
  xs <- seq(min(x), max(x), length.out = 100)
  kk <- min(8, length(unique(x)) - 1)
  g1 <- mgcv::gam(pos ~ s(x, k = kk), family = binomial, data = data.frame(x = x, pos = pos), method = "REML")
  p1 <- predict(g1, data.frame(x = xs), se.fit = TRUE)
  sm1 <- rbind(sm1, data.frame(covariate = lab, x = xs, fit = plogis(p1$fit), lo = plogis(p1$fit - 1.96 * p1$se.fit),
                               hi = plogis(p1$fit + 1.96 * p1$se.fit)))
  xp <- x[pos]; kk2 <- min(8, length(unique(xp)) - 1)
  g2 <- mgcv::gam(cbind(y, n - y) ~ s(x, k = kk2), family = quasibinomial(link = "cloglog"),
                  data = data.frame(x = xp, y = d$y[pos], n = d$n[pos]), method = "REML")
  p2 <- predict(g2, data.frame(x = xs[xs >= min(xp) & xs <= max(xp)]), se.fit = TRUE)
  sm2 <- rbind(sm2, data.frame(covariate = lab, x = xs[xs >= min(xp) & xs <= max(xp)], fit = p2$fit,
                               lo = p2$fit - 1.96 * p2$se.fit, hi = p2$fit + 1.96 * p2$se.fit))
  # straight-line slopes per standard deviation of the village average
  z <- (x - mean(x)) / sd(x)
  s1 <- coef(summary(glm(pos ~ z, family = binomial)))[2, 1:2]
  s2 <- coef(summary(glm(cbind(d$y[pos], d$n[pos] - d$y[pos]) ~ z[pos], family = quasibinomial(link = "cloglog"))))[2, 1:2]
  slopes <- rbind(slopes, data.frame(covariate = k, zone_logit = s1[1], zone_se = s1[2],
                                     intensity_cloglog = s2[1], intensity_se = s2[2], ratio = s1[1] / s2[1]))
}
for (dd in c("pts", "bins", "sm1", "sm2")) { z <- get(dd); z$covariate <- factor(z$covariate, labels); assign(dd, z) }
print(slopes, digits = 2, row.names = FALSE)

p_zone <- ggplot() +
  geom_jitter(data = pts, aes(x, pos), width = 0, height = 0.03, size = 0.5, alpha = 0.25, colour = "grey30") +
  geom_ribbon(data = sm1, aes(x, ymin = lo, ymax = hi), fill = "#1b7837", alpha = 0.2) +
  geom_line(data = sm1, aes(x, fit), colour = "#1b7837", linewidth = 0.9) +
  geom_pointrange(data = bins, aes(x, p, ymin = lo, ymax = hi), size = 0.2, colour = "black") +
  facet_wrap(~ covariate, scales = "free_x", ncol = 4) + theme_minimal(base_size = 9) +
  labs(x = "village average over the rivers within 40 km", y = "probability of at least one positive")
ggsave(file.path(fig_dir, "cov_part_zone.png"), p_zone, width = 11, height = 10, dpi = 180, bg = "white")

p_int <- ggplot() +
  geom_point(data = subset(pts, pos == 1), aes(x, cll), size = 0.6, alpha = 0.35, colour = "grey30") +
  geom_ribbon(data = sm2, aes(x, ymin = lo, ymax = hi), fill = "#2166ac", alpha = 0.2) +
  geom_line(data = sm2, aes(x, fit), colour = "#2166ac", linewidth = 0.9) +
  facet_wrap(~ covariate, scales = "free_x", ncol = 4) + theme_minimal(base_size = 9) +
  labs(x = "village average over the rivers within 40 km",
       y = "prevalence among villages with positives, log(-log(1 - p))")
ggsave(file.path(fig_dir, "cov_part_intensity.png"), p_int, width = 11, height = 10, dpi = 180, bg = "white")
saveRDS(list(slopes = slopes, bins = bins, sm_zone = sm1, sm_intensity = sm2), file.path(out_dir, "cov_two_parts.rds"))
