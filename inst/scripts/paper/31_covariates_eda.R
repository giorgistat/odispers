###############################################################################
# Model-free exploratory analysis of the RiverATLAS covariates as candidate
# breeding covariates in the exposure E.
#
# If log w(s) = theta log Q(s) + beta x(s), then for small beta
#   log E_i(beta) ~ log E_i(0) + beta * xbar_i,
#   xbar_i = sum_j K_ij w_j x_j / sum_j K_ij w_j,
# the average of x over the rivers within reach of village i, weighted by the
# dispersal kernel (and by breeding from discharge). xbar_i is therefore the
# village-level quantity through which a river covariate acts on prevalence,
# and it is computed here without fitting any model (theta = 0, i.e. weights
# proportional to river length times the kernel; theta = -0.9 as a check).
#
# For each covariate: one figure with (a) xbar_i at the villages, (b) the
# share of villages with at least one positive (zone of transmission) and
# (c) the pooled prevalence among villages with positives (intensity), by
# quintile of xbar_i, overall and within terciles of aridity. A summary table
# gives Spearman correlations with 95% spatial block-bootstrap intervals
# (1-degree cells), the sign within aridity terciles, the sensitivity to the
# kernel range (20, 40, 60 km) and theta, and the correlations between the
# xbar_i of different covariates. These feed the inclusion criterion.
###############################################################################
library(odispers)
library(sf)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
inp <- load_result("inputs.rds"); v <- inp$villages
net <- load_result("inputs_full_network.rds"); nodes <- net$nodes
data(sudan_river_covariates); data(sudan_borders)
rc <- as.data.frame(sudan_river_covariates)[match(nodes$HYRIV_ID, sudan_river_covariates$HYRIV_ID), ]

# covariates on the river (per node); expected sign of the effect on breeding
hard <- c(2, 4, 7, 8, 9, 10, 12, 13)                 # GLiM basement / volcanic classes
cv <- data.frame(
  log_discharge = log(nodes$DIS_AV_CMS),
  log_slope     = log1p(rc$sgr_dk_rav),
  steep_share   = as.numeric(rc$sgr_dk_rav >= 20),   # slope >= 2 m/km
  velocity      = 0.27 * log(nodes$DIS_AV_CMS) + 0.5 * log((rc$sgr_dk_rav + 0.5) / 1e4),
  permanence    = pmin(rc$dis_m3_pmn / pmax(rc$dis_m3_pyr, 1e-6), 1),
  inundation    = log1p(rc$inu_pc_ult),
  hard_rock     = as.numeric(rc$lit_cl_cmj %in% hard),
  forest        = rc$for_pc_use,
  lakes         = log1p(rc$lka_pc_use),
  regulation    = log1p(rc$dor_pc_pva),
  elevation     = rc$ele_mt_cav / 100,
  runoff        = log1p(rc$run_mm_cyr))
info <- data.frame(
  covariate = names(cv),
  label = c("log discharge (already in the model)", "log(1 + channel slope, dm/km)", "share of river with slope >= 2 m/km",
            "log velocity proxy, 0.27 log Q + 0.5 log S", "flow permanence (min month / mean)",
            "log(1 + upstream inundation %)", "share of river on hard rock", "upstream forest cover %",
            "log(1 + upstream lake %)", "log(1 + degree of regulation %)", "catchment elevation (100 m)",
            "log(1 + runoff, mm/yr)"),
  expected = c("?", "+", "+", "+", "+", "-", "+", "+", "?", "?", "none", "none"),
  rationale = c("primary", "primary", "primary", "primary", "primary", "primary", "secondary", "secondary",
                "weak", "weak", "weak", "weak"))

kern <- function(t) ifelse(t < 1, (1 - t)^4 * (1 + 4 * t), 0)
xbar_of <- function(pairs, phi, theta) {
  keep <- pairs$r < phi
  i <- pairs$i[keep]; j <- pairs$j[keep]
  k <- nodes$len[j] * kern(pairs$r[keep] / phi) * nodes$DIS_AV_CMS[j]^theta
  den <- tapply(k, factor(i, levels = seq_len(pairs$m)), sum)
  out <- sapply(cv, function(x) tapply(k * x[j], factor(i, levels = seq_len(pairs$m)), sum) / den)
  out[!is.finite(den) | den == 0, ] <- NA
  out
}
pairs60 <- river_pairs(v, nodes, max_range = 60)
X <- list(main = xbar_of(net$pairs, 40, 0), phi20 = xbar_of(net$pairs, 20, 0),
          phi60 = xbar_of(pairs60, 60, 0), theta = xbar_of(net$pairs, 40, -0.9))
ok <- complete.cases(X$main)
cat("villages with river within 40 km:", sum(ok), "of", nrow(v), "; positives among the others:",
    sum(v$y[!ok] > 0), "\n")

d <- st_drop_geometry(v)[ok, c("long", "lat", "n", "y", "AI")]
d$pos <- as.numeric(d$y > 0)
d$aridity <- cut(d$AI, quantile(d$AI, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
                 labels = c("dry tercile", "middle tercile", "humid tercile"))
d$cell <- paste(floor(d$long), floor(d$lat))                # 1-degree cells for the block bootstrap
d$elogit <- log((d$y + 0.5) / (d$n - d$y + 0.5))
cat("villages per 1-degree cell: median", median(table(d$cell)), "; cells:", length(unique(d$cell)), "\n")

# Spearman correlations: part 1 with the indicator of any positive (all villages),
# part 2 with the empirical logit among villages with positives
rho <- function(x, dd) c(zone = suppressWarnings(cor(x, dd$pos, method = "spearman", use = "complete.obs")),
                         intensity = suppressWarnings(cor(x[dd$pos == 1], dd$elogit[dd$pos == 1], method = "spearman",
                                                          use = "complete.obs")))
set.seed(2026)
cells <- unique(d$cell); B <- 2000
boot_idx <- replicate(B, unlist(split(seq_len(nrow(d)), d$cell)[sample(cells, replace = TRUE)], use.names = FALSE),
                      simplify = FALSE)
summ <- NULL
for (k in names(cv)) {
  x <- X$main[ok, k]
  r0 <- rho(x, d)
  bs <- sapply(boot_idx, function(id) rho(x[id], d[id, ]))
  ci <- apply(bs, 1, quantile, c(0.025, 0.975), na.rm = TRUE)
  strat <- sapply(levels(d$aridity), function(a) { s <- d$aridity == a; rho(x[s], d[s, ]) })
  sens <- sapply(c("phi20", "phi60", "theta"), function(s) rho(X[[s]][ok, k], d))
  for (part in c("zone", "intensity"))
    summ <- rbind(summ, data.frame(covariate = k, part = part, rho = r0[[part]],
      lower = ci[1, part], upper = ci[2, part],
      same_sign_strata = sum(sign(strat[part, ]) == sign(r0[[part]])),
      rho_dry = strat[part, 1], rho_mid = strat[part, 2], rho_humid = strat[part, 3],
      rho_phi20 = sens[part, "phi20"], rho_phi60 = sens[part, "phi60"], rho_theta = sens[part, "theta"],
      cor_main_phi20 = cor(x, X$phi20[ok, k], use = "complete.obs"), cor_main_phi60 = cor(x, X$phi60[ok, k])))
}
summ <- merge(summ, info[, c("covariate", "expected", "rationale")], sort = FALSE)
print(summ[, 1:6], digits = 2, row.names = FALSE)
print(summ[, c(1, 2, 7:9)], digits = 2, row.names = FALSE)
print(summ[, c(1, 2, 10:16)], digits = 2, row.names = FALSE)
cat("villages with river within 20 km:", sum(complete.cases(X$phi20[ok, ])), "\n")
cmat <- cor(X$main[ok, ], method = "spearman")
cat("\nSpearman correlations between the village averages:\n"); print(round(cmat, 2))

# one figure per covariate
ci_wilson <- function(s, n) { p <- s / n; z <- 1.96; c0 <- p + z^2 / (2 * n); h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  cbind((c0 - h) / (1 + z^2 / n), (c0 + h) / (1 + z^2 / n)) }
binned <- function(x, dd, nb = 5) {
  br <- unique(quantile(x, seq(0, 1, length.out = nb + 1)))
  b <- if (length(br) > 2) cut(x, br, include.lowest = TRUE) else factor(x)
  do.call(rbind, lapply(split(seq_along(x), b), function(id) {
    if (!length(id)) return(NULL)
    np <- sum(dd$pos[id]); w1 <- ci_wilson(np, length(id))
    pp <- id[dd$pos[id] == 1]
    sy <- sum(dd$y[pp]); sn <- sum(dd$n[pp])
    w2 <- if (sn > 0) ci_wilson(sy, sn) else cbind(NA, NA)
    data.frame(x = mean(x[id]), zone = np / length(id), zl = w1[1], zu = w1[2], nv = length(id),
               prev = if (sn > 0) sy / sn else NA, pl = w2[1], pu = w2[2], npos = length(pp))
  }))
}
dir.create(file.path(fig_dir, "eda"), showWarnings = FALSE)
for (k in names(cv)) {
  x <- X$main[ok, k]; lab <- info$label[info$covariate == k]
  b_all <- cbind(binned(x, d), stratum = "all villages")
  b_str <- do.call(rbind, lapply(levels(d$aridity), function(a) { s <- d$aridity == a
    cbind(binned(x[s], d[s, ], nb = 4), stratum = a) }))
  bb <- rbind(b_all, b_str)
  bb$stratum <- factor(bb$stratum, c("all villages", levels(d$aridity)))
  cols <- c("all villages" = "black", "dry tercile" = "#d6604d", "middle tercile" = "#999999", "humid tercile" = "#4393c3")
  pa <- ggplot() + geom_sf(data = sudan_borders, fill = "grey97", colour = "grey40", linewidth = 0.3) +
    geom_point(data = cbind(d, x = x), aes(long, lat, colour = x, shape = factor(pos)), size = 1.1) +
    scale_shape_manual(values = c(`0` = 1, `1` = 16), labels = c("no positives", "positives"), name = NULL) +
    scale_colour_viridis_c(name = NULL) + coord_sf(expand = FALSE) + theme_minimal(base_size = 9) +
    labs(x = NULL, y = NULL, title = "(a) village average")
  pb <- ggplot(bb, aes(x, zone, colour = stratum)) +
    geom_line(aes(linewidth = stratum == "all villages")) +
    geom_pointrange(data = subset(bb, stratum == "all villages"), aes(ymin = zl, ymax = zu), size = 0.25) +
    scale_linewidth_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.4), guide = "none") +
    scale_colour_manual(values = cols, name = NULL) + ylim(0, 1) + theme_minimal(base_size = 9) +
    labs(x = "village average (quintile mean)", y = "share of villages with positives", title = "(b) zone of transmission")
  pc <- ggplot(subset(bb, !is.na(prev)), aes(x, prev, colour = stratum)) +
    geom_line(aes(linewidth = stratum == "all villages")) +
    geom_pointrange(data = subset(bb, stratum == "all villages"), aes(ymin = pl, ymax = pu), size = 0.25) +
    scale_linewidth_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.4), guide = "none") +
    scale_colour_manual(values = cols, name = NULL) + theme_minimal(base_size = 9) +
    labs(x = "village average (quintile mean)", y = "prevalence among villages with positives", title = "(c) intensity")
  s <- subset(summ, covariate == k)
  sub <- sprintf("expected sign: %s (%s rationale).   Spearman, zone: %.2f (%.2f, %.2f);   intensity: %.2f (%.2f, %.2f)",
                 info$expected[info$covariate == k], info$rationale[info$covariate == k],
                 s$rho[1], s$lower[1], s$upper[1], s$rho[2], s$lower[2], s$upper[2])
  png(file.path(fig_dir, "eda", sprintf("eda_%s.png", k)), width = 13, height = 4.5, units = "in", res = 170)
  grid::grid.newpage()
  lay <- grid::grid.layout(2, 3, heights = grid::unit(c(0.6, 3.9), "in"), widths = grid::unit(c(4.6, 4.2, 4.2), "in"))
  grid::pushViewport(grid::viewport(layout = lay))
  grid::grid.text(lab, x = 0.01, y = 0.68, hjust = 0, gp = grid::gpar(fontsize = 12, fontface = "bold"),
                  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  grid::grid.text(sub, x = 0.01, y = 0.25, hjust = 0, gp = grid::gpar(fontsize = 9),
                  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  for (q in 1:3) print(list(pa, pb, pc)[[q]], vp = grid::viewport(layout.pos.row = 2, layout.pos.col = q))
  dev.off()
}
saveRDS(list(summary = summ, cor = cmat, xbar = X, ok = ok, villages = d), file.path(out_dir, "cov_eda.rds"))
