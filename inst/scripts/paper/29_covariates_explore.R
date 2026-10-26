###############################################################################
# Exploratory analysis of the RiverATLAS covariates (data(sudan_river_covariates))
# as candidate breeding covariates for the tail-down model
# (directed_taildown_joint.rds from 23_directed_profiles.R).
#
# step "describe": joins the covariates to the quadrature nodes, maps them
#   along the rivers, computes their correlations at the nodes that enter the
#   likelihood (within 40 km of a village), and plots the posterior mode of
#   the breeding field u of the tail-down model against each covariate. If a
#   covariate captures breeding habitat that log discharge misses, u should
#   vary with it.
# step "screen": adds each covariate (standardised) to the breeding part
#   log w = theta log Q + beta x + u, one at a time and then jointly, with the
#   hyperparameters (sigma_u, sigma_z, a_E, gamma, range) fixed at the
#   tail-down estimates, and reports the coefficient, its standard error and
#   the change in the Laplace log marginal likelihood. Since sigma_u is not
#   re-estimated, the gains are conservative.
# Usage: Rscript 29_covariates_explore.R <describe|screen>
###############################################################################
library(odispers)
library(sf)
library(Matrix)
library(ggplot2)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
step <- commandArgs(trailingOnly = TRUE)[1]
net <- load_result("inputs_full_network.rds"); inp <- load_result("inputs.rds")
nodes <- net$nodes
data(sudan_river_covariates); data(sudan_borders)
rc <- as.data.frame(sudan_river_covariates)
rc <- rc[match(nodes$HYRIV_ID, rc$HYRIV_ID), ]
stopifnot(!anyNA(rc$HYRIV_ID))

# GLiM lithology classes as numbered in HydroATLAS (11 = water bodies, which
# matches the lake reaches); crystalline basement and volcanic rocks are
# grouped as "hard rock", where rapids and stable substrate are common
hard <- c(2, 4, 7, 8, 9, 10, 12, 13)
cov <- data.frame(
  log_Q        = log(nodes$DIS_AV_CMS),
  log_slope    = log1p(rc$sgr_dk_rav),               # dm per km
  permanence   = pmin(rc$dis_m3_pmn / pmax(rc$dis_m3_pyr, 1e-6), 1),   # min month / mean
  inundation   = log1p(rc$inu_pc_ult),
  lakes        = log1p(rc$lka_pc_use),
  regulation   = log1p(rc$dor_pc_pva),
  forest       = rc$for_pc_use,
  elevation    = rc$ele_mt_cav / 100,
  runoff       = log1p(rc$run_mm_cyr),
  hard_rock    = as.numeric(rc$lit_cl_cmj %in% hard))
cov$log_Q2 <- cov$log_Q^2
cov$velocity <- 0.27 * cov$log_Q + 0.5 * log((rc$sgr_dk_rav + 0.5) / 1e4)  # log v, Manning + hydraulic geometry, up to a constant
labels <- c(log_Q = "log discharge", log_slope = "log(1 + slope, dm/km)", permanence = "flow permanence (min/mean)",
            inundation = "log(1 + inundation %)", lakes = "log(1 + lake %)", regulation = "log(1 + regulation %)",
            forest = "forest cover %", elevation = "elevation (100 m)", runoff = "log(1 + runoff mm)",
            hard_rock = "hard rock", log_Q2 = "log discharge squared", velocity = "log velocity proxy")

fit <- readRDS(file.path(out_dir, "directed_taildown_joint.rds"))
S <- fit$setup; used <- S$used
g <- directed_graph(river_graph(nodes, vertex_spacing = 5), "taildown")
S$A <- g$A[used, , drop = FALSE]
cand <- setdiff(names(cov), "log_Q")

if (step == "describe") {
  cu <- cov[used, ]
  u_used <- as.vector(S$A %*% fit$u)
  # weight of each node: kernel weight summed over the villages (how much the
  # data see that node). Weighting also by the fitted breeding weight would let
  # a few hotspots dominate (effective sample size about 40).
  infl <- colSums(S$K)
  cat("effective number of nodes under the kernel weights:", round(sum(infl)^2 / sum(infl^2)), "\n")
  cat("nodes in the likelihood:", length(used), "of", nrow(nodes), "\n")
  print(round(sapply(cu, function(x) c(mean = mean(x), sd = sd(x), q05 = quantile(x, 0.05), q95 = quantile(x, 0.95))), 3))
  cat("share of hard rock at used nodes:", round(mean(cu$hard_rock), 3), "\n")
  cat("lithology classes at used nodes:\n"); print(sort(table(rc$lit_cl_cmj[used]), decreasing = TRUE))

  # correlations among covariates (Spearman) at the used nodes
  cm <- cor(cu[, setdiff(names(cu), "log_Q2")], method = "spearman")
  print(round(cm, 2))
  cd <- as.data.frame(as.table(cm)); names(cd) <- c("a", "b", "r")
  p_cor <- ggplot(cd, aes(a, b, fill = r)) + geom_tile() + geom_text(aes(label = sprintf("%.2f", r)), size = 2.6) +
    scale_fill_gradient2(low = "#2166ac", high = "#b2182b", limits = c(-1, 1)) +
    theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x = NULL, y = NULL)
  ggsave(file.path(fig_dir, "cov_correlations.png"), p_cor, width = 7.5, height = 6.5, dpi = 200, bg = "white")

  # u against each covariate: influence-weighted binned means (deciles) with 95% intervals
  bins <- NULL; assoc <- NULL
  for (k in cand) {
    x <- cu[[k]]
    # deciles; for covariates with many zeros, a zero class plus quintiles of the positive values
    if (length(unique(x)) <= 2) b <- factor(x) else if (mean(x == 0) > 0.1) {
      br <- unique(quantile(x[x > 0], seq(0, 1, 0.2)))
      b <- factor(ifelse(x == 0, "0", as.character(cut(x, br, include.lowest = TRUE))))
    } else b <- cut(x, unique(quantile(x, seq(0, 1, 0.1))), include.lowest = TRUE)
    for (lv in levels(b)) {
      s <- b == lv; ww <- infl[s] / sum(infl[s])
      mu <- sum(ww * u_used[s]); neff <- 1 / sum(ww^2)
      se <- sqrt(sum(ww * (u_used[s] - mu)^2) / neff)
      bins <- rbind(bins, data.frame(covariate = labels[[k]], x = sum(ww * x[s]), u = mu, se = se))
    }
    assoc <- rbind(assoc, data.frame(covariate = k,
      spearman_all = cor(x, u_used, method = "spearman"),
      weighted_r = cov.wt(cbind(x, u_used), wt = infl / sum(infl), cor = TRUE)$cor[1, 2]))
  }
  print(assoc, digits = 3)
  p_u <- ggplot(bins, aes(x, u)) + geom_hline(yintercept = 0, colour = "grey60") +
    geom_pointrange(aes(ymin = u - 1.96 * se, ymax = u + 1.96 * se), size = 0.2, colour = "#1b7837") + geom_line(colour = "#1b7837") +
    facet_wrap(~ covariate, scales = "free_x", ncol = 4) + theme_minimal(base_size = 9) +
    labs(x = "covariate (decile mean)", y = "posterior mode of u (kernel-weighted mean)")
  ggsave(file.path(fig_dir, "cov_u_association.png"), p_u, width = 10, height = 6.5, dpi = 200, bg = "white")

  # maps along the rivers (every third node)
  ll <- st_coordinates(st_transform(st_as_sf(as.data.frame(nodes)[, c("x", "y")], coords = c("x", "y"),
                                             crs = attr(nodes, "crs")), 4326))
  sel <- seq(1, nrow(nodes), by = 3)
  mp <- function(k) {
    z <- cov[[k]][sel]; lim <- quantile(z, c(0.01, 0.99))
    data.frame(covariate = labels[[k]], lon = ll[sel, 1], lat = ll[sel, 2], value = (pmin(pmax(z, lim[1]), lim[2]) - lim[1]) / diff(lim))
  }
  md <- do.call(rbind, lapply(c("log_slope", "permanence", "inundation", "forest", "hard_rock", "velocity"), mp))
  md <- md[order(md$value), ]
  vil <- st_drop_geometry(inp$villages)
  p_map <- ggplot() + geom_sf(data = sudan_borders, fill = "grey97", colour = "grey40", linewidth = 0.3) +
    geom_point(data = md, aes(lon, lat, colour = value), size = 0.06) +
    geom_point(data = vil, aes(long, lat), shape = 1, size = 0.5, colour = "black", stroke = 0.2) +
    scale_colour_viridis_c(name = "relative value\n(1st-99th pct)") + facet_wrap(~ covariate, ncol = 3) +
    coord_sf(expand = FALSE) + theme_minimal(base_size = 9) + labs(x = NULL, y = NULL)
  ggsave(file.path(fig_dir, "cov_maps.png"), p_map, width = 11, height = 8, dpi = 200, bg = "white")
  saveRDS(list(assoc = assoc, bins = bins, cor = cm), file.path(out_dir, "cov_describe.rds"))
}

if (step == "screen") {
  R <- wm_precision(g, fit$hyper[["range_u"]], 1)
  fld <- list(field = "u", R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE),
              sig_u = fit$hyper[["sigma_u"]], sig_z = fit$hyper[["sigma_z"]],
              aE = fit$hyper[["a_E"]], gamma = fit$hyper[["gamma"]])
  cf <- fit$coefficients
  x0 <- list(beta = cf[1:2], c = cf[3:4], theta = cf[5], u = fit$u, z = fit$z)
  sc <- lapply(cov[used, cand], function(x) (x - mean(x)) / sd(x))
  run <- function(ks) {
    S2 <- S; S2$B <- cbind(S$B, do.call(cbind, sc[ks])); S2$m_theta <- c(S$m_theta, rep(0, length(ks)))
    x <- x0; x$theta <- c(x0$theta, rep(0, length(ks)))
    t0 <- Sys.time(); r <- odispers:::odisp_zil_laplace(fld, S2, x)
    q <- length(S2$m_theta); it <- 4 + seq_len(q)
    est <- r$x$theta; se <- sqrt(diag(r$Vfix)[it])
    cat(sprintf("%-40s logLik %.2f  (%.1f min)\n", paste(ks, collapse = "+"), r$logLik, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    list(ks = ks, logLik = r$logLik, logpost = r$logpost, theta = est, se = se)
  }
  base <- run(character(0))
  one <- lapply(cand, run); names(one) <- cand
  tab <- do.call(rbind, lapply(one, function(o) data.frame(covariate = o$ks, estimate = tail(o$theta, 1), se = tail(o$se, 1),
                                                         d_logLik = o$logLik - base$logLik, d_logpost = o$logpost - base$logpost)))
  tab$z <- tab$estimate / tab$se
  tab <- tab[order(-abs(tab$z)), ]
  print(tab, digits = 3, row.names = FALSE)
  # joint model with the habitat covariates (velocity and log_Q2 left out: velocity is a
  # linear combination of log Q and log slope)
  joint_ks <- c("log_slope", "permanence", "inundation", "lakes", "regulation", "forest", "elevation", "runoff", "hard_rock")
  joint <- run(joint_ks)
  jt <- data.frame(covariate = c("log_Q", joint_ks), estimate = joint$theta, se = joint$se)
  jt$z <- jt$estimate / jt$se
  cat("joint model: d_logLik =", round(joint$logLik - base$logLik, 2), "\n"); print(jt, digits = 3, row.names = FALSE)
  saveRDS(list(base = base, one = tab, joint = jt, joint_d_logLik = joint$logLik - base$logLik),
          file.path(out_dir, "cov_screen.rds"))
}
