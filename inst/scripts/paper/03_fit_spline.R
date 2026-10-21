###############################################################################
# Main model: range fixed (40 km, and 53 km as sensitivity), gamma = 0,
# breeding quadratic in log discharge, aridity in the intensity as a linear
# spline in log AI whose change point is estimated by profile likelihood.
# Writes paper_results/spline_fits.rds.
###############################################################################

library(odispers)
library(sf)
source("00_settings.R")

inp <- readRDS(file.path(out_dir, "inputs.rds"))
v <- inp$villages; nodes <- inp$nodes; pairs <- inp$pairs
la <- log(v$AI)

fit_cp <- function(cp, phi) {
  v$h <- pmax(la - cp, 0)
  suppressWarnings(odisp_fit(y ~ log(AI) + h, breeding, data = v, nodes = nodes, pairs = pairs,
                             max_range = max_range, fixed = list(gamma = 0, eta0 = log(phi))))
}

out <- list()
for (phi in c(phi_main, phi_sens)) {
  grid_cp <- seq(quantile(la, 0.05), quantile(la, 0.95), length.out = 31)
  ll <- sapply(grid_cp, function(cp) fit_cp(cp, phi)$logLik)
  b <- which.max(ll)                                     # refine around the best grid point
  fine <- seq(grid_cp[max(1, b - 1)], grid_cp[min(length(grid_cp), b + 1)], length.out = 21)
  llf <- sapply(fine, function(cp) fit_cp(cp, phi)$logLik)
  cp_hat <- fine[which.max(llf)]
  f <- fit_cp(cp_hat, phi)
  prof <- data.frame(cp = c(grid_cp, fine), logLik = c(ll, llf)); prof <- prof[order(prof$cp), ]
  ci <- range(prof$cp[prof$logLik > max(prof$logLik) - qchisq(0.95, 1) / 2])
  f_lin <- suppressWarnings(odisp_fit(y ~ log(AI), breeding, data = v, nodes = nodes, pairs = pairs,
                                      max_range = max_range, fixed = list(gamma = 0, eta0 = log(phi))))

  cat(sprintf("\nphi = %d km: change point AI = %.3f (95%% profile interval %.3f-%.3f)\n",
              phi, exp(cp_hat), exp(ci[1]), exp(ci[2])))
  cat(sprintf("logLik: linear in log AI %.2f, spline %.2f\n", f_lin$logLik, f$logLik))
  print(summary(f)$table[, 1:2])
  cf <- coef(f); V <- vcov(f)
  above <- cf["beta.log(AI)"] + cf["beta.h"]
  se_above <- sqrt(V["beta.log(AI)", "beta.log(AI)"] + V["beta.h", "beta.h"] + 2 * V["beta.log(AI)", "beta.h"])
  cat(sprintf("slope below change point %.2f (%.2f); above %.2f (%.2f)\n", cf["beta.log(AI)"],
              sqrt(V["beta.log(AI)", "beta.log(AI)"]), above, se_above))
  cls <- aridity_classes(v$AI)
  print(round(cbind(observed = tapply(v$y, cls, sum) / tapply(v$n, cls, sum),
                    fitted = tapply(f$fitted * v$n, cls, sum) / tapply(v$n, cls, sum)), 3))
  out[[as.character(phi)]] <- list(fit = f, fit_lin = f_lin, cp = cp_hat, ci = ci, prof = prof)
}
saveRDS(out, file.path(out_dir, "spline_fits.rds"))
