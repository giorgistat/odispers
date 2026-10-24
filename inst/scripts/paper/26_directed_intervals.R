###############################################################################
# 95% confidence intervals for the joint fits of the tail-down and tail-up
# models (23_directed_profiles.R, 24_tailup_joint_warm.R). Regression
# coefficients: Wald intervals from the Laplace approximation, conditional on
# the other parameters. sigma_u, sigma_z, a_E, gamma and the range: Wald
# intervals from the observed information of the approximate marginal
# likelihood, computed by central finite differences with the latent effects
# started at the mode of the fit, on the log scale for the standard
# deviations and the range.
# Usage: Rscript 26_directed_intervals.R <taildown|tailup>
###############################################################################
library(odispers)
library(Matrix)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))
type <- commandArgs(trailingOnly = TRUE)[1]
net <- load_result("inputs_full_network.rds")
g <- directed_graph(river_graph(net$nodes, vertex_spacing = 5), type)
fit <- readRDS(file.path(out_dir, sprintf("directed_%s_joint.rds", type)))
S <- fit$setup; S$A <- g$A[S$used, , drop = FALSE]
cf <- fit$coefficients
x_mode <- list(beta = cf[1:2], c = cf[3:4], theta = cf[5], u = fit$u, z = fit$z)
p0 <- c(log(fit$hyper[["sigma_u"]]), log(fit$hyper[["sigma_z"]]), fit$hyper[["a_E"]], fit$hyper[["gamma"]],
        log(fit$hyper[["range_u"]]))
names(p0) <- c("log_sigma_u", "log_sigma_z", "a_E", "gamma", "log_range")
ll_at <- function(par) {
  R <- wm_precision(g, exp(par[5]), 1)
  fld <- list(field = "u", R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE),
              sig_u = exp(par[1]), sig_z = exp(par[2]), aE = par[3], gamma = par[4])
  r <- odispers:::odisp_zil_laplace(fld, S, x_mode)
  cat(sprintf("  %s  logLik = %.4f\n", paste(sprintf("%.4f", par), collapse = " "), r$logLik))
  r$logLik
}
k <- length(p0); h <- c(0.05, 0.05, 0.02, 0.02, 0.05)
f0 <- ll_at(p0)
cat(sprintf("logLik at the estimate %.4f (fit: %.4f)\n", f0, fit$logLik))
H <- matrix(0, k, k); fp <- fm <- numeric(k)
for (i in 1:k) {
  e <- replace(numeric(k), i, h[i]); fp[i] <- ll_at(p0 + e); fm[i] <- ll_at(p0 - e)
  H[i, i] <- (fp[i] - 2 * f0 + fm[i]) / h[i]^2
}
for (i in 1:(k - 1)) for (j in (i + 1):k) {
  ei <- replace(numeric(k), i, h[i]); ej <- replace(numeric(k), j, h[j])
  H[i, j] <- H[j, i] <- (ll_at(p0 + ei + ej) - ll_at(p0 + ei - ej) - ll_at(p0 - ei + ej) + ll_at(p0 - ei - ej)) /
    (4 * h[i] * h[j])
}
cat("gradient at the estimate:", round((fp - fm) / (2 * h), 3), "\n")
ev <- eigen(-H, symmetric = TRUE)$values; cat("eigenvalues of the observed information:", signif(ev, 3), "\n")
V <- solve(-H); se <- sqrt(pmax(diag(V), 0)); z <- stats::qnorm(0.975)
logsc <- c(TRUE, TRUE, FALSE, FALSE, TRUE)
hyp <- data.frame(parameter = c("sigma_u", "sigma_z", "a_E", "gamma", "range_u"),
                  estimate = ifelse(logsc, exp(p0), p0),
                  lower = ifelse(logsc, exp(p0 - z * se), p0 - z * se),
                  upper = ifelse(logsc, exp(p0 + z * se), p0 + z * se))
sec <- sqrt(diag(fit$vcov))
reg <- data.frame(parameter = names(cf), estimate = cf, lower = cf - z * sec, upper = cf + z * sec)
tab <- rbind(reg, hyp); rownames(tab) <- NULL
print(tab, digits = 3)
saveRDS(list(table = tab, hessian = H, vcov_hyper = V), file.path(out_dir, sprintf("directed_%s_intervals.rds", type)))
