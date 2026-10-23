###############################################################################
# 95% confidence intervals for the parameters of the zero-inflated model of
# 17_fit_zil.R. Regression coefficients (a0, a1, b0, b1, theta): Wald
# intervals from the Laplace approximation at the posterior mode, conditional
# on the other parameters. sigma_u, sigma_z, a_E and gamma: Wald intervals
# from the observed information of the approximate marginal likelihood,
# computed by central finite differences, on the log scale for the standard
# deviations (so that their intervals are positive).
###############################################################################

library(odispers)
library(sf)
library(Matrix)
source(system.file("scripts/paper/00_settings.R", package = "odispers"))   # shared settings

fit <- load_result("zil_fit.rds")$fit
inp <- load_result("inputs.rds"); v <- inp$villages
net <- load_result("inputs_full_network.rds")
nodes <- net$nodes
graph <- river_graph(nodes, vertex_spacing = 5)
S <- fit$setup
R <- wm_precision(graph, fit$hyper[["range_u"]], 1)
base <- list(field = "u", R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE))
x_mode <- list(beta = S$m_beta, c = S$m_c, theta = S$m_theta, u = fit$u, z = fit$z)
x_mode$beta[] <- fit$coefficients[grep("^lambda\\.", names(fit$coefficients))]
x_mode$c[] <- fit$coefficients[grep("^pi\\.", names(fit$coefficients))]
x_mode$theta[] <- fit$coefficients[grep("^theta\\.", names(fit$coefficients))]

ll_at <- function(par) {   # par = (log sigma_u, log sigma_z, a_E, gamma)
  fld <- c(base, list(sig_u = exp(par[1]), sig_z = exp(par[2]), aE = par[3], gamma = par[4]))
  r <- odispers:::odisp_zil_laplace(fld, S, x_mode)
  cat(sprintf("  %s  logLik = %.4f\n", paste(sprintf("%.4f", par), collapse = " "), r$logLik))
  r$logLik
}
p0 <- fit$par; k <- length(p0)
h <- c(0.05, 0.05, 0.02, 0.02)
f0 <- ll_at(p0)
cat(sprintf("logLik at the estimate %.4f (fit: %.4f)\n", f0, fit$logLik))
H <- matrix(0, k, k)
fp <- fm <- numeric(k)
for (i in 1:k) {
  e <- replace(numeric(k), i, h[i])
  fp[i] <- ll_at(p0 + e); fm[i] <- ll_at(p0 - e)
  H[i, i] <- (fp[i] - 2 * f0 + fm[i]) / h[i]^2
}
for (i in 1:(k - 1)) for (j in (i + 1):k) {
  ei <- replace(numeric(k), i, h[i]); ej <- replace(numeric(k), j, h[j])
  H[i, j] <- H[j, i] <- (ll_at(p0 + ei + ej) - ll_at(p0 + ei - ej) - ll_at(p0 - ei + ej) + ll_at(p0 - ei - ej)) /
    (4 * h[i] * h[j])
}
cat("Gradient at the estimate (central differences):", round((fp - fm) / (2 * h), 3), "\n")
cat("Eigenvalues of -H:", signif(eigen(-H)$values, 3), "\n")
V <- solve(-H); se <- sqrt(diag(V))
z <- stats::qnorm(0.975)
hyp <- data.frame(parameter = c("sigma_u", "sigma_z", "a_E", "gamma"),
                  estimate = c(exp(p0[1:2]), p0[3:4]),
                  se_working = se, scale = c("log", "log", "identity", "identity"),
                  lower = c(exp(p0[1:2] - z * se[1:2]), p0[3:4] - z * se[3:4]),
                  upper = c(exp(p0[1:2] + z * se[1:2]), p0[3:4] + z * se[3:4]))
cf <- fit$coefficients; sec <- sqrt(diag(fit$vcov))
reg <- data.frame(parameter = names(cf), estimate = cf, se_working = sec, scale = "identity",
                  lower = cf - z * sec, upper = cf + z * sec)
tab <- rbind(reg, hyp); rownames(tab) <- NULL
print(tab, digits = 3)
saveRDS(list(table = tab, hessian = H, vcov_hyper = V), file.path(out_dir, "zil_intervals.rds"))
