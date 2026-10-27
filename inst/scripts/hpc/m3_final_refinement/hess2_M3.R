# Observed information of M3 at the estimate in whitened coordinates: finite differences along
# the eigenvectors of a first Hessian, each scaled to about half a standard error, so that every
# step changes the log-likelihood by a similar amount; all inner fits start from the mode at the
# estimate. The Hessian is then transformed back to the original parameters.
txt <- readLines("/home/claude/odispers_repo/inst/scripts/hpc/taildown_covariates.R")
stop_at <- grep("^# ---- helpers: cross-validation", txt)
commandArgs <- function(trailingOnly = TRUE) c("/home/claude/hpc_ml_hess", "M3")
eval(parse(text = txt[1:(stop_at - 1)]))
res <- readRDS("/home/claude/hpc_ml/fit_M3.rds"); f <- res$fit
H1 <- readRDS("/home/claude/hpc_ml_hess/hessian_M3.rds")
S <- f$setup; mode <- list(u = f$u, z = f$z)
phi0 <- c(f$coefficients, f$par); d <- length(phi0)
Rcache <- list()
ll <- function(ph) {
  cf <- ph[1:6]; ps <- ph[7:11]
  S2 <- S; S2$m_beta <- cf[1:2]; S2$m_c <- cf[3:4]; S2$m_theta <- cf[5:6]
  k <- sprintf("%.10f", ps[5])
  if (is.null(Rcache[[k]])) { R <- wm_precision(graph, exp(ps[5]), 1)
    Rcache[[k]] <<- list(R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE)) }
  fl <- list(field = "u", R = Rcache[[k]]$R, Rch = Rcache[[k]]$Rch, sig_u = exp(ps[1]), sig_z = exp(ps[2]),
             aE = ps[3], gamma = ps[4], range_u = exp(ps[5]))
  x <- list(beta = cf[1:2], c = cf[3:4], theta = cf[5:6], u = mode$u, z = mode$z)
  odispers:::odisp_zil_laplace(fl, S2, x)$logLik
}
f0 <- ll(phi0); log_msg("logLik at the estimate %.4f", f0)
Hold <- -solve(H1$vcov)
e <- eigen(-Hold, symmetric = TRUE)
A <- e$vectors %*% diag(1 / sqrt(pmax(abs(e$values), 1e-2)))       # phi = phi0 + A psi
h <- 0.5
fp <- fm <- numeric(d); Hp <- matrix(0, d, d); n_ev <- 0
lp <- function(psi) { n_ev <<- n_ev + 1; if (n_ev %% 20 == 0) log_msg("%d evaluations", n_ev); ll(phi0 + as.vector(A %*% psi)) }
for (i in 1:d) { ei <- replace(numeric(d), i, h); fp[i] <- lp(ei); fm[i] <- lp(-ei); Hp[i, i] <- (fp[i] - 2 * f0 + fm[i]) / h^2 }
for (i in 1:(d - 1)) for (j in (i + 1):d) {
  ei <- replace(numeric(d), i, h); ej <- replace(numeric(d), j, h)
  Hp[i, j] <- Hp[j, i] <- (lp(ei + ej) - lp(ei - ej) - lp(-ei + ej) + lp(-ei - ej)) / (4 * h^2)
}
gp <- (fp - fm) / (2 * h)
Ai <- solve(A); H <- t(Ai) %*% Hp %*% Ai; g <- as.vector(t(Ai) %*% gp)
dimnames(H) <- list(names(phi0), names(phi0))
ev <- eigen(-H, symmetric = TRUE, only.values = TRUE)$values
log_msg("whitened Hessian: eigenvalues %s; gradient (whitened) max |g| %.4f", paste(signif(range(eigen(-Hp, TRUE, TRUE)$values), 3), collapse = " to "), max(abs(gp)))
log_msg("eigenvalues of the observed information %s", paste(signif(range(ev), 3), collapse = " to "))
V <- solve(-H); se <- sqrt(pmax(diag(V), 0)); z <- qnorm(0.975)
logsc <- c(rep(FALSE, 6), TRUE, TRUE, FALSE, FALSE, TRUE)
tab <- data.frame(parameter = c(names(f$coefficients), "sigma_u", "sigma_z", "a_E", "gamma", "range_u (km)"),
                  estimate = ifelse(logsc, exp(phi0), phi0), se = ifelse(logsc, NA, se),
                  lower = ifelse(logsc, exp(phi0 - z * se), phi0 - z * se), upper = ifelse(logsc, exp(phi0 + z * se), phi0 + z * se),
                  se_log = ifelse(logsc, se, NA))
print(tab, digits = 4)
saveRDS(list(intervals = tab, vcov = V, H = H, eigen = ev, gradient = g, Hpsi = Hp, gpsi = gp, A = A, logLik = f0),
        "/home/claude/hpc_ml_hess/hessian2_M3.rds")
write.csv(tab, "/home/claude/hpc_ml_hess/intervals_M3.csv", row.names = FALSE)
log_msg("done")
