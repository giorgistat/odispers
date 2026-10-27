# Newton refinement of the M3 maximum likelihood estimate with finite-difference Hessians in
# whitened coordinates (steps of half a standard error along the eigenvectors of the previous
# Hessian); inner fits start from the mode at the current point; evaluations on 2 cores.
suppressPackageStartupMessages(library(parallel))
txt <- readLines("/home/claude/odispers_repo/inst/scripts/hpc/taildown_covariates.R")
stop_at <- grep("^# ---- helpers: cross-validation", txt)
commandArgs <- function(trailingOnly = TRUE) c("/home/claude/hpc_ml_hess", "M3")
eval(parse(text = txt[1:(stop_at - 1)]))
st_file <- "/home/claude/hpc_ml_hess/newton_state.rds"
res <- readRDS("/home/claude/hpc_ml/fit_M3.rds"); f <- res$fit; S <- f$setup
nm <- c(names(f$coefficients), names(f$par)); d <- 11
Rcache <- list()
lap <- function(ph, x0) {
  cf <- ph[1:6]; ps <- ph[7:11]
  S2 <- S; S2$m_beta <- cf[1:2]; S2$m_c <- cf[3:4]; S2$m_theta <- cf[5:6]
  k <- sprintf("%.10f", ps[5])
  if (is.null(Rcache[[k]])) { R <- wm_precision(graph, exp(ps[5]), 1)
    Rcache[[k]] <<- list(R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE)) }
  fl <- list(field = "u", R = Rcache[[k]]$R, Rch = Rcache[[k]]$Rch, sig_u = exp(ps[1]), sig_z = exp(ps[2]),
             aE = ps[3], gamma = ps[4], range_u = exp(ps[5]))
  x <- list(beta = cf[1:2], c = cf[3:4], theta = cf[5:6], u = x0$u, z = x0$z)
  o <- tryCatch(odispers:::odisp_zil_laplace(fl, S2, x), error = function(e) list(logLik = -Inf))
  list(logLik = o$logLik, x = if (is.finite(o$logLik)) list(u = o$x$u, z = o$x$z) else NULL)
}
whess <- function(phi, xc, f0, A, h = 0.5) {
  pts <- list(); for (i in 1:d) { ei <- replace(numeric(d), i, h); pts <- c(pts, list(ei, -ei)) }
  pr <- list(); for (i in 1:(d - 1)) for (j in (i + 1):d) { ei <- replace(numeric(d), i, h); ej <- replace(numeric(d), j, h)
    pts <- c(pts, list(ei + ej, ei - ej, -ei + ej, -ei - ej)); pr <- c(pr, list(c(i, j))) }
  vals <- unlist(mclapply(pts, function(p) lap(phi + as.vector(A %*% p), xc)$logLik, mc.cores = 2, mc.preschedule = FALSE))
  fp <- vals[seq(1, 2 * d, 2)]; fm <- vals[seq(2, 2 * d, 2)]
  Hp <- diag((fp - 2 * f0 + fm) / h^2, d)
  for (k in seq_along(pr)) { v4 <- vals[2 * d + 4 * (k - 1) + 1:4]; i <- pr[[k]][1]; j <- pr[[k]][2]
    Hp[i, j] <- Hp[j, i] <- (v4[1] - v4[2] - v4[3] + v4[4]) / (4 * h^2) }
  list(Hp = Hp, gp = (fp - fm) / (2 * h))
}
whiten <- function(H) { e <- eigen(-H, symmetric = TRUE); e$vectors %*% diag(1 / sqrt(pmax(abs(e$values), 1e-2))) }
if (file.exists(st_file)) { st <- readRDS(st_file); log_msg("resuming at iteration %d, logLik %.4f", st$iter, st$f0) } else {
  h2 <- readRDS("/home/claude/hpc_ml_hess/hessian2_M3.rds")
  st <- list(iter = 0, phi = c(f$coefficients, f$par), x = list(u = f$u, z = f$z), f0 = h2$logLik,
             Hp = h2$Hpsi, gp = h2$gpsi, A = h2$A, history = data.frame(iter = 0, logLik = h2$logLik, pred_gain = NA))
}
names(st$phi) <- nm
for (it in (st$iter + 1):8) {
  e <- eigen(-st$Hp, symmetric = TRUE); lam <- pmax(e$values, 1e-3 * max(e$values))
  stp <- as.vector(e$vectors %*% (crossprod(e$vectors, st$gp) / lam))
  pred <- 0.5 * sum(st$gp * stp)
  log_msg("iteration %d: logLik %.4f, predicted gain %.4f", it, st$f0, pred)
  if (pred < 0.02) { log_msg("converged"); break }
  if (sqrt(sum(stp^2)) > 4) stp <- stp * 4 / sqrt(sum(stp^2))       # trust region, 4 standard errors
  ok <- FALSE
  for (t in c(1, 0.5, 0.25, 0.125)) {
    r <- lap(st$phi + as.vector(st$A %*% (t * stp)), st$x)
    log_msg("  step %.3f: logLik %.4f", t, r$logLik)
    if (is.finite(r$logLik) && r$logLik > st$f0 + 1e-3) { ok <- TRUE; break }
  }
  if (!ok) { log_msg("no ascent along the Newton direction; stopping"); break }
  st$phi <- st$phi + as.vector(st$A %*% (t * stp)); st$x <- r$x; st$f0 <- r$logLik
  # re-converge the inner mode at the new point from itself (centre of the next Hessian)
  r2 <- lap(st$phi, st$x); st$x <- r2$x; st$f0 <- r2$logLik
  H <- t(solve(st$A)) %*% st$Hp %*% solve(st$A); st$A <- whiten(H)
  t0 <- Sys.time(); hw <- whess(st$phi, st$x, st$f0, st$A)
  st$Hp <- hw$Hp; st$gp <- hw$gp; st$iter <- it
  st$history <- rbind(st$history, data.frame(iter = it, logLik = st$f0, pred_gain = pred))
  saveRDS(st, st_file)
  log_msg("  Hessian in %.1f min; whitened eigenvalues %s, max |g| %.3f", as.numeric(difftime(Sys.time(), t0, units = "mins")),
          paste(signif(range(eigen(-st$Hp, TRUE, TRUE)$values), 3), collapse = " to "), max(abs(st$gp)))
}
Ai <- solve(st$A); H <- t(Ai) %*% st$Hp %*% Ai; dimnames(H) <- list(nm, nm)
V <- solve(-H); se <- sqrt(pmax(diag(V), 0)); z <- qnorm(0.975); ph <- st$phi
logsc <- c(rep(FALSE, 6), TRUE, TRUE, FALSE, FALSE, TRUE)
tab <- data.frame(parameter = c(nm[1:6], "sigma_u", "sigma_z", "a_E", "gamma", "range_u (km)"),
                  estimate = ifelse(logsc, exp(ph), ph), se = ifelse(logsc, NA, se),
                  lower = ifelse(logsc, exp(ph - z * se), ph - z * se), upper = ifelse(logsc, exp(ph + z * se), ph + z * se))
print(tab, digits = 4); print(st$history)
log_msg("final logLik %.4f; eigenvalues of the observed information %s", st$f0,
        paste(signif(range(eigen(-H, TRUE, TRUE)$values), 3), collapse = " to "))
saveRDS(list(phi = ph, x = st$x, logLik = st$f0, H = H, vcov = V, intervals = tab, history = st$history),
        "/home/claude/hpc_ml_hess/newton_final_M3.rds")
log_msg("done")
