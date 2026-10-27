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

st <- readRDS(st_file)   # state of the Newton iterations: point, mode, whitened Hessian
log_msg("start: logLik %.4f", st$f0)
bf_file <- "/home/claude/hpc_ml_hess/bfgs_state.rds"
phi <- st$phi; x <- st$x; f0 <- st$f0
A <- whiten(t(solve(st$A)) %*% st$Hp %*% solve(st$A))          # scaling from the latest Hessian
B <- diag(d); hist <- data.frame()
if (file.exists(bf_file)) { b <- readRDS(bf_file); phi <- b$phi; x <- b$x; f0 <- b$f0; A <- b$A; B <- b$B; hist <- b$hist
  log_msg("resuming BFGS at logLik %.4f", f0) }
grad <- function(phi, x, h = 0.5) {
  pts <- list(); for (i in 1:d) { e <- replace(numeric(d), i, h); pts <- c(pts, list(e, -e)) }
  v <- unlist(mclapply(pts, function(p) lap(phi + as.vector(A %*% p), x)$logLik, mc.cores = 2, mc.preschedule = FALSE))
  (v[seq(1, 2 * d, 2)] - v[seq(2, 2 * d, 2)]) / (2 * h)
}
g <- grad(phi, x)
for (it in 1:40) {
  p <- as.vector(solve(B, g)); pred <- 0.5 * sum(g * p)
  log_msg("BFGS %d: logLik %.4f, |g| %.3f, predicted gain %.4f", it, f0, sqrt(sum(g^2)), pred)
  hist <- rbind(hist, data.frame(it = it, logLik = f0, grad = sqrt(sum(g^2)), pred = pred, t(phi)))
  if (pred < 0.01 && max(abs(g)) < 0.1) { log_msg("converged"); break }
  if (sqrt(sum(p^2)) > 2) p <- p * 2 / sqrt(sum(p^2))
  ts <- c(1, 0.5, 0.25, 0.1)
  rs <- mclapply(ts, function(t) lap(phi + as.vector(A %*% (t * p)), x), mc.cores = 2, mc.preschedule = FALSE)
  fs <- sapply(rs, `[[`, "logLik"); k <- which.max(fs)
  log_msg("  line search: %s", paste(sprintf("t=%.2f: %.4f", ts, fs), collapse = ", "))
  if (!(fs[k] > f0 + 1e-4)) {
    if (!isTRUE(all.equal(B, diag(d)))) { log_msg("  no ascent: reset BFGS matrix"); B <- diag(d); next }
    log_msg("no ascent along the gradient; stopping"); break }
  s <- ts[k] * p; phi <- phi + as.vector(A %*% s); x <- rs[[k]]$x
  r2 <- lap(phi, x); x <- r2$x; f0 <- r2$logLik
  g_new <- grad(phi, x); yv <- -(g_new - g)                    # B approximates minus the Hessian
  if (sum(s * yv) > 1e-8) B <- B - (B %*% s %*% t(s) %*% B) / as.numeric(t(s) %*% B %*% s) + (yv %*% t(yv)) / sum(yv * s)
  g <- g_new
  saveRDS(list(phi = phi, x = x, f0 = f0, A = A, B = B, hist = hist), bf_file)
}
saveRDS(list(phi = phi, x = x, f0 = f0, A = A, B = B, hist = hist, g = g), bf_file)
print(round(phi, 4)); log_msg("done")
