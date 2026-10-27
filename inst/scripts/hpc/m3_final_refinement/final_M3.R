# Final M3: observed information at the BFGS estimate (finite differences along the eigenvectors
# of the BFGS approximation, half a standard error, two cores), Wald intervals, then the
# cross-validation and grid predictions of taildown_covariates.R at the new estimate.
suppressPackageStartupMessages(library(parallel))
txt <- readLines("/home/claude/odispers_repo/inst/scripts/hpc/taildown_covariates.R")
stop_at <- grep("^# ---- run the models", txt)
commandArgs <- function(trailingOnly = TRUE) c("/home/claude/hpc_ml_hess", "M3")
eval(parse(text = txt[1:(stop_at - 1)])); eval(parse(text = txt[stop_at + 1]))   # slim()
res <- readRDS("/home/claude/hpc_ml/fit_M3.rds"); f <- res$fit; S <- f$setup
b <- readRDS("/home/claude/hpc_ml_hess/bfgs_state.rds")
nm <- c(names(f$coefficients), names(f$par)); d <- 11
Rcache <- list()
mkfld <- function(ps) { k <- sprintf("%.10f", ps[5])
  if (is.null(Rcache[[k]])) { R <- wm_precision(graph, exp(ps[5]), 1)
    Rcache[[k]] <<- list(R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE)) }
  list(field = "u", R = Rcache[[k]]$R, Rch = Rcache[[k]]$Rch, sig_u = exp(ps[1]), sig_z = exp(ps[2]),
       aE = ps[3], gamma = ps[4], range_u = exp(ps[5])) }
mkS <- function(cf) { S2 <- S; S2$m_beta <- cf[1:2]; S2$m_c <- cf[3:4]; S2$m_theta <- cf[5:6]; S2 }
lap <- function(ph, x0) { cf <- ph[1:6]; ps <- ph[7:11]
  x <- list(beta = cf[1:2], c = cf[3:4], theta = cf[5:6], u = x0$u, z = x0$z)
  tryCatch(odispers:::odisp_zil_laplace(mkfld(ps), mkS(cf), x), error = function(e) list(logLik = -Inf)) }
phi <- b$phi; names(phi) <- nm
r0 <- lap(phi, b$x); x0 <- list(u = r0$x$u, z = r0$x$z); f0 <- r0$logLik
log_msg("logLik at the estimate %.4f", f0)
# scaling: minus the Hessian approximated by A^-T B A^-1
Hb <- -t(solve(b$A)) %*% b$B %*% solve(b$A)
e <- eigen(-Hb, symmetric = TRUE); A <- e$vectors %*% diag(1 / sqrt(pmax(e$values, 1e-2)))
h <- 0.5
pts <- list(); for (i in 1:d) { ei <- replace(numeric(d), i, h); pts <- c(pts, list(ei, -ei)) }
pr <- list(); for (i in 1:(d - 1)) for (j in (i + 1):d) { ei <- replace(numeric(d), i, h); ej <- replace(numeric(d), j, h)
  pts <- c(pts, list(ei + ej, ei - ej, -ei + ej, -ei - ej)); pr <- c(pr, list(c(i, j))) }
t0 <- Sys.time()
vals <- unlist(mclapply(pts, function(p) lap(phi + as.vector(A %*% p), x0)$logLik, mc.cores = 2, mc.preschedule = FALSE))
fp <- vals[seq(1, 2 * d, 2)]; fm <- vals[seq(2, 2 * d, 2)]
Hp <- diag((fp - 2 * f0 + fm) / h^2, d)
for (k in seq_along(pr)) { v4 <- vals[2 * d + 4 * (k - 1) + 1:4]; i <- pr[[k]][1]; j <- pr[[k]][2]
  Hp[i, j] <- Hp[j, i] <- (v4[1] - v4[2] - v4[3] + v4[4]) / (4 * h^2) }
gp <- (fp - fm) / (2 * h); Ai <- solve(A)
H <- t(Ai) %*% Hp %*% Ai; dimnames(H) <- list(nm, nm)
ev <- eigen(-H, TRUE, TRUE)$values
st <- solve(-Hp, gp)
log_msg("Hessian in %.1f min; eigenvalues (whitened) %s; max |g| (whitened) %.3f; Newton predicted gain %.4f",
        as.numeric(difftime(Sys.time(), t0, units = "mins")), paste(signif(range(eigen(-Hp, TRUE, TRUE)$values), 3), collapse = " to "),
        max(abs(gp)), 0.5 * sum(gp * st))
V <- solve(-H); se <- sqrt(pmax(diag(V), 0)); z <- qnorm(0.975)
logsc <- c(rep(FALSE, 6), TRUE, TRUE, FALSE, FALSE, TRUE)
tab <- data.frame(parameter = c(nm[1:6], "sigma_u", "sigma_z", "a_E", "gamma", "range_u (km)"),
                  estimate = ifelse(logsc, exp(phi), phi), se = ifelse(logsc, NA, se),
                  lower = ifelse(logsc, exp(phi - z * se), phi - z * se), upper = ifelse(logsc, exp(phi + z * se), phi + z * se))
print(tab, digits = 4)
# updated fit object
cf <- phi[1:6]; ps <- phi[7:11]; fl <- mkfld(ps)
fit <- f; fit$coefficients <- cf; fit$par <- ps; fit$vcov <- V; fit$intervals <- tab; fit$eigen <- ev
fit$gradient <- as.vector(crossprod(Ai, gp)); fit$logLik <- f0; fit$u <- r0$x$u; fit$z <- r0$x$z
fit$hyper <- c(sigma_u = exp(ps[[1]]), sigma_z = exp(ps[[2]]), a_E = ps[[3]], gamma = ps[[4]], range_u = exp(ps[[5]]))
fit$setup <- mkS(cf); fit$fld <- fl
fit$E <- exp(r0$a); fit$pi <- stats::plogis(r0$d); fit$p <- -expm1(-exp(r0$e))
out <- list(fit = fit)
saveRDS(list(fit = slim(fit)), file.path(out_dir, "fit_M3.rds"))
write.csv(tab, file.path(out_dir, "intervals_M3.csv"), row.names = FALSE)
t0 <- Sys.time(); out$cv <- cv_zil(fit, folds); out$cv_scores <- cv_scores(out$cv, v$AI[fit$setup$has])
log_msg("cross-validation in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins")))
print(out$cv_scores$table, digits = 4)
pg <- predict(fit, grid, nodes, graph, pairs = grid_pairs); xy <- st_coordinates(grid)
out$grid <- data.frame(lon = xy[, 1], lat = xy[, 2], AI = grid$AI, pg)
out$fit <- slim(fit)
saveRDS(out, file.path(out_dir, "fit_M3.rds"))
log_msg("done")
