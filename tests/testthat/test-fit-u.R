toy_u_setup <- function(m = 150, seed = 2) {
  set.seed(seed)
  # three reaches drawn from upstream to downstream; A and B join C
  rivers <- sf::st_sf(
    HYRIV_ID = 1:3, NEXT_DOWN = c(3, 3, 0), DIS_AV_CMS = c(5, 0.5, 50),
    geometry = sf::st_sfc(
      sf::st_linestring(rbind(c(30.5, 2.0), c(30.2, 1.2))),
      sf::st_linestring(rbind(c(29.2, 1.4), c(29.6, 1.5), c(30.2, 1.2))),
      sf::st_linestring(rbind(c(30.2, 1.2), c(30.3, 0.6), c(30.0, 0.0))),
      crs = 4326))
  lon <- runif(m, 29.3, 30.8); lat <- runif(m, 0.1, 1.9)
  v <- sf::st_as_sf(data.frame(lon = lon, lat = lat, n = 50L, AI = exp(-2 + 0.8 * lat)),
                    coords = c("lon", "lat"), crs = 4326)
  nodes <- river_nodes(rivers, spacing = 2)
  pairs <- river_pairs(v, nodes, max_range = 40)
  # simulate from the model without u, with a village effect
  d <- odisp_design(n ~ log(AI), ~ log(DIS_AV_CMS), v, nodes, pairs, max_range = 40)
  ev <- odispers:::odisp_eval(c(-3, 0.5, 0.3, log(30), 0), d, FALSE)
  lam <- ev$lambda * exp(rnorm(m, 0, 0.5))
  v$y <- stats::rbinom(m, v$n, 1 - exp(-lam))
  list(v = v, nodes = nodes, pairs = pairs)
}

test_that("the Laplace mode is a stationary point of the log posterior", {
  tt <- toy_u_setup()
  v <- tt$v[sapply(seq_len(nrow(tt$v)), function(i) any(tt$pairs$i == i & tt$pairs$r < 30)), ]
  pairs <- river_pairs(v, tt$nodes, 40)
  g <- river_graph(tt$nodes, vertex_spacing = 5)
  S <- odispers:::odisp_u_setup(y ~ log(AI), ~ log(DIS_AV_CMS), v, tt$nodes, g, pairs,
                                30, "n", "AI", 10)
  R <- wm_precision(g, 50); Rch <- Matrix::Cholesky(R, LDL = FALSE, perm = TRUE)
  fld <- list(field = "u", R = R, Rch = Rch, sig_u = 0.7, sig_z = 0.4)
  x0 <- list(beta = S$m_beta, theta = S$m_theta, u = rep(0, g$nv), s = NULL, z = rep(0, length(S$y)))
  res <- odispers:::odisp_u_laplace(fld, S, x0)
  expect_true(is.finite(res$logLik))
  f <- function(v) {
    x <- odispers:::odisp_u_unvec(v, x0)
    odispers:::odisp_u_lik(odispers:::odisp_u_eta(x, S)$eta, S)$ll - 0.5 * odispers:::odisp_u_prior(x, S, fld)
  }
  xv <- odispers:::odisp_u_vec(res$x)
  idx <- c(1:4, 4 + c(1, 5, 10), length(xv) - 0:2)   # beta, theta, some u, some z
  h <- 1e-5
  num <- sapply(idx, function(k) { e <- replace(numeric(length(xv)), k, h); (f(xv + e) - f(xv - e)) / (2 * h) })
  expect_lt(max(abs(num)), 1e-3)
})

test_that("odisp_fit_u runs and the field improves on the village effect alone", {
  tt <- toy_u_setup()
  v <- tt$v[sapply(seq_len(nrow(tt$v)), function(i) any(tt$pairs$i == i & tt$pairs$r < 30)), ]
  pairs <- river_pairs(v, tt$nodes, 40)
  g <- river_graph(tt$nodes, vertex_spacing = 5)
  fz <- odisp_fit_u(y ~ log(AI), ~ log(DIS_AV_CMS), v, tt$nodes, g, pairs, phi = 30, u = FALSE, trace = FALSE)
  fu <- odisp_fit_u(y ~ log(AI), ~ log(DIS_AV_CMS), v, tt$nodes, g, pairs, phi = 30, range_u = 50, trace = FALSE)
  expect_s3_class(fu, "odisp_ufit")
  expect_true(is.finite(fu$logLik) && is.finite(fz$logLik))
  expect_gte(fu$logLik, fz$logLik - 0.5)    # u can shrink to zero
  pr <- predict(fu, v, tt$nodes, g, pairs = pairs)
  expect_equal(nrow(pr), nrow(v))
  expect_true(all(pr$prevalence >= 0 & pr$prevalence <= 1))
})

test_that("the Laplace mode is stationary with a Gaussian process S and no village effect", {
  tt <- toy_u_setup()
  v <- tt$v[sapply(seq_len(nrow(tt$v)), function(i) any(tt$pairs$i == i & tt$pairs$r < 30)), ]
  pairs <- river_pairs(v, tt$nodes, 40)
  g <- river_graph(tt$nodes, vertex_spacing = 5)
  S <- odispers:::odisp_u_setup(y ~ log(AI), ~ log(DIS_AV_CMS), v, tt$nodes, g, pairs,
                                30, "n", "AI", 10)
  n <- length(S$y)
  fld <- list(field = "S", LC = chol(odispers:::odisp_s_corr(S$xy, 60) + diag(1e-8, n)), sig_s = 0.8)
  x0 <- list(beta = S$m_beta, theta = S$m_theta, u = NULL, s = rep(0, n), z = NULL)
  res <- odispers:::odisp_u_laplace(fld, S, x0)
  expect_true(is.finite(res$logLik))
  f <- function(v) {
    x <- odispers:::odisp_u_unvec(v, x0)
    odispers:::odisp_u_lik(odispers:::odisp_u_eta(x, S)$eta, S)$ll - 0.5 * odispers:::odisp_u_prior(x, S, fld)
  }
  xv <- odispers:::odisp_u_vec(res$x)
  idx <- c(1:3, 3 + c(1, 7, 20), length(xv) - 0:1)
  h <- 1e-5
  num <- sapply(idx, function(k) { e <- replace(numeric(length(xv)), k, h); (f(xv + e) - f(xv - e)) / (2 * h) })
  expect_lt(max(abs(num)), 1e-3)
  fs <- odisp_fit_u(y ~ log(AI), ~ log(DIS_AV_CMS), v, tt$nodes, g, pairs, phi = 30,
                    field = "S", range = 60, z = FALSE, trace = FALSE)
  expect_equal(fs$field, "S")
  pr <- predict(fs, v, tt$nodes, g, pairs = pairs)
  expect_true(all(is.finite(pr$S)))
})

test_that("u and S together: stationary mode and a fit", {
  tt <- toy_u_setup()
  v <- tt$v[sapply(seq_len(nrow(tt$v)), function(i) any(tt$pairs$i == i & tt$pairs$r < 30)), ]
  pairs <- river_pairs(v, tt$nodes, 40)
  g <- river_graph(tt$nodes, vertex_spacing = 5)
  S <- odispers:::odisp_u_setup(y ~ log(AI), ~ log(DIS_AV_CMS), v, tt$nodes, g, pairs,
                                30, "n", "AI", 10)
  n <- length(S$y)
  R <- wm_precision(g, 50)
  fld <- list(field = "uS", R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE),
              LC = chol(odispers:::odisp_s_corr(S$xy, 60) + diag(1e-8, n)),
              sig_u = 0.6, sig_s = 0.5, sig_z = 0.3)
  x0 <- list(beta = S$m_beta, theta = S$m_theta, u = rep(0, g$nv), s = rep(0, n), z = rep(0, n))
  res <- odispers:::odisp_u_laplace(fld, S, x0)
  f <- function(v) {
    x <- odispers:::odisp_u_unvec(v, x0)
    odispers:::odisp_u_lik(odispers:::odisp_u_eta(x, S)$eta, S)$ll - 0.5 * odispers:::odisp_u_prior(x, S, fld)
  }
  xv <- odispers:::odisp_u_vec(res$x)
  nu <- g$nv
  idx <- c(1:3, 3 + c(1, 5), 3 + nu + c(1, 9), length(xv) - 0:1)
  h <- 1e-5
  num <- sapply(idx, function(k) { e <- replace(numeric(length(xv)), k, h); (f(xv + e) - f(xv - e)) / (2 * h) })
  # the analytic gradient agrees with the numerical one, and is small at the mode
  ev <- odispers:::odisp_u_eta(res$x, S); lk <- odispers:::odisp_u_lik(ev$eta, S)
  ops <- odispers:::odisp_u_ops(res$x, ev, lk, S, fld)
  m <- odispers:::odisp_u_vec(list(beta = S$m_beta, theta = S$m_theta, u = 0 * res$x$u, s = 0 * res$x$s, z = 0 * res$x$z))
  ga <- ops$Jtw(lk$g) - ops$Pv(xv - m)
  expect_lt(max(abs(ga[idx] - num)), 1e-5)
  expect_lt(max(abs(num)), 5e-3)
  fus <- odisp_fit_u(y ~ log(AI), ~ log(DIS_AV_CMS), v, tt$nodes, g, pairs, phi = 30,
                     field = "uS", range = c(50, 60), z = TRUE, trace = FALSE)
  expect_true(all(c("sigma_u", "sigma_s", "sigma_z") %in% names(fus$hyper)))
  expect_true(all(is.finite(predict(fus, v, tt$nodes, g, pairs = pairs)$prevalence)))
})

test_that("zero-inflated likelihood: derivatives agree with finite differences", {
  set.seed(3)
  S <- list(n = c(50, 40, 30, 60), y = c(0, 3, 0, 12))
  a <- c(-1, 0.5, 1.2, 0.3); b <- c(-2, -1.5, -4, -1)
  f <- function(a, b) odispers:::odisp_zi_lik(a, b, S, 0.6, 0.8)$ll
  d <- odispers:::odisp_zi_lik(a, b, S, 0.6, 0.8)
  h <- 1e-5
  for (i in seq_along(a)) {
    ea <- replace(numeric(4), i, h)
    num_a <- (f(a + ea, b) - f(a - ea, b)) / (2 * h); num_b <- (f(a, b + ea) - f(a, b - ea)) / (2 * h)
    expect_equal(d$ga[i], num_a, tolerance = 1e-6); expect_equal(d$gb[i], num_b, tolerance = 1e-6)
    ga <- function(a, b) odispers:::odisp_zi_lik(a, b, S, 0.6, 0.8)$ga[i]
    gb <- function(a, b) odispers:::odisp_zi_lik(a, b, S, 0.6, 0.8)$gb[i]
    expect_equal(d$haa[i], (ga(a + ea, b) - ga(a - ea, b)) / (2 * h), tolerance = 1e-5)
    expect_equal(d$hab[i], (ga(a, b + ea) - ga(a, b - ea)) / (2 * h), tolerance = 1e-5)
    expect_equal(d$hbb[i], (gb(a, b + ea) - gb(a, b - ea)) / (2 * h), tolerance = 1e-5)
  }
})

test_that("zero-inflated model: Laplace mode is stationary and the fit runs", {
  tt <- toy_u_setup()
  v <- tt$v[sapply(seq_len(nrow(tt$v)), function(i) any(tt$pairs$i == i & tt$pairs$r < 30)), ]
  pairs <- river_pairs(v, tt$nodes, 40)
  g <- river_graph(tt$nodes, vertex_spacing = 5)
  S <- odispers:::odisp_u_setup(y ~ log(AI), ~ log(DIS_AV_CMS), v, tt$nodes, g, pairs, 30, "n", "AI", 10)
  R <- wm_precision(g, 50)
  fld <- list(field = "u", R = R, Rch = Matrix::Cholesky(R, LDL = FALSE, perm = TRUE),
              sig_u = 0.6, sig_z = 0.4, alpha = 2, gamma = 0.5)
  x0 <- list(beta = S$m_beta, theta = S$m_theta, u = rep(0, g$nv), s = NULL, z = rep(0, length(S$y)))
  res <- odispers:::odisp_zi_laplace(fld, S, x0)
  expect_true(is.finite(res$logLik))
  f <- function(v) {
    x <- odispers:::odisp_u_unvec(v, x0); ev <- odispers:::odisp_zi_eval(x, S)
    odispers:::odisp_zi_lik(ev$a, ev$b, S, fld$gamma, fld$alpha)$ll - 0.5 * odispers:::odisp_u_prior(x, S, fld)
  }
  xv <- odispers:::odisp_u_vec(res$x)
  idx <- c(1:3, 3 + c(1, 5), length(xv) - 0:1)
  h <- 1e-5
  num <- sapply(idx, function(k) { e <- replace(numeric(length(xv)), k, h); (f(xv + e) - f(xv - e)) / (2 * h) })
  expect_lt(max(abs(num)), 5e-3)
  fz <- odisp_fit_zi(y ~ log(AI), ~ log(DIS_AV_CMS), v, tt$nodes, g, pairs, phi = 30, range_u = 50,
                     fixed = list(sigma_u = 0.6, sigma_z = 0.4), trace = FALSE)
  expect_s3_class(fz, "odisp_zifit")
  expect_true(all(fz$pi >= 0 & fz$pi <= 1))
})
