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
  sig <- c(0.7, 0.4)
  x0 <- list(beta = S$m_beta, theta = S$m_theta, u = rep(0, g$nv), z = rep(0, length(S$y)))
  res <- odispers:::odisp_u_laplace(sig, S, R, Rch, x0)
  expect_true(is.finite(res$logLik))
  f <- function(v) {
    x <- odispers:::odisp_u_unvec(v, x0)
    odispers:::odisp_u_lik(odispers:::odisp_u_eta(x, S)$eta, S)$ll - 0.5 * odispers:::odisp_u_prior(x, S, R, sig)
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
