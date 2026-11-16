test_that("Gaussian kernels stay finite at representable extreme widths", {
  # A broad truncated Gaussian tends to a uniform density on its support.
  for (width in c(1e6, 1e100, 1e200)) {
    expect_equal(lgcp_gaussian_kernel(c(0, 1, 4, 5), width, 5, "line"),
                 c(.1, .1, .1, 0), tolerance = 1e-10)
    expect_equal(lgcp_gaussian_kernel(c(0, 1, 4, 5), width, 5, "area"),
                 c(rep(1/(25*pi), 3), 0), tolerance = 1e-10)
  }
  # Rescale by width to check narrow peaks without integrating a tiny interval.
  width <- 1e-200
  expect_equal(lgcp_gaussian_kernel(width*c(0, 1, 2), width, 5, "line") * width,
               dnorm(c(0, 1, 2)), tolerance = 1e-12)
  expect_identical(lgcp_gaussian_kernel(c(1, 5), width, 5, "area"), c(0, 0))
  # Just below the broad-kernel series branch, ordinary quadrature is accurate.
  for (scaling in c("line", "area")) {
    integrand <- function(r) {
      weight <- lgcp_gaussian_kernel(r, 501, 5, scaling)
      if (scaling == "line") 2*weight else 2*pi*r*weight
    }
    expect_equal(integrate(integrand, 0, 5)$value, 1, tolerance = 1e-12)
  }
})

test_that("feasible initial counts reuse shared sites and ignore explicit zero weights", {
  K <- Matrix::sparseMatrix(i = c(1, 2, 2, 3, 3), j = c(2, 2, 3, 1, 4),
                            x = c(.5, .1, 1, 0, .8), dims = c(4, 5))
  expect_identical(lgcp_initial_counts(K, c(1, 1, 1, 0)), c(0L, 1L, 0L, 1L, 0L))
  expect_identical(lgcp_initial_counts(K, numeric(4)), integer(5))
  expect_error(lgcp_initial_counts(K, c(0, 0, 0, 1)), "no positive kernel weight")
  K@x[] <- 0
  expect_error(lgcp_initial_counts(K, c(1, 0, 0, 0)), "kernel width")
})

test_that("unrepresentable kernels fail with a useful message before sampling", {
  a <- lgcp_test_data()
  fit <- function(width, scaling = "line") {
    odisp_fit_lgcp(y ~ 0, ~ 0, a$data, a$nodes, a$graph, phi = 5,
                   kernel = "gaussian", width = width, scaling = scaling,
                   n_iter = 4, burnin = 2)
  }
  # The second positive village is half a kilometre from the only midpoint.
  expect_error(fit(1e-200), "no positive kernel weight")
  expect_error(fit(1e-200, "area"), "kernel weights must be finite")
})

test_that("masked occurrence predictors do not enter count or exposure likelihoods", {
  K <- Matrix::Matrix(rbind(c(.8, .2), c(.3, .7)), sparse = TRUE)
  S <- list(K = K, kernel_tree = lgcp_tree_layout(K), y = c(1, 0), n = c(3, 0),
            link = "logit", gate_base = c(.2, -1000), gate_slope = 2,
            activity_sign = c(1L, 0L))
  counts <- c(1L, 2L)
  exposure <- as.numeric(K %*% counts)
  low_likelihood <- lgcp_exposure_loglik(exposure, c(0, 0), S)
  low_sweep <- withr::with_seed(112,
    lgcp_count_sweep(counts, c(1, 1), c(0, 0), S, 1L))
  S$gate_base[2] <- 1000
  expect_identical(lgcp_exposure_loglik(exposure, c(0, 0), S), low_likelihood)
  high_sweep <- withr::with_seed(112,
    lgcp_count_sweep(counts, c(1, 1), c(0, 0), S, 1L))
  expect_identical(high_sweep, low_sweep)
})
