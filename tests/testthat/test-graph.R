test_that("Whittle-Matern precision has exponential covariance on a chain", {
  g <- list(from = 1:400, to = 2:401, length = rep(0.5, 400), nv = 401)
  S <- as.matrix(solve(wm_precision(g, range = 20, sigma = 1.5)))
  expect_equal(S[201, 201], 1.5^2, tolerance = 1e-8)
  expect_equal(S[201, 221] / S[201, 201], exp(-0.1 * 10), tolerance = 1e-6)
  expect_equal(S[1, 1], 2 * 1.5^2, tolerance = 1e-6)     # doubled at a river end
})

test_that("subdividing an edge leaves the field unchanged", {
  g2 <- list(from = c(1, 2, 2), to = c(2, 3, 4), length = c(3, 7, 5), nv = 4)
  g3 <- list(from = c(1, 5, 2, 2), to = c(5, 2, 3, 4), length = c(1, 2, 7, 5), nv = 5)
  S2 <- as.matrix(solve(wm_precision(g2, 10)))
  S3 <- as.matrix(solve(wm_precision(g3, 10)))
  expect_equal(S2, S3[1:4, 1:4], tolerance = 1e-10)
})

test_that("river_graph joins reaches at confluences and interpolates nodes", {
  # two tributaries (ids 1, 2) flowing into reach 3, which is an outlet
  nodes <- data.frame(line = rep(1:3, c(3, 2, 4)), x = 0, y = 0,
                      len = rep(c(2, 1.5, 3), c(3, 2, 4)),
                      HYRIV_ID = rep(1:3, c(3, 2, 4)), NEXT_DOWN = rep(c(3, 3, 0), c(3, 2, 4)))
  class(nodes) <- c("odisp_network", "data.frame")
  g <- river_graph(nodes, vertex_spacing = 5)
  # reach 1 (6 km) has 2 edges, reach 2 (3 km) 1, reach 3 (12 km) 3
  expect_equal(length(g$from), 6)
  expect_equal(sum(g$length), 6 + 3 + 12)
  expect_equal(g$nv, 3 + 1 + 2 + 1)              # 3 upstream ends, interior, 1 outlet
  expect_equal(as.vector(Matrix::rowSums(g$A)), rep(1, 9))
  # the downstream end of reaches 1 and 2 is the upstream vertex of reach 3
  end1 <- g$to[2]; end2 <- g$to[3]
  expect_equal(end1, 3L); expect_equal(end2, 3L)
})
