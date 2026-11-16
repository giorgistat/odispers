test_that("truncated Gaussian kernels have the requested normalization and support", {
  for (width in c(.2, 2, 20)) {
    line <- integrate(function(r) lgcp_gaussian_kernel(abs(r), width, 5, "line"), -5, 5)$value
    area <- integrate(function(r) 2*pi*r*lgcp_gaussian_kernel(r, width, 5, "area"), 0, 5)$value
    expect_equal(line, 1, tolerance=1e-8)
    expect_equal(area, 1, tolerance=1e-8)
    expect_identical(lgcp_gaussian_kernel(c(5,6),width,5,"line"), c(0,0))
  }
  a <- lgcp_test_data()
  S <- lgcp_setup(y~0,~0,a$data,a$nodes,a$graph,NULL,5,2,.55,.8,"n","line",2,
                  kernel="gaussian",width=1.4)
  expected <- c(dnorm(c(0,.5),sd=1.4)/(2*pnorm(5/1.4)-1),0)
  expect_equal(as.numeric(S$K),expected)
  expect_equal(S$kernel_distance,c(0,.5))
  # Per-site weights do not acquire a second river-cell length factor.
  expect_equal(as.numeric(S$K %*% 2),2*expected)
  expect_equal(exp(S$log_mass),1.6)
})

test_that("width proposals use marked mass and both parts of the village likelihood", {
  distances <- c(.2,1.5,2.5,.4)
  S <- list(K=Matrix::Matrix(matrix(lgcp_gaussian_kernel(distances,1,3,"line"),2),sparse=TRUE),
            kernel_distance=distances,kernel_cutoff=3,kernel_scaling="line",
            y=c(1,0),n=c(4,0),link="logit",gate_base=c(-.2,.8),
            gate_slope=15,activity_sign=c(1L,-1L))
  counts <- c(2L,1L); high <- c(1L,1L); spread <- .7
  predictor <- c(-.3,.1)
  mass <- c(2,1.7)
  exposure <- as.numeric(S$K %*% mass)
  reference <- function(width, with_gate=TRUE) {
    K <- matrix(dnorm(distances,sd=width)/(pnorm(3/width)-pnorm(-3/width)),2)
    E <- as.numeric(K %*% mass)
    risk <- plogis(log(E)+predictor)
    # The sampler omits binomial coefficients, which cancel in MH ratios.
    value <- sum(dbinom(S$y,S$n,risk,log=TRUE))-sum(lchoose(S$n,S$y))
    if (with_gate) {
      eta <- S$gate_base+15*E
      # Evaluate the inactive tail directly; subtracting a rounded q from one
      # would make the independent reference less accurate than the sampler.
      probability_of_state <- c(plogis(eta[1]),plogis(eta[2],lower.tail=FALSE))
      value <- value + sum(log(probability_of_state))
    }
    value
  }
  expect_equal(lgcp_exposure_loglik(exposure,predictor,S),reference(1),tolerance=1e-12)
  different_decision <- FALSE
  original <- S$K
  for (seed in seq_len(30)) {
    candidate <- withr::with_seed(seed, exp(rnorm(1,sd=.5)))
    stopifnot(candidate>.01,candidate<100)
    log_u <- withr::with_seed(seed, {rnorm(1,sd=.5); log(runif(1))})
    ratio <- reference(candidate)-reference(1)
    without_gate <- reference(candidate,FALSE)-reference(1,FALSE)
    different_decision <- different_decision || ((log_u<ratio)!=(log_u<without_gate))
    update <- withr::with_seed(seed,
      lgcp_width_step(1,counts,high,spread,exposure,predictor,S,c(.01,100),.5))
    expect_identical(update$accepted,log_u<ratio)
    expect_equal(update$width,if(log_u<ratio) candidate else 1)
    expect_equal(update$exposure,as.numeric(update$K %*% mass),tolerance=1e-12)
    expect_identical(update$K@i,original@i)
    expect_identical(update$K@p,original@p)
    expect_identical(S$K,original)
  }
  expect_true(different_decision)
})

test_that("sampled widths reconstruct marked exposures with a fixed cutoff", {
  a <- lgcp_test_data()
  for (field in c("wm","taildown")) {
    graph <- if(field=="wm") a$graph else directed_graph(a$graph,"taildown")
    fit <- withr::with_seed(71107,
      odisp_fit_lgcp(y~0,~0,a$data,a$nodes,graph,phi=5,range=2,sigma_u=.55,
        source_rate=.8,n_iter=180,burnin=40,keep_counts=TRUE,keep_prevalence=TRUE,
        marks="two_point",link="logit",occurrence=~1,occurrence_exposure=TRUE,
        relocate=TRUE,field_update="block_slice",kernel="gaussian",width=1.2,
        estimate_width=TRUE,width_bounds=c(.25,10)))
    widths <- fit$draws[,"kernel_width"]
    spread <- fit$draws[,"mark_spread"]
    high <- as.numeric(fit$high_draws)
    mass <- (as.numeric(fit$count_draws)-high)*(1-spread)+high*(1+spread)
    E <- t(vapply(seq_along(widths),function(k) {
      mass[k]*c(dnorm(c(0,.5),sd=widths[k])/(2*pnorm(5/widths[k])-1),0)
    }, numeric(3)))
    expect_gt(length(unique(widths)),5)
    expect_true(all(widths>=.25 & widths<=10))
    expect_equal(fit$fitted$exposure,colMeans(E),tolerance=1e-12)
    expect_equal(fit$prevalence_draws[,1:2],E[,1:2]/(1+E[,1:2]),tolerance=1e-12)
    expect_true(all(fit$prevalence_draws[,3]==0))
    expect_true(is.finite(fit$acceptance[["width"]]))
    expect_true(fit$settings$estimate_width)
    expect_equal(fit$settings$phi,5)
    expect_null(fit$settings$width)
  }
})

test_that("Gaussian width controls reject inconsistent inputs", {
  a <- lgcp_test_data()
  fit <- function(...) odisp_fit_lgcp(y~0,~0,a$data,a$nodes,a$graph,n_iter=10,burnin=5,...)
  expect_error(fit(estimate_width=TRUE),"requires the Gaussian")
  expect_error(fit(kernel="gaussian",width=0),"finite positive")
  expect_error(fit(kernel="gaussian",estimate_width=TRUE,width_bounds=c(10,1)),"positive increasing")
  expect_error(fit(kernel="gaussian",estimate_width=TRUE,width=3,width_bounds=c(4,10)),"containing")
  fixed <- withr::with_seed(93,fit(kernel="gaussian",width=2))
  expect_true(all(fixed$draws[,"kernel_width"]==2))
  expect_equal(fixed$settings$width,2)
})
