# Extracted from test-fit-u.R:252

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "odispers", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
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

# test -------------------------------------------------------------------------
tt <- toy_u_setup(m=40,seed=94)
g <- directed_graph(river_graph(tt$nodes,vertex_spacing=5),"taildown")
S <- odisp_u_setup(y~AI,~log(DIS_AV_CMS),tt$v,tt$nodes,g,tt$pairs,30,"n","AI",10)
S$m_c <- c(0,0)
R <- wm_precision(g,50)
fld <- list(R=R,Rch=Matrix::Cholesky(R,LDL=FALSE),sig_u=.6,sig_z=.4,aE=.7,gamma=-.25)
x <- list(beta=c(-3,.2),c=c(-.4,.1),theta=.1,u=rep(.08,g$nv),z=rep(.05,length(S$y)))
xv <- odisp_zil_vec(x)
delta <- sin(seq_along(xv))/sqrt(length(xv))
evaluate <- function(vv) {
    state <- odisp_zil_unvec(vv,x)
    ev <- odisp_zil_eval(state,S,fld); lk <- odisp_zil_lik(ev$d,ev$e,S)
    op <- odisp_zil_ops(state,ev,lk,S,fld)
    centre <- odisp_zil_vec(list(beta=S$m_beta,c=S$m_c,theta=S$m_theta,u=0*x$u,z=0*x$z))
    list(value=lk$ll-.5*odisp_zil_prior(state,S,fld),
         gradient=op$Jtw(lk$gd,lk$ge)-op$Pv(vv-centre),hv=op$Hv(delta))
  }
for (q in c("linear","log")) for (p in c("linear","log")) {
    fld$occurrence_scale <- q; fld$risk_scale <- p
    h <- 1e-5; center <- evaluate(xv); plus <- evaluate(xv+h*delta); minus <- evaluate(xv-h*delta)
    expect_equal(sum(center$gradient*delta),(plus$value-minus$value)/(2*h),tolerance=2e-6)
    expect_equal(center$hv,(minus$gradient-plus$gradient)/(2*h),tolerance=2e-6)
  }
