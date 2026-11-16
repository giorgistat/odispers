# Sum over the marked counts and integrate the Gaussian effects in a one-cell example.
library(odispers)
library(sf)
crs <- "+proj=aeqd +lat_0=0 +lon_0=30 +datum=WGS84 +units=km +no_defs"
river <- st_sf(HYRIV_ID=1L,NEXT_DOWN=0L,geometry=st_sfc(st_linestring(rbind(c(0,0),c(2,0))),crs=crs))
nodes <- river_nodes(river,spacing=3,crs=crs)
graph <- river_graph(nodes)
v <- st_as_sf(data.frame(x=c(1,1.5),coord_y=0,n=c(14,10),y=c(6,0)),coords=c("x","coord_y"),crs=crs)
K <- c(1,.9^4*1.4)/(5*2/3)
var_u <- as.numeric(graph$A %*% solve(as.matrix(wm_precision(graph,2,.55))) %*% Matrix::t(graph$A))
n_values <- 0:80
mu <- function(z) 1.6*exp(sqrt(var_u)*z)
pN <- vapply(n_values,function(n) integrate(function(z) dnorm(z)*dpois(n,mu(z)),-10,10,abs.tol=1e-15)$value,0)
tail_bound <- integrate(function(z) dnorm(z)*ppois(80,mu(z),lower.tail=FALSE),-10,10,abs.tol=1e-15)$value + 2*pnorm(-10)
# First village must be active. The second contributes either g*(1-g)
# or g^2 times its binomial zero likelihood. Integrate these two terms in gamma.
moment <- function(active,power) integrate(function(gamma) {
  q <- plogis(gamma)
  dnorm(gamma,sd=1.5)*gamma^power*q*if(active) q else (1-q)
},-15,15,abs.tol=1e-15)$value
inactive <- vapply(0:1,function(k) moment(FALSE,k),0)
active <- vapply(0:1,function(k) moment(TRUE,k),0)
reference <- function(steps) {
  d <- (seq_len(steps)-.5)/steps
  total <- numeric(5)
  for (n in n_values) for(h in 0:n) {
    mass <- (n-h)*(1-d)+h*(1+d)
    odds1 <- K[1]*mass
    odds2 <- K[2]*mass
    w <- pN[n+1L]*dbinom(h,n,.5)*dbinom(6,14,odds1/(1+odds1))
    z <- (1+odds2)^-10
    wt <- w*(inactive[1]+active[1]*z)
    total <- total + c(sum(wt),n*sum(wt),sum(d*wt),sum(w*(inactive[2]+active[2]*z)),sum(w*active[1]*z))
  }
  setNames(total[-1]/total[1],c("total_sites","mark_spread","occurrence.(Intercept)","zero_village_active"))
}
truth <- reference(512)
stopifnot(max(abs(truth-reference(256)))<1e-5,tail_bound<1e-9)
batches <- lapply(c(1001,1002,1003,1004),function(seed) {
  set.seed(seed)
  fit <- odisp_fit_lgcp(y~0,~0,v,nodes,graph,phi=5,range=2,sigma_u=.55,source_rate=.8,
    prior_sd=1.5,field_update="block_slice",marks="two_point",occurrence=~1,link="logit",
    n_iter=5000L,burnin=1000L,keep_prevalence=TRUE)
  # The first positive village requires at least one site. Both kernels are positive.
  # Thus the second village has exactly zero risk iff its state is inactive.
  samples <- cbind(fit$draws[,names(truth)[1:3]],zero_village_active=as.integer(fit$prevalence_draws[,2]>0))
  rowsum(samples,rep(1:40,each=100))/100
})
x <- do.call(rbind,batches)
check <- data.frame(reference=truth,sampled=colMeans(x),mcse=apply(x,2,sd)/sqrt(nrow(x)))
print(check)
stopifnot(all(abs(check$sampled-check$reference)<4.5*check$mcse))
