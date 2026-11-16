# Reconstruct every model/fold's first replicate, and audit the held-out rows and survey totals.
library(odispers)
library(Matrix)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) %in% 1:2)
out <- args[1]; a <- readRDS(file.path(out,"inputs.rds"))
v <- a$villages
folds <- read.csv(file.path(out,"folds.csv"))
stopifnot(identical(a$spatial_fold,folds$fold),all(table(folds$fold)==180),
          all(tapply(folds$fold,folds$block,function(x) length(unique(x)))==1L),
          sum(v$y==0)==550L,sum(v$y)==2386, sum(v$n)==37831)
factor <- Cholesky(wm_precision(a$graph,25,.5),LDL=FALSE)
models <- c("baseline","river_climate","occurrence","both")
if (length(args)==2L && args[2] %in% c("flexible","design")) models <- c(models,"flexible_occurrence")
if (length(args)==2L && args[2]=="design") models <- c(models,"survey_design")
if (length(args)==2L && args[2]=="unit") models <- "unit_occurrence"
if (length(args)==2L && args[2]=="exposure") models <- "exposure_occurrence"
for (m in models) for (f in if(m=="unit_occurrence") 0L else 0:5) {
  path <- file.path(out,m,paste0("fold",f))
  expected_test <- if (f) which(a$spatial_fold==f) else seq_len(900)
  expected_training <- if(f) which(a$spatial_fold!=f) else seq_len(900)
  y <- readRDS(file.path(path,"replications.rds"))
  data <- read.csv(file.path(path,"predictions.csv"))
  first <- readRDS(file.path(path,"first_landscape.rds"))
  files <- if(m=="baseline" && f==0L) file.path("jonas_local/lgcp_full_20261106/logit/observed",
    paste0("two_point_",c(8106,8206,8306,8406),".rds")) else list.files(path,"^chain_[0-9]+[.]rds$",full.names=TRUE)
  expected_reps <- if(m=="exposure_occurrence") { if(f) 500L else 1000L } else { if(f) 400L else 600L }
  stopifnot(if(f) length(files) %in% c(2L,4L) else length(files)==4L,identical(y$training,expected_training),
            identical(y$test,expected_test),identical(data$village,expected_test),identical(first$test,expected_test),
            length(y$seeds)==expected_reps,all(y$y>=0),all(sweep(y$y,2,v$n[expected_test],"-")<=0),
            max(abs(data$prob_zero-data$prob_inactive-data$prob_zero_exposure_active-data$prob_binomial_zero_active_exposed))<1e-12,
            max(abs(data$simulated_prob_zero-colMeans(y$y==0)))<1e-12)
  pooled <- do.call(rbind,lapply(files,function(file) {
    fit <- readRDS(file)
    stopifnot(isTRUE(fit$settings$occurrence_exposure)==(m=="exposure_occurrence"))
    if(!is.null(fit$training)) stopifnot(identical(fit$training,expected_training),fit$fold==f)
    stopifnot(fit$settings$link=="logit",fit$settings$source_rate==.04,fit$marks==if(m=="unit_occurrence") "unit" else "two_point",
              fit$settings$range==25,fit$settings$phi==40,fit$settings$sigma_u==.5,
              (is.null(fit$settings$occurrence))==(m %in% c("baseline","river_climate")))
    fit$draws
  }))
  set.seed(y$seeds[1]); parameters <- pooled[sample.int(nrow(pooled),1L),]
  z <- rnorm(a$graph$nv)
  u <- as.numeric(solve(factor,solve(factor,z,system="Lt"),system="Pt"))
  terms <- a$source_covariates
  if (m %in% c("river_climate","both")) terms <- c(terms,"river_dry")
  design <- as.matrix(as.data.frame(a$nodes)[,terms])
  rate <- exp(log(.04)+as.numeric(a$graph$A%*%u)+as.numeric(design%*%parameters[paste0("theta.",terms)]))
  stopifnot(identical(parameters,first$parameters),max(abs(u-first$u))<1e-12,
            max(abs(rate-first$rate)/pmax(1,rate))<1e-12)
  sites <- rpois(nrow(a$nodes),a$nodes$len*rate)
  high <- if(m=="unit_occurrence") integer(length(sites)) else rbinom(length(sites),sites,.5)
  inside <- a$pairs$r<40 & a$pairs$i %in% expected_test
  d <- if(m=="unit_occurrence") 0 else unname(parameters["mark_spread"])
  distance <- a$pairs$r[inside]/40
  k <- 3/80*(1-distance)^4*(1+4*distance)
  mass <- (sites-high)*(1-d)+high*(1+d)
  groups <- rowsum(k*mass[a$pairs$j[inside]],a$pairs$i[inside])
  E <- numeric(length(expected_test))
  E[match(as.integer(rownames(groups)),expected_test)] <- groups[,1]
  eta <- parameters["beta.(Intercept)"]+parameters["beta.z_logAI"]*v$z_logAI[expected_test]
  odds <- E*exp(eta)
  p <- odds/(1+odds)
  g <- if(m %in% c("occurrence","both","flexible_occurrence","survey_design","unit_occurrence","exposure_occurrence")) {
    climate <- (v$AI[expected_test]-.4)/.2
    log_odds <- parameters["occurrence.(Intercept)"]+parameters["occurrence.climate"]*climate
    if (m=="exposure_occurrence") log_odds <- log_odds+parameters["occurrence.exposure"]*E
    if (m %in% c("flexible_occurrence","survey_design")) log_odds <- log_odds+parameters["occurrence.I(pmin(climate + 0.75, 0))"]*pmin(climate+.75,0)
    if (m=="survey_design") log_odds <- log_odds+parameters["occurrence.examined30"]*(v$n[expected_test]==30)+parameters["occurrence.examined50"]*(v$n[expected_test]==50)
    plogis(log_odds)
  } else rep(1,length(expected_test))
  active <- as.integer(runif(length(expected_test))<g)
  simulated <- rbinom(length(expected_test),v$n[expected_test],active*p)
  stopifnot(identical(sites,first$N),identical(high,first$H),max(abs(E-first$E))<1e-10,
    max(abs(p-first$p))<1e-12,max(abs(g-first$g))<1e-12,identical(active,first$active),
    identical(simulated,as.integer(first$y)),identical(simulated,as.integer(y$y[1,])))
  cat(m,f,length(expected_training),"training rows,",length(expected_test),"predicted rows,",
      nrow(y$y),"fresh landscapes and seed reconstruction verified\n")
}
