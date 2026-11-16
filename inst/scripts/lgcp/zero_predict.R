# Fresh full-graph replicated landscapes and training-landscape predictions.
library(odispers)
library(Matrix)
library(posterior)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==4L)
out <- args[1]
model <- args[2]
fold <- as.integer(args[3])
replicates <- as.integer(args[4])
a <- readRDS(file.path(out,"inputs.rds"))
v <- a$villages
v$examined30 <- as.integer(v$n==30L)
v$examined50 <- as.integer(v$n==50L)
test <- if (fold==0L) seq_len(nrow(v)) else which(a$spatial_fold==fold)
folder <- file.path(out,model,paste0("fold",fold))
dir.create(folder,showWarnings=FALSE,recursive=TRUE)
if (model=="baseline" && fold==0L) {
  path <- file.path("jonas_local/lgcp_full_20261106/logit/observed",paste0("two_point_",c(8106,8206,8306,8406),".rds"))
} else {
  path <- list.files(folder,"^chain_[0-9]+[.]rds$",full.names=TRUE)
}
stopifnot(if (fold==0L) length(path)==4L else length(path) %in% c(2L,4L))
selected <- lapply(path,function(f) {
  fit <- readRDS(f)
  stopifnot(fit$settings$link=="logit",fit$settings$phi==40,fit$settings$range==25,
            fit$settings$sigma_u==.5,fit$settings$source_rate==.04,fit$marks==if(model=="unit_occurrence") "unit" else "two_point",
            ncol(fit$prevalence_draws)==900)
  if (!is.null(fit$training)) stopifnot(identical(fit$training,if(fold) which(a$spatial_fold!=fold) else seq_len(nrow(v))))
  if (fold>0L) stopifnot(identical(fit$fold,fold))
  stopifnot(isTRUE(fit$settings$occurrence_exposure)==(model=="exposure_occurrence"))
  p <- fit$prevalence_draws[,test,drop=FALSE]
  # Conditional on training surveys, use the real test denominator, including when fitting masked n=0.
  zero <- colMeans(exp(sweep(log1p(-p),2,v$n[test],"*")))
  list(draws=fit$draws,conditional_zero=zero,conditional_mean=colMeans(p),seed=fit$seed)
})
pooled <- do.call(rbind,lapply(selected,"[[","draws"))
diagnostics <- do.call(rbind,lapply(colnames(pooled),function(k) {
  x <- do.call(cbind,lapply(selected,function(s) s$draws[,k]))
  q <- unname(quantile(x,c(.025,.5,.975)))
  data.frame(parameter=k,lo=q[1],median=q[2],hi=q[3],mean=mean(x),rhat=rhat(x),ess=ess_bulk(x),mcse=mcse_mean(x))
}))
geo_prior <- Cholesky(wm_precision(a$graph,25,.5),LDL=FALSE)
s <- a$pairs
inside <- s$r<40 & s$i %in% test
kernel <- sparseMatrix(i=match(s$i[inside],test),j=s$j[inside],x=3/80*wendland(s$r[inside]/40),
                       dims=c(length(test),nrow(a$nodes)))
theta_names <- grep("^theta[.]",colnames(pooled),value=TRUE)
Z <- as.matrix(as.data.frame(a$nodes)[,sub("^theta[.]","",theta_names),drop=FALSE])
gamma_names <- grep("^occurrence[.]",colnames(pooled),value=TRUE)
has_gate <- length(gamma_names)>0
if (has_gate) {
  gate_formula <- switch(model, flexible_occurrence=~climate + I(pmin(climate + .75, 0)),
     survey_design=~climate + I(pmin(climate + .75, 0)) + examined30 + examined50, ~climate)
  gate_matrix <- model.matrix(gate_formula,sf::st_drop_geometry(v)[test,,drop=FALSE])
  if (model=="exposure_occurrence") gate_matrix <- cbind(gate_matrix,exposure=numeric(length(test)))
  stopifnot(identical(gamma_names,paste0("occurrence.",colnames(gate_matrix))))
}
B <- replicates
M <- length(test)
yrep <- matrix(0L,B,M)
absent <- empty <- sampling_zero <- risk <- likelihood <- matrix(0,B,M)
# Seed shared across candidate models for a given fold. Each replicate has its own seed
# so a different number of gate draws in one model cannot shift subsequent replicates.
seeds <- 171106L + fold*10000L + seq_len(B)
draw_seed <- numeric(B)
if (model=="exposure_occurrence") gate_effect <- matrix(0,B,5,dimnames=list(NULL,
    c("mean_exposure","q25_exposure","q75_exposure","fresh_expected_zeros","same_draw_slope_zero")))
pooled_seed <- unlist(lapply(selected,function(x) rep(x$seed,nrow(x$draws))),use.names=FALSE)
for (b in seq_len(B)) {
  set.seed(seeds[b])
  index <- sample.int(nrow(pooled),1L)
  parameters <- pooled[index,]
  draw_seed[b] <- pooled_seed[index]
  z <- rnorm(a$graph$nv)
  u <- as.numeric(solve(geo_prior,solve(geo_prior,z,system="Lt"),system="Pt"))
  rate <- exp(log(.04) + as.numeric(a$graph$A %*% u) + as.numeric(Z %*% parameters[theta_names]))
  N <- rpois(nrow(a$nodes),a$nodes$len*rate)
  H <- if (model=="unit_occurrence") integer(length(N)) else rbinom(length(N),N,.5)
  d <- if (model=="unit_occurrence") 0 else unname(parameters["mark_spread"])
  mass <- (N-H)*(1-d) + H*(1+d)
  E <- as.numeric(kernel %*% mass)
  eta <- parameters["beta.(Intercept)"] + v$z_logAI[test]*parameters["beta.z_logAI"]
  p <- plogis(log(E)+eta)
  if (model=="exposure_occurrence") gate_matrix[,"exposure"] <- E
  g <- if (has_gate) plogis(as.numeric(gate_matrix %*% parameters[gamma_names])) else rep(1,M)
  if (model=="exposure_occurrence") {
    no_exposure_term <- plogis(parameters["occurrence.(Intercept)"]+
                               v$climate[test]*parameters["occurrence.climate"])
    prob_detection <- -expm1(v$n[test]*log1p(-p))
    gate_effect[b,] <- c(mean(E),quantile(E,c(.25,.75)),sum(1-g*prob_detection),
                         sum(1-no_exposure_term*prob_detection))
  }
  active <- as.integer(runif(M)<g)
  yrep[b,] <- rbinom(M,v$n[test],active*p)
  absent[b,] <- 1-g
  empty[b,] <- g*(E==0)
  sampling_zero[b,] <- g*(E>0)*dbinom(0,v$n[test],p)
  risk[b,] <- g*p
  likelihood[b,] <- (1-g)*(v$y[test]==0) + g*dbinom(v$y[test],v$n[test],p)
  if (b==1L) first <- list(parameters=parameters,seed=seeds[b],u=u,rate=rate,N=N,H=H,
                           E=E,p=p,g=g,active=active,y=yrep[b,],test=test)
}
prob_zero <- colMeans(absent + empty + sampling_zero)
counts <- v$n[test]
actual <- v$y[test]
crps <- vapply(seq_along(test),function(i) {
  x <- yrep[,i]/counts[i]
  mean(abs(x-actual[i]/counts[i])) - sum((2*seq_len(B)-B-1)*sort(x))/(B*(B-1))
},numeric(1))
# Count likelihood is an independent-landscape Monte Carlo mixture. Flag any zero Monte Carlo masses.
predictive_mass <- colMeans(likelihood)
predictions <- data.frame(village=test,fold=fold,model=model,y=actual,n=counts,AI=v$AI[test],
   prob_zero=prob_zero,prob_inactive=colMeans(absent),prob_zero_exposure_active=colMeans(empty),
   prob_binomial_zero_active_exposed=colMeans(sampling_zero),mean_prevalence=colMeans(risk),
   simulated_prob_zero=colMeans(yrep==0),crps=crps,neglog_score=-log(predictive_mass),
   conditional_zero=rowMeans(do.call(cbind,lapply(selected,"[[","conditional_zero"))),
   conditional_mean=rowMeans(do.call(cbind,lapply(selected,"[[","conditional_mean"))))
stopifnot(all(is.finite(predictions$prob_zero)),all(predictions$prob_zero>=0 & predictions$prob_zero<=1+1e-12),
          all(rowSums(yrep<0L)==0),all(sweep(yrep,2,counts,"-")<=0))
write.csv(diagnostics,file.path(folder,"diagnostics.csv"),row.names=FALSE)
if (model=="exposure_occurrence") write.csv(gate_effect,file.path(folder,"gate_effect.csv"),row.names=FALSE)
chain_replications <- data.frame(seed=draw_seed,zeros=rowSums(absent+empty+sampling_zero),
                                 expected_cases=as.numeric(risk %*% counts))
write.csv(aggregate(chain_replications[c("zeros","expected_cases")],chain_replications["seed"],mean),
          file.path(folder,"per_chain_predictions.csv"),row.names=FALSE)
write.csv(predictions,file.path(folder,"predictions.csv"),row.names=FALSE)
saveRDS(first,file.path(folder,"first_landscape.rds"))
temp <- tempfile(tmpdir=folder)
saveRDS(list(y=yrep,test=test,seeds=seeds,training=if (fold==0L) seq_len(nrow(v)) else which(a$spatial_fold!=fold),
             model=model,fold=fold),temp)
stopifnot(file.rename(temp,file.path(folder,"replications.rds")))
cat(model,fold,B,"landscapes; expected zeros",sum(prob_zero),"observed",sum(actual==0),
    "max Rhat",max(diagnostics$rhat),"min ESS",min(diagnostics$ess),"\n")
