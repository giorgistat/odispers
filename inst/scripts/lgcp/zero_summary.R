# Aggregate spatial holdout predictions; retain model/fold and aridity-specific calibration.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) %in% 1:2)
root <- args[1]
v <- read.csv(file.path(root,"folds.csv"))
models <- c("baseline","river_climate","occurrence","both")
if (length(args)==2L && args[2] %in% c("flexible","design")) models <- c(models,"flexible_occurrence")
if (length(args)==2L && args[2]=="design") models <- c(models,"survey_design")
all_rows <- list(); groups <- list(); diagnostic <- list()
raw <- list()
for (model in models) {
  for (source in c("full","heldout")) {
    folds <- if(source=="full") 0L else seq_len(5)
    folders <- file.path(root,model,paste0("fold",folds))
    predictions <- do.call(rbind,lapply(folders,function(d) read.csv(file.path(d,"predictions.csv"))))
    reps <- lapply(folders,function(d) readRDS(file.path(d,"replications.rds")))
    agreement <- all(vapply(reps,function(r) r$fold==0L || length(intersect(r$test,r$training))==0L,TRUE))
    stopifnot(agreement,length(unique(predictions$village))==900L,
              identical(unname(unlist(lapply(reps,"[[","test"))),predictions$village))
    order_id <- order(predictions$village)
    pred <- predictions[order_id,]
    Y <- do.call(cbind,lapply(reps,"[[","y"))[,order_id,drop=FALSE]
    stopifnot(identical(pred$y,v$y),identical(pred$n,v$n),identical(pred$village,v$village))
    aggregate_group <- function(ids,key,group) {
      q <- function(x) unname(quantile(x,c(.025,.5,.975)))
      zeros <- q(rowSums(Y[,ids,drop=FALSE]==0))
      prevalence <- q(rowSums(Y[,ids,drop=FALSE])/sum(v$n[ids]))
      high <- q(rowSums(sweep(Y[,ids,drop=FALSE],2,ceiling(v$n[ids]/2),">=")))
      observed <- v$y[ids]==0
      predicted_positive <- sum(1-pred$prob_zero[ids])
      data.frame(model=model,assessment=source,partition=key,group=group,n=length(ids),
                 observed_zeros=sum(observed),expected_zeros=sum(pred$prob_zero[ids]),
                 lo=zeros[1],median=zeros[2],hi=zeros[3],
                 inactive=sum(pred$prob_inactive[ids]),active_no_sites=sum(pred$prob_zero_exposure_active[ids]),
                 active_binomial=sum(pred$prob_binomial_zero_active_exposed[ids]),
                 brier=mean((pred$prob_zero[ids]-observed)^2),crps=mean(pred$crps[ids]),
                 mean_negative_log_mass=mean(pred$neglog_score[ids]),
                 observed_fraction=sum(v$y[ids])/sum(v$n[ids]),
                 observed_cases_per_positive=if(any(!observed)) sum(v$y[ids])/sum(!observed) else NA_real_,
                 expected_cases_per_positive=if(predicted_positive>0) sum(pred$n[ids]*pred$mean_prevalence[ids])/predicted_positive else NA_real_,
                 frac_lo=prevalence[1],frac_median=prevalence[2],frac_hi=prevalence[3],
                 observed_high=sum(v$y[ids]>=ceiling(v$n[ids]/2)),
                 high_lo=high[1],high_median=high[2],high_hi=high[3],
                 conditional_zeros=sum(pred$conditional_zero[ids]),
                 conditional_brier=mean((pred$conditional_zero[ids]-observed)^2),
                 simulation_mcse=sd(rowSums(Y[,ids,drop=FALSE]==0))/sqrt(nrow(Y)))
    }
    all_rows[[length(all_rows)+1L]] <- aggregate_group(1:900,"all","all")
    partitions <- list(aridity=cut(v$AI,quantile(v$AI,0:4/4),include.lowest=TRUE),
                       examined=ifelse(v$n==30,"30",ifelse(v$n==50,"50","other")),
                       spatial_fold=v$fold)
    for(key in names(partitions)) for(label in unique(as.character(partitions[[key]]))) {
      ids <- which(as.character(partitions[[key]])==label)
      groups[[length(groups)+1L]] <- aggregate_group(ids,key,label)
    }
    for (i in seq_along(folds)) {
      d <- read.csv(file.path(folders[i],"diagnostics.csv"))
      diagnostic[[length(diagnostic)+1L]] <- data.frame(model=model,fold=folds[i],rhat=max(d$rhat),minimum_ess=min(d$ess))
    }
    if(source=="heldout") raw[[model]] <- pred
  }
}
tab <- do.call(rbind,all_rows)
write.csv(tab,file.path(root,"comparison.csv"),row.names=FALSE)
write.csv(do.call(rbind,groups),file.path(root,"groups.csv"),row.names=FALSE)
write.csv(do.call(rbind,diagnostic),file.path(root,"diagnostics.csv"),row.names=FALSE)
# Paired resampling of 100 km blocks; descriptive uncertainty conditional on the fitted folds.
set.seed(12037)
block <- as.integer(factor(v$block))
bs <- replicate(1500,sample.int(max(block),max(block),replace=TRUE))
diffs <- list()
for (model in setdiff(models,"baseline")) for(comparator in setdiff(intersect(models,c("baseline","river_climate","occurrence","flexible_occurrence")),model)) {
  x <- raw[[model]]; y <- raw[[comparator]]
  stats <- list(brier=(x$prob_zero-(v$y==0))^2-(y$prob_zero-(v$y==0))^2,crps=x$crps-y$crps)
  for (metric in names(stats)) {
    value <- stats[[metric]]
    by_block <- tapply(value,block,sum)
    size <- tabulate(block)
    delta <- apply(bs,2,function(ids) sum(by_block[ids])/sum(size[ids]))
    interval <- unname(quantile(delta,c(.025,.975)))
    diffs[[length(diffs)+1L]] <- data.frame(model=model,comparator=comparator,metric=metric,delta=mean(value),lo=interval[1],hi=interval[2])
  }
}
write.csv(do.call(rbind,diffs),file.path(root,"paired_scores.csv"),row.names=FALSE)
print(tab[,c("model","assessment","expected_zeros","median","lo","hi","brier","crps","frac_median","high_median")],row.names=FALSE)

# Show the count deficit as well as geographic variation; common scale across the four panels.
palette <- c("#176D73","#BF6A3C","#293F69","#714E77","#527047","#996C18")[seq_along(models)]
png(file.path(root,"zero_calibration.png"),width=1740,height=850,res=165)
layout(matrix(c(1,1,2,3,4,4),3,byrow=TRUE),heights=c(.65,4.2,.55))
par(mar=rep(0,4),family="sans",fg="#233D41")
plot.new();text(.5,.72,"Predicting zero-positive village surveys",font=2,cex=1.55)
text(.5,.18,"Full river network  ·  900 villages  ·  Five spatial holdouts",cex=.97,col="#4B6568")
lab <- c("Original LGCP","River aridity","Village occurrence","Both","Flexible occurrence","Survey design")[seq_along(models)]
par(mar=c(3.9,8.8,2.2,1),las=1)
nmodels <- length(models)
plot.new();plot.window(c(190,700),c(.65,nmodels+.65))
abline(v=seq(200,700,100),col="#E0E7E6")
abline(v=550,lty=2,col="#545454");text(550,nmodels+.51,"550 observed",cex=.82)
a <- tab[tab$assessment=="full",]
for(i in seq_along(models)) {
  segments(a$lo[i],nmodels+1-i,a$hi[i],nmodels+1-i,col=palette[i],lwd=3)
  points(a$median[i],nmodels+1-i,col=palette[i],pch=19,cex=1.2)
}
axis(1,at=seq(200,700,100));axis(2,at=nmodels:1,labels=lab,tick=FALSE,cex.axis=.87)
mtext("Zero surveys in new full-river landscapes",side=1,line=2.25,cex=.84)
mtext("Full-data predictive intervals",3,line=.7,font=2)
par(mar=c(3.9,4.2,2.2,.8))
observed <- as.numeric(tapply(v$y==0,partitions$aridity,sum))
plot.new();plot.window(c(.75,4.25),c(-4,230));abline(h=seq(0,225,50),col="#E0E7E6")
for(i in seq_along(models)) {
  by_q <- tapply(raw[[models[i]]]$prob_zero,partitions$aridity,sum)
  lines(1:4,by_q,col=palette[i],lwd=2,pch=i+1,type="b")
}
lines(1:4,observed,col="#252E31",lwd=2.4,pch=16,type="b",lty=2)
axis(1,at=1:4,labels=c("Driest","2","3","Wettest"));axis(2,at=seq(0,225,50))
mtext("Village aridity quarter",1,line=2.25,cex=.86)
mtext("Zero surveys",2,line=2.6,las=0,cex=.84)
mtext("Held-out villages by aridity",3,line=.7,font=2)
par(mar=rep(0,4));plot.new()
legend("center",c(lab,"Observed"),col=c(palette,"#252E31"),lty=c(rep(1,nmodels),2),pch=c(seq_along(models)+1,16),horiz=TRUE,bty="n",cex=.78)
dev.off()
