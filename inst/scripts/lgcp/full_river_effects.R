# Compare environmental effects and zero-survey counts by observed aridity.
# Run after both summaries and forward simulations in a covariate analysis.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L,3L)) stop("supply the covariate analysis directory and optional fit/marked simulation subdirectories")
fit_dir <- if (length(args)==3L) args[2] else "observed"
marked_dir <- if (length(args)==3L) args[3] else "marked_predictive"
out <- args[1]
effects <- lapply(c("unit","two_point"), function(m)
  read.csv(file.path(out,fit_dir,paste0(m,"_covariate_effects.csv"))))
label_map <- c(theta.z_flow="Log discharge", theta.z_slope="River slope", theta.z_flood="Inundation",
               theta.z_lon="River longitude", theta.z_lat="River latitude", beta.z_logAI="Log village aridity")
labels <- unname(label_map[effects[[1]]$parameter])
rows <- rev(seq_along(labels))
stopifnot(identical(effects[[1]]$parameter,effects[[2]]$parameter), !anyNA(labels))
rng <- range(c(1,unlist(lapply(effects,function(e) c(e$lo,e$hi)))))
ticks <- 2^(floor(log2(rng[1])):ceiling(log2(rng[2])))
colour <- c("#247A83","#AB5938")
png(file.path(out,"covariate_effects.png"), width=1500, height=820+100*(length(rows)-4), res=165)
par(mar=c(4.3,8.5,4.5,2),family="sans", fg="#233D41",col.axis="#364C50",las=1,xaxs="i")
plot(NA, xlim=range(ticks)*c(.93,1.07), ylim=c(.6,max(rows)+.65), log="x",
     xlab="", ylab="", axes=FALSE)
for(y in rows) abline(h=y,col="#E3E9E8",lwd=1)
abline(v=1,col="#879798",lty=2)
axis(1,at=ticks,labels=format(ticks,trim=TRUE,digits=3),tck=-.015)
axis(2,at=rows,labels=labels,tick=FALSE,cex.axis=.96)
for(k in 1:2){
  e<-effects[[k]]
  y<-rows + c(.12,-.12)[k]
  segments(e$lo,y,e$hi,y,col=colour[k],lwd=2.5)
  points(e$median,y,pch=19,col=colour[k],cex=1.2)
}
title("Environmental effects in the full river models",col.main="#20383B",cex.main=1.35,line=2.5)
mtext("Posterior medians and 95% credible intervals for a one-SD increase",side=3,line=1.1,
      col="#4B6264",cex=.86)
survey_scale <- if ("village odds ratio" %in% effects[[1]]$scale) "odds" else "hazard"
mtext(sprintf("Multiplier on site density (river predictors) or infection %s (village aridity)",survey_scale),side=1,
      line=2.6,cex=.82)
legend("topright", inset=c(.01,0), legend=c("Unit strengths","Two strengths"),pch=19,
       col=colour,bty="n",cex=.83,horiz=TRUE)
dev.off()

a <- readRDS(file.path(out,"inputs.rds"))
v <- a$villages
band <- cut(v$AI,quantile(v$AI,seq(0,1,.25)),include.lowest=TRUE)
summaries <- list()
for (mode in c("unit","two_point")){
  sims <- if(mode=="unit") readRDS(file.path(out,"null","null_simulations.rds")) else
    readRDS(file.path(out,marked_dir,"marked_simulations.rds"))
  summaries[[mode]] <- do.call(rbind,lapply(split(seq_len(nrow(v)),band),function(i){
    zeros <- rowSums(sims$y[,i,drop=FALSE]==0)
    frac <- rowSums(sims$y[,i,drop=FALSE])/sum(v$n[i])
    z <- unname(quantile(zeros,c(.025,.5,.975)))
    f <- unname(quantile(frac,c(.025,.5,.975)))
    data.frame(model=mode,n=length(i),observed_zeros=sum(v$y[i]==0),zero_lo=z[1],zero_median=z[2],
               zero_hi=z[3],observed_fraction=sum(v$y[i])/sum(v$n[i]),
               fraction_lo=f[1],fraction_median=f[2],fraction_hi=f[3])
  }))
  summaries[[mode]]$aridity_band <- rownames(summaries[[mode]])
}
write.csv(do.call(rbind,summaries),file.path(out,"aridity_prediction_check.csv"),row.names=FALSE)
