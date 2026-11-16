# Aggregate the 3x3 full-data comparison. Input is the analysis directory
# containing coordinates/ (40,25) and range_grid/phi{r}_range{h}/.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("supply the full LGCP analysis directory")
root <- args[1]
grid <- file.path(root,"range_grid")
ref <- readRDS(file.path(root,"coordinates","inputs.rds"))
ref_nodes <- sort(unique(ref$pairs$j[ref$pairs$r<40]))
radii <- c(40,60,80)
ranges <- c(10,25,100)

case_directory <- function(radius, range, mode) {
  if (radius == 40 && range == 25) return(file.path(root, "coordinates"))
  directory <- file.path(grid, sprintf("phi%d_range%d", radius, range))
  # A longer fit replaces its original only after its independent output check succeeds.
  refined <- file.path(directory, "refined")
  if (file.exists(file.path(refined, paste0("verified_", mode, ".txt")))) return(refined)
  directory
}
# Fair finite-ensemble CRPS on prevalence, averaging the 900 villages equally.
# The pair term uses distinct replicate pairs and can be computed from the order statistics.
score_replications <- function(yrep, survey) {
  ensemble <- sweep(yrep, 2, survey$n, "/")
  actual <- survey$y / survey$n
  score <- function(x) {
    B <- nrow(x)
    weights <- 2*seq_len(B)-B-1
    mean(vapply(seq_len(ncol(x)), function(j) {
      mean(abs(x[,j]-actual[j])) - sum(weights*sort(x[,j])) / (B*(B-1))
    },0))
  }
  batches <- split(seq_len(nrow(ensemble)), rep_len(1:10, nrow(ensemble)))
  batch_scores <- vapply(batches,function(i) score(ensemble[i,,drop=FALSE]),0)
  c(value=score(ensemble),mcse=sd(batch_scores)/sqrt(length(batches)))
}
rows <- list()
effect_rows <- list()
for (radius in radii) for (range in ranges) {
  for(mode in c("unit","two_point")){
    case <- case_directory(radius, range, mode)
    pred_dir <- file.path(case, if(mode=="unit") "null" else "marked_predictive")
    env <- read.csv(file.path(pred_dir, if(mode=="unit") "null_envelope.csv" else "marked_envelope.csv"))
    sim_file <- if(mode=="unit") "null_simulations.rds" else "marked_simulations.rds"
    crps <- score_replications(readRDS(file.path(pred_dir,sim_file))$y, ref$villages)
    g <- read.csv(file.path(case,"observed",paste0(mode,"_parameters.csv")))
    v <- read.csv(file.path(case,"observed",paste0(mode,"_village_risk.csv")))
    e <- read.csv(file.path(case,"observed",paste0(mode,"_covariate_effects.csv")))
    effect_rows[[length(effect_rows)+1L]] <- cbind(radius_km=radius,field_range_km=range,model=mode,e)
    cells <- readRDS(file.path(case,"observed",paste0(mode,"_cells.rds")))
    get <- function(name,col) g[[col]][match(name,g$parameter)]
    total <- match("total_sites",g$parameter)
    zero <- match("zero_surveys",env$statistic)
    frac <- match("positive_fraction",env$statistic)
    high <- match("at_least_half_positive",env$statistic)
    runs <- read.csv(file.path(case,"observed",paste0(mode,"_runs.csv")))
    stopifnot(nrow(v)==900L,nrow(runs)==4L,all(ref_nodes %in% cells$node))
    rows[[length(rows)+1L]] <- data.frame(radius_km=radius, field_range_km=range, model=mode,
      fit_directory=substring(case,nchar(root)+2L),
      active_cells=nrow(cells),active_km=sum(ref$nodes$len[cells$node]),mean_sites=g$mean[total],
      sites_lo=g$lo[total],sites_hi=g$hi[total],sites_rhat=g$rhat[total],sites_mcse=g$mcse[total],
      mean_sites_on_common_40km_footprint=sum(cells$mean_count[cells$node %in% ref_nodes]),
      mark_spread=if(mode=="unit") 0 else get("mark_spread","mean"),
      mark_lo=if(mode=="unit") 0 else get("mark_spread","lo"),
      mark_hi=if(mode=="unit") 0 else get("mark_spread","hi"),
      mark_rhat=if(mode=="unit") NA_real_ else get("mark_spread","rhat"),
      mark_mcse=if(mode=="unit") NA_real_ else get("mark_spread","mcse"),
      predicted_zero=env$median[zero],zero_lo=env$lo[zero],zero_hi=env$hi[zero],
      predicted_positive_fraction=env$median[frac],predicted_high_villages=env$median[high],
      replicate_crps=unname(crps["value"]),crps_simulation_mcse=unname(crps["mcse"]),
      longitude_coef=get("theta.z_lon","mean"),latitude_coef=get("theta.z_lat","mean"),
      aridity_coef=get("beta.z_logAI","mean"),maximum_global_rhat=max(g$rhat),
      minimum_global_bulk_ess=min(g$bulk_ess),maximum_finite_village_rhat=max(v$rhat[is.finite(v$rhat)]),
      undefined_village_rhat=sum(!is.finite(v$rhat)),villages_above_101=sum(v$rhat>1.01,na.rm=TRUE),
      max_elapsed_seconds=max(runs$elapsed_s))
  }
}
tab <- do.call(rbind,rows)
tab$global_diagnostic <- ifelse(tab$maximum_global_rhat>1.05 | tab$minimum_global_bulk_ess<100,
                                "weak",ifelse(tab$maximum_global_rhat>1.01 | tab$minimum_global_bulk_ess<400,
                                              "caution","ok"))
write.csv(tab,file.path(grid,"grid_comparison.csv"),row.names=FALSE)
write.csv(do.call(rbind,effect_rows),file.path(grid,"grid_covariate_effects.csv"),row.names=FALSE)
# Keep other outcomes in the table. The plot compares one diagnostic across
# models against its observed value, without turning it into a model score.
png(file.path(grid,"range_comparison.png"),width=1740,height=790,res=160)
layout(matrix(c(1,1,1,2,3,4,5,5,5),3,byrow=TRUE),heights=c(.55,3.7,.35))
par(mar=rep(0,4),family="sans",fg="#223C40")
plot.new();text(.5,.74,"How the two ranges change the full-river fit",font=2,cex=1.55)
text(.5,.20,"All panels include river covariates, source coordinates and village aridity",cex=.97,col="#456065")
zero_col <- colorRampPalette(c("#ECF4F1","#78B4AB","#276D69","#124A4B"))(100)
mark_col <- colorRampPalette(c("#FAF0DC","#D7AC72","#A56843","#7C422D"))(100)
for(panel in 1:3){
  dat<-tab[tab$model==c("unit","two_point","two_point")[panel], ]
  value<-if(panel==3) dat$mark_spread else dat$predicted_zero
  scale<-if(panel==3) 1 else 550
  colours<-if(panel==3) mark_col else zero_col
  par(mar=c(3.2,3.3,2.6,.65));plot.new();plot.window(c(.4,3.6),c(.4,3.6),asp=1)
  for(k in seq_len(nrow(dat))){
    h<-match(dat$field_range_km[k],ranges);r<-match(dat$radius_km[k],radii)
    i<-pmax(1,pmin(100,round(value[k]/scale*99)+1))
    rect(r-.48,h-.48,r+.48,h+.48,col=colours[i],border=NA)
    label <- if(panel==3) sprintf("%.4f",value[k]) else format(value[k],trim=TRUE)
    if(dat$global_diagnostic[k]=="weak") label <- paste0(label,"*")
    text(r,h,label,
         cex=1.25,font=2,col=if(i>62) "white" else "#1D353B")
  }
  axis(1,at=1:3,labels=radii,tick=FALSE,cex.axis=.95)
  axis(2,at=1:3,labels=ranges,tick=FALSE,las=1,cex.axis=.95)
  mtext("Dispersal radius (km)",1,line=1.9,cex=.87)
  mtext("Field range (km)",2,line=1.9,cex=.87)
  mtext(c("Unit marks: zero surveys","Two marks: zero surveys","Estimated mark spread")[panel],
        3,line=1.1,font=2,cex=1.04)
}
par(mar=rep(0,4));plot.new()
caption <- "The first two panels show medians from 250 new landscapes per model. The surveys recorded zeros at 550 of 900 villages."
if(any(tab$global_diagnostic=="weak")) caption <- paste(caption, "An asterisk marks weak mixing.")
text(.5,.5,caption,cex=.86,col="#456065")
dev.off()

# Absolute weights show that the radius affects both peak and support.
distance <- seq(0,max(radii),length.out=801)
kernel <- vapply(radii,function(phi) 3/(2*phi)*odispers::wendland(distance/phi), numeric(length(distance)))
cols <- c("#267783","#B75B35","#595083")
png(file.path(grid,"kernel_profiles.png"),width=1440,height=780,res=170)
par(mar=c(4.5,5,4.5,1.3),family="sans",fg="#243C40",las=1)
plot(distance,kernel[,1],type="n",xlim=range(distance),ylim=c(0,.04),xaxs="i",yaxs="i",
     xlab="Straight-line distance from river site to village (km)",ylab="",axes=FALSE)
abline(h=seq(0,.04,.01),col="#DDE7E6",lwd=1)
for(j in seq_along(radii))lines(distance,kernel[,j],col=cols[j],lwd=3)
axis(1,at=seq(0,max(radii),10));axis(2,at=seq(0,.04,.01))
mtext("Exposure per unit-strength site",side=2,line=3.6,las=0)
title("Dispersal kernels for the three radii",font.main=2,col.main="#20383B",line=2.5)
mtext("Wendland weight with the package's line normalization",side=3,line=1.1,col="#48646A",cex=.88)
legend("topright",legend=paste(radii,"km radius"),col=cols,lwd=3,bty="n",cex=.9)
dev.off()

# Spatial pattern of replicated survey zeros for each unit-strength model.
# Each point uses the same 250-replicate Monte Carlo scale and village sample size.
data(sudan_borders,package="odispers")
land <- sf::st_geometry(sudan_borders)
pcols <- c("#E7EAE7","#AACCCC","#75A9A4","#397F7E","#185357")
png(file.path(grid,"null_spatial_comparison.png"),width=1710,height=1660,res=165)
layout(matrix(c(rep(1,3),2:10,rep(11,3)),5,byrow=TRUE),heights=c(.46,2.9,2.9,2.9,.43))
par(mar=rep(0,4),family="sans",fg="#223C40")
plot.new();text(.5,.72,"Where the unit-strength model predicts zero-positive surveys",font=2,cex=1.44)
text(.5,.20,"Estimated probability from 250 new landscapes per setting, at the 900 actual survey sites",cex=.93,col="#4C6167")
for(h in ranges) for(r in radii) {
  case <- case_directory(r, h, "unit")
  sim <- readRDS(file.path(case,"null","null_simulations.rds"))
  prob <- colMeans(sim$y==0)
  par(mar=c(.05,.2,1.5,.2));plot.new();plot.window(c(21.5,38.8),c(3.2,22.5),asp=1/cos(12*pi/180))
  plot(land,add=TRUE,col="#F5F6F4",border="#A6B2B2",lwd=.5)
  i <- order(prob)
  points(ref$villages$long[i],ref$villages$lat[i],col=pcols[1L+findInterval(prob[i],c(.1,.25,.5,.75))],
         pch=16,cex=.58)
  mtext(sprintf("Field %d km   ·   Kernel %d km",h,r),side=3,line=.35,font=2,cex=.95)
}
par(mar=rep(0,4));plot.new()
legend("center",legend=c("<10%","10–<25%","25–<50%","50–<75%","75% or higher"),
       col=pcols,pch=16,horiz=TRUE,bty="n",cex=.95)
dev.off()
