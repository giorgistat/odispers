# Map actual outcomes, spatial holdout predictions and the first fresh landscape.
# Run with the results directory of the simple aridity-dependent occurrence model.
library(sf)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
out <- args[1]
visits <- read.csv(file.path(out,"folds.csv"))
dir <- file.path(out,"occurrence")
pred <- do.call(rbind,lapply(seq_len(5),function(f) {
  x <- read.csv(file.path(dir,paste0("fold",f),"predictions.csv"))
  stopifnot(all(x$fold==f),setequal(x$village,which(visits$fold==f)))
  x
}))
pred <- pred[order(pred$village),]
replicates <- readRDS(file.path(dir,"fold0","replications.rds"))
stopifnot(identical(pred$village,visits$village),identical(pred$y,visits$y),
          identical(replicates$test,visits$village), nrow(replicates$y)==600L)
zeros <- rowSums(replicates$y==0)
interval <- round(quantile(zeros,c(.025,.975)))
data(sudan_borders,package="odispers")
data(sudan_rivers,package="odispers")
projection <- "+proj=aeqd +lat_0=13 +lon_0=30 +datum=WGS84 +units=km +no_defs"
land <- st_transform(st_geometry(sudan_borders),projection)
rivers <- st_transform(st_geometry(sudan_rivers[sudan_rivers$DIS_AV_CMS>=450,]),projection)
rivers <- st_intersection(rivers,st_union(land))
points <- st_as_sf(visits,coords=c("long","lat"),crs=4326)
xy <- st_coordinates(st_transform(points,projection))
bbox <- st_bbox(land)
ink <- "#203937"; muted <- "#506762"; paper <- "#FFFFFF"
zero_col <- "#086664"; positive_col <- "#B9633B"
prob_cols <- c("#D9DFCE","#AECAB1","#79AD96","#3D897B",zero_col)
file <- file.path(out,"zero_positive_villages.png")
png(file,width=2380,height=1340,res=195,bg=paper,type="cairo")
layout(matrix(c(1,1,1,2,3,4,5,6,7,8,8,8,9,9,9),5,byrow=TRUE),
       heights=c(.94,4.13,.42,.77,.46))
par(mar=rep(0,4),family="sans",fg=ink)
new_canvas <- function() {
  plot.new()
  plot.window(c(0,1),c(0,1),xaxs="i",yaxs="i")
}
new_canvas()
text(.5,.72,"Zero-positive village surveys on the river network",font=2,col=ink,cex=2.05)
text(.5,.26,"Sudan and South Sudan  |  900 survey locations  |  Aridity-dependent occurrence model",col=muted,cex=1.13)

panel <- function(title,subtitle,kind) {
  par(mar=c(.1,.15,2.7,.15),xaxs="i",yaxs="i")
  plot.new()
  plot.window(c(bbox["xmin"]-30,bbox["xmax"]+30),c(bbox["ymin"]-35,bbox["ymax"]+35),asp=1)
  plot(land,add=TRUE,col="#F6F4EC",border="#B7BDB5",lwd=.7)
  plot(rivers,add=TRUE,col="#A6C6C4",lwd=.56)
  if(kind=="heldout") {
    bin <- 1L+findInterval(pred$prob_zero,c(.2,.4,.6,.8))
    ord <- order(bin)
    points(xy[ord,,drop=FALSE],pch=21,cex=.64,bg=prob_cols[bin[ord]],col="#667D71",lwd=.26)
  } else {
    z <- if(kind=="actual") visits$y==0 else replicates$y[1,]==0
    points(xy[z,,drop=FALSE],pch=21,cex=.59,bg=zero_col,col="white",lwd=.25)
    points(xy[!z,,drop=FALSE],pch=23,cex=.62,bg=positive_col,col="white",lwd=.26)
  }
  mtext(title,3,line=1.12,font=2,cex=1.22,col=ink)
  mtext(subtitle,3,line=-.02,cex=.96,col=muted)
}
panel("Observed surveys",sprintf("%d villages recorded no positives",sum(visits$y==0)),"actual")
panel("Held-out predictions",sprintf("%.0f zero-positive villages expected",sum(pred$prob_zero)),"heldout")
panel("One simulated survey",sprintf("%d villages recorded no positives",zeros[1]),"sim")

par(mar=rep(0,4))
outcome_legend <- function() {
  new_canvas()
  legend("center",c("Zero positives","One or more"),pch=c(21,23),pt.bg=c(zero_col,positive_col),
         col=c(zero_col,positive_col),pt.cex=1.0,horiz=TRUE,bty="n",cex=.9,x.intersp=.55,text.col=muted)
}
outcome_legend()
new_canvas()
width <- .156
for(i in seq_along(prob_cols)) {
  left <- .11+(i-1)*width
  rect(left,.55,left+width,.88,col=prob_cols[i],border=NA)
}
text(.11,.26,"0%",adj=0,cex=.83,col=muted);text(.89,.26,"100%",adj=1,cex=.83,col=muted)
text(.5,.26,"Probability of zero positives",cex=.78,col=muted)
outcome_legend()

new_canvas()
segments(.06,.95,.94,.95,col="#D7DFDA",lwd=1)
text(.5,.59,sprintf("Across 600 fresh landscapes, the median was %d zero-positive villages  |  95%% interval %d–%d",
                     as.integer(median(zeros)),interval[1],interval[2]),cex=1.22,font=2,col=ink)
text(.5,.14,"The simulation shown is the first draw. Each landscape redraws the river sites, site strengths, village states and survey counts.",
     cex=.93,col=muted)
new_canvas()
text(.5,.44,"Held-out probabilities come from five spatial folds. Rivers shown have mean discharge of at least 450 m³/s.",cex=.84,col=muted)
dev.off()
cat(normalizePath(file),"\n")
