# Maps from the full fits and simulations. Run after full_river_summary.R
# and full_river_null.R (unit in null/, two_point in marked_predictive/).
library(odispers)
library(sf)
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 4L))
  stop("supply inputs.rds and the full analysis directory, optionally also fit and marked predictive subdirectories")
fit_dir <- if (length(args) == 4L) args[3] else "observed"
marked_dir <- if (length(args) == 4L) args[4] else "marked_predictive"
a <- readRDS(args[1])
out <- args[2]
phi <- if (is.null(a$phi)) 40 else a$phi
v <- a$villages
null <- readRDS(file.path(out, "null", "null_simulations.rds"))
marked <- readRDS(file.path(out, marked_dir, "marked_simulations.rds"))
data(sudan_borders, package = "odispers")
data(sudan_rivers, package = "odispers")
major <- st_geometry(sudan_rivers[sudan_rivers$DIS_AV_CMS >= 250, ])
ink <- "#20383B"
muted <- "#48646A"
colours <- c("#B6BFC1", "#82BCC4", "#29828A", "#CC823A", "#913E36")
point_colour <- function(p) {
  where <- ifelse(p == 0, 1L, 2L + (p >= .1) + (p >= .25) + (p >= .5))
  colours[where]
}

png(file.path(out, "village_prevalence.png"), width = 1820, height = 950, res = 155)
layout(matrix(c(1,1,1, 2,3,4, 5,5,5), 3, byrow = TRUE), heights = c(.63,4.6,.55))
par(mar = rep(0,4), family = "sans", fg = ink)
plot.new()
text(.5,.72, "Village prevalence across Sudan and South Sudan", cex=1.63, font=2, col=ink)
text(.5,.24,"Observed surveys and the first simulated survey from each fitted model", cex=1.05, col=muted)
outcomes <- list(v$y, null$y[1, ], marked$y[1, ])
titles <- c("Observed", "Unit-strength simulation", "Two-strength simulation")
for (j in seq_along(outcomes)) {
  par(mar = c(.1,.5,2.4,.5))
  plot.new()
  plot.window(xlim=c(21.5,38.8), ylim=c(3.2,22.5), asp=1/cos(12*pi/180))
  plot(st_geometry(sudan_borders), add=TRUE, col="#F3F5F4", border="#A1AFAF", lwd=.8)
  plot(major, add=TRUE, col="#C2D8DE", lwd=.55)
  y <- outcomes[[j]]
  prop <- y / v$n
  # Positives on top so zeros cannot hide clustered high-prevalence villages.
  o <- order(prop)
  points(v$long[o], v$lat[o], pch=21, cex=.62, bg=point_colour(prop[o]),
         col=grDevices::adjustcolor("#344F51", .48), lwd=.3)
  mtext(titles[j], line=1.25, font=2, cex=1.13, col=ink)
  mtext(sprintf("%d of 900 villages recorded zero positives",sum(y==0)), line=.02, cex=.85, col=muted)
}
par(mar=rep(0,4)); plot.new()
legend("center", legend=c("Zero","Below 10%","10–<25%","25–<50%","50% or higher"),
       pt.bg=colours, col="#527071", pch=21, pt.cex=1.4, horiz=TRUE, bty="n", cex=.99)
dev.off()

# Pair the two marks on the *same* first null site pattern at the first
# survey. Include sites outside the displayed window up to one kernel radius.
land <- readRDS(file.path(out, "null", "first_landscape.rds"))
params <- read.csv(file.path(out,fit_dir,"two_point_parameters.csv"))
d <- params$median[params$parameter=="mark_spread"]
set.seed(801106)
H <- rbinom(length(land$N), land$N, .5)
mass <- (land$N-H)*(1-d) + H*(1+d)
xy <- st_coordinates(st_transform(v, attr(a$nodes,"crs")))
center <- xy[1, ]
xs <- seq(center[1]-40, center[1]+40, length.out=101)
ys <- seq(center[2]-40, center[2]+40, length.out=101)
sites <- which(land$N>0 & abs(a$nodes$x-center[1]) <= 40+phi & abs(a$nodes$y-center[2]) <= 40+phi)
grid <- expand.grid(x=xs,y=ys)
E <- matrix(0, nrow(grid), 2)
for (i in sites) {
  r <- sqrt((grid$x-a$nodes$x[i])^2 + (grid$y-a$nodes$y[i])^2)
  k <- wendland(r/phi)/(phi*2/3)
  E[,1] <- E[,1] + land$N[i]*k
  E[,2] <- E[,2] + mass[i]*k
}
near <- abs(a$nodes$x-center[1]) <= 44 & abs(a$nodes$y-center[2]) <= 44
rivers <- st_geometry(st_transform(sudan_rivers[sudan_rivers$HYRIV_ID %in% a$nodes$HYRIV_ID[near], ],
                                    attr(a$nodes,"crs")))
inside <- sites[abs(a$nodes$x[sites]-center[1]) <= 40 & abs(a$nodes$y[sites]-center[2]) <= 40]
palette <- colorRampPalette(c("#FAFAEF","#F0DAAA","#D5AF79","#C07C51","#964B3F"))(80)
lim <- c(0, max(E))
png(file.path(out,"kernel_zoom.png"), width=1560, height=900, res=155)
layout(matrix(c(1,1,2,3,4,4), 3, byrow=TRUE), heights=c(.55,4,.65))
par(mar=rep(0,4), family="sans",fg=ink);plot.new()
text(.5,.72, "How marks change exposure around the same breeding sites", font=2,cex=1.45,col=ink)
text(.5,.21, sprintf("80 km window around the first survey  ·  %g km dispersal radius",phi),cex=.99,col=muted)
for (j in 1:2) {
  par(mar=c(.5, .5, 2, .5), pty="s")
  plot.new();plot.window(range(xs),range(ys),asp=1,xaxs="i",yaxs="i")
  image(xs,ys,matrix(E[,j],length(xs)),col=palette,zlim=lim,add=TRUE,useRaster=TRUE)
  plot(rivers,add=TRUE,col="#54828C",lwd=.6)
  strength <- if (j==1) land$N else mass
  points(a$nodes$x[inside],a$nodes$y[inside], cex=.2+.73*sqrt(strength[inside]),
         bg=adjustcolor("#243C46", .78),col="#FFFFFF",pch=21,lwd=.45)
  points(xy[,1],xy[,2],pch=3,cex=.8,col="#182A30",lwd=1.2)
  segments(xs[5],ys[5],xs[5]+10,ys[5],lwd=2)
  text(xs[5]+5,ys[5]+2.2,"10 km",cex=.7)
  box(col="#94A2A4")
  mtext(if(j==1) "Each site has strength 1" else sprintf("Site strengths %.4f and %.4f",1-d,1+d),
        font=2,cex=1.05,line=.5,col=ink)
}
par(mar=rep(0,4),pty="m");plot.new();plot.window(c(0,1),c(0,1))
rect(seq(.12,.44,length.out=81)[1:80],.34,seq(.12,.44,length.out=81)[-1],.59,col=palette,border=NA)
text(.12,.14,"0",cex=.73);text(.44,.14,formatC(lim[2],format="f",digits=2),cex=.73)
text(.28,.82,"Kernel exposure",cex=.87,col=ink)
points(.57,.49,pch=21,cex=1,bg="#243C46",col="white")
text(.595,.49,"Site (size = strength)",adj=0,cex=.83)
points(.82,.49,pch=3);text(.84,.49,"Survey",adj=0,cex=.83)
dev.off()

# The complete null realization, including cells outside the survey footprint.
# Draw in the same projected kilometre coordinates used by the kernels.
borders <- st_geometry(st_transform(sudan_borders, attr(a$nodes,"crs")))
limits <- list(x=range(a$nodes$x),y=range(a$nodes$y))
intensity_breaks <- c(-Inf,.01,.025,.05,.1,.25,Inf)
rate_col <- c("#DAE4E3","#A3C8C3","#61A59B","#2F807C","#C58A43","#974F39")
bin <- findInterval(land$rate_per_km, intensity_breaks)
source <- which(land$N>0)
png(file.path(out,"null_sites.png"),width=1620,height=970,res=160)
layout(matrix(c(1,1,2,3,4,4),3,byrow=TRUE),heights=c(.55,4.6,.48))
par(mar=rep(0,4),family="sans",fg=ink)
plot.new();text(.5,.74,"One null realization on the full river graph",font=2,cex=1.5,col=ink)
text(.5,.22,"The intensity determines the expected number of sites in each river cell",cex=.98,col=muted)
for(j in 1:2) {
  par(mar=c(.1,.15,2.2,.15))
  plot.new(); plot.window(limits$x, limits$y, asp=1)
  plot(borders, add=TRUE, col="#FBFAF6", border="#AAB6B6", lwd=.65)
  if(j==1) {
    points(a$nodes$x,a$nodes$y,pch=".",col=rate_col[bin])
  } else {
    points(a$nodes$x,a$nodes$y,pch=".",col="#CFDDDB")
    points(a$nodes$x[source],a$nodes$y[source],pch=16,cex=.17+.10*sqrt(land$N[source]),
           col=adjustcolor("#254F50", .76))
  }
  mtext(c("Conditional site intensity","Realized breeding sites")[j],line=1.15,font=2,cex=1.08,col=ink)
  mtext(c("Sites per kilometre",sprintf("%s sites on %s occupied cells",format(sum(land$N),big.mark=","),
                                        format(length(source),big.mark=",")))[j],
        line=.03,cex=.79,col=muted)
}
par(mar=rep(0,4));plot.new()
legend("center",legend=c("<0.01","0.01–0.025","0.025–0.05","0.05–0.10","0.10–0.25","≥0.25 sites/km"),
       fill=rate_col,border=NA,horiz=TRUE,bty="n",cex=.8)
dev.off()
