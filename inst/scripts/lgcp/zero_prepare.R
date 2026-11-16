# Freeze the spatial folds and climate predictors before fitting the zero-survey comparison.
library(odispers)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
a <- readRDS(args[1])
out <- args[2]
dir.create(out, showWarnings=FALSE, recursive=TRUE)
v <- a$villages
xy <- sf::st_coordinates(sf::st_transform(v, attr(a$nodes, "crs")))
block <- paste(floor(xy[,1]/100), floor(xy[,2]/100), sep="_")
# Assign whole 100-km squares to five equally sized folds. Use no outcomes.
sizes <- sort(table(block), decreasing=TRUE)
load <- integer(5)
block_fold <- setNames(integer(length(sizes)), names(sizes))
for (b in names(sizes)) {
  f <- which.min(load)
  block_fold[b] <- f
  load[f] <- load[f] + sizes[b]
}
a$spatial_fold <- unname(block_fold[block])
v$climate <- (v$AI - .4)/.2
a$villages <- v
data(sudan_aridity, package="odispers")
r <- terra::unwrap(sudan_aridity)
points <- sf::st_as_sf(as.data.frame(a$nodes)[,c("x","y")],coords=c("x","y"),crs=attr(a$nodes,"crs"))
lonlat <- sf::st_coordinates(sf::st_transform(points,4326))
river_ai <- terra::extract(r,lonlat,method="bilinear")[,1]/10000
active_nodes <- unique(a$pairs$j[a$pairs$r<a$phi])
stopifnot(!anyNA(river_ai[active_nodes]))
# Unmeasured river cells outside every kernel get the mean measured dry-side effect.
# These cells cannot affect a fitted or simulated village outcome.
dry <- pmin(river_ai - .45, 0)/.2
missing <- !is.finite(dry)
dry[missing] <- weighted.mean(dry[!missing],a$nodes$len[!missing])
a$nodes$river_dry <- dry
saveRDS(a,file.path(out,"inputs.rds"),compress=FALSE)
folds <- data.frame(village=seq_len(nrow(v)),long=v$long,lat=v$lat,AI=v$AI,n=v$n,y=v$y,
                    block=block,fold=a$spatial_fold)
write.csv(folds,file.path(out,"folds.csv"),row.names=FALSE)
cat("Cells lacking aridity:",sum(missing),"; within reach:",sum(missing[active_nodes]),"\n")
print(table(a$spatial_fold))
