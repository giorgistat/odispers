# Compare the location of observed and predicted zero surveys in the held-out villages.
library(sf)
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
root <- args[1]
folds <- read.csv(file.path(root,"folds.csv"))
models <- c("baseline","river_climate","occurrence","both","flexible_occurrence","survey_design")
held <- lapply(models,function(model) {
  rows <- do.call(rbind,lapply(1:5,function(f) {
    read.csv(file.path(root,model,paste0("fold",f),"predictions.csv"))
  }))
  rows[order(rows$village),]
})
names(held) <- models
reliability <- list()
by_block <- list()
for (model in models) {
  p <- held[[model]]$prob_zero
  # Equal-size probability groups show whether the predicted event rate calibrates.
  bin <- ceiling(rank(p,ties.method="first")/90)
  reliability[[model]] <- do.call(rbind,lapply(split(seq_len(900),bin),function(i) {
    data.frame(model=model,bin=bin[i[1]],n=length(i),mean_prediction=mean(p[i]),observed=mean(folds$y[i]==0))
  }))
  by_block[[model]] <- do.call(rbind,lapply(split(seq_len(900),folds$block),function(i) {
    data.frame(model=model,block=folds$block[i[1]],n=length(i),long=mean(folds$long[i]),lat=mean(folds$lat[i]),
               expected=sum(p[i]),observed=sum(folds$y[i]==0))
  }))
}
write.csv(do.call(rbind,reliability),file.path(root,"heldout_reliability.csv"),row.names=FALSE)
write.csv(do.call(rbind,by_block),file.path(root,"heldout_spatial_blocks.csv"),row.names=FALSE)
data(sudan_borders,package="odispers")
land <- st_geometry(sudan_borders)
cols <- c("#DBE3E0","#ABC8C3","#74A49F","#367B79","#174D50")
png(file.path(root,"heldout_locations.png"),width=1740,height=880,res=155)
layout(matrix(c(1,1,1,2,3,4,5,5,5),3,byrow=TRUE),heights=c(.6,4.8,.55))
par(mar=rep(0,4),family="sans",fg="#233D41")
plot.new();text(.5,.74,"Where the zero-positive surveys occur",cex=1.53,font=2)
text(.5,.18,"Five spatial holdouts  ·  Same full river network and survey denominators",cex=.95,col="#4D6568")
for(model in c("observed","baseline","flexible_occurrence")) {
  par(mar=c(.1,.2,1.7,.2));plot.new();plot.window(c(21.1,39.2),c(3.2,22.5),asp=1/cos(12*pi/180))
  plot(land,add=TRUE,col="#F7F7F3",border="#AFBCB9",lwd=.6)
  p <- if(model=="observed") as.numeric(folds$y==0) else held[[model]]$prob_zero
  order_p <- order(p)
  col <- if (model=="observed") ifelse(p==1,cols[5],"#B8704D") else cols[1+findInterval(p,c(.2,.4,.6,.8))]
  points(folds$long[order_p],folds$lat[order_p],pch=16,col=col[order_p],cex=.65)
  title <- switch(model,observed="Survey outcomes",baseline="Original LGCP",flexible_occurrence="Aridity-dependent occurrence")
  mtext(title,3,line=.8,font=2,cex=1.05)
}
par(mar=rep(0,4));plot.new()
legend("center",c("Survey had positives","Survey had zero positives","Predicted <20%","20–40%","40–60%","60–80%","≥80%"),
       col=c("#B8704D",cols[5],cols),pch=16,bty="n",horiz=TRUE,cex=.8)
dev.off()
