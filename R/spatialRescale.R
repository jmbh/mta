##
# interpolates n equally distant time points on each trajectory
##


spatialRescale1Traj <- function(x, y, nResc) {
  
  # get distance
  x2 <- c(0,x[-length(x)])
  y2 <- c(0,y[-length(x)])
  xdiff <- sqrt((x2-x)^2 + (y2-y)^2)
  cs <- cumsum(xdiff) #cumulative distance
  
  #get equally spaced x steps:
  xsteps <- seq(0,max(cs), length=nResc)
  x_rescaled <- y_rescaled <- numeric(nResc)
  
  for(i in 1:nResc) {
    #get two adjacent x
    check <- xsteps[i]>cs
    sum_c <- sum(check)
    sum_c[sum_c==0]<- 1
    #weigted sum of the x
    w1 <- abs(xsteps[i]-cs[sum_c])
    w2 <- abs(xsteps[i]-cs[sum_c+1])
    #save
    x_rescaled[i] <- x[sum_c] * w2/(w1+w2) + x[sum_c+1] * w1/(w1+w2)
    y_rescaled[i] <- y[sum_c] * w2/(w1+w2) + y[sum_c+1] * w1/(w1+w2)
    if(i==1 & is.na(x_rescaled[i])==TRUE) { x_rescaled[i]<-0 }
    if(i==1 & is.na(y_rescaled[i])==TRUE) { y_rescaled[i]<-0 }
  }
  out <- cbind(x_rescaled, y_rescaled, 1:nResc)
  colnames(out) <- c("x", "y", "t")
  return(out)
} 



spatialRescale <- function(data, i.id, i.xyt, nResc) {
  
  # input check
  mov_check_d <- ddply(data, c(i.id), function(x) { !(sd(x$x)>0 & sd(x$y)>0)  })
  mov_check <- mov_check_d[, ncol(mov_check_d)]
  if(sum(mov_check)>0) {
    print(mov_check_d[mov_check,])
    stop("The above trials have zero movement. Movement is necessary for spatial Rescaling.")
  }
  
  for(ptp in unique(data$id.ptp)) {
    for(trial in unique(data$id.trial[data$id.ptp==ptp])) {
      x <- subset(data, id.ptp==ptp & id.trial == trial)
      spatialRescale1Traj(x[,i.xyt[1]], x[,i.xyt[2]], nResc)
    }
  }
  
  data_resc <- ddply(data, i.id, function(x) { 
    out <- spatialRescale1Traj(x[,i.xyt[1]], x[,i.xyt[2]], nResc) 
    return(out)
    })
  
  return(data_resc)
  
}




