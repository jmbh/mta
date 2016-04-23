##
# interpolates n equally distant time points on each trajectory
##


spatialRescale <- function(data, i.id, i.xyt, nResc) {
  
  # input check
  mov_check_d <- ddply(data, c(i.id), function(x) { !(sd(x$x)>0 & sd(x$y)>0)  })
  mov_check   <- mov_check_d[, ncol(mov_check_d)]
  
  if(sum(mov_check)>0) {
    print(mov_check_d[mov_check,])
    stop("The above trials have zero movement. Movement is necessary for spatial Rescaling.")
    }
  
  dat   = dlply(data,i.id,function(x) as.matrix(x[,i.xyt]))
  
  dat_r = spatialRescale_c(dat, nResc)
  front = ddply(data[,i.id],i.id,function(x) matrix(1:nResc,ncol=1))
  data_resc = data.frame(front, dat_r[,1], dat_r[,2])
  data_resc = data_resc[,c(1,2,4,5,3)]
  names(data_resc)[3:5] = i.xyt
  
  return(data_resc)


  }

