

prepr <- function(data, # data frame with x,y,t and flagging variables
                  i.xyt = c('x','y','t'), # names of columns for x,y,t in this order
                  i.id  = c('ptp','trial'), # names of variables if id variables
                  type = "time", # also: spatial 
                  steps = 101, 
                  start2zero = TRUE, 
                  fliponeside = TRUE, 
                  stretch = list('start' = c(0,0),'left' = c(-1,1.5)),
                  takeAllvar = FALSE) 
  
{
  
  # +++ input checks +++
  stopifnot(class(data[,i.xyt[1]])=='numeric' | class(data[,i.xyt[1]])=='integer')
  stopifnot(class(data[,i.xyt[2]])=='numeric' | class(data[,i.xyt[2]])=='integer')
  stopifnot(class(data[,i.xyt[3]])=='numeric' | class(data[,i.xyt[3]])=='integer')
  if(sum(is.na(data[,c(i.xyt,i.id)])>0)) {  stop("No missing values permitted.")  }
  
  # +++ prepare data +++
  dat <- data.frame(data[,c(i.id,i.xyt)])
  
  # +++ set starting point of each trajectory to (0,0) +++
  if(start2zero==TRUE) {
    dat <- ddply(dat, i.id, function(x) {
      x[,i.xyt[1]] <- x[,i.xyt[1]] - x[1,i.xyt[1]] # set x-start to zero
      x[,i.xyt[2]] <- x[,i.xyt[2]] - x[1,i.xyt[2]] # set y-start to zero
      x
    })
  }
  
  # +++ calculate side of chosen box+++
  dat$choice <- getside(dat, i.id)
  
  # +++ flip trajectories to left side ++++
  if(fliponeside & start2zero) dat$x = ifelse(dat$choice == 1, dat$x*-1, dat$x)
  if(fliponeside & !start2zero) stop('Flipping should only be done for start2zero == T')
  
  # +++ normalize wrt time +++
  if(type=='time') {
    n_dat <- ddply(dat, i.id, function(traj) {
      rnorm <- (traj$t-traj$t[1]) / max((traj$t-traj$t[1])) * steps
      a.x <- approx(rnorm, traj$x,  xout = 0:(steps-1), method = "linear") 
      a.y <- approx(rnorm, traj$y,  xout = 0:(steps-1), method = "linear") 
      cbind(a.x$y, a.y$y, 0:(steps-1))
      })
    colnames(n_dat) <- c(i.id, i.xyt)
    dat <- n_dat[,c( i.id, i.xyt)]
    }
  
  # +++ normalize wrt space +++
  if(type=='spatial') {
    if(is.na(steps)) {stop("Please specify the number of points to be interpolated on each trajectory (steps)")}
    dat <- spatialRescale(data, i.id, i.xyt, steps)
    }
  
  # +++ stretch +++
    if(is.na(stretch)[1]==FALSE) {
    dat_str <- ddply(dat, i.id, function(traj) {
      # starting point
      X <- traj$x - traj$x[1]; X <- X + stretch$start[1]
      Y <- traj$y - traj$y[1]; Y <- Y + stretch$start[2]
      # end point
      X <- (X / X[length(X)]) * stretch$left[1]     
      Y <- (Y / Y[length(Y)]) * stretch$left[2]     
      t <- traj$t
      cbind(X,Y,t)
      })
    colnames(dat_str) <- c(i.id, i.xyt)
    dat <- dat_str[,c(i.id, i.xyt)]
    }
  
  # +++ add aux variables if specified +++
  namesv <- names(data)[!names(data) %in% c(i.id, i.xyt)]
    if(takeAllvar==TRUE & !is.null(namesv)) {
    aux_vars <- ddply(data, i.id, function(x) {
      namesv <- names(x)[!names(x) %in% c(i.id, i.xyt)]
      nv <- length(namesv)
      dat_aux <- x[,!names(x) %in% c(i.id, i.xyt)]
      if(nv==1) {
        dat_aux_1r <- dat_aux[1]  
        } else {
        dat_aux_1r <- unlist(dat_aux[1,])
        }
      m <- as.data.frame(matrix(rep(dat_aux_1r, times=steps), steps, nv, byrow = TRUE))
      names(m) <- namesv
      return(m)
      })
    dat[,namesv] <- aux_vars[, ! names(aux_vars) %in% i.id]
    dat <- dat[,c(i.id, i.xyt, other, namesv)]
    }
  
  # output
  call <- list('i.xyt'=i.xyt, 'i.id'=i.id,  'type'=type, 
               'steps'=steps, 'start2zero'=start2zero, 'stretch'=stretch)
  outlist <- list('call'=call, 'data'=dat)
  
  class(outlist) <- 'mta'
  
  return(outlist) 
}







