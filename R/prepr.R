
prepr <- function(data, # data frame with x,y,t and flagging variables
                  i.xyt, # names of columns for x,y,t in this order
                  i.id, # names of variables if id variables
                  type = "time", # also: spatial 
                  steps = 101, 
                  start2zero = TRUE, #
                  stretch = NA,
                  takeAllvar = FALSE) # same form as layout
  
{
  
  # +++ input checks +++
  stopifnot(class(data[,i.xyt[1]])=='numeric' | class(data[,i.xyt[1]])=='integer')
  stopifnot(class(data[,i.xyt[2]])=='numeric' | class(data[,i.xyt[2]])=='integer')
  stopifnot(class(data[,i.xyt[3]])=='numeric' | class(data[,i.xyt[3]])=='integer')
  if(sum(is.na(data)>0)) {  stop("No missing values permitted.")  }
  
  # +++ prepare data +++
  dat <- data.frame(data[,c(i.xyt, i.id)])
  colnames(dat)[1:3] <- cn <-c('x', 'y', 't')
  
  # +++ set starting point of each trajectory to (0,0) +++
  if(start2zero==TRUE) {
    dat <- ddply(dat, colnames(dat)[-(1:3)], function(x) {
      x[,1] <- x[1,1] - x[,1] # set x-start to zero
      x[,2] <- x[1,2] - x[,2] # set y-start to zero
      x
    })
  }
  
  # +++ normalize wrt time +++
  if(type=='time') {
    
    n_dat <- ddply(dat, i.id, function(traj) {
      rnorm <- (traj$t-traj$t[1]) / max((traj$t-traj$t[1])) * steps
      a.x <- approx(rnorm, traj$x,  xout = 0:(steps-1), method = "linear") 
      a.y <- approx(rnorm, traj$y,  xout = 0:(steps-1), method = "linear") 
      cbind(a.x$y, a.y$y, 0:(steps-1))
    })
    colnames(n_dat) <- c(i.id, cn)
    dat <- n_dat[,c(cn, i.id)]
    
  }
  
  # +++ normalize wrt space +++
  if(type=='spatial') {
    
    stop("not impemented yet ...")
    
  }
  
  # +++ stretch +++
  
  if(is.na(stretch)[1]==FALSE) {
    
    dat_str <- ddply(dat, i.id, function(traj) {
      # starting point
      X <- traj$x - traj$x[1]; X <- X + stretch$start[1]
      Y <- traj$y - traj$y[1]; Y <- Y + stretch$start[2]
      # end point
      X <- (X / abs(X[length(X)])) * abs(stretch$left[1])     
      Y <- (Y / abs(Y[length(Y)])) * abs(stretch$left[2])     
      t <- traj$t
      cbind(X,Y,t)
    })
    colnames(dat_str) <- c(i.id, cn)
    dat <- dat_str[,c(cn, i.id)]
    
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
    
  }
  
  # reorder columns
  dat <- dat[,c(i.id, i.xyt, namesv)]
  
  # output
  call <- list('i.xyt'=i.xyt, 'i.id'=i.id,  'type'=type, 
               'steps'=steps, 'start2zero'=start2zero, 'stretch'=stretch)
  outlist <- list('call'=call, 'data'=dat)
  
  class(outlist) <- 'mta'
  
  return(outlist) 
}










