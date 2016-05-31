
# calculates chosen option from data

getside <- function(dat, i.id = c('ptp','trial')) {
  
  dat$side <- 0
  mean.start <- mean(ddply(dat[,c(i.id,'x')], c(i.id), function(x) { 
    return(x$x[1])
    })$V1)
  
  side <- ddply(dat[,c(i.id,'side','x')], c(i.id), function(x) {
    if(x$x[nrow(x)]>mean.start) { x$side <- 1  } 
    return(x)
  })
  
  return(side$side)
  }

