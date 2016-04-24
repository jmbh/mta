

derivatives <- function(data,
                        i.xyt = c('x','y','t'), 
                        i.id  = c('ptp','trial'))
  
{
 dist = function(x,y) abs(-1.5*x-y)/sqrt(3.25)  
 dists = dlply(data,c(i.id),function(x) dist(x$x,x$y))
 MAD = sapply(dists,max)
 AAD = sapply(dists,mean)
 return(list('MAD'=MAD,'AAD'=AAD,'dists'=dists))  
}