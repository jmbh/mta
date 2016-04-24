

derivatives <- function(data,
                        i.xyt, 
                        i.id)
  
{
 dist = function(x,y) abs(-1.5*x-y)/sqrt(3.25)  
 dists = dlply(data,c(i.id),function(x) dist(x$x,x$y))
 MAD = sapply(dists,max)
 AAD = sapply(dists,mean)
  
  
}