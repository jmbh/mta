
# computes within cluster and between cluster variance/distance & average trajectory in each cluster

clusterstat <- function(data_rescaled, #data
                       i.id, 
                       i.xyt
) {
  
  cl <- data_rescaled$cl

  # eucl distance function
  edist <- function(x1,y1,x2,y2) { sum(sqrt((x1-x2)^2+(y1-y2)^2)) }
  
  # calc mean traj
  l_mtraj <- list()
  for(i in unique(cl)) {
    l_mtraj[[i]] <- ddply(data_rescaled[cl==i,], "t", function(x) {
      out <- cbind(mean(x$x),
                   mean(x$y))
      return(out)
    })[,-1]
    colnames(l_mtraj[[i]]) <- c('x', 'y')
  }
  
  # calc between cluster distance
  dims <- length(l_mtraj)
  m <- matrix(NA,dims,dims)
  for(i in 1:dims) {
    for(j in i:dims) {
      m[i,j] <-   edist(l_mtraj[[i]][,1], l_mtraj[[i]][,2], 
                        l_mtraj[[j]][,1], l_mtraj[[j]][,2])
    }
  }
  md_between <- mean(m[upper.tri(m)])
  
  # calc within cluster distance
  d_trial <- ddply(data_rescaled, c(i.id,"cl"), function(xw) {
    meantraj <- l_mtraj[[xw$cl[1]]]
    ed <- edist(meantraj[,1], meantraj[,2], xw$x, xw$y)
    return(ed)
  })
  md_within <- mean(aggregate(V1~cl, data=d_trial, mean)$V1)
  
  outlist <- list("Mdist.between"=md_between, "Mdist.within"=md_within, "M_traj"=l_mtraj)
  
  return(outlist)
}

