

gap_statistic <- function(dist, # p x p distance matrix of
                          kseq, #sequence of ks to be checked 
                          steps, # number of points on synthetic curves
                          lambda = .7, 
                          bezier = NULL,
                          xcor = c(0,1,-1), # x-coordinates of start ,nonselected, selected box inthatorder
                          ycor = c(0,1.5,1.5) # y-coordinates
)
  
{
  
  # ----- aux function -----
  
  rBezier = function(x, y, w = 1, resol = 100) { 
    
    t = seq(0,1,length = resol)           #Resolution
    B = cbind( (1-t)^2, 2*t*(1-t),t^2 )   #Bernstein
    
    bx = ( B[,1] * x[1] + w * B[,2] * x[2] + B[,3] * x[3] ) / 
      ( B[,1]        + w * B[,2]        + B[,3]) 
    by = ( B[,1] * y[1] + w * B[,2] * y[2] + B[,3] * y[3] ) / 
      ( B[,1]        + w * B[,2]        + B[,3]) 
    
    return(data.frame('x'=bx,'y'=by))
  }
  
  # ----- generate synthetic data -----
  
  n <- dim(dist)[1]
  
  # create curves
  l_curves <- vector('list', length=n)
  set.seed(1)
  l_curves <- lapply(l_curves, function(s) {
    bz <-  rBezier(x,y,w=rexp(1,.7), resol=steps)
    return(bz)
  })
  # restruct data in dataframe
  syn_data <- do.call(rbind, l_curves)
  colnames(syn_data) <- c('x', 'y')
  syn_data$id <- rep(1:n, each=steps)
  
  # ----- run hierarchical clustering -----
  
  empty_list <- vector('list', length=length(kseq))
  l_real <- list('clusters'=empty_list, 'WCD'=empty_list)
  l_syn <- list('clusters'=empty_list, 'WCD'=empty_list, 'SynData'=syn_data)
  l_all <- list(l_real, l_syn)
  
  for(type in 1:2) { # 1=real, 2=synthetic
    
    if(type==1) { # data has to come first, because we then overwrite 'dist'
      distobj <- dist(dist)
    } else {
      dist <- mta:::distmat(1:n, syn_data$x , syn_data$y, steps)
      distobj <- dist(dist)
    }
    
    hc <- hclust(distobj, method = 'complete')
    
    for(k in kseq) {
      
      l_all[[type]][[1]][[k]] <- cl <- cutree(hc, k) # cut tree & define cluster
      # compute mean within cluster dissimilarity
      l_diss <- list()
      for(i in 1:k) l_diss[[i]] <- mean(dist[cl==i, cl==i])
      l_all[[type]][[2]][[k]] <- sum(unlist(l_diss)*table(cl)/sum(table(cl))) #weighting of dissimilarities
      
    }
    
  } # end of loop type
  
  outlist <- list('WCD_data'=unlist(l_all[[1]][[2]]), 
                  'WCD_syn'=unlist(l_all[[2]][[2]]),
                  'Gaps'= unlist(l_all[[2]][[2]])-unlist(l_all[[1]][[2]]),
                  'cl_data'=l_all[[1]],
                  'cl_syn'=l_all[[2]])
  
  return(outlist)
  
} # EoF
