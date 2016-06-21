
gap_statistic <- function(dist, # p x p distance matrix of
                          kseq, #sequence of ks to be checked 
                          steps=10, # number of points on synthetic curves
                          method='gauss', # 'bezier' for bezier curves, 'gauss' for a multivariate gaussian with cov=0, matching the dimension
                          dim = NULL,
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
  
  if(method=='bezier') {
    
    # create curves
    l_curves <- vector('list', length=n)
    set.seed(1)
    l_curves <- lapply(l_curves, function(s) {
      bz <-  rBezier(xcor,ycor,w=rexp(1,lambda), resol=steps)
      bz <- matrix(c(bz$x, bz$y), nrow=1)
      return(bz)
    })
    # restruct data in dataframe
    syn_data <- do.call(rbind, l_curves)
    #colnames(syn_data) <- c('x', 'y')
    #syn_data$id <- rep(1:n, each=steps)
    
  } else {
    
    Sigma <- matrix(0, dim, dim)
    diag(Sigma) <- 1
    mu <- rep(0, dim)
    syn_data <- mvrnorm(n, mu, Sigma)
    
  }
  
  # ----- run hierarchical clustering -----
  
  empty_list <- vector('list', length=length(kseq))
  l_real <- list('clusters'=empty_list, 'WCD'=empty_list)
  l_syn <- list('clusters'=empty_list, 'WCD'=empty_list, 'SynData'=syn_data)
  l_all <- list(l_real, l_syn)
  l_silhouette <- list()
  
  for(type in 1:2) { # 1=real, 2=synthetic
    
    if(type==1) { # data has to come first, because we then overwrite 'dist'
      distobj <- as.dist(dist)
    } else {
      distobj <- dist(syn_data)
      dist <- as.matrix(distobj)
    }
    
    hc <- hclust(distobj, method = 'complete')
    
    for(k in kseq) {
      
      l_all[[type]][[1]][[k]] <- cl <- cutree(hc, k) # cut tree & define cluster
      # compute mean within cluster dissimilarity
      l_diss <- list()
      
      for(i in 1:k) l_diss[[i]] <- mean(dist[cl==i, cl==i])
      
      l_all[[type]][[2]][[k]] <- sum(unlist(l_diss)*table(cl)/sum(table(cl))) #weighting of dissimilarities
      
      # compute silhouette statistic
      if(type==1){
        if(k==1) { 
          l_silhouette[[k]] <- NA
        } else {
          l_silhouette[[k]] <- mean(silhouette(cl, dist)[,3])
          }
        
      }  
    }
    
  } # end of loop type
  
  # Prepare output Gap Statistic
  WCD_data = unlist(l_all[[1]][[2]]) # Within cluster D for Data
  WCD_syn = unlist(l_all[[2]][[2]])
  
  # take log
  WCD_data_log <- log(WCD_data)
  WCD_syn_log <- log(WCD_syn)
  
  WCD_data_log <- WCD_data_log - WCD_data_log[1]
  WCD_syn_log <- WCD_syn_log - WCD_syn_log[1]
  
  Gaps <- WCD_syn-WCD_data
  
  # Also Compute JUMP Statistic
  WCD_transf <- WCD_data^(- dim/2)
  jump <- (WCD_transf - c(0,WCD_transf[-length(WCD_transf)]))[-1]
  k_jump <- which.max(jump) + 1
  
  # Compute Slope Statistic
  sil <- unlist(l_silhouette)
  p <- 1
  slope <- -(sil[-1] - sil[-length(sil)]) * sil[-1]^p
  
  
  outlist <- list('WCD_data'=WCD_data, 
                  'WCD_syn'=unlist(l_all[[2]][[2]]),
                  'Gaps'= Gaps,
                  'kopt'=which.max(Gaps),
                  'kopt_jump'=k_jump,
                  'silhouettes'=sil,
                  'kopt_sil'=which.max(slope),
                  'cl_data'=l_all[[1]],
                  'cl_syn'=l_all[[2]])
  
  return(outlist)
  
} # EoF
