
gap_statistic2 <- function(x, # n x p data matrix
                           kseq, #sequence of ks to be checked 
                           kmIter) # restarts of k means algorithm
{ 
  
  # ----- aux functionS -----

  getWCFandSil <- function(x, k, type='data') {
    
    # clustering
    if(k==1) {
      cl <- rep(1,nrow(x))
    } else {
      check_k <- FALSE
      counter <- 0
      while(check_k == FALSE) {
        
        # re-start k means algorithm a couple of times
        l_km <- list()
        for(km in 1:kmIter) {
          l_km[[km]] <- flexclust::kcca(x, k=k, kccaFamily("kmeans")) #save whole model  
        }
        
        WCD <- unlist(lapply(l_km, function(x) mean(x@clusinfo$av_dist)))
        km_model <- l_km[[which.min(WCD)]] # pick k-means clustering with smallest WCD
        
        cl <- km_model@cluster
        
        if(length(unique(cl))==k) check_k <- TRUE
        counter <- counter + 1
        if(counter>100) stop(paste0('k means solution with ', k, 'centers always converges to solution with empty cluster.'))
      }

    }
    
    # calc distance matric
    dmat <- as.matrix(dist(x))
    
    # calc within cluster dissimilarity
    l_diss <- list()
    for(i in 1:k) l_diss[[i]] <- mean(dmat[cl==i, cl==i])
    
    tb <- table(cl)
    WCD <- sum(unlist(l_diss) * tb/sum(tb))
    
    if(k>1) {
      if(type=='data') {
        Sil <- mean(silhouette(cl, dmat)[,3])
      } else {
        Sil <- NULL
      }
    } else {
      Sil <- NULL
    }
    
    # calc MSE (for jump statistic)
    v_mse <- numeric(nrow(x))
    if(k>1) {
    for(i in 1:nrow(x)) {
      diffs <- (x[i,] - km_model@centers[km_model@cluster[i],])
      mse <- t(diffs) %*% diffs
      v_mse[i] <- mse
    }
    MSE <- sum(v_mse) / ncol(x)
    } else {
      SE <- apply(x, 1, function(inst){
        diffs <- inst - colMeans(x)
        return(t(diffs) %*% diffs)
      })
      MSE <- sum(SE) / ncol(x)
    }
    
    outlist <- list('WCD' = WCD, 'Sil'=Sil, 'MSE'=MSE)
    return(outlist)
  }
  
  
  UniformData <- function(x) {
    
    # get dimensions of data
    unifdims <- t(apply(x, 2, range))
    diffs <- abs(unifdims[,1] - unifdims[,2])
    unifdims2 <- cbind(0,diffs)  
    
    # sample data and combine to data matrix
    unif_dims <- list()
    for(i in 1:dims) unif_dims[[i]] <- runif(n, unifdims2[i,1], unifdims2[i,2])
    data_synt <- do.call(cbind, unif_dims)
    
    return(data_synt)
  }
  
  
  # ----- run k-means clustering on real data and generated data -----
  
  n <- nrow(x)
  dims <- ncol(x)
  
  # First: Real Data
  l_WCD_data <- l_Sil <- l_MSE <-  list()
  for(k in c(1, kseq)) {
    obj <- getWCFandSil(x, k, type = 'data')
    l_WCD_data[[k]] <- obj$WCD
    l_Sil[[k]] <- obj$Sil
    l_MSE[[k]] <- obj$MSE
  }
  l_WCD_data <- unlist(l_WCD_data)
  l_Sil <- unlist(l_Sil)
  l_MSE <- unlist(l_MSE)
  
  # Then: Synthetic Data
  l_WCD_syndata_runs <- list()
  for(runs in 1:10) {
    l_WCD_syndata1 <- list()
    data_synt <- UniformData(x)
    for(k in c(1,kseq)) {
      obj <- getWCFandSil(data_synt, k, type = 'simulated')
      l_WCD_syndata1[[k]] <- obj$WCD
    }
    l_WCD_syndata_runs[[runs]] <- unlist(l_WCD_syndata1)
  }
  l_WCD_syntetic <- colMeans(do.call(rbind, l_WCD_syndata_runs))
  
  # Compute Gap Statistic
  WCD_data_log <- log(l_WCD_data)
  WCD_syn_log <- log(l_WCD_syntetic)
  
  WCD_data_log <- WCD_data_log - WCD_data_log[1]
  WCD_syn_log <- WCD_syn_log - WCD_syn_log[1]
  
  Gaps <- WCD_syn_log-WCD_data_log
  
  # Compute Slope Statistic
  sil <- l_Sil
  p <- 1
  slope <- -(sil[-1] - sil[-length(sil)]) * sil[-1]^p
  kopt_sil <- which.max(slope)+1 #add one because we have sequence 2:...
  ## Also Compute JUMP Statistic
  MSE_transf <- l_MSE^(- dims/2)
  jump <- (MSE_transf - c(0, MSE_transf[-length(MSE_transf)]))[-1]
  k_jump <- which.max(jump) + 1 #add one because we have sequence 2:...
  
  outlist <- list('WCD_data'=l_WCD_data, 
                  'WCD_syn'=l_WCD_syntetic,
                  'Gaps'= Gaps,
                  'kopt'=which.max(Gaps),
                  'kopt_jump'=k_jump,
                  'silhouettes'=sil,
                  'jumps'=jump,
                  'kopt_sil'=kopt_sil) 
  
  return(outlist)
  
} # EoF


