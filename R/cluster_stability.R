
cluster_stability <- function(dist, # p x p distance matrix 
                              kseq, # sequence of ks tested
                              B = 10, # number of bootstrap samples 
                              norm = FALSE, # norm over pw equal assign,FALSE=as in Wang etal
                              ...) # other arguments passed to hclust
{
  # input checks
  p <- ncol(dist)
  m_instab <- matrix(NA, B, length(kseq)) 
  
  # for bootstrap
  for(b in 1:B) {
    
    # take two bootstrap samples
    l_dist <- list()
    l_ind <- list()
    l_clust <- list()
    for(samp in 1:2) {
      l_ind[[samp]] <- sample(1:p, p, replace=T)
      l_dist[[samp]] <- dist[ l_ind[[samp]],  l_ind[[samp]]]
      l_clust[[samp]] <- fastcluster::hclust(as.dist(l_dist[[samp]]), method='complete')
    }
    
    # for k (most efficient nesting)
    for(k in kseq) {
      
      l_cl <- list()
      l_pairwise <- list() # save pairwise assignment similarities
      
      for(samp in 1:2) {
        # cut dendogram
        l_cl[[samp]]  <- cutree(l_clust[[samp]], k)
        
        # get indicator matrix, entry 1 = same cluster
        dm <- (as.matrix(dist(l_cl[[samp]]))==0)*1 
        diag(dm) <- 0
        colnames(dm) <- rownames(dm) <- l_ind[[samp]]
        dm <- dm[order(l_ind[[samp]]), order(l_ind[[samp]])]
        l_pairwise[[samp]] <- dm
      }
      
      # take intersection of objects in bi,bj
      no_both <- (1:p)[(1:p %in% l_ind[[1]]) & (1:p %in% l_ind[[2]])]
      
      l_pairwise_short <- list()
      for(samp in 1:2) {
        # subset
        ind <- colnames(l_pairwise[[samp]]) %in% no_both
        l_pairwise_short[[samp]] <- l_pairwise[[samp]][ind, ind]
        # remove duplicates
        dupl <- !duplicated(colnames(l_pairwise_short[[samp]]))
        l_pairwise_short[[samp]] <- l_pairwise_short[[samp]][dupl, dupl]
      }
      
      # compute instability
      if(norm) {
        norm_val <- max(1,(sum(l_pairwise_short[[1]]) + sum(l_pairwise_short[[1]]) ) / 2)
        m_instab[b,which(kseq==k)] <- sum(abs(l_pairwise_short[[1]] - l_pairwise_short[[2]])) / norm_val 
      } else {
        m_instab[b,which(kseq==k)] <- mean(abs(l_pairwise_short[[1]] - l_pairwise_short[[2]])) 
      }
      
    } # end for k
    
  } # end for B
  
  return(colMeans(m_instab))
  
} # EoF
