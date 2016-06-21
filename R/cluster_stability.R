

cluster_stability <- function(dist, # p x p distance matrix 
                              kseq, # sequence of ks tested
                              Bcomp = 10, # number of bootstrap comparisons
                              norm = FALSE, # norm over pw equal assign,FALSE=as in Wang etal
                              pbar = TRUE,
                              ...) # other arguments passed to hclust
{
  
  # input checks
  p <- ncol(dist)
  m_instab <- matrix(NA, Bcomp, length(kseq)) 
  
  # compute necesarry bootstrap samples to get Bcomp comparisons
  Bsamp <- ceiling(.5 * (sqrt(8*Bcomp+1) + 1)) # solution to equation Bcomp = Bsamp(Bsamp-1)/2
  
  # draw bootstrap samples
  l_ind <- list()
  l_dist <- list()
  l_clust <- list()
  for(b in 1:Bsamp) {
    l_ind[[b]] <-  sample(1:p, p, replace=T)
    l_dist[[b]] <- dist[ l_ind[[b]],  l_ind[[b]]]
    l_clust[[b]] <- fastcluster::hclust(as.dist(l_dist[[b]]), method='complete')
  }
  
  # all possible combinations
  combs <- combn(1:Bsamp ,2)
  
  if(pbar)  pb <- txtProgressBar(min=0, max=Bcomp, style = 2)

  # loop over first Bcomp comparisons
  for(bc in 1:Bcomp) {
    
    # for k (most efficient nesting)
    for(k in kseq) {
      
      l_cl <- list()
      l_pairwise <- list() # save pairwise assignment similarities

      for(samp in combs[,bc]) {
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
      no_both <- (1:p)[(1:p %in% l_ind[[combs[1,bc]]]) & (1:p %in% l_ind[[combs[2,bc]]])]
      
      l_pairwise_short <- list()
      for(samp in combs[,bc]) {
        # subset
        ind <- colnames(l_pairwise[[samp]]) %in% no_both
        l_pairwise_short[[samp]] <- l_pairwise[[samp]][ind, ind]
        # remove duplicates
        dupl <- !duplicated(colnames(l_pairwise_short[[samp]]))
        l_pairwise_short[[samp]] <- l_pairwise_short[[samp]][dupl, dupl]
      }
      
      # compute instability
      if(norm) {
        #norm_val <- max(1,(sum(l_pairwise_short[[1]]) + sum(l_pairwise_short[[1]]) ) / 2)
        tab1 <- as.numeric(table(l_cl[[combs[1,bc]]]))
        tab2 <- as.numeric(table(l_cl[[combs[2,bc]]]))
        
        norm_val <- instab(tab1, tab2, 100)
        
        ind_mat1 <- l_pairwise_short[[combs[1,bc]]]
        ind_mat2 <- l_pairwise_short[[combs[2,bc]]]
        
        scale_dim <- (dim(ind_mat1)[1] * (dim(ind_mat1)[1]-1)) / 2
        
        m_instab[bc,which(kseq==k)] <- sum(abs(ind_mat1 - ind_mat2)) / (norm_val*scale_dim) 
        } else {
        m_instab[bc,which(kseq==k)] <- mean(abs(l_pairwise_short[[combs[1,bc]]] - l_pairwise_short[[combs[1,bc]]])) 
      }
        
    } # end for k
  
    if(pbar) setTxtProgressBar(pb, bc)
  
    
  } # end for B
  
  stab <- colMeans(m_instab)
  kopt <- which.min(stab[-1])+1

  
  outlist <- list('stabilities'=stab, 'kopt'=kopt)
  
  return(outlist)
  
} # EoF
