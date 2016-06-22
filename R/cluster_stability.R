
cluster_stability <- function(dist, # p x p distance matrix 
                              kseq, # sequence of ks tested
                              Bcomp = 10, # number of bootstrap comparisons
                              norm = FALSE, # norm over pw equal assign,FALSE=as in Wang etal
                              pbar = TRUE,
                              linkage='complete',
                              ...) # other arguments passed to hclust
{
  
  # Input checks
  p <- ncol(dist)
  m_instab <- matrix(NA, Bcomp, length(kseq)) 
  
  # Compute necesarry bootstrap samples to get Bcomp comparisons
  Bsamp <- ceiling(.5 * (sqrt(8*Bcomp+1) + 1)) # solution to equation Bcomp = Bsamp(Bsamp-1)/2
  
  # ----- Draw bootstrap samples -----
  
  l_ind <- list()
  l_dist <- list()
  l_clust <- list()
  for(b in 1:Bsamp) {
    l_ind[[b]] <-  sample(1:p, p, replace=T)
    l_ind[[b]] <- l_ind[[b]][order(l_ind[[b]])] # order
    l_dist[[b]] <- dist[ l_ind[[b]],  l_ind[[b]] ]
    l_clust[[b]] <- fastcluster::hclust(as.dist(l_dist[[b]]), method=linkage)
  }
  
  
  # ----- compute indices for each pair: which objects are in both samples? -----
  
  combs <- combn(1:Bsamp ,2) # All possible combinations
  l_indices <- list()
  for(bc in 1:Bcomp) {
    # overlap
    IS <- intersect(l_ind[[combs[1,bc]]],l_ind[[combs[2,bc]]])
    l_pair <- list()
    count <- 1
    for(r in combs[,bc]) {
      
      ind_is_r1 <- l_ind[[r]] %in% IS # indicator: indice in both?
      ind2_r1 <- !duplicated(l_ind[[r]])
      ind3_r1 <- ind_is_r1 & ind2_r1
      
      #ind_is_r2 <- l_ind[[2]] %in% IS # indicator: indice in both?
      #ind2_r2 <- !duplicated(l_ind[[2]])
      #ind3_r2 <- ind_is_r2 & ind2_r2
      
      l_pair[[count]] <- ind3_r1
      count<-count+1
    }
    l_indices[[bc]] <- l_pair
  }
  
  #l_ind[[1]][l_indices[[1]][[1]]]
  #l_ind[[2]][l_indices[[1]][[2]]]
  
  
  # ----- Loop over first Bcomp comparisons -----
  
  if(pbar)  pb <- txtProgressBar(min=0, max=Bcomp, style = 2)
  
  for(bc in 1:Bcomp) {
    
    for(k in kseq) {

      l_cl <- list()
      l_pairind <- list() # are two objects in same cluster (1 yes, 0 no)
      count <- 1
      for(r in combs[,bc]) {
        cl_long  <- cutree(l_clust[[r]], k) # cut dendogram
        l_cl[[count]] <- cl_long[l_indices[[bc]][[count]]] # only take the ones in the intersection set
        l_pairind[[count]] <- (as.numeric(dist(l_cl[[count]]))==0)*1
        count <- count+1
      }
    
      InStab <- mean(abs(l_pairind[[1]] - l_pairind[[2]])) 
      
      # Normalize
      if(norm==FALSE) {
        m_instab[bc,which(kseq==k)] <- InStab
      } else {
        tb1 <- table(l_cl[[1]])
        tb2 <- table(l_cl[[2]])
        norm_val <- instab(tb1, tb2, 100)
        m_instab[bc,which(kseq==k)] <- InStab/norm_val
      }
        
    } # end for k
  
    if(pbar) setTxtProgressBar(pb, bc)
  
  } # end for B
  
  
  
  instabM <- colMeans(m_instab)
  kopt <- which.min(instabM)+(min(kseq)-1)

  
  outlist <- list('instabilities'=instabM, 'kopt'=kopt)
  
  return(outlist)
  
} # EoF
