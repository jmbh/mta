

# dummy input
#x <- dat
#Bcomp <- 10
#kseq <- 2:10
#norm <- TRUE

# testing function

#p0 <- cluster_stability2(x, 2:10, norm=FALSE, prediction=FALSE)
#p1 <- cluster_stability2(x, 2:10, norm=FALSE, prediction=TRUE)

#plot(p0$instabilities, type='l', col='red', ylim=c(0,.4))
#lines(p1$instabilities)

# runs with k-means only

cluster_stability2 <- function(x, # n x p data matrix 
                               kseq, # sequence of ks tested
                               Bcomp = 10, # number of bootstrap comparisons
                               norm = FALSE, # norm over pw equal assign,FALSE=as in Wang etal
                               prediction = TRUE, # use prediction approach, if FALSE, use brute pair in equal cluster approach
                               pbar = TRUE,
                               ...) # other arguments passed to hclust
{
  
  # Input checks
  p <- nrow(x)
  m_instab <- matrix(NA, Bcomp, length(kseq)) 
  
  # Compute necesarry bootstrap samples to get Bcomp comparisons
  Bsamp <- ceiling(.5 * (sqrt(8*Bcomp+1) + 1)) # solution to equation Bcomp = Bsamp(Bsamp-1)/2
  
  # ----- Draw bootstrap samples -----
  
  l_ind <- list()
  for(b in 1:Bsamp) {
    l_ind[[b]] <-  sample(1:p, p, replace=T)
    l_ind[[b]] <- l_ind[[b]][order(l_ind[[b]])] # order
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
      
      if(prediction) {
        
        l_km_models <- list()
        count <- 1
        for(r in combs[,bc]) {
          l_km_models[[count]] <- kcca(x[l_ind[[r]],], k=k, kccaFamily("kmeans")) #save whole mode
          count <- count+1
        }
        
        # make predictions
        cl_a <- predict <- predict(l_km_models[[1]], newdata=x)
        cl_b <- predict <- predict(l_km_models[[2]], newdata=x)
        
        # count pairwise equal assignments
        same_a <- as.numeric(dist(cl_a)==0)*1
        same_b <- as.numeric(dist(cl_b)==0)*1
        
        # compute Instability
        InStab <- mean(abs(same_a - same_b)) 
        
        # Normalize
        if(norm==FALSE) {
          m_instab[bc,which(kseq==k)] <- InStab
        } else {
          tb1 <- table(l_cl[[1]])
          tb2 <- table(l_cl[[2]])
          norm_val <- instab(tb1, tb2, 100)
          m_instab[bc,which(kseq==k)] <- InStab/norm_val
        }
        
        
      } else {
        
        l_cl <- list()
        l_pairind <- list() # are two objects in same cluster (1 yes, 0 no)
        count <- 1
        for(r in combs[,bc]) {
          cl_long <- kcca(x[l_ind[[r]],], k=k, kccaFamily("kmeans"))@second
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
        
      } # end if: prediction TRUE/FALSE
      
    } # end for k
    
    if(pbar) setTxtProgressBar(pb, bc)
    
  } # end for B
  
  
  
  instabM <- colMeans(m_instab)
  kopt <- which.min(instabM)+(min(kseq)-1)
  
  
  outlist <- list('instabilities'=instabM, 'kopt'=kopt)
  
  return(outlist)
  
} # EoF

