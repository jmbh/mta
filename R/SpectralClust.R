#from: http://stackoverflow.com/questions/31074113/cluster-unseen-points-using-spectral-clustering



S <- SpectralClust(train, k=k, test=test) #make predictions


SpectralClust <- function(train, k, test) {
  
  ##The paper gives great instructions on how to perform spectral clustering
  ##so I'll be following the steps
  ##Ng, A. Y., Jordan, M. I., & Weiss, Y. (2002). On spectral clustering: Analysis and an algorithm. Advances in neural information processing systems, 2, 849-856.
  ##Pg.2 http://papers.nips.cc/paper/2092-on-spectral-clustering-analysis-and-an-algorithm.pdf
  #1. Form the affinity matrix
  #k = 2L #This is the number ofo clusters we will train
  K = rbfdot(sigma = 300) #Our kernel
  A = kernelMatrix(K, train) #Caution choosing your kernel product function, some have higher numerical imprecision
  diag(A) = 0
  #2. Define the diagonal matrix D and the laplacean matrix L
  D = diag(rowSums(A))
  L = diag(1/sqrt(diag(D))) %*% A %*% diag(1/sqrt(diag(D)))
  #3. Find the eigenvectors of L
  X = eigen(L, symmetric = TRUE)$vectors[,1:k]
  #4. Form Y from X
  Y = X/sqrt(rowSums(X^2))
  #5. Cluster (k-means)
  kM = kmeans(Y, centers = k, iter.max = 100L, nstart = 1000L)
  #6. This is the cluster assignment of the original data
  cl = fitted(kM, "classes")
  ##Projection on eigen vectors, see the ranges, it shows how there's a single preferential direction
  #plot(jitter(Y, .1), ylab = "2nd eigenfunction", xlab = "1st eigenfunction", col = adjustcolor(rainbow(3)[2*cl-1], alpha = .5))
  
  if(is.null(test)) {
    
    return(cl)
    
  } else {
    
    ##LET'S TRY TEST DATA NOW
    B = kernelMatrix(K, test, train) #The kernel product between train and test data
    
    ##We project on the learned eigenfunctions
    f = tcrossprod(B, t(Y))
    #This part is described in Bengio, Y., Vincent, P., Paiement, J. F., Delalleau, O., Ouimet, M., & Le Roux, N. (2003). Spectral clustering and kernel PCA are learning eigenfunctions (Vol. 1239). CIRANO.
    #Pg.12 http://www.cirano.qc.ca/pdf/publication/2003s-19.pdf
    
    ##And assign clusters based on the centers in that space
    new.cl = apply(as.matrix(f), 1, function(x) { which.max(tcrossprod(x,kM$centers)) } ) #This computes the distance to the k-means centers on the transformed space
    
    return(new.cl)
    
  }
  
}
