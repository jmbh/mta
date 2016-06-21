
library(MASS)
library(cluster)

s <- .3
n <- 300
kseq <- 1:20
l_stability <- l_gap <- l_jump <- l_slope <- list()

l_res <- list()

for(method in 1:2) {
  
for(j in 1:10) {
  
  # ---------- Generate Data ----------
  set.seed(j)
  
  data <- rbind(cbind(rnorm(n, 1, s),rnorm(n, 1, s)),
                cbind(rnorm(n, 1, s),rnorm(n, 2.5, s)),
                cbind(rnorm(n, 1, s),rnorm(n, 3, s)),
                cbind(rnorm(n, 2, s),rnorm(n, 1, s)),
                cbind(rnorm(n, 2, s),rnorm(n, 1.5, s)),
                cbind(rnorm(n, 2, s),rnorm(n, 3, s)),
                cbind(rnorm(n, 3, s),rnorm(n, 1, s)),
                cbind(rnorm(n, 3, s),rnorm(n, 1.5, s)),
                cbind(rnorm(n, 3, s),rnorm(n, 3, s)))
  
  plot(data) # show data
  #dim(data)
  dist <- distm <- as.matrix(dist(data)) # calc distance matrix
  
  # ---------- Call different Functions ----------
  
  kseq <- 2:20
  C2 <- cluster_stability2(dist, kseq, Bcomp = 30, norm=TRUE)
  C2_n <- cluster_stability2(dist, kseq, Bcomp = 30, norm=FALSE)
 
  plot(kseq, C2$instabilities, ylim=c(0,.8), type='l')
  lines(kseq, C2_n$instabilities, col='red')
  abline(v=9, col='blue')
  
  
  
  
  cluster_stability2
  # speed testing
  tt<-proc.time()[3]
  set.seed(2)
  C2 <- cluster_stability2(dist=distm, kseq = kseq, Bcomp = 10, norm = TRUE, pbar=FALSE)
  C2
  C2$instabilities
  plot(C2$instabilities)
  proc.time()[3]-tt
  
  # a) Cluster stability
  cl_a <- cluster_stability(dist=distm, kseq = kseq, B = 20, norm = FALSE)
  l_stability[[j]] <- cl_a$kopt
  
  # b,c,d) Gap, Jump, Slope Statistic
  cl_b <- gap_statistic(dist = distm, kseq = kseq, method = 'gauss', dim = 2)
  l_gap[[j]] <- cl_b$kopt
  l_jump[[j]] <- cl_b$kopt_jump
  l_slope[[j]] <- cl_b$kopt_sil

  print(j)  
}

stability <- unlist(l_stability)
gap <- unlist(l_gap)
jump <- unlist(l_jump)
slope <- unlist(l_slope)

l_res[[method]] <- cbind(stability, gap, jump, slope)

}

par(mfrow=c(2,1), mar=c(2,2,2,2))
boxplot(l_res[[1]], ylim=c(1,10))
boxplot(l_res[[2]], ylim=c(1,10))



