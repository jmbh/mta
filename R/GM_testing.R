
library(MASS)
library(cluster)

s <- 1
n <- 150
kseq <- 1:20
plot(data)
l_stability <- l_gap <- l_jump <- l_slope <- list()

l_res <- list()

for(method in 1:2) {
  
for(j in 1:10) {
  
  # ---------- Generate Data ----------
  set.seed(j)
  data <- rbind(cbind(rnorm(n, 1, s),rnorm(n, 0, s)),
                cbind(rnorm(n, 1, s),rnorm(n, 1, s)),
                cbind(rnorm(n, 0, s),rnorm(n, 1, s)),
                cbind(rnorm(n, 0, s),rnorm(n, 0, s)))
  
  #plot(data) # show data
  #dim(data)
  dist <- distm <- as.matrix(dist(data)) # calc distance matrix
  
  # ---------- Call different Functions ----------
  
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



