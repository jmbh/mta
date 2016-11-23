

## TESTING 

#library(plyr)
#library(Rcpp)
#library(fastcluster)

# get functions
#setwd("G:\\_projects_ongoing\\MPI\\mta\\R")
#source("getside.R")
#source("spatialRescale.R")
#source("prepr.R")
#setwd("G:\\_projects_ongoing\\MPI\\mta\\src")
#sourceCpp("distmat.cpp")

# data
#setwd("G:\\_projects_ongoing\\MPI\\mta\\data")
#load("prototypes.RData")
#load("data_sp2015.RData")

# function input for testing
#i.xyt <- c('x', 'y', 't')
#i.id <- c('id.ptp', 'id.trial')
#stretch <- list("start"=c(0,0), "left"=c(-1,1.5), "right"=c(1,1.5))

#dat_norm <- prepr(data_sp2015, i.xyt, i.id, stretch = stretch)
#data <- dat_norm$data

#clus_obj <- trajcluster(data, i.xyt, i.id, type=c("hierarchical", "prototypes"),
#                        nclust=4, nResc = 10, prototypes = prototypes) 



trajcluster <- function(data, 
                        i.xyt, 
                        i.id, 
                        type=c("hierarchical", "prototypes"),
                        nclust = 4, # numer of desired clusters if no prototypes
                        nResc = 10, # number of data points for spatial rescaling
                        prototypes = NA, # list of prototypical trajectories
                        subsampN = NA, # maximal N number of trajectories; to render analysis possible in the case of huge number of trajectories
                        stretch = list('start' = c(0,0),'left' = c(-1,1.5)), #
                        method = 'complete')
  
{
  
  
  # +++ subsampling +++
  n_trials = nrow(unique(data[,i.id]))
  if(is.na(subsampN)==FALSE) {
    if(n_trials>subsampN) {
      ids = apply(data[,i.id],1,paste0,collapse='-')
      sel = sample(unique(ids),subsampN)
      data <- data[ids %in% sel,]  
    }
  }

  # +++ stretch +++
  dat_str <- ddply(data, i.id, function(traj) {
      # starting point
      X <- traj$x - traj$x[1]; X <- X + stretch$start[1]
      Y <- traj$y - traj$y[1]; Y <- Y + stretch$start[2]
      # end point
      X <- (X / abs(X[length(X)])) * abs(stretch$left[1])     
      Y <- (Y / abs(Y[length(Y)])) * abs(stretch$left[2])     
      t <- traj$t
      res = data.frame(cbind(X,Y,t))
      names(res) = i.xyt
      return(res)
      })

  # +++ compute some derivatives form data +++
  side <- getside(dat_str, i.id)
  dat_str[,i.xyt[1]][side==1] <- dat_str[,i.xyt[1]][side==1]*-1 #flip all trjectories to the left
    
  # +++ check input +++
  if(sum(is.na(dat_str))>0) {stop('No missing values permitted.')}
  
  # +++ rescale trajectories +++
  dat_resc <- spatialRescale(dat_str, i.id, i.xyt, nResc)
  
  
  # ++++++++++++++++++ hierarchical clustering ++++++++++++++++++
  
  if('hierarchical' %in% type) {
  
  # +++ calculate distance matrix +++
  id <- 1:nrow(unique(dat_resc[,i.id]))
  distm <- distmat(id, dat_resc$x , dat_resc$y, nResc)

  # +++ clustering +++
  md <- as.dist(distm)
  clust_obj <- fastcluster::hclust(md, method = method, )
  cl <- rep(cutree(clust_obj, nclust), each=nResc)
  dat_resc$cl <- cl
  
  # +++ clustering measures: within/between variance, #cluster suggestion (gap stat),  +++
  
  # clusterstat() function/ TO BE DONE
  
  # +++ define output +++
  
  l_hier <- list('distmat'=distm, 'clust_obj'=clust_obj, "cluster"=cl)
  
  } else {  l_hier =NULL }
  
  # ++++++++++++++++++ prototype clustering ++++++++++++++++++
  
  if('prototypes' %in% type) {
    
    # +++ rescale prototypes +++
    l_proto_sc <- list()
    n_proto <- length(prototypes)
    for(i in 1:n_proto){
      pr_df <- as.data.frame(prototypes[[i]])
      colnames(pr_df) <- c("x", "y")
      pr_df$id <- i
      l_proto_sc[[i]] <- f_rescale_c(pr_df[,1], pr_df[,2], nResc)   
    }
    
    # +++ calculate distances & assign lable +++
    proto_dists <- ddply(dat_resc, i.id, function(x) {
      dists <- rep(NA,n_proto)
      for(i in 1:n_proto) {
        dat_pr_i <- l_proto_sc[[i]]
        dists[i] <- sqrt(sum(c((x[,i.xyt[1]] - dat_pr_i[,1])^2,
                               (x[,i.xyt[2]] - dat_pr_i[,2])^2))) 
      }
      return(dists)
    })
    colnames(proto_dists)[-(1:length(i.id))] <- paste0("DistProto",1:n_proto)
    proto_dists$cproto <- apply(proto_dists[,-(1:length(i.id))], 1, which.min)
    
    # +++ define output +++
    l_proto <- list('prototypes_resc'=l_proto_sc,
                    'protoclust'=proto_dists)    
    
  } else {  l_proto = NULL }

  # +++ output +++
  
  l_call <- list('i.xyt'=i.xyt, 'i.id'=i.id, 'nclust'=nclust,  'prototypes'=prototypes, 
                 'side'=side, 'subsampN'=subsampN, 'nResc'=nResc)
  
  outlist <- list('call'=l_call,
                  'data_res'=dat_resc,
                  'hierarchical'=l_hier,
                  'prototypes'=l_proto) 

  return(outlist)
  
}

