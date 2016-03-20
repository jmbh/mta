

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
                        subsampN = NA # maximal N number of trajectories; to render analysis possible in the case of huge number of trajectories
)
  
{
  
  # +++ check input +++
  if(sum(is.na(data))>0) {stop('No missing values permitted.')}
  rowcheck <- ddply(data, i.id, nrow)
  steps <- rowcheck[,ncol(rowcheck)]
  if(var(steps)>0) {stop('All trajectories must have the same length!')}
  
  # +++ compute some derivatives form data +++
  n_trial <- nrow(rowcheck)
  side <- getside(data, i.id)
  step <- steps[1]
  data$x[side==1] <- data$x[side==1]*-1 #flip all trjectories to the left
  
  # +++ subsampling +++
  if(is.na(subsampN)==FALSE) {
    if(n_trials>subsampN) {
      id <- rep(1:n_trials, each=step)
      s_id <- sample(1:n_trials, subsampN, replace = FALSE)
      data <- dat[id %in% s_id,]  
    }
  }
  
  # +++ rescale trajectories +++
  dat_resc <- spatialRescale(data, i.id, i.xyt, nResc)
  
  # ++++++++++++++++++ hierarchical clustering ++++++++++++++++++
  
  if('hierarchical' %in% type) {
  
  # +++ calculate distance matrix +++
  id <- 1:n_trial
  distm <- distmat(id, dat_resc$x , dat_resc$y, nResc)

  # +++ clustering +++
  md <- as.dist(distm)
  clust_obj <- fastcluster::hclust(md, method = "ward.D")
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
      l_proto_sc[[i]] <- spatialRescale(pr_df, "id", i.xyt, nResc)   
    }
    
    # +++ calculate distances & assign lable +++
    proto_dists <- ddply(dat_resc, i.id, function(x) {
      dists <- rep(NA,n_proto)
      for(i in 1:n_proto) {
        dat_pr_i <- l_proto_sc[[i]]
        dists[i] <- sqrt(sum(c((x[,i.xyt[1]] - dat_pr_i[,2])^2,
                               (x[,i.xyt[2]] - dat_pr_i[,3])^2))) 
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

