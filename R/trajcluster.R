


trajcluster <- function(data, 
                        i.xyt, 
                        i.id, 
                        side, 
                        nclust = 3, #numer of desired clusters if no prototypes
                        prototypes = NA,  # list of prototypical trajectories
                        ... # some other things that we pass to the cluster function
)
  
{
  
  # output: good question! but definitely:
    # 1) all input (except data)
    # 2) if unsupervised: the fast-clust object 
    # 3) if supervised (protos): distance matrix 
    # 4) flagging of trajectories
    # ???
  
  
}