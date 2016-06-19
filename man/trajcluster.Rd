\name{trajcluster}

\alias{trajcluster}

\title{
Clustering mouse tracking trajectories
}
\description{
Function to cluster mouse trajectory using hierarchical cluster analysis or distance to given prototypical trajectories}
\usage{
trajcluster(data, i.xyt, i.id, type=c("hierarchical", "prototypes"), nclust, nResc = 10, 
      prototypes = NA, subsampN = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A dataframe with x,y values of time-normalizes trajectories (all trajectories have the same length) and id variable(s).}
  \item{i.xyt}{
Vector containing the column names of the x, y, t variables in that order.
}
  \item{i.id}{
Vector containing the column names of the indicator variables that uniquely identify single trajectories (e.g. c('experiment1', 'trial')).
}
  \item{type}{\code{type='hierarchical'} performs hierarchical cluster analysis on the trajectories. The distance measure is the sum of the euclidean distances between points 1,2,...,n of the two compared trajectories. The number of extracted clusters is specified with the parameter \code{nclust}.
  \code{type='hierarchical'} calculates the distance of each trajectory to a number of prototypical trajectories that are provided via \code{prototypes}. Using these distances, trajectories are classified to the prototype trajectories. \code{type=c("hierarchical", "prototypes")} performs both analyses.}
  \item{nclust}{The number of clusters that should be extracted in the hierarchical clustering method}
  \item{nResc}{Before calculating the distance matrix, all trajectories are spatially normalized, i.e. we distribute \code{nResc} equally (spatially) spaced points on each trajectory. This improves clustering performance, as most points in movement trajectories are at the start- and endpoint of the trajectory, which are relatively uninformative with regards to the shape of the trajectory.}
  \item{prototypes}{A list containing prototypical trajectories, to which the method  \code{type=c("prototypes")} classifies all trajectories. Each prototype has to be specified as a n x 2 - matrix, with x values in the first and y values in the second column.}
  \item{subsampN}{Takes a random subsample from the original data and performs the clustering analysis on this subsample. This is useful for datasets with a large number of trajectories, which could render hierarchical clustering computationally infeasible.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a list containing:
\item{call}{The function input except the data.}
\item{data_res}{A data frame containing the spatially normalized data.}
\item{hierarchical}{The results of the hierarchical cluster analysis. \code{distmat} contains the distance matrix used for clustering. \code{clust_obj} contains the \code{hclust} object from the \code{fastcluster} package. Using this object, you can compute clustering solutions with a different number of clusters, without rerunning the whole function. \code{cluster} is a vector containing the classification of each trajectory.}
\item{prototypes}{Contains the results of the prototype classification. \code{prototypes_resc} contains the spatially normalized prototypes (Both the empirical trajectories and the prototypes have to have the same number of data points to compute a difference). \code{protoclust} contains a data frame that specifies the distance of each trajectory to each prototypical trajectory. Also, every trajectory is calssified to the closest prototype.}

}
\references{
Spivey, M. J., Grosjean, M., & Knoblich, G. (2005). Continuous attraction toward phonological competitors. Proceedings of the National Academy of Sciences of the United States of America, 102(29), 10393-10398.
}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

\dontrun{

# THIS EXAMPLE DOES NOT RUN ANYMORE

# we use a part of the following example dataset

head(data_sp2015)

# we use prepr() to time normalize all trajectories to 101 time-steps and strech them to a norm display

layout_stretch <- list("start"=c(0,0), "left"=c(-1,1.5), "right"=c(1,1.5))

output <- prepr(data = data_sp2015[1:1000,], 
                i.xyt = c('x', 'y', 't'), 
                i.id =  c('id.ptp', 'id.trial'), 
                type = "time", 
                steps = 101, 
                start2zero = TRUE, 
                stretch = layout_stretch)

head(output$data)

# clustering; prototypes from example datasets "prototypes"

out_clust <- trajcluster(output$data, 
                        i.xyt =  c('x', 'y', 't'), 
                        i.id = c('id.ptp', 'id.trial'),
                        type=c("hierarchical", "prototypes"),
                        nclust = 4, 
                        nResc = 10, 
                        prototypes = prototypes, 
                        subsampN = NA)

}

}

