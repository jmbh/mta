\name{gap_statistic}
\alias{gap_statistic}

\title{
Computing the Gap Statistic
}
\description{
Computing within cluster dissimilarity for different numbers of cluster k and the Gap statistic.}
\usage{
gap_statistic(dist, kseq, steps, lambda = .7, bezier = NULL, 
              xcor = c(0,1,-1),  ycor = c(0,1.5,1.5))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist}{p x p distance matrix, where p is the number of trajectories.}
  \item{kseq}{A sequence of cluster sized, for which the within cluster dissimilarity should be computed}
  \item{steps}{The number of data points on the generated trajectories.}
  \item{lambda}{Lambda parameter of the exponential distribution, from which the parameters for the Bezier curve are drawn, which we use to generate trajectories. To get a feel for this parameter, try out values and look at the generated data.}
  \item{bezier}{Provide a custom sequence of Bezier curve parameters.}
\item{xcor}{x-coordinates of the start point, the non-selected box and the end-point in that order.}
\item{ycor}{y-coordinates of the start point, the non-selected box and the end-point in that order.}

  }
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a list:
\item{WCD_data}{Within cluster dissimilarity for each k for the real data.}
\item{WCD_syn}{Within cluster dissimilarity for each k for the generated data.}
\item{Gaps}{The difference values of \code{WCD_syn} - \code{WCD_data}}
\item{cl_data}{Cluster assignments for each k}
\item{cl_syn}{Cluster assignments for each k and the generated data}
}
\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}

\examples{

\dontrun{

# Example data
head(data_sp2015)

# Preprocess Data
stre <- list("start"=c(0,0), "left"=c(-1,1.5), "right"=c(1,1.5))
steps <- 5
data_pr <- prepr(data_sp2015, i.xyt = c('x', 'y', 't'), i.id=c('id.ptp', 'id.trial'), 
                  type='spatial', stretch = stre, steps = steps)$data
data_pr$x[data_pr$chosen.box==1] <- data_pr$x[data_pr$chosen.box==1] * -1 # all endpoints to one side

# Compute distance matrix
id <- 1:1344 # no of unique trials
dist <- mta:::distmat(id, data_pr$x , data_pr$y, steps)

# Call Gap statistic function
gapobj <- gap_statistic(dist, kseq=1:10, steps=5)



}

}

