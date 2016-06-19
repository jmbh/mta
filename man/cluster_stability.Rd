\name{cluster_stability}
\alias{cluster_stability}

\title{
Bootstrap Cluster Stability
}
\description{
Computes cluster the stability for different values of k via the non-parametric bootstrap.}
\usage{
cluster_stability(dist, kseq, B, norm = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dist}{p x p distance matrix, where p is the number of objects.}
  \item{kseq}{A sequence of cluster sizes, for which cluster stability should be computed.}
  \item{B}{Number of bootstrap samples}
  \item{norm}{Default is FALSE and corresponds to the method by Fang & Wang, 2012.}
  \item{...}{Additional arguments passed to hclust.}
  }
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a vector of cluster instability indices, one for each k in \code{kseq}.
}

\references{

Fang, Y., & Wang, J. (2012). Selection of the number of clusters via the bootstrap method. Computational Statistics & Data Analysis, 56(3), 468-477. 

}


\author{
Jonas Haslbeck <jonashaslbeck@gmail.com>
}

\examples{

\dontrun{

# simple Gaussian mixture
data <- c(rnorm(100,0,1), 
          rnorm(100,5,1), 
          rnorm(100,10,1))
hist(data, breaks=40) # look at mixture

# compute distance matrix
dist <- as.matrix(dist(data))

kseq <- 2:10 # define k sequence of interest
set.seed(1) # make reproducible
instobj <- cluster_stability(dist, kseq, B=25)

# visualize instability as a function of k:
plot(kseq, instobj, ylim=c(0,.15), type='l', 
     xlab='k', ylab='Cluster Instability')

# correctly identifies k=3!

}

}

