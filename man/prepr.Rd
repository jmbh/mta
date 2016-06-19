
\name{prepr}

\alias{prepr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Preprocessing data from 2D mousetracking data in dual choice tasks.
}
\description{
Function normalizes trajectories with respect to time or space.}
\usage{
prepr(data, i.xyt, i.id, type, steps = 101, start2zero = TRUE, stretch = NA, takeAllvar = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing the x, y, t variables and indicator variables which uniquely identify single trajectories}
  \item{i.xyt}{
Vector containing the column names of the x, y, t variables in that order.
}
  \item{i.id}{
Vector containing the column names of the indicator variables that jointly uniquely identify single trajectories (e.g. \code{c('participant_id', 'trial_id')}).
}
  \item{type} { \code{type='time'} normalizes trajectory with respect to time (all trajectories have the same number of steps of time steps). This is the standard data preprocessing procedure for mousetracking. \code{type='spatial'} normalizes trajectories with respect to space, i.e. we place a specified number of equally spaced (with respect to space) points on the trajectory. This is useful for comparing trajectories using a distance measure.}
  \item{steps}{The number of steps that should intrapolated (with respect to time or space, depending on \code{type} )}
  \item{start2zero}{If TRUE, the startpoint of each trajectory is set to (0,0)}
  \item{stretch}{Each trajectory is 'streched' into a specified norm-display. The norm display is specified as a list with the coordinates of the start point and the two possible end points of the trajectories. See the example below. This is useful for comparing data from various studies with different layout.}
  \item{takeAllvar}{Default \code{FALSE}. If \code{TRUE}, the returned dataframe contains all columns of the input dataframe. Otherwise only the columns x,y,t and the id-variables are returned.}
  
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function returns a list:

\item{call}{Returns all function input except the data.}
\item{data}{Returns a data frame with normalized trajectories.}
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

# THIS EXAMPLE DOES NOT WORK ANYMORE !!!


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

}

}

