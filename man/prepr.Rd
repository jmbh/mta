
\name{prepr}

\alias{prepr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Preprocessing computer mouse trajectories
}
\description{
Function to normalize trajectories across time or space and to strech trajectory to a norm display.}
}
\usage{
prepr(data, i.xyt, i.id, layout, type, steps = 101, start2zero = TRUE, stretch = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame containing the x, y, t variables and indicator variables which uniquely identify single trajectories}
  \item{i.xyt}{
Vector with the column names of the x, y, t variables in that order.
}
  \item{i.id}{
Vector with the column names of the indicator variables that uniquely identify single trajectories (e.g. c('experiment1', 'trial')).
}
  \item{layout}{
A list that specifyies the start- and end points in the experimental setup. See examples for the exact specification.
}
  \item{steps}{The number of steps that should intrapolated}
  \item{start2zero}{If TRUE, the startpoint of each trajectory is set to (0,0)}
  \item{stretch}{Each trajectory is 'streched' into a specified norm-display. The norm display is specified as \code{layout}. This is useful for comparing data from various studies with different layout.}
  
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Spivey, M. J., Grosjean, M., & Knoblich, G. (2005). Continuous attraction toward phonological competitors. Proceedings of the National Academy of Sciences of the United States of America, 102(29), 10393-10398.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# we use the following example data:  
head(dataraw)

# we use prepr() to time normalize all trajectories to 101 time-steps and strech them to a norm display

layout <- list("start"=c(960,230), "left"=c(130,905), "right"=c(1830,905)) 
layout_strech <- list("start"=c(0,0), "left"=c(-1,1.5), "right"=c(1,1.5))

output <- prepr(data = dataraw, 
                i.xyt = c('x', 'y', 't'), 
                i.id = 'trial', 
                layout, 
                type = "time", 
                steps = 101, 
                start2zero = TRUE, 
                stretch = layout_strech)

head(output$data)

}

