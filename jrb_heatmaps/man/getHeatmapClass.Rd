\name{getHeatmapClass}
\alias{getHeatmapClass}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getHeatmapClass(resultsFromHeatmap, classIndex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultsFromHeatmap}{
%%     ~~Describe \code{resultsFromHeatmap} here~~
}
  \item{classIndex}{
%%     ~~Describe \code{classIndex} here~~
}
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
%% ~put references to the literature/web site here ~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (resultsFromHeatmap, classIndex) 
{
    res = resultsFromHeatmap
    sizes = res[[1]]
    o = res[[2]]
    asPlotted = res[[3]]
    begin = 1
    end = sizes[1]
    if (classIndex > 1) {
        begin = sum(sizes[1:(classIndex - 1)]) + 1
        end = sum(sizes[1:(classIndex)])
    }
    return(rownames(asPlotted)[begin:end])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
