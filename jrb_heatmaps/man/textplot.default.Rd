\name{textplot.default}
\alias{textplot.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
textplot.default(object, halign = c("center", "left", "right"), valign = c("center", "top", "bottom"), cex, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{halign}{
%%     ~~Describe \code{halign} here~~
}
  \item{valign}{
%%     ~~Describe \code{valign} here~~
}
  \item{cex}{
%%     ~~Describe \code{cex} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (object, halign = c("center", "left", "right"), valign = c("center", 
    "top", "bottom"), cex, ...) 
{
    if (is.matrix(object) || (is.vector(object) && length(object) > 
        1)) 
        return(textplot.matrix(object, halign, valign, cex, ...))
    halign <- match.arg(halign)
    valign <- match.arg(valign)
    textplot.character(object, halign, valign, cex, ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
