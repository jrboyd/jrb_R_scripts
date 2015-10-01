\name{replaceTabs.inner}
\alias{replaceTabs.inner}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
replaceTabs.inner(text, width = 8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{text}{
%%     ~~Describe \code{text} here~~
}
  \item{width}{
%%     ~~Describe \code{width} here~~
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
function (text, width = 8) 
{
    spaces <- "        "
    if (nchar(text) < 1) 
        return(text)
    text.split <- strsplit(text, "\t")[[1]]
    if (length(text.split) == 1) 
        return(text)
    else {
        nSpaceAdd <- 8 - nchar(text.split)\%\%8
        nSpaceAdd[length(nSpaceAdd)] <- 0
        nSpaceAdd[nSpaceAdd == 8] <- 0
        retval <- ""
        for (i in 1:length(text.split)) {
            tmp.text <- chartr("\t", " ", text.split[i])
            retval <- paste(retval, tmp.text, substr(spaces, 
                0, nSpaceAdd[i] - 1), sep = "")
        }
        return(retval)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
