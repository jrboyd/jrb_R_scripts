\name{plot.HeatmapLists}
\alias{plot.HeatmapLists}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.HeatmapLists(heatmap3_output)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{heatmap3_output}{
%%     ~~Describe \code{heatmap3_output} here~~
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
function (heatmap3_output) 
{
    res = heatmap3_output
    classSizes = res[[1]]
    asPlotted = res[[3]]
    classColors = res[[4]]
    data_out = matrix(rownames(asPlotted), ncol = 1)
    colors_out = matrix("black", nrow = sum(classSizes), ncol = 1)
    soFar = 0
    for (j in 1:length(classSizes)) {
        size = classSizes[j]
        start = soFar + 1
        end = soFar + size
        soFar = soFar + size
        colors_out[start:end, ] = classColors[j]
    }
    colSize = 50
    while (nrow(data_out)\%\%colSize != 0) {
        data_out = rbind(data_out, "")
        colors_out = rbind(colors_out, "black")
    }
    a = matrix(colors_out, nrow = colSize)
    b = matrix(data_out, nrow = colSize)
    step = 6
    while (ncol(a)\%\%step != 0) {
        a = cbind(a, rep("black", nrow(a)))
        b = cbind(b, rep("", nrow(a)))
    }
    for (i in 1:(ncol(a)/6)) {
        start = (i - 1) * 6 + 1
        end = i * 6
        print(start)
        print(end)
        textplot(b[, start:end], col.data = a[, start:end], show.rownames = F, 
            show.colnames = F, halign = "left", hadj = 0, cex = 0.5)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
