\name{textplot.character}
\alias{textplot.character}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
textplot.character(object, halign = c("center", "left", "right"), valign = c("center", "top", "bottom"), cex, fixed.width = TRUE, cspace = 1, lspace = 1, mar = c(0, 0, 3, 0) + 0.1, tab.width = 8, ...)
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
  \item{fixed.width}{
%%     ~~Describe \code{fixed.width} here~~
}
  \item{cspace}{
%%     ~~Describe \code{cspace} here~~
}
  \item{lspace}{
%%     ~~Describe \code{lspace} here~~
}
  \item{mar}{
%%     ~~Describe \code{mar} here~~
}
  \item{tab.width}{
%%     ~~Describe \code{tab.width} here~~
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
    "top", "bottom"), cex, fixed.width = TRUE, cspace = 1, lspace = 1, 
    mar = c(0, 0, 3, 0) + 0.1, tab.width = 8, ...) 
{
    object <- paste(object, collapse = "\n", sep = "")
    object <- replaceTabs(object, width = tab.width)
    halign = match.arg(halign)
    valign = match.arg(valign)
    plot.new()
    opar <- par()[c("mar", "xpd", "cex", "family")]
    on.exit(par(opar))
    par(mar = mar, xpd = FALSE)
    if (fixed.width) 
        par(family = "mono")
    plot.window(xlim = c(0, 1), ylim = c(0, 1), log = "", asp = NA)
    slist <- unlist(lapply(object, function(x) strsplit(x, "\n")))
    slist <- lapply(slist, function(x) unlist(strsplit(x, "")))
    slen <- sapply(slist, length)
    slines <- length(slist)
    if (missing(cex)) {
        lastloop <- FALSE
        cex <- 1
    }
    else lastloop <- TRUE
    for (i in 1:20) {
        oldcex <- cex
        cwidth <- max(sapply(unlist(slist), strwidth, cex = cex)) * 
            cspace
        cheight <- max(sapply(unlist(slist), strheight, cex = cex)) * 
            (lspace + 0.5)
        width <- strwidth(object, cex = cex)
        height <- strheight(object, cex = cex)
        if (lastloop) 
            break
        cex <- cex/max(width, height)
        if (abs(oldcex - cex) < 0.001) {
            lastloop <- TRUE
        }
    }
    if (halign == "left") 
        xpos <- 0
    else if (halign == "center") 
        xpos <- 0 + (1 - width)/2
    else xpos <- 0 + (1 - width)
    if (valign == "top") 
        ypos <- 1
    else if (valign == "center") 
        ypos <- 1 - (1 - height)/2
    else ypos <- 1 - (1 - height)
    text(x = xpos, y = ypos, labels = object, adj = c(0, 1), 
        cex = cex, ...)
    par(opar)
    invisible(cex)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
