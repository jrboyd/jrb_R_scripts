\name{textplot.matrix}
\alias{textplot.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
textplot.matrix(object, halign = c("center", "left", "right"), valign = c("center", "top", "bottom"), cex, cmar = 2, rmar = 0.5, show.rownames = TRUE, show.colnames = TRUE, hadj = 1, vadj = 1, mar = c(1, 1, 4, 1) + 0.1, col.data = par("col"), col.rownames = par("col"), col.colnames = par("col"), ...)
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
  \item{cmar}{
%%     ~~Describe \code{cmar} here~~
}
  \item{rmar}{
%%     ~~Describe \code{rmar} here~~
}
  \item{show.rownames}{
%%     ~~Describe \code{show.rownames} here~~
}
  \item{show.colnames}{
%%     ~~Describe \code{show.colnames} here~~
}
  \item{hadj}{
%%     ~~Describe \code{hadj} here~~
}
  \item{vadj}{
%%     ~~Describe \code{vadj} here~~
}
  \item{mar}{
%%     ~~Describe \code{mar} here~~
}
  \item{col.data}{
%%     ~~Describe \code{col.data} here~~
}
  \item{col.rownames}{
%%     ~~Describe \code{col.rownames} here~~
}
  \item{col.colnames}{
%%     ~~Describe \code{col.colnames} here~~
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
    "top", "bottom"), cex, cmar = 2, rmar = 0.5, show.rownames = TRUE, 
    show.colnames = TRUE, hadj = 1, vadj = 1, mar = c(1, 1, 4, 
        1) + 0.1, col.data = par("col"), col.rownames = par("col"), 
    col.colnames = par("col"), ...) 
{
    if (is.vector(object)) 
        object <- t(as.matrix(object))
    else object <- as.matrix(object)
    if (length(col.data) == 1) 
        col.data <- matrix(col.data, nrow = nrow(object), ncol = ncol(object))
    else if (nrow(col.data) != nrow(object) || ncol(col.data) != 
        ncol(object)) 
        stop("Dimensions of 'col.data' do not match dimensions of 'object'.")
    if (length(col.rownames) == 1) 
        col.rownames <- rep(col.rownames, nrow(object))
    if (length(col.colnames) == 1) 
        if (show.rownames) 
            col.colnames <- rep(col.colnames, ncol(object) + 
                1)
        else col.colnames <- rep(col.colnames, ncol(object))
    halign = match.arg(halign)
    valign = match.arg(valign)
    opar <- par()[c("mar", "xpd", "cex")]
    on.exit(par(opar))
    par(mar = mar, xpd = FALSE)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1), log = "", asp = NA)
    if (is.null(colnames(object))) 
        colnames(object) <- paste("[,", 1:ncol(object), "]", 
            sep = "")
    if (is.null(rownames(object))) 
        rownames(object) <- paste("[", 1:nrow(object), ",]", 
            sep = "")
    if (show.rownames) {
        object <- cbind(rownames(object), object)
        col.data <- cbind(col.rownames, col.data)
    }
    if (show.colnames) {
        object <- rbind(colnames(object), object)
        col.data <- rbind(col.colnames, col.data)
    }
    if (missing(cex)) {
        cex <- 1
        lastloop <- FALSE
    }
    else {
        lastloop <- TRUE
    }
    for (i in 1:20) {
        oldcex <- cex
        width <- sum(apply(object, 2, function(x) max(strwidth(x, 
            cex = cex)))) + strwidth("M", cex = cex) * cmar * 
            ncol(object)
        height <- strheight("M", cex = cex) * nrow(object) * 
            (1 + rmar)
        if (lastloop) 
            break
        cex <- cex/max(width, height)
        if (abs(oldcex - cex) < 0.001) {
            lastloop <- TRUE
        }
    }
    rowheight <- strheight("W", cex = cex) * (1 + rmar)
    colwidth <- apply(object, 2, function(XX) max(strwidth(XX, 
        cex = cex))) + strwidth("W") * cmar
    width <- sum(colwidth)
    height <- rowheight * nrow(object)
    if (halign == "left") 
        xpos <- 0
    else if (halign == "center") 
        xpos <- 0 + (1 - width)/2
    else xpos <- 0 + (1 - width)
    if (valign == "top") 
        ypos <- 1
    else if (valign == "center") 
        ypos <- 1 - (1 - height)/2
    else ypos <- 0 + height
    x <- xpos
    y <- ypos
    xpos <- x
    for (i in 1:ncol(object)) {
        xpos <- xpos + hadj * colwidth[i]
        for (j in 1:nrow(object)) {
            ypos <- y - (j - 1) * rowheight
            if ((show.rownames && i == 1) || (show.colnames && 
                j == 1)) 
                text(xpos, ypos, object[j, i], adj = c(hadj, 
                  vadj), cex = cex, font = 2, col = col.data[j, 
                  i], ...)
            else text(xpos, ypos, object[j, i], adj = c(hadj, 
                vadj), cex = cex, font = 1, col = col.data[j, 
                i], ...)
        }
        xpos <- xpos + (1 - hadj) * colwidth[i]
    }
    par(opar)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
