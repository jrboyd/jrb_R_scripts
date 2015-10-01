\name{heatmap.3_kmeans_wrapper}
\alias{heatmap.3_kmeans_wrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
heatmap.3_kmeans_wrapper(dat, nclust = 4, hmap_colors = c("blue", "white", "red"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{nclust}{
%%     ~~Describe \code{nclust} here~~
}
  \item{hmap_colors}{
%%     ~~Describe \code{hmap_colors} here~~
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
function (dat, nclust = 4, hmap_colors = c("blue", "white", "red"), 
    ...) 
{
    set.seed(1)
    kclust = kmeans(dat, centers = nclust)
    o = order(kclust$cluster)
    kclust$cluster = kclust$cluster[o]
    dat = dat[o, ]
    for (i in 1:nclust) {
        keep = kclust$cluster == i
        subset = dat[keep, ]
        o = order(apply(subset, 1, function(x) {
            return(max(x) - min(x))
        }), decreasing = T)
        rownames(dat)[keep] = rownames(subset)[o]
        dat[keep, ] = subset[o, ]
    }
    o = order(rowSums(kclust$centers), decreasing = T)
    kclust$centers = kclust$centers[o, ]
    plot_dat = matrix(0, nrow = 0, ncol = ncol(dat))
    new_cluster = integer()
    colnames(plot_dat) = colnames(dat)
    j = 1
    for (i in o) {
        keep = kclust$cluster == i
        subset = dat[keep, ]
        plot_dat = rbind(plot_dat, subset)
        new_cluster = c(new_cluster, rep(j, sum(keep)))
        j = j + 1
    }
    clust_sizes = sapply(o, function(x) {
        sum(kclust$cluster == x)
    })
    cr = colorRamp(c("blue", "white", "red"))
    colors = rgb(cr(0:100/100)/255)
    seps = cumsum(sapply(1:nclust, function(x) sum(x == new_cluster)))
    override_o = cbind(1:nclust, clust_sizes)
    heatmap.3(plot_dat, trace = "n", Rowv = F, Colv = F, scale = "n", 
        cexCol = 1.6, cexRow = 0.4, col = colors, density.info = "n", 
        key.xlab = "z-score of log gene norm-counts", key.title = "", 
        labRow = "", override_o = override_o, nsplits = 1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
