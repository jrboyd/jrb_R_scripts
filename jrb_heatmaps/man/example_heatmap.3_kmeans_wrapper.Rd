\name{example_heatmap.3_kmeans_wrapper}
\alias{example_heatmap.3_kmeans_wrapper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
example_heatmap.3_kmeans_wrapper(n_row = 2000, n_col = 6, n_spikes = 4, n_clust = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n_row}{
%%     ~~Describe \code{n_row} here~~
}
  \item{n_col}{
%%     ~~Describe \code{n_col} here~~
}
  \item{n_spikes}{
%%     ~~Describe \code{n_spikes} here~~
}
  \item{n_clust}{
%%     ~~Describe \code{n_clust} here~~
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
function (n_row = 2000, n_col = 6, n_spikes = 4, n_clust = 5) 
{
    MIN = 95
    MAX = 105
    set.seed(1)
    dat = matrix(runif(n_row * n_col, min = MIN, max = MAX), 
        nrow = n_row, ncol = n_col)
    rand = function(n, of) {
        order(runif(of))[1:n]
    }
    for (i in 1:n_spikes) {
        n_choose = (1/n_spikes) * n_row
        rows_key = rand(n_choose, n_row)
        cols_key = rand(2, of = n_col)
        bump = runif(1, min = 10, max = 30)
        print(round(bump))
        if (runif(1) > 0.5) {
            bump = -bump
        }
        dat[rows_key, cols_key] = dat[rows_key, cols_key] + bump
    }
    heatmap.3_kmeans_wrapper(dat = dat, nclust = n_clust)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
