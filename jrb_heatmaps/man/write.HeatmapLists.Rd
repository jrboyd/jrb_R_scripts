\name{write.HeatmapLists}
\alias{write.HeatmapLists}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
write.HeatmapLists(heatmap3_output, name = "heatmap_lists", output_directory = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{heatmap3_output}{
%%     ~~Describe \code{heatmap3_output} here~~
}
  \item{name}{
%%     ~~Describe \code{name} here~~
}
  \item{output_directory}{
%%     ~~Describe \code{output_directory} here~~
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
function (heatmap3_output, name = "heatmap_lists", output_directory = "") 
{
    res = heatmap3_output
    print(output_directory)
    if (!file.exists(output_directory)) 
        dir.create(output_directory)
    root = paste(output_directory, name, sep = "/")
    print(root)
    classSizes = res[[1]]
    asPlotted = res[[3]]
    classColors = res[[4]]
    data_out = matrix(0, nrow = sum(classSizes), ncol = 3)
    data_out[, 1:2] = as.matrix(indexDict[rownames(asPlotted), 
        c(1, 4)])
    colors_out = matrix("black", nrow = sum(classSizes), ncol = 3)
    soFar = 0
    for (j in 1:length(classSizes)) {
        size = classSizes[j]
        start = soFar + 1
        end = soFar + size
        soFar = soFar + size
        data_out[start:end, 3] = paste("Class", j)
        colors_out[start:end, 1:3] = classColors[j]
    }
    pdf(paste(root, ".pdf", sep = ""))
    textplot(data_out[, 1:2], col.data = colors_out[, 1:2], show.rownames = F, 
        show.colnames = F, halign = "left", hadj = 0)
    dev.off()
    png(paste(root, ".png", sep = ""))
    textplot(as.matrix(data_out[, 1]), col.data = as.matrix(colors_out[, 
        3]), show.rownames = F, show.colnames = F)
    dev.off()
    write.table(x = sub(data_out, pattern = ",", replacement = " "), 
        file = paste(root, ".csv", sep = ""), quote = F, row.names = F, 
        col.names = F, sep = ",")
    rownames(asPlotted) = indexDict[rownames(asPlotted), 1]
    asPlotted = cbind(rownames(asPlotted), asPlotted)
    colnames(asPlotted)[1] = "Symbol"
    asPlotted = unlist(strsplit(asPlotted[, 1], split = ","))
    write.table(x = asPlotted, file = paste(root, "_ipa.txt", 
        sep = ""), col.names = F, row.names = F, quote = F, sep = "\t")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
