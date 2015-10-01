\name{heatmap.3}
\alias{heatmap.3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
heatmap.3(x, nsplits = 1, override_o = NA, clusterSort = F, Rowv = TRUE, Colv = if (symm) "Rowv" else FALSE, distfun = dist, hclustfun = hclust, dendrogram = c("both", "row", "column", "none"), reorderfun = function(d, w) reorder(d, w), symm = FALSE, classCount = 6, scale = c("none", "row", "column"), minVal = min(x), maxVal = max(x), na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || scale != "none", col = "heat.colors", rowsep, sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan", na.color = par("bg"), trace = c("none", "column", "row", "both"), tracecol = "cyan", hline = median(breaks), vline = median(breaks), linecol = tracecol, margins = c(0, 0), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5, key = TRUE, keysize = 1.5, density.info = c("histogram", "density", "none"), denscol = tracecol, symkey = min(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, key.title = NULL, key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, key.ytickfun = NULL, key.par = list(), main = NULL, xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{nsplits}{
%%     ~~Describe \code{nsplits} here~~
}
  \item{override_o}{
%%     ~~Describe \code{override_o} here~~
}
  \item{clusterSort}{
%%     ~~Describe \code{clusterSort} here~~
}
  \item{Rowv}{
%%     ~~Describe \code{Rowv} here~~
}
  \item{Colv}{
%%     ~~Describe \code{Colv} here~~
}
  \item{distfun}{
%%     ~~Describe \code{distfun} here~~
}
  \item{hclustfun}{
%%     ~~Describe \code{hclustfun} here~~
}
  \item{dendrogram}{
%%     ~~Describe \code{dendrogram} here~~
}
  \item{reorderfun}{
%%     ~~Describe \code{reorderfun} here~~
}
  \item{symm}{
%%     ~~Describe \code{symm} here~~
}
  \item{classCount}{
%%     ~~Describe \code{classCount} here~~
}
  \item{scale}{
%%     ~~Describe \code{scale} here~~
}
  \item{minVal}{
%%     ~~Describe \code{minVal} here~~
}
  \item{maxVal}{
%%     ~~Describe \code{maxVal} here~~
}
  \item{na.rm}{
%%     ~~Describe \code{na.rm} here~~
}
  \item{revC}{
%%     ~~Describe \code{revC} here~~
}
  \item{add.expr}{
%%     ~~Describe \code{add.expr} here~~
}
  \item{breaks}{
%%     ~~Describe \code{breaks} here~~
}
  \item{symbreaks}{
%%     ~~Describe \code{symbreaks} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{rowsep}{
%%     ~~Describe \code{rowsep} here~~
}
  \item{sepcolor}{
%%     ~~Describe \code{sepcolor} here~~
}
  \item{sepwidth}{
%%     ~~Describe \code{sepwidth} here~~
}
  \item{cellnote}{
%%     ~~Describe \code{cellnote} here~~
}
  \item{notecex}{
%%     ~~Describe \code{notecex} here~~
}
  \item{notecol}{
%%     ~~Describe \code{notecol} here~~
}
  \item{na.color}{
%%     ~~Describe \code{na.color} here~~
}
  \item{trace}{
%%     ~~Describe \code{trace} here~~
}
  \item{tracecol}{
%%     ~~Describe \code{tracecol} here~~
}
  \item{hline}{
%%     ~~Describe \code{hline} here~~
}
  \item{vline}{
%%     ~~Describe \code{vline} here~~
}
  \item{linecol}{
%%     ~~Describe \code{linecol} here~~
}
  \item{margins}{
%%     ~~Describe \code{margins} here~~
}
  \item{ColSideColors}{
%%     ~~Describe \code{ColSideColors} here~~
}
  \item{RowSideColors}{
%%     ~~Describe \code{RowSideColors} here~~
}
  \item{cexRow}{
%%     ~~Describe \code{cexRow} here~~
}
  \item{cexCol}{
%%     ~~Describe \code{cexCol} here~~
}
  \item{labRow}{
%%     ~~Describe \code{labRow} here~~
}
  \item{labCol}{
%%     ~~Describe \code{labCol} here~~
}
  \item{srtRow}{
%%     ~~Describe \code{srtRow} here~~
}
  \item{srtCol}{
%%     ~~Describe \code{srtCol} here~~
}
  \item{adjRow}{
%%     ~~Describe \code{adjRow} here~~
}
  \item{adjCol}{
%%     ~~Describe \code{adjCol} here~~
}
  \item{offsetRow}{
%%     ~~Describe \code{offsetRow} here~~
}
  \item{offsetCol}{
%%     ~~Describe \code{offsetCol} here~~
}
  \item{key}{
%%     ~~Describe \code{key} here~~
}
  \item{keysize}{
%%     ~~Describe \code{keysize} here~~
}
  \item{density.info}{
%%     ~~Describe \code{density.info} here~~
}
  \item{denscol}{
%%     ~~Describe \code{denscol} here~~
}
  \item{symkey}{
%%     ~~Describe \code{symkey} here~~
}
  \item{densadj}{
%%     ~~Describe \code{densadj} here~~
}
  \item{key.title}{
%%     ~~Describe \code{key.title} here~~
}
  \item{key.xlab}{
%%     ~~Describe \code{key.xlab} here~~
}
  \item{key.ylab}{
%%     ~~Describe \code{key.ylab} here~~
}
  \item{key.xtickfun}{
%%     ~~Describe \code{key.xtickfun} here~~
}
  \item{key.ytickfun}{
%%     ~~Describe \code{key.ytickfun} here~~
}
  \item{key.par}{
%%     ~~Describe \code{key.par} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
  \item{xlab}{
%%     ~~Describe \code{xlab} here~~
}
  \item{ylab}{
%%     ~~Describe \code{ylab} here~~
}
  \item{lmat}{
%%     ~~Describe \code{lmat} here~~
}
  \item{lhei}{
%%     ~~Describe \code{lhei} here~~
}
  \item{lwid}{
%%     ~~Describe \code{lwid} here~~
}
  \item{extrafun}{
%%     ~~Describe \code{extrafun} here~~
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
function (x, nsplits = 1, override_o = NA, clusterSort = F, Rowv = TRUE, 
    Colv = if (symm) "Rowv" else FALSE, distfun = dist, hclustfun = hclust, 
    dendrogram = c("both", "row", "column", "none"), reorderfun = function(d, 
        w) reorder(d, w), symm = FALSE, classCount = 6, scale = c("none", 
        "row", "column"), minVal = min(x), maxVal = max(x), na.rm = TRUE, 
    revC = identical(Colv, "Rowv"), add.expr, breaks, symbreaks = min(x < 
        0, na.rm = TRUE) || scale != "none", col = "heat.colors", 
    rowsep, sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, 
    notecex = 1, notecol = "cyan", na.color = par("bg"), trace = c("none", 
        "column", "row", "both"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(0, 
        0), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    srtRow = NULL, srtCol = NULL, adjRow = c(0, NA), adjCol = c(NA, 
        0), offsetRow = 0.5, offsetCol = 0.5, key = TRUE, keysize = 1.5, 
    density.info = c("histogram", "density", "none"), denscol = tracecol, 
    symkey = min(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, 
    key.title = NULL, key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, 
    key.ytickfun = NULL, key.par = list(), main = NULL, xlab = NULL, 
    ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, 
    ...) 
{
    scale01 <- function(x, low = minVal, high = maxVal) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    if (clusterSort) 
        dendrogram = "none"
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram \%in\% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dendrogram <- "none"
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram \%in\% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
        if (length(rowInd) > nr || any(rowInd < 1 | rowInd > 
            nr)) 
            stop("Rowv dendrogram doesn't match size of x")
    }
    else if (is.integer(Rowv)) {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (length(override_o) == 1 && is.na(override_o)) {
        dataClasses = cutree(tree = hcr, k = classCount)
        orderOfClasses = unique(dataClasses[rowInd])
        tmp = classCount:1
        tmp = tmp[order(orderOfClasses)]
        require("RColorBrewer")
        colorClasses = brewer.pal(n = classCount, name = "Dark2")
        colorClasses = colorClasses[tmp]
        RowSideColors = colorClasses[dataClasses]
    }
    else {
        dataClasses = integer()
        hidden = sapply(1:nrow(override_o), function(x) {
            dataClasses <<- c(dataClasses, rep(override_o[x, 
                1], override_o[x, 2]))
        })
        classCount = nrow(override_o)
        colorClasses = brewer.pal(n = classCount, name = "Dark2")
        RowSideColors = colorClasses[dataClasses]
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
        if (length(colInd) > nc || any(colInd < 1 | colInd > 
            nc)) 
            stop("Colv dendrogram doesn't match size of x")
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    if (clusterSort) {
        for (i in 1:classCount) {
            keep = dataClasses[rowInd] == i
            subset = x[rowInd, , drop = F][keep, , drop = F]
            in_cluster_o = order(rowSums(subset))
            x[rowInd, ][keep, ] = subset[in_cluster_o, ]
        }
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(minVal, na.rm = na.rm), max(maxVal, 
                na.rm = na.rm), length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 2, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 2, c(rep(NA, nrow(lmat) - 
                1), 2))
            lwid <- c(lwid[1], 0.2, lwid[2], 0.2)
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    linePlotRegion = matrix(c(0, rep(x = 1, classCount) + max(lmat), 
        0, 2:(classCount + 1) + max(lmat)), nrow = classCount + 
        1, ncol = 2, byrow = FALSE)
    combinedRegions = matrix(0, ncol = (ncol(lmat) + ncol(linePlotRegion)), 
        nrow = nrow(linePlotRegion))
    combinedRegions[1, 1:ncol(lmat)] = lmat[1, ]
    for (i in 2:nrow(combinedRegions)) combinedRegions[i, 1:ncol(lmat)] = lmat[2, 
        ]
    combinedRegions[, (ncol(lmat) + 1):ncol(combinedRegions)] = linePlotRegion
    lwid = c(lwid, 1, 1)
    lhei = c(1.5, rep(4/classCount, classCount))
    if (margins[1] > 0) {
        combinedRegions = cbind(combinedRegions, rep(0, nrow(combinedRegions)))
        lwid = c(lwid, margins[1])
    }
    if (margins[2] > 0) {
        combinedRegions = rbind(combinedRegions, rep(0, ncol(combinedRegions)))
        lhei = c(lhei, margins[2])
    }
    nf = layout(combinedRegions, widths = lwid, heights = lhei, 
        respect = TRUE)
    if (!missing(RowSideColors)) {
        par(mar = c(0, 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
        par(mar = c(0, 0.5, 0, 0))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, 0))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(0, 0, 0, 0))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    if (is.null(srtCol)) 
        axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
            offsetCol, tick = 0, cex.axis = 0.45 * cexCol, hadj = adjCol[1], 
            padj = adjCol[2])
    else {
        if (is.numeric(srtCol)) {
            if (missing(adjCol) || is.null(adjCol)) 
                adjCol = c(1, NA)
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
                strheight("M"), labels = labCol, adj = adjCol, 
                cex = cexCol * 0.45, srt = srtCol)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtCol ignored.")
    }
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = 0 - 1.25)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = 0 - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if ((ncol(t(x))\%\%nsplits) != 0) {
        print("ncol not divisible by num of splits! ignoring splits.")
        nsplits = 1
    }
    win = ncol(t(x))/nsplits
    if (nsplits > 1) {
        colsep = 1:(nsplits - 1)
        colsep = colsep * win
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    }
    sepcolor = "black"
    o = matrix(0, nrow = classCount, ncol = 2)
    colnames(o) = c("class index", "class size")
    ordered = dataClasses[rev(rowInd)]
    for (i in 1:classCount) {
        totalLength = length(ordered)
        classIndex = ordered[1]
        keep = ordered != ordered[1]
        ordered = ordered[keep]
        classSize = totalLength - length(ordered)
        o[i, ] = c(classIndex, classSize)
    }
    if (length(override_o) == 1 && !is.na(override_o)) {
        colnames(override_o) = colnames(o)
        o = override_o
    }
    rowsep = 1:(classCount - 1)
    for (i in 1:(classCount - 1)) {
        rowsep[i] = sum(o[1:i, 2])
    }
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace \%in\% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace \%in\% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i - 0.5 + hline.vals, col = linecol, 
                  lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(0, 0, 0, 0))
    if (dendrogram \%in\% c("both", "row")) {
        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, 0))
    if (dendrogram \%in\% c("both", "column")) {
        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.1 * op[["cex.main"]])
    if (T) {
        mar <- c(3.5, 3.5, 5, 2)
        if (!is.null(key.xlab) && is.na(key.xlab)) 
            mar[1] <- 2
        if (!is.null(key.title) && is.na(key.title)) 
            mar[3] <- 1
        par(mar = mar, cex = 0.5, mgp = c(2.1, 1, 0))
        if (length(key.par) > 0) 
            do.call(par, key.par)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(minVal, na.rm = TRUE)
            max.raw <- max(maxVal, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        if (is.null(key.xtickfun)) {
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            xargs <- list(at = xv, labels = lv)
        }
        else {
            xargs <- key.xtickfun()
        }
        xargs$side <- 1
        do.call(axis, xargs)
        if (is.null(key.xlab)) {
            if (scale == "row") 
                key.xlab <- "Row Z-Score"
            else if (scale == "column") 
                key.xlab <- "Column Z-Score"
            else key.xlab <- "Value"
        }
        if (!is.na(key.xlab)) {
            mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5, 
                cex = 0.8)
        }
        title("Color Key")
    }
    else plot.new()
    revx = t(x)
    revx = revx[rev(1:nrow(revx)), ]
    print(o)
    fetchClass = function(n) {
        start = 1
        if (n > 1) 
            start = sum(o[1:(n - 1), 2]) + 1
        end = sum(o[1:n, 2])
        c1 = as.matrix(revx[start:end, , drop = F])
        nc1 = c1[1, ]
        if (nrow(c1) > 1) 
            nc1 = t(apply(c1, MARGIN = 1, FUN = function(row) return(row/(maxVal - 
                minVal))))
        mc1 = c1
        if (nrow(c1) > 1) 
            mc1 = colMeans(c1)
        return(list(c1, nc1, mc1))
    }
    avgA = matrix(0, classCount, ncol(t(x)))
    for (i in 1:nrow(avgA)) {
        avgA[i, ] = fetchClass(i)[[3]]
    }
    par(mai = rep(0, 4))
    plot(x = c(0, 1), frame.plot = FALSE, y = c(0, 1), xaxt = "n", 
        yaxt = "n", type = "n", xlab = "", ylab = "", xlim = c(0, 
            1), ylim = c(0, 1))
    min = -0.04
    max = 1.04
    rng = max - min
    par(mai = rep(0, 4))
    lines(c(min, max), c(max, max), lty = 2)
    for (i in 1:classCount) {
        heatmapFraction = sum(o[1:i, 2])/sum(o[, 2])
        lineplotFraction = i/classCount
        lines(c(min, max), c(max - rng * heatmapFraction, max - 
            rng * lineplotFraction), lty = 2)
    }
    par(mar = c(0, 0, 0, 0.5))
    for (i in 1:classCount) {
        xrange <- as.numeric(range(1:ncol(avgA)))
        yrange <- range(avgA[i, ])
        xspace = (xrange[2] - xrange[1]) * 0.15
        yspace = (yrange[2] - yrange[1]) * 0.15
        plot(xrange, yrange, xaxt = "n", yaxt = "n", type = "n", 
            xlab = "", ylab = "", ylim = c(min(yrange) - yspace, 
                max(yrange) + yspace), xlim = c(min(xrange) - 
                xspace, max(xrange) + xspace))
        colors <- colorClasses[o[, 1]]
        linetype <- c(1:classCount)
        plotchar <- seq(18, 18 + classCount, 1)
        vals <- avgA[i, ]
        for (s in 1:nsplits) {
            start = (s - 1) * win + 1
            end = s * win
            xs = start:end
            lines(xs, vals[xs], type = "b", lwd = 2.5, lty = 1, 
                col = colors[i], pch = plotchar[3])
        }
    }
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    if (!is.null(extrafun)) 
        extrafun()
    x = t(x)
    x = x[rev(1:nrow(x)), ]
    rowInd = rev(rowInd)
    return(list(o[, 2], rowInd, x, colors))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
