## $Id: heatmap.2.R 1823 2014-06-30 19:26:21Z warnes $
library(gplots)
library(gtools)
getHeatmapClass = function(resultsFromHeatmap, classIndex){
  res = resultsFromHeatmap
  sizes = res[[1]]
  o = res[[2]]
  asPlotted = res[[3]]
  begin = 1
  end = sizes[1]
  if(classIndex > 1){
    begin = sum(sizes[1:(classIndex-1)]) + 1
    end = sum(sizes[1:(classIndex)])
  }
  #asPlotted = asPlotted[o,]
  return(rownames(asPlotted)[begin:end])
}


heatmap.3 <- function (x,
                       #returns list of classsizes, row order vector for input x, the data as plotted in heatmap, all with row 1 at top of heatmap and last row at bottom of heatmap
                       ## set number of splits
                       nsplits = 1,
                       ## dendrogram control
                       clusterSort = F,
                       Rowv = TRUE,
                       Colv=if(symm)"Rowv" else FALSE,
                       distfun = dist,
                       hclustfun = hclust,
                       dendrogram = c("both","row","column","none"),
                       reorderfun = function(d, w) reorder(d, w),
                       symm = FALSE,
                       classCount = 6,
                       
                       ## data scaling
                       scale = c("none","row", "column"),
                       minVal = min(x),
                       maxVal = max(x),
                       na.rm=TRUE,
                       
                       ## image plot
                       revC = identical(Colv, "Rowv"),
                       add.expr,
                       
                       ## mapping data to colors
                       breaks,
                       symbreaks=min(x < 0, na.rm=TRUE) || scale!="none",
                       
                       ## colors
                       col="heat.colors",
                       
                       ## block sepration
                       #colsep,
                       rowsep,
                       sepcolor="white",
                       sepwidth=c(0.05,0.05),
                       
                       ## cell labeling
                       cellnote,
                       notecex=1.0,
                       notecol="cyan",
                       na.color=par("bg"),
                       
                       ## level trace
                       trace=c("none","column","row","both"),
                       tracecol="cyan",
                       hline=median(breaks),
                       vline=median(breaks),
                       linecol=tracecol,
                       
                       ## Row/Column Labeling
                       margins = c(0, 0),
                       ColSideColors,
                       RowSideColors,
                       cexRow = 0.2 + 1/log10(nr),
                       cexCol = 0.2 + 1/log10(nc),
                       labRow = NULL,
                       labCol = NULL,
                       srtRow = NULL,
                       srtCol = NULL,
                       adjRow = c(0,NA),
                       adjCol = c(NA,0),
                       offsetRow = 0.5,
                       offsetCol = 0.5,
                       
                       ## color key + density info
                       key = TRUE,
                       keysize = 1.5,
                       density.info=c("histogram","density","none"),
                       denscol=tracecol,
                       symkey = min(x < 0, na.rm=TRUE) || symbreaks,
                       densadj = 0.25,
                       key.title = NULL,
                       key.xlab = NULL,
                       key.ylab = NULL,
                       key.xtickfun = NULL,
                       key.ytickfun = NULL,
                       key.par=list(),
                       
                       ## plot labels
                       main = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       
                       ## plot layout
                       lmat = NULL,
                       lhei = NULL,
                       lwid = NULL,
                       
                       ## extras
                       extrafun=NULL,
                       ...
)
{
  scale01 <- function(x, low=minVal, high=maxVal )
  {
    x <- (x-low)/(high - low)
    x
  }
  
  retval <- list()
  
  scale <- if(symm && missing(scale)) "none" else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  if(clusterSort) dendrogram = 'none'
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  
  if(length(col)==1 && is.character(col) )
    col <- get(col, mode="function")
  
  if(!missing(breaks) && (scale!="none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.",
            "Please consider using only one or the other.")
  
  if ( is.null(Rowv) || is.na(Rowv) )
    Rowv <- FALSE
  if ( is.null(Colv) || is.na(Colv) )
    Colv <- FALSE
  else if( Colv=="Rowv" && !isTRUE(Rowv) )
    Colv <- FALSE
  
  
  if(length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  
  nr <- di[1]
  nc <- di[2]
  
  if(nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  
  if(!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  
  if(missing(cellnote))
    cellnote <- matrix("", ncol=ncol(x), nrow=nrow(x))
  
  if(!inherits(Rowv, "dendrogram")) {
    ## Check if Rowv and dendrogram arguments are consistent
    if ( ( (!isTRUE(Rowv)) || (is.null(Rowv))) &&
           (dendrogram %in% c("both","row") ) )
    {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else
        dendrogram <- "none"
      
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
      
    }
  }
  
  if(!inherits(Colv, "dendrogram")) {
    ## Check if Colv and dendrogram arguments are consistent
    if ( ( (!isTRUE(Colv)) || (is.null(Colv)))
         && (dendrogram %in% c("both","column")) )
    {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else
        dendrogram <- "none"
      #jrb this warning makes no sense to me - gone
      #         warning("Discrepancy: Colv is FALSE, while dendrogram is `",
      #                 dendrogram, "'. Omitting column dendogram.")
    }
  }
  
  
  ## by default order by row/col mean
  ## if(is.null(Rowv)) Rowv <- rowMeans(x, na.rm = na.rm)
  ## if(is.null(Colv)) Colv <- colMeans(x, na.rm = na.rm)
  
  ## get the dendrograms and reordering indices
  
  ## if( dendrogram %in% c("both","row") )
  ##  { ## dendrogram option is used *only* for display purposes
  if(inherits(Rowv, "dendrogram"))
  {
    ddr <- Rowv ## use Rowv 'as-is', when it is dendrogram
    rowInd <- order.dendrogram(ddr)
    if(length(rowInd)>nr || any(rowInd<1 | rowInd > nr ))
      stop("Rowv dendrogram doesn't match size of x")
  }
  else if (is.integer(Rowv))
  { ## Compute dendrogram and do reordering based on given vector
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <-  reorderfun(ddr, Rowv)
    
    rowInd <- order.dendrogram(ddr)
    
    
    
    if(nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv))
  { ## If TRUE, compute dendrogram and do reordering based on rowMeans
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)
    
    rowInd <- order.dendrogram(ddr)
    if(nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  } else {
    rowInd <- nr:1
  }
  #jrb row side colors is default and get data classes
  dataClasses = cutree(tree=hcr, k=classCount)
  orderOfClasses = unique(dataClasses[rowInd])
  tmp = classCount:1
  tmp = tmp[order(orderOfClasses)]
  require('RColorBrewer')
  colorClasses = brewer.pal(n = classCount,name = "Dark2")
  colorClasses = colorClasses[tmp]
  RowSideColors = colorClasses[dataClasses]
  ## if( dendrogram %in% c("both","column") )
  ##   {
  if(inherits(Colv, "dendrogram"))
  {
    ddc <- Colv ## use Colv 'as-is', when it is dendrogram
    colInd <- order.dendrogram(ddc)
    if(length(colInd)>nc || any(colInd<1 | colInd > nc ))
      stop("Colv dendrogram doesn't match size of x")
  }
  else if(identical(Colv, "Rowv")) {
    if(nr != nc)
      stop('Colv = "Rowv" but nrow(x) != ncol(x)')
    if(exists("ddr"))
    {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    } else
      colInd <- rowInd
  } else if(is.integer(Colv))
  {## Compute dendrogram and do reordering based on given vector
    hcc <- hclustfun(distfun(if(symm)x else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    
    colInd <- order.dendrogram(ddc)
    if(nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv))
  {## If TRUE, compute dendrogram and do reordering based on rowMeans
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if(symm)x else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)
    
    colInd <- order.dendrogram(ddc)
    if(nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else
  {
    colInd <- 1:nc
  }
  
  
#   o=matrix(0,nrow=classCount,ncol=2)
#   colnames(o)=c('class index', 'class size')
#   ordered=dataClasses[rev(rowInd)]
#   for(i in 1:classCount){
#     totalLength=length(ordered)
#     classIndex = ordered[1]
#     keep=ordered!=ordered[1]
#     ordered=ordered[keep]
#     classSize = totalLength - length(ordered)
#     o[i,] = c(classIndex, classSize)
#   }
  
  if(clusterSort){#sort within each cluster, max rowSums at top
    for(i in 1:classCount){
      keep = dataClasses[rowInd] == i 
      subset = x[rowInd,,drop = F][keep,,drop = F]
      in_cluster_o = order(rowSums(subset))
      x[rowInd,][keep,] = subset[in_cluster_o,]
    }
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  
  ## reorder x & cellnote
  x <- x[rowInd, colInd]
  
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  
  if(is.null(labRow))
    labRow <- if(is.null(rownames(x))) (1:nr)[rowInd] else rownames(x)
  else
    labRow <- labRow[rowInd]
  
  if(is.null(labCol))
    labCol <- if(is.null(colnames(x))) (1:nc)[colInd] else colnames(x)
  else
    labCol <- labCol[colInd]
  
  if(scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <-  sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if(scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <-  sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  
  ## Set up breaks and force values outside the range into the endmost bins
  if(missing(breaks) || is.null(breaks) || length(breaks)<1 )
  {
    if( missing(col) ||  is.function(col) )
      breaks <- 16
    else
      breaks <- length(col)+1
  }
  
  if(length(breaks)==1)
  {
    if(!symbreaks)
      breaks <- seq( min(minVal, na.rm=na.rm), max(maxVal,na.rm=na.rm), length=breaks)
    else
    {
      extreme <- max(abs(x), na.rm=TRUE)
      breaks <- seq( -extreme, extreme, length=breaks )
    }
  }
  
  nbr <- length(breaks)
  ncol <- length(breaks)-1
  
  if(class(col)=="function")
    col <- col(ncol)
  
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  
  x[x<min.breaks] <- min.breaks
  x[x>max.breaks] <- max.breaks
  
  
  ## Calculate the plot layout
  if( missing(lhei) || is.null(lhei) )
    lhei <- c(keysize, 4)
  
  if( missing(lwid) || is.null(lwid) )
    lwid <- c(keysize, 4)
  
  if( missing(lmat) || is.null(lmat) )
  {
    lmat <- rbind(4:3, 2:1)
    
    if(!missing(ColSideColors)) { ## add middle row to layout
      if(!is.character(ColSideColors) || length(ColSideColors) != nc)
        stop("'ColSideColors' must be a character vector of length ncol(x)")
      lmat <- rbind(lmat[1,]+1, c(NA,1), lmat[2,]+1)
      lhei <- c(lhei[1], 0.2, lhei[2])
    }
    
    if(!missing(RowSideColors)) { ## add middle column to layout
      if(!is.character(RowSideColors) || length(RowSideColors) != nr)
        stop("'RowSideColors' must be a character vector of length nrow(x)")
      lmat <- cbind(lmat[,1]+2, c(rep(NA, nrow(lmat)-1), 1), lmat[,2]+2,c(rep(NA, nrow(lmat)-1), 2))
      lwid <- c(lwid[1], 0.2, lwid[2], .2)
    }
    
    lmat[is.na(lmat)] <- 0
  }
  
  if(length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  
  if(length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  
  ## Graphics `output' -----------------------
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  #jrb modify layout to accomadate line plots and margins
  linePlotRegion = matrix(c(0,rep(x=1,classCount)+max(lmat),0,2:(classCount+1)+max(lmat)),nrow=classCount+1,ncol=2,byrow = FALSE)
  combinedRegions = matrix(0,ncol = (ncol(lmat) + ncol(linePlotRegion)), nrow = nrow(linePlotRegion))
  combinedRegions[1,1:ncol(lmat)] = lmat[1,]
  for(i in 2:nrow(combinedRegions))
    combinedRegions[i,1:ncol(lmat)] = lmat[2,]
  combinedRegions[,(ncol(lmat)+1):ncol(combinedRegions)] = linePlotRegion
  #print(combinedRegions)
  lwid=c(lwid,1,1)
  #print(lwid)
  lhei=c(1.5,rep(4/classCount,classCount))
  if(margins[1] > 0){
    combinedRegions=cbind(combinedRegions, rep(0, nrow(combinedRegions)))
    lwid=c(lwid, margins[1])
  }
  if(margins[2] > 0){
    combinedRegions=rbind(combinedRegions, rep(0, ncol(combinedRegions)))
    lhei=c(lhei, margins[2])
  }
  #print(lhei)
  nf=layout(combinedRegions, widths = lwid, heights = lhei, respect = TRUE)
  ## draw the side bars
  if(!missing(RowSideColors)) {
    par(mar = c(0,0, 0,0.5))
    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    par(mar = c(0,0.5, 0,0))
    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
  }
  if(!missing(ColSideColors)) {
    par(mar = c(0.5,0, 0,0))
    image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
  }
  ## draw the main carpet
  par(mar = c(0, 0, 0, 0))
  #if(scale != "none" || !symm)
  #  {
  x <- t(x)
  cellnote <- t(cellnote)
  #  }
  if(revC)
  { ## x columns reversed
    iy <- nr:1
    if(exists("ddr"))
      ddr <- rev(ddr)
    x <- x[,iy]
    cellnote <- cellnote[,iy]
  }
  else iy <- 1:nr
  
  ## display the main carpet
  image(1:nc, 1:nr, x, xlim = 0.5+ c(0, nc), ylim = 0.5+ c(0, nr),
        axes = FALSE, xlab = "", ylab = "", col=col, breaks=breaks,
        ...)
  retval$carpet <- x
  if(exists("ddr"))
    retval$rowDendrogram <- ddr
  if(exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  
  ## fill 'na' positions with na.color
  if(!invalid(na.color) & any(is.na(x)))
  {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col=na.color, add=TRUE)
  }
  
  ## add column labels
  if(is.null(srtCol))
    axis(1,
         1:nc,
         labels= labCol,
         las= 2,
         line= -0.5 + offsetCol,
         tick= 0,
         cex.axis= .45 * cexCol,
         hadj=adjCol[1],
         padj=adjCol[2]
    )
  else
  {
    if(is.numeric(srtCol))
    {
      if(missing(adjCol) || is.null(adjCol))
        adjCol=c(1,NA)
      
      xpd.orig <- par("xpd")
      par(xpd=NA)
      xpos <- axis(1, 1:nc, labels=rep("", nc), las=2, tick=0)
      text(x=xpos,
           y=par("usr")[3] - (1.0 + offsetCol) * strheight("M"),
           labels=labCol,
           ##pos=1,
           adj=adjCol,
           cex=cexCol * .45,
           srt=srtCol)
      par(xpd=xpd.orig)
    }
    else
      warning("Invalid value for srtCol ignored.")
  }
  
  #jrb i don't ever want row labels
  ## add row labels
  #   if(is.null(srtRow))
  #     {
  #       axis(4,
  #            iy,
  #            labels=labRow,
  #            las=2,
  #            line=-0.5+offsetRow,
  #            tick=0,
  #            cex.axis=cexRow,
  #            hadj=adjRow[1],
  #            padj=adjRow[2]
  #            )
  #     }
  #   else
  #     {
  #       if(is.numeric(srtRow))
  #         {
  #           xpd.orig <- par("xpd")
  #           par(xpd=NA)
  #           ypos <- axis(4, iy, labels=rep("", nr), las=2, line= -0.5, tick=0)
  #           text(x=par("usr")[2] + (1.0 + offsetRow) * strwidth("M"),
  #                y=ypos,
  #                labels=labRow,
  #                adj=adjRow,
  #                cex=cexRow,
  #                srt=srtRow
  #                )
  #           par(xpd=xpd.orig)
  #         }
  #       else
  #         warning("Invalid value for srtRow ignored.")
  #     }
  
  
  
  ## add row and column headings (xlab, ylab)
  if(!is.null(xlab)) mtext(xlab, side = 1, line = 0 - 1.25)
  if(!is.null(ylab)) mtext(ylab, side = 4, line = 0 - 1.25)
  
  ## perform user-specified function
  if (!missing(add.expr))
    eval(substitute(add.expr))
  
  ## add 'background' colored spaces to visually separate sections
  
  if((ncol(t(x)) %% nsplits) != 0){
    print('ncol not divisible by num of splits! ignoring splits.')
    nsplits = 1
  }
  win = ncol(t(x)) / nsplits
  if(nsplits > 1){
    colsep = 1:(nsplits-1)
    colsep = colsep * win
    for(csep in colsep)
      rect(xleft =csep+0.5,               ybottom=0,
           xright=csep+0.5+sepwidth[1],   ytop=ncol(x)+1,
           lty=1, lwd=1, col=sepcolor, border=sepcolor)
  }
  
  
  #jrb set rowsep to distinguish classes
  sepcolor='black'
  o=matrix(0,nrow=classCount,ncol=2)
  colnames(o)=c('class index', 'class size')
  ordered=dataClasses[rev(rowInd)]
  for(i in 1:classCount){
    totalLength=length(ordered)
    classIndex = ordered[1]
    keep=ordered!=ordered[1]
    ordered=ordered[keep]
    classSize = totalLength - length(ordered)
    o[i,] = c(classIndex, classSize)
  }
  rowsep=1:(classCount-1)
  for(i in 1:(classCount-1)){
    #calculate number of rows covered so far as fraction of total
    rowsep[i] = sum(o[1:i,2])
  }
  if(!missing(rowsep))
    for(rsep in rowsep)
      rect(xleft =0,          ybottom= (ncol(x)+1-rsep)-0.5,
           xright=nrow(x)+1,  ytop   = (ncol(x)+1-rsep)-0.5 - sepwidth[2],
           lty=1, lwd=1, col=sepcolor, border=sepcolor)
  
  
  ## show traces
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled  <- scale01(t(x), min.scale, max.scale)
  
  if(trace %in% c("both","column") )
  {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for( i in colInd )
    {
      if(!is.null(vline))
      {
        abline(v=i-0.5 + vline.vals, col=linecol, lty=2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[,i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv)-0.5
      lines(x=xv, y=yv, lwd=1, col=tracecol, type="s")
    }
  }
  
  
  if(trace %in% c("both","row") )
  {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for( i in rowInd )
    {
      if(!is.null(hline))
      {
        abline(h=i - 0.5 + hline.vals, col=linecol, lty=2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i,] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1-0.5
      lines(x=xv, y=yv, lwd=1, col=tracecol, type="s")
    }
  }
  
  
  
  if(!missing(cellnote))
    text(x=c(row(cellnote)),
         y=c(col(cellnote)),
         labels=c(cellnote),
         col=notecol,
         cex=notecex)
  
  ## the two dendrograms :
  par(mar = c(0, 0, 0, 0))
  if( dendrogram %in% c("both","row") )
  {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else
    plot.new()
  
  par(mar = c(0, 0, if(!is.null(main)) 5 else 0, 0))
  
  if( dendrogram %in% c("both","column") )
  {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else
    plot.new()
  
  ## title
  if(!is.null(main)) title(main, cex.main = 1.1*op[["cex.main"]])
  
  ## Add the color-key
  if(T)
  {
    mar <- c(3.5, 3.5, 5, 2)
    if (!is.null(key.xlab) && is.na(key.xlab))
      mar[1] <- 2
    #       if (!is.null(key.ylab) && is.na(key.ylab))
    #           mar[2] <- 2
    if (!is.null(key.title) && is.na(key.title))
      mar[3] <- 1
    par(mar = mar, cex=0.5, mgp=c(2.1, 1, 0))
    if (length(key.par) > 0)
      do.call(par, key.par)
    tmpbreaks <- breaks
    
    if(symkey)
    {
      max.raw <- max(abs(c(x,breaks)),na.rm=TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm=TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm=TRUE)
    }
    else
    {
      min.raw <- min(minVal, na.rm=TRUE) ## Again, modified to use scaled
      max.raw <- max(maxVal, na.rm=TRUE) ## or unscaled (SD 12/2/03)
    }
    
    z <- seq(min.raw, max.raw, length=length(col))
    image(z=matrix(z, ncol=1),
          col=col, breaks=tmpbreaks,
          xaxt="n", yaxt="n")
    
    par(usr=c(0,1,0,1))
    if (is.null(key.xtickfun)) {
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      xargs <- list(at=xv, labels=lv)
    } else {
      xargs <- key.xtickfun()
    }
    xargs$side <- 1
    do.call(axis, xargs)
    if (is.null(key.xlab)) {
      if(scale=="row")
        key.xlab <- "Row Z-Score"
      else if(scale=="column")
        key.xlab <- "Column Z-Score"
      else
        key.xlab <- "Value"
    }
    if (!is.na(key.xlab)) {
      mtext(side=1, key.xlab, line=par("mgp")[1], padj=0.5,cex = .8)
    }
    
    title("Color Key")
  }else
    plot.new()
  
  revx=t(x)
  revx=revx[rev(1:nrow(revx)),]
  #jrb draw line plots
  fetchClass=function(n){
    # data: original data
    # n: the class number from line plot(1 at top, last at bottom)
    # classDesc: [[1]] returned by makeFigure
    # classAssign: [[2]] returned by makeFigure
    # return list of values for all class members, values relative to column 1, and mean of each column's values relative to column 1
    start=1
    if(n>1)
      start = sum(o[1:(n-1),2])+1
    end=sum(o[1:n,2])
    c1 = as.matrix(revx[start:end,,drop = F])
    #print(c1)
    nc1 = c1[1,]
    if(nrow(c1) > 1)
      nc1 = t(apply(c1,MARGIN = 1,FUN = function(row) return(row/(maxVal-minVal))))
    #nc1=c1
    mc1 = c1
    if(nrow(c1) > 1) mc1 = colMeans(c1)
    #   plot(x = range(1:3),y = range(nc1), col = rgb(0,0,0,.3))
    #   for(i in 1:nrow(nc1)){
    #     lines(1:3,nc1[i,])
    #   }
    return(list(c1,nc1,mc1))
  }
  
  avgA = matrix(0,classCount, ncol(t(x)))
  for(i in 1:nrow(avgA)){
    avgA[i,] = fetchClass(i)[[3]]
  }
  
  #need to reverse order
  #o=o[classCount:1,]
  par(mai=rep(0,4))
  plot(x=c(0,1),frame.plot=FALSE, y=c(0,1), xaxt='n',yaxt='n', type="n", xlab="",
       ylab="",xlim=c(0,1),ylim=c(0,1))
  #no idea how these 'actual' margins are selected
  min=-.04
  max=1.04
  rng=max-min
  
  ### draw the dotted lines connecting line plots to heatmap
  par('mai'=rep(0,4))
  #line at top
  lines(c(min,max),c(max,max),lty=2)
  for(i in 1:classCount){
    #calculate number of rows covered so far as fraction of total
    heatmapFraction = sum(o[1:i,2]) / sum(o[,2])
    lineplotFraction = i/classCount
    #x is always from min to max
    #y goes from variable fraction on heatmap to constant fraction of plotting column
    lines(c(min,max),c(max - rng*heatmapFraction,max - rng*lineplotFraction),lty=2)
    
  }
  
  par(mar = c(0,0, 0,0.5))
  #draw line chart represented each class from clustering
  for (i in 1:classCount) { 
    xrange <- as.numeric(range(1:ncol(avgA)))
    yrange <- range(avgA)
    #print(avgA)
    
    
    # set up the plot 
    #op <- par(mar = rep(.01, 4))
    plot(xrange, yrange, xaxt='n',yaxt='n', type="n", xlab="",
         ylab="", ylim = c(minVal - .1*abs(minVal - maxVal),maxVal + .1*abs(minVal - maxVal))) 
    colors <- colorClasses[o[,1]] 
    linetype <- c(1:classCount) 
    plotchar <- seq(18,18+classCount,1)
    
    #   axis(side=1,tick=TRUE,at=days)
    tree <- avgA[i,]
    #print(tree)
    for(s in 1:nsplits){
      start = (s - 1) * win + 1
      end = s * win
      xs = start:end#first half of profile
      #print(xs)
      lines(xs, tree[xs], type="b", lwd=2.5,
            lty=1, col=colors[i], pch=plotchar[3]) 
    }
  } 
  
  
  ## Create a table showing how colors match to (transformed) data ranges
  retval$colorTable <- data.frame(
    low=retval$breaks[-length(retval$breaks)],
    high=retval$breaks[-1],
    color=retval$col
  )
  
  ## If user has provided an extra function, call it.
  if(!is.null(extrafun))
    extrafun()
  x=t(x)
  x=x[rev(1:nrow(x)),]
  rowInd=rev(rowInd)
  return(list(o[,2],rowInd,x,colors))
}


# $Id: textplot.R 1213 2007-11-01 20:20:10Z warnes $

textplot <- function(object, halign="center", valign="center", cex, ... )
  UseMethod('textplot')


textplot.default <- function(object,
                             halign=c("center","left","right"),
                             valign=c("center","top","bottom"),
                             cex, ... )
{
  
  if (is.matrix(object) || (is.vector(object) && length(object)>1) )
    return(textplot.matrix(object, halign, valign, cex, ... ))
  
  halign <- match.arg(halign)
  valign <- match.arg(valign)
  
  textplot.character(object, halign,  valign, cex, ...)
}


textplot.data.frame <- function(object,
                                halign=c("center","left","right"),
                                valign=c("center","top","bottom"),
                                cex, ... )
  textplot.matrix(object, halign, valign, cex, ... )


textplot.matrix <- function(object,
                            halign=c("center","left","right"),
                            valign=c("center","top","bottom"),
                            cex, cmar=2, rmar=0.5,
                            show.rownames=TRUE, show.colnames=TRUE,
                            hadj=1,
                            vadj=1,
                            mar= c(1,1,4,1)+0.1,
                            col.data=par("col"),
                            col.rownames=par("col"),
                            col.colnames=par("col"),
                            ... )
{
  
  if(is.vector(object))
    object <- t(as.matrix(object))
  else
    object <- as.matrix(object)
  
  # check dimensions of col.data, col.rownames, col.colnames
  if(length(col.data)==1)
    col.data <- matrix(col.data, nrow=nrow(object), ncol=ncol(object))
  else
    if( nrow(col.data)!=nrow(object) || ncol(col.data)!=ncol(object) )
      stop("Dimensions of 'col.data' do not match dimensions of 'object'.")
  
  if(length(col.rownames)==1)
    col.rownames <- rep(col.rownames, nrow(object))      
  
  if(length(col.colnames)==1)
    if(show.rownames)
      col.colnames <- rep(col.colnames, ncol(object)+1)
  else
    col.colnames <- rep(col.colnames, ncol(object))
  
  halign=match.arg(halign)
  valign=match.arg(valign)
  
  opar <- par()[c("mar","xpd","cex")]
  on.exit( par(opar) )
  par(mar=mar, xpd=FALSE )
  
  # setup plot area
  plot.new()
  plot.window(xlim=c(0,1),ylim=c(0,1), log = "", asp=NA)
  
  
  
  # add 'r-style' row and column labels if not present
  if( is.null(colnames(object) ) )
    colnames(object) <- paste( "[,", 1:ncol(object), "]", sep="" )
  if( is.null(rownames(object)) )
    rownames(object) <- paste( "[", 1:nrow(object), ",]", sep="")
  
  
  # extend the matrix to include row and column labels
  if( show.rownames )
  {
    object <- cbind( rownames(object), object )
    col.data <- cbind( col.rownames, col.data )
    
  }
  if( show.colnames )
  {
    object <- rbind( colnames(object), object )
    col.data <- rbind( col.colnames, col.data )
  }
  
  # set the character size
  if( missing(cex) )
  {
    cex <- 1.0
    lastloop <- FALSE
  }
  else
  {
    lastloop <- TRUE
  }
  
  for (i in 1:20)
  {
    oldcex <- cex
    
    width  <- sum(
      apply( object, 2,
             function(x) max(strwidth(x,cex=cex) ) )
    ) +
      strwidth('M', cex=cex) * cmar * ncol(object)
    
    height <- strheight('M', cex=cex) * nrow(object) * (1 + rmar)
    
    if(lastloop) break
    
    cex <- cex / max(width,height)
    
    if (abs(oldcex - cex) < 0.001)
    {
      lastloop <- TRUE
    }
  }
  
  
  # compute the individual row and column heights
  rowheight<-strheight("W",cex=cex) * (1 + rmar)
  colwidth<- apply( object, 2, function(XX) max(strwidth(XX, cex=cex)) ) +
    strwidth("W")*cmar
  
  
  width  <- sum(colwidth)
  height <- rowheight*nrow(object)
  
  # setup x alignment
  if(halign=="left")
    xpos <- 0
  else if(halign=="center")
    xpos <- 0 + (1-width)/2
  else #if(halign=="right")
    xpos <- 0 + (1-width)
  
  # setup y alignment
  if(valign=="top")
    ypos <- 1
  else if (valign=="center")
    ypos <- 1 - (1-height)/2
  else #if (valign=="bottom")
    ypos <- 0 + height
  #xpos = xpos - .1
  x <- xpos
  y <- ypos
  
  # iterate across elements, plotting them
  xpos<-x
  for(i in 1:ncol(object)) {
    xpos <- xpos + hadj*colwidth[i]
    for(j in 1:nrow(object)) {
      ypos<-y-(j-1)*rowheight
      if( (show.rownames && i==1) || (show.colnames && j==1) )
        text(xpos, ypos, object[j,i], adj=c(hadj,vadj), cex=cex, font=2,
             col=col.data[j,i], ... )
      else
        text(xpos, ypos, object[j,i], adj=c(hadj,vadj), cex=cex, font=1,
             col=col.data[j,i], ... )
    }
    xpos <- xpos + (1-hadj)*colwidth[i]
  }
  
  par(opar)
}

textplot.character <- function (object,
                                halign = c("center", "left", "right"),
                                valign = c("center", "top", "bottom"),
                                cex, fixed.width=TRUE,
                                cspace=1,
                                lspace=1,
                                mar=c(0,0,3,0)+0.1,
                                tab.width=8,
                                ...)
{
  object <- paste(object,collapse="\n",sep="")
  object <- replaceTabs(object, width=tab.width)
  
  halign = match.arg(halign)
  valign = match.arg(valign)
  plot.new()
  
  opar <- par()[c("mar","xpd","cex","family")]
  on.exit( par(opar) )
  
  par(mar=mar,xpd=FALSE )
  if(fixed.width)
    par(family="mono")
  
  plot.window(xlim = c(0, 1), ylim = c(0, 1), log = "", asp = NA)
  
  slist   <- unlist(lapply(object, function(x) strsplit(x,'\n')))
  slist   <- lapply(slist, function(x) unlist(strsplit(x,'')))
  
  slen    <- sapply(slist, length)
  slines  <- length(slist)
  
  if (missing(cex))
  {
    lastloop <- FALSE
    cex <- 1
  }
  else
    lastloop <- TRUE
  
  
  for (i in 1:20)
  {
    oldcex <- cex
    #cat("cex=",cex,"\n")
    #cat("i=",i,"\n")
    #cat("calculating width...")
    cwidth  <- max(sapply(unlist(slist), strwidth,  cex=cex)) * cspace
    #cat("done.\n")
    #cat("calculating height...")
    cheight <- max(sapply(unlist(slist), strheight, cex=cex)) * ( lspace + 0.5 )
    #cat("done.\n")
    
    width <- strwidth(object, cex=cex)
    height <- strheight(object, cex=cex)
    
    if(lastloop) break
    
    cex <- cex  / max(width, height)
    
    if (abs(oldcex - cex) < 0.001)
    {
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
  
  text(x=xpos, y=ypos, labels=object, adj=c(0,1),
       cex=cex, ...)
  
  par(opar)
  invisible(cex)
}

write.HeatmapLists = function(heatmap3_output, name = "heatmap_lists", output_directory = ""){
  res = heatmap3_output
  print(output_directory)
  if(!file.exists(output_directory))
    dir.create(output_directory)
  root = paste(output_directory, name, sep = "/")
  print(root)
  classSizes = res[[1]]
  asPlotted = res[[3]]
  classColors = res[[4]]
  data_out = matrix(0, nrow = sum(classSizes), ncol = 3)
  data_out[,1:2] = as.matrix(indexDict[rownames(asPlotted),c(1,4)])
  colors_out = matrix('black', nrow = sum(classSizes), ncol = 3)
  soFar = 0
  for(j in 1:length(classSizes)){
    size = classSizes[j]
    start = soFar + 1
    end = soFar + size
    soFar = soFar + size
    data_out[start:end,3] = paste('Class',j)
    colors_out[start:end,1:3] = classColors[j]
  }
  pdf(paste(root,'.pdf', sep = ''))#, width = 8.5, height = .2 * nrow(data_out))
  textplot(data_out[,1:2], col.data = colors_out[,1:2], show.rownames = F, show.colnames = F, halign = "left", hadj = 0)
  dev.off()
  
  png(paste(root,'.png', sep = ''))
  textplot(as.matrix(data_out[,1]), col.data = as.matrix(colors_out[,3]), show.rownames = F, show.colnames = F)
  dev.off()
  
  write.table(x = sub(data_out,pattern = ",",replacement = " "),file = paste(root,'.csv', sep = ''), quote = F, row.names = F, col.names = F, sep = ',')
  
  rownames(asPlotted) = indexDict[rownames(asPlotted),1]
  asPlotted = cbind(rownames(asPlotted),asPlotted)
  colnames(asPlotted)[1] = 'Symbol'
  asPlotted = unlist(strsplit(asPlotted[,1], split = ','))
  write.table(x = asPlotted,file = paste(root,'_ipa.txt', sep = ''), col.names = F, row.names = F, quote = F, sep = '\t')
  
}

plot.HeatmapLists = function(heatmap3_output){
  res = heatmap3_output
  classSizes = res[[1]]
  asPlotted = res[[3]]
  classColors = res[[4]]
  data_out = matrix(rownames(asPlotted), ncol = 1)
  colors_out = matrix('black', nrow = sum(classSizes), ncol = 1)
  soFar = 0
  for(j in 1:length(classSizes)){
    size = classSizes[j]
    start = soFar + 1
    end = soFar + size
    soFar = soFar + size
    colors_out[start:end,] = classColors[j]
  }
  a = matrix("black", nrow = 100, ncol = 10)
  b = matrix("", nrow = 100, ncol = 10)
  a[1:nrow(colors_out)] = colors_out[1:nrow(colors_out)]
  b[1:nrow(data_out)] = data_out[1:nrow(data_out)]
  textplot(b, col.data = a, show.rownames = F, show.colnames = F, halign = "left", hadj = 0)
}