heatmap.3 <-
function (x,
                       #returns list of classsizes, row order vector for input x, the data as plotted in heatmap, all with row 1 at top of heatmap and last row at bottom of heatmap
                       ## set number of splits
                       nsplits = 1,
                       ## dendrogram control
                       override_o = NA,
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
      
#       warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
#               dendrogram, "'. Omitting row dendogram.")
      
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
  if(length(override_o) == 1 && is.na(override_o)){
    #jrb row side colors is default and get data classes
    dataClasses = cutree(tree=hcr, k=classCount)
    orderOfClasses = unique(dataClasses[rowInd])
    tmp = classCount:1
    tmp = tmp[order(orderOfClasses)]
    require('RColorBrewer')
    colorClasses = brewer.pal(n = classCount,name = "Dark2")
    colorClasses = colorClasses[tmp]
    RowSideColors = colorClasses[dataClasses]
  }else{
    
    dataClasses = integer()
    hidden = sapply(1:nrow(override_o), function(x){
      dataClasses <<- c(dataClasses, rep(override_o[x,1], override_o[x,2]))
    })
    
    classCount = nrow(override_o)
    colorClasses = brewer.pal(n = classCount,name = "Dark2")
    RowSideColors = colorClasses[dataClasses]
  }
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
  lwid=c(lwid,1,1)
  lhei=c(1.5,rep(4/classCount,classCount))
  if(margins[1] > 0){
    combinedRegions=cbind(combinedRegions, rep(0, nrow(combinedRegions)))
    lwid=c(lwid, margins[1])
  }
  if(margins[2] > 0){
    combinedRegions=rbind(combinedRegions, rep(0, ncol(combinedRegions)))
    lhei=c(lhei, margins[2])
  }
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
  if(length(override_o) == 1 && !is.na(override_o)){
    colnames(override_o) = colnames(o)
    o = override_o
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
  print(o)
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
    yrange <- range(avgA[i,])
    xspace = (xrange[2] - xrange[1]) * .15
    yspace = (yrange[2] - yrange[1]) * .15
    # set up the plot 
    #op <- par(mar = rep(.01, 4))
    plot(xrange, yrange, xaxt='n',yaxt='n', type="n", xlab="",
         ylab="", ylim = c(min(yrange)-yspace, max(yrange)+yspace), xlim = c(min(xrange)-xspace, max(xrange)+xspace))  #c(minVal - .1*abs(minVal - maxVal),maxVal + .1*abs(minVal - maxVal))) 
    colors <- colorClasses[o[,1]] 
    linetype <- c(1:classCount) 
    plotchar <- seq(18,18+classCount,1)
    
    #   axis(side=1,tick=TRUE,at=days)
    vals <- avgA[i,]
    for(s in 1:nsplits){
      start = (s - 1) * win + 1
      end = s * win
      xs = start:end#first half of profile
      lines(xs, vals[xs], type="b", lwd=2.5,
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
