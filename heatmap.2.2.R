#dummy out tracing, dummy out clustering
heatmap.2.2 = function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
                        distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                                           "row", "column", "none"), reorderfun = function(d, w) reorder(d, 
                                                                                                                                         w), symm = FALSE, scale = c("none", "row", "column"), 
                        na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
                        symbreaks = any(x < 0, na.rm = TRUE) || scale != "none", 
                        col = "heat.colors", colsep.minor = -1, colsep.major = -1, rowsep.minor = -1, rowsep.major = -1, sepcolor = "white", 
                        sepwidth.minor = 1, sepwidth.major = 3, cellnote, notecex = 1, notecol = "cyan", 
                        na.color = par("bg"), trace = c("column", "row", "both", 
                                                        "none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks), 
                        linecol = tracecol, mar.top = 2, mar.left = 2, mar.bot = 5, mar.right = 5, ColSideColors, RowSideColors = NULL, 
                        cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL,
                        labCol = NULL, labColAbove = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, 
                                                                                                    NA), adjCol = c(.5, .5), offsetRow = 0.5, offsetCol = 0.5, 
                        colRow = NULL, colCol = NULL, key = TRUE, keysize = .1, 
                        density.info = c("histogram", "density", "none"), denscol = tracecol, 
                        symkey = any(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, 
                        key.title = NULL, key.lab = NA, key.ylab = NULL, key.xtickfun = NULL, 
                        key.ytickfun = NULL, key.par = list(), main = NULL, xlab = NULL, 
                        ylab = NULL, lmat = NULL, lhei = 1, lwid = NULL, extrafun = NULL, rowsep.major.labels = NULL,
                        ...) 
{
  x.original = x
  rowsep.minor.original = rowsep.minor
  rowsep.major.original = rowsep.major
  
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  
  plot0 = function(width = 1, height = 1){
    fudge = 0.037037037037
    plot(c(0+fudge*width, width-fudge * width), c(0+fudge*height, height-fudge * height), type = 'n', xlab = '', ylab = '', axes = F, )
  }
  
  dev.width = par('din')[1]
  dev.height = par('din')[2]
  
  lmat = 1
  body.height = dev.height
  body.i = 1
  lwid = 1
  label.height = 1
  main.height = 1.5
  lhei = body.height
  
  if(!is.null(labColAbove)){
    lmat = rbind(max(lmat) + 1, lmat)
    lhei[body.i] = lhei[body.i] - label.height
    lhei = c(label.height, lhei)
    body.i = body.i + 1
  }
  if(!is.null(main)){
    lmat = rbind(max(lmat) + 1, lmat)
    lhei[body.i] = lhei[body.i] - main.height
    lhei = c(main.height, lhei)
    body.i = body.i + 1
  }
  if(!is.null(labCol)){
    lmat = rbind(lmat, max(lmat) + 1)
    lhei[body.i] = lhei[body.i] - label.height
    lhei = c(lhei, label.height)
  }
  if(key){
    key.height = 1
    lmat = rbind(lmat, max(lmat) + 1)
    lhei[body.i] = lhei[body.i] - key.height
    lhei = c(lhei, key.height)
  }
  RowSideColors.size = 0
  if (!is.null(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) != nrow(x)) 
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    RowSideColors.size = .8
    lmat <- cbind(lmat, c(rep(0, body.i-1), max(lmat)+1, rep(0, nrow(lmat)-body.i)))
    lwid <- c(dev.width - RowSideColors.size, RowSideColors.size)
  }
  if (!is.null(rowsep.major.labels)) {
    labRowSize = 1
    lmat <- cbind(lmat, c(rep(0, body.i-1), max(lmat)+1, rep(0, nrow(lmat)-body.i)))
    lwid <- c(dev.width - RowSideColors.size - labRowSize, labRowSize)
  }
  
  
  print(lmat)
  #print(lhei)
  
  
  
  layout(lmat, heights = lhei, widths = lwid)
  
  
  
  x = x[nrow(x):1,]#reverse so heatmap top row is top row of input matrix
  RowSideColors = rev(RowSideColors)
  nr <- nrow(x)
  nc <- ncol(x)
  
  
  insert_col = function(n, i){#insert n column following index i
    x <<- cbind(x[,1:i], matrix(NA, nrow = nrow(x), ncol = n),  x[,(i+1):ncol(x)])
    affected = colsep.minor > i
    colsep.minor[affected] <<- colsep.minor[affected] + n
    affected = colsep.major > i
    colsep.major[affected] <<- colsep.major[affected] + n
    if(!is.null(labCol)) labCol <<- c(labCol[1:i], rep('', n), labCol[(i + 1):length(labCol)])
    if(!is.null(labColAbove)) labColAbove <<- c(labColAbove[1:i], rep('', n), labColAbove[(i + 1):length(labColAbove)])
  }
  
  insert_row = function(n, i){#insert n rows following index i
    x <<- rbind(x[1:i,], matrix(NA, ncol = ncol(x), nrow = n),  x[(i+1):nrow(x),])
    affected = rowsep.minor > i
    rowsep.minor[affected] <<- rowsep.minor[affected] + n
    affected = rowsep.major > i
    rowsep.major[affected] <<- rowsep.major[affected] + n
    RowSideColors <<- c(RowSideColors[1:i], rep(NA, n), RowSideColors[(i+1):length(RowSideColors)])
  }
  
  #print(dim(x))
  do_col_major = colsep.major[1] > -1
  if (colsep.minor[1] > -1){ #insert  columns in x and color
    for (i in 1:length(colsep.minor)){ 
      csep = colsep.minor[i]
      insert_col(sepwidth.minor, csep)
    }
  }
  if (do_col_major){ #insert  columns in x and color
    for (i in 1:length(colsep.major)){ 
      csep = colsep.major[i]
      insert_col(sepwidth.major, csep)
    }
  }
  #print(dim(x))
  x = x[nrow(x):1,]#reverse rows temporarily
  RowSideColors = rev(RowSideColors)
  do_row_major = rowsep.major[1] > -1
  if (rowsep.minor[1] > -1){
    for (i in 1:length(rowsep.minor)){ 
      rsep = rowsep.minor[i]
      insert_row(sepwidth.minor, rsep)
    }
  }
  if (do_row_major){
    for (i in 1:length(rowsep.major)){ 
      rsep = rowsep.major[i]
      insert_row(sepwidth.major, rsep)
    }
  }
  x = x[nrow(x):1,]#unreverse rows
  RowSideColors = rev(RowSideColors)
  #print(dim(x))
  nc = ncol(x)
  nr = nrow(x)
  par(mai = rep(0,4))
  #plots the body of the heatap
  image(0:nc, 0:nr, t(x), xlim = c(0, nc), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
  
  #print(dim(x))
  #plots column labels above column
  apply_col_labels = function(col_labels, yshift = .5, yadj = .5){
    plot0(ncol(x))
    if (is.null(srtCol) && is.null(colCol)) 
      axis(1, 1:nc-.5, labels = col_labels, las = 2, line = -0.5 + 
             offsetCol, tick = 0, cex.axis = cexCol, padj = 1)
    else {
      if (is.null(srtCol) || is.numeric(srtCol)) {
        xpd.orig <- par("xpd")
        par(xpd = NA)
        xpos <- 1:nc-.5
        text(x = xpos, y = yshift, labels = col_labels, adj = c(.5, yadj),
             cex = cexCol, srt = srtCol, col = colCol)
        #print(colCol)
        par(xpd = xpd.orig)
      }
      else warning("Invalid value for srtCol ignored.")
    }
  }
  if(!is.null(labColAbove)){
    apply_col_labels(labColAbove, .1, 0)
  }
  
  #draw title
  if(!is.null(main)){
    plot.new()
    text(.5,.5, main, cex = 4)
  }
  
  #draws column labels below column
  if(!is.null(labCol)){
    apply_col_labels(labCol, .9, 1)
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
    if (missing(col) || is.function(col)) 
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                    length = breaks)
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
  if (key) {
    mai <- c(.5, .2, .2, .2)
    par(mai = mai, cex = 0.75, mgp = c(2, 1, 0))
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
      min.raw <- min.breaks
      max.raw <- max.breaks
    }
    z <- seq(min.raw, max.raw, by = min(diff(breaks)/4))
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
    if (!is.na(key.lab)) {
      mtext(side = 1, key.lab, line = par("mgp")[1], padj = 0.5, 
            cex = par("cex") * par("cex.lab"))
    }
  }
  if(!is.null(RowSideColors)){
    #print(length(RowSideColors))
    #print(nr)
    tmp = as.factor(RowSideColors)
    vals = matrix(as.numeric(tmp), ncol = 1)
    lev = levels(tmp)
    #print(dim(vals))
    #print(length(0:nr))
    par(mai = c(0,.1,0,.1))
    image(0:1, 0:nr, t(vals), xlim = c(0, 1), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "", col = lev)
  }
  if(!is.null(rowsep.major.labels)){
    par(mai = rep(0,4))
    plot0()
    for(i in 1:length(rowsep.major.labels)){
      rsep_prev = 0
      if(i > 1){
        rsep_prev = rowsep.major[i - 1] + sepwidth.major
      }
      rsep_curr = rowsep.major.original[i]
      if(i > 1){
        rsep_curr = rsep_prev + rowsep.major.original[i] - rowsep.major.original[i - 1]
      }
      if(i > length(rowsep.major.original)){
        rsep_curr = nrow(x)
      }
#     else{
#         rsep_curr = nrow(x)
#       }
      
      #print(nrow(x))
      rsep_mid = (mean(c(rsep_prev, rsep_curr)))
      
      ypos = 1 - (rsep_mid / nrow(x))
      text(.5, ypos, rowsep.major.labels[i], adj = c(.5,.5) )
#       text(.5, -.5/nrow(x), 'test', adj = c(.5,.5) )
#       text(.5, .5/nrow(x), 'test', adj = c(.5,.5) )
#       text(.5, 1.5/nrow(x), 'test', adj = c(.5,.5) )
#       text(.5, 1-2.5/nrow(x), 'test', adj = c(.5,.5) )
    }
    
  }
  
  #print(par('usr'))
}