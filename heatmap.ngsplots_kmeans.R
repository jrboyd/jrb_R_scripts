#plots K27me3 and K4me3 combinations of heatmaps
#clustering is kmeans then sorted within clusters

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
                        ylab = NULL, lmat = NULL, lhei = 1, lwid = NULL, extrafun = NULL, rowsep.major.labels = NULL, left_mai = .7, cex.main = 2,
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
  label.height = .5
  main.height = 1.2
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
  labRowSize = 0
  if (!is.null(rowsep.major.labels)) {
    labRowSize = 1
    lmat <- cbind(lmat, c(rep(0, body.i-1), max(lmat)+1, rep(0, nrow(lmat)-body.i)))
    lwid <- c(dev.width - RowSideColors.size - labRowSize, labRowSize)
  }
  if(left_mai > 0){
    lmat <- cbind(rep(0, nrow(lmat)), lmat)
    lwid <- c(left_mai, dev.width - RowSideColors.size - labRowSize - left_mai)
  }
  
  
  #   print(lmat)
  #   print(lhei)
  #   print(lwid)
  
  
  
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
    text(.5,.5, main, cex = cex.main)
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
      mtext(side = 3, key.lab, line = 1, padj = 1, 
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

#source('jrb_R_scripts//heatmap.3-split.R')
#source('heatmap.2.2.R')
library('reshape')
library('ggplot2')

heatmap.ngsplots = function(ngs_profiles, main_title = NULL, profiles_to_plot, nclust = 6, labels_below = NA, labels_above = NA, fg_toPlot = character(), fg_label = NULL, sortClustersByTotal = F, hmapColors = c('royalblue1', 'black', 'yellow'), ...){
  if(length(labels_below) < 2 && is.na(labels_below)){
    labels_below = profiles_to_plot
  }
  if(length(profiles_to_plot) %% length(labels_below) != 0){
    stop('length of labels_below must divide profiles_to_plot evenly!')
  }
  if(length(profiles_to_plot) %% length(labels_above) != 0){
    stop('length of labels_above must divide profiles_to_plot evenly!')
  }
  
  
  prof = matrix(0, ncol = 0, nrow = nrow(ngs_profiles[[1]]))
  hidden = lapply(profiles_to_plot, function(x){
    prof <<- cbind(prof, ngs_profiles[[x]])
  })
  
  #HARDCODED DEBUG STUFF
  #   prof = prof[1:100,]
  #   fg_toPlot = rownames(prof)[1]
  
  
  cr = colorRamp(hmapColors)
  tmp = 0:100 / 100
  colors = rgb(cr(tmp)/255)
  
  len = length(profiles_to_plot)
  
  cseps = 1:(len-1) * ncol(prof) / len
  
  labColAbove = rep('', ncol(prof))
  if(!is.na(labels_above[1])){
    labels_loc = round((1:length(labels_above) -.499) * ncol(prof) / length(labels_above))
    labColAbove[labels_loc] = labels_above
  }
  
  
  labels = rep('', ncol(prof))
  labels_loc = round((1:length(labels_below)-.499) * ncol(prof) / length(labels_below))
  labels[labels_loc] = labels_below
  set.seed(1)
  kmclust = kmeans(prof, centers = nclust, iter.max = 10)
  o = order(kmclust$cluster)#sort prof to be in cluster number order
  prof = prof[o,, drop = F]
  kmclust$cluster = kmclust$cluster[o]#keep cluster assignment sorted in same order as prof
  
  for(i in 1:nclust){#sort within each cluster by rowSum
    size = kmclust$size
    start = 1
    if(i > 1){
      start = sum(kmclust$size[1:(i-1)]) + 1
    }
    end = sum(kmclust$size[1:(i)])
    #print(paste(start, end))
    o = order(rowSums(prof[start:end,,drop = F]), decreasing = T)
    prof[start:end,] = prof[start:end,,drop = F][o,, drop = F]
  }
  
  prof_ordered = matrix(0, ncol = ncol(prof), nrow = 0)
  
  
  
  get_kmclust = function(i){
    start = 1
    if(i > 1){
      start = sum(kmclust$size[1:(i-1)]) + 1
    }
    end = sum(kmclust$size[1:(i)])
    return(prof[start:end,])
  }
  if(sortClustersByTotal){#sort clusters by total of their centers
    kmclust_order = order(apply(kmclust$centers, 1, sum), decreasing = T)
    for(i in kmclust_order){#sort k means clusters by their center hierarchical clutering
      new_cluster = get_kmclust(i)
      start = nrow(prof_ordered) + 1
      prof_ordered = rbind(prof_ordered, new_cluster)
      end = nrow(prof_ordered)
      #kmclust$cluster[start:end] = i
    }
  }else{
    hiclust = hclust(dist(kmclust$centers))
    kmclust_order = hiclust$order
    for(i in kmclust_order){#sort k means clusters by their center hierarchical clutering
      new_cluster = get_kmclust(i)
      start = nrow(prof_ordered) + 1
      prof_ordered = rbind(prof_ordered, new_cluster)
      end = nrow(prof_ordered)
      #kmclust$cluster[start:end] = i
    }
  }
  
  
  rseps = 1:(nclust-1)#cacluated row seperations
  kmclust$size = kmclust$size[kmclust_order]
  
  for(i in 1:nclust){
    size = kmclust$size
    start = 1
    if(i > 1){
      start = sum(kmclust$size[1:(i-1)]) + 1
    }
    end = sum(kmclust$size[1:(i)])
    if(i < nclust){
      rseps[i] = end
    }
  }
  
  rColorChoices = RColorBrewer::brewer.pal(nclust, 'Dark2')
  rColors = rColorChoices[kmclust$cluster[rownames(prof_ordered)]]
  names(rColors) = rownames(prof_ordered)
  #fg_toPlot = intersect(fg_toPlot, rownames(prof_ordered))
  if(length(fg_toPlot) > 0){#extract fg_plot as special cluster
    rseps = c(0, rseps)
    rColors[fg_toPlot] = 'white'
    #move each fg_toPlot to top of profile, reduce rseps accordingly
    fg_indexes = 1:nrow(prof_ordered)  
    names(fg_indexes) = rownames(prof_ordered)
    keep = rep(F, length(fg_indexes))
    names(keep) = names(fg_indexes)
    keep[fg_toPlot] = T
    fg_indexes = fg_indexes[keep]
    #fg_indexes = fg_indexes[fg_toPlot]
    fg_indexes = sort(fg_indexes)#sort is necessary, as long as indexes don't change order relative to one another, they can be moved independently
    #insert_at = 1
    
    for(i in fg_indexes){
      #print(i)
      affected = rseps < i
      rseps[affected] = rseps[affected] + 1
      a = integer()
      b = i
      c = integer()
      if(i > 1) a = 1:(i-1)
      if(i < nrow(prof_ordered)) c = (i+1):nrow(prof_ordered)
      new_order = c(b, a, c)
      prof_ordered = prof_ordered[new_order ,]
      #rColors = rColors[new_order]
      
    }
    
    #sort the new cluster
    o = order(rowSums(prof_ordered[1:(rseps[1]), ]), decreasing = T)
    prof_ordered[1:(rseps[1]), ] = prof_ordered[1:(rseps[1]), ][o,]
    
    #may be sufficient to reverse new cluster
    #prof_ordered[1:(rseps[1]), ] = prof_ordered[(rseps[1]):1, ]
  }
  
  
  hidden = heatmap.2.2(prof_ordered, RowSideColors = rColors[rownames(prof_ordered)], Rowv = F, Colv = F, dendrogram = 'n', trace = 'n', scale = 'n', labRow = '', labCol = labels, srtCol = 0, adjCol = c(.5,1), cexCol = 1, col = colors, symm = T, rowsep.major = rseps, colsep.minor = cseps, sepwidth.minor = 2, sepwidth.major = max(nrow(prof_ordered)/100,1), density.info = 'none', labColAbove = labColAbove, na.color = 'red', rowsep.major.labels = fg_label, main = main_title, ...)
  cluster_members = list()
  #print(nrow(prof_ordered))
  for(i in 0:length(rseps)){
    start = 1
    if(i > 0) start = rseps[i] + 1
    end = nrow(prof_ordered)
    if(i < length(rseps)) end = rseps[i + 1]
    
    #print(start)
    #print(end)
    cluster_members[[i+1]] = rownames(prof_ordered)[start:end]
  }
  ngs_res = list(cluster_members = cluster_members, colors = unique(rColors[rownames(prof_ordered)]), as_plotted = prof_ordered)
  return(ngs_res)
}

cl = 'MCF10A'
hm = c('H3K27AC','H3K27ME3','H3K4AC') #, 'H3K4ME3', 'H4K20ME3')
toPlot = paste(cl, hm, sep = '_')
ngs_profiles = list()
nr = 10
nc = 49
for(tp in toPlot){
  dat = matrix(runif(nr * nc, 0, 100), nrow = nr, ncol = nc)
  rownames(dat) = 1:nrow(dat)
  reduce = runif(nr/2, 1, nr)
  dat[reduce,] = dat[reduce,] / 5
  ngs_profiles[[tp]] = dat
}



#png(paste0('test8sidecol2_', cl, '.png'), width = 4, height = 8, units = 'in', res = 450)
nclust = 3
pdf(file = NULL)
res = heatmap.ngsplots(ngs_profiles,# hmapColors = c('#2166ac', 'black', '#b2182b'),
                       main_title = 'title', toPlot, nclust = nclust, labels_below = hm, labels_above = cl, 
                       fg_toPlot = c(9,2), fg_label = c('selected', 1:nclust), sortClustersByTotal = T)
dev.off()
res = heatmap.ngsplots(ngs_profiles,# hmapColors = c('#2166ac', 'black', '#b2182b'),
                       main_title = 'title', toPlot, nclust = nclust, labels_below = hm, labels_above = cl, 
                       fg_toPlot = c(9,2), fg_label = paste(sapply(res[[1]], length), 'genes'), sortClustersByTotal = T)


#dev.off()


