#plots K27me3 and K4me3 combinations of heatmaps
#clustering is kmeans then sorted within clusters



#source('jrb_R_scripts//heatmap.3-split.R')
#source('heatmap.2.2.R')
library('png')
library('reshape')
library('ggplot2')

plot0 = function(width = 1, height = 1){
  fudge = 0.037037037037
  plot(c(0+fudge*width, width-fudge * width), c(0+fudge*height, height-fudge * height), type = 'n', xlab = '', ylab = '', axes = F, )
}

startRasterMode = function(width = 1, height = 1){
  png('tmp.png', units = 'in', width = width, height = height, res = 600)
  par(mai = rep(0,4))
}

stopRasterMode = function(mai = NA){
  dev.off()
  if(length(mai) > 1 && !is.na(mai)){
    par(mai = mai)
  }
  plot0()
  rasterImage(readPNG("tmp.png", native = FALSE), xleft = 0, xright = 1, ytop = 1, ybottom = 0,
              interpolate = FALSE)
}

heatmap.ngsplots = function(ngs_profiles, 
                            main_title = NULL, 
                            profiles_to_plot = NA, 
                            nclust = 6, 
                            labels_below = NA, 
                            labels_above = NA, 
                            fg_toPlot = character(), 
                            labels_right = NA, 
                            sortClustersByTotal = F, 
                            hmapColors = c('royalblue1', 'royalblue1', 'black', 'yellow'), 
                            labelWithCounts = F, 
                            fg_label = NA, 
                            label_clusters = T,
                            key.lab = 'log2 FE',
                            key.cex = 1.3,
                            cex.main = 2.5,
                            cex.row = 1.3,
                            cex.col = 2){
  if(length(profiles_to_plot) == 1 && is.na(profiles_to_plot)){
    profiles_to_plot = names(ngs_profiles)
  }
  if(length(labels_below) < 2 && is.na(labels_below)){
    labels_below = profiles_to_plot
  }
  if(length(profiles_to_plot) %% length(labels_below) != 0){
    stop('length of labels_below must divide profiles_to_plot evenly!')
  }
  if(length(profiles_to_plot) %% length(labels_above) != 0){
    stop('length of labels_above must divide profiles_to_plot evenly!')
  }
  
  
  prof = matrix(0, ncol = 0, nrow = nrow(ngs_profiles[[1]]))#assemble single matrix by joining selected profiles (matrices)
  hidden = lapply(profiles_to_plot, function(x){
    prof <<- cbind(prof, ngs_profiles[[x]])
  })
  
  cr = colorRamp(hmapColors)
  tmp = 0:100 / 100
  colors = rgb(cr(tmp)/255)
  
  len = length(labels_below)
  len_major = length(labels_above)
  cseps = 1:(len-1) * ncol(prof) / len #column separators are between joined profiles
  if(len == 1) cseps = -1
  cseps_major = 1:(len_major-1) * ncol(prof) / len_major
  if(len_major == 1) cseps_major = -1
  
  labels.above = rep('', ncol(prof))#column labels are centered above each profile
  if(!is.na(labels_above[1])){
    labels_loc = round((1:length(labels_above) -.499) * ncol(prof) / length(labels_above))
    labels.above[labels_loc] = labels_above
  }
  
  
  labels.below = rep('', ncol(prof))
  labels_loc = round((1:length(labels_below)-.499) * ncol(prof) / length(labels_below))
  labels.below[labels_loc] = labels_below
  
  set.seed(1)
  kmclust = kmeans(prof, centers = nclust, iter.max = 10)
  
  
  
  if(length(fg_toPlot) > 0){#extract fg_plot as special cluster
    #rseps = c(0, rseps)#add an
    kmclust$cluster[fg_toPlot] = 0#change cluster of fg
    kmclust$cluster = kmclust$cluster + 1
    
    
    nclust = nclust + 1
    kmclust$size = sapply(1:nclust, function(x){
      return(sum(kmclust$cluster == x))
    })
    kmclust$centers = rbind(colMeans(prof[fg_toPlot,]), kmclust$centers)
    rownames(kmclust$centers) = 1:nclust
  }
  
  o = order(kmclust$cluster)#sort prof to be in cluster number order
  prof = prof[o,, drop = F]
  kmclust$cluster = kmclust$cluster[o]#keep cluster assignment sorted in same order as prof
  
  for(i in 1:nclust){#sort within each cluster by rowSum
    #size = kmclust$size
    start = 1
    if(i > 1){
      start = sum(kmclust$size[1:(i-1)]) + 1
    }
    end = sum(kmclust$size[1:(i)])
    #print(paste(start, end))
    o = order(rowSums(prof[start:end,,drop = F]), decreasing = T)
    prof[start:end,] = prof[start:end,,drop = F][o,, drop = F]
  }
  
  prof_ordered = matrix(0, ncol = ncol(prof), nrow = 0)#prof_ordered contains profiles rearranged at cluster level
  get_kmclust = function(i){
    start = 1
    if(i > 1){
      start = sum(kmclust$size[1:(i-1)]) + 1
    }
    end = sum(kmclust$size[1:(i)])
    return(prof[start:end,,drop = F])
  }
  
  kmclust_order = 1:nclust
  if(length(fg_toPlot) > 0){
    sortByDist = F
    if(sortByDist){
      km_dist = as.matrix(dist(kmclust$centers))
      o = order(km_dist[,1])
      kmclust_order = o
    }else{
      kmclust_order = order(apply(kmclust$centers, 1, sum), decreasing = T)
      not_1 = kmclust_order != 1#move 1 to beginning
      kmclust_order = c(1, kmclust_order[not_1])
    }
    
    #    
  }else{
    if(sortClustersByTotal){#sort clusters by total of their centers
      kmclust_order = order(apply(kmclust$centers, 1, sum), decreasing = T)
    }else{
      hiclust = hclust(dist(kmclust$centers))#use hierarchical clustering on centers
      kmclust_order = hiclust$order
    }
  }
  
  rColorChoices = RColorBrewer::brewer.pal(nclust, 'Dark2')#set cluster id colors
  if(length(fg_toPlot) > 0){#correct for special selection color
    rColorChoices = c('white', rColorChoices[2:length(rColorChoices)-1])
  }
  rColors = character()
    
  for(i in 1:length(kmclust_order)){#sort k means clusters 
    new_cluster = get_kmclust(kmclust_order[i])
    prof_ordered = rbind(prof_ordered, new_cluster)
    
    cluster_color = rColorChoices[i]
    rColors = c(rColors, rep(cluster_color, nrow(new_cluster)))
  }
  names(rColors) = rownames(prof_ordered)
  
  
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
  
  if(labelWithCounts){
    if(is.null(labels_right)){
      labels_right = 1:nclust
    } else {
      labels_right = paste(kmclust$size, labels_right)
    }
  }else{
    if(length(labels_right) == 1){
      labels_right = rep(labels_right, nclust)
    }else if(length(labels_right) != nclust){
      stop('length of labels_right must be 1 or = to # of clusters (+1 if fg_toPlot is specified)')
    }
  }
  if(length(fg_toPlot) == 0 && !is.na(fg_label)){
    stop('cannot add fg_label, fg_toPlot is empty')
  }
  RowSideLabels = character()
  if(label_clusters){
    RowSideLabels = 1:nclust#paste('cluster', 1:nclust)
    if(length(fg_toPlot > 0)){
      if(!is.na(fg_label)){
        RowSideLabels = c(fg_label, 1:(nclust-1))# paste("cluster", 1:(nclust-1)))#special label for fg_toPlot
      }
    }
  }else if(!is.na(fg_label)){
    RowSideLabels = c(fg_label, rep('', nclust-1))#special label for fg_toPlot
  }
  extra_spacer = 0
  if(length(fg_toPlot) > 0){
    extra_spacer = 2
  }
  args = list(
    x = prof_ordered,key.lab = key.lab, key.cex = key.cex, 
    RowSideLabels = RowSideLabels,
    RowSideColors = rColors[rownames(prof_ordered)],
    labels.below = labels.below, cexRow = cex.row, cexCol = cex.col, col = colors, 
    rowsep.major = c(rep(rseps[1], extra_spacer), rseps), 
    colsep.minor = cseps, 
    colsep.major = cseps_major,
    sepwidth.minor = .01, 
    sepwidth.major = .05, 
    labels.above = labels.above, 
    na.color = 'red', 
    labels.rowsep = c(labels_right[1], rep('', extra_spacer), labels_right[2:length(labels_right)]), 
    main = main_title, cex.main = cex.main
    )
  
  hidden = heatmap.2.2(args$x, key.lab = args$key.lab, args$key.cex, #key.par = list(mai = c(.5,0,0,0)),
                       RowSideLabels = args$RowSideLabels,
                       RowSideColors = args$RowSideColors,
                       labels.below = args$labels.below, cexRow = args$cexRow, cexCol = args$cexCol, col = args$col, 
                       rowsep.major = args$rowsep.major, 
                       colsep.minor = args$colsep.minor, 
                       colsep.major = args$colsep.major,
                       sepwidth.minor = args$sepwidth.minor, 
                       sepwidth.major = args$sepwidth.major, 
                       labels.above = args$labels.above, 
                       na.color = args$na.colors, 
                       labels.rowsep = args$labels.rowsep,
                       main = args$main,
                       cex.main = args$cex.main)
  cluster_members = list()
  for(i in 0:length(rseps)){
    start = 1
    if(i > 0) start = rseps[i] + 1
    end = nrow(prof_ordered)
    if(i < length(rseps)) end = rseps[i + 1]
    
    cluster_members[[i+1]] = rownames(prof_ordered)[start:end]
  }
  ngs_res = list(cluster_members = cluster_members, colors = unique(rColors[rownames(prof_ordered)]), as_plotted = prof_ordered, args = args)
  return(ngs_res)
}

heatmap.replot_ngsplots = function(ngs_profiles, hmap_res, labels_above, labels_below){
  
  args = hmap_res$args
  
  prof = matrix(0, ncol = 0, nrow = nrow(ngs_profiles[[1]]))#assemble single matrix by joining selected profiles (matrices)
  hidden = lapply(ngs_profiles, function(x){
    prof <<- cbind(prof, x)
  })
  
  len = length(ngs_profiles)
  
  cseps = 1:(len-1) * ncol(prof) / len #column separators are between joined profiles
  
  labels.above = rep('', ncol(prof))#column labels are centered above each profile
  if(!is.na(labels_above[1])){
    labels_loc = round((1:length(labels_above) -.499) * ncol(prof) / length(labels_above))
    labels.above[labels_loc] = labels_above
  }
  
  
  labels.below = rep('', ncol(prof))
  labels_loc = round((1:length(labels_below)-.499) * ncol(prof) / length(labels_below))
  labels.below[labels_loc] = labels_below
  
  prof = prof[rownames(args$x),]
  
  
  hidden = heatmap.2.2(x = prof, key.lab = args$key.lab, key.cex = args$key.cex, #key.par = list(mai = c(.5,0,0,0)),
                       RowSideLabels = args$RowSideLabels,
                       RowSideColors = args$RowSideColors,
                       labels.below = labels.below,
                       labels.above = labels.above,
                       cexRow = args$cexRow,
                       cexCol = args$cexCol,
                       col = args$col, 
                       rowsep.major = args$rowsep.major, 
                       colsep.minor = cseps, 
                       sepwidth.minor = args$sepwidth.minor, 
                       sepwidth.major = args$sepwidth.major, 
                        
                       na.color = args$na.color, 
                       labels.rowsep = args$labels.rowsep, 
                       main = args$main,
                       cex.main = args$cex.main)
}

#dev.off()
#dummy out tracing, dummy out clustering
heatmap.2.2 = function (x,
                        col, 
                        
                        #dividing up the plot
                        colsep.minor = -5, 
                        colsep.major = -1, 
                        rowsep.minor = -1, 
                        rowsep.major = -1, 
                        sepwidth.minor = .01, 
                        sepwidth.major = .05, 
                        na.color = par("bg"), 
                        
                        #color and label for clusters
                        RowSideLabels = NA,
                        RowSideColors = NULL, 
                        
                        #title
                        main = NULL,
                        cex.main = 2,
                        
                        #labelling rows and columns
                        cexRow = 0.2 + 1/log10(nr), 
                        cexCol = 0.2 + 1/log10(nc),
                        labels.below = NULL, 
                        labels.above = NULL, 
                        labels.rowsep = NULL,
                        
                        #left margin size
                        left_mai = .7,
                        
                        #color key
                        key = T, 
                        key.height = 1.5,
                        key.lab = 'Color Key', 
                        key.cex = 1,
                        key.xtickfun = NULL, 
                        key.par = list()) 
{
  x.original = x
  rowsep.minor = -1
  rowsep.minor.original = rowsep.minor
  rowsep.major.original = rowsep.major
  colsep.minor.original = colsep.minor
  colsep.major.original = colsep.major
  #print(list(rowsep.minor.original, rowsep.major.original, colsep.minor.original, colsep.major.original))
  #sepwidth.minor = round(sepwidth.minor)
  #sepwidth.major = round(sepwidth.major)
  
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  
  
  
  dev.width = par('din')[1]
  dev.height = par('din')[2]
  
  lmat = 1
  body.height = dev.height
  body.width = dev.width
  body.iy = 1
  body.ix = 1
  lwid = dev.width
  label.height = .5
  main.height = .8
  lhei = body.height
  
  if(!is.null(labels.above)){
    lmat = rbind(max(lmat) + 1, lmat)
    lhei[body.iy] = lhei[body.iy] - label.height
    lhei = c(label.height, lhei)
    body.iy = body.iy + 1
  }
  if(!is.null(main)){
    lmat = rbind(max(lmat) + 1, lmat)
    lhei[body.iy] = lhei[body.iy] - main.height
    lhei = c(main.height, lhei)
    body.iy = body.iy + 1
  }
  if(!is.null(labels.below)){
    lmat = rbind(lmat, max(lmat) + 1)
    lhei[body.iy] = lhei[body.iy] - label.height
    lhei = c(lhei, label.height)
  }
  if(key){
    
    lmat = rbind(lmat, max(lmat) + 1)
    lhei[body.iy] = lhei[body.iy] - key.height
    lhei = c(lhei, key.height)
  }
  RowSideColors.size = 0
  if (!is.null(RowSideColors)) {
    if (!is.character(RowSideColors) || length(RowSideColors) != nrow(x)) 
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    RowSideColors.size = 1.5
    lmat <- cbind(lmat, c(rep(0, body.iy-1), max(lmat)+1, rep(0, nrow(lmat)-body.iy)))
    lwid[body.ix] = lwid[body.ix] - RowSideColors.size
    lwid <- c(lwid, RowSideColors.size)
  }
  labRowSize = 0
  if (!is.null(labels.rowsep)) {
    labRowSize = 1
    lmat <- cbind(lmat, c(rep(0, body.iy-1), max(lmat)+1, rep(0, nrow(lmat)-body.iy)))
    lwid[body.ix] = lwid[body.ix] - labRowSize
    lwid <- c(lwid, labRowSize)
  }
  if(left_mai > 0){
    lmat <- cbind(rep(0, nrow(lmat)), lmat)
    lwid[body.ix] = lwid[body.ix] - left_mai
    lwid <- c(left_mai, lwid)
    body.ix = body.ix + 1
  }
  
#   print(lmat)
#   print(lhei)
#   print(lwid)
  
  layout(lmat, heights = lhei, widths = lwid)
  
  breaks <- length(col) + 1
  symbreaks = any(x < 0, na.rm = TRUE)
  symkey = F#any(x < 0, na.rm = TRUE) || symbreaks
  
  if (length(breaks) == 1) {
    if (!symbreaks) 
      breaks <- seq(min(x, na.rm = T), max(x, na.rm = T), 
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
    if(!is.null(labels.below)) labels.below <<- c(labels.below[1:i], rep('', n), labels.below[(i + 1):length(labels.below)])
    if(!is.null(labels.above)) labels.above <<- c(labels.above[1:i], rep('', n), labels.above[(i + 1):length(labels.above)])
  }
  
  insert_row = function(n, i){#insert n rows following index i
    x <<- rbind(x[1:i,], matrix(NA, ncol = ncol(x), nrow = n),  x[(i+1):nrow(x),])
    affected = rowsep.minor >= i
    rowsep.minor[affected] <<- rowsep.minor[affected] + n
    affected = rowsep.major >= i
    rowsep.major[affected] <<- rowsep.major[affected] + n
    RowSideColors <<- c(RowSideColors[1:i], rep(NA, n), RowSideColors[(i+1):length(RowSideColors)])
  }
  
  do_col_major = colsep.major.original[1] > -1#test before value are changed
  if (colsep.minor.original[1] > -1){ #insert  columns in x and color
    for (i in 1:length(colsep.minor)){ 
      csep = colsep.minor[i]
      insert_col(round(sepwidth.minor*nc), csep)
    }
  }
  if (do_col_major){ #insert  columns in x and color
    for (i in 1:length(colsep.major)){ 
      csep = colsep.major[i]
      insert_col(round(sepwidth.major*nc), csep)
    }
  }
  x = x[nrow(x):1,]#reverse rows temporarily
  RowSideColors = rev(RowSideColors)
  do_row_major = rowsep.major[1] > -1#test before values are changed
  if (rowsep.minor.original[1] > -1){
    for (i in 1:length(rowsep.minor)){ 
      rsep = rowsep.minor[i]
      insert_row(round(sepwidth.minor*nr), rsep)
    }
  }
  if (do_row_major){
    for (i in 1:length(rowsep.major)){ 
      rsep = rowsep.major[i]
      insert_row(round(sepwidth.major*nr), rsep)
    }
  }
  
  #unreverse rows
  x = x[nrow(x):1,]
  RowSideColors = rev(RowSideColors)
  nc = ncol(x)
  nr = nrow(x)
  
  #plots the body of the heatap
  startRasterMode(width = lwid[body.ix], height = lhei[body.iy])
  image(0:nc, 0:nr, t(x), xlim = c(0, nc), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks)
  stopRasterMode(mai = rep(0,4))
  
  #plots column labels above column
  srtCol = 0
  adjCol = c(.5,1)
  colCol = NULL
  
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
        par(xpd = xpd.orig)
      }
      else warning("Invalid value for srtCol ignored.")
    }
  }
  if(!is.null(labels.above)){
    apply_col_labels(labels.above, .1, 0)
  }
  
  #draw title
  if(!is.null(main)){
    plot.new()
    par(xpd = NA)
    text(.5,.5, main, cex = cex.main)
  }
  
  #draws column labels below column
  if(!is.null(labels.below)){
    apply_col_labels(labels.below, .9, 1)
  }
  
  
  if (key) {
    mai <- c(.5, 0,.5,0)
    par(mai = mai, cex = 0.75, mgp = c(2, 1, 0))
    if (length(key.par) > 0) 
      do.call(par, key.par)
    tmpbreaks <- breaks
    if (F) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)-.001
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)+.001
    }
    else {
      min.raw <- min(x, na.rm = T)
      max.raw <- max(x, na.rm = T)
    }
    z <- seq(min.raw, max.raw, by = min(diff(breaks)/4))
    
    startRasterMode(height = lhei[body.iy])
    #draw the key color gradient
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
          xaxt = "n", yaxt = "n")
    stopRasterMode(mai = mai)
    if (is.null(key.xtickfun)) {
      lv <- -100:100 * 2
      keep = lv > min.raw & lv < max.raw
      lv = lv[keep]
#       print(lv)
      xmin = par('usr')[1]
      xmax = par('usr')[2]
      lvmin = min(lv)
      lvmax = max(lv)
      xv <- (lv - min.raw) / (max.raw - min.raw) * (xmax - xmin) + xmin
#       print(xv)
      xargs <- list(at = xv, labels = lv)
    }
    else {
      xargs <- key.xtickfun()
    }
    xargs$side <- 1
    par(xpd = NA)
    do.call(axis, xargs)
#     MAX = max(abs(c(min(x), max(x))))
#     par(usr = c(-MAX, MAX, 0, 1))
#     axis(side = 1, )
    if (!is.na(key.lab)) {
      mtext(side = 3, key.lab, line = 2, padj = 1, 
            cex = key.cex)
    }
  }
  if(!is.null(RowSideColors)){#plot colored block for each cluster
    tmp = as.factor(RowSideColors)
    vals = matrix(as.numeric(tmp), ncol = 1)
    
    lev = levels(tmp)
    startRasterMode()
    #par(mai = c(0,.1,0,.1))
    image(0:1, 0:nr, t(vals), xlim = c(0, 1), ylim = c(0, nr), axes = FALSE, xlab = "", ylab = "", col = lev)
    stopRasterMode(mai = c(0,.1,0,.1))
    cluster_levels = unique(vals)
    cluster_levels = cluster_levels[!is.na(cluster_levels)]
    cluster_ids = vals[,1]
    names(cluster_ids) = 1:length(cluster_ids)
    cluster_ids = cluster_ids[!is.na(cluster_ids)]
    par(xpd = NA)
    for(i in 1:length(RowSideLabels)){
      keep = (cluster_ids == cluster_levels[i])
      center = mean(as.numeric(names(keep[keep])))
      rowLab = rev(RowSideLabels)[i]
      text(.5,center/nr, rowLab, adj = c(.5,.5), cex = cexRow)
      
    }
    #par(xpd = NA)
    
  }
  if(!is.null(labels.rowsep)){
    par(mai = rep(0,4))
    plot0(height = nr)
    for(i in 1:length(labels.rowsep)){
      rsep_prev = 1
      if(i > 1){
        rsep_prev = rowsep.major[i - 1] + sepwidth.major
      }
      rsep_curr = rowsep.major.original[i]
      if(i > 1){
        rsep_curr = rsep_prev - 1 + rowsep.major.original[i] - rowsep.major.original[i - 1]
      }
      if(i > length(rowsep.major.original)){
        rsep_curr = nrow(x)
      }
      rsep_mid = (mean(c(rsep_prev, rsep_curr)))
      ypos = nr - rsep_mid + 1
      
      text(.5, ypos, labels.rowsep[i], adj = c(.5,.5), cex = cexRow )
    }
  }
}

do_example = F
if(do_example){
  cl = 'MCF10A'
  hm = c('H3K4AC', 'H3K4ME3')#c('H3K27AC','H3K27ME3','H3K4AC', 'H3K4ME3', 'H4K20ME3')
  toPlot = paste(cl, hm, sep = '_')
  ngs_profiles = list()
  nr = 50
  nc = 49
  for(tp in toPlot){
    dat = matrix(runif(nr * nc, -3, 7), nrow = nr, ncol = nc)
    rownames(dat) = 1:nrow(dat)
    reduce = runif(nr/2, 1, nr)
    dat[reduce,] = dat[reduce,] - 4
    dat = ifelse(dat < -3, -3, dat)
    ngs_profiles[[tp]] = dat
  }
  
  
  
  nclust = 4
  res = heatmap.ngsplots(ngs_profiles = ngs_profiles,profiles_to_plot = toPlot, label_clusters = T, fg_label = 'selected', 
                         main_title = 'title', nclust = nclust, labels_below = hm, labels_above = cl,
                         fg_toPlot = 1:10, labels_right = 'genes', sortClustersByTotal = T, labelWithCounts = T)
}