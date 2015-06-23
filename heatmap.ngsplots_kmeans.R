#plots K27me3 and K4me3 combinations of heatmaps
#clustering is kmeans then sorted within clusters

#source('jrb_R_scripts//heatmap.3-split.R')
source('heatmap.2.2.R')
library('reshape')
library('ggplot2')

heatmap.ngsplots = function(ngs_profiles, profiles_to_plot, nclust = 6, labels_below = NA, labels_above = NA, fg_toPlot = character(), fg_label = NULL){
  if(length(labels_below) < 2 && is.na(labels_below)){
    labels_below = profiles_to_plot
  }
  if(length(labels_below) !=length(profiles_to_plot)){
    stop('length of labels_below below must be the same as profiles_to_plot!')
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
  
  
  cr = colorRamp(c('blue', 'black', 'yellow'))
  tmp = 0:100 / 100
  colors = rgb(cr(tmp)/255)
  
  len = length(profiles_to_plot)
  
  cseps = 1:(len-1) * ncol(prof) / len
  
  labColAbove = rep('', ncol(prof))
  if(!is.na(labels_above)){
    labels_loc = round((1:length(labels_above) -.499) * ncol(prof) / length(labels_above))
    labColAbove[labels_loc] = labels_above
  }
  
  
  labels = rep('', ncol(prof))
  labels_loc = round((1:len-.499) * ncol(prof) / len)
  labels[labels_loc] = labels_below
  set.seed(1)
  kmclust = kmeans(prof, centers = nclust, iter.max = 10)
  o = order(kmclust$cluster)
  prof = prof[o,, drop = F]
  kmclust$cluster = kmclust$cluster[o]
  
  for(i in 1:nclust){
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
  hiclust = hclust(dist(kmclust$centers))
  kmclust_order = hiclust$order
  prof_ordered = matrix(0, ncol = ncol(prof), nrow = 0)
  
  
  
  get_kmclust = function(i){
    start = 1
    if(i > 1){
      start = sum(kmclust$size[1:(i-1)]) + 1
    }
    end = sum(kmclust$size[1:(i)])
    return(prof[start:end,])
  }
  for(i in kmclust_order){#sort k means clusters by their center hierarchical clutering
    new_cluster = get_kmclust(i)
    start = nrow(prof_ordered) + 1
    prof_ordered = rbind(prof_ordered, new_cluster)
    end = nrow(prof_ordered)
    #kmclust$cluster[start:end] = i
    
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
  if(length(fg_toPlot) > 0){
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
      print(i)
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
  }
  
  
  res = heatmap.2.2(prof_ordered, RowSideColors = rColors[rownames(prof_ordered)], Rowv = F, Colv = F, dendrogram = 'n', trace = 'n', scale = 'n', labRow = '', labCol = labels, srtCol = 0, adjCol = c(.5,1), cexCol = 1, col = colors, symm = T, rowsep.major = rseps, colsep.minor = cseps, sepwidth.minor = 2, sepwidth.major = max(nrow(prof_ordered)/100,1), density.info = 'none', main = NULL, labColAbove = labColAbove, na.color = 'red', rowsep.major.labels = fg_label)
  return(res)
}

cl = 'MCF10A'
hm = c('H3K27AC','H3K27ME3','H3K4AC', 'H3K4ME3', 'H4K20ME3')
toPlot = paste(cl, hm, sep = '_')
ngs_profiles = list()
nr = 20
nc = 5
for(tp in toPlot){
  dat = matrix(runif(nr * nc, 0, 100), nrow = nr, ncol = nc)
  rownames(dat) = 1:nrow(dat)
  ngs_profiles[[tp]] = dat
}



#png(paste0('test8sidecol2_', cl, '.png'), width = 4, height = 8, units = 'in', res = 450)
nclust = 3
  res = heatmap.ngsplots(ngs_profiles, toPlot, nclust = 3, labels_below = hm, labels_above = cl, fg_toPlot = 1:3, fg_label = c('selected', 1:nclust))
#dev.off()