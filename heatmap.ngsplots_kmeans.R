#plots K27me3 and K4me3 combinations of heatmaps
#clustering is kmeans then sorted within clusters

#source('jrb_R_scripts//heatmap.3-split.R')
source('heatmap.2.1.R')
library('reshape')
library('ggplot2')

heatmap.ngsplots = function(ngs_profiles, profiles_to_plot, nclust = 6, labels_below = NA, labels_above = NA){
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
  
  
  
  cr = colorRamp(c('blue', 'black', 'yellow'))
  tmp = 0:100 / 100
  colors = rgb(cr(tmp)/255)
  
  len = length(profiles_to_plot)
  
  cseps = 1:(len-1) * ncol(prof) / len
  
  labColAbove = rep('', ncol(prof))
  if(!is.na(labels_above)){
    labels_loc = round((1:length(labels_above) -.5) * ncol(prof) / length(labels_above))
    labColAbove[labels_loc] = labels_above
  }
  
  
  labels = rep('', ncol(prof))
  labels_loc = round((1:len-.5) * ncol(prof) / len)
  labels[labels_loc] = labels_below
  
  kmclust = kmeans(prof, centers = nclust, iter.max = 10)
  o = order(kmclust$cluster)
  prof = prof[o,]
  kmclust$cluster = kmclust$cluster[o]
  
  for(i in 1:nclust){
    size = kmclust$size
    start = 1
    if(i > 1){
      start = sum(kmclust$size[1:(i-1)]) + 1
    }
    end = sum(kmclust$size[1:(i)])
    #print(paste(start, end))
    o = order(rowSums(prof[start:end,]), decreasing = T)
    prof[start:end,] = prof[start:end,][o,]
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
  for(i in kmclust_order){
    prof_ordered = rbind(prof_ordered, get_kmclust(i))
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
      rseps[i] = end-1
    }
  }
  
  res = heatmap.2.2(prof_ordered, Rowv = F, Colv = F, dendrogram = 'n', trace = 'n', scale = 'n', labRow = '', labCol = labels, srtCol = 0, adjCol = c(.5,1), cexCol = 1, col = colors, symm = T, rowsep.major = rseps, colsep.minor = cseps, sepwidth.minor = 2, sepwidth.major = 150, density.info = 'none', main = NULL, labColAbove = labColAbove, na.color = 'red')
  return(res)
}

# cl = 'MCF10A'
# hm = c('H3K27AC','H3K27ME3','H3K4AC', 'H3K4ME3', 'H4K20ME3')
# toPlot = paste(cl, hm, sep = '_')
# 
# png(paste0('test8_', cl, '.png'), width = 4, height = 8, units = 'in', res = 450)
#   res = heatmap.ngsplots(ngs_profiles, toPlot, nclust = 8, labels_below = hm, labels_above = cl)
# dev.off()