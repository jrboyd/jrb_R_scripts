write.HeatmapLists <-
function(heatmap3_output, name = "heatmap_lists", output_directory = ""){
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
