plot.HeatmapLists <-
function(heatmap3_output){
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
  colSize = 50
  while(nrow(data_out) %% colSize != 0){
    data_out = rbind(data_out, '')
    colors_out = rbind(colors_out, 'black')
  }
  
  a = matrix(colors_out, nrow = colSize)
  b = matrix(data_out, nrow = colSize)
  step = 6
  while(ncol(a) %% step != 0){
    a = cbind(a, rep('black', nrow(a)))
    b = cbind(b, rep('', nrow(a)))
  }
  
  for(i in 1:(ncol(a)/6)){
    start = (i - 1) * 6 + 1
    end = i * 6
    print(start)
    print(end)
    textplot(b[,start:end], col.data = a[,start:end], show.rownames = F, show.colnames = F, halign = "left", hadj = 0, cex = .5)
  }
}
