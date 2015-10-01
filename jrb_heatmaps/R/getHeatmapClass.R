getHeatmapClass <-
function(resultsFromHeatmap, classIndex){
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
