plot.hmap_res_lists = function(hmap_res, nper_column = 50, nper_row = 10, cex = 1, col_spacing = 1, row_spacing = 1){
  res = hmap_res
  class_members = res$cluster_members
  class_colors = res$colors
  class_sizes = sapply(class_members, length)
  #   classSizes = res[[1]]
  #   asPlotted = res[[3]]
  #   class_colors = res[[4]]
  data_out = matrix(unlist(class_members), ncol = 1)
  colors_out = matrix('black', nrow = sum(class_sizes), ncol = 1)
  rownames(colors_out) = unlist(class_members)
  for(i in 1:length(class_members)){
    colors_out[class_members[[i]],] = class_colors[i]
  }  
  #   soFar = 0
  #   for(j in 1:length(classSizes)){
  #     size = classSizes[j]
  #     start = soFar + 1
  #     end = soFar + size
  #     soFar = soFar + size
  #     colors_out[start:end,] = class_colors[j]
  #   }
  #grow nrow of matrices to nper_column
  while(nrow(data_out) %% nper_column != 0){
    data_out = rbind(data_out, '')
    colors_out = rbind(colors_out, 'black')
  }
  
  color_data = matrix(colors_out, nrow = nper_column)
  names_data = matrix(data_out, nrow = nper_column)
  #grow ncol to nper_row
  while(ncol(color_data) %% nper_row != 0){ #grow matrices to nper_row
    color_data = cbind(color_data, rep('black', nrow(color_data)))
    names_data = cbind(names_data, rep('', nrow(color_data)))
  }
  #rev rows for plotting
#   color_data = color_data[nrow(color_data):1,]
#   names_data = names_data[nrow(names_data):1,]
  
  for(i in 1:(ncol(color_data)/nper_row)){
    
    start = (i - 1) * nper_row + 1
    end = i * nper_row
    #     print(start)
    #     print(end)
    names_tmp = names_data[,start:end, drop = F]
    color_tmp = color_data[,start:end, drop = F]
#     print(names_tmp)
#     print(color_tmp)
    plot(1, xlim = c(0,(ncol(names_tmp) + 1)), ylim = c((nrow(names_tmp) + 1), 0), type = 'n', xlab = '', ylab = '', axes = F)
    for(x in 1:ncol(names_tmp)) for(y in 1:nrow(names_tmp)){
      if(color_tmp[y,x] == 'white'){
        rect((x-.5) * col_spacing,
             (y-.5) * row_spacing,
             (x+.5) * col_spacing,
             (y+.5) * row_spacing,
             bty = 'n', col = 'black',border = "transparent")
      }
    }
    for(x in 1:ncol(names_tmp)) for(y in 1:nrow(names_tmp)){
      
      text(x * col_spacing,y * row_spacing,names_tmp[y,x], col = color_tmp[y,x], cex = cex, adj = c(.5,.5))
    }
    
    
  }
}