library(gplots)
library(xlsx)
library(RColorBrewer)
source('H:/R_workspace/jrb_R_scripts/heatmap.3-split.R')

#a wrapper for heatmap.2, primarily to add kmeans clustering(much faster for larger datasets)
#changes default colors
#outputs clusters to excel file
example_heatmap.3_kmeans_wrapper = function(n_row = 2000, n_col = 6, n_spikes = 4, n_clust = 5){
  #creates ranomized data with dimensions n_row and n_col
  #n_spikes is the number of simulated clusters added to dataset
  #n_clust is the number of clusters used by kmeans (should capture the n_spikes clusters)
  MIN = 95
  MAX = 105
  set.seed(1)
  dat = matrix(runif(n_row * n_col, min = MIN, max = MAX), nrow = n_row, ncol = n_col)
  rand =  function(n, of){
    order(runif(of))[1:n]
  }
  for(i in 1:n_spikes){
    n_choose = (1 / n_spikes)*n_row
    rows_key = rand(n_choose, n_row )
    cols_key = rand(2, of = n_col)
    bump = runif(1, min = 10, max = 30)
    print(round(bump))
    if(runif(1) > .5){
      bump = -bump
    }
    dat[rows_key, cols_key] = dat[rows_key, cols_key] + bump
  }
  heatmap.3_kmeans_wrapper(dat = dat, nclust = n_clust)
}

heatmap.3_kmeans_wrapper = function(dat, nclust = 4, hmap_colors = c("blue", 'white', 'red'), ...){
  set.seed(1)
  kclust = kmeans(dat, centers = nclust)
  o = order(kclust$cluster)
  kclust$cluster = kclust$cluster[o]
  dat = dat[o,]
  
  for(i in 1:nclust){
    keep = kclust$cluster == i
    subset = dat[keep,]
    o = order(apply(subset, 1, function(x){
      return(max(x)-min(x))
    }), decreasing = T)
    rownames(dat)[keep] = rownames(subset)[o]
    dat[keep,] = subset[o,]
  }
  
  o = order(rowSums(kclust$centers), decreasing = T)
  kclust$centers = kclust$centers[o,]
  plot_dat = matrix(0, nrow = 0, ncol = ncol(dat))
  new_cluster = integer()
  colnames(plot_dat) = colnames(dat)
  j = 1
  for(i in o){
    keep = kclust$cluster == i
    subset = dat[keep,]
    plot_dat = rbind(plot_dat, subset)
    new_cluster = c(new_cluster, rep(j, sum(keep)))
    j = j + 1
  }
  
  clust_sizes = sapply(o, function(x){
    sum(kclust$cluster == x)
  })
  
  cr = colorRamp(c('blue', 'white', 'red'))
  colors = rgb(cr(0:100/100)/255)
  # pdf("diff_mRNA_mMSCs_heatmap.pdf", width = 6, height = 8)
  seps = cumsum(sapply(1:nclust, function(x)sum(x == new_cluster)))
  
  
  override_o = cbind(1:nclust, clust_sizes)
  heatmap.3(plot_dat, trace = 'n', Rowv = F, Colv = F, scale = 'n',  cexCol = 1.6, cexRow = .4, col = colors, density.info = 'n', key.xlab = "z-score of log gene norm-counts", key.title = "", labRow = "", override_o = override_o, nsplits = 1)
  
#   library(xlsx)
#   clust_memb_wb = createWorkbook()
#   sheets = list()
#   for(i in 1:nclust){
#     new_sheet = createSheet(clust_memb_wb, paste("cluster", i))
#     keep = new_cluster == i
#     members = rownames(plot_dat)[keep]
#     df = as.data.frame(my_fpkm[members,])
#     df = cbind(sym[members], df)
#     addDataFrame(df, new_sheet)
#     
#   }
#   saveWorkbook(clust_memb_wb, "diff_mRNA_mMSCs_cluster_members.xlsx")
}

example_heatmap.3_kmeans_wrapper()
