example_heatmap.3_kmeans_wrapper <-
function(n_row = 2000, n_col = 6, n_spikes = 4, n_clust = 5){
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
