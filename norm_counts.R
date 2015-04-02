norm.read_depth = function(counts, reads_per = 10^6){
  #performs cpm reads type normalization by sample
  #counts : feature x sample matrix
  #eachs column of counts is assumed to be one sample
  #defaults to reads_per million (10^6)
  counts = apply(counts, 2, function(x)return(x / sum(x) * reads_per))
  return(counts)
}

norm.feature_length = function(counts, feature_lengths, reads_per = 10^3){
  #performs rpk type normalization by feature
  #counts : feature x sample matrix, rows must have names for each feature to cross to feature_lengths
  #feature_lengths : named vector of feature lengths, must have entry matching each named row in counts
  #reads_per : defaults to reads per kilobase (10^3)
  
  #sort feature_lengths
  feature_lengths = feature_lengths[rownames(counts),,drop = F]
  
  tmp = cbind(counts, feature_lengths)
  
  counts = t(apply(tmp, 1, function(x)return(x[1:(length(x)-1)] / x[length(x)] * reads_per)))
  
  return(counts)
}

norm.remove_unaligned = function(counts, n = 5){
  #removes last n rows of counts
  #htseq-count returns 5 rows of unaligned read counts
  return(counts[1:(nrow(counts)-n),])
}

norm.fe = function(counts){
  #assumption that colnames have the form 'cell_line histone_mod'
  #each non-input histone_mod will be divided by cell_line specific input
  feData = matrix(0, nrow = nrow(counts), ncol = 0)
  for(i in 1:ncol(counts)){
    cname = colnames(counts)[i]
    l = strsplit(cname, ' ')[[1]][1]
    m = strsplit(cname, ' ')[[1]][2]
    if(m == 'input'){
      next
    }
    fe_dat = counts[,i] / counts[,paste(l, 'input')]
    feData = cbind(feData, fe_dat)
    colnames(feData)[ncol(feData)] = cname
  }
  return(feData)
}

