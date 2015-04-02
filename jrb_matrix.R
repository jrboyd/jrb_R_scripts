#matrix related functions

jrb.split_colnames = function(mat, split = '_', trim = '\\.'){
  #returns matrix of elements of colnames, defaults are ideal for filenames
  #mat is a matrix with colnames with descriptive elements separated by split(default '_'), and with unwanted elements following trim default('\\.')
  #split : separator of elements
  #trim : terminates desired elements
  #example : cellLine_histoneMod.asdf.fastq returns cellLine and histoneMod elements as 2 rows of matrix
  #number of elements must be the same for each colname
  cnames = colnames(mat)
  cnames = unlist(lapply(strsplit(cnames, trim), function(x)return(x[1])))
  cnames = strsplit(cnames, split)
  nrow = length(cnames[[1]])
  if(!all(lapply(cnames, length) == nrow))
    stop('number of elements in colnames not uniform!')
  out = matrix(unlist(cnames), nrow = nrow)
  return(out)
}