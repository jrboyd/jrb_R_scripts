
countDir.read = function(countDir){
  #given a directory containing counts files, return assembled matrix containing all data
  #count files must have same number of rows in the same order (no changing feature annotation)
  #count_unaligned : htseq-count outputs 5 rows of unaligned read types, should they be counted?
  #NOTE : keep as true for counts not from htseq-count
  files = dir(countDir, full.names = T)
  
  dat = NULL
  rnames = NULL
  cnames = dir(countDir)
  
  for(f in files){
    tmp = read.table(f, row.names = 1)
    if(is.null(dat)){
      dat = matrix(0, nrow = nrow(tmp), ncol = 0)
      rnames = rownames(tmp)
    }
    dat = cbind(dat, tmp)
  }
  
  colnames(dat) = cnames
  rownames(dat) = rnames
  return(dat)
}



countDir.example = function(){
  d = 'h:/projects/Terri/paper_v4/chrm_counts'
  tmp = countDir.read(d)
  return(tmp)
}


if(F){
  indexDict = read.table("ref/hg38.index2symbol.1kb_ext_promoters.txt", row.names = 1, stringsAsFactors = F)
  ref2sym = as.matrix(indexDict[, 1])
  rownames(ref2sym) = rownames(indexDict)
  
  ref2lengths = indexDict[, 2]
  names(ref2lengths) = rownames(indexDict)
  
  countsDir = 'counts'
  countFiles = dir(countsDir, full.names = T)
  data = read.table(countFiles[1])[,2, drop = F]
  rownames(data) = read.table(countFiles[1])[,1]
  for(cf in countFiles[2:length(countFiles)]){
    dat = read.table(cf)[,2]
    data = cbind(data, dat)
  }
  
  indi_header = countFiles
  indi_header = matrix(unlist(strsplit(indi_header, split = "[_/.]")), ncol = ncol(data))
  
  # fix mcf7, remove high count sequencing rep for mcf7 4ac
  
  l_row = 2
  m_row = 3
  r_row = 4
  
  lines = unique(indi_header[l_row, ])
  mods = unique(indi_header[m_row, ])
  reps = unique(indi_header[r_row, ])
  
  
  
  print("aggregating counts...")
  agg_data = matrix(0, nrow = nrow(data), ncol = length(lines) * length(mods))
  colnames(agg_data) = 1:ncol(agg_data)
  rownames(agg_data) = rownames(data)
  i = 1
  for (l in lines) {
    for (m in mods) {
      keep = (indi_header[l_row, ] == l & indi_header[m_row, ] == m)
      agg_data[, i] = rowSums(data[, keep, drop = F])
      colnames(agg_data)[i] = paste(l, m)
      i = i + 1
    }
  }
  
  
  
  # calc cpm
  print("normalizing per million reads...")
  mod_data = agg_data[1:(nrow(agg_data) - 4), ]
  cpm_factors = 1e+06/colSums(mod_data)
  mod_data = t(t(mod_data) * cpm_factors)
  # mod_data = mod_data[1:(nrow(mod_data)-5),] norm for feature length
  mod_data = mod_data[1:(nrow(mod_data) - 1), ]
  
  featureLengths = ref2lengths[rownames(mod_data)]
  print("normalizing for kb feature length...")
  mod_data = mod_data/featureLengths * 1000
  
  print("applying log2 transform...")
  rpkm_data = log2(mod_data + 1)
  fe_data = matrix(0, nrow = nrow(rpkm_data), ncol = length(lines) * (length(mods) - 1))
  keep = !grepl("input", colnames(rpkm_data))
  colnames(fe_data) = colnames(rpkm_data)[keep]
  rownames(fe_data) = rownames(rpkm_data)
  
  print("calculating FE over input...")
  for (i in 1:ncol(fe_data)) {
    key = colnames(fe_data)[i]
    l = strsplit(key, " ")[[1]][1]
    key_input = paste(l, "input")
    # print(paste(key, 'and', key_input))
    fe_data[, i] = rpkm_data[, key] - rpkm_data[, key_input]
  }
  
  print("flooring data at 1 FE, log2(1) = 0...")
  fe_data = ifelse(fe_data < 0, 0, fe_data)
  print("removing undetected entries...")
  keep = apply(fe_data, 1, max) > 0
  fe_data = fe_data[keep, ]
  print("sorting columns by histone mark (4ac 4me3) then by line (mcf10a, mcf7, mda231)...")
  markData_4me3_4ac = fe_data[, c(1, 3, 5, 2, 4, 6)]
  save(markData_4me3_4ac, file = "oldcounts_data.save")
  print("markData_4me3_4ac is set and saved!") 
}
