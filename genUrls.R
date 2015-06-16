
genUrls = function(){
  base = 'track name=LINE_MOD_logFE description="LINE MOD" visibility=2 autoScale=off, maxHeightPixels=16:128:128 viewLimits=0:3 color=RGB bigDataUrl=http://www.uvm.edu/~jrboyd/KZ_hESC/LINE_MOD_treat_pileup_logFE.bw type=bigWig'
  lines = c("MCF7", "MDA231")
  lines = 'H1-hESC'
  mods = c('H3K4ME3', 'H3K4AC', 'H3K27ME3', 'H3K27AC', 'H4K20ME3', 'H4K12AC')
  mods = c('H3K4me2', 'H3K4me3', 'H3K27me3', 'c-Myc')
  colors = c('252,141,98')#, '141,160,203')
  names(colors) = lines
  final = character()
  
  for(m in mods){
    for(l in lines){
      str = gsub('LINE', l, base)
      str = gsub('MOD', m , str)
      str = gsub('RGB', colors[l], str)
      print(str)
      final = paste(final, str, sep = '\n')
    }
  }
  return(final)
}

writeUrls = function(genUrls_output, fileName = 'UCSC_urls.txt'){
  write.table(genUrls_output, fileName, col.names = F, row.names = F, quote = F)
}

writeUrls(genUrls())
