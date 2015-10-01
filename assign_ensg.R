assign_ensg = function(test_res, ref_dict, overlap = 0) {
    if(overlap > 1 | overlap < 0){
      stop('overlap must be between 0 and 1')
    }
    all_chrms = union(unique(test_res[, 1]), unique(ref_dict$chrm))
    output = character()  #matrix(0, ncol = ncol(test_res), nrow = 0)
    for (chrm in all_chrms) {
        print(chrm)
        output = c(output, assign_ensg_chrm(test_res, ref_dict, chrm, overlap))  #rbind(output, assign_ensg_chrm(test_res, ref_dict, chrm))
        # break
    }
    return(output)
}

assign_ensg_chrm = function(test_res, ref_dict, chrm, overlap) {
    # columns must be chrm, start, end, ref_dict rownames are ids sort all by start then by chrm returns
    # test_res filtered by intersection with ref_dict, rownames assigned by ref_dict
    test_res = test_res[test_res[, 1] == chrm, ]
    ref_dict = ref_dict[ref_dict$chrm == chrm, ]
    
    if (nrow(test_res) == 0 || nrow(ref_dict) == 0) {
        print(paste(chrm, "is emtpy"))
        return(matrix(0, ncol = ncol(test_res), nrow = 0))
    }
    
    o = order(ref_dict$start)
    ref_dict = ref_dict[o, ]
    
    o = order(test_res[, 2])
    test_res = test_res[o, ]
    
    all_ensg = character()
    # intersect ma_res with ensgs, assigning ensg or removing peaks
    for (i in 1:nrow(test_res)) {
        # iterate test_res
        test_start = test_res[i, 2]
        test_end = test_res[i, 3]
        
        ref_ensg = rownames(ref_dict)
        ref_starts = ref_dict$start
        ref_ends = ref_dict$end
        
        ends_after_start = test_end >= ref_starts
        starts_before_end = test_start <= ref_ends
        
        is_overlap = ends_after_start & starts_before_end
        intersecting_ensg = ref_ensg[is_overlap]
        #need to test for overlap
        if(F){#debug values for overlap
          test_start = 0; test_end = 100
          ref_starts = c(10, -50, 100, 70, -80, 110, -500); ref_ends = c(20, 10, 150, 120, -70, 120, 200)
          answers = c(10, 10, 0, 30, 0, 0, 100)#n bp
        }
        if(overlap > 0 & length(intersecting_ensg) > 0){
          max_starts = sapply(ref_starts[is_overlap], function(x){
            max(x, test_start)
          })
          min_ends = sapply(ref_ends[is_overlap], function(x){
            min(x, test_end)
          })
          bp_overlap = min_ends - max_starts
          bp_overlap = min_in_place(bp_overlap, 0)
          test_fraction = bp_overlap / (test_end - test_start)
          ref_fraction = bp_overlap / (ref_ends[is_overlap] - ref_starts[is_overlap])
          keep = (test_fraction >= overlap) | (ref_fraction >= overlap)
          # print(sum(keep))
          intersecting_ensg = intersecting_ensg[keep]
        }
        
        
        if(length(intersecting_ensg) == 0){
          intersecting_ensg = paste0("no_match_", rownames(test_res)[i])
        }
        all_ensg = c(all_ensg, paste(intersecting_ensg, collapse = ';'))
    }
    return(all_ensg)
} 

min_in_place <- function(x, a) {
  x[x <= a] <- a
  return(x)
}

resolve_overlaps_in_ref = function(ref_dict){
  overlap_perc = .8
  new_dict = ref_dict[integer(),]
  for(chrm in unique(ref_dict$chrm)){
    print(chrm)
    keep = ref_dict$chrm == chrm
    chrm_dict = ref_dict[keep,]
    
    i = 1
    if(nrow(chrm_dict) == 1){
      new_dict = rbind(new_dict, chrm_dict)
    }else while(i <= nrow(chrm_dict)){
      i_start = chrm_dict$start[i]
      i_end = chrm_dict$end[i]
      after_i = chrm_dict[(i+1):nrow(chrm_dict),]
      is_overlap = i_start < after_i$end & i_end > after_i$start
      if(i == nrow(chrm_dict) || sum(is_overlap) == 0){
        new_dict = rbind(new_dict, chrm_dict[i,])
        i = i + 1#do nothing, not overlap
      }else{
        start = i
        end = i + sum(is_overlap)
        overlap_dict = chrm_dict[start:end,]  
        #print(overlap_dict)
        new_dict = rbind(new_dict, overlap_dict[nrow(overlap_dict)/2,])
        #print(chrm_dict[start:end,])
        i = i + sum(is_overlap)+1#is overlap, find block
      }
    }
  }
  return(new_dict)
}
