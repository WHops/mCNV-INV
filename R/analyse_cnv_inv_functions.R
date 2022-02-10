# # Whoeps, 20th Sep 2021
# df = data.frame(Category = c('First','First','Second'), x=c(1,2,5))
# 
# df %>% group_by(Category) %>% mutate(sum = sum(x))
# blub


#' add_SD_colnames
#' @description Append colnames to SD. A bit silly to have a separate function
#' for this but it is what it is :) 
#'
#' @param SD_f [df] SD pair dataframe
#' @return SD_f, the same datamframe but with colnames.
#'
#' @author Wolfram Höps
#' @export
add_SD_colnames <- function(SD_df){
  
  SD_colnames = c('source_chr', 'source_start', 'source_end', 'orientation',
                  'target_chr','target_start', 'target_end', 'pairID',
                  'inv_chr', 'inv_start', 'inv_end', 'overlap_size')
  
  colnames(SD_df) = SD_colnames
  
  return(SD_df)
}


#' SD_keep_only_pairs_above_len
#' @description keep only pairs in which source_end is at least 'min_sd_len'-bps away from start.  
#'
#' @param SD_f [df] SD pair dataframe
#' @param min_sd_len [numeric] minimal SDpair-len [bp]
#' @return SD_f, the same datamframe but filtered. 
#'
#' @author Wolfram Höps
#' @export
SD_keep_only_pairs_above_len <- function(SD_f, min_sd_len){
  # Keep only pairs in which source_end is at least 'min_sd_len'-bps away from start.  
  SD_f = SD_f %>% group_by(pairID) %>% filter(all((source_end - source_start) > min_sd_len))

  return(SD_f)
}

#' add_sd_stats
#' @description enrich SD dataframe with additional information. 
#' Useful for easier data processing downstream. 
#'
#' @param SD_f [df] SD pair dataframe
#' @return SD_f, the same datamframe but with more columns
#'
#' @author Wolfram Höps
#' @export
add_sd_stats <- function(SD_f){
  
  SD_f$sdlen = SD_f$source_end - SD_f$source_start
  SD_f$pairlen = abs(SD_f$source_start - SD_f$target_end)
  SD_f$invlen = SD_f$inv_end - SD_f$inv_start
  
  return(SD_f)
}


#' find_SDs_with_inv_interference_potential
#' @description determine for all pairs if they have any overlap, and if all overlaps are the same.
#' The interesting pairs are the ones where this is both F.
#'
#' @param SD_f [df] SD pair dataframe
#' @return SD_f, the same datamframe but filtered
#'
#' @author Wolfram Höps
#' @export
find_SDs_with_inv_interference_potential <- function(SD_f){
  
  # Noverlap: pairs which dont overlap any invs.
  # sameverlap: pairs where both pairs overlap same inv (i.e. are nested in the inv)
  SD2 = SD_f %>% group_by(pairID) %>% 
    mutate(noverlap = all(inv_chr=='.'),
           sameverlap = var(inv_start) == 0)
  
  # the 'i' is for interesting. 
  SDi = SD2[(SD2$noverlap == F) & (SD2$sameverlap == F),]
  
  SDi$filled = SDi$sdlen / (SDi$inv_end - SDi$inv_start) 
  
  return(SDi)
}

#' prep_df_bed
#' @description # Prepare a bedfile that IGV will color 
#'
#' @param SD_f [df] SD pair dataframe
#' @return df_bed, a dataframe in bed format. 
#'
#' @author Wolfram Höps
#' @export
prep_df_bed <- function(SDi_f){
  
  df_bed = unique(SDi_f[,c('source_chr', 'source_start', 'target_end', 'orientation')])
  df_bed = df_bed[df_bed$source_start < df_bed$target_end,]
  df_bed$V5 = 0
  df_bed$V6 = df_bed$orientation
  df_bed$V7 = df_bed$source_start
  df_bed$V8 = df_bed$target_end
  df_bed$V9 = 'col'
  df_bed[df_bed$orientation == '+',]$V9 = "201,50,91"
  df_bed[df_bed$orientation == '-',]$V9 = "50,194,50"
  
  return(df_bed)
}

#' prep_df_bed2
#' @description Prepare a bedfile that IGV will color
#'
#' @param SDi_f [df] SD pair dataframe
#' @return df_bed, a dataframe in bed format. 
#'
#' @author Wolfram Höps
#' @export
prep_df_bed2 <- function(SDi_f){
  df_bed = unique(SDi_f[,c('source_chr', 'source_start', 'target_end', 'orientation')])
  #df_bed = df_bed[df_bed$source_start < df_bed$target_end,]
  df_bed$realstart = df_bed$source_start
  df_bed[df_bed$source_start > df_bed$target_end,]$realstart = df_bed[df_bed$source_start > df_bed$target_end,]$target_end
  df_bed$realend = df_bed$target_end
  df_bed[df_bed$source_start > df_bed$target_end,]$realend = df_bed[df_bed$source_start > df_bed$target_end,]$source_start
  
  df_bed$source_start = df_bed$realstart
  df_bed$target_end = df_bed$realend
  
  df_bed$realstart = NULL
  df_bed$realend = NULL
  df_bed$V5 = 0
  df_bed$V6 = df_bed$orientation
  df_bed$V7 = df_bed$source_start
  df_bed$V8 = df_bed$target_end
  df_bed$V9 = 'col'
  df_bed[df_bed$orientation == '+',]$V9 = "201,50,91"
  df_bed[df_bed$orientation == '-',]$V9 = "50,194,50"
  
  return(df_bed)
}


#' filter_down_SD_INV_list
#' @description 
#'
#' @param SDi_with_GT [df] SD pair dataframe
#' @return df_bed, a dataframe in bed format. 
#'
#' @author Wolfram Höps
#' @export
filter_down_SD_INV_list <- function(SDi_with_GT){
  
  
  SDi_simp = SDi_with_GT[SDi_with_GT$n_switches == 1,]
  
  invcenter = SDi_simp %>% group_by(inv_start) %>% slice(1)
  invcenter = invcenter[invcenter$nnoreads < 10,]
  invcenter = invcenter[invcenter$nlowconf < 10,]
  return(invcenter)
}

#' filter_down_SD_INV_list_nofilter_but_turn
#' @description 
#'
#' @param SDi_with_GT [df] SD pair dataframe
#' @return  
#'
#' @author Wolfram Höps
#' @export
filter_down_SD_INV_list_nofilter_but_turn <- function(SDi_with_GT){
  
  SDi_simp = SDi_with_GT[SDi_with_GT$n_switches == 1,]
  
  invcenter = SDi_simp %>% group_by(inv_start) %>% slice(1)
  invcenter = invcenter[invcenter$nnoreads < 10,]
  invcenter = invcenter[invcenter$nlowconf < 10,]
  return(invcenter)
}

#' mark_protective_risk_mixed_invs
#' @description 
#'
#' @param invcenter
#' @param n1
#' @param n2
#' @return  
#'
#' @author Wolfram Höps
#' @export
mark_protective_risk_mixed_invs <- function(invcenter, n1, n2){
  
  invcenter$mix = 'Mixed'
  invcenter[invcenter$pctpos < n1,]$mix = 'SV risk factor'
  invcenter[invcenter$pctpos > n2,]$mix = 'Protective'
  
  # While we are at it, also note the length from SD1 to SD2
  invcenter$SDpairlen = abs(invcenter$source_start - invcenter$target_start)
  
  return(invcenter)
}

#' make_plot
#' @description 
#'
#' @param invcenter
#' @param limit
#' @return  
#'
#' @author Wolfram Höps
#' @export
make_plot <- function(invcenter, limit){
  p = ggplot(invcenter[(invcenter$inv_end-invcenter$inv_start < limit) & (invcenter$mix != 0),]) +
    geom_boxplot(aes(group=mix, x=mix, y= (nhet+(2*nhom)) / (43*2), fill=mix)) + scale_color_viridis_b() +
    geom_beeswarm(aes(x=mix, y=(nhet+(2*nhom))/(43*2))) + 
    theme_bw() + 
    labs(x='', y='% Inverted alleles', title='Inversions < 200kb')
  return(p)
}


#' make_wilcox
#' @description run t.test
#'
#' @param df
#' @param group1
#' @param group2
#' @return  
#'
#' @author Wolfram Höps
#' @export
make_wilcox <- function(df, group1, group2){
  df$inv_alleles = df$nhet + 2*df$nhom
  x = df[df$mix == group1,]$inv_alleles
  y = df[df$mix == group2,]$inv_alleles
  wilcox.test(x,y, alternative='greater', exact=F)
}

#' moveme
#' @description A helperfunction
#'
#' @param df
#' @param group1
#' @param group2
#' @return  
#'
#' @author Ananda Mahto,
#' https://gist.github.com/mrdwab/7183841
#' @export
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

#' determine_directionality_of_affected_SDs_below_inv
#' @description A helperfunction
#'
#' @param df
#' @param group1
#' @param group2
#' @return  
#'
#' @author Wolfram Höps
#' @export
determine_directionality_of_affected_SDs_below_inv <- function(SD_df){
  SDi_chromcenter = SD_df %>% group_by(inv_start) %>% mutate(
    inv_n_alterations = length(inv_chr),
    allp = all(orientation == '+'),
    alln = all(orientation == '-'),
    pctpos = sum((abs(
      target_end - source_start
    ))[orientation == '+']) / sum(abs(target_end - source_start)),
    mix =
      any(allp) + 2 * any(alln)
  )
  
  return(SDi_chromcenter)
}

#' SDi_chromcenter_nswitch_gts
#' @description Add genotypes to SD dataframe
#'
#' @param SD_df
#' @param gt_link
#' @return  SD_df but with GT fields
#'
#' @author Wolfram Höps
#' @export
add_gts <- function(SD_df, gt_link){
  
  gts_cols_to_keep = c(
    'seqnames',
    'start',
    'end',
    'nref',
    'nhet',
    'nhom',
    'ncomplex',
    'ninvdup',
    'nnoreads',
    'nlowconf'
  )
  
  # Bring GTs into the game.
  gts = read.table(gt_link, sep = '\t', header = T)[, gts_cols_to_keep]
  colnames(gts)[1:3] = c('inv_chr', 'inv_start', 'inv_end')
  
  # Join...
  SD_df_gts = left_join(SD_df,
                        gts,
                        by = c('inv_chr', 'inv_start', 'inv_end'))
  
  return(SD_df_gts)
  
}

#' mark_inv_class
#' @description mark_inv_class
#'
#' @param inv_df
#' @param inv_noSD_link
#' @return inv_df but with class info
#'
#' @author Wolfram Höps
#' @export
mark_inv_class <- function(inv_df, inv_noSD_link){
  invnosd = read.table(inv_noSD_link, header = F, sep = '\t')
  colnames(invnosd) = c('inv_chr' , 'inv_start'  , 'inv_end', 'class')
  inv_df = left_join(inv_df, invnosd, by = c('inv_chr' , 'inv_start'  , 'inv_end'))
  
  return(inv_df)
}

#' part1_find_all_SDpairs_affected_by_inv
#' @description # Part 1: find all SD pairs which are affected by an inversion.
#' This produces two tables that can be viewed in IGV
#' @param
#' @return
#'
#' @author Wolfram Höps
#' 
#' @export
part1_find_all_SDpairs_affected_by_inv <- function(SD_link,
                                                   params,
                                                   outdir) {
  # Load SD data table
  SD = add_SD_colnames(read.table(SD_link))
  
  # Cut down to min len
  SD_minlen = SD_keep_only_pairs_above_len(SD, params$min_sd_len)
  
  # Add some information on sdlen, pairlen, invlen
  SD_minlen_stats = add_sd_stats(SD_minlen)
  
  # Exclude extremely long CNVs (>10Mb), which are not expect to appear.
  SD_minlen_stats_pairlen = SD_minlen_stats[SD_minlen_stats$pairlen < params$max_pairlen,]
  
  # Keep only SDs where an inv has a potential for orientation change.
  # This means one of them, but not both, is overlapped by the same inversion.
  # These are the 'interesting' SDs (SDi).
  SDi = find_SDs_with_inv_interference_potential(SD_minlen_stats_pairlen)
  
  SDi_save = SDi; SDi_save[,c('noverlap','sameverlap','filled')] = NULL
  
  # Made two bedfiles for IGV.
  df_bed = prep_df_bed(SDi)
  SDi_bed = SDi[, c('source_chr', 'source_start', 'source_end')]
  
  
  
  # Save results
  if (params$save) {
    dir.create(outdir, showWarnings = F)
    write.table(
      SDi_save,
      file = paste0(outdir, '/sheetE_SDi.tsv'),
      col.names = T,
      row.names = F,
      quote = F,
      sep = '\t'
    )
    # write.table(
    #   file = paste0(outdir, "/all_flipped_SD_pairs.bed"),
    #   na.omit(df_bed),
    #   col.names = F,
    #   row.names = F,
    #   quote = F,
    #   sep = '\t'
    # )
    # write.table(
    #   file = paste0(outdir, "/all_SDs_in_flipped_pairs.bed"),
    #   na.omit(SDi_bed),
    #   col.names = F,
    #   row.names = F,
    #   quote = F,
    #   sep = '\t'
    # )
  }
  
  return(SDi)
  
}




#' part2_find_protective_and_risk_inversions
#' @description # Part 2: Find 'protective' and 'risk' Inversions:
#' Those are inversions that lead to creation of only inverted or only direct SDs.
#' I.e. non-messy sites :)
#' @param
#' @return
#'
#' @author Wolfram Höps
#' @export
part2_find_protective_and_risk_inversions <-
  function(SDi, params, outdir) {
    
    SD_cols_to_keep = c(
      'source_chr',
      'source_start',
      'source_end',
      'target_end',
      'target_start',
      'orientation',
      'inv_chr',
      'inv_start',
      'inv_end',
      'pairID'
    )
    
    colnames_mcnv_list = c('chr',
                           'start',
                           'end',
                           'mCNV name',
                           'genome',
                           'gt',
                           'del/dup',
                           'len')
    
    # Filter to only SDs being touched by Inversions (As opposed to keeping the pairs.)
    SDi_touch_inv = SDi[SDi$inv_chr != '.',]
    
    SDi_regions = SDi_touch_inv[, SD_cols_to_keep]
    
    # Inversion-wise statistics: are all affected SD pairs in pos or neg directions?
    # How many SD pairs does this inversion affect?
    SDi_chromcenter = determine_directionality_of_affected_SDs_below_inv(SDi_regions)
    
    # Another caveat: There are SD pairs that are affected by >1 inversions.
    # The 'nswitch' parameter explains per SD pair, how many inversions can affect its orientation.
    SDi_chromcenter_nswitch = SDi_chromcenter %>% group_by(pairID) %>% mutate(n_switches = length(inv_chr))
    
    if (params$use_only_high_confidence_invs) {
      SDi_chromcenter_nswitch_gts = add_gts(SDi_chromcenter_nswitch, gt_link)
      
      # This is a big step here. We want to filter our list of affected SDs by:
      # - SD pair should only be affected by ONE inversion (otherwise interpretation gets tricky)
      # - Inversion should have <10 GTs with noreads or lowconf reads.
      SDi_chromcenter_nswitch_gts_goodsamples = filter_down_SD_INV_list(SDi_chromcenter_nswitch_gts)
    } else {
      SDi_chromcenter_nswitch_gts_goodsamples =  SDi_chromcenter_nswitch %>% group_by(inv_chr, inv_start, inv_end) %>% slice(1)
    }
    
    
    
    # Now, mark which inversions are likely protective, risky and mixed.
    # We consider 'protective' if  < sds_protect_ro_risk_th of SD-pair-basepairs switch from protect to risk
    # We consider 'risky' if > 1-x% of SD-pair-basepairs switch from protect to risk
    inv_df = mark_protective_risk_mixed_invs(
      invcenter = SDi_chromcenter_nswitch_gts_goodsamples,
      n1 = params$sds_protect_ro_risk_th,
      n2 = 1 - params$sds_protect_ro_risk_th
    )
    
    # More info still. We also want to know the inversion class. SD-mediated? nonSD-mediated?
    inv_df_class = mark_inv_class(inv_df, inv_noSD_link)
    
    # Save that thing
    cols_sheetA_keep = c(
      'inv_chr','inv_start','inv_end','inv_n_alterations','pctpos','mix'
    )
    write.table(
      inv_df_class[,cols_sheetA_keep],
      file = paste0(outdir, '/sheetA_flipping_inversions.tsv'),
      sep = '\t',
      quote = F,
      col.names = T,
      row.names = F
    )
    
    
    
    # Add mCNVs
    mcnv_list = read.table(mcnv_link, header = F, sep = '\t')
    colnames(mcnv_list) = colnames_mcnv_list
    mcnv_list$chr = paste0('chr', mcnv_list$chr)
    
    # Bedtools
    inv_df_class = inv_df_class[moveme(names(inv_df_class),
                                       "inv_chr first; inv_start after inv_chr; inv_end after inv_start")]
    inv_df_protect_risk = inv_df_class[inv_df_class$mix != 'Mixed',]
    
    # Do mCNVs better...
    SDs_flipped = left_join(SDi,
                            inv_df_protect_risk,
                            by = c('inv_chr', 'inv_start', 'inv_end'))
    SDs_flipped2 = SDs_flipped[!is.na(SDs_flipped$SDpairlen),]
    
    # Add recurrence
    rec = read.table(recurrence_link, header = T, sep = '\t')
    colnames(rec)[1:3] = c('inv_chr', 'inv_start', 'inv_end')
    rec2 = rec[, c(
      'inv_chr',
      'inv_start',
      'inv_end',
      'verdictRecurrence_hufsah',
      'verdictRecurrence_benson'
    )]
    
    # Sheet A: 
    test = left_join(inv_df_class, rec2, by = c('inv_chr', 'inv_start', 'inv_end'))
    
    
    inv_df_protect_risk = left_join(inv_df_protect_risk,
                                    rec2,
                                    by = c('inv_chr', 'inv_start', 'inv_end'))
    
    # SAVE
    write.table(
      inv_df_protect_risk[order(inv_df_protect_risk$inv_chr),],
      file = paste0(outdir, '/sheetB_protect_risk_loci.tsv'),
      sep = '\t',
      quote = F,
      col.names = T,
      row.names = F
    )
    
    
    overlap_main = bedtoolsr::bt.intersect(inv_df_protect_risk, mcnv_list, wao =
                                             T)
    colnames(overlap_main) = c(colnames(inv_df_protect_risk), colnames(mcnv_list))
    
    SDs_flipped2$min = rowMins(as.matrix(SDs_flipped2[, c('source_start.x',
                                                          'target_start.x',
                                                          'source_end.x',
                                                          'target_end.x')]))
    SDs_flipped2$max = rowMaxs(as.matrix(SDs_flipped2[, c('source_start.x',
                                                          'target_start.x',
                                                          'source_end.x',
                                                          'target_end.x')]))
    
    SDs_flipped3 = SDs_flipped2[moveme(names(SDs_flipped2),
                                       "source_chr.x first; min after source_chr.x; max after min")]
    
    overlap_sec = bedtoolsr::bt.intersect(SDs_flipped3, mcnv_list, wao = T)
    colnames(overlap_sec) = c(colnames(SDs_flipped3), colnames(mcnv_list))
    
    overlap_sec_to_save_cols = c(
      'source_chr.x',
      'min',
      'max',
      'orientation.x',
      'pairID.x',
      'inv_chr',
      'inv_start',
      'inv_end',
      'source_start.y',
      'mix',
      'n_switches',
      'chr',
      'start',
      'end',
      'mCNV name',
      'gt',
      'del/dup',
      'len'
    )
    
    overlap_save = overlap_sec[(overlap_sec$n_switches == 1), overlap_sec_to_save_cols]
    overlap_save2 = overlap_save[(overlap_save$`mCNV name` != '.'),]
    write.table(
      overlap_save2,
      file = paste0(outdir, '/sheetD_flipped_SD_pairs_overlapping_mCNVs.tsv'),
      row.names = F,
      col.names = T,
      sep = '\t',
      quote = F
    )
    overlap_main2 = dplyr::left_join(overlap_main, rec2, by = c('inv_chr', 'inv_start', 'inv_end'))
    
    
    # # To_save
    cols_to_save = c(
      'inv_chr',
      'inv_start',
      'inv_end',
      'orientation',
      'inv_n_alterations',
      'mix',
      'nref',
      'nhet',
      'nhom',
      'ncomplex',
      'ninvdup',
      'nnoreads',
      'nlowconf',
      'class',
      'mCNV name'
    )
    
    colnames(overlap_save) = stringr::word(colnames(overlap_save),1,sep = "\\.")
    overlap_save = overlap_main2[, cols_to_save]
    overlap_save = overlap_save[overlap_save$`mCNV name` != '.',]
    colnames(overlap_save)[colnames(overlap_save) == 'orientation'] = 'SDpairs_orientation'
    colnames(overlap_save)[colnames(overlap_save) == 'inv_n_alterations'] = 'SDpairs_flipped'
    colnames(overlap_save)[colnames(overlap_save) == 'class'] = 'Inversion_class'
    colnames(overlap_save)[colnames(overlap_save) == 'mix'] = 'inv_role'
    
    write.table(
      overlap_save,
      file = paste0(outdir, '/sheetC_invs_affecting_sds_affecting_mcnvs_strict.tsv'),
      row.names = F,
      col.names = T,
      quote = F,
      sep = '\t'
    )
    ### Now
    # Made a bedfile for igv
    # inv_df_bed = prep_df_bed2(inv_df[inv_df$mix != 'Mixed',])
    # write.table(
    #   inv_df_bed,
    #   file = paste0(outdir, '/bedversion_sheetB_protect_risk.bed'),
    #   sep = '\t',
    #   col.names = F,
    #   row.names = F,
    #   quote = F
    # )
    # Make another bedfile for igv^^
    SDi_bed = SDi[, c('source_chr', 'source_start', 'source_end')]
    
  }

