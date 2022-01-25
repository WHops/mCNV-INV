# # Whoeps, 20th Sep 2021
# df = data.frame(Category = c('First','First','Second'), x=c(1,2,5))
# 
# df %>% group_by(Category) %>% mutate(sum = sum(x))
# blub


# Functions for mCNV-INV-SD analysis.
# To be documented better...

cut_down_len <- function(SD_f, params_f){
  #SD_f = tibble(SD[1:10,])
  #SD_f$lendissatisfy = (SD_f$source_end - SD_f$source_start) < params_f$min_sd_len
  SD_f = SD_f %>% group_by(pairID) %>% mutate(valid = all(source_end - source_start > params_f$min_sd_len))
  #SD_f$lds = as.data.frame(xtabs(lendissatisfy ~ pairID, SD_f))$Freq
  SD_f = SD_f[SD_f$valid > 0,]
  return(SD_f)
}

add_sd_stats <- function(SD_f){
  SD_f$sdlen = SD_f$source_end - SD_f$source_start
  SD_f$pairlen = abs(SD_f$source_start - SD_f$target_end)
  SD_f$invlen = SD_f$inv_end - SD_f$inv_start
  return(SD_f)
}

# Determine for all pairs if they have any overlap, and if all overlaps are the same.
# The interesting pairs are the ones where this is both F.
find_SDs_with_inv_interference_potential <- function(SD_f){
  SD2 = SD_f %>% group_by(pairID) %>% 
    mutate(noverlap = all(inv_chr=='.'),
           sameverlap = var(inv_start) == 0)
  SDi = SD2[(SD2$noverlap == F) & (SD2$sameverlap == F),]
  SDi$filled = SDi$sdlen/ (SDi$inv_end - SDi$inv_start) 
  return(SDi)
}

# Prepare a bedfile that IGV will color 
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

# Prepare a bedfile that IGV will color 
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


filter_down_SD_INV_list <- function(SDi_with_GT){
  
  SDi_simp = SDi_with_GT[SDi_with_GT$n_switches == 1,]
  
  invcenter = SDi_simp %>% group_by(inv_start) %>% slice(1)
  invcenter = invcenter[invcenter$nnoreads < 10,]
  invcenter = invcenter[invcenter$nlowconf < 10,]
  return(invcenter)
}

filter_down_SD_INV_list_nofilter_but_turn <- function(SDi_with_GT){
  
  SDi_simp = SDi_with_GT[SDi_with_GT$n_switches == 1,]
  
  invcenter = SDi_simp %>% group_by(inv_start) %>% slice(1)
  invcenter = invcenter[invcenter$nnoreads < 10,]
  invcenter = invcenter[invcenter$nlowconf < 10,]
  return(invcenter)
}

mark_protective_risk_mixed_invs <- function(invcenter, n1, n2){
  
  invcenter$mix = 'Mixed'
  invcenter[invcenter$pctpos < n1,]$mix = 'SV risk factor'
  invcenter[invcenter$pctpos > n2,]$mix = 'Protective'
  
  # While we are at it, also note the length from SD1 to SD2
  invcenter$SDpairlen = abs(invcenter$source_start - invcenter$target_start)
  
  return(invcenter)
}

make_plot <- function(invcenter, limit){
  p = ggplot(invcenter[(invcenter$inv_end-invcenter$inv_start < limit) & (invcenter$mix != 0),]) +
    geom_boxplot(aes(group=mix, x=mix, y= (nhet+(2*nhom)) / (43*2), fill=mix)) + scale_color_viridis_b() +
    geom_beeswarm(aes(x=mix, y=(nhet+(2*nhom))/(43*2))) + 
    theme_bw() + 
    labs(x='', y='% Inverted alleles', title='Inversions < 200kb')
  return(p)
}


# Do t-test
make_wilcox <- function(df, group1, group2){
  df$inv_alleles = df$nhet + 2*df$nhom
  x = df[df$mix == group1,]$inv_alleles
  y = df[df$mix == group2,]$inv_alleles
  wilcox.test(x,y, alternative='greater', exact=F)
}


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
