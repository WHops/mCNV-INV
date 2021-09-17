# Whoeps, 17th Sep 2021

# Libraries and functions #

library(dplyr)
library(ggplot2)

cut_down_len <- function(SD_f, params_f){
  SD_f = SD_f %>% group_by(pairID) %>% 
    mutate(valid = all(source_end - source_start > params_f$min_sd_len))
  SD_f = SD_f[SD_f$valid == T,]
  return(SD_f)
}

add_sd_stats <- function(SD_f){
  SD_f$sdlen = SD_f$source_end - SD_f$source_start
  SD_f$pairlen = abs(SD_f$source_end - SD_f$target_end)
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

######### Set parameters ##########################################
params = list()
params$min_sd_len = 10000 # 10kb
params$max_pairlen = 10 * 1000 * 1000 # 10Mb. 

SD_link = '../data/SD/SDs_with_inv.bed'

######### GO ######################################################

# Load data table
SD = read.table(SD_link)
colnames(SD) = c('source_chr', 'source_start', 'source_end', 'orientation',
                'target_chr','target_start', 'target_end', 'pairID',
                'inv_chr', 'inv_start', 'inv_end', 'overlap_size')

# Cut down to min len
SD_cut = cut_down_len(SD, params)

# Add some information on sdlen, pairlen, invlen
SD_cut_stats = add_sd_stats(SD_cut)

# This is dubious. I want to exclude CNVs > 10Mb. Because I feel like they don't happen.
SD_cut_stats_pairlen = SD_cut_stats[SD_cut_stats$pairlen < params$max_pairlen,]

# Only SDs where an inv has a potential for orientation change
SDi = find_SDs_with_inv_interference_potential(SD_cut_stats_pairlen)

# Made a bedfile for igv
df_bed = prep_df_bed(SDi)
# Make another bedfile for igv^^
SDi_bed = SDi[,c('source_chr', 'source_start', 'source_end')]

# Save ##############################################################
write.table(file = "~/Desktop/df_bed.bed", na.omit(df_bed), col.names = F, row.names = F, quote=F, sep='\t')
write.table(file = "~/Desktop/sdi_bed.bed", na.omit(SDi_bed), col.names = F, row.names = F, quote=F, sep='\t')



