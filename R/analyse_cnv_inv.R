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


prep_SD_oneswitch <- function(SDi_with_GT, n1, n2){
  
  SDi_simp = SDi_with_GT[SDi_with_GT$n_switches == 1,]
  
  # SDi_simp = SDi_with_GT_s1 %>% group_by(inv_start) %>% mutate(inv_n_alterations = length(inv_chr),
  #                                                                  allp = all(orientation=='+'),
  #                                                                  alln = all(orientation=='-'),
  #                                                                  pctpos = sum((abs(target_end - source_start))[orientation=='+']) / sum(abs(target_end-source_start)),
  #                                                                  mix=any(allp) + 2*any(alln))
  invcenter = SDi_simp %>% group_by(inv_start) %>% slice(1)
  invcenter = invcenter[invcenter$nnoreads < 10,]
  invcenter = invcenter[invcenter$nlowconf < 10,]
  
  
  invcenter$catpos = 'c Mixed'
  invcenter[invcenter$pctpos < n1,]$catpos = 'b SV risk factor'
  invcenter[invcenter$pctpos > n2,]$catpos = 'a Protective'
  invcenter$mix = invcenter$catpos
  
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
# write.table(file = "~/Desktop/df_bed.bed", na.omit(df_bed), col.names = F, row.names = F, quote=F, sep='\t')
# write.table(file = "~/Desktop/sdi_bed.bed", na.omit(SDi_bed), col.names = F, row.names = F, quote=F, sep='\t')


# Part 2. Find 'foci'

SDi_touch_inv = SDi[SDi$inv_chr != '.',]
#SDi_plus = SDi_touch_inv[SDi_touch_inv$source_start < SDi_touch_inv$target_start,]
SDi_regions = SDi_touch_inv[,c('source_chr','source_start', 'target_end', 'orientation', 'inv_chr', 'inv_start', 'inv_end', 'pairID')]

SDi_chromcenter = SDi_regions %>% group_by(inv_start) %>% mutate(inv_n_alterations = length(inv_chr),
                                                                 allp = all(orientation=='+'),
                                                                 alln = all(orientation=='-'),
                                                                 pctpos = sum((abs(target_end - source_start))[orientation=='+']) / sum(abs(target_end-source_start)),
                                                                 mix=any(allp) + 2*any(alln))

# good, but we need to think of double switches. 

SDi_chromcenter_nswitch = SDi_chromcenter %>% group_by(pairID) %>% mutate(n_switches = length(inv_chr))

gt_link = "/Users/hoeps/PhD/projects/huminvs/analyses_paper/inv_gts/polished_table.tsv"
gts = read.table(gt_link, sep='\t', header=T)
gts = gts[, c('seqnames','start','end','nref','nhet','nhom','ncomplex','ninvdup','nnoreads','nlowconf')]
colnames(gts)[1:3] = c('inv_chr','inv_start','inv_end')

SDi_with_GT = left_join(SDi_chromcenter_nswitch, gts, by=c('inv_chr','inv_start','inv_end'))

invcenter2 = prep_SD_oneswitch(SDi_with_GT, 0.2, 0.8)


# Combine with more information
SD_noSD_link = '/Users/hoeps/PhD/projects/huminvs/analyses_paper/inv_refinement/mastertable_jul1/for_overlap/sorted/4_horsemen/tailored/lab/all_sorted.txt'
sdnosd = read.table(SD_noSD_link, header=F, sep='\t')
colnames(sdnosd) = c('inv_chr' ,'inv_start'  ,'inv_end','class')
invcenter = left_join(invcenter2, sdnosd, by=c('inv_chr' ,'inv_start'  ,'inv_end'))
invcenter$SDpairlen = abs(invcenter$source_start - invcenter$target_end)
#invcenter = invcenter[invcenter$class != 'invs_unprocessed.bed',]


limit = 1e10 #200000
g = make_plot(invcenter, limit=limit)
g
make_wilcox(invcenter[invcenter$inv_end-invcenter$inv_start < limit,], 'a Protective', 'b SV risk factor')
make_wilcox(invcenter[invcenter$inv_end-invcenter$inv_start < limit,], 'a Protective', 'c Mixed')
make_wilcox(invcenter[invcenter$inv_end-invcenter$inv_start < limit,], 'c Mixed', 'b SV risk factor')

ablub = (invcenter[invcenter$inv_end-invcenter$inv_start < limit,])
print(blub$class)


ggplot(invcenter) + geom_histogram(aes(x=log10(abs(inv_end - inv_start)), fill=mix), breaks=5)

# Now off to new adventures.
# SDi_with_GT_s2 = SDi_with_GT[SDi_with_GT$n_switches > 1,]
# 
# SDi_simp = SDi_with_GT_s2 %>% group_by(inv_start) %>% mutate(inv_n_alterations = length(inv_chr),
#                                                              allp = all(orientation=='+'),
#                                                              alln = all(orientation=='-'),
#                                                              pctpos = sum((abs(target_end - source_start))[orientation=='+']) / sum(abs(target_end-source_start)),
#                                                              mix=any(allp) + 2*any(alln))
# invcenter = SDi_simp %>% group_by(inv_start) %>% slice(1)

# 
# invcenter$catpos = 'c Mixed'
# invcenter[invcenter$pctpos < 0.1,]$catpos = 'b SV risk factor'
# invcenter[invcenter$pctpos > 0.9,]$catpos = 'a Protective'
# invcenter$mix = invcenter$catpos
# 
# g = make_plot(invcenter, limit=limit)
# g
# 
# blub2 = (invcenter)


#for (lim in c(100000,125000,150000, 500000)){
#}
# ggplot(invcenter) + geom_beeswarm(aes(x=log10(inv_end - inv_start), y=mix))  
# blub = blub[blub$nlowconf < 5,]
# 
# ggplot(blub) + 
#   geom_beeswarm(aes(x=pctpos, y=ncomplex, color=log10(inv_end-inv_start))) + scale_color_viridis_b()
