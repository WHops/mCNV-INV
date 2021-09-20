# Whoeps, 17th Sep 2021

# Libraries and functions #

library(dplyr)
library(ggplot2)

source('./analyse_cnv_inv_functions.R')


######### Set parameters ##########################################
params = list()
params$min_sd_len = 10000 # 10kb
params$max_pairlen = 10 * 1000 * 1000 # 10Mb. 

SD_link = '../data/SD/SDs_with_inv.bed'
gt_link = "../data/INV/polished_GT_table.tsv"
SD_noSD_link = '../data/INV/inversion_class_refinement.txt'

######### GO ######################################################

# Part 1: find all SD pairs which are affected by an inversion. 
# This produces two tables that can be viewed in IGV

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

# Save #
write.table(file = "~/Desktop/df_bed.bed", na.omit(df_bed), col.names = F, row.names = F, quote=F, sep='\t')
write.table(file = "~/Desktop/sdi_bed.bed", na.omit(SDi_bed), col.names = F, row.names = F, quote=F, sep='\t')

#############################################################

# Part 2: Find 'protective' and 'risk' Inversions: 
# Those are inversions that lead to creation of only inverted or only direct SDs. 
# I.e. non-messy sites :) 
# This part is in general a bit more experimental than part 1.

# Filter to only SDs being touched by Inversions
SDi_touch_inv = SDi[SDi$inv_chr != '.',]

#SDi_plus = SDi_touch_inv[SDi_touch_inv$source_start < SDi_touch_inv$target_start,]
SDi_regions = SDi_touch_inv[,c('source_chr','source_start', 'target_end', 'orientation', 'inv_chr', 'inv_start', 'inv_end', 'pairID')]

# Inversion-wise statistics: are all affected SD pairs in pos or neg directions? 
# How many SD pairs does this inversion affect?
SDi_chromcenter = SDi_regions %>% group_by(inv_start) %>% mutate(inv_n_alterations = length(inv_chr),
                                                                 allp = all(orientation=='+'),
                                                                 alln = all(orientation=='-'),
                                                                 pctpos = sum((abs(target_end - source_start))[orientation=='+']) / sum(abs(target_end-source_start)),
                                                                 mix=any(allp) + 2*any(alln))

# Another caveat: There are SD pairs that are affected by >1 inversions. 
# The 'nswitch' parameter explains per SD pair, how many inversions can affect its orientation. 
SDi_chromcenter_nswitch = SDi_chromcenter %>% group_by(pairID) %>% mutate(n_switches = length(inv_chr))

# Bring GTs into the game. 
gts = read.table(gt_link, sep='\t', header=T)
gts = gts[, c('seqnames','start','end','nref','nhet','nhom','ncomplex','ninvdup','nnoreads','nlowconf')]
colnames(gts)[1:3] = c('inv_chr','inv_start','inv_end')
SDi_with_GT = left_join(SDi_chromcenter_nswitch, gts, by=c('inv_chr','inv_start','inv_end'))

# This is a big step here. We want to filter our list of affected SDs by:
# - SD pair should only be affected by ONE inversion (otherwise interpretation gets tricky)
# - Inversion should have <10 GTs with noreads or lowconf reads.
SDi_with_GT_goodsamples = filter_down_SD_INV_list(SDi_with_GT,)

# Now, mark which inversions are likely protective, risky and mixed. 
# We consider 'protective' if <20% of SD-pair-basepairs switch from protect to risk
# We consider 'risky' if >80% of SD-pair-basepairs switch from protect to risk
invcenter = mark_protective_risk_mixed_invs(SDi_with_GT_goodsamples, 0.2,0.8)


# More info still. We also want to know the inversion class. SD-mediated? nonSD-mediated?
sdnosd = read.table(SD_noSD_link, header=F, sep='\t')
colnames(sdnosd) = c('inv_chr' ,'inv_start'  ,'inv_end','class')
invcenter = left_join(invcenter2, sdnosd, by=c('inv_chr' ,'inv_start'  ,'inv_end'))
#invcenter = invcenter[invcenter$class != 'invs_unprocessed.bed',]

# Make plots and stat.tests
limit = 200000
g = make_plot(invcenter, limit=limit)
g
make_wilcox(invcenter[invcenter$inv_end-invcenter$inv_start < limit,], 'a Protective', 'b SV risk factor')
make_wilcox(invcenter[invcenter$inv_end-invcenter$inv_start < limit,], 'a Protective', 'c Mixed')
make_wilcox(invcenter[invcenter$inv_end-invcenter$inv_start < limit,], 'c Mixed', 'b SV risk factor')

ablub = (invcenter[invcenter$inv_end-invcenter$inv_start < limit,])
print(blub$class)


ggplot(invcenter) + geom_histogram(aes(x=log10(abs(inv_end - inv_start)), fill=mix), bins=5)



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
