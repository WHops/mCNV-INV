# Whoeps, 17th Sep 2021

# Libraries and functions #

library(dplyr)
library(ggplot2)
library(ggbeeswarm)
source('./analyse_cnv_inv_functions.R')


######### Set parameters ##########################################
params = list()
params$min_sd_len = 10000 # 10kb
params$max_pairlen = 1e10 # Unlimited
params$sds_protect_ro_risk_th = 0.1#000000000001#0.00000000001#0.10 #1e-10 # 0 + tiny
params$save = T


SD_link = '../data/SD/SDs_with_inv.bed'
gt_link = "../data/INV/polished_GT_table.tsv"
SD_noSD_link = '../data/INV/inversion_class_refinement.txt'
mcnv_link = '../data/mCNV/as_bed/decipher.bed'

outdir = '../res/'
######### GO ######################################################


#' part1_find_all_SDpairs_affected_by_inv
#' @description # Part 1: find all SD pairs which are affected by an inversion.
#' This produces two tables that can be viewed in IGV
#' @param
#' @return
#'
#' @author Wolfram Höps
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
  SD_minlen_stats_pairlen = SD_minlen_stats[SD_minlen_stats$pairlen < params$max_pairlen, ]
  
  # Keep only SDs where an inv has a potential for orientation change.
  # This means one of them, but not both, is overlapped by the same inversion.
  # These are the 'interesting' SDs (SDi).
  SDi = find_SDs_with_inv_interference_potential(SD_minlen_stats_pairlen)
  
  # Filter to only SDs being touched by Inversions (As opposed to keeping the pairs.)
  SDi_touch_inv = SDi[SDi$inv_chr != '.', ]
  
  # Made two bedfiles for IGV.
  df_bed = prep_df_bed(SDi)
  SDi_bed = SDi[, c('source_chr', 'source_start', 'source_end')]
  
  # Save results
  if (params$save) {
    dir.create(outdir, showWarnings = F)
    write.table(
      SDi,
      file = paste0(outdir, '/SDi.tsv'),
      col.names = T,
      row.names = F,
      quote = F,
      sep = '\t'
    )
    write.table(
      file = paste0(outdir, "bedfile_1.bed"),
      na.omit(df_bed),
      col.names = F,
      row.names = F,
      quote = F,
      sep = '\t'
    )
    write.table(
      file = paste0("bedfile_2.bed"),
      na.omit(SDi_bed),
      col.names = F,
      row.names = F,
      quote = F,
      sep = '\t'
    )
  }
  
  return(SDi_touch_inv)
  
}


SDi_touch_inv = part1_find_all_SDpairs_affected_by_inv(SD_link, params, outdir)


#' part1_find_all_SDpairs_affected_by_inv
#' @description # Part 2: Find 'protective' and 'risk' Inversions:
#' Those are inversions that lead to creation of only inverted or only direct SDs.
#' I.e. non-messy sites :)
#' This part is in general a bit more experimental than part 1.
#' @param
#' @return
#'
#' @author Wolfram Höps
#' @export
part2_find_protective_and_risk_inversions <-
  function(SDi_touch_inv) {
    cols_to_keep = c(
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
    SDi_regions = SDi_touch_inv[, cols_to_keep]
    
    # Inversion-wise statistics: are all affected SD pairs in pos or neg directions?
    # How many SD pairs does this inversion affect?
    
    SDi_chromcenter = determine_directionality_of_affected_SDs_below_inv(SDi_regions)

    # Another caveat: There are SD pairs that are affected by >1 inversions.
    # The 'nswitch' parameter explains per SD pair, how many inversions can affect its orientation.
    SDi_chromcenter_nswitch = SDi_chromcenter %>% group_by(pairID) %>% mutate(n_switches = length(inv_chr))
    
    # Bring GTs into the game.
    gts = read.table(gt_link, sep = '\t', header = T)
    gts = gts[, c(
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
    )]
    colnames(gts)[1:3] = c('inv_chr', 'inv_start', 'inv_end')
    SDi_with_GT = left_join(SDi_chromcenter_nswitch,
                            gts,
                            by = c('inv_chr', 'inv_start', 'inv_end'))
    
    # # For the record
    # length(unique(SDi_with_GT[SDi_with_GT$n_switches == 1,]))
    
    # This is a big step here. We want to filter our list of affected SDs by:
    # - SD pair should only be affected by ONE inversion (otherwise interpretation gets tricky)
    # - Inversion should have <10 GTs with noreads or lowconf reads.
    SDi_with_GT_goodsamples = filter_down_SD_INV_list(SDi_with_GT)
    
    # TURN THIS ON IF YOU NEED the full table!
    SDi_with_GT_goodsamples =  SDi_with_GT %>% group_by(inv_chr, inv_start, inv_end) %>% slice(1)
    # Now, mark which inversions are likely protective, risky and mixed.
    # We consider 'protective' if  < sds_protect_ro_risk_th of SD-pair-basepairs switch from protect to risk
    # We consider 'risky' if > 1-x% of SD-pair-basepairs switch from protect to risk
    invcenter = mark_protective_risk_mixed_invs(
      SDi_with_GT_goodsamples,
      params$sds_protect_ro_risk_th,
      1 - params$sds_protect_ro_risk_th
    )
    
    #write.table(invcenter, file='../final/invcenter.txt', row.names=F, col.names=T, sep='\t', quote = F)
    # More info still. We also want to know the inversion class. SD-mediated? nonSD-mediated?
    sdnosd = read.table(SD_noSD_link, header = F, sep = '\t')
    colnames(sdnosd) = c('inv_chr' , 'inv_start'  , 'inv_end', 'class')
    invcenter = left_join(invcenter, sdnosd, by = c('inv_chr' , 'inv_start'  , 'inv_end'))
    
    # Add mCNVs
    mcnv_list = read.table(mcnv_link, header = F, sep = '\t')
    colnames(mcnv_list) = c('chr',
                            'start',
                            'end',
                            'mCNV name',
                            'genome',
                            'gt',
                            'del/dup',
                            'len')
    mcnv_list$chr = paste0('chr', mcnv_list$chr)
    
    # Bedtools
    invcenter = invcenter[moveme(names(invcenter),
                                 "inv_chr first; inv_start after inv_chr; inv_end after inv_start")]
    invcenter_protect_risk = invcenter[invcenter$mix != 'Mixed', ]
    
    # Do mCNVs better...
    SDs_flipped = left_join(SDi,
                            invcenter_protect_risk,
                            by = c('inv_chr', 'inv_start', 'inv_end'))
    SDs_flipped2 = SDs_flipped[!is.na(SDs_flipped$SDpairlen), ]
    # Add recurrence
    recurrence_link = '/Users/hoeps/PhD/projects/huminvs/limix-inversion-project/limixploteR/data_new/recurrence_hufsah.tsv'
    rec = read.table(recurrence_link, header = T, sep = '\t')
    colnames(rec)[1:3] = c('inv_chr', 'inv_start', 'inv_end')
    rec2 = rec[, c(
      'inv_chr',
      'inv_start',
      'inv_end',
      'verdictRecurrence_hufsah',
      'verdictRecurrence_benson'
    )]
    
    invcenter_protect_risk = left_join(invcenter_protect_risk,
                                       rec2,
                                       by = c('inv_chr', 'inv_start', 'inv_end'))
    
    # SAVE
    write.table(
      invcenter_protect_risk[order(invcenter_protect_risk$inv_chr), ],
      file = '../final/protect_risk_loci.tsv',
      sep = '\t',
      quote = F,
      col.names = T,
      row.names = F
    )
    
    
    overlap_main = bedtoolsr::bt.intersect(invcenter_protect_risk, mcnv_list, wao =
                                             T)
    colnames(overlap_main) = c(colnames(invcenter_protect_risk), colnames(mcnv_list))
    head(overlap_main)
    
    library(matrixStats)
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
    head(overlap_sec)
    
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
    overlap_save2 = overlap_save[(overlap_save$`mCNV name` != '.'), ]
    write.table(
      overlap_save2,
      file = '../final/flipped_SD_pairs_overlapping_mCNVs.txt',
      row.names = F,
      col.names = T,
      sep = '\t',
      quote = F
    )
    overlap_main2 = left_join(overlap_main, rec2, by = c('chr', 'start', 'end'))
    
    
    # To_save
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
    overlap_save = overlap_main[, cols_to_save]
    overlap_save = overlap_save[overlap_save$`mCNV name` != '.', ]
    colnames(overlap_save)[colnames(overlap_save) == 'orientation'] = 'SDpairs_orientation'
    colnames(overlap_save)[colnames(overlap_save) == 'inv_n_alterations'] = 'SDpairs_flipped'
    colnames(overlap_save)[colnames(overlap_save) == 'class'] = 'Inversion_class'
    colnames(overlap_save)[colnames(overlap_save) == 'mix'] = 'inv_role'
    
    write.table(
      overlap_save,
      file = 'invs_affecting_sds_affecting_mcnvs_strict.tsv',
      row.names = F,
      col.names = T,
      quote = F,
      sep = '\t'
    )
    ### Now
    # Made a bedfile for igv
    invcenter_bed = prep_df_bed2(invcenter[invcenter$mix != 'Mixed', ])
    write.table(
      invcenter_bed,
      file = '~/Desktop/protect_risk.bed',
      sep = '\t',
      col.names = F,
      row.names = F,
      quote = F
    )
    # Make another bedfile for igv^^
    SDi_bed = SDi[, c('source_chr', 'source_start', 'source_end')]
    
  }
# Make plots and stat.tests #

# Plot length of all interesting SD pairs
ggplot(data = SDi[(SDi$target_start > SDi$source_start), ]) +
  geom_histogram(aes(x = abs(target_start - source_start), fill = orientation)) +
  scale_x_log10() +
  theme_bw()

# plot length of inversting SD pairs being flipped by one inversion only, in the same directions
ggplot(data = invcenter[invcenter$mix != 'Mixed', ]) +
  geom_histogram(aes(x = SDpairlen, fill = mix)) +
  scale_x_log10() +
  theme_bw()


# Max size of inversion to consider in plot
limit = 200000

# Plot that invcenter plot for genotypes
g = make_plot(invcenter, limit = limit)
g

# Make tests
make_wilcox(invcenter[invcenter$inv_end - invcenter$inv_start < limit, ], 'Protective', 'SV risk factor')
make_wilcox(invcenter[invcenter$inv_end - invcenter$inv_start < limit, ], 'Protective', 'Mixed')
make_wilcox(invcenter[invcenter$inv_end - invcenter$inv_start < limit, ], 'Mixed', 'SV risk factor')



### NUMBERS REPORTED IN THE PAPER
n_SDs = length(unique(SDi$pairID))
n_invs = length(unique(paste0(SDi$inv_end, SDi$inv_start)))

SDi2 = SDi[SDi$inv_start != '-1', ]
n_invs - sum(table(paste0(SDi2$inv_start, SDi2$inv_end)) == 1)
max(table(paste0(SDi2$inv_start, SDi2$inv_end)))


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
# invcenter$catpos = 'Mixed'
# invcenter[invcenter$pctpos < 0.1,]$catpos = 'SV risk factor'
# invcenter[invcenter$pctpos > 0.9,]$catpos = 'Protective'
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
