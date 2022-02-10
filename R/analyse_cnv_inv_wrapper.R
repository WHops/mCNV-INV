# Whoeps, 17th Sep 2021

#' Small helperfunction for the wrapper to 
#' source the functions file. 
#' Taken from Stackoverflow:
#' https://stackoverflow.com/questions/42815889/r-source-and-path-to-source-files/42817150#42817150
#' @author MrFlick, Stackoverflow
source_here <- function(x, ...) {
  dir <- "."
  if(sys.nframe()>0) {
    frame <- sys.frame(1)
    if (!is.null(frame$ofile)) {
      dir <- dirname(frame$ofile)
    }
  }
  source(file.path(dir, x), ...)
}

# Libraries and functions #

library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(matrixStats)
library(stringr)

source_here('R/analyse_cnv_inv_functions.R')


######### Set parameters ##########################################
params = list()
params$min_sd_len = 10000 # Consider only SDs longer than this..
params$max_pairlen = 1e10 # Unlimited
params$sds_protect_ro_risk_th = 0.1 # Fraction of SDs flipped in SAME direction to be considered non-mixed.
params$use_only_high_confidence_invs = T # Use only high confidence inversions?
params$save = T # Save output files

######### Set input files ##########################################

# Inversion data: position, 
gt_link = "data/INV/inversions_genotypes.tsv"
inv_noSD_link = 'data/INV/inversions.tsv'
recurrence_link = 'data/INV/recurrence.tsv'

SD_link = 'data/SD/SDs_with_inv.bed'
mcnv_link = 'data/mCNV/mcnvs.tsv'

outdir = 'res/'

######### Run analysis ######################################################


SDi = part1_find_all_SDpairs_affected_by_inv(SD_link, params, outdir)
part2_find_protective_and_risk_inversions(SDi, params, outdir)

######### Inform human ######################################################
print(paste0('Done! Output files written to: ', outdir))

