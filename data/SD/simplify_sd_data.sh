#!/bin/bash
set -euo pipefail

cat <(echo "chr_source\tstart_source\tend_source\torientation\tchr_target\tstart_target\tend_target") <(cut -f 1,2,3,7,8,9,6 GRCh38GenomicSuperDup.bed) > SDs_simplified.bed

#1/5: same chr. 
#2/6: different coordinates
awk '{FS=OFS="\t"} (($1==$5) && ($2 != $6)) {print $0}' SDs_simplified.bed > SDs_samechr.bed

awk '{FS=OFS="\t"} {print $0, int((NR+1)/2)}' SDs_samechr.bed > SDs_samechr_number.bed 
bedtools intersect -a SDs_samechr_number.bed -b ../INV/invs_final_inner.bed -wao -f 1 > SDs_with_inv.bed
