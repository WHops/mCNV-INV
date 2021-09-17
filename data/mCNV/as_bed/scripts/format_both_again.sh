# Decipher to bed
gawk '{FS=OFS="\t"} (match($3,/([0-9XY]*)?:([0-9]*)-([0-9]*)/, a)){print a[1], a[2],a[3],$1, $4,$5,$6,$7}' decipher_formatted.bed > decipher.bed

# Coe to bed
awk '{FS=OFS="\t"} {print $3,$4,$5,$1,$2}' coe_cnvs.bed > coe.bed
