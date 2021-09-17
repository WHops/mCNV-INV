awk '{FS=OFS="\t"} ($1==$5)' SDs_simplified.bed > SDs_samechr.bed

