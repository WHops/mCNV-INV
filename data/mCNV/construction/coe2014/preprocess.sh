
# variants_for_nstd_100.csv downloaded from: https://www.ncbi.nlm.nih.gov/dbvar/studies/nstd100/download/?type=v

# Preprocessing: 
cat <(echo "ID\tType\tchr\tstart\tend") <(awk '{FS=","; OFS="\t"} ((NR>1) && ($21="remapped")) {print $2,$3,$14,$17,$18}' variants_for_nstd100.csv | sed 's/"//g') > coe_cnvs.txt

