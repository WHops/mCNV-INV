
# decipher.txt was taken from https://www.deciphergenomics.org/disorders/syndromes/list

# Then, reformatting started here: 
cat decipher.txt | tr '\n' '\t' | sed 's/Mb/Mb\n/g' | sed 's/kb/kb\n/g' > decipher_formatted.txt

# Last, go into vim, first columns needs a bit of polishing
