
# mcnv-inv

## Overview of this workflow

mcnv-inv is a repo collecting functions and wrappers for the systematic analysis of the overlap interplays between morbid copy-number variants (mCNVs), polymorphic inversions (INVs) and Segmental Duplications (SDs). 

Scripts here were used for analyses reported in Porubsky et al., 2021 (https://doi.org/10.1101/2021.12.20.472354)


## Running the script


### Input Data 
Example input files are provided in data/

Input tables are to be stored in data/ and its subfolders. Required files:


	data/INV/inversions.tsv	# Table containing inversion coordinates [coord, start, end, SV_class]
	data/INV/inversions_genotypes.tsv # Table containing inversion genotypes 
	data/INV/recurrence.tsv # Table containing recurrence information	
	data/SD/SDs_with_inv.bed # Table containing SDs, including columns defining their overlapping inversions
	data/mCNV/mcnvs.tsv		# Table containing morbid CNVs, in table format obtained from decipher database 



### Instructions	

Optional: you can adjust parameters (such as minimal/maximal SD size) in the header of *analyse_cnv_inv_wrapper.R*.

Run the wrapper script using: 
 `Rscript R/analyse_cnv_inv_wrapper.R `


### Output

Output files are by default, stored in *res/*. 

The main wrapper script will produce five output tables, named after sheets A-E of table S11, reported in Porubsky et al. 2021.

	sheetA_flipping_inversions.tsv	Inversions affecting SD pairs
	sheetB_protect_risk_loci.tsv	Inversions 'flipping' SD-pairs in the same orientation
	sheetC_invs_affecting_sds_affecting_mcnvs_strict.tsv	Inversions 'flipping' SD-pairs overlapping with a morbid CNV. 
	sheetD_flipped_SD_pairs_overlapping_mCNVs.txt SD-pairs overlapping a morbid CNV and flipped by n=1 Inversion
	sheetE_SDi.tsv All SD pairs flipped by any inversions. Every SD of every pair is noted in this table, with a total of 1094 pairs affected by an inversion. 




## References

This code is associated to our publication on human inversion recurrence:
Porubsky et al., 2021 (https://doi.org/10.1101/2021.12.20.472354)

