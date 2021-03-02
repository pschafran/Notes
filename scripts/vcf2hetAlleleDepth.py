"""
------------------------------------------------------------------------------------------
Description: Converts the variant-calling files (vcf) from ipyrad output to a format required for the use with gbs2ploidy.

Usage: python vcf2hetAlleleDepth.py your_input.vcf
3 optional arguements: -N, -o, and -id. For more information: see bottom

File Name: vcf2hetAlleleDepth.py
Author: Carol Rowe
Date Created: 2018-08-23 using Python 3.6.2

gbs2ploidy reference: Gompert Z, Mock KE (2017) Detection of individual ploidy levels with genotyping-by-sequencing (GBS) analysis. Molecular Ecology Resources 17(6): 1156–1167.)
Gompert and Mock (2017) 
gbs2ploidy aspen data/files: https://datadryad.org//resource/doi:10.5061/dryad.5vs40?show=full

NOTES: 
The format for gbs2ploidy is: two columns for each sample and one row for each SNP.
Homozygous sites and missing data are NA.
Heterozygous sites have counts for each allele.

For example:
NA NA NA NA NA NA 67 1 NA NA
8 6 11 8 NA NA 3 3 NA NA
NA NA 11 8 NA NA NA NA NA NA
NA NA 11 8 NA NA NA NA NA NA
NA NA NA NA NA NA NA NA NA NA
13 1 NA NA NA NA NA NA NA NA
13 1 NA NA NA NA NA NA NA NA

This would represent 5 samples (two columns for each sample, for a total of 10 columns) and 7 SNPs (there are 7 rows).  
The first sample has a heterozygous SNP at the second SNP location (row 2) where there are counts of 8 and 6 for the heterozygous alleles.

gbs2ploidy requires an id file. This file is created for you in this script: "HAD_ID.csv". You have the option to rename this file when you run the script.
The id file is just a list of your sample names (ids). 
Note: The id file used by Gompert and Mock (2017) has the following headers (though, most columns were not used for gbs2ploidy):
pop,id,lab_id,ploidy_updated,prD7Mar16,prT7Mar16,Comments

-----------------------------------------------------------------------------------------
"""

import argparse
import pandas as pd
import copy

__author__ = "Carol Rowe"


def vcf2hetAlleleDepth(filename, skiplines, outfile, id_outfile):	
	# DF = DATAFRAME
	#  READ IN THE FILE AS A DF AND GET JUST THE NEEDED COLUMNS
	# read in the .vcf file from the ipyrad output into a df (dataframe)
	vcf = pd.read_csv(filename, sep="\t", skiprows=skiplines)
	# Slice df: keep "REF" and remaining sample cols
	# .vcf ipyrad files have consistent format, so know which columns to remove
	the_real_data = vcf.drop(vcf.columns[[0,1,2,4,5,6,7,8]], axis=1)

	# GET SAMPLE NAMES TO A LIST
	# get the column names (your sample names) as a list
	the_columns = the_real_data.columns
	# remove the "REF" name from the list
	the_columns = the_columns[1:]
	#print("Your input data contains {} samples." .format(len(the_columns)))

	# CREATE A TEMPORARY DF TO HOLD THE PARSED INFO THAT WE WILL NEED
	# create an empty, temporary df with a single column, "REF", for the reference allele
	temp_df = the_real_data[['REF']].copy()

	# FILL TEMP_DF WITH: CATG_INDEX, AND TWO COLUMNS FOR EACH SAMPLE (SAMPLE_CATG AND SAMPLE_GT): GT = genotype, CATG = counts for each bp
	# The CATG_index column:
	# Function to convert CATG to numeric value for indexing them from a list
	# order is: CATG
	def CATGtoIndex(row):
		if row["REF"] == 'C' :
			return 0
		if row["REF"] == 'A' :
			return 1
		if row["REF"] == 'T' :
			return 2
		if row["REF"] == 'G' :
			return 3
	# now add the CATG number conversion column to the temp_df
	temp_df["CATG_index"] = temp_df.apply(lambda row: CATGtoIndex(row), axis=1)
	# get the sample_CATG and sample_GT columns
	# Slice the CATG counts into a list:
	# i.e. from: 0/0:38:0,0,0,38 get: [0,0,0,38] which is the CATG counts, and 0/0 which is the genotype
	for col in the_columns:
		# each sample will be split into two columns; here we create the names for the cols
		GT = col + "_GT"
		CATG = col + "_CATG"
		temp_df[CATG] = the_real_data[col].str.split(":").str.get(2).str.split(',').tolist()
		temp_df[GT] = the_real_data[col].str.split(":").str.get(0)
	index_range = len(temp_df.index)
	#print("Your input data includes {} SNPs." .format(index_range))

	# REF AND CATG_INDEX COLUMNS TO LISTS:
	ref_list = temp_df['REF'].values.tolist()
	ref_index = temp_df['CATG_index'].values.tolist()

	# CREATE THE FINAL DF WHICH, AT THIS POINT IS EMPTY
	final_df = pd.DataFrame(index=range(index_range))


	# FILL OUR FINAL DF......
	# SPLIT THE TEMP_DF INTO TWO: HETEROZYGOUS AND HOMOZYGOUS
	# RETRIEVE COUNTS FOR THE HETEROZYGOUS. FILL IN NA FOR ALL HOMOZYGOUS
	# PUT HETEROZYGOUS AND HOMOZYGOUS TABELS BACK TOGETHER IN PROPER ORDER
	# AND THEN ADD TO THE FINAL DF
	# Loop through the list of sample names
	for col in the_columns:
		# create an empty, temporary storage df
		df2 = pd.DataFrame(index=range(len(temp_df.index)))
		# retrieving just the columns with our current sample from the temp_df
		df2 = temp_df.filter(regex=col)
		# REF and CATG_index lists to df2
		df2 = df2.assign(CATG_index = ref_index)
		df2 = df2.assign(REF = ref_list)
		#print(col)
		# Get the actual column names so we can reference them easily
		GT_name = col + "_GT" 
		CATG_name = col + "_CATG"   
		# get all the genotypes to a list
		GT_col = df2[GT_name].values.tolist()
		# need just the unique values
		GT_unique = set(GT_col)
		# as a list (GT unique list)
		GT_un_list = list(GT_unique)
		# get just heterozygous genotypes
		# heterozygous sites have: 1/0, 0/1, 2/0, 0/2, etc.
		het_list = []
		for item in GT_un_list:
			if item.split('/')[0] != item.split('/')[1]:
				het_list.append(item)
		# GETTING DF OF JUST HETEROZYGOUS SNPs
		hetz_df = pd.DataFrame(df2.loc[df2[GT_name].isin(het_list)])
		# need to keep the order, i.e. keep the index values.
		# Thus, I will make the index values a new column
		hetz_df['index1'] = hetz_df.index
		# Get the correct CATG value for our alleles
		# First, CATG count values to list
		CATG_values = hetz_df[CATG_name].values.tolist()
		# Make sure these are saved as integers
		CATG_val_int = [[int(j) for j in i] for i in CATG_values]
		# Get the heterozygous GATC index values to list
		CATG_index = hetz_df['CATG_index'].values.tolist()
		# new, empty list to which we will add the "a" column of allele counts
		new_list = []
		# loop for length equal to your number of heterozygous reads (or loci, depending on your terminology)
		for num in range(len(hetz_df.index)):    
			# col_val_int is a list of lists, so get the correct indexing:
			# [num] to get the "outer" list and
			# [CATG_index[num]] to get the correct C,A,T, or G index for the allele count 
			new_list.append( CATG_val_int[num][CATG_index[num]] )
			# now that we saved the needed number, let's delete it from the CATG_val_int list
			del CATG_val_int[num][CATG_index[num]]
			# Now, we will get the maximum value from the remaining list. That should be the 'b' allele
		b = [max(item) for item in CATG_val_int]
	# New column names for the "a" and "b" allele counts
		col_a = col + "_a"
		col_b = col + "_b"
		# Add the allele counts to the heterozygous df
		hetz_df[col_a] = new_list
		hetz_df[col_b] = b
		# Grab just columns we want:
		final_hetz_df = hetz_df[['index1', col_a, col_b]]
		# GETTING JUST THE HOMOZYGOUS DF    
		# Now, let's get the homozygous alleles by getting the index numbers not in the het.
		final_hetz_indeces = final_hetz_df['index1'].values.tolist()
		total_indeces = list(df2.index)
		final_homz_indeces = list(set(total_indeces) - set(final_hetz_indeces))
		# create the empty homozygous df
		final_homz_df = pd.DataFrame(index=range(len(final_homz_indeces)))
		# add columns to the homozygous df: index, col_a, and col_b 
		final_homz_df['index1'] = final_homz_indeces
		# col_a and col_b (allele counts) all have NA
		final_homz_df[col_a] = 'NA'
		final_homz_df[col_b] = 'NA'

		# Now stack the two df on top of each other
		sub_final_df = pd.concat([final_homz_df, final_hetz_df], axis=0)
		# Sort in the index1 column
		sub_final_df.sort_values(by=['index1'], ascending=True, inplace=True)
		# have to reset the index (makes life easier when adding to the final_df)
		sub_final_df.reset_index(drop=True, inplace=True)
		# add the allele count columns to the final_df
		final_df[[col_a, col_b]] = sub_final_df[[col_a, col_b]]

	# Save your final df to a file:
	# no headers, no index, and space as a separator
	final_df.to_csv(outfile, header=False, index=False, sep=" ")

	Count_Row=final_df.shape[0] #gives number of row count
	Count_Col=int((final_df.shape[1])/2) #gives number of col count
	print("Your final data file contains {} samples and {} SNPs." .format(Count_Col, Count_Row))


	# And now for the csv file which will just be a column of your sample names.
	id_index_range = len(the_columns)
	id_df = pd.DataFrame(index=range(id_index_range))
	id_df['id'] = the_columns
	# save with header, with default ',' separator:
	id_df.to_csv(id_outfile, index=False)


# if name in main so we can run vcf2hetAlleleDepth.py by itself (main)
# or it can still be used if you use it embedded (import vcf2hetAlleleDepth) within another script (name)
if __name__ in '__main__':
	# This allows the --help to show the docstring at the top of this script.
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	# Add arguments (4):
	# First argument is mandatory
	parser.add_argument('input', metavar='vcf_input_file', help='Enter the name of your vcf file.')
	# Next 3 arguments are optional
	parser.add_argument('-N','--skiplines', help='Number of lines to skip in vcf file until header of the data table. Default = 10', type=int, default=10, required=False)
	parser.add_argument('-o','--outfile', help='Output heterozygous allele counts file name. Default = hetAlleleDepth.txt', type=str, default='hetAlleleDepth.txt', required=False)
	parser.add_argument('-id','--id_outfile', help='Output file with list of samples in your vcf file. Default = HAD_ID.csv', type=str, default='HAD_ID.csv', required=False)
	# Array for all arguments passed to script:
	args = parser.parse_args()

	# Now, we can access the arguments input by the user (or use defaults), and apply to our function	
	vcf2hetAlleleDepth(args.input, args.skiplines, args.outfile, args.id_outfile)

