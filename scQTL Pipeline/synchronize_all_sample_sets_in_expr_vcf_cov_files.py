import numpy as np
import sys
import os
import random
from operator import itemgetter, attrgetter, methodcaller
# import matplotlib.pyplot as plt
# import os.path

'''
USAGE:
python synchronize_all_sample_sets_in_expr_vcf_cov_files.py o.vcf OPC__cmc_kellis.bed cmc_kellis_covariates_v1.bed output_vcf_file.vcf output_expr_file.bed output_cov_file.bed

'''


input_vcf_file_name = sys.argv[1]
input_expr_file_name = sys.argv[2]
input_cov_file_name = sys.argv[3]
output_vcf_file_name = sys.argv[4]
output_expr_file_name = sys.argv[5]
output_cov_file_name = sys.argv[6]


####  First build list_of_common_samples_in_vcf_expr_cov_files
list_of_INTERSECTED_sampleIDs_in_vcf_expr_cov_files = list()
list_of_sampleIDs_in_vcf_file = list()
input_vcf_file = open(input_vcf_file_name, "r")
for line in input_vcf_file:
	if "CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT" in line: ## ie, at header
		ln_elms = list()
		ln_elms = line.split("\t")
		col_index = 9
		while col_index < len(ln_elms):
			sampleID = ln_elms[col_index].strip()
			list_of_sampleIDs_in_vcf_file.append(sampleID)
			col_index += 1
input_vcf_file.close()
'''
for s in list_of_sampleIDs_in_vcf_file:
	print "|" + s + "|"
'''


list_of_sampleIDs_in_expr_file = list()
input_expr_file = open(input_expr_file_name, "r")
for line in input_expr_file:
	if "chr	start	end	gene	length	strand" in line: # ie, at header
		ln_elms = list()
		ln_elms = line.split("\t")
		col_index = 6
		while col_index < len(ln_elms):
			sampleID = ln_elms[col_index].strip()
			list_of_sampleIDs_in_expr_file.append(sampleID)
			col_index += 1
input_expr_file.close()
'''
for s in list_of_sampleIDs_in_expr_file:
	#print "|" + s + "|"
'''


list_of_sampleIDs_in_cov_file = list()
input_cov_file = open(input_cov_file_name, "r")
for line in input_cov_file:
	ln_elms = list()
	ln_elms = line.split("\t")
	if ln_elms[0] == "SampleID": # ie, at header
		col_index = 1
		while col_index < len(ln_elms):
			sampleID = ln_elms[col_index].strip()
			list_of_sampleIDs_in_cov_file.append(sampleID)
			col_index += 1
input_cov_file.close()
'''
for s in list_of_sampleIDs_in_expr_file:
	#print "|" + s + "|"
'''


for sampleID in list_of_sampleIDs_in_vcf_file:
	if sampleID in list_of_sampleIDs_in_expr_file:
		if sampleID in list_of_sampleIDs_in_cov_file:
			if sampleID not in list_of_INTERSECTED_sampleIDs_in_vcf_expr_cov_files:
				list_of_INTERSECTED_sampleIDs_in_vcf_expr_cov_files.append(sampleID)





####  Work on VCF file
### build dict: vcf_sample_ID___2___col_index
vcf_sample_ID___2___col_index = {}
input_vcf_file = open(input_vcf_file_name, "r")
for line in input_vcf_file:
	if "CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT" in line: ## ie, at header
		ln_elms = list()
		ln_elms = line.split("\t")
		col_index = 9
		while col_index < len(ln_elms):
			sampleID = ln_elms[col_index].strip()
			vcf_sample_ID___2___col_index[sampleID] = col_index
			col_index += 1
input_vcf_file.close()


### print output vcf file
output_vcf_file = open(output_vcf_file_name, "w")

# (print header material)
output_vcf_file.write("##fileformat=VCFv4.3" + "\n")
header_line_to_print = "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT"
for sampleID in list_of_INTERSECTED_sampleIDs_in_vcf_expr_cov_files:
	header_line_to_print = header_line_to_print + "\t" + sampleID 
output_vcf_file.write(header_line_to_print + "\n")

# (go through each line and print out relevant columns)
input_vcf_file = open(input_vcf_file_name, "r")
list_of_variant_IDs = list()
for line in input_vcf_file:
	if "#" not in line:  ## we've reached a variant (ie, non-header) line
		ln_elms = list()
		ln_elms = line.split("\t")
		variant_ID = ln_elms[2]
		if variant_ID not in list_of_variant_IDs:
			list_of_variant_IDs.append(variant_ID)
			line_to_print = ln_elms[0] + "\t" + ln_elms[1] + "\t" + ln_elms[2] + "\t" + ln_elms[3] + "\t" + ln_elms[4] + "\t" + ln_elms[5] + "\t" + ln_elms[6] + "\t" + ln_elms[7] + "\t" + ln_elms[8]
			for sampleID in list_of_INTERSECTED_sampleIDs_in_vcf_expr_cov_files:
				rel_col_index = vcf_sample_ID___2___col_index[sampleID]
				line_to_print = line_to_print + "\t" + ln_elms[rel_col_index].strip()
			output_vcf_file.write(line_to_print + "\n")
		else:
			print "\nWARNING: duplicate variant ID found: " + variant_ID
input_vcf_file.close()
output_vcf_file.close()






####  Work on expr file
### build dict: expr_sample_ID___2___col_index
expr_sample_ID___2___col_index = {}
input_expr_file = open(input_expr_file_name, "r")
list_of_genes = list()
for line in input_expr_file:
	if "chr	start	end	gene	length	strand	" in line: ## ie, at header
		ln_elms = list()
		ln_elms = line.split("\t")
		gene = ln_elms[3]
		if gene not in list_of_genes:
			list_of_genes.append(gene)
			col_index = 6
			while col_index < len(ln_elms):
				sampleID = ln_elms[col_index].strip()
				expr_sample_ID___2___col_index[sampleID] = col_index
				col_index += 1
		else:
			print "\nWARNING: duplicate gene found: " + gene
input_expr_file.close()

### print output expr file
output_expr_file = open(output_expr_file_name, "w")

# (print header material)
header_line_to_print = "#chr	start	end	gene	length	strand"
for sampleID in list_of_INTERSECTED_sampleIDs_in_vcf_expr_cov_files:
	header_line_to_print = header_line_to_print + "\t" + sampleID 
output_expr_file.write(header_line_to_print + "\n")

# (go through each line and print out relevant columns)
input_expr_file = open(input_expr_file_name, "r")
for line in input_expr_file:
	if "#" not in line:  ## we've reached gene line
		ln_elms = list()
		ln_elms = line.split("\t")
		line_to_print = ln_elms[0] + "\t" + ln_elms[1] + "\t" + ln_elms[2] + "\t" + ln_elms[3] + "\t" + ln_elms[4] + "\t" + ln_elms[5]
		for sampleID in list_of_INTERSECTED_sampleIDs_in_vcf_expr_cov_files:
			rel_col_index = expr_sample_ID___2___col_index[sampleID]
			line_to_print = line_to_print + "\t" + ln_elms[rel_col_index].strip()
		output_expr_file.write(line_to_print + "\n")
input_expr_file.close()
output_expr_file.close()






####  Work on cov file
### build dict: cov_sample_ID___2___col_index
cov_sample_ID___2___col_index = {}
input_cov_file = open(input_cov_file_name, "r")
for line in input_cov_file:
	if "SampleID" in line: ## ie, at header
		ln_elms = list()
		ln_elms = line.split("\t")
		col_index = 1
		while col_index < len(ln_elms):
			sampleID = ln_elms[col_index].strip()
			cov_sample_ID___2___col_index[sampleID] = col_index
			col_index += 1
input_cov_file.close()

### print output cov file
output_cov_file = open(output_cov_file_name, "w")

# (print header material)
header_line_to_print = "SampleID"
for sampleID in list_of_INTERSECTED_sampleIDs_in_vcf_expr_cov_files:
	header_line_to_print = header_line_to_print + "\t" + sampleID 
output_cov_file.write(header_line_to_print + "\n")

# (go through each line and print out relevant columns)
input_cov_file = open(input_cov_file_name, "r")
for line in input_cov_file:
	if "SampleID" not in line:  ## we've reached a covariate line
		ln_elms = list()
		ln_elms = line.split("\t")
		line_to_print = ln_elms[0]
		for sampleID in list_of_INTERSECTED_sampleIDs_in_vcf_expr_cov_files:
			rel_col_index = cov_sample_ID___2___col_index[sampleID]
			line_to_print = line_to_print + "\t" + ln_elms[rel_col_index].strip()
		output_cov_file.write(line_to_print + "\n")
input_cov_file.close()
output_cov_file.close()





print "\n\n -- fin -- \n"

