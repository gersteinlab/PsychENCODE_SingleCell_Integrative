import numpy as np
import sys
import os
import random
from operator import itemgetter, attrgetter, methodcaller
# import matplotlib.pyplot as plt
# import os.path

'''
USAGE:
python build_expr_matrix_frm_xz_matrix_v3.py Astro list_of_intersected_genes.dat list_retained_indivs_for__Astro.txt /ysm-gpfs/pi/gerstein/dc547/brain/v4_scqtls/ps_bulk_umi_matrs/cpm_vanilla_normalized/ gene_info__ordered.dat

'''


cell_type = sys.argv[1]
print "\nCurrently processing cell type " + cell_type	+ " ... \n" 



##  Build list_of_genes
list_of_genes = list()
file_w_list_of_genes = open(sys.argv[2], "r")
ordered_list_of_sampleIDs = list()
for line in file_w_list_of_genes:
	gene = line.strip()
	list_of_genes.append(gene)
file_w_list_of_genes.close()
'''
for e in list_of_genes:
	print e
'''



def fix_cell_type_name(curr_cell_type_pre):
	revised_cell_type_name = ""
	for ch in curr_cell_type_pre:
		if ch == "/":
			revised_ch = "."
		elif ch == " ":
			revised_ch = "."
		else:
			revised_ch = ch
		revised_cell_type_name = revised_cell_type_name + revised_ch
	return revised_cell_type_name



##  Build ordered_list_of_sampleIDs
file_w_ordered_list_of_sampleIDs = open(sys.argv[3], "r")
ordered_list_of_sampleIDs = list()
for line in file_w_ordered_list_of_sampleIDs:
	sampleID = line.strip()
	ordered_list_of_sampleIDs.append(sampleID)
file_w_ordered_list_of_sampleIDs.close()
'''
for e in ordered_list_of_sampleIDs:
	print e
'''



'''
##  Build ordered_list_of_sampleIDs_containing_cell_type
ordered_list_of_sampleIDs_containing_cell_type = list()
dir_w_input_matrices = sys.argv[2]
list_of_input_matrix_file_names = os.listdir(dir_w_input_matrices)
for file_name in list_of_input_matrix_file_names:
	file_name_elms = list()
	file_name_elms = file_name.split("\t")
	sampleID = file_name_elms[0]
	if sampleID in ordered_list_of_sampleIDs:
		expr_matrix_file_name = dir_w_input_matrices + "/" + sampleID + ".bed"
		expr_matrix_file = open(expr_matrix_file_name, "r")
		for line in expr_matrix_file:
			ln_elms = list()
			ln_elms = line.split("\t")
			if ln_elms[0] == "featurekey":  ## at header line
				col_index = 1
				while col_index < len(ln_elms):
					curr_cell_type_pre = ln_elms[col_index].strip()
					curr_cell_type = fix_cell_type_name(curr_cell_type_pre)
					if curr_cell_type == cell_type:
						if sampleID not in ordered_list_of_sampleIDs_containing_cell_type:
							ordered_list_of_sampleIDs_containing_cell_type.append(sampleID)
					col_index += 1
		expr_matrix_file.close()
"""
for sampleID__gene__tag in sampleID__gene__tag____2____expr_val:
	print sampleID__gene__tag + "\t" + sampleID__gene__tag____2____expr_val[sampleID__gene__tag]
"""
'''



##  Build dict sampleID__gene__tag____2____expr_val
sampleID__gene__tag____2____expr_val = {}
dir_w_input_matrices = sys.argv[4]
list_of_input_matrix_file_names = os.listdir(dir_w_input_matrices)
list_of_genes_to_keep = list()
for file_name in list_of_input_matrix_file_names:
	sampleID = file_name[13:-30]
	if sampleID in ordered_list_of_sampleIDs:
		#print "Currently looking at sampleID " + sampleID
		expr_matrix_file_name = dir_w_input_matrices + "/" + file_name
		expr_matrix_file = open(expr_matrix_file_name, "r")
		for line in expr_matrix_file:
			ln_elms = list()
			ln_elms = line.split("\t")
			####print "reached here"
			if ln_elms[0] == "featurekey":  ## at header line
				col_index = 1
				#print "       Reached here 0"
				while col_index < len(ln_elms):
					curr_cell_type_pre = ln_elms[col_index].strip()
					curr_cell_type = fix_cell_type_name(curr_cell_type_pre)
					#print curr_cell_type_pre + "\t\t" + curr_cell_type
					if curr_cell_type == cell_type:
						#print "       Reached here 1 -- w/in sample " + sampleID + "   curr_cell_type = " + curr_cell_type
						col_index_associated_w_needed_cell_type = col_index
					col_index += 1
			else:
				gene = ln_elms[0]
				sampleID__gene__tag = sampleID + "__" + gene
				#print "\n\nsampleID__gene__tag  =  " + sampleID__gene__tag
				#print "col_index_associated_w_needed_cell_type = " + str(col_index_associated_w_needed_cell_type)
				#print "ln_elms  =  " + str(ln_elms)
				expr_val_of_gene_in_rel_cell_type = ln_elms[col_index_associated_w_needed_cell_type].strip()
				sampleID__gene__tag____2____expr_val[sampleID__gene__tag] = expr_val_of_gene_in_rel_cell_type
				list_of_genes_to_keep.append(gene)
		expr_matrix_file.close()
'''
for sampleID__gene__tag in sampleID__gene__tag____2____expr_val:
	print sampleID__gene__tag + "\t" + sampleID__gene__tag____2____expr_val[sampleID__gene__tag]
'''



##  Build dictionaries w/info needed for each gene
genes_info_file = open(sys.argv[5])
gene___2___chr = {}
gene___2___start = {}
gene___2___end = {}
gene___2___length = {}
gene___2___strand = {}
ordered_list_of_genes_in_genome = list()
for line in genes_info_file:
	if ("start" not in line)  and  ("end" not in line):  ## ie, if not in header line
		ln_elms = list()
		ln_elms = line.split("\t")
		chrom = ln_elms[0]
		start = ln_elms[1]
		end = ln_elms[2]
		gene = ln_elms[3]
		length = ln_elms[4]
		strand = ln_elms[5]
		gene___2___chr[gene] = chrom
		gene___2___start[gene] = start
		gene___2___end[gene] = end
		gene___2___length[gene] = length
		gene___2___strand[gene] = strand
		ordered_list_of_genes_in_genome.append(gene)
genes_info_file.close()



## Print output file
out_file_name = cell_type + "__v5.bed"
out_file = open(out_file_name, "w")

## First print header line (using the ordered_list_of_sampleIDs)
header_line = "#chr" + "\t" + "start" + "\t" + "end" + "\t" + "gene" + "\t" + "length" + "\t" + "strand"
for sampleID in ordered_list_of_sampleIDs:
	header_line = header_line + "\t" + sampleID
out_file.write(header_line + "\n")
header_line_elms = list()
header_line_elms = header_line.split("\t")
num_cols_in_header = len(header_line_elms)



## Now print the expr values for each gene
list_of_genes_already_incl_in_output = list()
for gene in ordered_list_of_genes_in_genome:
	if gene in list_of_genes_to_keep:
		chrom = gene___2___chr[gene]
		start = gene___2___start[gene]
		end = gene___2___end[gene]
		length = gene___2___length[gene]
		strand = gene___2___strand[gene].strip()
		gene_has_zero_expr_accr_all_samples = True
		first_six_elms_of_line_to_print = chrom + "\t" + start + "\t" + end + "\t" + gene + "\t" + length + "\t" + strand
		if first_six_elms_of_line_to_print not in list_of_genes_already_incl_in_output:
			list_of_genes_already_incl_in_output.append(first_six_elms_of_line_to_print)
			line_to_print = first_six_elms_of_line_to_print.strip()
			#print "       Reached here 2"
			for sampleID in ordered_list_of_sampleIDs:
				sampleID__gene__tag = sampleID + "__" + gene
				if sampleID__gene__tag in sampleID__gene__tag____2____expr_val:
					expr_val = sampleID__gene__tag____2____expr_val[sampleID__gene__tag]
					#print "       Reached here 3"
					if float(expr_val) > 0.0:
						gene_has_zero_expr_accr_all_samples = False
						#print "       Reached here 5"
					line_to_print = line_to_print + "\t" + expr_val.strip()
			line_to_print_elms = list()
			line_to_print_elms = line_to_print.split("\t")
			num_cols_in_curr_line = len(line_to_print_elms)
			if num_cols_in_curr_line == num_cols_in_header:  ## we want to exclude genes in which the number of fields in the expr matrix does not match the number of fields in the header of the file
				#print "       Reached here 7"
				if gene_has_zero_expr_accr_all_samples == False: ## we don't want to call QTLs on genes in which all expr vals = 0!
					#print "       Reached here 6"
					out_file.write(line_to_print + "\n")
			else:
				print "\n\n# cols w/in this row [= " + str(num_cols_in_curr_line) + "]  does not match # of expected # of cols (ie, # cols in header) [= " + str(num_cols_in_header) + "]"
				print line_to_print + "\n\n"
				sys.exit()
		else:
			print "Duplicate gene entry in input expr matrix:"
			print first_six_elms_of_line_to_print + "\n\n"
			sys.exit()
out_file.close()



print "\n ... Just finished w/cell type " + cell_type + "\n\n"
