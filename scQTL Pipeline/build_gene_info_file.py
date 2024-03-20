import numpy as np
import sys
import os
import random
from operator import itemgetter, attrgetter, methodcaller
# import matplotlib.pyplot as plt
# import os.path

'''
USAGE:
python build_gene_info_file.py list_of_genes_in__CMC_and_SZBD.dat /Users/dc/Desktop/rsch/pec2/qtl_tools/genes.gtf gene_info__CMC_Kellis.dat

'''


## Extact data from input file
input_gene_list_file = open(sys.argv[1], "r")
list_of_genes_in_input_file = list()
for line in input_gene_list_file:
	ln_elms = list()
	ln_elms = line.split("\t")
	gene = ln_elms[0].strip()
	if gene not in list_of_genes_in_input_file:
		list_of_genes_in_input_file.append(gene)
input_gene_list_file.close()


## Build dictionaries w/info needed for each gene
genes_info_file = open(sys.argv[2], "r")
gene___2___chr = {}
gene___2___start = {}
gene___2___end = {}
gene___2___length = {}
gene___2___strand = {}
for line in genes_info_file:
	if line[0] != "#": # ie: if not at header line
		ln_elms = list()
		ln_elms = line.split("\t")
		if ln_elms[2] == "gene":
			gene_name_info_field = ln_elms[8]
			gene_name_info_field_elms = list()
			gene_name_info_field_elms = gene_name_info_field.split(";")
			gene_name_component = gene_name_info_field_elms[3] ## ex: gene_name "LINC01409"
			gene_name_component_elms = list()
			gene_name_component_elms = gene_name_component.split()
			if gene_name_component_elms[0] == "gene_name":
				gene_name_of_curr_row_pre = gene_name_component_elms[1]
				gene_name_of_curr_row = gene_name_of_curr_row_pre[1:-1]
			else:
				print "\n\nERROR: unrecognized gene_name_component sub-field: " + gene_name_component + "\n\n"
				sys.exit("\n\nERROR: unrecognized gene_name_component sub-field: " + gene_name_component + "\n\n")

			if gene_name_of_curr_row in list_of_genes_in_input_file:
				chrom = ln_elms[0]
				strand = ln_elms[6]
				if strand == "+":
					start_TSS = int(ln_elms[3])
					end_TSS = start_TSS + 1
					end = int(ln_elms[4])
				elif strand == "-":
					start_TSS = int(ln_elms[4])
					end_TSS = start_TSS + 1
					end = int(ln_elms[3])
				length = abs(end - start_TSS)
				gene___2___chr[gene_name_of_curr_row] = chrom
				gene___2___start[gene_name_of_curr_row] = start_TSS
				gene___2___end[gene_name_of_curr_row] = end_TSS
				gene___2___length[gene_name_of_curr_row] = length
				gene___2___strand[gene_name_of_curr_row] = strand
genes_info_file.close()


## Build the header line and print it to the output file
output_file = open(sys.argv[3], "w")
header_line_to_print = "#chr	start	end	gene	length	strand"
output_file.write(header_line_to_print + "\n")


## Finally -- print the gene info tooutput file
for gene in list_of_genes_in_input_file:
	if (gene in gene___2___chr)  and  (gene in gene___2___start)  and  (gene in gene___2___end)  and  (gene in gene___2___length)  and  (gene in gene___2___strand):
		line_to_print = gene___2___chr[gene] + "\t" + str(gene___2___start[gene]) + "\t" + str(gene___2___end[gene]) + "\t" + gene + "\t" + str(gene___2___length[gene]) + "\t" + gene___2___strand[gene]
		output_file.write(line_to_print + "\n")
output_file.close()


print "\n\n -- fin -- \n"

