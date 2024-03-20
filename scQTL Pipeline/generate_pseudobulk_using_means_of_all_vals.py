import numpy as np
import sys
import os
import random
from operator import itemgetter, attrgetter, methodcaller
# import matplotlib.pyplot as plt
# import os.path

'''
USAGE:
python generate_pseudobulk_using_means_of_all_vals.py list_of_21_cell_types.dat ./input_umi_dir/ ./ps_bulk_umi_matrs/ 300 50 PEC2_all_cohorts_sample_mapping_adults__wo_duplicates_filt.txt

'''


tot_nuclei_thresh = int(sys.argv[4])
num_nuclei_thresh = int(sys.argv[5])


input_file_w_cell_types = open(sys.argv[1], "r")
list_of_cell_types = list()
for line in input_file_w_cell_types:
	cell_type = line.strip()
	list_of_cell_types.append(cell_type)
input_file_w_cell_types.close()


list_of_expr_samples_to_process = list()
PEC2_all_cohorts_sample_mapping__file = open(sys.argv[6], "r")
for line in PEC2_all_cohorts_sample_mapping__file:
	if "Cohort	scRNA_ID	GenotypeID	Covariate_ID" not in line:  ## ie -- if not at header line
		ln_elms = list()
		ln_elms = line.split("\t")
		expr_sampleID = ln_elms[1]
		if expr_sampleID not in list_of_expr_samples_to_process:
			list_of_expr_samples_to_process.append(expr_sampleID)
PEC2_all_cohorts_sample_mapping__file.close()


dir_w_input_UMI_files = sys.argv[2]
dir_w_ps_bulk_UMI_files = sys.argv[3]
for file_name in os.listdir(dir_w_input_UMI_files):
	current_expr_sampleID = file_name[0:-21]
	if current_expr_sampleID in list_of_expr_samples_to_process:
		full_file_name = dir_w_input_UMI_files + "/" + file_name

		'''
		cell_type___2___list_of_retained_indivs = {}
		for cell_type in list_of_cell_types:
			list_of_retained_indivs = list()
			file_name_w_retained_indivs = "list_retained_indivs_for__" + cell_type + ".txt"
			file_w_retained_indivs = open(file_name_w_retained_indivs, "r")
			for line in file_w_retained_indivs:
				sampleID = line.strip()
				list_of_retained_indivs.append(sampleID)
			file_w_retained_indivs.close()
			cell_type___2___list_of_retained_indivs[cell_type] = list_of_retained_indivs
		#for cell_type in cell_type___2___list_of_retained_indivs:
		#	print cell_type + "\n" + str(cell_type___2___list_of_retained_indivs[cell_type])
		#	print "\n\n\n"
		'''

		##  Build dictionary (cell_type___2___list_of_col_indeces) to keep track of which rows correspond to which cell types
		cell_type___2___list_of_col_indeces = {}
		list_of_cell_types_in_indiv = list()
		for cell_type in list_of_cell_types:
			cell_type___2___list_of_col_indeces[cell_type] = list()
		input_file = open(full_file_name, "r")
		for line in input_file:
			ln_elms = list()
			ln_elms = line.split("\t")
			tot_num_nuclei = len(ln_elms) - 1
			if tot_num_nuclei >= tot_nuclei_thresh:
				if ln_elms[0] == "featurekey":   ## we're at the header line
					curr_col_index = 1
					while curr_col_index < len(ln_elms):
						nuc_name_pre = ln_elms[curr_col_index]
						nuc_name_elms = list()
						nuc_name_pre_elms = list()
						if ("Sst" in nuc_name_pre):
							cell_type = "Sst__Sst.Chodl"
						elif ("Endo" in nuc_name_pre)  or  ("VLMC" in nuc_name_pre):
							cell_type = "Endo__VLMC"
						elif ("Chandelier" in nuc_name_pre)  or  ("Pvalb" in nuc_name_pre):
							cell_type = "Chandelier__Pvalb"
						elif "L6 IT Car3" in nuc_name_pre:
							cell_type = "L6 IT Car3"
						elif "L6 IT" in nuc_name_pre:
							cell_type = "L6 IT"
						elif "Lamp5 Lhx6" in nuc_name_pre:
							cell_type = "Lamp5 Lhx6"
						elif "Lamp5" in nuc_name_pre:
							cell_type = "Lamp5"
						elif "Astro" in nuc_name_pre:
							cell_type = "Astro"
						elif "L2/3 IT" in nuc_name_pre:
							cell_type = "L2/3 IT"
						elif "L4 IT" in nuc_name_pre:
							cell_type = "L4 IT"
						elif "L5/6 NP" in nuc_name_pre:
							cell_type = "L5/6 NP"
						elif "L5 ET" in nuc_name_pre:
							cell_type = "L5 ET"
						elif "L5 IT" in nuc_name_pre:
							cell_type = "L5 IT"
						elif "L6 CT" in nuc_name_pre:
							cell_type = "L6 CT"
						elif "Micro/PVM" in nuc_name_pre:
							cell_type = "Micro/PVM"
						elif "Micro" in nuc_name_pre:
							cell_type = "Micro/PVM"
						elif "OPC" in nuc_name_pre:
							cell_type = "OPC"
						elif "Oligo" in nuc_name_pre:
							cell_type = "Oligo"
						elif "Sncg" in nuc_name_pre:
							cell_type = "Sncg"
						elif "Vip" in nuc_name_pre:
							cell_type = "Vip"
						elif "L6b" in nuc_name_pre:
							cell_type = "L6b"
						elif "Pax6" in nuc_name_pre:
							cell_type = "Pax6"
						elif "Immune" in nuc_name_pre:
							cell_type = "Immune"
						elif "PC" in nuc_name_pre:   #########################
							cell_type = "PC"   #########################
						elif "SMC" in nuc_name_pre:
							cell_type = "SMC"
						elif "RB" in nuc_name_pre:
							cell_type = "RB"
						else:
							print "\nUNRECOGNIZED CELL TYPE: " + nuc_name_pre + "\n"
							#sys.exit()
						#print "curr cell_type: " + cell_type + "\n"
						if cell_type not in list_of_cell_types_in_indiv:
							list_of_cell_types_in_indiv.append(cell_type)
						list_of_col_indeces = list()
						list_of_col_indeces = cell_type___2___list_of_col_indeces[cell_type]
						list_of_col_indeces.append(curr_col_index)   ####################################
						cell_type___2___list_of_col_indeces[cell_type] = list_of_col_indeces
						curr_col_index += 1
		input_file.close()
		#for cell_type in list_of_cell_types:
		#	print cell_type
		#	print str(cell_type___2___list_of_col_indeces[cell_type])
		#	print "\n\n\n\n"

		if tot_num_nuclei >= tot_nuclei_thresh:
			##  Within each cell type -- determine the mean expression for each gene
			input_file = open(full_file_name, "r")
			out_file_name = dir_w_ps_bulk_UMI_files + "/pseudo_bulk__" + file_name
			#out_file_name = dir_w_ps_bulk_UMI_files + "/" + file_name
			out_file = open(out_file_name, "w")
			for line in input_file:
				ln_elms = list()
				ln_elms = line.split("\t")
				line_to_print = ln_elms[0]
				if ln_elms[0] != "featurekey":   ## we're at a gene line (ie, not at the header line)
					gene = ln_elms[0]
					#print "curr at gene " + gene
					for cell_type in list_of_cell_types_in_indiv:
						summed_UMI_counts_for_curr_gene = 0.0
						num_nuclie_for_cell_type = len(cell_type___2___list_of_col_indeces[cell_type])
						if num_nuclie_for_cell_type >= num_nuclei_thresh:
							for col_index in cell_type___2___list_of_col_indeces[cell_type]:
								summed_UMI_counts_for_curr_gene += float(ln_elms[col_index])
							mean_exp_for_gene_in_cell_type = float(summed_UMI_counts_for_curr_gene) / float(num_nuclie_for_cell_type)
							#print "\t\t" + cell_type + "\t" + str(summed_UMI_counts_for_curr_gene) + "\t" + str(num_nuclie_for_cell_type) + "\t" + str(mean_exp_for_gene_in_cell_type)
							line_to_print = line_to_print + "\t" + str(mean_exp_for_gene_in_cell_type)
						#else:
						#	print "Within " + full_file_name + ":  Insufficient number of nuclei found for cell type " + cell_type + "\n"
				else:  ## header line
					for cell_type in list_of_cell_types_in_indiv:
						num_nuclie_for_cell_type = len(cell_type___2___list_of_col_indeces[cell_type])
						if num_nuclie_for_cell_type >= num_nuclei_thresh:
							line_to_print = line_to_print + "\t" + cell_type
				out_file.write(line_to_print + "\n")
			input_file.close()
			out_file.close()


print "\n\n -- fin -- \n"

