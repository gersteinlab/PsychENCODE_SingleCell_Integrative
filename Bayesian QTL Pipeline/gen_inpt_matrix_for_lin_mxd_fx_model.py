import numpy as np
import sys
import os
import random
from operator import itemgetter, attrgetter, methodcaller
# import matplotlib.pyplot as plt
# import os.path

'''
USAGE:
python gen_inpt_matrix_for_lin_mxd_fx_model.py sig_eGene_SNP_pairs.dat list_of_all_cell_types.dat /Users/dc/Desktop/rsch/pec2/bayesian_lin_mixed_FX_models/sample_input_expr_dir/ /Users/dc/Desktop/rsch/pec2/bayesian_lin_mixed_FX_models/sample_input_vcf_dir/ /Users/dc/Desktop/rsch/pec2/bayesian_lin_mixed_FX_models/sample_input_cov_dir/

'''

input_file_name_w_eGene_eSNP_pairs = sys.argv[1]
file_name_w_list_of_cellTypes = sys.argv[2]
dir_w_expr_matrices = sys.argv[3]
dir_w_vcf_matrices = sys.argv[4]
dir_w_cov_matrices = sys.argv[5]

input_file_w_eGene_eSNP_pairs = open(input_file_name_w_eGene_eSNP_pairs, "r")
for inpt_line in input_file_w_eGene_eSNP_pairs:
	inpt_line_elms = list()
	inpt_line_elms = inpt_line.split()
	eGene = inpt_line_elms[0]
	eSNP = inpt_line_elms[1].strip()


	out_file_name = eGene + "__" + eSNP + "__" + "input_for_mxd_fx_model.dat"
	out_file = open(out_file_name, "w")


	list_of_cellTypes = list()
	file_w_list_of_cellTypes = open(file_name_w_list_of_cellTypes, "r")
	for line in file_w_list_of_cellTypes:
		cellType = line.strip()
		list_of_cellTypes.append(cellType)
	file_w_list_of_cellTypes.close()


	cellType__sampleID_tag____2____expr = {}
	cellType___2___list_of_all_sampleIDs = {}
	list_of_cellTypes_w_this_gene = list()
	for cellType in list_of_cellTypes:
		#print "curr cellType  =  " + cellType
		list_of_all_sampleIDs_for_curr_cellType = list()
		input_exp_file_name = dir_w_expr_matrices + "/" + cellType + ".expr.bed"  ######### will need to add ".gz" at end
		expr_file = open(input_exp_file_name, "r")
		col_index___2___sampleID = {}
		col_index___2___sampleID.clear()
		for line in expr_file:
			#print "curr cellType  =  " + cellType
			ln_elms = list()
			ln_elms = line.split()
			if ln_elms[0] == "#chr":  ## at header line
				col_index = 6
				while col_index < len(ln_elms):
					sampleID = ln_elms[col_index].strip()
					col_index___2___sampleID[col_index] = sampleID
					list_of_all_sampleIDs_for_curr_cellType.append(sampleID)
					col_index += 1
			else:
				curr_gene = ln_elms[3]
				if curr_gene == eGene:
					list_of_cellTypes_w_this_gene.append(cellType)
					col_index = 6
					while col_index < len(ln_elms):
						expr_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr[cellType__sampleID_tag] = float(expr_val)
						col_index += 1
		cellType___2___list_of_all_sampleIDs[cellType] = list_of_all_sampleIDs_for_curr_cellType
		expr_file.close()
	'''
	for elm in cellType__sampleID_tag____2____expr:
		print elm + "\t" + str(cellType__sampleID_tag____2____expr[elm])
	'''


	cellType__sampleID_tag____2____genotypeDosage = {}
	for cellType in list_of_cellTypes:
		input_vcf_file_name = dir_w_vcf_matrices + cellType + ".vcf"  ######### will need to add ".gz" at end
		vcf_file = open(input_vcf_file_name, "r")
		col_index___2___sampleID = {}
		col_index___2___sampleID.clear()
		for line in vcf_file:
			if "##fileformat=VCFv4.3" not in line:
				ln_elms = list()
				ln_elms = line.split()
				if ln_elms[0] == "#CHROM":  ## at header line
					col_index = 9
					while col_index < len(ln_elms):
						sampleID = ln_elms[col_index].strip()
						col_index___2___sampleID[col_index] = sampleID
						col_index += 1
				else:
					curr_SNP_ID = ln_elms[2]
					if curr_SNP_ID == eSNP:
						col_index = 9
						while col_index < len(ln_elms):
							genotypeDosage_pre = ln_elms[col_index].strip()
							if genotypeDosage_pre == "0/0":
								genotypeDosage = 0.0
							elif (genotypeDosage_pre == "1/0")  or  (genotypeDosage_pre == "0/1"):
								genotypeDosage = 1.0
							elif (genotypeDosage_pre == "1/1"):
								genotypeDosage = 2.0
							else:
								print "\n\nERROR: unrecognized genotypeDosage_pre in vcf file : " + genotypeDosage_pre
								sys.exit()
							sampleID = col_index___2___sampleID[col_index]
							cellType__sampleID_tag = cellType + "__" + sampleID
							cellType__sampleID_tag____2____genotypeDosage[cellType__sampleID_tag] = genotypeDosage
							col_index += 1
		vcf_file.close()


	cellType__sampleID_tag____2____gender_is_female = {}
	cellType__sampleID_tag____2____ageDeath = {}
	cellType__sampleID_tag____2____diagnosis_is_ASD = {}
	cellType__sampleID_tag____2____diagnosis_is_Alzheimers_dementia = {}
	cellType__sampleID_tag____2____diagnosis_is_Bipolar_Disorder = {}
	cellType__sampleID_tag____2____diagnosis_is_Schizophrenia = {}
	cellType__sampleID_tag____2____diagnosis_is_cognitive_impairment = {}
	cellType__sampleID_tag____2____diagnosis_is_control = {}
	cellType__sampleID_tag____2____diagnosis_MDD = {}
	cellType__sampleID_tag____2____diagnosis_PTSD = {}
	cellType__sampleID_tag____2____diagnosis_Williams = {}
	cellType__sampleID_tag____2____cohort_is_CMC = {}
	cellType__sampleID_tag____2____cohort_is_DevBrain = {}
	cellType__sampleID_tag____2____cohort_is_IsoHuB = {}
	cellType__sampleID_tag____2____cohort_is_Kreigstein = {}
	cellType__sampleID_tag____2____cohort_is_ROSMAP = {}
	cellType__sampleID_tag____2____cohort_is_SZBD = {}
	cellType__sampleID_tag____2____cohort_is_UCLA_ASD = {}
	cellType__sampleID_tag____2____cohort_is_Urban = {}
	cellType__sampleID_tag____2____cohort_Girgenti = {}
	cellType__sampleID_tag____2____cohort_LIBD = {}
	cellType__sampleID_tag____2____cohort_Ma_Sestan = {}
	cellType__sampleID_tag____2____cohort_PTSDBrainomics = {}
	cellType__sampleID_tag____2____gt_PC_1 = {}
	cellType__sampleID_tag____2____gt_PC_2 = {}
	cellType__sampleID_tag____2____gt_PC_3 = {}
	cellType__sampleID_tag____2____gt_PC_4 = {}
	cellType__sampleID_tag____2____gt_PC_5 = {}
	cellType__sampleID_tag____2____expr_PC1 = {}
	cellType__sampleID_tag____2____expr_PC2 = {}
	cellType__sampleID_tag____2____expr_PC3 = {}
	cellType__sampleID_tag____2____expr_PC4 = {}
	cellType__sampleID_tag____2____expr_PC5 = {}
	cellType__sampleID_tag____2____expr_PC6 = {}
	cellType__sampleID_tag____2____expr_PC7 = {}
	cellType__sampleID_tag____2____expr_PC8 = {}
	cellType__sampleID_tag____2____expr_PC9 = {}
	cellType__sampleID_tag____2____expr_PC10 = {}
	cellType__sampleID_tag____2____expr_PC11 = {}
	cellType__sampleID_tag____2____expr_PC12 = {}
	cellType__sampleID_tag____2____expr_PC13 = {}
	cellType__sampleID_tag____2____expr_PC14 = {}
	cellType__sampleID_tag____2____expr_PC15 = {}
	cellType__sampleID_tag____2____expr_PC16 = {}
	cellType__sampleID_tag____2____expr_PC17 = {}
	cellType__sampleID_tag____2____expr_PC18 = {}
	cellType__sampleID_tag____2____expr_PC19 = {}
	cellType__sampleID_tag____2____expr_PC20 = {}
	cellType__sampleID_tag____2____expr_PC21 = {}
	cellType__sampleID_tag____2____expr_PC22 = {}
	cellType__sampleID_tag____2____expr_PC23 = {}
	cellType__sampleID_tag____2____expr_PC24 = {}
	cellType__sampleID_tag____2____expr_PC25 = {}
	cellType__sampleID_tag____2____expr_PC26 = {}
	cellType__sampleID_tag____2____expr_PC27 = {}
	cellType__sampleID_tag____2____expr_PC28 = {}
	cellType__sampleID_tag____2____expr_PC29 = {}
	cellType__sampleID_tag____2____expr_PC30 = {}
	cellType__sampleID_tag____2____expr_PC31 = {}
	cellType__sampleID_tag____2____expr_PC32 = {}
	cellType__sampleID_tag____2____expr_PC33 = {}
	cellType__sampleID_tag____2____expr_PC34 = {}
	cellType__sampleID_tag____2____expr_PC35 = {}
	cellType__sampleID_tag____2____expr_PC36 = {}
	cellType__sampleID_tag____2____expr_PC37 = {}
	cellType__sampleID_tag____2____expr_PC38 = {}
	cellType__sampleID_tag____2____expr_PC39 = {}
	cellType__sampleID_tag____2____expr_PC40 = {}
	cellType__sampleID_tag____2____expr_PC41 = {}
	cellType__sampleID_tag____2____expr_PC42 = {}
	cellType__sampleID_tag____2____expr_PC43 = {}
	cellType__sampleID_tag____2____expr_PC44 = {}
	cellType__sampleID_tag____2____expr_PC45 = {}
	cellType__sampleID_tag____2____expr_PC46 = {}
	cellType__sampleID_tag____2____expr_PC47 = {}
	cellType__sampleID_tag____2____expr_PC48 = {}
	cellType__sampleID_tag____2____expr_PC49 = {}
	cellType__sampleID_tag____2____expr_PC50 = {}
	cellType__sampleID_tag____2____expr_PC51 = {}
	cellType__sampleID_tag____2____expr_PC52 = {}
	cellType__sampleID_tag____2____expr_PC53 = {}
	cellType__sampleID_tag____2____expr_PC54 = {}
	cellType__sampleID_tag____2____expr_PC55 = {}
	cellType__sampleID_tag____2____expr_PC56 = {}
	cellType__sampleID_tag____2____expr_PC57 = {}
	cellType__sampleID_tag____2____expr_PC58 = {}
	cellType__sampleID_tag____2____expr_PC59 = {}
	cellType__sampleID_tag____2____expr_PC60 = {}
	cellType__sampleID_tag____2____expr_PC61 = {}
	cellType__sampleID_tag____2____expr_PC62 = {}
	cellType__sampleID_tag____2____expr_PC63 = {}
	cellType__sampleID_tag____2____expr_PC64 = {}
	cellType__sampleID_tag____2____expr_PC65 = {}
	cellType__sampleID_tag____2____expr_PC66 = {}
	cellType__sampleID_tag____2____expr_PC67 = {}
	cellType__sampleID_tag____2____expr_PC68 = {}
	cellType__sampleID_tag____2____expr_PC69 = {}
	cellType__sampleID_tag____2____expr_PC70 = {}
	cellType__sampleID_tag____2____expr_PC71 = {}
	cellType__sampleID_tag____2____expr_PC72 = {}
	cellType__sampleID_tag____2____expr_PC73 = {}
	cellType__sampleID_tag____2____expr_PC74 = {}
	cellType__sampleID_tag____2____expr_PC75 = {}
	cellType__sampleID_tag____2____expr_PC76 = {}
	cellType__sampleID_tag____2____expr_PC77 = {}
	cellType__sampleID_tag____2____expr_PC78 = {}
	cellType__sampleID_tag____2____expr_PC79 = {}
	cellType__sampleID_tag____2____expr_PC80 = {}
	cellType__sampleID_tag____2____expr_PC81 = {}
	cellType__sampleID_tag____2____expr_PC82 = {}
	cellType__sampleID_tag____2____expr_PC83 = {}
	cellType__sampleID_tag____2____expr_PC84 = {}
	cellType__sampleID_tag____2____expr_PC85 = {}
	cellType__sampleID_tag____2____expr_PC86 = {}
	cellType__sampleID_tag____2____expr_PC87 = {}
	cellType__sampleID_tag____2____expr_PC88 = {}
	cellType__sampleID_tag____2____expr_PC89 = {}
	cellType__sampleID_tag____2____expr_PC90 = {}
	cellType__sampleID_tag____2____expr_PC91 = {}
	cellType__sampleID_tag____2____expr_PC92 = {}
	cellType__sampleID_tag____2____expr_PC93 = {}
	cellType__sampleID_tag____2____expr_PC94 = {}
	cellType__sampleID_tag____2____expr_PC95 = {}
	cellType__sampleID_tag____2____expr_PC96 = {}
	cellType__sampleID_tag____2____expr_PC97 = {}
	cellType__sampleID_tag____2____expr_PC98 = {}
	cellType__sampleID_tag____2____expr_PC99 = {}
	cellType__sampleID_tag____2____expr_PC100 = {}
	for cellType in list_of_cellTypes:
		input_cov_file_name = dir_w_cov_matrices + cellType + ".cov.bed"  ######### will need to add ".gz" at end
		cov_file = open(input_cov_file_name, "r")
		col_index___2___sampleID = {}
		col_index___2___sampleID.clear()
		list_of_covariates = list()
		for line in cov_file:
			ln_elms = list()
			ln_elms = line.split()
			if ln_elms[0] == "SampleID":  ## at header line
				col_index = 1
				while col_index < len(ln_elms):
					sampleID = ln_elms[col_index].strip()
					col_index___2___sampleID[col_index] = sampleID
					col_index += 1
			else:
				covariate = ln_elms[0]
				list_of_covariates.append(covariate)
				if covariate == "gender_is_female":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____gender_is_female[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "ageDeath":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____ageDeath[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "diagnosis_is_ASD":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____diagnosis_is_ASD[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "diagnosis_is_Alzheimers_dementia":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____diagnosis_is_Alzheimers_dementia[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "diagnosis_is_Bipolar_Disorder":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____diagnosis_is_Bipolar_Disorder[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "diagnosis_is_Schizophrenia":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____diagnosis_is_Schizophrenia[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "diagnosis_is_cognitive_impairment":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____diagnosis_is_cognitive_impairment[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "diagnosis_MDD":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____diagnosis_MDD[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "diagnosis_PTSD":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____diagnosis_PTSD[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "diagnosis_Williams":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____diagnosis_Williams[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "diagnosis_is_control":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____diagnosis_is_control[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "cohort_is_CMC":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____cohort_is_CMC[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "cohort_is_DevBrain":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____cohort_is_DevBrain[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "cohort_is_IsoHuB":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____cohort_is_IsoHuB[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "cohort_is_Kreigstein":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____cohort_is_Kreigstein[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "cohort_is_ROSMAP":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____cohort_is_ROSMAP[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "cohort_is_SZBD":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____cohort_is_SZBD[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "cohort_is_UCLA-ASD":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____cohort_is_UCLA_ASD[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "cohort_is_Urban":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____cohort_is_Urban[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "cohort_Girgenti":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____cohort_Girgenti[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "cohort_LIBD":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____cohort_LIBD[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "cohort_Ma_Sestan":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____cohort_Ma_Sestan[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "cohort_PTSDBrainomics":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____cohort_PTSDBrainomics[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "gt_PC_1":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____gt_PC_1[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "gt_PC_2":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____gt_PC_2[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "gt_PC_3":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____gt_PC_3[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "gt_PC_4":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____gt_PC_4[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "gt_PC_5":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____gt_PC_5[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC1":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC1[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC2":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC2[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC3":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC3[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC4":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC4[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC5":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC5[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC6":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC6[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC7":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC7[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC8":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC8[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC9":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC9[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC10":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC10[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC11":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC11[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC12":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC12[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC13":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC13[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC14":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC14[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC15":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC15[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC16":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC16[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC17":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC17[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC18":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC18[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC19":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC19[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC20":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC20[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC21":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC21[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC22":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC22[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC23":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC23[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC24":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC24[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC25":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC25[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC26":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC26[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC27":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC27[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC28":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC28[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC29":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC29[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC30":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC30[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC31":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC31[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC32":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC32[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC33":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC33[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC34":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC34[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC35":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC35[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC36":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC36[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC37":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC37[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC38":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC38[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC39":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC39[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC40":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC40[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC41":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC41[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC42":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC42[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC43":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC43[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC44":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC44[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC45":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC45[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC46":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC46[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC47":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC47[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC48":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC48[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC49":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC49[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC50":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC50[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC51":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC51[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC52":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC52[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC53":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC53[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC54":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC54[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC55":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC55[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC56":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC56[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC57":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC57[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC58":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC58[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC59":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC59[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC60":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC60[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC61":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC61[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC62":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC62[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC63":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC63[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC64":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC64[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC65":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC65[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC66":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC66[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC67":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC67[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC68":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC68[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC69":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC69[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC70":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC70[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC71":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC71[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC72":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC72[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC73":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC73[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC74":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC74[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC75":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC75[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC76":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC76[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC77":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC77[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC78":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC78[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC79":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC79[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC80":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC80[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC81":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC81[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC82":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC82[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC83":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC83[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC84":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC84[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC85":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC85[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC86":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC86[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC87":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC87[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC88":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC88[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC89":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC89[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC90":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC90[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC91":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC91[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC92":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC92[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC93":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC93[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC94":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC94[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC95":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC95[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC96":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC96[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC97":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC97[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC98":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC98[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC99":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC99[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
				elif covariate == "expr_PC100":
					col_index = 1
					while col_index < len(ln_elms):
						cov_val = ln_elms[col_index].strip()
						sampleID = col_index___2___sampleID[col_index]
						cellType__sampleID_tag = cellType + "__" + sampleID
						cellType__sampleID_tag____2____expr_PC100[cellType__sampleID_tag] = float(cov_val)
						col_index += 1
		cov_file.close()


	header_line = "cellType" + "\t" + "genotypeDosage" + "\t" + "cov__gender_is_female" + "\t" + "cov__ageDeath" + "\t" + "cov__diagnosis_is_ASD" + "\t" + "cov__diagnosis_is_Alzheimers_dementia" + "\t" + "cov__diagnosis_is_Bipolar_Disorder" + "\t" + "cov__diagnosis_is_Schizophrenia" + "\t" + "cov__diagnosis_is_cognitive_impairment" + "\t" + "cov_diagnosis_MDD" + "\t" + "cov_diagnosis_PTSDr" + "\t" + "cov_diagnosis_Williams" + "\t" + "diagnosis_is_control" + "\t" + "cov__cohort_is_CMC" + "\t" + "cov__cohort_is_DevBrain" + "\t" + "cov__cohort_is_IsoHuB" + "\t" + "cov__cohort_is_Kreigstein" + "\t" + "cov__cohort_is_ROSMAP" + "\t" + "cov__cohort_is_SZBD" + "\t" + "cov__cohort_is_UCLA_ASD" + "\t" + "cov__cohort_is_Urban" + "\t" + "cov__cohort_Girgenti" + "\t" + "cov__cohort_LIBD" + "\t" + "cov__cohort_Ma_Sestan" + "\t" + "cov__cohort_PTSDBrainomics" + "\t" + "cov__gt_PC_1" + "\t" + "cov__gt_PC_2" + "\t" + "cov__gt_PC_3" + "\t" + "cov__gt_PC_4" + "\t" + "cov__gt_PC_5" + "\t" + "cov__expr_PC1" + "\t" + "cov__expr_PC2" + "\t" + "cov__expr_PC3" + "\t" + "cov__expr_PC4" + "\t" + "cov__expr_PC5" + "\t" + "cov__expr_PC6" + "\t" + "cov__expr_PC7" + "\t" + "cov__expr_PC8" + "\t" + "cov__expr_PC9" + "\t" + "cov__expr_PC10" + "\t" + "cov__expr_PC11" + "\t" + "cov__expr_PC12" + "\t" + "cov__expr_PC13" + "\t" + "cov__expr_PC14" + "\t" + "cov__expr_PC15" + "\t" + "cov__expr_PC16" + "\t" + "cov__expr_PC17" + "\t" + "cov__expr_PC18" + "\t" + "cov__expr_PC19" + "\t" + "cov__expr_PC20" + "\t" + "cov__expr_PC21" + "\t" + "cov__expr_PC22" + "\t" + "cov__expr_PC23" + "\t" + "cov__expr_PC24" + "\t" + "cov__expr_PC25" + "\t" + "cov__expr_PC26" + "\t" + "cov__expr_PC27" + "\t" + "cov__expr_PC28" + "\t" + "cov__expr_PC29" + "\t" + "cov__expr_PC30" + "\t" + "cov__expr_PC31" + "\t" + "cov__expr_PC32" + "\t" + "cov__expr_PC33" + "\t" + "cov__expr_PC34" + "\t" + "cov__expr_PC35" + "\t" + "cov__expr_PC36" + "\t" + "cov__expr_PC37" + "\t" + "cov__expr_PC38" + "\t" + "cov__expr_PC39" + "\t" + "cov__expr_PC40" + "\t" + "cov__expr_PC41" + "\t" + "cov__expr_PC42" + "\t" + "cov__expr_PC43" + "\t" + "cov__expr_PC44" + "\t" + "cov__expr_PC45" + "\t" + "cov__expr_PC46" + "\t" + "cov__expr_PC47" + "\t" + "cov__expr_PC48" + "\t" + "cov__expr_PC49" + "\t" + "cov__expr_PC50" + "\t" + "cov__expr_PC51" + "\t" + "cov__expr_PC52" + "\t" + "cov__expr_PC53" + "\t" + "cov__expr_PC54" + "\t" + "cov__expr_PC55" + "\t" + "cov__expr_PC56" + "\t" + "cov__expr_PC57" + "\t" + "cov__expr_PC58" + "\t" + "cov__expr_PC59" + "\t" + "cov__expr_PC60" + "\t" + "cov__expr_PC61" + "\t" + "cov__expr_PC62" + "\t" + "cov__expr_PC63" + "\t" + "cov__expr_PC64" + "\t" + "cov__expr_PC65" + "\t" + "cov__expr_PC66" + "\t" + "cov__expr_PC67" + "\t" + "cov__expr_PC68" + "\t" + "cov__expr_PC69" + "\t" + "cov__expr_PC70" + "\t" + "cov__expr_PC71" + "\t" + "cov__expr_PC72" + "\t" + "cov__expr_PC73" + "\t" + "cov__expr_PC74" + "\t" + "cov__expr_PC75" + "\t" + "cov__expr_PC76" + "\t" + "cov__expr_PC77" + "\t" + "cov__expr_PC78" + "\t" + "cov__expr_PC79" + "\t" + "cov__expr_PC80" + "\t" + "cov__expr_PC81" + "\t" + "cov__expr_PC82" + "\t" + "cov__expr_PC83" + "\t" + "cov__expr_PC84" + "\t" + "cov__expr_PC85" + "\t" + "cov__expr_PC86" + "\t" + "cov__expr_PC87" + "\t" + "cov__expr_PC88" + "\t" + "cov__expr_PC89" + "\t" + "cov__expr_PC90" + "\t" + "cov__expr_PC91" + "\t" + "cov__expr_PC92" + "\t" + "cov__expr_PC93" + "\t" + "cov__expr_PC94" + "\t" + "cov__expr_PC95" + "\t" + "cov__expr_PC96" + "\t" + "cov__expr_PC97" + "\t" + "cov__expr_PC98" + "\t" + "cov__expr_PC99" + "\t" + "cov__expr_PC100" + "\t" + "expr"
	out_file.write(header_line + "\n")

	'''
	for covariate in list_of_covariates:
		if (covariate != "diagnosis_is_control"):
			header_line = header_line + "\t" + covariate
	header_line = header_line + "\t" + "expr_val"
	'''
	#print "header_line  = \n" + header_line + "\n"


	for cellType in list_of_cellTypes:
		if cellType in list_of_cellTypes_w_this_gene:
			list_of_all_sampleIDs_for_curr_cellType = cellType___2___list_of_all_sampleIDs[cellType]
			if len(list_of_all_sampleIDs_for_curr_cellType) > 0:
				for sampleID in list_of_all_sampleIDs_for_curr_cellType:
					cellType__sampleID_tag = cellType + "__" + sampleID
					if (cellType__sampleID_tag in cellType__sampleID_tag____2____genotypeDosage):
						if (cellType__sampleID_tag in cellType__sampleID_tag____2____gender_is_female):
							if (cellType__sampleID_tag in cellType__sampleID_tag____2____expr):
								genotypeDosage = cellType__sampleID_tag____2____genotypeDosage[cellType__sampleID_tag]
								cov__gender_is_female = cellType__sampleID_tag____2____gender_is_female[cellType__sampleID_tag]
								cov__ageDeath = cellType__sampleID_tag____2____ageDeath[cellType__sampleID_tag]
								cov__diagnosis_is_ASD = cellType__sampleID_tag____2____diagnosis_is_ASD[cellType__sampleID_tag]
								cov__diagnosis_is_Alzheimers_dementia = cellType__sampleID_tag____2____diagnosis_is_Alzheimers_dementia[cellType__sampleID_tag]
								cov__diagnosis_is_Bipolar_Disorder = cellType__sampleID_tag____2____diagnosis_is_Bipolar_Disorder[cellType__sampleID_tag]
								cov__diagnosis_is_Schizophrenia = cellType__sampleID_tag____2____diagnosis_is_Schizophrenia[cellType__sampleID_tag]
								cov__diagnosis_is_cognitive_impairment = cellType__sampleID_tag____2____diagnosis_is_cognitive_impairment[cellType__sampleID_tag]
								cov_diagnosis_MDD = cellType__sampleID_tag____2____diagnosis_MDD[cellType__sampleID_tag]
								cov_diagnosis_PTSDr = cellType__sampleID_tag____2____diagnosis_PTSD[cellType__sampleID_tag]
								cov_diagnosis_Williams = cellType__sampleID_tag____2____diagnosis_Williams[cellType__sampleID_tag]
								diagnosis_is_control = cellType__sampleID_tag____2____diagnosis_is_control[cellType__sampleID_tag]
								cov__cohort_is_CMC = cellType__sampleID_tag____2____cohort_is_CMC[cellType__sampleID_tag]
								cov__cohort_is_DevBrain = cellType__sampleID_tag____2____cohort_is_DevBrain[cellType__sampleID_tag]
								cov__cohort_is_IsoHuB = cellType__sampleID_tag____2____cohort_is_IsoHuB[cellType__sampleID_tag]
								cov__cohort_is_Kreigstein = cellType__sampleID_tag____2____cohort_is_Kreigstein[cellType__sampleID_tag]
								cov__cohort_is_ROSMAP = cellType__sampleID_tag____2____cohort_is_ROSMAP[cellType__sampleID_tag]
								cov__cohort_is_SZBD = cellType__sampleID_tag____2____cohort_is_SZBD[cellType__sampleID_tag]
								cov__cohort_is_UCLA_ASD = cellType__sampleID_tag____2____cohort_is_UCLA_ASD[cellType__sampleID_tag]
								cov__cohort_is_Urban = cellType__sampleID_tag____2____cohort_is_Urban[cellType__sampleID_tag]
								cov__cohort_Girgenti = cellType__sampleID_tag____2____cohort_Girgenti[cellType__sampleID_tag]
								cov__cohort_LIBD = cellType__sampleID_tag____2____cohort_LIBD[cellType__sampleID_tag]
								cov__cohort_Ma_Sestan = cellType__sampleID_tag____2____cohort_Ma_Sestan[cellType__sampleID_tag]
								cov__cohort_PTSDBrainomics = cellType__sampleID_tag____2____cohort_PTSDBrainomics[cellType__sampleID_tag]
								cov__gt_PC_1 = cellType__sampleID_tag____2____gt_PC_1[cellType__sampleID_tag]
								cov__gt_PC_2 = cellType__sampleID_tag____2____gt_PC_2[cellType__sampleID_tag]
								cov__gt_PC_3 = cellType__sampleID_tag____2____gt_PC_3[cellType__sampleID_tag]
								cov__gt_PC_4 = cellType__sampleID_tag____2____gt_PC_4[cellType__sampleID_tag]
								cov__gt_PC_5 = cellType__sampleID_tag____2____gt_PC_5[cellType__sampleID_tag]
								cov__expr_PC1 = cellType__sampleID_tag____2____expr_PC1[cellType__sampleID_tag]
								cov__expr_PC2 = cellType__sampleID_tag____2____expr_PC2[cellType__sampleID_tag]
								cov__expr_PC3 = cellType__sampleID_tag____2____expr_PC3[cellType__sampleID_tag]
								cov__expr_PC4 = cellType__sampleID_tag____2____expr_PC4[cellType__sampleID_tag]
								cov__expr_PC5 = cellType__sampleID_tag____2____expr_PC5[cellType__sampleID_tag]
								cov__expr_PC6 = cellType__sampleID_tag____2____expr_PC6[cellType__sampleID_tag]
								cov__expr_PC7 = cellType__sampleID_tag____2____expr_PC7[cellType__sampleID_tag]
								cov__expr_PC8 = cellType__sampleID_tag____2____expr_PC8[cellType__sampleID_tag]
								cov__expr_PC9 = cellType__sampleID_tag____2____expr_PC9[cellType__sampleID_tag]
								cov__expr_PC10 = cellType__sampleID_tag____2____expr_PC10[cellType__sampleID_tag]
								cov__expr_PC11 = cellType__sampleID_tag____2____expr_PC11[cellType__sampleID_tag]
								cov__expr_PC12 = cellType__sampleID_tag____2____expr_PC12[cellType__sampleID_tag]
								cov__expr_PC13 = cellType__sampleID_tag____2____expr_PC13[cellType__sampleID_tag]
								cov__expr_PC14 = cellType__sampleID_tag____2____expr_PC14[cellType__sampleID_tag]
								cov__expr_PC15 = cellType__sampleID_tag____2____expr_PC15[cellType__sampleID_tag]
								cov__expr_PC16 = cellType__sampleID_tag____2____expr_PC16[cellType__sampleID_tag]
								cov__expr_PC17 = cellType__sampleID_tag____2____expr_PC17[cellType__sampleID_tag]
								cov__expr_PC18 = cellType__sampleID_tag____2____expr_PC18[cellType__sampleID_tag]
								cov__expr_PC19 = cellType__sampleID_tag____2____expr_PC19[cellType__sampleID_tag]
								cov__expr_PC20 = cellType__sampleID_tag____2____expr_PC20[cellType__sampleID_tag]
								cov__expr_PC21 = cellType__sampleID_tag____2____expr_PC21[cellType__sampleID_tag]
								cov__expr_PC22 = cellType__sampleID_tag____2____expr_PC22[cellType__sampleID_tag]
								cov__expr_PC23 = cellType__sampleID_tag____2____expr_PC23[cellType__sampleID_tag]
								cov__expr_PC24 = cellType__sampleID_tag____2____expr_PC24[cellType__sampleID_tag]
								cov__expr_PC25 = cellType__sampleID_tag____2____expr_PC25[cellType__sampleID_tag]
								cov__expr_PC26 = cellType__sampleID_tag____2____expr_PC26[cellType__sampleID_tag]
								cov__expr_PC27 = cellType__sampleID_tag____2____expr_PC27[cellType__sampleID_tag]
								cov__expr_PC28 = cellType__sampleID_tag____2____expr_PC28[cellType__sampleID_tag]
								cov__expr_PC29 = cellType__sampleID_tag____2____expr_PC29[cellType__sampleID_tag]
								cov__expr_PC30 = cellType__sampleID_tag____2____expr_PC30[cellType__sampleID_tag]
								cov__expr_PC31 = cellType__sampleID_tag____2____expr_PC31[cellType__sampleID_tag]
								cov__expr_PC32 = cellType__sampleID_tag____2____expr_PC32[cellType__sampleID_tag]
								cov__expr_PC33 = cellType__sampleID_tag____2____expr_PC33[cellType__sampleID_tag]
								cov__expr_PC34 = cellType__sampleID_tag____2____expr_PC34[cellType__sampleID_tag]
								cov__expr_PC35 = cellType__sampleID_tag____2____expr_PC35[cellType__sampleID_tag]
								cov__expr_PC36 = cellType__sampleID_tag____2____expr_PC36[cellType__sampleID_tag]
								cov__expr_PC37 = cellType__sampleID_tag____2____expr_PC37[cellType__sampleID_tag]
								cov__expr_PC38 = cellType__sampleID_tag____2____expr_PC38[cellType__sampleID_tag]
								cov__expr_PC39 = cellType__sampleID_tag____2____expr_PC39[cellType__sampleID_tag]
								cov__expr_PC40 = cellType__sampleID_tag____2____expr_PC40[cellType__sampleID_tag]
								cov__expr_PC41 = cellType__sampleID_tag____2____expr_PC41[cellType__sampleID_tag]
								cov__expr_PC42 = cellType__sampleID_tag____2____expr_PC42[cellType__sampleID_tag]
								cov__expr_PC43 = cellType__sampleID_tag____2____expr_PC43[cellType__sampleID_tag]
								cov__expr_PC44 = cellType__sampleID_tag____2____expr_PC44[cellType__sampleID_tag]
								cov__expr_PC45 = cellType__sampleID_tag____2____expr_PC45[cellType__sampleID_tag]
								cov__expr_PC46 = cellType__sampleID_tag____2____expr_PC46[cellType__sampleID_tag]
								cov__expr_PC47 = cellType__sampleID_tag____2____expr_PC47[cellType__sampleID_tag]
								cov__expr_PC48 = cellType__sampleID_tag____2____expr_PC48[cellType__sampleID_tag]
								cov__expr_PC49 = cellType__sampleID_tag____2____expr_PC49[cellType__sampleID_tag]
								cov__expr_PC50 = cellType__sampleID_tag____2____expr_PC50[cellType__sampleID_tag]
								cov__expr_PC51 = cellType__sampleID_tag____2____expr_PC51[cellType__sampleID_tag]
								cov__expr_PC52 = cellType__sampleID_tag____2____expr_PC52[cellType__sampleID_tag]
								cov__expr_PC53 = cellType__sampleID_tag____2____expr_PC53[cellType__sampleID_tag]
								cov__expr_PC54 = cellType__sampleID_tag____2____expr_PC54[cellType__sampleID_tag]
								cov__expr_PC55 = cellType__sampleID_tag____2____expr_PC55[cellType__sampleID_tag]
								cov__expr_PC56 = cellType__sampleID_tag____2____expr_PC56[cellType__sampleID_tag]
								cov__expr_PC57 = cellType__sampleID_tag____2____expr_PC57[cellType__sampleID_tag]
								cov__expr_PC58 = cellType__sampleID_tag____2____expr_PC58[cellType__sampleID_tag]
								cov__expr_PC59 = cellType__sampleID_tag____2____expr_PC59[cellType__sampleID_tag]
								cov__expr_PC60 = cellType__sampleID_tag____2____expr_PC60[cellType__sampleID_tag]
								cov__expr_PC61 = cellType__sampleID_tag____2____expr_PC61[cellType__sampleID_tag]
								cov__expr_PC62 = cellType__sampleID_tag____2____expr_PC62[cellType__sampleID_tag]
								cov__expr_PC63 = cellType__sampleID_tag____2____expr_PC63[cellType__sampleID_tag]
								cov__expr_PC64 = cellType__sampleID_tag____2____expr_PC64[cellType__sampleID_tag]
								cov__expr_PC65 = cellType__sampleID_tag____2____expr_PC65[cellType__sampleID_tag]
								cov__expr_PC66 = cellType__sampleID_tag____2____expr_PC66[cellType__sampleID_tag]
								cov__expr_PC67 = cellType__sampleID_tag____2____expr_PC67[cellType__sampleID_tag]
								cov__expr_PC68 = cellType__sampleID_tag____2____expr_PC68[cellType__sampleID_tag]
								cov__expr_PC69 = cellType__sampleID_tag____2____expr_PC69[cellType__sampleID_tag]
								cov__expr_PC70 = cellType__sampleID_tag____2____expr_PC70[cellType__sampleID_tag]
								cov__expr_PC71 = cellType__sampleID_tag____2____expr_PC71[cellType__sampleID_tag]
								cov__expr_PC72 = cellType__sampleID_tag____2____expr_PC72[cellType__sampleID_tag]
								cov__expr_PC73 = cellType__sampleID_tag____2____expr_PC73[cellType__sampleID_tag]
								cov__expr_PC74 = cellType__sampleID_tag____2____expr_PC74[cellType__sampleID_tag]
								cov__expr_PC75 = cellType__sampleID_tag____2____expr_PC75[cellType__sampleID_tag]
								cov__expr_PC76 = cellType__sampleID_tag____2____expr_PC76[cellType__sampleID_tag]
								cov__expr_PC77 = cellType__sampleID_tag____2____expr_PC77[cellType__sampleID_tag]
								cov__expr_PC78 = cellType__sampleID_tag____2____expr_PC78[cellType__sampleID_tag]
								cov__expr_PC79 = cellType__sampleID_tag____2____expr_PC79[cellType__sampleID_tag]
								cov__expr_PC80 = cellType__sampleID_tag____2____expr_PC80[cellType__sampleID_tag]
								cov__expr_PC81 = cellType__sampleID_tag____2____expr_PC81[cellType__sampleID_tag]
								cov__expr_PC82 = cellType__sampleID_tag____2____expr_PC82[cellType__sampleID_tag]
								cov__expr_PC83 = cellType__sampleID_tag____2____expr_PC83[cellType__sampleID_tag]
								cov__expr_PC84 = cellType__sampleID_tag____2____expr_PC84[cellType__sampleID_tag]
								cov__expr_PC85 = cellType__sampleID_tag____2____expr_PC85[cellType__sampleID_tag]
								cov__expr_PC86 = cellType__sampleID_tag____2____expr_PC86[cellType__sampleID_tag]
								cov__expr_PC87 = cellType__sampleID_tag____2____expr_PC87[cellType__sampleID_tag]
								cov__expr_PC88 = cellType__sampleID_tag____2____expr_PC88[cellType__sampleID_tag]
								cov__expr_PC89 = cellType__sampleID_tag____2____expr_PC89[cellType__sampleID_tag]
								cov__expr_PC90 = cellType__sampleID_tag____2____expr_PC90[cellType__sampleID_tag]
								cov__expr_PC91 = cellType__sampleID_tag____2____expr_PC91[cellType__sampleID_tag]
								cov__expr_PC92 = cellType__sampleID_tag____2____expr_PC92[cellType__sampleID_tag]
								cov__expr_PC93 = cellType__sampleID_tag____2____expr_PC93[cellType__sampleID_tag]
								cov__expr_PC94 = cellType__sampleID_tag____2____expr_PC94[cellType__sampleID_tag]
								cov__expr_PC95 = cellType__sampleID_tag____2____expr_PC95[cellType__sampleID_tag]
								cov__expr_PC96 = cellType__sampleID_tag____2____expr_PC96[cellType__sampleID_tag]
								cov__expr_PC97 = cellType__sampleID_tag____2____expr_PC97[cellType__sampleID_tag]
								cov__expr_PC98 = cellType__sampleID_tag____2____expr_PC98[cellType__sampleID_tag]
								cov__expr_PC99 = cellType__sampleID_tag____2____expr_PC99[cellType__sampleID_tag]
								cov__expr_PC100 = cellType__sampleID_tag____2____expr_PC100[cellType__sampleID_tag]

								expr = cellType__sampleID_tag____2____expr[cellType__sampleID_tag]
								line_to_print = cellType + "\t" + str(genotypeDosage) + "\t" + str(cov__gender_is_female) + "\t" + str(cov__ageDeath) + "\t" + str(cov__diagnosis_is_ASD) + "\t" + str(cov__diagnosis_is_Alzheimers_dementia) + "\t" + str(cov__diagnosis_is_Bipolar_Disorder) + "\t" + str(cov__diagnosis_is_Schizophrenia) + "\t" + str(cov__diagnosis_is_cognitive_impairment) + "\t" + str(cov_diagnosis_MDD) + "\t" + str(cov_diagnosis_PTSDr) + "\t" + str(cov_diagnosis_Williams) + "\t" + str(diagnosis_is_control) + "\t" + str(cov__cohort_is_CMC) + "\t" + str(cov__cohort_is_DevBrain) + "\t" + str(cov__cohort_is_IsoHuB) + "\t" + str(cov__cohort_is_Kreigstein) + "\t" + str(cov__cohort_is_ROSMAP) + "\t" + str(cov__cohort_is_SZBD) + "\t" + str(cov__cohort_is_UCLA_ASD) + "\t" + str(cov__cohort_is_Urban) + "\t" + str(cov__cohort_Girgenti) + "\t" + str(cov__cohort_LIBD) + "\t" + str(cov__cohort_Ma_Sestan) + "\t" + str(cov__cohort_PTSDBrainomics) + "\t" + str(cov__gt_PC_1) + "\t" + str(cov__gt_PC_2) + "\t" + str(cov__gt_PC_3) + "\t" + str(cov__gt_PC_4) + "\t" + str(cov__gt_PC_5) + "\t" + str(cov__expr_PC1) + "\t" + str(cov__expr_PC2) + "\t" + str(cov__expr_PC3) + "\t" + str(cov__expr_PC4) + "\t" + str(cov__expr_PC5) + "\t" + str(cov__expr_PC6) + "\t" + str(cov__expr_PC7) + "\t" + str(cov__expr_PC8) + "\t" + str(cov__expr_PC9) + "\t" + str(cov__expr_PC10) + "\t" + str(cov__expr_PC11) + "\t" + str(cov__expr_PC12) + "\t" + str(cov__expr_PC13) + "\t" + str(cov__expr_PC14) + "\t" + str(cov__expr_PC15) + "\t" + str(cov__expr_PC16) + "\t" + str(cov__expr_PC17) + "\t" + str(cov__expr_PC18) + "\t" + str(cov__expr_PC19) + "\t" + str(cov__expr_PC20) + "\t" + str(cov__expr_PC21) + "\t" + str(cov__expr_PC22) + "\t" + str(cov__expr_PC23) + "\t" + str(cov__expr_PC24) + "\t" + str(cov__expr_PC25) + "\t" + str(cov__expr_PC26) + "\t" + str(cov__expr_PC27) + "\t" + str(cov__expr_PC28) + "\t" + str(cov__expr_PC29) + "\t" + str(cov__expr_PC30) + "\t" + str(cov__expr_PC31) + "\t" + str(cov__expr_PC32) + "\t" + str(cov__expr_PC33) + "\t" + str(cov__expr_PC34) + "\t" + str(cov__expr_PC35) + "\t" + str(cov__expr_PC36) + "\t" + str(cov__expr_PC37) + "\t" + str(cov__expr_PC38) + "\t" + str(cov__expr_PC39) + "\t" + str(cov__expr_PC40) + "\t" + str(cov__expr_PC41) + "\t" + str(cov__expr_PC42) + "\t" + str(cov__expr_PC43) + "\t" + str(cov__expr_PC44) + "\t" + str(cov__expr_PC45) + "\t" + str(cov__expr_PC46) + "\t" + str(cov__expr_PC47) + "\t" + str(cov__expr_PC48) + "\t" + str(cov__expr_PC49) + "\t" + str(cov__expr_PC50) + "\t" + str(cov__expr_PC51) + "\t" + str(cov__expr_PC52) + "\t" + str(cov__expr_PC53) + "\t" + str(cov__expr_PC54) + "\t" + str(cov__expr_PC55) + "\t" + str(cov__expr_PC56) + "\t" + str(cov__expr_PC57) + "\t" + str(cov__expr_PC58) + "\t" + str(cov__expr_PC59) + "\t" + str(cov__expr_PC60) + "\t" + str(cov__expr_PC61) + "\t" + str(cov__expr_PC62) + "\t" + str(cov__expr_PC63) + "\t" + str(cov__expr_PC64) + "\t" + str(cov__expr_PC65) + "\t" + str(cov__expr_PC66) + "\t" + str(cov__expr_PC67) + "\t" + str(cov__expr_PC68) + "\t" + str(cov__expr_PC69) + "\t" + str(cov__expr_PC70) + "\t" + str(cov__expr_PC71) + "\t" + str(cov__expr_PC72) + "\t" + str(cov__expr_PC73) + "\t" + str(cov__expr_PC74) + "\t" + str(cov__expr_PC75) + "\t" + str(cov__expr_PC76) + "\t" + str(cov__expr_PC77) + "\t" + str(cov__expr_PC78) + "\t" + str(cov__expr_PC79) + "\t" + str(cov__expr_PC80) + "\t" + str(cov__expr_PC81) + "\t" + str(cov__expr_PC82) + "\t" + str(cov__expr_PC83) + "\t" + str(cov__expr_PC84) + "\t" + str(cov__expr_PC85) + "\t" + str(cov__expr_PC86) + "\t" + str(cov__expr_PC87) + "\t" + str(cov__expr_PC88) + "\t" + str(cov__expr_PC89) + "\t" + str(cov__expr_PC90) + "\t" + str(cov__expr_PC91) + "\t" + str(cov__expr_PC92) + "\t" + str(cov__expr_PC93) + "\t" + str(cov__expr_PC94) + "\t" + str(cov__expr_PC95) + "\t" + str(cov__expr_PC96) + "\t" + str(cov__expr_PC97) + "\t" + str(cov__expr_PC98) + "\t" + str(cov__expr_PC99) + "\t" + str(cov__expr_PC100) + "\t" + str(expr)
								out_file.write(line_to_print + "\n")
							else:
								print cellType__sampleID_tag + " not in cellType__sampleID_tag____2____expr for " + eGene + "__" + eSNP + "\n"
						else:
							print cellType__sampleID_tag + " not in cellType__sampleID_tag____2____gender_is_female for " + eGene + "__" + eSNP + "\n"
					else:
						print cellType__sampleID_tag + " not in cellType__sampleID_tag____2____genotypeDosage for " + eGene + "__" + eSNP + "\n"
	print "\n\n   -- fininshed w " + eGene + "__" + eSNP + " --   \n"
	out_file.close()
input_file_w_eGene_eSNP_pairs.close()

