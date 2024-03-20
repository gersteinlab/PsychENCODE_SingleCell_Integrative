import numpy as np
import sys
import os
import random
from operator import itemgetter, attrgetter, methodcaller
# import matplotlib.pyplot as plt
# import os.path

'''
USAGE:
python det_nominal_thresholds.py /Users/dc/Desktop/rsch/pec2/v5_scqtls/ld_pruning_plink/output_10000_permuts/sig_QTLs_fdr_0.05/

'''


permutations_dir = sys.argv[1]
permutations_dir_wheel = os.listdir(permutations_dir)
cell_type___2___pt = {}
for file_name in permutations_dir_wheel:
	## Get cell_type
	cell_type = "NULL"
	if "100_expr_PCs" in file_name:
		cell_type = file_name[13:-37]
	elif "20_expr_PCs" in file_name:
		cell_type = file_name[13:-36]
	#print "cell type is: " + cell_type


	## Open associated permutations file and determine pt val (ie, least sig p_adj val)
	full_file_name = permutations_dir + "/" + file_name
	permutations_file = open(full_file_name, "r")
	max_p_adj = -1000
	for line in permutations_file:
		ln_elms = list()
		ln_elms = line.split()
		p_adj = float(ln_elms[19])
		if p_adj > max_p_adj:
			max_p_adj = p_adj
	permutations_file.close()
	pt = max_p_adj
	cell_type___2___pt[cell_type] = pt
	print cell_type + "\t" + str(cell_type___2___pt[cell_type])

	## For each significant eGene -- print an R command to determine the nominal p-val threshold
	permutations_file = open(full_file_name, "r")
	out_file_name = cell_type + "__det_gene_level_nominal_p_val_threshs.r"
	out_file = open(out_file_name, "w")
	list_of_sig_eGenes = list()
	for line in permutations_file:
		ln_elms = list()
		ln_elms = line.split()
		eGene = ln_elms[0]
		list_of_sig_eGenes.append(eGene)
		alpha = float(ln_elms[13])
		beta = float(ln_elms[14])
		#print "alpha = " + str(alpha) + "\t\t beta = " + str(beta)
		#print file_name + "\t" + str(alpha) + "\t" + str(beta) + "\t" + str(pt) + "\t" + cell_type
		qbeta_line = "nominal_p_val_thresh = qbeta(" + str(pt) + ", " + str(alpha) + ", " + str(beta) + ", ncp = 0, lower.tail = TRUE, log.p = FALSE)"
		out_file.write("eGene = \"" + eGene + "\"" + "\n")
		out_file.write(qbeta_line + "\n")
		out_file.write("line_to_print = paste(eGene, nominal_p_val_thresh)" + "\n")
		out_file.write("print(line_to_print, quote=FALSE)" + "\n\n")
	permutations_file.close()
	out_file.close()




print "\n\n -- fin -- \n"

