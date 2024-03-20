import numpy as np
import sys
import os
import random
from operator import itemgetter, attrgetter, methodcaller
# import matplotlib.pyplot as plt
# import os.path

'''
USAGE:
python pull_sig_QTLs.py /ysm-gpfs/pi/gerstein/dc547/brain/v5_scqtls/nominal_pass_to_call_QTLs/nominals_out/ /ysm-gpfs/pi/gerstein/dc547/brain/v5_scqtls/sig_QTLs/

'''

nominals_dir = sys.argv[1]
out_dir_to_store_sig_QTLs = sys.argv[2]
nominals_dir_wheel = os.listdir(nominals_dir)
for file_name in nominals_dir_wheel:
	##  Get cell_type
	cell_type = file_name[9:-4]
	print file_name + "\t" + cell_type

	##  Build dict gene___2___nom_p_val_thresh
	gene___2___nom_p_val_thresh = {}
	file_name_w_gene_level_nominal_p_val_threshs = cell_type + "__gene_level_nominal_p_val_threshs.dat"
	file_w_gene_level_nominal_p_val_threshs = open(file_name_w_gene_level_nominal_p_val_threshs, "r")
	list_of_sig_eGenes = list()
	for line in file_w_gene_level_nominal_p_val_threshs:
		ln_elms = list()
		ln_elms = line.split()
		gene = ln_elms[1]
		list_of_sig_eGenes.append(gene)
		nominal_p_val_thresh = float(ln_elms[2])
		gene___2___nom_p_val_thresh[gene] = nominal_p_val_thresh

	##  Identify significant eQTLs
	full_file_name = nominals_dir + "/" + file_name
	nominals_file = open(full_file_name, "r")
	out_file_name = out_dir_to_store_sig_QTLs + "/" + cell_type + "_sig_QTLs.dat"
	out_file = open(out_file_name, "w")
	for line in nominals_file:
		ln_elms = list()
		ln_elms = line.split()
		gene = ln_elms[0]
		if gene in list_of_sig_eGenes:
			nominal_p_val = float(ln_elms[11])
			if nominal_p_val <= gene___2___nom_p_val_thresh[gene]:
				out_file.write(line)
	nominals_file.close()
	out_file.close()


print "\n\n -- fin -- \n"

