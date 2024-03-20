import numpy as np
import sys
import os
import random
from operator import itemgetter, attrgetter, methodcaller
# import matplotlib.pyplot as plt
# import os.path

'''
python retain_covariates.py ret_covariates___OR2T10__chr1:248678345:G:A.dat OR2T10__chr1:248678345:G:A__input_for_mxd_fx_model.dat > OR2T10__chr1:248678345:G:A__input_for_mxd_fx_model_sanit_COV.dat

'''

list_of_covariates_to_ret = list()
file_w_covariates_to_ret = open(sys.argv[1], "r")
for line in file_w_covariates_to_ret:
	covariate = line.strip()
	list_of_covariates_to_ret.append(int(covariate))
file_w_covariates_to_ret.close()


list_of_col_indeces_to_retain = list()
list_of_col_indeces_to_retain.append(0)  ## retain cellType
list_of_col_indeces_to_retain.append(1)  ## retain genotypeDosage
list_of_col_indeces_to_retain.append(130)  ## retain expression
orig_file = open(sys.argv[2], "r")
for line in orig_file:
	ln_elms = list()
	ln_elms = line.split()
	if ln_elms[0] == "cellType": ## at header line
		col_index = 2
		while col_index < len(ln_elms):
			curr_cov_index = col_index - 1
			if curr_cov_index in list_of_covariates_to_ret:
				list_of_col_indeces_to_retain.append(col_index)
			col_index += 1
orig_file.close()
'''
for i in list_of_col_indeces_to_retain:
	print str(i)
'''


orig_file = open(sys.argv[2], "r")
for line in orig_file:
	ln_elms = list()
	ln_elms = line.split()
	line_to_print = ln_elms[0]
	col_index = 1
	while col_index < len(ln_elms):
		if col_index in list_of_col_indeces_to_retain:
			line_to_print = line_to_print + "\t" + ln_elms[col_index].strip()
		col_index += 1
	print(line_to_print)
orig_file.close()

