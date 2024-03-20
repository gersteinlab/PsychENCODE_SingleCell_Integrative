import numpy as np
import sys
import os
import random
from operator import itemgetter, attrgetter, methodcaller
# import matplotlib.pyplot as plt
# import os.path

'''
python retain_cellTypes.py ret_cellTypes___OR2T10__chr1:248678345:G:A.dat OR2T10__chr1:248678345:G:A__input_for_mxd_fx_model_sanit_COV.dat > OR2T10__chr1:248678345:G:A__input_for_mxd_fx_model_sanit_COV_retCellTypes.dat

'''


list_of_cellTypes_to_ret = list()
file_w_cellTypes_to_ret = open(sys.argv[1], "r")
for line in file_w_cellTypes_to_ret:
	cellType = line.strip()
	list_of_cellTypes_to_ret.append(cellType)
file_w_cellTypes_to_ret.close()


orig_file = open(sys.argv[2], "r")
for line in orig_file:
	ln_elms = list()
	ln_elms = line.split()
	cellType = ln_elms[0]
	if (cellType in list_of_cellTypes_to_ret)   or   (ln_elms[0] == "cellType"):
		print(line.strip())
orig_file.close()

