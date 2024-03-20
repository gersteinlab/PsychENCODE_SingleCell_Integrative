import pandas as pd
import math
import numpy as np
import sys
import os
import random
from operator import itemgetter, attrgetter, methodcaller


df = pd.read_csv("ps_bulk_umi_matrs/pseudo_bulk__Br2743-annotated_matrix.txt", sep="\t", index_col="featurekey")
col_sum = df.sum(axis=0)
df = df.divide(col_sum)
df = df.multiply(1000000.0)
df = np.log2(df.add(1))
df.to_csv('ps_bulk_umi_matrs/cpm_vanilla_normalized/pseudo_bulk__Br2743-annotated_matrix_CPM_norm.txt', index=True, sep="\t")


