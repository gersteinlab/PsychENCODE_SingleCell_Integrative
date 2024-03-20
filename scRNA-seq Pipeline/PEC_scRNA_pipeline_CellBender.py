import sys
import os
import subprocess
import numpy as np
import pandas as pd
import math
import time
import random
from copy import deepcopy
import argparse
import multiprocessing as mp
from functools import partial
import json
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Part 2 of the PEC single-nucleus/single-cell RNA-seq processing pipeline: CellBender')
    parser.add_argument('-J','--inputjson', required=True, help='Input json file for the CellBender runs')
    args = parser.parse_args()

    inputjson = args.inputjson
    with open(inputjson,'r') as infile:
        inputdict = json.load(infile)


    currdir = inputdict['currdir']
    runprefix = inputdict['runprefix']
    transcriptome = inputdict['transcriptome']
    n_jobs_pg = inputdict['n_jobs_pg']
    mito_file = inputdict['mito_file']
    cb_partition = inputdict['cb_partition']
    cb_gpus = inputdict['cb_gpus']
    cb_cpus = inputdict['cb_cpus']
    cb_nodes = inputdict['cb_nodes']
    raw_samplefile_dict = inputdict['raw_samplefile_dict']
    sample_to_batch_map = inputdict['sample_to_batch_map']
    countsdict = inputdict['countsdict']
    hashing_dict = inputdict['hashing']
    hto_dict = inputdict['hto_count']
    if currdir[-1] != "/":
        currdir += "/"

    #### Cellbender remove-background (optional)
    print("checkpoint 2")

    cb_output_dict = {}
    for sampletag, samplefile in raw_samplefile_dict.items():
      samplecounts = countsdict[sampletag]

      cb_dir = f"{currdir}{sampletag}_cb_outputs"
      if not os.path.exists(cb_dir):
          os.mkdir(cb_dir)
      cb_sbatch = f"sbatch --job-name=cb_{sampletag} --time 15:00:00 --partition={cb_partition} --nodes {cb_nodes} --gpus {cb_gpus} --cpus-per-gpu={cb_cpus} --mem-per-cpu=40G cellbender remove-background --input={samplefile} --output={currdir}{sampletag}_cb_outputs/cellbender-output.h5 --expected-cells={samplecounts} --total-droplets-included {samplecounts+20000} --fpr 0.01 --epochs 150 --cuda"
      os.environ['MKL_THREADING_LAYER']='GNU'
      subprocess.call(cb_sbatch,shell=True)
      cb_output = f"{currdir}{sampletag}_cb_outputs/cellbender-output_filtered.h5"
      batchnum = sample_to_batch_map[sampletag]
      if batchnum not in cb_output_dict.keys():
          cb_output_dict[batchnum] = [[sampletag,cb_output]]
      else:
          cb_output_dict[batchnum].append([sampletag,cb_output])
    print("checkpoint 3")
    #### Write json file

    pegasus_input_dict = cb_output_dict

    json_dict = {}
    print("checkpoint 4")
    print(pegasus_input_dict)
    for batchnum, samplelist in pegasus_input_dict.items():
        hto_count = hto_dict[str(batchnum)]
        print(f"HTO_dict = {hto_dict}")
        print(f"HTO_count = {hto_count}")
        hto_count = [hto for hto in hto_count if hto != None]
        hashing = hashing_dict[str(batchnum)]
        if "True" in hashing:
            hashing = "True"
        else:
            hashing = "False"
        pg_dict = {"sampledetails":inputdict['sampledetails'],"matrix_directory": samplelist, "mito_file": mito_file, "transcriptome": transcriptome, "n_jobs_pg": n_jobs_pg, "hashing": hashing, "currdir": currdir}
        if hashing == "True":
            pg_dict["hto_file"] = hto_count[0]

        ###Set default parameters
        if "qc_min_umis" in inputdict.keys():
            pg_dict["qc_min_umis"] = inputdict["qc_min_umis"]
        if "qc_percent_mito" in inputdict.keys():
            pg_dict["qc_percent_mito"] = inputdict["qc_percent_mito"]
        if "qc_min_genes" in inputdict.keys():
            pg_dict["qc_min_genes"] = inputdict["qc_min_genes"]
        if "dd_bst_n_iters" in inputdict.keys():
            pg_dict["dd_bst_n_iters"] = inputdict["dd_bst_n_iters"]
        if "dd_bst_use_pheno" in inputdict.keys():
            pg_dict["dd_bst_use_pheno"] = inputdict["dd_bst_use_pheno"]
        if "dd_bst_std_scaling" in inputdict.keys():
            pg_dict["dd_bst_std_scaling"] = inputdict["dd_bst_std_scaling"]
        if "dd_pred_pthresh" in inputdict.keys():
            pg_dict["dd_pred_pthresh"] = inputdict["dd_pred_pthresh"]
        if "dd_pred_voterthresh" in inputdict.keys():
            pg_dict["dd_pred_voterthresh"] = inputdict["dd_pred_voterthresh"]
        if "hvg_n_top" in inputdict.keys():
            pg_dict["hvg_n_top"] = inputdict["hvg_n_top"]
        jsonString = json.dumps(pg_dict)
        jsonFileName = f"{runprefix}_Batch{batchnum}_pg_input.json"
        jsonFile = open(jsonFileName, "w")
        jsonFile.write(jsonString)
        jsonFile.close()
        json_dict[batchnum] = jsonFileName
    print("checkpoint 5")
    print(json_dict)
