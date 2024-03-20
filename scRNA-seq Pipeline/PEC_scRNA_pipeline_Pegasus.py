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
    parser = argparse.ArgumentParser(description='Part 3 of the PEC single-nucleus/single-cell RNA-seq processing pipeline: Pegasus')
    parser.add_argument('-J','--inputjson', required=True, help='Input json file for all the runs')
    args = parser.parse_args()

    inputjson = args.inputjson
    with open(inputjson,'r') as infile:
        inputdict = json.load(infile)

    currdir = inputdict['currdir']
    runprefix = inputdict['runprefix']
    transcriptome = inputdict['transcriptome']
    mito_file = inputdict['mito_file']
    n_jobs_pg = inputdict['n_jobs_pg']

    if currdir[-1] != "/":
        currdir += "/"

    batchdict = {}
    for sampledict in inputdict['sampledetails']:
        fastqfolder = sampledict['fastqfolder']
        samplefolder = sampledict['samplefoldername']
        sampletag = sampledict['sampletag']
        batchnum = sampledict['batchnum']
        chemistry = sampledict['chemistry']
        hashing = sampledict['hashing']
        if hashing == "True":
            hashingfastq1 = sampledict['hashingfastq1']
            hashingfastq2 = sampledict['hashingfastq2']
            if 'tags' in sampledict.keys():
                tags = sampledict['tags']
            if 'whitelist' in sampledict.keys():
                whitelist = sampledict['whitelist']
        if batchnum not in batchdict.keys():
            batchdict[batchnum] = [(sampletag,currdir+samplefolder,chemistry)]
        else:
            batchdict[batchnum].append((sampletag,currdir+samplefolder,chemistry))

    #### Call Pegasus script
    for batchnum in batchdict.keys():
        pegasus_call = f"python3 Pegasus-Pipeline.py -J={runprefix}_Batch{batchnum}_pg_input.json -S={currdir}{runprefix}_Batch{batchnum}_pg_outputs"
        subprocess.call(pegasus_call,shell=True)
    print("checkpoint 6")
