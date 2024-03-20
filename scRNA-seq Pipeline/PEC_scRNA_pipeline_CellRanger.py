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

def run_cellranger(transcriptome,numproc,sampledict):
    fastqfolder = sampledict['fastqfolder']
    samplefolder = sampledict['samplefoldername']
    sampletag = sampledict['sampletag']
    cellflag = 'False'
    if 'cells' in sampledict.keys():
        cellflag = sampledict['cells']

    #### CellRanger count run(s)
    #### 1. RNA-seq count matrix run
    print(f"Running CellRanger count for the RNA-seq count matrix for sample {sampletag}")
    if cellflag == 'True':
        count_call = f"cellranger count  --id={samplefolder} --localcores={numproc} --fastqs={fastqfolder} --sample={sampletag} --transcriptome={transcriptome} --jobmode=slurm.template"
    else:
        count_call = f"cellranger count  --id={samplefolder} --localcores={numproc} --include-introns --fastqs={fastqfolder} --sample={sampletag} --transcriptome={transcriptome} --jobmode=slurm.template"
    print(count_call)
    subprocess.call(count_call,shell=True)
    print(f"Finished submitting CellRanger count for the RNA-seq count matrix for sample {sampletag}")
    output = f"{samplefolder}"

    return output

#### 2. CellHashing and MULTI-seq: runs for antibody counts
def run_citeseqcount(currdir,counts_by_batch,numproc,sampledict):
    fastqfolder = sampledict['fastqfolder']
    samplefolder = sampledict['samplefoldername']
    sampletag = sampledict['sampletag']
    batchnum = sampledict['batchnum']
    chemistry = sampledict['chemistry']
    hashing = sampledict['hashing']

    if chemistry == "v2":
        umil = 26
    if chemistry == "v3":
        umil = 28
    if hashing == "True":
        hashingfastq1 = sampledict['hashingfastq1']
        hashingfastq2 = sampledict['hashingfastq2']
        tags = sampledict['tags']
        print(hashingfastq1)
        print(hashingfastq2)
        print(f"Running CellRanger count for the Antibody count matrix for sample {sampletag}")
        hashing_call = f"CITE-seq-Count -R1={hashingfastq1} -R2={hashingfastq2} -t={tags} -cbf 1 -cbl 16 -umif 17 -umil {umil} -cells {counts_by_batch[batchnum]} -T {numproc} -o {currdir+samplefolder}/CITE-seq-Results"
        print(hashing_call)
        #if not os.path.exists(f"{currdir+samplefolder}/CITE-seq-Results"):
        if "Set5" in samplefolder:
            subprocess.call(hashing_call,shell=True)
        hto_count = f"{currdir+samplefolder}/CITE-seq-Results/umi_count"
        print(f"Finished running CellRanger count for the Antibody count matrix for sample {sampletag}")
    else:
        hto_count = None

    return hto_count
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Part 1 of the PEC single-nucleus/single-cell RNA-seq processing pipeline: CellRanger count, aggr and CellHashing')
    parser.add_argument('-J','--inputjson', required=True, help='Input json file for all the runs')
    parser.add_argument('-N','--nprocesses', required=False, type=int, help='Number of processes for the wrapper script',default = 1)
    args = parser.parse_args()

    inputjson = args.inputjson
    ncellr = args.nprocesses
    with open(inputjson,'r') as infile:
        inputdict = json.load(infile)

    currdir = inputdict['currdir']
    runprefix = inputdict['runprefix']
    transcriptome = inputdict['transcriptome']
    mito_file = inputdict['mito_file']
    numproc = inputdict['numproc']
    n_jobs_pg = inputdict['n_jobs_pg']
    runcellbender = inputdict['runcellbender']
    if runcellbender=="True":
        if 'cb_partition' in inputdict.keys():
            cb_partition = inputdict['cb_partition']
        else:
            cb_partition = "general"
        if 'cb_gpus' in inputdict.keys():
            cb_gpus = inputdict['cb_gpus']
        else:
            cb_gpus = 1
        if 'cb_cpus' in inputdict.keys():
            cb_cpus = inputdict['cb_cpus']
        else:
            cb_cpus = 1
        if 'cb_nodes' in inputdict.keys():
            cb_nodes = inputdict['cb_nodes']
        else:
            cb_nodes = 1

    if currdir[-1] != "/":
        currdir += "/"

    pool = mp.Pool(processes=ncellr)
    fix_cellr=partial(run_cellranger,transcriptome,numproc)
    count_output_list = pool.map(fix_cellr, inputdict['sampledetails'])
    pool.close()
    pool.join()

    for output in count_output_list:
        move_outputs = f"mv {output} {currdir}"
        subprocess.call(move_outputs,shell=True)

    raw_samplefile_dict = {}
    batchdict = {}
    countsdict = {}
    hashingdict = {}
    htodict = {}
    sample_to_batch_map = {}
    counts_by_batch = {}
    for sampledict in inputdict['sampledetails']:
        fastqfolder = sampledict['fastqfolder']
        samplefolder = sampledict['samplefoldername']
        sampletag = sampledict['sampletag']
        batchnum = sampledict['batchnum']
        chemistry = sampledict['chemistry']
        hashing = sampledict['hashing']

        raw_samplefile = f"{currdir+samplefolder}/outs/raw_feature_bc_matrix.h5"


        if runcellbender=="True":
            cellb_path = f"{currdir}{sampletag}_cb_outputs"
            if not os.path.exists(cellb_path):
                os.mkdir(cellb_path)


        df_counts = pd.read_csv(f"{currdir+samplefolder}/outs/metrics_summary.csv")
        median_umis_per_cell = df_counts['Median UMI Counts per Cell'].iloc[0]
        median_umis_per_cell = str(median_umis_per_cell)
        median_umis_per_cell = median_umis_per_cell.replace(',',"")
        median_umis_per_cell = int(median_umis_per_cell)

        if median_umis_per_cell < 1000:
            print(f"Likely Low-quality sample detected! Sample {sampletag} has only {median_umis_per_cell} Median UMIs per cell. Eliminating from further analysis.")
            continue
        else:
            sample_to_batch_map[sampletag] = batchnum
            raw_samplefile_dict[sampletag] = raw_samplefile
            if batchnum not in batchdict.keys():
                batchdict[batchnum] = [(sampletag,currdir+samplefolder,chemistry)]
            else:
                batchdict[batchnum].append((sampletag,currdir+samplefolder,chemistry))
        expected_cells = df_counts.iloc[0,0]
        expected_cells = str(expected_cells)
        expected_cells = expected_cells.replace(',',"")
        expected_cells = int(expected_cells)
        countsdict[sampletag] = expected_cells
        if batchnum not in counts_by_batch.keys():
            counts_by_batch[batchnum] = expected_cells+10000
        else:
            counts_by_batch[batchnum] += expected_cells+10000
        if batchnum not in hashingdict.keys():
            hashingdict[batchnum] = [hashing]
        else:
            hashingdict[batchnum].append(hashing)

    pool = mp.Pool(processes=ncellr)
    fix_csc=partial(run_citeseqcount,currdir,counts_by_batch,numproc)
    hto_count_list = pool.map(fix_csc, inputdict['sampledetails'])
    pool.close()
    pool.join()

    for i, sampledict in enumerate(inputdict['sampledetails']):
        hto_count = hto_count_list[i]
        batchnum = sampledict['batchnum']
        if batchnum not in htodict.keys():
            htodict[batchnum] = [hto_count]
        else:
            htodict[batchnum].append(hto_count)
    print(countsdict)
    if runcellbender!="True":
        filtered_samplefile_dict = {}
        for batchnum, sampletuples in batchdict.items():
            if len(sampletuples) > 1:
                # #### Write input for CellRanger aggr
                batchid = f"batch_{batchnum}_aggr"
                batchmetafile = f"batch_{batchnum}_aggr.csv"
                with open(batchmetafile,"w") as outfile:
                    outfile.write("sample_id,molecule_h5,batch\n")
                    for samplecase in sampletuples:
                        outfile.write(",".join([samplecase[0],f"{samplecase[1]}/outs/molecule_info.h5",f"{samplecase[2]}_lib"])+"\n")
                #### CellRanger aggr
                print(f"Running CellRanger aggr for the RNA-seq count matrices for batch {batchnum}")
                aggr_call = f"cellranger aggr --id={batchid} --csv={batchmetafile} --normalize=none"
                subprocess.call(aggr_call,shell=True)
                print(f"Finished running CellRanger count for the RNA-seq count matrices for batch {batchnum}")

                filtered_samplefile = f"{currdir+batchid}/outs/count/filtered_feature_bc_matrix.h5"
            else:

                filtered_samplefile = f"{currdir+samplefolder}/outs/filtered_feature_bc_matrix.h5"

            filtered_samplefile_dict[batchnum] = filtered_samplefile


        pegasus_input_dict = filtered_samplefile_dict
        json_dict = {}

        for batchnum, sampleh5 in pegasus_input_dict.items():
            hto_count = htodict[batchnum]
            hto_count = [hto for hto in hto_count if hto != None]
            hashing = hashingdict[batchnum]
            if "True" in hashing:
                hashing = "True"
            else:
                hashing = "False"
            pg_dict = {"sampledetails":inputdict['sampledetails'],"matrix_directory": [sampleh5], "mito_file": mito_file, "transcriptome": transcriptome, "n_jobs_pg": n_jobs_pg, "hashing": hashing, "hto_file": hto_count[0], "currdir": currdir}
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
    else:
        cellbender_dict = {}
        cellbender_dict['currdir'] = currdir
        cellbender_dict['runprefix'] = runprefix
        cellbender_dict['sampledetails'] = inputdict['sampledetails']
        cellbender_dict['n_jobs_pg'] = n_jobs_pg
        cellbender_dict['mito_file'] = mito_file
        cellbender_dict['transcriptome'] = transcriptome
        cellbender_dict['cb_partition'] = cb_partition
        cellbender_dict['cb_gpus'] = cb_gpus
        cellbender_dict['cb_cpus'] = cb_cpus
        cellbender_dict['cb_nodes'] = cb_nodes
        cellbender_dict['raw_samplefile_dict'] = raw_samplefile_dict
        cellbender_dict['sample_to_batch_map'] = sample_to_batch_map
        cellbender_dict['countsdict'] = countsdict
        cellbender_dict['hashing'] = hashingdict
        cellbender_dict['hto_count'] = htodict
        jsonString = json.dumps(cellbender_dict)
        jsonFileName = f"{runprefix}_CellBender_input.json"
        jsonFile = open(jsonFileName, "w")
        jsonFile.write(jsonString)
        jsonFile.close()
