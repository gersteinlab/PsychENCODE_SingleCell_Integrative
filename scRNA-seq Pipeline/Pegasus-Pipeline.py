import numpy as np
import pandas as pd
import pegasus as pg
import math
import matplotlib
import seaborn as sns
import json
import csv
import re
import doubletdetection as dd
import scrublet as scr

from pegasusio import UnimodalData, MultimodalData

import sys
import os
import subprocess
import time
import random
import argparse
import multiprocessing as mp

#############################################################################################
if __name__=="__main__":
    ###Setting up program parameters annd optionns available to the user
    parser = argparse.ArgumentParser()
    parser.add_argument('-J','--jsonfile', required=True, help='')
    parser.add_argument('-S','--samplename', required=True, help='Name of the sample to be processed (will be designated as the folder name containing all files)')

    args = parser.parse_args()

    ###Create local variables
    jsonfile = args.jsonfile
    samplename = args.samplename

    batchname = jsonfile.split("_")[-3]
    ###Create directory for outputs
    if not os.path.exists(samplename):
        os.mkdir(samplename)

    ###Read in jsonfile
    f = open(jsonfile,)
    jdict = json.load(f)

    ###Set default parameters
    if "qc_min_umis" not in jdict.keys():
        jdict["qc_min_umis"] = 500
    if "qc_percent_mito" not in jdict.keys():
        jdict["qc_percent_mito"] = 10
    if "qc_min_genes" not in jdict.keys():
        jdict["qc_min_genes"] = 200
    if "dd_bst_n_iters" not in jdict.keys():
        jdict["dd_bst_n_iters"] = 25
    if "dd_bst_use_pheno" not in jdict.keys():
        jdict["dd_bst_use_pheno"] = False
    if "dd_bst_std_scaling" not in jdict.keys():
        jdict["dd_bst_std_scaling"] = True
    if "dd_pred_pthresh" not in jdict.keys():
        jdict["dd_pred_pthresh"] = 1e-16
    if "dd_pred_voterthresh" not in jdict.keys():
        jdict["dd_pred_voterthresh"] = 0.3
    if "hvg_n_top" not in jdict.keys():
        jdict["hvg_n_top"] = 5000
    if "n_jobs" not in jdict.keys():
        jdict["n_jobs"] = 10

    print(jdict)

    currdir = jdict["currdir"]

    ###Create summary stats text file
    summary_file = open(f"{samplename}/{batchname}_summary_stats.txt","w")
    summary_file.write("Parameters used:\n")
    summary_file.write(json.dumps(jdict))
    summary_file.write("\n")

    ###Write out pg.aggregate csv file
    header = ["Sample","Location"]
    filename = f"{batchname}_pg_aggregate.csv"
    csvfile = open(f"{samplename}/{batchname}_pg_aggregate.csv", 'w')
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(header)


    for dataset in jdict["matrix_directory"]:
        print("Importing count matrix")
        data = pg.read_input(dataset[1])
        summary_file.write(f"\nSize of count matrix {dataset[0]} (# of obs, # of genes):"+str(data.X.get_shape()))
        summary_file.write("\n")
        ###Read count matrix if only one sample
        # if len(jdict["matrix_directory"]) == 1:
        #     print("Importing count matrix")
        #     data = pg.read_input(jdict["matrix_directory"][0][1])
        #     summary_file.write("\nSize of count matrix (# of obs, # of genes):"+str(data.X.get_shape()))
        #     summary_file.write("\n")

        ###Read count matrix if aggregate needed (more than one sample)
        # else:




        ###QC Metrics and filtration and log-normalization
        print("Beginning QC metrics")
        pg.qc_metrics(data, min_umis=jdict["qc_min_umis"], percent_mito=jdict["qc_percent_mito"], min_genes=jdict["qc_min_genes"])
        df_qc = pg.get_filter_stats(data)
        print(data)
        summary_file.write("\nQC metrics stats:\n")
        summary_file.write(df_qc.to_string(header = True, index = True))
        pg.filter_data(data)
        print(data)

        ###Filter out mitochondrial genes
        print("Beginning mito gene filtration")
        mito_df = pd.read_csv(jdict["mito_file"])
        mito_df = mito_df.loc[0:1135,'HumanGeneID':'Symbol']
        mito_list = []
        for i in range(mito_df.shape[0]):
            mito_list.append(mito_df.loc[i,'Symbol'])
            mito_list = set(mito_list)
            mito_list = list(mito_list)
        non_mito_list = []
        for i in data.var_names:
            if i in mito_list:
                non_mito_list.append(False)
            else:
                non_mito_list.append(True)
        data_subset = data[:, non_mito_list].copy()
        data = data_subset
        data = MultimodalData(data)
        print(data)
        summary_file.write("\nSize of count matrix post mito gene filtration:"+str(data.X.get_shape()))
        summary_file.write("\n")

        pg.identify_robust_genes(data)
        pg.log_norm(data)
        summary_file.write("\n")

        ###Demultiplexing
        print("Beginning demultiplexing")
        if jdict["hashing"] == "True":
            print(f"HTO data = {jdict['hto_file']}")

            features_file = jdict["hto_file"]+"/features.tsv.gz"
            feature_metadata = pd.read_csv(features_file, sep="\t", header=None)
            print(f"Feature metadata = {feature_metadata}")
            print(f"Feature metadata shape = {feature_metadata.shape}")
            # feature_metadata.iloc[:, 0] = feature_metadata.iloc[:, 0].str.replace("_","")
            # feature_metadata.iloc[:, 0] = feature_metadata.iloc[:, 0].str.replace("-","_")
            # feature_metadata.iloc[~feature_metadata.iloc[:,0].str.contains("_"), 0] = feature_metadata[~feature_metadata.iloc[:,0].str.contains("_")].astype('str') + '_0'
            # print(feature_metadata)
            # features_updated = jdict["hto_file"]+"/features_clean.tsv.gz"
            # feature_metadata.to_csv(features_updated,sep="\t",index=False,compression="gzip", header = False)
            #
            # rm_call = f"rm {features_file}"
            # subprocess.call(rm_call,shell=True)
            #
            # rename_call = f"mv {features_updated} {features_file}"
            # subprocess.call(rename_call,shell=True)


            hto_data = pg.read_input(jdict["hto_file"]+"/matrix.mtx.gz", genome = 'hashing_HTO', modality="hashing")
            features = pd.read_csv(jdict["hto_file"]+"/features.tsv.gz", header=None)
            barcodes = pd.read_csv(jdict["hto_file"]+"/barcodes.tsv.gz", header=None)
            hto_data.var_names = features[0]
            hto_data.obs_names = barcodes[0]

            pg.estimate_background_probs(hto_data)
            print(hto_data.uns["background_probs"])
            pg.demultiplex(data, hto_data)
            data_subset = data[data.obs["demux_type"] == "singlet",:].copy()
            data = data_subset
            data = MultimodalData(data)
            print(data)

            summary_file.write("\nSize of count matrix post hashing:"+str(data.X.get_shape()))
            summary_file.write("\n")

        ###Doublet detection -- Scrublet
        print("Beginning doublet detection - scrublet")
        summary_file.write("\nDoublet detection and filtration – Scrublet:")
        data.select_matrix('raw.X')
        counts_matrix = data.X
        scrub = scr.Scrublet(counts_matrix)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        doublet = scrub.predicted_doublets_
        if doublet is not None:
            data.obs["doublet"] = doublet
            print(doublet)
            data_subset = data[data.obs["doublet"] == False,:].copy()
            data = data_subset
            data = MultimodalData(data)

            summary_file.write("\nSize of count matrix post doublet filtration:"+str(data.X.get_shape()))
            summary_file.write("\n")

        ###Doublet detection -- Doublet Detection
        print("Beginning doublet detection - DD")
        summary_file.write("\nDoublet detection and filtration – Doublet Detection:")
        clf = dd.BoostClassifier(n_iters=jdict["dd_bst_n_iters"], clustering_algorithm="phenograph", standard_scaling=jdict["dd_bst_std_scaling"])
        data.select_matrix('raw.X')
        print(f"after select raw = {np.max(data.X.T.todense())}")
        doublets = clf.fit(data.X).predict(p_thresh=jdict["dd_pred_pthresh"], voter_thresh=jdict["dd_pred_voterthresh"])
        doublet_score = clf.doublet_score()
        data.obs["doublet"] = doublets
        data.obs["doublet_score"] = doublet_score
        print(doublets)
        data_subset = data[data.obs["doublet"] == 0,:].copy()
        data = data_subset
        data = MultimodalData(data)
        data.select_matrix('X')
        #data_TPM_norm = data_TPM.copy()
        #data_TPM_norm.X = (10**6)*normalize(data_TPM.X,norm='l1',axis=1)
        #print(data_TPM_norm.var['featureid'])
        summary_file.write("\nSize of count matrix post doublet filtration:"+str(data.X.get_shape()))
        summary_file.write("\n")
        ###Aggregate count matrices
        dataset_anndata = f"{samplename}/{dataset[0]}.h5ad"
        pg.write_output(data,dataset_anndata)
        csvwriter.writerow([dataset[0],dataset_anndata])
        print(data)

    csvfile.close()
    ###Aggregate count matrices
    print("Aggregating count matrices")
    print(f"{samplename}/{batchname}_pg_aggregate.csv")
    data = pg.aggregate_matrices(f"{samplename}/{batchname}_pg_aggregate.csv")
    summary_file.write("\nSize of aggregated count matrix (# of obs, # of genes):"+str(data.X.get_shape()))
    summary_file.write("\n")
    print(data)
    pg.identify_robust_genes(data)
    ###UMAP pre-Harmony
    data_pre = data.copy()

    pg.highly_variable_features(data_pre, batch="Channel", n_top=jdict["hvg_n_top"])
    pg.pca(data_pre)
    pg.neighbors(data_pre, n_jobs = jdict["n_jobs"])
    pg.leiden(data_pre)
    pg.umap(data_pre, n_jobs = jdict["n_jobs"])
    ###save UMAP figure
    umap_fig_pre = pg.scatter(data_pre, attrs=['leiden_labels', 'Channel'], basis='umap', return_fig=True)
    umap_fig_pre.savefig(f"{samplename}/umap_fig_pre.png")

    ###HVG, PCA, Harmony, Neighbors, Leiden, UMAP
    pg.highly_variable_features(data, batch="Channel", n_top=jdict["hvg_n_top"])
    pg.pca(data)
    pca_key = pg.run_harmony(data, n_jobs = jdict["n_jobs"])
    pg.neighbors(data, rep=pca_key, n_jobs = jdict["n_jobs"])
    pg.leiden(data, rep=pca_key)
    pg.umap(data, rep=pca_key, n_jobs = jdict["n_jobs"])
    print(data)
    ###save UMAP figure
    umap_fig_post = pg.scatter(data, attrs=['leiden_labels', 'Channel'], basis='umap', return_fig=True)
    umap_fig_post.savefig(f"{samplename}/umap_fig_post.png")
    ###save UMAP coordinates
    umap_coord = data.obsm['X_umap']
    first_coord = []
    second_coord = []
    count = 0
    for i in umap_coord:
        for j in i:
            if count % 2 == 0:
                first_coord.append(j)
            count += 1
    count_2 = 0
    for i in umap_coord:
        for j in i:
            if count_2 % 2 == 1:
                second_coord.append(j)
            count_2 += 1
    first_coord = pd.Series(first_coord)
    second_coord = pd.Series(second_coord)
    d = {'barcodekey':data.obs_names,'first_coord':first_coord, 'second_coord':second_coord}
    umap_df = pd.DataFrame(d)
    umap_df.to_csv(f"{samplename}/umap_coords.csv", index=False)
    ###save summary stats on cluster sizes
    summary_file.write("\nCluster sizes:\n")
    summary_file.write(pd.DataFrame(data.obs[['leiden_labels']].value_counts()).to_string(header = False, index = True))
    summary_file.write("\n")

    ###Marker Gene Analysis
    pg.de_analysis(data, cluster='leiden_labels', t=True)
    print("de analysis done")
    marker_dict = pg.markers(data)
    print("marker gene analysis done")
    ###creating up and down regulated dataframes
    master_list_up = []
    master_list_down = []
    for keys in marker_dict:
        value = marker_dict[keys]
        for j in value:
            if j == 'up':
                df_value = marker_dict[keys]['up']
                master_list_up.append(df_value.index)
    for keys in marker_dict:
        value = marker_dict[keys]
        for j in value:
            if j == 'down':
                df_value_d = marker_dict[keys]['down']
                master_list_down.append(df_value_d.index)
    up_marker_df = pd.DataFrame(master_list_up)
    up_marker_df = up_marker_df.transpose()
    up_marker_df.to_csv(f"{samplename}/up_markers.csv", index=False)
    down_marker_df = pd.DataFrame(master_list_down)
    down_marker_df = down_marker_df.transpose()
    down_marker_df.to_csv(f"{samplename}/down_markers.csv", index=False)

    ###Infer Cell Types

    pg.infer_cell_types(data, markers='AHBA_PFC_filtered.json', output_file=f"{samplename}/infer_cell_types_AHBA_markers")
    pg.infer_cell_types(data, markers='Hybrid_subclass_markers.json', output_file=f"{samplename}/infer_cell_types_Hybrid_markers")
    print("infer cell types using the Bakken et al markers done")

    if jdict["hashing"] == "True":
        HTOnames = set(data.obs["assignment"].values)
        print(f"Sample names = {HTOnames}")
        for sample in HTOnames:
            data_sample = MultimodalData(data[data.obs["assignment"]==sample,:].copy())
            dataset_anndata = f"{samplename}/{batchname}_{sample}_Processed.h5ad"
            pg.write_output(data_sample,dataset_anndata)
    # data_TPM_norm.obs[['leiden_labels']] = data.obs[['leiden_labels']]
    # if not os.path.exists(f'{samplename}'):
    #     os.mkdir(f'{samplename}')
    # pg.write_output(data_TPM_norm,f'{samplename}/Filtered_TPM_Multimodal_object.h5ad')

    summary_file.close()
    print("done")
