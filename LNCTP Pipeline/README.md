# Cornerstone

Welcome, this repository contains the code and partial data for the LNCTP model used in our paper titled "Single-Cell Genomics Integration for Brain Studies", [under review] for publication in Science Journal.

## Overview
Single-cell genomics is a powerful tool, especially for heterogeneous tissues such as the brain. Yet, little is understood about how cell-level gene expression is influenced by individual genetic variants. Here, we integrate single-nucleus, multi-omics datasets from the PsychENCODE consortium to create a uniformly processed resource, comprising >2.8M nuclei from the prefrontal cortex across 388 individuals. The resource enables assessment of population-level variation in expression and chromatin for 28 cell types across gene families (e.g., for neurotransmitters and drug targets). It identifies >550K cell-type-specific regulatory elements and >1.4M single-cell eQTLs, which we use, in turn, to build cell-type regulatory networks and cell-to-cell communications networks. These networks provide a comprehensive picture of expression changes in aging and neuropsychiatric disorders. Finally, they enable the construction of an integrative model that accurately imputes cell-type gene expression. The model also prioritizes potential disease-risk genes and drug targets, with their associated cell types.

For a comprehensive understanding, the official publication provides detailed insights into our methodologies, results, and discussions. We recommend reading the paper alongside exploring this repository.

## Publication

- **Title**: Single-cell genomics & regulatory networks for 388 human brains
- **Authors**: [Name], [Co-Author`s Names], et al.
- **Journal**: Science Journal [?]
- **Link**: [Read the full paper here](PLACEHOLDER_FOR_PAPER_LINK)

## Key Highlights

- **Extensive Dataset**: More than 2.8 million nuclei from the prefrontal cortex across 388 individuals.
- **Cell-Type Specific Regulatory Elements**: Identification of over 550K unique regulatory elements.
- **Single-Cell eQTLs Exploration**: Unearthed more than 1.4 million eQTLs offering insights into cell-type regulatory networks and inter-cellular communication.
- **Aging and Neuropsychiatric Disorders**: Comprehensive network models highlighting expression changes associated with aging and various neuropsychiatric disorders.
- **Integrative Model Construction**: A model that accurately imputes cell-type gene expression, prioritizing potential disease-risk genes and drug targets, associating them with specific cell types.

## Repository Structure

- `*code_*/`: Contains all scripts, data and saved models LNCTP training, testing, and analysis for different types of phenotypes.
- `data/`: Processed datasets for the model.
- `models_saved/`: Optimal models selected for predicting disorders.
- `models/`: All trained models during the optimization process.

## Getting Started

1. Clone this repository: `git clone https://github.com/gersteinlab/Cornerstone.git`.
2. If modules are available, do `module load TensorFlow/2.7.1-foss-2020b-CUDA-11.3.1`
3. Do `python lnctp_*_models.py` for training and testing. Corresponding results will be saved (i.e. in `model_accuracy_summary.csv`).
4. Explore the scripts in the `*code_*/` directory to understand the LNCTP training and optimization pipeline.

## Dedicated Website

For more detailed information, interactive visualizations, and a dockerized version of the LNCTP for a quick start, please visit our [dedicated resource website](http://brainscope.psychencode.org/).

## Contributions

We welcome contributions! Please submit a pull request or raise an issue to discuss proposed changes.

## LNCTP Model Training

This section provides detailed instructions for training unaries and python models for the LNCTP model.

### I. Schizophrenia, Bipolar disorder, and ASD

#### Unary Training
Local predictors for each gene were trained using a Lasso loss function, to predict the z-score normalized expression from the eQTL SNPs associated with each gene.

```bash
module load MATLAB
matlab -nodisplay
```
- MATLAB Version: R2023a

```bash
run A_train_unariesZscore.m
```
- Unary models will be saved under 'unary_out/tag/unary_*.mat', where tag is asd, bip, etc.
- Change the corresponding variables (i.e. 'datDir', 'tag') to train for SCZ, ASD, or BIP.

#### GMRF Training
[To be completed by JW]

#### DNN Training
Optimized DNN via architectures and hyperparameters for best prediction power.

```bash
module load MATLAB
module load R/4.3.0-foss-2020b
run data2csv_wcgna_ros.m
```
- Take the [GMRF outputs] as inputs; and generating corresponding csv files for 'wgcna_modules_ros.R' to use

```bash
Rscript wgcna_modules_ros.R
```
- Take the [tarin-val-tests (from data2csv_wcgna_ros.m); geneIds_ens.txt; TableS5_Network_annotations.csv] as inputs, and generating training data for 'model_testing_10fold_updated_all_preds_weights_ros.py' to use

```bash
module purge
module load TensorFlow/2.7.1-foss-2020b-CUDA-11.3.1
module load UCX/1.9.0-GCCcore-10.2.0
python model_testing_10fold_updated_all_preds_weights_*.py
```
- Change * to asd, bip, or scz
- Modify the network parameters corresponding for different input features (i.e. 'params', 'n_params)

### II. ROSMAP dataset

All steps are the same as in [Schizophrenia, Bipolar disorder, and ASD](#schizophrenia-bipolar-disorder-and-asd), except using different source code:
- `A_train_unariesZscore_ros.m`
- [GMRF - JW]
- `wgcna_modules_ros.R`
- `model_testing_10fold_updated_all_preds_weights_ros.py`
