# Linear Network of Cell-Type Phenotypes (LNCTP) 

This folder contains the code and partial data for the LNCTP model.


## Subfolder Structure for each phenotype

- `*.py`:  Scripts
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
