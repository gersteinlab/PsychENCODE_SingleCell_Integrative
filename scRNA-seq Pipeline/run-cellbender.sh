#!/bin/bash
#SBATCH --job-name=cellranger-job
#SBATCH --out="slurm-%j.out"
#SBATCH --time=24:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=1
#SBATCH --mem-per-cpu=40G
#SBATCH --partition=pi_gerstein

python3 PEC_scRNA_pipeline_CellBender.py  -J CMC-Sets_1-5_CellBender_input.json
python3 PEC_scRNA_pipeline_CellBender.py  -J CMC-Sets_6-8_CellBender_input.json 
