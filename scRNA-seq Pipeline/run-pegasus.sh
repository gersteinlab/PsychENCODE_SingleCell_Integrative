#!/bin/bash
#SBATCH --job-name=cellranger-job
#SBATCH --out="slurm-%j.out"
#SBATCH --time=48:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=12
#SBATCH --mem-per-cpu=40G
#SBATCH --partition=pi_gerstein

python3 PEC_scRNA_pipeline_Pegasus.py -J /ysm-gpfs/pi/gerstein/ech43/PEC-pipeline-08-13/NeMO_A46_samples_1.json

