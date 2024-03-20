#!/bin/bash
#SBATCH --job-name=cellranger-job
#SBATCH --out="slurm-%j.out"
#SBATCH --time=24:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=12
#SBATCH --mem-per-cpu=20G
#SBATCH --partition=pi_gerstein

python3 PEC_scRNA_pipeline_CellRanger.py  -N 12 -J /ysm-gpfs/pi/gerstein/ech43/PEC-pipeline-08-13/jsons-CMC/cmc-1-5.json
python3 PEC_scRNA_pipeline_CellRanger.py  -N 12 -J /ysm-gpfs/pi/gerstein/ech43/PEC-pipeline-08-13/jsons-CMC/cmc-6-8.json
