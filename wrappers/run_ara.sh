#!/bin/bash
#SBATCH --qos long
#SBATCH --mem 4G
#SBATCH -c 2
#SBATCH -t 10-00:00:00
#SBATCH -o logs/%x_%j.log
#SBATCH -e logs/%x_%j.log
#SBATCH -J SQevi_ara

module load anaconda
source activate snakemake

/home/patienza/oscars/LR_annotation/sqanti_evidence -c configs/arabidopsis_config.yaml -j 20 -s -se "--default-resources slurm_account=gge" -R -e "--rerun-triggers input params code --max-jobs-per-second 5 --max-status-checks-per-second 5"
