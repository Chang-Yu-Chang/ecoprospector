#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=community_selection
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu=6000

# Activate anaconda dependent packages
module load miniconda
source activate py37_dev
pip install -e test/community-simulator/

python simulate_algorithm.py 4