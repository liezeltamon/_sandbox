#!/bin/bash

#SBATCH --partition=short
#SBATCH --mem=20G
#SBATCH --ntasks=2
#SBATCH --time=1-00:00:00
#SBATCH --output=%j_%x.log.out
#SBATCH --error=%j_%x.log.err

load_mamba
mamba activate sandbox

python3 test.py
