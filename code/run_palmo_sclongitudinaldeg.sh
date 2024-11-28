#!/bin/bash

#SBATCH --partition=short
#SBATCH --mem=80G
#SBATCH --ntasks=1
#SBATCH --time=1-00:00:00
#SBATCH --output=%j_%x.log.out
#SBATCH --error=%j_%x.log.err

module load R-base/4.3.0
module load R-cbrg/current

Rscript --vanilla palmo_sclongitudinaldeg.R
