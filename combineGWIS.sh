#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=combineGWIS
#SBATCH --ntasks=1
#SBATCH --time=4:00:00
#SBATCH --mem=18
#SBATCH --output=combineGWIS.%j.out
#SBATCH --error=combineGWIS.%j.err
#SBATCH --mail-user=ahc87874@uga.edu
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

module load R/4.1.3-foss-2020b
R CMD BATCH combineGWIS.R
