#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=combineGWIS
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=20:00:00
#SBATCH --mem=1000
#SBATCH --output=combineGWIS.%j.out
#SBATCH --error=combineGWIS.%j.err
#SBATCH --mail-user=ahc87874@uga.edu
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

module load R/4.1.3-foss-2020b
R CMD BATCH combineGWIS.R
