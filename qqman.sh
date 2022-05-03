#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=qqman
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=4:00:00
#SBATCH --mem=180000
#SBATCH --output=qqman.%j.out
#SBATCH --error=qqman.%j.err
#SBATCH --mail-user=ahc87874@uga.edu
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

module load R/4.1.3-foss-2020b
R CMD BATCH qqman.R
