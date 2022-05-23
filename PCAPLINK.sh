#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=PCAPLINK
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=160:00:00
#SBATCH --mem=200000
#SBATCH --output=PCAPLINK.%j.out
#SBATCH --error=PCAPLINK.%j.err
#SBATCH --array=1-22
#SBATCH --mail-user=ahc87874@uga.edu
#SBATCH --mail-type=ALL

i=$SLURM_ARRAY_TASK_ID

cd /scratch/ahc87874/Check

ml PLINK/2.00-alpha2.3-x86_64-20210920-dev

indir=("/scratch/ahc87874/Check/genoQC")
outdir=("/scratch/ahc87874/Check/PCAPLINK")

mkdir -p $outdir

plink2 \
--bgen "$indir"/chr"$i".bgen ref-first \
--sample "$indir"/chr"$i".sample \
--maf 0.01 \
--pca approx \
--out "$outdir"/chr"$i"_PCA
