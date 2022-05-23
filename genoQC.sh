#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=genoQC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=180000
#SBATCH --output=genoQC.%j.out
#SBATCH --error=genoQC.%j.err
#SBATCH --array=1-22

i=$SLURM_ARRAY_TASK_ID

cd /scratch/ahc87874/Check

ml PLINK/2.00-alpha2.3-x86_64-20210920-dev

###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###STEP 1. GENOTYPE QC PLINK-=-=-=-=-=-=-=-=
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
echo "-=-=-=-=-=-=-=-STEP 1-=-=-=-=-=-=-=-\n\n"

genoindir=("/scratch/ahc87874/bgen_v1.2_UKBsource")
mfiscoredir=("/scratch/ahc87874/Fall2021Practice/UKBpgen/mfi/info0.5")
outdir=("/scratch/ahc87874/Check/genoQC")
mkdir -p $outdir

plink2 \
--bgen $genoindir/ukb_imp_chr"$i"_v3.bgen ref-first \
--sample $genoindir/ukb_imp_v3.sample \
--extract $mfiscoredir/ukb_mfi_chr"$i"_v3_0.5.txt \
--mind 0.05 \
--geno 0.02 \
--hwe 1e-06 \
--maf 0.01 \
--autosome \
--maj-ref \
--keep /scratch/ahc87874/Check/GWAS_phenoQC_IDS_M1_Veg.txt \
--max-alleles 2 \
--freq
--export bgen-1.2 bits=8 \
--out "$outdir"/chr"$i"
