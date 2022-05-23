#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=GEM
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=144:00:00
#SBATCH --mem=30000
#SBATCH --output=GEM.%j.out
#SBATCH --error=GEM.%j.err
#SBATCH --mail-user=ahc87874@uga.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-22

i=$SLURM_ARRAY_TASK_ID

cd /scratch/ahc87874/Check

#ml PLINK/2.00-alpha2.3-x86_64-20210920-dev
ml GEM/1.4.1-foss-2019b


genoindir=("/scratch/ahc87874/Check/genoQC")
phenodir=("/scratch/ahc87874/Check/PhenoEigen")
outdir=("/scratch/ahc87874/Check/GEM")

phenotypes=("TotalCholesterol" "LDLCholesterol" "HDLCholesterol" "Triglycerides")

exposures=("Consistent_Self_Reported_Vegetarian_across_all_24hr" "Self_Reported_Vegetarian_plus_strict_initial_and24")


for j in ${phenotypes[@]} 
        do

for e in ${exposures[@]} 
        do

mkdir -p $outdir/$j

echo running "$j" and "$e"


GEM \
--bgen $genoindir/chr"$i".bgen \
--sample $genoindir/chr"$i".sample \
--pheno-file $phenodir/GWAS_pheno_M1_Veg_eigen"$i".csv \
--sampleid-name IID \
--pheno-name $j \
--covar-names Age Age2 Sex Geno_batch BMI statins \
center1 center2 center3 center4 center5 \
center6 center7 center8 center9 center10 \
center11 center12 center13 center14 center15 \
center16 center17 center18 center20 \
PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 \
--robust 1 \
--exposure-names "$e" \
--thread 16 \
--out $outdir/$j/"$j"x"$e"-chr"$i"

done
done
