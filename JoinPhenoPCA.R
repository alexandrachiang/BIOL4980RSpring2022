setwd("/scratch/ahc87874/Check")

library(plyr)
library(dplyr)
library(tidyverse)

pheno <- as_tibble(read.table("/scratch/ahc87874/Check/GWAS_pheno_M1_Veg.txt", header=TRUE))

outdir = "/scratch/ahc87874/Check/PhenoEigen"

for (i in 1:22) {
eigenPCA <- as_tibble(read.table(paste("/scratch/ahc87874/Check/PCAPLINK/chr", i, "_PCA.eigenvec", sep=""), header=FALSE))

colnames(eigenPCA) <- c("bFID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

phenoeigen <- left_join(pheno, eigenPCA,  by=c("IID"))

write.table(phenoeigen, 
	paste(outdir, "/GWAS_pheno_M1_Veg_eigen", i, ".txt", sep=""),
	row.names=FALSE, quote=FALSE)

write.csv(phenoeigen, paste(outdir, "/GWAS_pheno_M1_Veg_eigen", i, ".csv", sep=""),row.names=FALSE, quote=FALSE)
}

#x<-as_tibble(read.csv("GWAS_pheno_M1_Veg_eigen22.csv"))
