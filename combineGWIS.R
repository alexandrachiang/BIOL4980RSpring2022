library(plyr)
library(tidyverse)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#TotalxCSRV
indir="/scratch/ahc87874/Check/GEM/TotalCholesterol"

for (i in 1:22) {
  infile<-as_tibble(read.table(paste(indir, paste("TotalCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE))

  #Add to input
  if (i == 1) {
    TotalxCSRV<-infile
  } else {
    TotalxCSRV<-rbind(TotalxCSRV, infile)
  }
}

TotalxCSRV <- TotalxCSRV %>% mutate(pheno = "Total", exposure = "CSRV") %>% select(pheno, exposure, everything())

outdir="/scratch/ahc87874/Check/SNPsfull"
#Make table of top 10 SNPs
attach(TotalxCSRV)
newdata <- TotalxCSRV[order(robust_P_Value_Interaction),]

sigSNPs<-newdata%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/TotalxCSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

TotalxCSRVtopSNPs <- newdata[1:10,]
write.table(TotalxCSRVtopSNPs, 
	paste(outdir, "/TotalxCSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

outdir="/scratch/ahc87874/Check/FUMA"
#Save table of all chr for pheno x exposure
write.table(TotalxCSRV, 
	paste(outdir, "/TotalxCSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
  
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#TotalxSSRV
indir="/scratch/ahc87874/Check/GEM/TotalCholesterol"

for (i in 1:22) {
  infile<-as_tibble(read.table(paste(indir, paste("TotalCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE))

  #Add to input
  if (i == 1) {
    TotalxSSRV<-infile
  } else {
    TotalxSSRV<-rbind(TotalxSSRV, infile)
  }
}

TotalxSSRV <- TotalxSSRV %>% mutate(pheno = "Total", exposure = "SSRV") %>% select(pheno, exposure, everything())

outdir="/scratch/ahc87874/Check/SNPsfull"
#Make table of top 10 SNPs
attach(TotalxSSRV)
newdata <- TotalxSSRV[order(robust_P_Value_Interaction),]

sigSNPs<-newdata%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/TotalxSSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

TotalxSSRVtopSNPs <- newdata[1:10,]
write.table(TotalxSSRVtopSNPs, 
	paste(outdir, "/TotalxSSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

outdir="/scratch/ahc87874/Check/FUMA"
#Save table of all chr for pheno x exposure
write.table(TotalxSSRV, 
	paste(outdir, "/TotalxSSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
  
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#LDLxCSRV
indir="/scratch/ahc87874/Check/GEM/LDLCholesterol"

for (i in 1:22) {
  infile<-as_tibble(read.table(paste(indir, paste("LDLCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE))

  #Add to input
  if (i == 1) {
    LDLxCSRV<-infile
  } else {
    LDLxCSRV<-rbind(LDLxCSRV, infile)
  }
}

LDLxCSRV <- LDLxCSRV %>% mutate(pheno = "LDL", exposure = "CSRV") %>% select(pheno, exposure, everything())

outdir="/scratch/ahc87874/Check/SNPsfull"
#Make table of top 10 SNPs
attach(LDLxCSRV)
newdata <- LDLxCSRV[order(robust_P_Value_Interaction),]

sigSNPs<-newdata%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/LDLxCSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

LDLxCSRVtopSNPs <- newdata[1:10,]
write.table(LDLxCSRVtopSNPs, 
	paste(outdir, "/LDLxCSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

outdir="/scratch/ahc87874/Check/FUMA"
#Save table of all chr for pheno x exposure
write.table(LDLxCSRV, 
	paste(outdir, "/LDLxCSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
  
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#LDLxSSRV
indir="/scratch/ahc87874/Check/GEM/LDLCholesterol"

for (i in 1:22) {
  infile<-as_tibble(read.table(paste(indir, paste("LDLCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE))

  #Add to input
  if (i == 1) {
    LDLxSSRV<-infile
  } else {
    LDLxSSRV<-rbind(LDLxSSRV, infile)
  }
}

LDLxSSRV <- LDLxSSRV %>% mutate(pheno = "LDL", exposure = "SSRV") %>% select(pheno, exposure, everything())

outdir="/scratch/ahc87874/Check/SNPsfull"
#Make table of top 10 SNPs
attach(LDLxSSRV)
newdata <- LDLxSSRV[order(robust_P_Value_Interaction),]

sigSNPs<-newdata%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/LDLxSSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

LDLxSSRVtopSNPs <- newdata[1:10,]
write.table(LDLxSSRVtopSNPs, 
	paste(outdir, "/LDLxSSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

outdir="/scratch/ahc87874/Check/FUMA"
#Save table of all chr for pheno x exposure
write.table(LDLxSSRV, 
	paste(outdir, "/LDLxSSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
  
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#HDLxCSRV
indir="/scratch/ahc87874/Check/GEM/HDLCholesterol"

for (i in 1:22) {
  infile<-as_tibble(read.table(paste(indir, paste("HDLCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE))

  #Add to input
  if (i == 1) {
    HDLxCSRV<-infile
  } else {
    HDLxCSRV<-rbind(HDLxCSRV, infile)
  }
}

HDLxCSRV <- HDLxCSRV %>% mutate(pheno = "HDL", exposure = "CSRV") %>% select(pheno, exposure, everything())

outdir="/scratch/ahc87874/Check/SNPsfull"
#Make table of top 10 SNPs
attach(HDLxCSRV)
newdata <- HDLxCSRV[order(robust_P_Value_Interaction),]

sigSNPs<-newdata%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/HDLxCSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

HDLxCSRVtopSNPs <- newdata[1:10,]
write.table(HDLxCSRVtopSNPs, 
	paste(outdir, "/HDLxCSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

outdir="/scratch/ahc87874/Check/FUMA"
#Save table of all chr for pheno x exposure
write.table(HDLxCSRV, 
	paste(outdir, "/HDLxCSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
  
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#HDLxSSRV
indir="/scratch/ahc87874/Check/GEM/HDLCholesterol"

for (i in 1:22) {
  infile<-as_tibble(read.table(paste(indir, paste("HDLCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE))

  #Add to input
  if (i == 1) {
    HDLxSSRV<-infile
  } else {
    HDLxSSRV<-rbind(HDLxSSRV, infile)
  }
}

HDLxSSRV <- HDLxSSRV %>% mutate(pheno = "HDL", exposure = "SSRV") %>% select(pheno, exposure, everything())

outdir="/scratch/ahc87874/Check/SNPsfull"
#Make table of top 10 SNPs
attach(HDLxSSRV)
newdata <- HDLxSSRV[order(robust_P_Value_Interaction),]

sigSNPs<-newdata%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/HDLxSSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

HDLxSSRVtopSNPs <- newdata[1:10,]
write.table(HDLxSSRVtopSNPs, 
	paste(outdir, "/HDLxSSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

outdir="/scratch/ahc87874/Check/FUMA"
#Save table of all chr for pheno x exposure
write.table(HDLxSSRV, 
	paste(outdir, "/HDLxSSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
  
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#TAGxCSRV
indir="/scratch/ahc87874/Check/GEM/Triglycerides"

for (i in 1:22) {
  infile<-as_tibble(read.table(paste(indir, paste("TriglyceridesxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE))

  #Add to input
  if (i == 1) {
    TAGxCSRV<-infile
  } else {
    TAGxCSRV<-rbind(TAGxCSRV, infile)
  }
}

TAGxCSRV <- TAGxCSRV %>% mutate(pheno = "TAG", exposure = "CSRV") %>% select(pheno, exposure, everything())

outdir="/scratch/ahc87874/Check/SNPsfull"
#Make table of top 10 SNPs
attach(TAGxCSRV)
newdata <- TAGxCSRV[order(robust_P_Value_Interaction),]

sigSNPs<-newdata%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/TAGxCSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

TAGxCSRVtopSNPs <- newdata[1:10,]
write.table(TAGxCSRVtopSNPs, 
	paste(outdir, "/TAGxCSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

outdir="/scratch/ahc87874/Check/FUMA"
#Save table of all chr for pheno x exposure
write.table(TAGxCSRV, 
	paste(outdir, "/TAGxCSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
  
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#TAGxSSRV
indir="/scratch/ahc87874/Check/GEM/Triglycerides"

for (i in 1:22) {
  infile<-as_tibble(read.table(paste(indir, paste("TriglyceridesxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE))

  #Add to input
  if (i == 1) {
    TAGxSSRV<-infile
  } else {
    TAGxSSRV<-rbind(TAGxSSRV, infile)
  }
}

TAGxSSRV <- TAGxSSRV %>% mutate(pheno = "TAG", exposure = "SSRV") %>% select(pheno, exposure, everything())

outdir="/scratch/ahc87874/Check/SNPsfull"
#Make table of top 10 SNPs
attach(TAGxSSRV)
newdata <- TAGxSSRV[order(robust_P_Value_Interaction),]

sigSNPs<-newdata%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/TAGxSSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

TAGxSSRVtopSNPs <- newdata[1:10,]
write.table(TAGxSSRVtopSNPs, 
	paste(outdir, "/TAGxSSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)

outdir="/scratch/ahc87874/Check/FUMA"
#Save table of all chr for pheno x exposure
write.table(TAGxSSRV, 
	paste(outdir, "/TAGxSSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
  
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

setwd("/scratch/ahc87874/Check/FUMA")
TotalxCSRV <- as_tibble(read.table("TotalxCSRVall.txt", header = TRUE))
TotalxSSRV <- as_tibble(read.table("TotalxSSRVall.txt", header = TRUE))
LDLxCSRV <- as_tibble(read.table("LDLxCSRVall.txt", header = TRUE))
LDLxSSRV <- as_tibble(read.table("LDLxSSRVall.txt", header = TRUE))
HDLxCSRV <- as_tibble(read.table("HDLxCSRVall.txt", header = TRUE))
HDLxSSRV <- as_tibble(read.table("HDLxSSRVall.txt", header = TRUE))
TAGxCSRV <- as_tibble(read.table("TAGxCSRVall.txt", header = TRUE))
TAGxSSRV <- as_tibble(read.table("TAGxSSRVall.txt", header = TRUE))

x <- c("pheno", "exposure", "SNPID", "RSID", "CHR", "POS", "Non_Effect_Allele", "Effect_Allele", "N_Samples", "AF", 
       "N_1", "AF_1", "N_0", "AF_0", 
       "Beta_Marginal", "robust_SE_Beta_Marginal", 
       "Beta_G.", "robust_SE_Beta_G.", 
       "robust_P_Value_Marginal", "robust_P_Value_Interaction", "robust_P_Value_Joint"
)
colnames(TotalxCSRVtopSNPs) <- x
colnames(TotalxSSRVtopSNPs) <- x
colnames(LDLxCSRVtopSNPs) <- x
colnames(LDLxSSRVtopSNPs) <- x
colnames(HDLxCSRVtopSNPs) <- x
colnames(HDLxSSRVtopSNPs) <- x
colnames(TAGxCSRVtopSNPs) <- x
colnames(TAGxSSRVtopSNPs) <- x

topSNPs <- rbind(TotalxCSRVtopSNPs, TotalxSSRVtopSNPs, 
		 LDLxCSRVtopSNPs, LDLxSSRVtopSNPs,
		 HDLxCSRVtopSNPs, HDLxSSRVtopSNPs,
		 TAGxCSRVtopSNPs, TAGxSSRVtopSNPs)
outdir="/scratch/ahc87874/Check/SNPsfull"
attach(topSNPs)
topSNPs <- topSNPs[order(robust_P_Value_Interaction),]
write.table(topSNPs, 
	paste(outdir, "/topSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
write.csv(topSNPs, paste(outdir, "/topSNPs.csv", sep=""), 
	  row.names=FALSE, quote=FALSE)

setwd("/scratch/ahc87874/Check/FUMA")
outdir="/scratch/ahc87874/Check/FUMA"
colnames(TotalxCSRV) <- x
colnames(TotalxSSRV) <- x
colnames(LDLxCSRV) <- x
colnames(LDLxSSRV) <- x
colnames(HDLxCSRV) <- x
colnames(HDLxSSRV) <- x
colnames(TAGxCSRV) <- x
colnames(TAGxSSRV) <- x

#Group by Exposure
CSRVall <- rbind(TotalxCSRV, LDLxCSRV, HDLxCSRV, TAGxCSRV)
SSRVall <- rbind(TotalxSSRV, LDLxSSRV, HDLxSSRV, TAGxSSRV)
write.table(CSRVall, 
	paste(outdir, "/CSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
write.table(SSRVall, 
	paste(outdir, "/SSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
gzip("CSRVall.txt")
gzip("SSRVall.txt")

#Group by Pheno
Totalall <- rbind(TotalxCSRV, TotalxSSRV)
LDLall <- rbind(LDLxCSRV, LDLxSSRV)
HDLall <- rbind(HDLxCSRV, HDLxSSRV)
TAGall <- rbind(TAGxCSRV, TAGxSSRV)
write.table(Totalall, 
	paste(outdir, "/Totalall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
write.table(LDLall, 
	paste(outdir, "/LDLall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
write.table(HDLall, 
	paste(outdir, "/HDLall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
write.table(TAGall, 
	paste(outdir, "/TAGall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
gzip("Totalall.txt")
gzip("LDLall.txt")
gzip("HDLall.txt")
gzip("TAGall.txt")

#Group All
All <- rbind(TotalxCSRV, TotalxSSRV, LDLxCSRV, LDLxSSRV, HDLxCSRV, HDLxSSRV, TAGxCSRV, TAGxSSRV)
write.table(All, 
	paste(outdir, "/All.txt", sep=""),
	row.names=FALSE, quote=FALSE)
gzip("All.txt")
