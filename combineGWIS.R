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

outdir="/scratch/ahc87874/Check/FUMA"
#Save table of all chr for pheno x exposure
write.table(TAGxSSRV, 
	paste(outdir, "/TAGxSSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
  
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

outdir="/scratch/ahc87874/Check/FUMA"

#Group by Exposure
CSRVall <- rbind(TotalxCSRV, LDLxCSRV, HDLxCSRV, TAGxCSRV)
SSRVall <- rbind(TotalxSSRV, LDLxSSRV, HDLxSSRV, TAGxSSRV)
write.table(CSRVall, 
	paste(outdir, "/CSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
write.table(SSRVall, 
	paste(outdir, "/SSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)

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

#Group All
All <- rbind(TotalxCSRV, TotalxSSRV, LDLxCSRV, LDLxSSRV, HDLxCSRV, HDLxSSRV, TAGxCSRV, TAGxSSRV)
write.table(All, 
	paste(outdir, "/All.txt", sep=""),
	row.names=FALSE, quote=FALSE)

