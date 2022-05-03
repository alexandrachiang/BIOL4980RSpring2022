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

outdir="/scratch/ahc87874/Check/FUMA"
#Save table of all chr for pheno x exposure
write.table(TAGxSSRV, 
	paste(outdir, "/TAGxSSRVall.txt", sep=""),
	row.names=FALSE, quote=FALSE)
  
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#Group by Exposure
#Group by Pheno
#Group All
