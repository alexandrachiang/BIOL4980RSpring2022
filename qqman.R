library(tidyverse)
library(qqman)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#TotalxCSRV
indir="/scratch/ahc87874/Check/GEM/TotalCholesterol"
#"/scratch/ahc87874/Check/GEM/TotalCholesterol" "/scratch/ahc87874/Check/GEM/LDLCholesterol" "/scratch/ahc87874/Check/GEM/Triglycerides"
#"TotalCholesterol" "LDLCholesterol" "HDLCholesterol" "Triglycerides"

for (i in 1:22) {
infile<-read.table(paste(indir, paste("TotalCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE)
#"TotalCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr" "TotalCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr"
#"LDLCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr" "LDLCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr"
#"TriglyceridesxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr"	
	
infile<-as_tibble(infile) 
  
#names(infile)<-c("CHROM", "POS", "ID", "REF", "ALT", "A1", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE")
#[1] "SNPID"
#[2] "RSID"
#[3] "CHR"
#[4] "POS"
#[5] "Non_Effect_Allele"
#[6] "Effect_Allele"
#[7] "N_Samples"
#[8] "AF"
#[9] "N_Consistent_Self_Reported_Vegetarian_across_all_24hr_1"
#[10] "AF_Consistent_Self_Reported_Vegetarian_across_all_24hr_1"
#[11] "N_Consistent_Self_Reported_Vegetarian_across_all_24hr_0"
#[12] "AF_Consistent_Self_Reported_Vegetarian_across_all_24hr_0"
#[13] "Beta_Marginal"
#[14] "robust_SE_Beta_Marginal"
#[15] "Beta_G.Consistent_Self_Reported_Vegetarian_across_all_24hr"
#[16] "robust_SE_Beta_G.Consistent_Self_Reported_Vegetarian_across_all_24hr"
#[17] "robust_P_Value_Marginal"
#[18] "robust_P_Value_Interaction"
#[19] "robust_P_Value_Joint"
  
#Subset data
infile1<-infile%>%select(CHR, POS, robust_P_Value_Interaction, RSID)

#Get qqman format
colnames(infile1)<-c("CHR", "BP", "P", "SNP")

#Add to input
if (i == 1) {
	infileall<-infile1
} else {
	infileall<-rbind(infileall, infile1)
}

}

outdir="/scratch/ahc87874/Check/SNPs"
#Make table of sig SNPs (P < 1e-5)
sigSNPs<-infileall%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/TotalxCSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
#Make table of top 10 SNPs
attach(infileall)
newdata <- infileall[order(P),]
newdata <- newdata[1:10,]
write.table(newdata, 
	paste(outdir, "/TotalxCSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
pvalue<-newdata[10,3]

#Make manhattan plot
outdir="/scratch/ahc87874/Check/manplots"
plotoutputfile<-paste(outdir, "/TotalCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr.png", sep="")

png(filename=plotoutputfile, type="cairo", width=600, height=300)
manhattan(infileall,ylim = c(0, 8.15), col = c("firebrick1", "black"), suggestiveline = T, genomewideline = T, main = "Manhattan Plot of CSRVxTotal GWIS", annotatePval = 1e-5)
#highlight = newdata
#firebrick1 deepskyblue1
dev.off()

#Make qq plot
outdir="/scratch/ahc87874/Check/qqplots"
plotoutputfile<-paste(outdir, "/TotalCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr.png", sep="")

png(filename=plotoutputfile, type="cairo")
qq(infileall$P, main = "Q-Q plot of CSRVxTotal GWIS p-values")
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#TotalxSSRV
indir="/scratch/ahc87874/Check/GEM/TotalCholesterol"
#"/scratch/ahc87874/Check/GEM/TotalCholesterol" "/scratch/ahc87874/Check/GEM/LDLCholesterol" "/scratch/ahc87874/Check/GEMsingle/SSRVvsTri/Triglycerides"
#"TotalCholesterol" "LDLCholesterol" "HDLCholesterol" "Triglycerides"

for (i in 1:22) {
infile<-read.table(paste(indir, paste("TotalCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE)
#"TotalCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr" "TotalCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr"
#"LDLCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr" "LDLCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr"
#"TriglyceridesxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr"	
	
infile<-as_tibble(infile) 
  
#Subset data
infile1<-infile%>%select(CHR, POS, robust_P_Value_Interaction, RSID)

#Get qqman format
colnames(infile1)<-c("CHR", "BP", "P", "SNP")

#Add to input
if (i == 1) {
	infileall<-infile1
} else {
	infileall<-rbind(infileall, infile1)
}

}

outdir="/scratch/ahc87874/Check/SNPs"
#Make table of sig SNPs (P < 1e-5)
sigSNPs<-infileall%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/TotalxSSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
#Make table of top 10 SNPs
attach(infileall)
newdata <- infileall[order(P),]
newdata <- newdata[1:10,]
write.table(newdata, 
	paste(outdir, "/TotalxSSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
pvalue<-newdata[10,3]

#Make manhattan plot
outdir="/scratch/ahc87874/Check/manplots"
plotoutputfile<-paste(outdir, "/TotalCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24.png", sep="")

png(filename=plotoutputfile, type="cairo", width=600, height=300)
manhattan(infileall,ylim = c(0, 8.15), col = c("deepskyblue1", "black"), suggestiveline = T, genomewideline = T, main = "Manhattan Plot of SSRVxTotal GWIS", annotatePval = 1e-5)
#highlight = newdata
#firebrick1 deepskyblue1
dev.off()

#Make qq plot
outdir="/scratch/ahc87874/Check/qqplots"
plotoutputfile<-paste(outdir, "/TotalCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24.png", sep="")

png(filename=plotoutputfile, type="cairo")
qq(infileall$P, main = "Q-Q plot of SSRVxTotal GWIS p-values")
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#LDLxCSRV
indir="/scratch/ahc87874/Check/GEM/LDLCholesterol"

for (i in 1:22) {
infile<-read.table(paste(indir, paste("LDLCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE)
	
infile<-as_tibble(infile) 

#Subset data
infile1<-infile%>%select(CHR, POS, robust_P_Value_Interaction, RSID)

#Get qqman format
colnames(infile1)<-c("CHR", "BP", "P", "SNP")

#Add to input
if (i == 1) {
	infileall<-infile1
} else {
	infileall<-rbind(infileall, infile1)
}

}

outdir="/scratch/ahc87874/Check/SNPs"
#Make table of sig SNPs (P < 1e-5)
sigSNPs<-infileall%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/LDLxCSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
#Make table of top 10 SNPs
attach(infileall)
newdata <- infileall[order(P),]
newdata <- newdata[1:10,]
write.table(newdata, 
	paste(outdir, "/LDLxCSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
pvalue<-newdata[10,3]

#Make manhattan plot
outdir="/scratch/ahc87874/Check/manplots"
plotoutputfile<-paste(outdir, "/LDLCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr.png", sep="")

png(filename=plotoutputfile, type="cairo", width=600, height=300)
manhattan(infileall,ylim = c(0, 8.15), col = c("firebrick1", "black"), suggestiveline = T, genomewideline = T, main = "Manhattan Plot of CSRVxLDL GWIS", annotatePval = 1e-5)
dev.off()

#Make qq plot
outdir="/scratch/ahc87874/Check/qqplots"
plotoutputfile<-paste(outdir, "/LDLCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr.png", sep="")

png(filename=plotoutputfile, type="cairo")
qq(infileall$P, main = "Q-Q plot of CSRVxLDL GWIS p-values")
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#LDLxSSRV
indir="/scratch/ahc87874/Check/GEM/LDLCholesterol"

for (i in 1:22) {
infile<-read.table(paste(indir, paste("LDLCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE)
	
infile<-as_tibble(infile) 
  
#Subset data
infile1<-infile%>%select(CHR, POS, robust_P_Value_Interaction, RSID)

#Get qqman format
colnames(infile1)<-c("CHR", "BP", "P", "SNP")

#Add to input
if (i == 1) {
	infileall<-infile1
} else {
	infileall<-rbind(infileall, infile1)
}

}

outdir="/scratch/ahc87874/Check/SNPs"
#Make table of sig SNPs (P < 1e-5)
sigSNPs<-infileall%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/LDLxSSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
#Make table of top 10 SNPs
attach(infileall)
newdata <- infileall[order(P),]
newdata <- newdata[1:10,]
write.table(newdata, 
	paste(outdir, "/LDLxSSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
pvalue<-newdata[10,3]

#Make manhattan plot
outdir="/scratch/ahc87874/Check/manplots"
plotoutputfile<-paste(outdir, "/LDLCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24.png", sep="")

png(filename=plotoutputfile, type="cairo", width=600, height=300)
manhattan(infileall,ylim = c(0, 8.15), col = c("deepskyblue1", "black"), suggestiveline = T, genomewideline = T, main = "Manhattan Plot of SSRVxLDL GWIS", annotatePval = 1e-5)
dev.off()

#Make qq plot
outdir="/scratch/ahc87874/Check/qqplots"
plotoutputfile<-paste(outdir, "/LDLCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24.png", sep="")

png(filename=plotoutputfile, type="cairo")
qq(infileall$P, main = "Q-Q plot of SSRVxLDL GWIS p-values")
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#HDLxCSRV
indir="/scratch/ahc87874/Check/GEM/HDLCholesterol"

for (i in 1:22) {
infile<-read.table(paste(indir, paste("HDLCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE)
	
infile<-as_tibble(infile) 

#Subset data
infile1<-infile%>%select(CHR, POS, robust_P_Value_Interaction, RSID)

#Get qqman format
colnames(infile1)<-c("CHR", "BP", "P", "SNP")

#Add to input
if (i == 1) {
	infileall<-infile1
} else {
	infileall<-rbind(infileall, infile1)
}

}

outdir="/scratch/ahc87874/Check/SNPs"
#Make table of sig SNPs (P < 1e-5)
sigSNPs<-infileall%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/HDLxCSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
#Make table of top 10 SNPs
attach(infileall)
newdata <- infileall[order(P),]
newdata <- newdata[1:10,]
write.table(newdata, 
	paste(outdir, "/HDLxCSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
pvalue<-newdata[10,3]

#Make manhattan plot
outdir="/scratch/ahc87874/Check/manplots"
plotoutputfile<-paste(outdir, "/HDLCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr.png", sep="")

png(filename=plotoutputfile, type="cairo", width=600, height=300)
manhattan(infileall,ylim = c(0, 8.15), col = c("firebrick1", "black"), suggestiveline = T, genomewideline = T, main = "Manhattan Plot of CSRVxHDL GWIS", annotatePval = 1e-5)
dev.off()

#Make qq plot
outdir="/scratch/ahc87874/Check/qqplots"
plotoutputfile<-paste(outdir, "/HDLCholesterolxConsistent_Self_Reported_Vegetarian_across_all_24hr.png", sep="")

png(filename=plotoutputfile, type="cairo")
qq(infileall$P, main = "Q-Q plot of CSRVxHDL GWIS p-values")
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#HDLxSSRV
indir="/scratch/ahc87874/Check/GEM/HDLCholesterol"

for (i in 1:22) {
infile<-read.table(paste(indir, paste("HDLCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE)
	
infile<-as_tibble(infile) 
  
#Subset data
infile1<-infile%>%select(CHR, POS, robust_P_Value_Interaction, RSID)

#Get qqman format
colnames(infile1)<-c("CHR", "BP", "P", "SNP")

#Add to input
if (i == 1) {
	infileall<-infile1
} else {
	infileall<-rbind(infileall, infile1)
}

}

outdir="/scratch/ahc87874/Check/SNPs"
#Make table of sig SNPs (P < 1e-5)
sigSNPs<-infileall%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/HDLxSSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
#Make table of top 10 SNPs
attach(infileall)
newdata <- infileall[order(P),]
newdata <- newdata[1:10,]
write.table(newdata, 
	paste(outdir, "/HDLxSSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
pvalue<-newdata[10,3]

#Make manhattan plot
outdir="/scratch/ahc87874/Check/manplots"
plotoutputfile<-paste(outdir, "/HDLCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24.png", sep="")

png(filename=plotoutputfile, type="cairo", width=600, height=300)
manhattan(infileall,ylim = c(0, 8.15), col = c("deepskyblue1", "black"), suggestiveline = T, genomewideline = T, main = "Manhattan Plot of SSRVxHDL GWIS", annotatePval = 1e-5)
dev.off()

#Make qq plot
outdir="/scratch/ahc87874/Check/qqplots"
plotoutputfile<-paste(outdir, "/HDLCholesterolxSelf_Reported_Vegetarian_plus_strict_initial_and24.png", sep="")

png(filename=plotoutputfile, type="cairo")
qq(infileall$P, main = "Q-Q plot of SSRVxHDL GWIS p-values")
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#TrixCSRV
indir="/scratch/ahc87874/Check/GEM/Triglycerides"

for (i in 1:22) {
infile<-read.table(paste(indir, paste("TriglyceridesxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE)
	
infile<-as_tibble(infile) 
  
#Subset data
infile1<-infile%>%select(CHR, POS, robust_P_Value_Interaction, RSID)

#Get qqman format
colnames(infile1)<-c("CHR", "BP", "P", "SNP")

#Add to input
if (i == 1) {
	infileall<-infile1
} else {
	infileall<-rbind(infileall, infile1)
}

}

outdir="/scratch/ahc87874/Check/SNPs"
#Make table of sig SNPs (P < 1e-5)
sigSNPs<-infileall%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/TAGxCSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
#Make table of top 10 SNPs
attach(infileall)
newdata <- infileall[order(P),]
newdata <- newdata[1:10,]
write.table(newdata, 
	paste(outdir, "/TAGxCSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
pvalue<-newdata[10,3]

#Make manhattan plot
outdir="/scratch/ahc87874/Check/manplots"
plotoutputfile<-paste(outdir, "/TriglyceridesxConsistent_Self_Reported_Vegetarian_across_all_24hr.png", sep="")

png(filename=plotoutputfile, type="cairo", width=600, height=300)
manhattan(infileall,ylim = c(0, 8.15), col = c("firebrick1", "black"), suggestiveline = T, genomewideline = T, main = "Manhattan Plot of CSRVxTAG GWIS", annotatePval = 1e-5)
dev.off()

#Make qq plot
outdir="/scratch/ahc87874/Check/qqplots"
plotoutputfile<-paste(outdir, "/TriglyceridesxConsistent_Self_Reported_Vegetarian_across_all_24hr-chr.png", sep="")

png(filename=plotoutputfile, type="cairo")
qq(infileall$P, main = "Q-Q plot of CSRVxTAG GWIS p-values")
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#TrixSSRV
indir="/scratch/ahc87874/Check/GEM/Triglycerides"

for (i in 1:22) {
infile<-read.table(paste(indir, paste("TriglyceridesxSelf_Reported_Vegetarian_plus_strict_initial_and24-chr", i, sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE)
	
infile<-as_tibble(infile) 
  
#Subset data
infile1<-infile%>%select(CHR, POS, robust_P_Value_Interaction, RSID)

#Get qqman format
colnames(infile1)<-c("CHR", "BP", "P", "SNP")

#Add to input
if (i == 1) {
	infileall<-infile1
} else {
	infileall<-rbind(infileall, infile1)
}

}

outdir="/scratch/ahc87874/Check/SNPs"
#Make table of sig SNPs (P < 1e-5)
sigSNPs<-infileall%>%filter(P<=1e-5)
write.table(sigSNPs, 
	paste(outdir, "/TAGxSSRVsigSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
#Make table of top 10 SNPs
attach(infileall)
newdata <- infileall[order(P),]
newdata <- newdata[1:10,]
write.table(newdata, 
	paste(outdir, "/TAGxSSRVtopSNPs.txt", sep=""),
	row.names=FALSE, quote=FALSE)
	
pvalue<-newdata[10,3]

#Make manhattan plot
outdir="/scratch/ahc87874/Check/manplots"
plotoutputfile<-paste(outdir, "/TriglyceridesxSelf_Reported_Vegetarian_plus_strict_initial_and24.png", sep="")

png(filename=plotoutputfile, type="cairo", width=600, height=300)
manhattan(infileall,ylim = c(0, 8.15), col = c("deepskyblue1", "black"), suggestiveline = T, genomewideline = T, main = "Manhattan Plot of SSRVxTAG GWIS", annotatePval = 1e-5)
dev.off()

#Make qq plot
outdir="/scratch/ahc87874/Check/qqplots"
plotoutputfile<-paste(outdir, "/TriglyceridesxSelf_Reported_Vegetarian_plus_strict_initial_and24.png", sep="")

png(filename=plotoutputfile, type="cairo")
qq(infileall$P, main = "Q-Q plot of SSRVxTAG GWIS p-values")
dev.off()
