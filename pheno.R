# combine load_UKBphenotables.R, UKB_participantQC.R, VegetarianQC1.R, VegetarianQC2.R, and pheno.R w/o PCs

qlogin3
Rload
R
setwd("/scratch/ahc87874/Check")

library(plyr)
library(dplyr)
library(tidyverse)

#=-=-=-=-=-=-=-=-=-=-+-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#UKB_participantQC.R
### This is a script for generating lists of UK Biobank participants which meet QC/filtering criteria.
source("/scratch/ahc87874/Fall2021Practice/Project/load_UKBphenotables.R")

pan<-read_tsv("/scratch/ahc87874/Fall2021Practice/Project/all_pops_non_eur_pruned_within_pop_pc_covs.tsv")
pan<-as_tibble(pan)
pan$s<-as.integer(pan$s)
table(pan$pop, useNA = "always")

bridge<-read.table("/scratch/ahc87874/Fall2021Practice/Project/ukb48818bridge31063.txt")
bridge<-as_tibble(bridge)
colnames(bridge)<-c("IID", "panID")

pan2<-pan%>%select(s, pop)%>%
    left_join(bridge, by=c("s"="panID"))

pan2

##===============================PART A=================================================
#Generate a list of participants who pass the following QC criteria:
#1. Genetic ethnicity = Caucasian VIA PAN UKBB
#2. Not an outlier for heterogeneity and missing genotype rate (poor quality genotype)
#3. No Sex chromosome aneuploidy
#4. Self-reported sex matches genetic sex
#5. Do not have high degree of genetic kinship (Ten or more third-degree relatives identified)
#6. Does not appear in "maximum_set_of_unrelated_individuals.MF.pl"

bd_QC<- bd %>% select(f.eid, f.31.0.0, f.22001.0.0, f.21000.0.0,
                      f.22027.0.0, f.22019.0.0,
                      f.22021.0.0)

colnames(bd_QC)<-c("IID", "Sex", "Genetic_Sex", "Race",
                   "Outliers_for_het_or_missing", "SexchrAneuploidy",
                   "Genetic_kinship")

#1. Genetic ethnicity = Caucasian VIA PAN UKBB
#Join UKB cols with with Pan UKBB
bd_QC<-as_tibble(bd_QC) #502,527
bd_QC<-bd_QC%>%inner_join(pan2, by="IID") #448193

#Filter by Genetic ethnicity = Caucasian VIA PAN UKBB
bd_QC<-bd_QC[bd_QC$pop=="EUR",] #426881

#2. Not an outlier for heterogeneity and missing genotype rate (poor quality genotype)
bd_QC<-bd_QC%>%
    filter(is.na(Outliers_for_het_or_missing) | Outliers_for_het_or_missing !="Yes") #426433

#3. No Sex chromosome aneuploidy
bd_QC<-bd_QC%>%
    filter(is.na(SexchrAneuploidy) | SexchrAneuploidy != "Yes") #425854

#4. Self-reported sex matches genetic sex
#If Sex does not equal genetic sex, exclude participant
bd_QC<-bd_QC[bd_QC$Sex == bd_QC$Genetic_Sex,] #425683

#5. Do not have high degree of genetic kinship (Ten or more third-degree relatives identified)
bd_QC<- bd_QC%>%
    filter(is.na(Genetic_kinship) |
               Genetic_kinship != "Ten or more third-degree relatives identified") #425510
               
#6. Does not appear in "maximum_set_of_unrelated_individuals.MF.pl"
#Filter related file by those in QC
relatives<-read.table("ukb48818_rel_s488282.dat", header=T)

#From maximum_set_of_unrelated_individuals.MF.pl output:
max_unrelated<-read.table("ukb48818_rel_s488282_output.dat")
max_unrelated<-as.integer(unlist(max_unrelated))
bd_QC<-bd_QC%>%filter(!IID %in% max_unrelated) #356980

QCkeepparticipants<-bd_QC%>%select(IID)

write.table(QCkeepparticipants, file = "/scratch/ahc87874/Check/bd_QC-keep.txt",
            row.names = FALSE, quote = FALSE)
            
#Start with 502527 participants
#End with 356980 participants, removed 145547

#=-=-=-=-=-=-=-=-=-=-+-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#VegetarianQC1.R

#The key of this script is solving the problem that if someone has an
#NA in field 20086 (24hr-recall: special diet followed), it could be 
#because they did not answer this question, or because they didn't 
#take the 24h recall. This parses out the answer to that question.

#Load UK Biobank datasets-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#UKB data for people who took the 24hour recall survey
bd24<-as_tibble(read.table("/scratch/ahc87874/Fall2021Practice/Project/UKB_24hourRecall-ParticipantInstancesTaken.txt",
                 header=TRUE, stringsAsFactors = FALSE))
#sum(bd24$ever_took_24HR)#[1] 211018

#24hr special diet columns
SpecialDiet<-bd%>%select(f.eid, f.20086.0.0, f.20086.1.0,
                            f.20086.2.0, f.20086.3.0, f.20086.4.0)

#Column names
colnames(SpecialDiet)<- c("IID", "Veg0","Veg1", "Veg2", "Veg3", "Veg4")
table(SpecialDiet$Veg0) 


#***Use sapply here to perform repeated functions over columns***
#Code columns so vegetarian or vegan = 1, other = 0
SpecialDiet[,2:6]<-sapply(SpecialDiet[,2:6], as.character)
SpecialDiet[,2:6]<-sapply(SpecialDiet[,2:6], 
                         mapvalues, c(NA, "Low calorie", "Gluten-free", 
                                      "Lactose-free", "Other",
                                      "Vegetarian", "Vegan"), 
                         c(0,0,0,0,0,1,1))

###This section creates a table that has both answers to special diet
###columns (Veg0-4) and columns indicating whether people took the
###24hr survey (took0-4)

TooksurveyAndVeg <- left_join(SpecialDiet, bd24,  by=c("IID"))
TooksurveyAndVeg<-as_tibble(TooksurveyAndVeg)
TooksurveyAndVeg<-TooksurveyAndVeg[TooksurveyAndVeg$ever_took_24HR==1,]
TooksurveyAndVeg<-TooksurveyAndVeg%>%select(-ever_took_24HR)

colnames(TooksurveyAndVeg)[7:11]<-c("took0", "took1", "took2", "took3", "took4")
#nrow(TooksurveyAndVeg) #[1] 211018

#TooksurveyAndVeg[!complete.cases(TooksurveyAndVeg),] #0

#initialize answer columns
TooksurveyAndVeg$Answer0<-"d" #column 12
TooksurveyAndVeg$Answer1<-"d" #column 13
TooksurveyAndVeg$Answer2<-"d" #column 14
TooksurveyAndVeg$Answer3<-"d" #column 15
TooksurveyAndVeg$Answer4<-"d" #column 16

#There are 4 conditions here:
#a) took survey, was vegetarian
#b) took survey, wasn't vegetarian
#c) didn't take survey in this instance
#d) code failed to update answer

#This loop takes between 10 and 20 minutes

for(instance in 7:11){ #columns 7 through 11 are the "took" columns
    for(row in 1:nrow(TooksurveyAndVeg)){ #loop through every row in the table
        if(TooksurveyAndVeg[row, instance] == 1){ #if they took the survey in this instance
            if(TooksurveyAndVeg[row, instance-5]==1){ #if they had 1 in the same instance of veg
                TooksurveyAndVeg[row, instance+5]<-"a" 
            }
            else TooksurveyAndVeg[row, instance+5]<-"b"
        }
        else(TooksurveyAndVeg[row, instance+5]<-"c")
    }
}


TooksurveyAndVeg
#adding column for diet 
TooksurveyAndVeg$Diet0<-"0"
TooksurveyAndVeg$Diet1<-"0"
TooksurveyAndVeg$Diet2<-"0"
TooksurveyAndVeg$Diet3<-"0"
TooksurveyAndVeg$Diet4<-"0"

#********************************************************************
#set b's to -100 because they once answered that they were not 
#vegetarian in the survey; we are only looking to designate people 
#with 100% of answers as yes vegetarian

TooksurveyAndVeg


TooksurveyAndVeg$Diet0[TooksurveyAndVeg$Answer0=="a"]<-1
TooksurveyAndVeg$Diet0[TooksurveyAndVeg$Answer0=="b"]<-(-100)
TooksurveyAndVeg$Diet0[TooksurveyAndVeg$Answer0=="c"]<-0
TooksurveyAndVeg$Diet1[TooksurveyAndVeg$Answer1=="a"]<-1
TooksurveyAndVeg$Diet1[TooksurveyAndVeg$Answer1=="b"]<-(-100)
TooksurveyAndVeg$Diet1[TooksurveyAndVeg$Answer1=="c"]<-0
TooksurveyAndVeg$Diet2[TooksurveyAndVeg$Answer2=="a"]<-1
TooksurveyAndVeg$Diet2[TooksurveyAndVeg$Answer2=="b"]<-(-100)
TooksurveyAndVeg$Diet2[TooksurveyAndVeg$Answer2=="c"]<-0
TooksurveyAndVeg$Diet3[TooksurveyAndVeg$Answer3=="a"]<-1
TooksurveyAndVeg$Diet3[TooksurveyAndVeg$Answer3=="b"]<-(-100)
TooksurveyAndVeg$Diet3[TooksurveyAndVeg$Answer3=="c"]<-0
TooksurveyAndVeg$Diet4[TooksurveyAndVeg$Answer4=="a"]<-1
TooksurveyAndVeg$Diet4[TooksurveyAndVeg$Answer4=="b"]<-(-100)
TooksurveyAndVeg$Diet4[TooksurveyAndVeg$Answer4=="c"]<-0



#********************************************************************

Diet<-TooksurveyAndVeg%>%select(IID, Diet0, Diet1, Diet2, Diet3, Diet4) #subsetting information

Diet

Diet<-Diet%>%mutate_if(is.character, as.numeric) #changing Diet to numeric 
#this above line was necessary because these rows were initialized 
#with "o"'s (a character) instead of "0"'s (a number)
Diet

Diet$row_sum = rowSums(Diet[,c(2,3,4,5,6)]) #sum of the rows

Diet
Diet%>%filter(row_sum== 0) #nobody. good.

Diet$Veg<-"a" #creating column for final answer

Diet

Diet$Veg[Diet$row_sum < 0]<-"Non-vegetarian"
Diet$Veg[Diet$row_sum > 0]<-"Vegetarian"

Diet
Diet%>%filter(Veg== "Vegetarian") #5,738
Diet%>%filter(Veg== "Non-vegetarian") #205,280
5738/205280#[1] 0.02795207


Diet$Vegb<-0 #creating column for final answer

Diet$Vegb[Diet$Veg == "Non-vegetarian" ]<-0
Diet$Vegb[Diet$Veg ==  "Vegetarian" ]<-1

sum(Diet$Vegb) #[1] 3321

veg<-Diet%>%select("IID", "Vegb")

colnames(veg)<-c("IID", "Consistent_Self-Reported_Vegetarian_across_all_24hr")

write.table(veg, file = "/scratch/ahc87874/Check/vegQC1_04032021.txt", 
            sep = "\t", col.names = TRUE, quote = FALSE,
            row.names = FALSE)

#=-=-=-=-=-=-=-=-=-=-+-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#VegetarianQC2.R

vegqc1 <- veg

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
#First check the 24hr columns for "Meat consumer" and "fish consumer"
#Across all instances where the 24hr was taken.
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=

Meatfish24<-bd%>%select(f.eid, 
                        #Category 100106: Meat/fish yesterday. 24HR.
                        f.103000.0.0, f.103000.1.0, f.103000.2.0, 
                        f.103000.3.0, f.103000.4.0, 
                        
                        f.103140.0.0, f.103140.1.0, f.103140.2.0,
                        f.103140.3.0, f.103140.4.0)


colnames(Meatfish24)<-c("IID", "Meat0", "Meat1","Meat2","Meat3","Meat4",
                        "Fish0","Fish1","Fish2","Fish3","Fish4")

TooksurveyAndMeat<- left_join(Meatfish24, bd24,  by=c("IID"))
TooksurveyAndMeat<-as_tibble(TooksurveyAndMeat)
TooksurveyAndMeat<-TooksurveyAndMeat[TooksurveyAndMeat$ever_took_24HR==1,]
TooksurveyAndMeat<-TooksurveyAndMeat%>%select(-ever_took_24HR)

TooksurveyAndMeat[,2:11]<-sapply(TooksurveyAndMeat[,2:11], 
                                 mapvalues, c("No", "Yes"), c(0, 1))

TooksurveyAndMeat$Meatfish0<-0 #column 12
TooksurveyAndMeat$Meatfish1<-0 #column 13
TooksurveyAndMeat$Meatfish2<-0 #column 14
TooksurveyAndMeat$Meatfish3<-0 #column 15
TooksurveyAndMeat$Meatfish4<-0 #column 16

#There are 3 conditions here:
#1 = took survey, didn't say they eat meat or fish
#-100 = took survey, did say they eat meat or fish
#0 didn't take survey in this instance


#This loop takes between 10 and 20 minutes start 5:14
for(instance in 12:16){ #columns 12 through 16 are the "took" columns
    for(row in 1:nrow(TooksurveyAndMeat)){ #loop through every row in the table
        if(TooksurveyAndMeat[row, instance] == 1){ #if they took the survey in this instance
            if((TooksurveyAndMeat[row, instance-10]==0) & (TooksurveyAndMeat[row, instance-5]==0))
                { #if they had 0 in the same instance of Meat AND in Fish
                TooksurveyAndMeat[row, instance+5]<-1
                }
            else TooksurveyAndMeat[row, instance+5]<-(-100)
        }
        else(TooksurveyAndMeat[row, instance+5]<-0)
    }
}

Diet<-TooksurveyAndMeat%>%select(IID, Meatfish0, Meatfish1, Meatfish2, 
                                 Meatfish3, Meatfish4)

Diet$row_sum = rowSums(Diet[,c(2,3,4,5,6)]) #sum of the rows
Diet$Meatfish[Diet$row_sum > 0]<-0 #didn't eat meat or fish
Diet$Meatfish[Diet$row_sum < 0]<-1 #did eat meat or fish

#negative number = ate meat or fish in these columns.
MF<-Diet%>%select(IID, Meatfish)
vegqc2<-left_join(vegqc1, MF, by="IID")
vegqc2[vegqc2$'Consistent_Self-Reported_Vegetarian_across_all_24hr'==1 & vegqc2$Meatfish==1,] #563

#Get number of people identified by this QC check:: 563


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
#Second check with data from initial assessment for those who ate meat
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=

#selecting data to tibble
new<-bd%>%select(f.eid, 
                    #Category 100052: Diet. Initial assessment
                    f.1329.0.0, f.1339.0.0, 
                    f.1349.0.0, f.1359.0.0, 
                    f.1369.0.0, f.1379.0.0, 
                    f.1389.0.0
                    )

#column names for tibble
colnames(new)<- c("IID", 
                  #Category 100052: Diet. Initial assessment
                  "Oily_fish_intake", "Non_oily_fish_intake",
                  "Processed_meat_intake","Poultry_intake",
                  "Beef_intake","Lamb.mutton_intake",
                  "Pork_intake"
)


#Category 100052: keep only those who answered "No"
#Category 100106: keep only those who answered "No"

#Combine diet data with result from VegetarianQC1.R
(new<-left_join(vegqc1, new, by="IID"))

new[,3:9]<-sapply(new[,3:9], as.character)
new
nrow(new) #[1] 211018
sapply(new, function(x) sum(is.na(x)))
#85 people did not answer to the Initial assessment diet questions.
#Remove them here:
new<-new[!is.na(new$Pork_intake),]
nrow(new) #[1] 210933
new$strict_initial<-0

#QC for diet data
new$strict_initial[(new$Oily_fish_intake =="Never") &
    (new$Non_oily_fish_intake =="Never") &
    (new$Processed_meat_intake =="Never") &
    (new$Poultry_intake =="Never") &
    (new$Beef_intake=="Never") &
    (new$Lamb.mutton_intake =="Never") &
    (new$Pork_intake =="Never")] <-1

sum(new$strict_initial) #[1] 4749
strictinitaltable<-new%>%select(IID, strict_initial)
colnames(vegqc2)[colnames(vegqc2)=="Meatfish"]<-"Meatfish24"
vegqc2<-inner_join(vegqc2, strictinitaltable, by="IID")
vegqc2$strict_initial_and24<-0
vegqc2$strict_initial_and24[(vegqc2$'Consistent_Self-Reported_Vegetarian_across_all_24hr'==1)&
                             (vegqc2$strict_initial==1)&
                             (vegqc2$Meatfish24==0)]<-1

sum(vegqc2$strict_initial_and24) #3784

vegqc3<-vegqc2%>%select(IID, 'Consistent_Self-Reported_Vegetarian_across_all_24hr',
                        strict_initial_and24)


colnames(vegqc3)<-c("IID", "Consistent_Self_Reported_Vegetarian_across_all_24hr", 
                    "Self_Reported_Vegetarian_plus_strict_initial_and24")
table(vegqc3$Consistent_Self_Reported_Vegetarian_across_all_24hr, useNA = "always")
#      0      1   <NA> 
#   205200   5733      0

table(vegqc3$Self_Reported_Vegetarian_plus_strict_initial_and24, useNA = "always")
#     0      1   <NA> 
#   207149   3784      0 


write.table(vegqc3, file = "vegQC2_04032021.txt", 
            sep = "\t", col.names = TRUE, quote = FALSE,
            row.names = FALSE)

#=-=-=-=-=-=-=-=-=-=-+-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#pheno.R

suppressMessages(library(ggpubr))
suppressMessages(library(RNOmni))

source("/scratch/ahc87874/Fall2021Practice/Project/manyColsToDummy.R")

###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#Load data=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

withdrawn<-read.csv("/scratch/ahc87874/Fall2021Practice/Project/w48818_20210809.csv", header = FALSE)

QCids<-read.table("bd_QC-keep.txt",header=TRUE)

#Phenotypes  ------------------------------------------------------------------------------
#Covariates 
#Model 1: sex, age, age squared, genotyping array, and assessment center indicators (sites of recruitment); 
#Cholesterol;  


pheno<-bd%>%select(f.eid, f.21003.0.0, f.31.0.0, 
                   f.189.0.0,
                   f.54.0.0, f.22000.0.0
                    )

colnames(pheno)<-c("IID", "Age", "Sex",  
                   "Townsend",
                   "Assessment_center", "Geno_batch"
                    )

pheno2<-bd_join4%>%select(f.eid,
                          f.21001.0.0, f.30690.0.0, f.30780.0.0,
			  f.30760.0.0, f.30870.0.0
                          )

colnames(pheno2)<-c("IID",
                    "BMI", "TotalCholesterol", "LDLCholesterol",
		    "HDLCholesterol", "Triglycerides"
		   )

new<-left_join(pheno, pheno2, by="IID")
new<-as_tibble(new)

#Remove withdrawn participants------------------------------------
new<-new[!(new$IID %in% withdrawn$V1), ]

#QC participants via output of UKB_participantQC.R----------------

new<-new[(new$IID %in% QCids$IID),]

#Age squared----------------------------
new$Age2<-new$Age^2

#Make dummy 0/1 cols for each assessment center----------------------
#table(pheno$Assessment_center)
centers<-unique(new$Assessment_center)
centercols<-paste("center", 1:22, sep="")
new[centercols]<-0

for (i in 1:length(centers)){
    new[new$Assessment_center==centers[i],][centercols[i]]<-1
}

new<-new%>%select(-Assessment_center)
new

#Genotype batch
new$Geno_batch1<-0
new$Geno_batch1[new$Geno_batch>0]<-1
#sum(pheno$Geno_batch1) #[1] 438313
new$Geno_batch<-new$Geno_batch1
new<-new%>%select(-Geno_batch1)
#table(new$Geno_batch) #it worked


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#Switch sex values to numeric
new$Sex<-mapvalues(as.character(new$Sex), 
                     c("Male", "Female"), c(0,1))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#Statins
statincols<-c(sprintf("f.20003.0.%s", 0:47))
statincodes<-c(1141146234,1141192414,1140910632,1140888594,1140864592,
	1141146138,1140861970,1140888648,1141192410,
	1141188146,1140861958,1140881748,1141200040)

manyColsToDummy(statincodes, bd_join4[,statincols], "statinoutput")
statinoutput$statins<-rowSums(statinoutput) 
statinoutput$statins[statinoutput$statins>1]<-1

statinoutput$IID<-bd_join4$f.eid

statinoutput<-statinoutput%>%select(IID, statins)

new<-left_join(new, statinoutput, by="IID")

participants1<-new%>%select(IID)
participants1$FID<-participants1$IID
participants1<-participants1%>%select(FID, IID)

###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###WRITE OUTPUT=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
###=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

outdir="/scratch/ahc87874/Fall2021Practice/Project/pheno"

#Model 1
write.table(participants1, 
	paste(outdir, "/GWAS_phenoQC_IDS_M1.txt",sep=""), 
	row.names=FALSE, quote=FALSE)

write.table(new, 
	paste(outdir, "/GWAS_pheno_M1.txt", sep=""),
	row.names=FALSE, quote=FALSE)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#Add vegetarian columns
vegqc3<-read.delim("vegQC2_04032021.txt", header=TRUE, sep="\t")

#Inner join
new2<-as_tibble(inner_join(new, vegqc3, by="IID"))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#PCs
#pan<-read_tsv("/scratch/ahc87874/Fall2021Practice/Project/all_pops_non_eur_pruned_within_pop_pc_covs.tsv")
#pan<-as_tibble(pan)
#pan$s<-as.integer(pan$s)
#table(pan$pop, useNA = "always")

#bridge<-read.table("/scratch/ahc87874/Fall2021Practice/Project/ukb48818bridge31063.txt")
#bridge<-as_tibble(bridge)
#colnames(bridge)<-c("IID", "panID")

#pan2<-pan%>%select(s, pop, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)%>%
#    left_join(bridge, by=c("s"="panID"))

#Inner join
#new3<-as_tibble(inner_join(new2, pan2, by="IID"))

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#Complete cases
#new4<-new3[complete.cases(new3),]
#removes 26,173
#163,712

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
FID<-new2$IID
new2<-bind_cols(FID,new2)
colnames(new2)[1] <- "FID"

participants <- new2 %>% select("FID", "IID")
participants$FID <- participants$IID
participants <- participants%>%select(FID, IID)

#Output
outdir="/scratch/ahc87874/Check"

write.table(participants, 
	paste(outdir, "/GWAS_phenoQC_IDS_M1_Veg.txt", sep=""), 
	row.names=FALSE, quote=FALSE)

write.table(new2, 
	paste(outdir, "/GWAS_pheno_M1_Veg.txt", sep=""),
	row.names=FALSE, quote=FALSE)

write.csv(new2, paste(outdir, "/GWAS_pheno_M1_Veg.csv", sep=""),row.names=FALSE, quote=FALSE)

#x<-as_tibble(read.csv("GWAS_pheno_M1_Veg.csv"))
#pheno <- as_tibble(read.table("/scratch/ahc87874/Check/GWAS_pheno_M1_Veg.txt", header=TRUE))       
       
#=-=-=-=-=-=-=-=-=-=-+-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=       
       
Means <- apply(new2, 2, mean, rm.na = TRUE)
SDs <- apply(new2, 2, sd, rm.na = TRUE)
cbind(Means, SDs)

       
pheno2 <- pheno %>%
       select(LDLCholesterol, HDLCholesterol, TotalCholesterol, Triglycerides)
       
apply(pheno2[pheno$CSRV=="Non-Vegetarian", ], 2, median, na.rm = TRUE)
apply(pheno2[pheno$CSRV=="Vegetarian", ], 2, median, na.rm = TRUE)
apply(pheno2[pheno$SSRV=="Non-Vegetarian", ], 2, median, na.rm = TRUE)
apply(pheno2[pheno$SSRV=="Vegetarian", ], 2, median, na.rm = TRUE)
