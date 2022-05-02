library(tidyverse)

setwd("/scratch/ahc87874/Check")

pheno<-as_tibble(read.csv("GWAS_pheno_M1_Veg.csv"))

colnames(pheno)[36]<-"CSRV"
colnames(pheno)[37]<-"SSRV"

pheno$Sex[pheno$Sex==1] <- "Female"
pheno$Sex[pheno$Sex==0] <- "Male"
pheno$Sex <- as.factor(pheno$Sex)
 
pheno$CSRV[pheno$CSRV==0] <- "Non-Vegetarian"
pheno$CSRV[pheno$CSRV==1] <- "Vegetarian"
pheno$CSRV <- as.factor(pheno$CSRV)

pheno$SSRV[pheno$SSRV==0] <- "Non-Vegetarian"
pheno$SSRV[pheno$SSRV==1] <- "Vegetarian"
pheno$SSRV <- as.factor(pheno$SSRV)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

##Total Cholesterol
##Covariates Age, Sex, BMI, CSRV/SSRV
Model1.lm <- lm(TotalCholesterol ~ Age + Sex + BMI + CSRV, data = pheno)
Model2.lm <- lm(TotalCholesterol ~ Age + Sex + BMI + SSRV, data = pheno)

##LDL Cholesterol
##Covariates Age, Sex, BMI, CSRV/SSRV
Model3.lm <- lm(LDLCholesterol ~ Age + Sex + BMI + CSRV, data = pheno)
Model4.lm <- lm(LDLCholesterol ~ Age + Sex + BMI + SSRV, data = pheno)

##HDL Cholesterol
##Covariates Age, Sex, BMI, CSRV/SSRV
Model5.lm <- lm(HDLCholesterol ~ Age + Sex + BMI + CSRV, data = pheno)
Model6.lm <- lm(HDLCholesterol ~ Age + Sex + BMI + SSRV, data = pheno)

##Triglycerides
##Covariates Age, Sex, BMI, CSRV/SSRV
Model7.lm <- lm(Triglycerides ~ Age + Sex + BMI + CSRV, data = pheno)
Model8.lm <- lm(Triglycerides ~ Age + Sex + BMI + SSRV, data = pheno)

summary(Model1.lm)
summary(Model2.lm)
summary(Model3.lm)
summary(Model4.lm)
summary(Model5.lm)
summary(Model6.lm)
summary(Model7.lm)
summary(Model8.lm)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

##Total Cholesterol
##Covariates Age, Sex, CSRV/SSRV
Model1.lm <- lm(TotalCholesterol ~ Age + Sex + CSRV, data = pheno)
Model2.lm <- lm(TotalCholesterol ~ Age + Sex + SSRV, data = pheno)

##LDL Cholesterol
##Covariates Age, Sex, CSRV/SSRV
Model3.lm <- lm(LDLCholesterol ~ Age + Sex + CSRV, data = pheno)
Model4.lm <- lm(LDLCholesterol ~ Age + Sex + SSRV, data = pheno)

##HDL Cholesterol
##Covariates Age, Sex, CSRV/SSRV
Model5.lm <- lm(HDLCholesterol ~ Age + Sex + CSRV, data = pheno)
Model6.lm <- lm(HDLCholesterol ~ Age + Sex + SSRV, data = pheno)

##Triglycerides
##Covariates Age, Sex, CSRV/SSRV
Model7.lm <- lm(Triglycerides ~ Age + Sex + CSRV, data = pheno)
Model8.lm <- lm(Triglycerides ~ Age + Sex + SSRV, data = pheno)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#Males Only
##Total Cholesterol
##Covariates Age, Sex, BMI, CSRV/SSRV
Model1.lm <- lm(TotalCholesterol ~ Age + BMI + CSRV, data = pheno[pheno$Sex=="Male",])
Model2.lm <- lm(TotalCholesterol ~ Age + BMI + SSRV, data = pheno[pheno$Sex=="Male",])

##LDL Cholesterol
##Covariates Age, Sex, BMI, CSRV/SSRV
Model3.lm <- lm(LDLCholesterol ~ Age + BMI + CSRV, data = pheno[pheno$Sex=="Male",])
Model4.lm <- lm(LDLCholesterol ~ Age + BMI + SSRV, data = pheno[pheno$Sex=="Male",])

##HDL Cholesterol
##Covariates Age, Sex, BMI, CSRV/SSRV
Model5.lm <- lm(HDLCholesterol ~ Age + BMI + CSRV, data = pheno[pheno$Sex=="Male",])
Model6.lm <- lm(HDLCholesterol ~ Age + BMI + SSRV, data = pheno[pheno$Sex=="Male",])

##Triglycerides
##Covariates Age, Sex, BMI, CSRV/SSRV
Model7.lm <- lm(Triglycerides ~ Age + BMI + CSRV, data = pheno[pheno$Sex=="Male",])
Model8.lm <- lm(Triglycerides ~ Age + BMI + SSRV, data = pheno[pheno$Sex=="Male",])

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#Females Only
##Total Cholesterol
##Covariates Age, Sex, BMI, CSRV/SSRV
Model1.lm <- lm(TotalCholesterol ~ Age + BMI + CSRV, data = pheno[pheno$Sex=="Female",])
Model2.lm <- lm(TotalCholesterol ~ Age + BMI + SSRV, data = pheno[pheno$Sex=="Female",])

##LDL Cholesterol
##Covariates Age, Sex, BMI, CSRV/SSRV
Model3.lm <- lm(LDLCholesterol ~ Age + BMI + CSRV, data = pheno[pheno$Sex=="Female",])
Model4.lm <- lm(LDLCholesterol ~ Age + BMI + SSRV, data = pheno[pheno$Sex=="Female",])

##HDL Cholesterol
##Covariates Age, Sex, BMI, CSRV/SSRV
Model5.lm <- lm(HDLCholesterol ~ Age + BMI + CSRV, data = pheno[pheno$Sex=="Female",])
Model6.lm <- lm(HDLCholesterol ~ Age + BMI + SSRV, data = pheno[pheno$Sex=="Female",])

##Triglycerides
##Covariates Age, Sex, BMI, CSRV/SSRV
Model7.lm <- lm(Triglycerides ~ Age + BMI + CSRV, data = pheno[pheno$Sex=="Female",])
Model8.lm <- lm(Triglycerides ~ Age + BMI + SSRV, data = pheno[pheno$Sex=="Female",])
