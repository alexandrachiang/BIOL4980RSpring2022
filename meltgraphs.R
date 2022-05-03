library(tidyverse)
library(plyr)
library(dplyr)
library(ggnewscale)
library(utils)
library(reshape2)

setwd("/scratch/ahc87874/Check")

pheno <- as_tibble(read.csv("GWAS_pheno_M1_Veg.csv"))
pvals <- as_tibble(read.csv("modelpvals.csv"))

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

phenomelt <- melt(pheno, id.vars = c("CSRV", "SSRV", "Sex"), 
                         measure.vars = c("TotalCholesterol", "LDLCholesterol", "HDLCholesterol", "Triglycerides"))

phenomeltsex <- phenomelt
phenomeltsex$Sex <- "Both"
phenomelt<-rbind(phenomelt, phenomeltsex)

phenomelt<-as_tibble(phenomelt)

phenoCSRV<-phenomelt%>%select(CSRV, Sex, variable, value)
phenoSSRV<-phenomelt%>%select(SSRV, Sex, variable, value)

CSRVTotal<-phenoCSRV[phenoCSRV$variable=="TotalCholesterol",]
CSRVLDL<-phenoCSRV[phenoCSRV$variable=="LDLCholesterol",]
CSRVHDL<-phenoCSRV[phenoCSRV$variable=="HDLCholesterol",]
CSRVTri<-phenoCSRV[phenoCSRV$variable=="Triglycerides",]

SSRVTotal<-phenoSSRV[phenoSSRV$variable=="TotalCholesterol",]
SSRVLDL<-phenoSSRV[phenoSSRV$variable=="LDLCholesterol",]
SSRVHDL<-phenoSSRV[phenoSSRV$variable=="HDLCholesterol",]
SSRVTri<-phenoSSRV[phenoSSRV$variable=="Triglycerides",]

setwd("/scratch/ahc87874/Check/graphs")
               
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

##Relationship of Total Cholesterol
##Relationship of CSRV & Total Cholesterol
meannonveg1 <- signif(mean(CSRVTotal$value[CSRVTotal$CSRV=="Non-Vegetarian" & CSRVTotal$Sex=="Both"], na.rm=TRUE), digits = 5)
meanveg1 <- signif(mean(CSRVTotal$value[CSRVTotal$CSRV=="Vegetarian" & CSRVTotal$Sex=="Both"], na.rm=TRUE), digits = 5)
meannonvegF1 <- signif(mean(CSRVTotal$value[CSRVTotal$CSRV=="Non-Vegetarian" & CSRVTotal$Sex=="Female"], na.rm=TRUE), digits = 5)
meanvegF1 <- signif(mean(CSRVTotal$value[CSRVTotal$CSRV=="Vegetarian" & CSRVTotal$Sex=="Female"], na.rm=TRUE), digits = 5)
meannonvegM1 <- signif(mean(CSRVTotal$value[CSRVTotal$CSRV=="Non-Vegetarian" & CSRVTotal$Sex=="Male"], na.rm=TRUE), digits = 5)
meanvegM1 <- signif(mean(CSRVTotal$value[CSRVTotal$CSRV=="Vegetarian" & CSRVTotal$Sex=="Male"], na.rm=TRUE), digits = 5)
    
xlab1 <-paste(levels(CSRVTotal$Sex),"\n(N=",table(CSRVTotal$Sex),")","\n(p=",pvals$CSRVvsTotal,")",sep="")

graphCSRVTotal<- ggplot(data = CSRVTotal, aes(x = Sex, y = value, fill = CSRV, color = CSRV)) + 
  geom_boxplot(alpha = 0.3, position = "identity") + 
    ylim(0, 17) + 
  labs(fill = "CSRV",
       x = "Sex",
       y = "Total Cholesterol (mmol/L)",
       title = "Relationship of CSRV & Total Cholesterol") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
  annotate("text", x = 2.5, y = meannonveg1+1, label = meannonveg1, color = "#F8766D") +
  annotate("text", x = 1.5, y = meannonvegM1+1, label = meannonvegM1, color = "#F8766D") +
  annotate("text", x = 0.5, y = meannonvegF1+1, label = meannonvegF1, color = "#F8766D") +
  annotate("text", x = 2.5, y = meanveg1-1, label = meanveg1, color = "#00BA38") +
  annotate("text", x = 1.5, y = meanvegM1-1, label = meanvegM1, color = "#00BA38") +
  annotate("text", x = 0.5, y = meanvegF1-1, label = meanvegF1, color = "#00BA38") +
  scale_fill_manual(values = c("#F8766D", "#00BA38")) +
  scale_color_manual(values = c("#F8766D", "#00BA38")) +
  scale_x_discrete(labels = xlab1) + 
  coord_flip()

png(filename="CSRVvsTotal.png", type="cairo", width=600, height=300)
graphCSRVTotal
dev.off()

##Relationship of SSRV & Total Cholesterol
meannonveg1 <- signif(mean(SSRVTotal$value[SSRVTotal$SSRV=="Non-Vegetarian" & SSRVTotal$Sex=="Both"], na.rm=TRUE), digits = 5)
meanveg1 <- signif(mean(SSRVTotal$value[SSRVTotal$SSRV=="Vegetarian" & SSRVTotal$Sex=="Both"], na.rm=TRUE), digits = 5)
meannonvegF1 <- signif(mean(SSRVTotal$value[SSRVTotal$SSRV=="Non-Vegetarian" & SSRVTotal$Sex=="Female"], na.rm=TRUE), digits = 5)
meanvegF1 <- signif(mean(SSRVTotal$value[SSRVTotal$SSRV=="Vegetarian" & SSRVTotal$Sex=="Female"], na.rm=TRUE), digits = 5)
meannonvegM1 <- signif(mean(SSRVTotal$value[SSRVTotal$SSRV=="Non-Vegetarian" & SSRVTotal$Sex=="Male"], na.rm=TRUE), digits = 5)
meanvegM1 <- signif(mean(SSRVTotal$value[SSRVTotal$SSRV=="Vegetarian" & SSRVTotal$Sex=="Male"], na.rm=TRUE), digits = 5)
    
xlab1 <-paste(levels(SSRVTotal$Sex),"\n(N=",table(SSRVTotal$Sex),")","\n(p=",pvals$SSRVvsTotal,")",sep="")

graphSSRVTotal<- ggplot(data = SSRVTotal, aes(x = Sex, y = value, fill = SSRV, color = SSRV)) + 
  geom_boxplot(alpha = 0.3, position = "identity") + 
    ylim(0, 17) + 
  labs(fill = "SSRV",
       x = "Sex",
       y = "Total Cholesterol (mmol/L)",
       title = "Relationship of SSRV & Total Cholesterol") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
  annotate("text", x = 2.5, y = meannonveg1+1, label = meannonveg1, color = "#F8766D") +
  annotate("text", x = 1.5, y = meannonvegM1+1, label = meannonvegM1, color = "#F8766D") +
  annotate("text", x = 0.5, y = meannonvegF1+1, label = meannonvegF1, color = "#F8766D") +
  annotate("text", x = 2.5, y = meanveg1-1, label = meanveg1, color = "#00BA38") +
  annotate("text", x = 1.5, y = meanvegM1-1, label = meanvegM1, color = "#00BA38") +
  annotate("text", x = 0.5, y = meanvegF1-1, label = meanvegF1, color = "#00BA38") +
  scale_fill_manual(values = c("#F8766D", "#00BA38")) +
  scale_color_manual(values = c("#F8766D", "#00BA38")) +
  scale_x_discrete(labels = xlab1) + 
  coord_flip()

png(filename="SSRVvsTotal.png", type="cairo", width=600, height=300)
graphSSRVTotal
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

##Relationship of LDL Cholesterol
##Relationship of CSRV & LDL Cholesterol
meannonveg1 <- signif(mean(CSRVLDL$value[CSRVLDL$CSRV=="Non-Vegetarian" & CSRVLDL$Sex=="Both"], na.rm=TRUE), digits = 5)
meanveg1 <- signif(mean(CSRVLDL$value[CSRVLDL$CSRV=="Vegetarian" & CSRVLDL$Sex=="Both"], na.rm=TRUE), digits = 5)
meannonvegF1 <- signif(mean(CSRVLDL$value[CSRVLDL$CSRV=="Non-Vegetarian" & CSRVLDL$Sex=="Female"], na.rm=TRUE), digits = 5)
meanvegF1 <- signif(mean(CSRVLDL$value[CSRVLDL$CSRV=="Vegetarian" & CSRVLDL$Sex=="Female"], na.rm=TRUE), digits = 5)
meannonvegM1 <- signif(mean(CSRVLDL$value[CSRVLDL$CSRV=="Non-Vegetarian" & CSRVLDL$Sex=="Male"], na.rm=TRUE), digits = 5)
meanvegM1 <- signif(mean(CSRVLDL$value[CSRVLDL$CSRV=="Vegetarian" & CSRVLDL$Sex=="Male"], na.rm=TRUE), digits = 5)
    
xlab1 <-paste(levels(CSRVLDL$Sex),"\n(N=",table(CSRVLDL$Sex),")","\n(p=",pvals$CSRVvsLDL,")",sep="")

graphCSRVLDL<- ggplot(data = CSRVLDL, aes(x = Sex, y = value, fill = CSRV, color = CSRV)) + 
  geom_boxplot(alpha = 0.3, position = "identity") + 
    ylim(0, 10) + 
  labs(fill = "CSRV",
       x = "Sex",
       y = "LDL Cholesterol (mmol/L)",
       title = "Relationship of CSRV & LDL Cholesterol") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
  annotate("text", x = 2.5, y = meannonveg1+.6, label = meannonveg1, color = "#F8766D") +
  annotate("text", x = 1.5, y = meannonvegM1+.6, label = meannonvegM1, color = "#F8766D") +
  annotate("text", x = 0.5, y = meannonvegF1+.6, label = meannonvegF1, color = "#F8766D") +
  annotate("text", x = 2.5, y = meanveg1-.6, label = meanveg1, color = "#00BA38") +
  annotate("text", x = 1.5, y = meanvegM1-.6, label = meanvegM1, color = "#00BA38") +
  annotate("text", x = 0.5, y = meanvegF1-.6, label = meanvegF1, color = "#00BA38") +
  scale_fill_manual(values = c("#F8766D", "#00BA38")) +
  scale_color_manual(values = c("#F8766D", "#00BA38")) +
  scale_x_discrete(labels = xlab1) + 
  coord_flip()

png(filename="CSRVvsLDL.png", type="cairo", width=600, height=300)
graphCSRVLDL
dev.off()

##Relationship of SSRV & LDL Cholesterol
meannonveg1 <- signif(mean(SSRVLDL$value[SSRVLDL$SSRV=="Non-Vegetarian" & SSRVLDL$Sex=="Both"], na.rm=TRUE), digits = 5)
meanveg1 <- signif(mean(SSRVLDL$value[SSRVLDL$SSRV=="Vegetarian" & SSRVLDL$Sex=="Both"], na.rm=TRUE), digits = 5)
meannonvegF1 <- signif(mean(SSRVLDL$value[SSRVLDL$SSRV=="Non-Vegetarian" & SSRVLDL$Sex=="Female"], na.rm=TRUE), digits = 5)
meanvegF1 <- signif(mean(SSRVLDL$value[SSRVLDL$SSRV=="Vegetarian" & SSRVLDL$Sex=="Female"], na.rm=TRUE), digits = 5)
meannonvegM1 <- signif(mean(SSRVLDL$value[SSRVLDL$SSRV=="Non-Vegetarian" & SSRVLDL$Sex=="Male"], na.rm=TRUE), digits = 5)
meanvegM1 <- signif(mean(SSRVLDL$value[SSRVLDL$SSRV=="Vegetarian" & SSRVLDL$Sex=="Male"], na.rm=TRUE), digits = 5)
    
xlab1 <-paste(levels(SSRVLDL$Sex),"\n(N=",table(SSRVLDL$Sex),")","\n(p=",pvals$SSRVvsLDL,")",sep="")

graphSSRVLDL<- ggplot(data = SSRVLDL, aes(x = Sex, y = value, fill = SSRV, color = SSRV)) + 
  geom_boxplot(alpha = 0.3, position = "identity") + 
    ylim(0, 10) + 
  labs(fill = "SSRV",
       x = "Sex",
       y = "LDL Cholesterol (mmol/L)",
       title = "Relationship of SSRV & LDL Cholesterol") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
  annotate("text", x = 2.5, y = meannonveg1+1, label = meannonveg1, color = "#F8766D") +
  annotate("text", x = 1.5, y = meannonvegM1+1, label = meannonvegM1, color = "#F8766D") +
  annotate("text", x = 0.5, y = meannonvegF1+1, label = meannonvegF1, color = "#F8766D") +
  annotate("text", x = 2.5, y = meanveg1-1, label = meanveg1, color = "#00BA38") +
  annotate("text", x = 1.5, y = meanvegM1-1, label = meanvegM1, color = "#00BA38") +
  annotate("text", x = 0.5, y = meanvegF1-1, label = meanvegF1, color = "#00BA38") +
  scale_fill_manual(values = c("#F8766D", "#00BA38")) +
  scale_color_manual(values = c("#F8766D", "#00BA38")) +
  scale_x_discrete(labels = xlab1) + 
  coord_flip()

png(filename="SSRVvsLDL.png", type="cairo", width=600, height=300)
graphSSRVLDL
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

##Relationship of HDL Cholesterol
##Relationship of CSRV & HDL Cholesterol
meannonveg1 <- signif(mean(CSRVHDL$value[CSRVHDL$CSRV=="Non-Vegetarian" & CSRVHDL$Sex=="Both"], na.rm=TRUE), digits = 5)
meanveg1 <- signif(mean(CSRVHDL$value[CSRVHDL$CSRV=="Vegetarian" & CSRVHDL$Sex=="Both"], na.rm=TRUE), digits = 5)
meannonvegF1 <- signif(mean(CSRVHDL$value[CSRVHDL$CSRV=="Non-Vegetarian" & CSRVHDL$Sex=="Female"], na.rm=TRUE), digits = 5)
meanvegF1 <- signif(mean(CSRVHDL$value[CSRVHDL$CSRV=="Vegetarian" & CSRVHDL$Sex=="Female"], na.rm=TRUE), digits = 5)
meannonvegM1 <- signif(mean(CSRVHDL$value[CSRVHDL$CSRV=="Non-Vegetarian" & CSRVHDL$Sex=="Male"], na.rm=TRUE), digits = 5)
meanvegM1 <- signif(mean(CSRVHDL$value[CSRVHDL$CSRV=="Vegetarian" & CSRVHDL$Sex=="Male"], na.rm=TRUE), digits = 5)
    
xlab1 <-paste(levels(CSRVHDL$Sex),"\n(N=",table(CSRVHDL$Sex),")","\n(p=",pvals$CSRVvsHDL,")",sep="")

graphCSRVHDL<- ggplot(data = CSRVHDL, aes(x = Sex, y = value, fill = CSRV, color = CSRV)) + 
  geom_boxplot(alpha = 0.3, position = "identity") + 
    ylim(0, 4.3) + 
  labs(fill = "CSRV",
       x = "Sex",
       y = "HDL Cholesterol (mmol/L)",
       title = "Relationship of CSRV & HDL Cholesterol") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
  annotate("text", x = 2.5, y = meannonveg1-.3, label = meannonveg1, color = "#F8766D") +
  annotate("text", x = 1.5, y = meannonvegM1+.3, label = meannonvegM1, color = "#F8766D") +
  annotate("text", x = 0.5, y = meannonvegF1+.3, label = meannonvegF1, color = "#F8766D") +
  annotate("text", x = 2.5, y = meanveg1+.3, label = meanveg1, color = "#00BA38") +
  annotate("text", x = 1.5, y = meanvegM1-.3, label = meanvegM1, color = "#00BA38") +
  annotate("text", x = 0.5, y = meanvegF1-.3, label = meanvegF1, color = "#00BA38") +
  scale_fill_manual(values = c("#F8766D", "#00BA38")) +
  scale_color_manual(values = c("#F8766D", "#00BA38")) +
  scale_x_discrete(labels = xlab1) + 
  coord_flip()

png(filename="CSRVvsHDL.png", type="cairo", width=600, height=300)
graphCSRVHDL
dev.off()

##Relationship of SSRV & HDL Cholesterol
meannonveg1 <- signif(mean(SSRVHDL$value[SSRVHDL$SSRV=="Non-Vegetarian" & SSRVHDL$Sex=="Both"], na.rm=TRUE), digits = 5)
meanveg1 <- signif(mean(SSRVHDL$value[SSRVHDL$SSRV=="Vegetarian" & SSRVHDL$Sex=="Both"], na.rm=TRUE), digits = 5)
meannonvegF1 <- signif(mean(SSRVHDL$value[SSRVHDL$SSRV=="Non-Vegetarian" & SSRVHDL$Sex=="Female"], na.rm=TRUE), digits = 5)
meanvegF1 <- signif(mean(SSRVHDL$value[SSRVHDL$SSRV=="Vegetarian" & SSRVHDL$Sex=="Female"], na.rm=TRUE), digits = 5)
meannonvegM1 <- signif(mean(SSRVHDL$value[SSRVHDL$SSRV=="Non-Vegetarian" & SSRVHDL$Sex=="Male"], na.rm=TRUE), digits = 5)
meanvegM1 <- signif(mean(SSRVHDL$value[SSRVHDL$SSRV=="Vegetarian" & SSRVHDL$Sex=="Male"], na.rm=TRUE), digits = 5)
    
xlab1 <-paste(levels(SSRVHDL$Sex),"\n(N=",table(SSRVHDL$Sex),")","\n(p=",pvals$SSRVvsHDL,")",sep="")

graphSSRVHDL<- ggplot(data = SSRVHDL, aes(x = Sex, y = value, fill = SSRV, color = SSRV)) + 
  geom_boxplot(alpha = 0.3, position = "identity") + 
    ylim(0, 4.3) + 
  labs(fill = "SSRV",
       x = "Sex",
       y = "HDL Cholesterol (mmol/L)",
       title = "Relationship of SSRV & HDL Cholesterol") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
  annotate("text", x = 2.5, y = meannonveg1-.3, label = meannonveg1, color = "#F8766D") +
  annotate("text", x = 1.5, y = meannonvegM1+.3, label = meannonvegM1, color = "#F8766D") +
  annotate("text", x = 0.5, y = meannonvegF1+.3, label = meannonvegF1, color = "#F8766D") +
  annotate("text", x = 2.5, y = meanveg1+.3, label = meanveg1, color = "#00BA38") +
  annotate("text", x = 1.5, y = meanvegM1-.3, label = meanvegM1, color = "#00BA38") +
  annotate("text", x = 0.5, y = meanvegF1-.3, label = meanvegF1, color = "#00BA38") +
  scale_fill_manual(values = c("#F8766D", "#00BA38")) +
  scale_color_manual(values = c("#F8766D", "#00BA38")) +
  scale_x_discrete(labels = xlab1) + 
  coord_flip()

png(filename="SSRVvsHDL.png", type="cairo", width=600, height=300)
graphSSRVHDL
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

##Relationship of Triglyceride
##Relationship of CSRV & Triglyceride
meannonveg1 <- signif(mean(CSRVTri$value[CSRVTri$CSRV=="Non-Vegetarian" & CSRVTri$Sex=="Both"], na.rm=TRUE), digits = 5)
meanveg1 <- signif(mean(CSRVTri$value[CSRVTri$CSRV=="Vegetarian" & CSRVTri$Sex=="Both"], na.rm=TRUE), digits = 5)
meannonvegF1 <- signif(mean(CSRVTri$value[CSRVTri$CSRV=="Non-Vegetarian" & CSRVTri$Sex=="Female"], na.rm=TRUE), digits = 5)
meanvegF1 <- signif(mean(CSRVTri$value[CSRVTri$CSRV=="Vegetarian" & CSRVTri$Sex=="Female"], na.rm=TRUE), digits = 5)
meannonvegM1 <- signif(mean(CSRVTri$value[CSRVTri$CSRV=="Non-Vegetarian" & CSRVTri$Sex=="Male"], na.rm=TRUE), digits = 5)
meanvegM1 <- signif(mean(CSRVTri$value[CSRVTri$CSRV=="Vegetarian" & CSRVTri$Sex=="Male"], na.rm=TRUE), digits = 5)
    
xlab1 <-paste(levels(CSRVTri$Sex),"\n(N=",table(CSRVTri$Sex),")","\n(p=",pvals$CSRVvsTri,")",sep="")

graphCSRVTri<- ggplot(data = CSRVTri, aes(x = Sex, y = value, fill = CSRV, color = CSRV)) + 
  geom_boxplot(alpha = 0.3, position = "identity") + 
    ylim(0, 12) + 
  labs(fill = "CSRV",
       x = "Sex",
       y = "Triglyceride (mmol/L)",
       title = "Relationship of CSRV & Triglyceride") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
  annotate("text", x = 2.5, y = meannonveg1+1, label = meannonveg1, color = "#F8766D") +
  annotate("text", x = 1.5, y = meannonvegM1-1, label = meannonvegM1, color = "#F8766D") +
  annotate("text", x = 0.5, y = meannonvegF1+1, label = meannonvegF1, color = "#F8766D") +
  annotate("text", x = 2.5, y = meanveg1-1, label = meanveg1, color = "#00BA38") +
  annotate("text", x = 1.5, y = meanvegM1+1, label = meanvegM1, color = "#00BA38") +
  annotate("text", x = 0.5, y = meanvegF1-1, label = meanvegF1, color = "#00BA38") +
  scale_fill_manual(values = c("#F8766D", "#00BA38")) +
  scale_color_manual(values = c("#F8766D", "#00BA38")) +
  scale_x_discrete(labels = xlab1) + 
  coord_flip()

png(filename="CSRVvsTri.png", type="cairo", width=600, height=300)
graphCSRVTri
dev.off()

##Relationship of SSRV & Triglyceride
meannonveg1 <- signif(mean(SSRVTri$value[SSRVTri$SSRV=="Non-Vegetarian" & SSRVTri$Sex=="Both"], na.rm=TRUE), digits = 5)
meanveg1 <- signif(mean(SSRVTri$value[SSRVTri$SSRV=="Vegetarian" & SSRVTri$Sex=="Both"], na.rm=TRUE), digits = 5)
meannonvegF1 <- signif(mean(SSRVTri$value[SSRVTri$SSRV=="Non-Vegetarian" & SSRVTri$Sex=="Female"], na.rm=TRUE), digits = 5)
meanvegF1 <- signif(mean(SSRVTri$value[SSRVTri$SSRV=="Vegetarian" & SSRVTri$Sex=="Female"], na.rm=TRUE), digits = 5)
meannonvegM1 <- signif(mean(SSRVTri$value[SSRVTri$SSRV=="Non-Vegetarian" & SSRVTri$Sex=="Male"], na.rm=TRUE), digits = 5)
meanvegM1 <- signif(mean(SSRVTri$value[SSRVTri$SSRV=="Vegetarian" & SSRVTri$Sex=="Male"], na.rm=TRUE), digits = 5)
    
xlab1 <-paste(levels(SSRVTri$Sex),"\n(N=",table(SSRVTri$Sex),")","\n(p=",pvals$SSRVvsTri,")",sep="")

graphSSRVTri<- ggplot(data = SSRVTri, aes(x = Sex, y = value, fill = SSRV, color = SSRV)) + 
  geom_boxplot(alpha = 0.3, position = "identity") + 
    ylim(0, 12) + 
  labs(fill = "SSRV",
       x = "Sex",
       y = "Triglyceride (mmol/L)",
       title = "Relationship of SSRV & Triglyceride") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3) +
  annotate("text", x = 2.5, y = meannonveg1+1, label = meannonveg1, color = "#F8766D") +
  annotate("text", x = 1.5, y = meannonvegM1-1, label = meannonvegM1, color = "#F8766D") +
  annotate("text", x = 0.5, y = meannonvegF1+1, label = meannonvegF1, color = "#F8766D") +
  annotate("text", x = 2.5, y = meanveg1-1, label = meanveg1, color = "#00BA38") +
  annotate("text", x = 1.5, y = meanvegM1+1, label = meanvegM1, color = "#00BA38") +
  annotate("text", x = 0.5, y = meanvegF1-1, label = meanvegF1, color = "#00BA38") +
  scale_fill_manual(values = c("#F8766D", "#00BA38")) +
  scale_color_manual(values = c("#F8766D", "#00BA38")) +
  scale_x_discrete(labels = xlab1) + 
  coord_flip()

png(filename="SSRVvsTri.png", type="cairo", width=600, height=300)
graphSSRVTri
dev.off()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
