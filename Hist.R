#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#Make histograms of Age, Sex, and BMI
AgeHist <- ggplot(pheno, aes(x = Age, fill = Sex)) +
  geom_histogram(breaks = seq(35, 75, 5), position = "dodge", color = "white") +
  scale_x_continuous(breaks = seq(35, 75, 5)) +
  labs(fill = "Sex",
       x = "Age",
       y = "Count",
       title = "Distribution of Age") + 
  theme(plot.title = element_text(hjust = 0.5))
  
BMIHist <- ggplot(pheno, aes(x = BMI, fill = Sex)) +
  geom_histogram(breaks = seq(10, 75, 5), position = "dodge", color = "white") +
  scale_x_continuous(breaks = seq(10, 75, 5)) +
  labs(fill = "Sex",
       x = "BMI",
       y = "Count",
       title = "Distribution of BMI") + 
  theme(plot.title = element_text(hjust = 0.5))
   
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

setwd("/scratch/ahc87874/Check/graphs")

png(filename="AgeHist.png", type="cairo", width=500, height=300)
AgeHist
dev.off()

png(filename="BMIHist.png", type="cairo", width=500, height=300)
BMIHist
dev.off()
