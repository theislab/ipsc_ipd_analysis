###Script to generate plots for depicting primary cilia distribution + statistical tests in Fig. 5,7,8,9 and Supplementary Fig. 7,8,10


setwd("")



library(ggplot2)
library(readxl)
library(sm)
library(devEMF)
library(car)
library(lme4)

##Import data
Ctrl <- c("UKERi1JF-R1", "UKERiG3G-R1", "UKERi1E4-R1", "UKERiO3H-R1", "UKERi82A-R1")
sPD <- c("UKERiJ2C-R1", "UKERiM89-R1", "UKERiC99-R1", "UKERiR66-R1", "UKERiAY6-R1", "UKERiPX7-R1", "UKERi88H-R1")

#Length Table from the Source Data file contain in the 
#first column "Disease_state": information about disease state (Ctrl vs sPD)
#second column "Patient ID": patient identifier (e.g. see above or in column) or "mouse ID"; 
#third column "PC Length [µm]": each cilia length measured per cell; 

length <- read_excel("length.xlsx", col_types = c("text", "text", "numeric"))


##assign variables
length_Ctrl <- as.numeric(na.exclude(length[length$`Patient ID` %in% Ctrl,]$`PC Length [µm]`))
length_sPD <- as.numeric(na.exclude(length[length$`Patient ID` %in% sPD,]$`PC Length [µm]`))


##Density plot and bootstrap hypothesis test
data <- c(length_Ctrl, length_sPD)
group.index <- rep(1:2, c(length(length_Ctrl), length(length_sPD)))

sm.density.compare(data, group = group.index, 
                   lwd = 2, 
                   lty = c(1,1), 
                   col=c("darkblue", "red"), 
                   xlab= "PC length [?m]",
                   size =8,
                   
) #or model = "equal"


##KS test
ks.test(length_Ctrl, length_sPD, alternative = "two.sided")


##Distribution per patient
#Order to table for plot
length$`Patient ID` <- factor(length$`Patient ID`, levels = c("UKERi1JF-R1", "UKERiG3G-R1", "UKERi1E4-R1", "UKERiO3H-R1", "UKERi82A-R1",
                                                         "UKERiJ2C-R1", "UKERiM89-R1", "UKERiC99-R1", "UKERiR66-R1", "UKERiAY6-R1", "UKERiPX7-R1", "UKERi88H-R1"))
#Calculate average primary cilia length of Ctrl individuals
mean <- mean(length[length$Disease_state == "Ctrl",]$`PC Length [µm]`)

#Plot primary cilia length per patient with horizontal line as average Ctrl length
plot <-ggplot(length,aes(x= `Patient ID`, y= `PC Length [µm]`, group = `Patient ID`, color = Disease_state) ) + 
  geom_boxplot() + scale_color_manual(values=c("darkblue" , "red")) + 
  theme_classic() + geom_hline(yintercept=mean) ; plot


##Linear mixed effects model
mod0 <- lme4::lmer(`PC Length [µm]` ~  Disease_state + (1 | Disease_state:`Patient ID`), data = length, REML = FALSE); summary(mod0)
Anova(mod0)
