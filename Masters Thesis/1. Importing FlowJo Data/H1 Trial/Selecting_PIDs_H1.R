# Clear environment
rm(list=ls())

# Clear all plots
if(!is.null(dev.list())) dev.off()


#---------------------
# Load Packages
#---------------------
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(viridis)
library("readxl")
library(ggstatsplot)
library(gridExtra)

#---------------------
# Load Data
#---------------------

cd4_counts <- read.csv("cd4_counts.csv", check.names = F)
cd4_freq   <- read.csv("cd4_frequency.csv", check.names = F)
cd8_counts <- read.csv("cd8_counts.csv", check.names = F)
cd8_freq   <- read.csv("cd8_frequency.csv", check.names = F)



# Set Groups as factors
cd4_counts$Group <- as.factor(cd4_counts$Group)
cd4_freq$Group   <- as.factor(cd4_freq$Group)
cd8_counts$Group <- as.factor(cd8_counts$Group)
cd8_freq$Group   <- as.factor(cd8_freq$Group)



# import gender data
gender_dat <- read_excel("THYB04_Demographics05102022.xlsx")
head(gender_dat)
gender_dat <- as.data.frame(gender_dat[,c(1,4)])
gender_dat$ParticipantNo <- as.numeric(gender_dat$ParticipantNo)
head(gender_dat)
gender_dat$ParticipantNo
dim(gender_dat)
max(gender_dat$ParticipantNo)  # so particpant numbers == row numbers
gender_dat$Gender <- as.factor(gender_dat$Gender)


#---------------------
# Random Sampling
#---------------------


# We want 15 people in each group*QFT stratum
# And we want to have approximately 50:50 male and female in each stratum
# "Group 1 QFTneg"
table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[1])])])
# "Group 1 QFTpos"
table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[2])])]) # use whole dataset
# "Group 2 QFTneg"
table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[3])])]) # use all men
# "Group 2 QFTpos" 
table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[4])])])
# "Group 3 QFTneg"
table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[5])])])
# "Group 3 QFTpos"
table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[6])])]) # use all men
# "Group 4 QFTneg"
table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[7])])]) # use all men
# "Group 4 QFTpos"
table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[8])])]) # use all men


# Only three groups for which we need to sample men
# Sample women of all groups


set.seed(1234)
n <- sample(c(7,8), 1)  # randomly choose the number of women to pick
f <- unique(gender_dat$ParticipantNo[which(gender_dat$Gender=="Female")])           # everyone who is female
m <- unique(gender_dat$ParticipantNo[which(gender_dat$Gender=="Male")])             # everyone who is male
a <- unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[1])])   # everyone in stratum 1 
Fgrp1_QFTneg_samp <- sample(intersect(a,f),n)
Mgrp1_QFTneg_samp <- sample(intersect(a,m),(15-n))
grp1_QFTneg_samp <- c(Fgrp1_QFTneg_samp, Mgrp1_QFTneg_samp)

grp1_QFTpos_samp <- unique(cd4_counts$PID[which(cd4_counts$Group == "Group 1 QFTpos")]) ## not enough pids to sample

a <- unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[3])])  # everyone in stratum 3
n <- 15-length(intersect(a,m))     # number of women (using all men)
Fgrp2_QFTneg_samp <- sample(intersect(a,f),n)
Mgrp2_QFTneg_samp <- intersect(a,m)
grp2_QFTneg_samp  <- c(Fgrp2_QFTneg_samp, Mgrp2_QFTneg_samp)

n <- sample(c(7,8), 1)  # randomly choose the number of women to pick
a <- unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[4])])
Fgrp2_QFTpos_samp <- sample(intersect(a,f),n)
Mgrp2_QFTpos_samp <- sample(intersect(a,m),(15-n))
grp2_QFTpos_samp  <- c(Fgrp2_QFTpos_samp, Mgrp2_QFTpos_samp)

n <- sample(c(7,8), 1)  # randomly choose the number of women to pick
a <- unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[5])])
Fgrp3_QFTneg_samp <- sample(intersect(a,f),n)
Mgrp3_QFTneg_samp <- sample(intersect(a,m),(15-n))
grp3_QFTneg_samp  <- c(Fgrp3_QFTneg_samp, Mgrp3_QFTneg_samp)

a <- unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[6])])  
n <- 15-length(intersect(a,m))     # number of women (using all men)
Fgrp3_QFTpos_samp <- sample(intersect(a,f),n)
Mgrp3_QFTpos_samp <- intersect(a,m)
grp3_QFTpos_samp  <- c(Fgrp3_QFTpos_samp, Mgrp3_QFTpos_samp)

a <- unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[7])])  
n <- 15-length(intersect(a,m))     # number of women (using all men)
Fgrp4_QFTneg_samp <- sample(intersect(a,f),n)
Mgrp4_QFTneg_samp <- intersect(a,m)
grp4_QFTneg_samp  <- c(Fgrp4_QFTneg_samp, Mgrp4_QFTneg_samp)

a <- unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[8])])  
n <- 15-length(intersect(a,m))     # number of women (using all men)
Fgrp4_QFTpos_samp <- sample(intersect(a,f),n)
Mgrp4_QFTpos_samp <- intersect(a,m)
grp4_QFTpos_samp  <- c(Fgrp4_QFTpos_samp, Mgrp4_QFTpos_samp)



samp <- c(grp1_QFTneg_samp, grp1_QFTpos_samp, grp2_QFTneg_samp,grp2_QFTpos_samp,
          grp3_QFTneg_samp, grp3_QFTpos_samp, grp4_QFTneg_samp, grp4_QFTpos_samp)
length(unique(samp)) == (15*7 + 13)                  # check


cd4_counts_samp <- cd4_counts[which(apply(data.frame(cd4_counts$PID), 1, function(r) any(r %in% samp))),] 
length(unique(cd4_counts_samp$PID)) == (15*7 + 13)   # check

cd4_freq_samp <- cd4_freq[which(apply(data.frame(cd4_freq$PID), 1, function(r) any(r %in% samp))),] 
length(unique(cd4_freq_samp$PID)) == (15*7 + 13)     # check

cd8_counts_samp <- cd8_counts[which(apply(data.frame(cd8_counts$PID), 1, function(r) any(r %in% samp))),] 
length(unique(cd8_counts_samp$PID)) == (15*7 + 13)    # check

cd8_freq_samp <- cd8_freq[which(apply(data.frame(cd8_freq$PID), 1, function(r) any(r %in% samp))),] 
length(unique(cd8_freq_samp$PID)) == (15*7 + 13)      # check




#---------------------
# Check Sampling
#---------------------
# We need to check that the distributions of the samples are the same as that of the full datasets
# Need to compare distributions of Net Total Response

#--------------------------------------------------------- Net Total Response
# Net Total Response is expressed as a frequency.
# We calculate this by summing all cells with at least one marker - so sum along rows but ignore last column (the 
# column counting the number of cells that express none of the markers).
# We then divide this by the total number of CD4 counts and express as a percentage.
# Then subtract UNS from Stim.
# Note that we are only interested in the distributions of the peak visit, i.e. Day 70 or 14


# Extract only Day 70 and Day 14 observations
cd4_samp_ana <- cd4_counts_samp[c(which(cd4_counts_samp$Day == 14),which(cd4_counts_samp$Day == 70)),]   # Sampled data
cd8_samp_ana <- cd8_counts_samp[c(which(cd8_counts_samp$Day == 14),which(cd8_counts_samp$Day == 70)),]
cd4_ana <- cd4_counts[c(which(cd4_counts$Day == 14),which(cd4_counts$Day == 70)),]                  # Full data
cd8_ana <- cd8_counts[c(which(cd8_counts$Day == 14),which(cd8_counts$Day == 70)),]


# Only keep counts for cells expressing at least one of INFg, IL2,TNF or L17
cd4_samp_ana <- cd4_samp_ana[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd4_samp_ana), fixed = TRUE)==T)]
cd8_samp_ana <- cd8_samp_ana[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd8_samp_ana), fixed = TRUE)==T)]
cd4_ana <- cd4_ana[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd4_ana), fixed = TRUE)==T)]
cd8_ana <- cd8_ana[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd8_ana), fixed = TRUE)==T)]


# Sum across rows
cd4_samp_ana <- cd4_samp_ana %>%
  mutate(Total = select(., 6:65) %>% rowSums(na.rm = TRUE))
cd8_samp_ana <- cd8_samp_ana %>%
  mutate(Total = select(., 6:65) %>% rowSums(na.rm = TRUE))
cd4_ana <- cd4_ana %>%
  mutate(Total = select(., 6:65) %>% rowSums(na.rm = TRUE))
cd8_ana <- cd8_ana %>%
  mutate(Total = select(., 6:65) %>% rowSums(na.rm = TRUE))


# Convert to frequency
cd4_samp_ana$Total <- cd4_samp_ana$Total/cd4_samp_ana$CD4count
cd8_samp_ana$Total <- cd8_samp_ana$Total/cd8_samp_ana$CD8count
cd4_ana$Total <- cd4_ana$Total/cd4_ana$CD4count
cd8_ana$Total <- cd8_ana$Total/cd8_ana$CD8count


# Remove unnecessary rows
cd4_samp_ana <- cd4_samp_ana[,c(1,2,3,4,66)] 
cd8_samp_ana <- cd8_samp_ana[,c(1,2,3,4,66)] 
cd4_ana <- cd4_ana[,c(1,2,3,4,66)] 
cd8_ana <- cd8_ana[,c(1,2,3,4,66)] 


# Subtract UNS from Stim
cd4_samp_ana_w <- reshape(cd4_samp_ana, idvar = c("Group", "PID", "Day"), timevar = "Stim", direction = "wide")
head(cd4_samp_ana_w)
cd4_samp_ana_w_bs <- cd4_samp_ana_w   # background subtract sample
cd4_samp_ana_w_bs[,4:6] <- cd4_samp_ana_w_bs[,4:6] - cd4_samp_ana_w_bs$Total.UNS
cd4_samp_ana_w_bs <- cd4_samp_ana_w_bs[,-7]
names(cd4_samp_ana_w_bs) <- c("Group","PID","Day","Ag85B_bs", "Esat6_bs", "PHA_bs")
head(cd4_samp_ana_w_bs)

cd8_samp_ana_w <- reshape(cd8_samp_ana, idvar = c("Group", "PID", "Day"), timevar = "Stim", direction = "wide")
head(cd8_samp_ana_w)
cd8_samp_ana_w_bs <- cd8_samp_ana_w   # background subtract sample
cd8_samp_ana_w_bs[,4:6] <- cd8_samp_ana_w_bs[,4:6] - cd8_samp_ana_w_bs$Total.UNS
cd8_samp_ana_w_bs <- cd8_samp_ana_w_bs[,-7]
names(cd8_samp_ana_w_bs) <- c("Group","PID","Day","Ag85B_bs", "Esat6_bs", "PHA_bs")
head(cd8_samp_ana_w_bs)

cd4_ana_w <- reshape(cd4_ana, idvar = c("Group", "PID", "Day"), timevar = "Stim", direction = "wide")
head(cd4_ana_w)
cd4_ana_w_bs <- cd4_ana_w   # background subtract sample
cd4_ana_w_bs[,4:6] <- cd4_ana_w_bs[,4:6] - cd4_ana_w_bs$Total.UNS
cd4_ana_w_bs <- cd4_ana_w_bs[,-7]
names(cd4_ana_w_bs) <- c("Group","PID","Day","Ag85B_bs", "Esat6_bs", "PHA_bs")
head(cd4_ana_w_bs)

cd8_ana_w <- reshape(cd8_ana, idvar = c("Group", "PID", "Day"), timevar = "Stim", direction = "wide")
head(cd8_ana_w)
cd8_ana_w_bs <- cd8_ana_w   # background subtract sample
cd8_ana_w_bs[,4:6] <- cd8_ana_w_bs[,4:6] - cd8_ana_w_bs$Total.UNS
cd8_ana_w_bs <- cd8_ana_w_bs[,-7]
names(cd8_ana_w_bs) <- c("Group","PID","Day","Ag85B_bs", "Esat6_bs", "PHA_bs")
head(cd8_ana_w_bs)


# Convert proportions to frequencies (i.e. percentages)
cd4_samp_ana_w_bs[,4:6] <- cd4_samp_ana_w_bs[,4:6]*100
cd8_samp_ana_w_bs[,4:6] <- cd8_samp_ana_w_bs[,4:6]*100
cd4_ana_w_bs[,4:6] <- cd4_ana_w_bs[,4:6]*100
cd8_ana_w_bs[,4:6] <- cd8_ana_w_bs[,4:6]*100






##########################################################################################
#                                CD4
###########################################################################################


#--------------------------------------------------------- BOX PLOTS


# Group 1 QFT-
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[1]),])), 
          rep("Ag85B Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[1]),])), 
          rep("Esat6 Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[1]),])), 
          rep("Esat6 Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[1]),]))),
  value=c(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[1]),4], 
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[1]),4], 
          cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[1]),5],
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[1]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD4 Group 1 QFT-") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))




#----------------------------------------------------------------------------------------------


# Group 2 QFT-
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[3]),])), 
          rep("Ag85B Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[3]),])), 
          rep("Esat6 Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[3]),])), 
          rep("Esat6 Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[3]),]))),
  value=c(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[3]),4], 
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[3]),4], 
          cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[3]),5],
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[3]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD4 Group 2 QFT-") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))



#----------------------------------------------------------------------------------------------


# Group 2 QFT+
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[4]),])), 
          rep("Ag85B Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[4]),])), 
          rep("Esat6 Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[4]),])), 
          rep("Esat6 Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[4]),]))),
  value=c(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[4]),4], 
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[4]),4], 
          cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[4]),5],
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[4]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD4 Group 2 QFT+") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))



#----------------------------------------------------------------------------------------------


# Group 3 QFT-
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[5]),])), 
          rep("Ag85B Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[5]),])), 
          rep("Esat6 Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[5]),])), 
          rep("Esat6 Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[5]),]))),
  value=c(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[5]),4], 
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[5]),4], 
          cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[5]),5],
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[5]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD4 Group 3 QFT-") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))



#----------------------------------------------------------------------------------------------


# Group 3 QFT+
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[6]),])), 
          rep("Ag85B Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[6]),])), 
          rep("Esat6 Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[6]),])), 
          rep("Esat6 Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[6]),]))),
  value=c(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[6]),4], 
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[6]),4], 
          cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[6]),5],
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[6]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD4 Group 3 QFT+") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))




#----------------------------------------------------------------------------------------------


# Group 4 QFT-
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[7]),])), 
          rep("Ag85B Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[7]),])), 
          rep("Esat6 Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[7]),])), 
          rep("Esat6 Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[7]),]))),
  value=c(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[7]),4], 
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[7]),4], 
          cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[7]),5],
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[7]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD4 Group 4 QFT-") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))



#----------------------------------------------------------------------------------------------


# Group 4 QFT+
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[8]),])), 
          rep("Ag85B Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[8]),])), 
          rep("Esat6 Full",nrow(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[8]),])), 
          rep("Esat6 Sample",nrow(cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[8]),]))),
  value=c(cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[8]),4], 
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[8]),4], 
          cd4_ana_w_bs[which(cd4_ana_w_bs$Group == levels(cd4_ana_w_bs$Group)[8]),5],
          cd4_samp_ana_w_bs[which(cd4_samp_ana_w_bs$Group == levels(cd4_samp_ana_w_bs$Group)[8]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD4 Group 4 QFT+") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))






##########################################################################################
#                                CD8
###########################################################################################



# Group 1 QFT-
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[1]),])), 
          rep("Ag85B Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[1]),])), 
          rep("Esat6 Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[1]),])), 
          rep("Esat6 Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[1]),]))),
  value=c(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[1]),4], 
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[1]),4], 
          cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[1]),5],
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[1]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD8 Group 1 QFT-") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))



#----------------------------------------------------------------------------------------------


# Group 2 QFT-
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[3]),])), 
          rep("Ag85B Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[3]),])), 
          rep("Esat6 Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[3]),])), 
          rep("Esat6 Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[3]),]))),
  value=c(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[3]),4], 
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[3]),4], 
          cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[3]),5],
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[3]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD8 Group 2 QFT-") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))



#----------------------------------------------------------------------------------------------


# Group 2 QFT+
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[4]),])), 
          rep("Ag85B Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[4]),])), 
          rep("Esat6 Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[4]),])), 
          rep("Esat6 Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[4]),]))),
  value=c(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[4]),4], 
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[4]),4], 
          cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[4]),5],
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[4]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD8 Group 2 QFT+") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))



#----------------------------------------------------------------------------------------------


# Group 3 QFT-
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[5]),])), 
          rep("Ag85B Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[5]),])), 
          rep("Esat6 Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[5]),])), 
          rep("Esat6 Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[5]),]))),
  value=c(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[5]),4], 
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[5]),4], 
          cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[5]),5],
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[5]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD8 Group 3 QFT-") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))




#----------------------------------------------------------------------------------------------


# Group 3 QFT+
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[6]),])), 
          rep("Ag85B Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[6]),])), 
          rep("Esat6 Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[6]),])), 
          rep("Esat6 Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[6]),]))),
  value=c(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[6]),4], 
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[6]),4], 
          cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[6]),5],
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[6]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD8 Group 3 QFT+") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))



#----------------------------------------------------------------------------------------------


# Group 4 QFT-
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[7]),])), 
          rep("Ag85B Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[7]),])), 
          rep("Esat6 Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[7]),])), 
          rep("Esat6 Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[7]),]))),
  value=c(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[7]),4], 
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[7]),4], 
          cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[7]),5],
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[7]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD8 Group 4 QFT-") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))



#----------------------------------------------------------------------------------------------


# Group 4 QFT+
data <- data.frame(
  name=c( rep("Ag85B Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[8]),])), 
          rep("Ag85B Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[8]),])), 
          rep("Esat6 Full",nrow(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[8]),])), 
          rep("Esat6 Sample",nrow(cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[8]),]))),
  value=c(cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[8]),4], 
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[8]),4], 
          cd8_ana_w_bs[which(cd8_ana_w_bs$Group == levels(cd8_ana_w_bs$Group)[8]),5],
          cd8_samp_ana_w_bs[which(cd8_samp_ana_w_bs$Group == levels(cd8_samp_ana_w_bs$Group)[8]),5]))
# Plot
data %>%
  ggplot( aes(x=name, y=value, fill=name)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_minimal()+
  theme(
    legend.position="none",
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)) +
  ggtitle("CD8 Group 4 QFT+") +
  ylab("Frequency (%)")+
  xlab("")+
  scale_fill_manual(values=c("#FCBBA1", "#FC9272", "#CCECE6", "#99D8C9"))





#--------------------------
# CHECK GENDER
#-------------------------


df <- as.data.frame(table(gender_dat$Gender[grp1_QFTneg_samp]))
levels(df$Var1) <- c("F", "M")
piechart1 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 1 QFT- (Sample)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40)) +
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
df <- as.data.frame(table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[1])])]))
levels(df$Var1) <- c("F", "M")
piechart2 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 1 QFT- (Full)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40))+
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
grid.arrange(piechart1, piechart2, ncol=2)




df <- as.data.frame(table(gender_dat$Gender[grp2_QFTneg_samp]))
levels(df$Var1) <- c("F", "M")
piechart1 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 2 QFT- (Sample)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40)) +
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
df <- as.data.frame(table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[3])])]))
levels(df$Var1) <- c("F", "M")
piechart2 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 2 QFT- (Full)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40))+
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
grid.arrange(piechart1, piechart2, ncol=2)




df <- as.data.frame(table(gender_dat$Gender[grp2_QFTpos_samp]))
levels(df$Var1) <- c("F", "M")
piechart1 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 2 QFT+ (Sample)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40)) +
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
df <- as.data.frame(table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[4])])]))
levels(df$Var1) <- c("F", "M")
piechart2 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 2 QFT+ (Full)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40))+
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
grid.arrange(piechart1, piechart2, ncol=2)



df <- as.data.frame(table(gender_dat$Gender[grp3_QFTneg_samp]))
levels(df$Var1) <- c("F", "M")
piechart1 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 3 QFT- (Sample)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40)) +
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
df <- as.data.frame(table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[5])])]))
levels(df$Var1) <- c("F", "M")
piechart2 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 3 QFT- (Full)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40))+
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
grid.arrange(piechart1, piechart2, ncol=2)



df <- as.data.frame(table(gender_dat$Gender[grp3_QFTpos_samp]))
levels(df$Var1) <- c("F", "M")
piechart1 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 3 QFT+ (Sample)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40)) +
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
df <- as.data.frame(table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[6])])]))
levels(df$Var1) <- c("F", "M")
piechart2 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 3 QFT+ (Full)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40))+
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
grid.arrange(piechart1, piechart2, ncol=2)




df <- as.data.frame(table(gender_dat$Gender[grp4_QFTneg_samp]))
levels(df$Var1) <- c("F", "M")
piechart1 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 4 QFT- (Sample)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40)) +
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
df <- as.data.frame(table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[7])])]))
levels(df$Var1) <- c("F", "M")
piechart2 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 4 QFT- (Full)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40))+
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
grid.arrange(piechart1, piechart2, ncol=2)



df <- as.data.frame(table(gender_dat$Gender[grp4_QFTpos_samp]))
levels(df$Var1) <- c("F", "M")
piechart1 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 4 QFT+ (Sample)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40)) +
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
df <- as.data.frame(table(gender_dat$Gender[unique(cd4_counts$PID[which(cd4_counts$Group==levels(cd4_counts$Group)[8])])]))
levels(df$Var1) <- c("F", "M")
piechart2 <- ggplot(df, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(width=1, stat="identity") +
  coord_polar("y", start=0) +
  xlab("") +
  ylab("Group 4 QFT+ (Full)") +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(), text = element_text(size = 40))+
  geom_text(aes(label = paste0(Freq,
                               " (",
                               scales::percent(Freq / sum(Freq)),
                               ")")),
            position = position_stack(vjust = 0.5),  size=15)+
  scale_fill_manual(values=c("#FC9272", "#99D8C9"), name="")
grid.arrange(piechart1, piechart2, ncol=2)













#--------------------------
# TABLE OF SAMPLE PIDS
#-------------------------

samp <- c(sort(grp1_QFTneg_samp), sort(grp1_QFTpos_samp), sort(grp2_QFTneg_samp), sort(grp2_QFTpos_samp),
          sort(grp3_QFTneg_samp), sort(grp3_QFTpos_samp), sort(grp4_QFTneg_samp), sort(grp4_QFTpos_samp))
tab <- data.frame(PID = samp,
                  Group = c(rep("Group 1 QFT-",length(grp1_QFTneg_samp)), rep("Group 1 QFT+",length(grp1_QFTpos_samp)),
                            rep("Group 2 QFT-",length(grp2_QFTneg_samp)), rep("Group 2 QFT+",length(grp2_QFTpos_samp)),
                            rep("Group 3 QFT-",length(grp3_QFTneg_samp)), rep("Group 3 QFT+",length(grp3_QFTpos_samp)),
                            rep("Group 4 QFT-",length(grp4_QFTneg_samp)), rep("Group 4 QFT+",length(grp4_QFTpos_samp))))
tab




#----------------------
# Write .csv files
#----------------------

write.csv(tab, "Sampled_PIDs_per_group_H1.csv", row.names = FALSE) 

write.csv(cd4_counts_samp, "cd4_counts_sample_H1.csv", row.names = FALSE) 
write.csv(cd4_freq_samp, "cd4_frequency_sample_H1.csv", row.names = FALSE) 

write.csv(cd8_counts_samp, "cd8_counts_sample_H1.csv", row.names = FALSE)  
write.csv(cd8_freq_samp, "cd8_frequency_sample_H1.csv", row.names = FALSE) 






