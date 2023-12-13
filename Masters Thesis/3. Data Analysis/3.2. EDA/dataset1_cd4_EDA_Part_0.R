# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load Libraries
library(dplyr)
library(ggplot2)
library(viridis)
library(tidyr)
library(RColorBrewer)
library(ungeviz)
library(DescTools)
library(boot)
library(ggpubr)
library(rstatix)
library(ggprism)
library(grid)
library(gridExtra)
library(patchwork)

#---------------------
# Import Data
#---------------------

cd4_counts <- read.csv("clean_cd4_counts_dataset_1.csv", check.names = F)
str(cd4_counts)

#---------------------
# Add Time Category
#---------------------
cd4_counts <- data.frame(timepnt = NA, cd4_counts)

cd4_counts$timepnt[which(cd4_counts$Day==0)] <- "baseline"
cd4_counts$timepnt[which(cd4_counts$Day==14)] <- "peak"
cd4_counts$timepnt[which(cd4_counts$Day==70)] <- "peak"
cd4_counts$timepnt[which(cd4_counts$Day==126)] <- "peak"
cd4_counts$timepnt[which(cd4_counts$Day==224)] <- "memory"
cd4_counts$timepnt[which(cd4_counts$Day==292)] <- "memory"
cd4_counts$timepnt[which(cd4_counts$Day==210)] <- "memory"

which(is.na(cd4_counts$timepnt) == T)

cd4_counts$timepnt <- as.factor(cd4_counts$timepnt)
cd4_counts$timepnt <- factor(cd4_counts$timepnt, levels=c("baseline", "peak", "memory"))
levels(cd4_counts$timepnt)


#---------------------
# Add Dose
#---------------------
# Need to add in dose and schedule information
cd4_counts <- data.frame(Dose = NA, cd4_counts)

cd4_counts[intersect(which(cd4_counts$Study == "H56-035"),which(cd4_counts$Group == 1)),1] <- 5
cd4_counts[intersect(which(cd4_counts$Study == "H56-035"),which(cd4_counts$Group == 2)),1] <- 5
cd4_counts[intersect(which(cd4_counts$Study == "H56-035"),which(cd4_counts$Group == 3)),1] <- 5
cd4_counts[intersect(which(cd4_counts$Study == "H56-035"),which(cd4_counts$Group == 4)),1] <- 5

cd4_counts[intersect(which(cd4_counts$Study == "H1"),which(cd4_counts$Group == 1)),1] <- 15
cd4_counts[intersect(which(cd4_counts$Study == "H1"),which(cd4_counts$Group == 2)),1] <- 50
cd4_counts[intersect(which(cd4_counts$Study == "H1"),which(cd4_counts$Group == 3)),1] <- 15
cd4_counts[intersect(which(cd4_counts$Study == "H1"),which(cd4_counts$Group == 4)),1] <- 50

cd4_counts[intersect(which(cd4_counts$Study == "H56-032"),which(cd4_counts$Group == 1)),1] <- 50
cd4_counts[intersect(which(cd4_counts$Study == "H56-032"),which(cd4_counts$Group == 2)),1] <- 15
cd4_counts[intersect(which(cd4_counts$Study == "H56-032"),which(cd4_counts$Group == 3)),1] <- 50

# check
cd4_counts[which(is.na(cd4_counts$Dose)==T),]
# convert to factor
cd4_counts$Dose <- as.factor(cd4_counts$Dose)


#---------------------
# Add Schedule
#---------------------
# Need to add in dose and schedule information
cd4_counts <- data.frame(Schedule = NA, cd4_counts)

cd4_counts[intersect(which(cd4_counts$Study == "H56-035"),which(cd4_counts$Group == 1)),1] <- 2
cd4_counts[intersect(which(cd4_counts$Study == "H56-035"),which(cd4_counts$Group == 2)),1] <- 3
cd4_counts[intersect(which(cd4_counts$Study == "H56-035"),which(cd4_counts$Group == 3)),1] <- 2
cd4_counts[intersect(which(cd4_counts$Study == "H56-035"),which(cd4_counts$Group == 4)),1] <- 3

cd4_counts[intersect(which(cd4_counts$Study == "H1"),which(cd4_counts$Group == 1)),1] <- 2
cd4_counts[intersect(which(cd4_counts$Study == "H1"),which(cd4_counts$Group == 2)),1] <- 2
cd4_counts[intersect(which(cd4_counts$Study == "H1"),which(cd4_counts$Group == 3)),1] <- 1
cd4_counts[intersect(which(cd4_counts$Study == "H1"),which(cd4_counts$Group == 4)),1] <- 1

cd4_counts[intersect(which(cd4_counts$Study == "H56-032"),which(cd4_counts$Group == 1)),1] <- 3
cd4_counts[intersect(which(cd4_counts$Study == "H56-032"),which(cd4_counts$Group == 2)),1] <- 3
cd4_counts[intersect(which(cd4_counts$Study == "H56-032"),which(cd4_counts$Group == 3)),1] <- 3

# check 
cd4_counts[which(is.na(cd4_counts$Schedule)==T),]
# convert to factor
cd4_counts$Schedule <- as.factor(cd4_counts$Schedule)
head(cd4_counts)


#---------------------
# Tables
#---------------------
# How many PIDs do we have in each conc*schedule combo per Antigen?
# How many observations do we have in each conc*schedule combo per Antigen?

#---------------------------------------------------------------------------------- Number of PIDs:

# Ag85B
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==5 & 
                           cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==5 & 
                           cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT pos","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==5 & 
                           cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==5 & 
                           cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT pos","PID"]))

length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==15 & 
                           cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==15 & 
                           cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT pos","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==15 & 
                           cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==15 & 
                           cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT pos","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==15 & 
                           cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==15 & 
                           cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT pos","PID"]))

length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==50 & 
                           cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==50 & 
                           cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT pos","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==50 & 
                           cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==50 & 
                           cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT pos","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==50 & 
                           cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==50 & 
                           cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT pos","PID"]))




# ESAT6
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==5 & 
                           cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==5 & 
                           cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT pos","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==5 & 
                           cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==5 & 
                           cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT pos","PID"]))

length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==15 & 
                           cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==15 & 
                           cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT pos","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==15 & 
                           cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==15 & 
                           cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT pos","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==15 & 
                           cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==15 & 
                           cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT pos","PID"]))

length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==50 & 
                           cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==50 & 
                           cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT pos","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==50 & 
                           cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==50 & 
                           cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT pos","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==50 & 
                           cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT neg","PID"]))
length(unique(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==50 & 
                           cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT pos","PID"]))






#---------------------------------------------------------------------------------- Number of timepoints:

# Ag85B
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==5 & 
                   cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==5 & 
                   cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT pos","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==5 & 
                   cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==5 & 
                   cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT pos","timepnt"])

table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==15 & 
                   cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==15 & 
                   cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT pos","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==15 & 
                   cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==15 & 
                   cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT pos","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==15 & 
                   cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==15 & 
                   cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT pos","timepnt"])

table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==50 & 
                   cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==50 & 
                   cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT pos","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==50 & 
                   cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==50 & 
                   cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT pos","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==50 & 
                   cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["Dose"]==50 & 
                   cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT pos","timepnt"])



# ESAT6
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==5 & 
                   cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==5 & 
                   cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT pos","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==5 & 
                   cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==5 & 
                   cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT pos","timepnt"])

table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==15 & 
                   cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==15 & 
                   cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT pos","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==15 & 
                   cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==15 & 
                   cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT pos","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==15 & 
                   cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==15 & 
                   cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT pos","timepnt"])

table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==50 & 
                   cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==50 & 
                   cd4_counts["Schedule"]==1  & cd4_counts["QFT"]=="QFT pos","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==50 & 
                   cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==50 & 
                   cd4_counts["Schedule"]==2  & cd4_counts["QFT"]=="QFT pos","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==50 & 
                   cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT neg","timepnt"])
table(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["Dose"]==50 & 
                   cd4_counts["Schedule"]==3  & cd4_counts["QFT"]=="QFT pos","timepnt"])



#------------------------------------------------------- Number of observations for each timepoint
nrow(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["timepnt"]=="baseline",])
nrow(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["timepnt"]=="peak",])
nrow(cd4_counts[cd4_counts["Stim"]=="Ag85B" & cd4_counts["timepnt"]=="memory",])


nrow(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["timepnt"]=="baseline",])
nrow(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["timepnt"]=="peak",])
nrow(cd4_counts[cd4_counts["Stim"]=="ESAT6" & cd4_counts["timepnt"]=="memory",])

