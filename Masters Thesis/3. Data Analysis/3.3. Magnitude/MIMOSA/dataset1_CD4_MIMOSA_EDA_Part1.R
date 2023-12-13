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
library(DescTools)
library(boot)
library(ggpubr)
library(rstatix)
library(ggprism)
library(grid)
library(gridExtra)
library(patchwork)
library(Rmisc) # to get CI



#---------------------
# Import Data
#---------------------

cd4_counts <- read.csv("Filtered_cd4_counts_dataset_1.csv", check.names = F)
str(cd4_counts)


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


#---------------------
# Import Old Data
#---------------------

cd4_counts_old <- read.csv("clean_cd4_counts_dataset_1.csv", check.names = F)
str(cd4_counts_old)

#---------------------
# Add Time Category
#---------------------
cd4_counts_old <- data.frame(timepnt = NA, cd4_counts_old)

cd4_counts_old$timepnt[which(cd4_counts_old$Day==0)] <- "baseline"
cd4_counts_old$timepnt[which(cd4_counts_old$Day==14)] <- "peak"
cd4_counts_old$timepnt[which(cd4_counts_old$Day==70)] <- "peak"
cd4_counts_old$timepnt[which(cd4_counts_old$Day==126)] <- "peak"
cd4_counts_old$timepnt[which(cd4_counts_old$Day==224)] <- "memory"
cd4_counts_old$timepnt[which(cd4_counts_old$Day==292)] <- "memory"
cd4_counts_old$timepnt[which(cd4_counts_old$Day==210)] <- "memory"

which(is.na(cd4_counts_old$timepnt) == T)

cd4_counts_old$timepnt <- as.factor(cd4_counts_old$timepnt)
cd4_counts_old$timepnt <- factor(cd4_counts_old$timepnt, levels=c("baseline", "peak", "memory"))
levels(cd4_counts_old$timepnt)


#---------------------
# Add Dose
#---------------------
# Need to add in dose and schedule information
cd4_counts_old <- data.frame(Dose = NA,cd4_counts_old)

cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-035"),which(cd4_counts_old$Group == 1)),1] <- 5
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-035"),which(cd4_counts_old$Group == 2)),1] <- 5
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-035"),which(cd4_counts_old$Group == 3)),1] <- 5
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-035"),which(cd4_counts_old$Group == 4)),1] <- 5

cd4_counts_old[intersect(which(cd4_counts_old$Study == "H1"),which(cd4_counts_old$Group == 1)),1] <- 15
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H1"),which(cd4_counts_old$Group == 2)),1] <- 50
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H1"),which(cd4_counts_old$Group == 3)),1] <- 15
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H1"),which(cd4_counts_old$Group == 4)),1] <- 50

cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-032"),which(cd4_counts_old$Group == 1)),1] <- 50
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-032"),which(cd4_counts_old$Group == 2)),1] <- 15
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-032"),which(cd4_counts_old$Group == 3)),1] <- 50

# check
cd4_counts_old[which(is.na(cd4_counts_old$Dose)==T),]
# convert to factor
cd4_counts_old$Dose <- as.factor(cd4_counts_old$Dose)


#---------------------
# Add Schedule
#---------------------
# Need to add in dose and schedule information
cd4_counts_old <- data.frame(Schedule = NA, cd4_counts_old)

cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-035"),which(cd4_counts_old$Group == 1)),1] <- 2
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-035"),which(cd4_counts_old$Group == 2)),1] <- 3
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-035"),which(cd4_counts_old$Group == 3)),1] <- 2
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-035"),which(cd4_counts_old$Group == 4)),1] <- 3

cd4_counts_old[intersect(which(cd4_counts_old$Study == "H1"),which(cd4_counts_old$Group == 1)),1] <- 2
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H1"),which(cd4_counts_old$Group == 2)),1] <- 2
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H1"),which(cd4_counts_old$Group == 3)),1] <- 1
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H1"),which(cd4_counts_old$Group == 4)),1] <- 1

cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-032"),which(cd4_counts_old$Group == 1)),1] <- 3
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-032"),which(cd4_counts_old$Group == 2)),1] <- 3
cd4_counts_old[intersect(which(cd4_counts_old$Study == "H56-032"),which(cd4_counts_old$Group == 3)),1] <- 3

# check 
cd4_counts_old[which(is.na(cd4_counts_old$Schedule)==T),]
# convert to factor
cd4_counts_old$Schedule <- as.factor(cd4_counts_old$Schedule)
head(cd4_counts_old)


#------------------------
# Plot % of Responders
#------------------------

# first make dataframe
(df_ag85b <- data.frame(timpnt = rep(c(rep("baseline", 3), rep("peak", 3), rep("memory", 3)), 6), 
                 admin = rep(c(1,2,3), 3*6), 
                 QFT = rep(c(rep("QFT neg", 9), rep("QFT pos", 9)), 3), 
                 Conc = c(rep(5, 9*2), rep(15, 9*2), rep(50, 9*2)), 
                 Tot = NA, 
                 Responders= NA))
df_esat6 <- df_ag85b

# Calculate Tot and # Responders:
for(i in 1:nrow(df_ag85b)){
  df_ag85b$Tot[i] <- nrow(cd4_counts_old[cd4_counts_old[,"Stim"]== "Ag85B" &
                                           cd4_counts_old[,"timepnt"]==df_ag85b$timpnt[i] &
                                           cd4_counts_old[,"Schedule"]==df_ag85b$admin[i] &
                                           cd4_counts_old[,"QFT"]==df_ag85b$QFT[i] &
                                           cd4_counts_old[,"Dose"]==df_ag85b$Conc[i],])
  df_ag85b$Responders[i] <- nrow(cd4_counts[cd4_counts[,"Stim"]== "Ag85B" &
                                              cd4_counts[,"timepnt"]==df_ag85b$timpnt[i] &
                                              cd4_counts[,"Schedule"]==df_ag85b$admin[i] &
                                              cd4_counts[,"QFT"]==df_ag85b$QFT[i] &
                                              cd4_counts[,"Dose"]==df_ag85b$Conc[i],])
  
}

for(i in 1:nrow(df_esat6)){
  df_esat6$Tot[i] <- nrow(cd4_counts_old[cd4_counts_old[,"Stim"]== "ESAT6" &
                                           cd4_counts_old[,"timepnt"]==df_esat6$timpnt[i] &
                                           cd4_counts_old[,"Schedule"]==df_esat6$admin[i] &
                                           cd4_counts_old[,"QFT"]==df_esat6$QFT[i] &
                                           cd4_counts_old[,"Dose"]==df_esat6$Conc[i],])
  df_esat6$Responders[i] <- nrow(cd4_counts[cd4_counts[,"Stim"]== "ESAT6" &
                                              cd4_counts[,"timepnt"]==df_esat6$timpnt[i] &
                                              cd4_counts[,"Schedule"]==df_esat6$admin[i] &
                                              cd4_counts[,"QFT"]==df_esat6$QFT[i] &
                                              cd4_counts[,"Dose"]==df_esat6$Conc[i],])
  
}
df_ag85b
df_esat6




# Remove zeros:
df_ag85b <- df_ag85b[-which(df_ag85b$Tot == 0),]
df_esat6 <- df_esat6[-which(df_esat6$Tot == 0),]


# Calculate % of responders:
df_ag85b$percent <- (df_ag85b$Responders/df_ag85b$Tot)*100
df_esat6$percent <- (df_esat6$Responders/df_esat6$Tot)*100


# Convert to factors:
df_ag85b$timpnt <- as.factor(df_ag85b$timpnt)
df_ag85b$timpnt <- factor(df_ag85b$timpnt , levels=c("baseline", "peak", "memory"))
df_ag85b$admin <- as.factor(df_ag85b$admin)
df_ag85b$QFT <- as.factor(df_ag85b$QFT)
df_ag85b$Conc <- as.factor(df_ag85b$Conc)

df_esat6$timpnt <- as.factor(df_esat6$timpnt)
df_esat6$timpnt <- factor(df_esat6$timpnt , levels=c("baseline", "peak", "memory"))
df_esat6$admin <- as.factor(df_esat6$admin)
df_esat6$QFT <- as.factor(df_esat6$QFT)
df_esat6$Conc <- as.factor(df_esat6$Conc)


#------------------------------------------------------------------------------------------ Ag85B Plot:
ggplot(df_ag85b[df_ag85b[,"QFT"]=="QFT neg" & df_ag85b[,"Conc"]==5,], aes(timpnt, percent, fill=admin)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  ylab("Response Rate (%)") +
  ylim(0,100)+
  xlab("")+
  geom_text(aes(y = percent, label = paste(round(percent), "%", sep = "")), size = 10, position = position_dodge(0.9), vjust=1.5, color="white") +
  scale_fill_manual(values=c("#FF0000","#660000"), name = "Admin.") 

ggplot(df_ag85b[df_ag85b[,"QFT"]=="QFT pos" & df_ag85b[,"Conc"]==5,], aes(timpnt, percent, fill=admin)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  ylab("Response Rate (%)") +
  ylim(0,100)+
  xlab("")+
  geom_text(aes(y = percent, label = paste(round(percent), "%", sep = "")), size = 10, position = position_dodge(0.9), vjust=1.5, color="white") +
  scale_fill_manual(values=c("blue", "#000066"), name = "Admin.") 


ggplot(df_ag85b[df_ag85b[,"QFT"]=="QFT neg" & df_ag85b[,"Conc"]==15,], aes(timpnt, percent, fill=admin)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  ylab("Response Rate (%)") +
  ylim(0,100)+
  xlab("")+
  geom_text(aes(y = percent, label = paste(round(percent), "%", sep = "")), size = 10, position = position_dodge(0.9), vjust=1.5, color="white") +
  scale_fill_manual(values=c("#FF9999","#FF0000"), name = "Admin.") 

ggplot(df_ag85b[df_ag85b[,"QFT"]=="QFT pos" & df_ag85b[,"Conc"]==15,], aes(timpnt, percent, fill=admin)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  ylab("Response Rate (%)") +
  ylim(0,100)+
  xlab("")+
  geom_text(aes(y = percent, label = paste(round(percent),"%", sep = "")), size = 8, position = position_dodge(0.9), vjust=1.5, color="white") +
  scale_fill_manual(values=c("#99CCFF", "blue", "#000066"), name = "Admin.") 



ggplot(df_ag85b[df_ag85b[,"QFT"]=="QFT neg" & df_ag85b[,"Conc"]==50,], aes(timpnt, percent, fill=admin)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  ylab("Response Rate (%)") +
  ylim(0,100)+
  xlab("")+
  geom_text(aes(y = percent, label = paste(round(percent),"%", sep = "")), size = 10, position = position_dodge(0.9), vjust=1.5, color="white") +
  scale_fill_manual(values=c("#FF9999","#FF0000","#660000"), name = "Admin.") 

ggplot(df_ag85b[df_ag85b[,"QFT"]=="QFT pos" & df_ag85b[,"Conc"]==50,], aes(timpnt, percent, fill=admin)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  ylab("Response Rate (%)") +
  ylim(0,100)+
  xlab("")+
  geom_text(aes(y = percent, label = paste(round(percent),"%", sep = "")), size = 10, position = position_dodge(0.9), vjust=1.5, color="white") +
  scale_fill_manual(values=c("#99CCFF", "blue", "#000066"), name = "Admin.") 





#------------------------------------------------------------------------------------------ ESAT6 Plot:
ggplot(df_esat6[df_esat6[,"QFT"]=="QFT neg" & df_esat6[,"Conc"]==5,], aes(timpnt, percent, fill=admin)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  ylab("Response Rate (%)") +
  ylim(0,100)+
  xlab("")+
  geom_text(aes(y = percent, label = paste(round(percent), "%", sep = "")), size = 10, position = position_dodge(0.9), vjust=1.5, color="white") +
  scale_fill_manual(values=c("#FF0000","#660000"), name = "Admin.") 

ggplot(df_esat6[df_esat6[,"QFT"]=="QFT pos" & df_esat6[,"Conc"]==5,], aes(timpnt, percent, fill=admin)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  ylab("Response Rate (%)") +
  ylim(0,100)+
  xlab("")+
  geom_text(aes(y = percent, label = paste(round(percent), "%", sep = "")), size = 10, position = position_dodge(0.9), vjust=1.5, color="white") +
  scale_fill_manual(values=c("blue", "#000066"), name = "Admin.") 


ggplot(df_esat6[df_esat6[,"QFT"]=="QFT neg" & df_esat6[,"Conc"]==15,], aes(timpnt, percent, fill=admin)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  ylab("Response Rate (%)") +
  ylim(0,100)+
  xlab("")+
  geom_text(aes(y = percent, label = paste(round(percent), "%", sep = "")), size = 10, position = position_dodge(0.9), vjust=1.5, color="white") +
  scale_fill_manual(values=c("#FF9999","#FF0000"), name = "Admin.") 

ggplot(df_esat6[df_esat6[,"QFT"]=="QFT pos" & df_esat6[,"Conc"]==15,], aes(timpnt, percent, fill=admin)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  ylab("Response Rate (%)") +
  ylim(0,100)+
  xlab("")+
  geom_text(aes(y = percent, label = paste(round(percent),"%", sep = "")), size = 8, position = position_dodge(0.9), vjust=1.5, color="white") +
  scale_fill_manual(values=c("#99CCFF", "blue", "#000066"), name = "Admin.") 



ggplot(df_esat6[df_esat6[,"QFT"]=="QFT neg" & df_esat6[,"Conc"]==50,], aes(timpnt, percent, fill=admin)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  ylab("Response Rate (%)") +
  ylim(0,100)+
  xlab("")+
  geom_text(aes(y = percent, label = paste(round(percent),"%", sep = "")), size = 10, position = position_dodge(0.9), vjust=1.5, color="white") +
  scale_fill_manual(values=c("#FF9999","#FF0000","#660000"), name = "Admin.") 

ggplot(df_esat6[df_esat6[,"QFT"]=="QFT pos" & df_esat6[,"Conc"]==50,], aes(timpnt, percent, fill=admin)) +
  geom_bar(stat='identity', position='dodge') +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  ylab("Response Rate (%)") +
  ylim(0,100)+
  xlab("")+
  geom_text(aes(y = percent, label = paste(round(percent),"%", sep = "")), size = 10, position = position_dodge(0.9), vjust=1.5, color="white") +
  scale_fill_manual(values=c("#99CCFF", "blue", "#000066"), name = "Admin.") 











