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
library(Rmisc) # to get CI

# Suppress Scientific Notation
options(scipen=999)

#---------------------
# Import Data
#---------------------

cd4_counts <- read.csv("clean_cd4_counts_dataset_1.csv", check.names = F)
str(cd4_counts)

#---------------------
# Calculate TRF
#---------------------

# Only keep counts for cells expressing at least one of INFg, IL2,TNF or L17
#cd4_trf <- cd4_counts[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd4_counts), fixed = TRUE)==T)]
#str(cd4_trf)

# Convert to Frequencies
freq <- cd4_counts
freq[,8:71] <- (cd4_counts[,8:71]/cd4_counts[,7])*100

# Background Subtract
cd4_bs <- data.frame()
for(i in 1:length(unique(cd4_counts$PID))){
  pid <- unique(cd4_counts$PID)[i]
  dat <- freq[which(freq$PID == pid),]
  days <- unique(dat$Day)
  for (j in 1:length(days)) {
    day <- days[j]
    df  <- dat[which(dat$Day == day),]
    df[-which(df$Stim == "UNS"),8:71]  <- df[-which(df$Stim == "UNS"),8:71] - df[rep(which(df$Stim == "UNS"), (length(df$Stim)-1)),8:71]
    cd4_bs <- rbind(cd4_bs,df)
  }
}
summary(cd4_bs)
which(is.na(cd4_bs)==T)


# Deal with Negative Frequencies:
for (i in 8:ncol(cd4_bs)) {
  # select one column of marker combinations
  column <- cd4_bs[,i]
  # extract negative values
  negval <- column[which(column<0)]
  
  # only proceed if there were negative values found
  if(length(negval)>1){ # have to have more than one value to get a CI
    # find upper and lower bounds of negative value
    ci  <- CI(negval,ci = 0.8)
    # calculate threshold
    cut <- abs(ci["upper"]-ci["lower"])
    # find indices of the values below the cut-off
    ind <- which(column<cut)
    # set all values less than the cut-off to 0 (including negative numbers)
    cd4_bs[ind,i] <- 0 }
  
  if(length(negval)==1){ # If only 1 value, abandon CI and set neg value to 0
    # find index of the value below the 0
    ind <- which(column<0)
    # set value to 0 
    cd4_bs[ind,i] <- 0 }
  
}
# Check:
sapply(sign(cd4_bs[,-(1:8)]), table)





#---------------------
# Add Time Category
#---------------------
cd4_bs <- data.frame(timepnt = NA, cd4_bs)

cd4_bs$timepnt[which(cd4_bs$Day==0)] <- "baseline"
cd4_bs$timepnt[which(cd4_bs$Day==14)] <- "peak"
cd4_bs$timepnt[which(cd4_bs$Day==70)] <- "peak"
cd4_bs$timepnt[which(cd4_bs$Day==126)] <- "peak"
cd4_bs$timepnt[which(cd4_bs$Day==224)] <- "memory"
cd4_bs$timepnt[which(cd4_bs$Day==292)] <- "memory"
cd4_bs$timepnt[which(cd4_bs$Day==210)] <- "memory"

which(is.na(cd4_bs$timepnt) == T)

cd4_bs$timepnt <- as.factor(cd4_bs$timepnt)
cd4_bs$timepnt <- factor(cd4_bs$timepnt, levels=c("baseline", "peak", "memory"))
levels(cd4_bs$timepnt)


#---------------------
# Add Dose
#---------------------
# Need to add in dose and schedule information
cd4_bs <- data.frame(Dose = NA, cd4_bs)

cd4_bs[intersect(which(cd4_bs$Study == "H56-035"),which(cd4_bs$Group == 1)),1] <- 5
cd4_bs[intersect(which(cd4_bs$Study == "H56-035"),which(cd4_bs$Group == 2)),1] <- 5
cd4_bs[intersect(which(cd4_bs$Study == "H56-035"),which(cd4_bs$Group == 3)),1] <- 5
cd4_bs[intersect(which(cd4_bs$Study == "H56-035"),which(cd4_bs$Group == 4)),1] <- 5

cd4_bs[intersect(which(cd4_bs$Study == "H1"),which(cd4_bs$Group == 1)),1] <- 15
cd4_bs[intersect(which(cd4_bs$Study == "H1"),which(cd4_bs$Group == 2)),1] <- 50
cd4_bs[intersect(which(cd4_bs$Study == "H1"),which(cd4_bs$Group == 3)),1] <- 15
cd4_bs[intersect(which(cd4_bs$Study == "H1"),which(cd4_bs$Group == 4)),1] <- 50

cd4_bs[intersect(which(cd4_bs$Study == "H56-032"),which(cd4_bs$Group == 1)),1] <- 50
cd4_bs[intersect(which(cd4_bs$Study == "H56-032"),which(cd4_bs$Group == 2)),1] <- 15
cd4_bs[intersect(which(cd4_bs$Study == "H56-032"),which(cd4_bs$Group == 3)),1] <- 50

# check
cd4_bs[which(is.na(cd4_bs$Dose)==T),]
# convert to factor
cd4_bs$Dose <- as.factor(cd4_bs$Dose)


#---------------------
# Add Schedule
#---------------------
# Need to add in dose and schedule information
cd4_bs <- data.frame(Schedule = NA, cd4_bs)

cd4_bs[intersect(which(cd4_bs$Study == "H56-035"),which(cd4_bs$Group == 1)),1] <- 2
cd4_bs[intersect(which(cd4_bs$Study == "H56-035"),which(cd4_bs$Group == 2)),1] <- 3
cd4_bs[intersect(which(cd4_bs$Study == "H56-035"),which(cd4_bs$Group == 3)),1] <- 2
cd4_bs[intersect(which(cd4_bs$Study == "H56-035"),which(cd4_bs$Group == 4)),1] <- 3

cd4_bs[intersect(which(cd4_bs$Study == "H1"),which(cd4_bs$Group == 1)),1] <- 2
cd4_bs[intersect(which(cd4_bs$Study == "H1"),which(cd4_bs$Group == 2)),1] <- 2
cd4_bs[intersect(which(cd4_bs$Study == "H1"),which(cd4_bs$Group == 3)),1] <- 1
cd4_bs[intersect(which(cd4_bs$Study == "H1"),which(cd4_bs$Group == 4)),1] <- 1

cd4_bs[intersect(which(cd4_bs$Study == "H56-032"),which(cd4_bs$Group == 1)),1] <- 3
cd4_bs[intersect(which(cd4_bs$Study == "H56-032"),which(cd4_bs$Group == 2)),1] <- 3
cd4_bs[intersect(which(cd4_bs$Study == "H56-032"),which(cd4_bs$Group == 3)),1] <- 3

# check 
cd4_bs[which(is.na(cd4_bs$Schedule)==T),]
# convert to factor
cd4_bs$Schedule <- as.factor(cd4_bs$Schedule)
head(cd4_bs)

str(cd4_bs)




#---------------------
# Keep only Il2+
#---------------------

il2 <- cd4_bs[,-which(grepl("L2c..", names(cd4_bs), fixed = TRUE)==T)]
str(il2)

il2 <- il2[il2[,"Schedule"]==2 & il2[, "Dose"]==5,]

# Calculate frequency of il2 expressing cells 
il2$Freq <- rowSums(il2[,11:42]) # bs frequency of CD4


#----------------------------------------------------------------- Peak Ag85B:
dftest      <- il2[il2[,"timepnt"]=="peak" & il2[,"Stim"]=="Ag85B",]
# Wilcoxon Test:
stat.test   <- dftest %>%
  wilcox_test(Freq ~ QFT, p.adjust.method = "BH") 
stat.test 
stat.test$p <- round(stat.test$p, 3)
stat.test   <- stat.test %>%
  add_xy_position(x = "QFT", dodge = 0.8)
# Log Transform Response:
dataf      <- dftest
dataf$Freq <- log(dataf$Freq+0.00001)
ggplot(dataf, aes(y=Freq, x=QFT)) +
  geom_boxplot(alpha=1, color=c("#660000", "#000066"), size =1, aes(fill = QFT), show.legend = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 28)) +
  ylab("Logged % of CD4 cells expressing IL-2") +
  xlab("")+
  ylim(-11.6,1)+
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  scale_fill_manual(values=c("#FF0000", "blue"))+
  stat_pvalue_manual(stat.test, 
                     label = "p", 
                     tip.length = 0.02,
                     size = 10, 
                     bracket.size = 0.8) 

#----------------------------------------------------------------- Peak ESAT6:
dftest <- il2[il2[,"timepnt"]=="peak" & il2[,"Stim"]=="ESAT6",]
# Wilcoxon Test:
stat.test <- dftest %>%
  wilcox_test(Freq ~ QFT, p.adjust.method = "BH") 
stat.test 
stat.test$p <- round(stat.test$p, 3)
stat.test <- stat.test %>%
  add_xy_position(x = "QFT", dodge = 0.8)
# Log Transform Response:
dataf <- dftest
dataf$Freq <- log(dataf$Freq+0.00001)
ggplot(dataf, aes(y=Freq, x=QFT)) +
  geom_boxplot(alpha=1, color=c("#660000", "#000066"), size =1, aes(fill = QFT), show.legend = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 28)) +
  ylab("Logged % of CD4 cells expressing IL-2") +
  xlab("")+
  ylim(-6.4,1)+
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  scale_fill_manual(values=c("#FF0000", "blue"))+
  stat_pvalue_manual(stat.test, 
                     label = "p", 
                     tip.length = 0.02,
                     bracket.nudge.y = 0.2, 
                     size = 10, 
                     bracket.size = 0.8) 


#----------------------------------------------------------------- Ag85B Memory:
dftest <- il2[il2[,"timepnt"]=="memory"& il2[,"Stim"]=="Ag85B",]
# Wilcoxon Test:
stat.test <- dftest %>%
  wilcox_test(Freq ~ QFT, p.adjust.method = "BH") 
stat.test 
stat.test <- stat.test %>%
  add_xy_position(x = "QFT", dodge = 0.8)
# Log Transform Response:
dataf <- dftest
dataf$Freq <- log(dataf$Freq+0.00001)
ggplot(dataf, aes(y=Freq, x=QFT)) +
  geom_boxplot(alpha=1, color=c("#660000", "#000066"), size =1, aes(fill = QFT), show.legend = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 28)) +
  ylab("Logged % of CD4 cells expressing IL-2") +
  xlab("")+
  ylim(-11.6,1)+
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  scale_fill_manual(values=c("#FF0000", "blue"))+
  stat_pvalue_manual(stat.test, 
                     label = "p", 
                     tip.length = 0.02,
                     bracket.nudge.y = 0.2, 
                     size = 10, 
                     bracket.size = 0.8) 


#----------------------------------------------------------------- ESAT6 Memory:
dftest <- il2[il2[,"timepnt"]=="memory"& il2[,"Stim"]=="ESAT6",]
# Wilcoxon Test:
stat.test <- dftest %>%
  wilcox_test(Freq ~ QFT, p.adjust.method = "BH") 
stat.test 
stat.test$p <- round(stat.test$p, 6)
stat.test <- stat.test %>%
  add_xy_position(x = "QFT", dodge = 0.8)
# Log Transform Response:
dataf <- dftest
dataf$Freq <- log(dataf$Freq+0.00001)
ggplot(dataf, aes(y=Freq, x=QFT)) +
  geom_boxplot(alpha=1, color=c("#660000", "#000066"), size =1, aes(fill = QFT), show.legend = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 28)) +
  ylab("Logged % of CD4 cells expressing IL-2") +
  xlab("")+
  ylim(-6.4,1)+
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  scale_fill_manual(values=c("#FF0000", "blue"))+
  stat_pvalue_manual(stat.test, 
                     label = "p", 
                     tip.length = 0.02,
                     bracket.nudge.y = 0.2, 
                     size = 10, 
                     bracket.size = 0.8) 



#---------------------
# Keep only Il17+
#---------------------

il17 <- cd4_bs[,-which(grepl("L17c..", names(cd4_bs), fixed = TRUE)==T)]
str(il17)

il17 <- il17[il17[,"Schedule"]==2 & il17[, "Dose"]==5,]

# Calculate frequency of il17 expressing cells 
il17$Freq <- rowSums(il17[,11:42]) # bs frequency of CD4


#----------------------------------------------------------------- Peak Ag85B:
dftest <- il17[il17[,"timepnt"]=="peak"& il17[,"Stim"]=="Ag85B",]
# Wilcoxon Test:
stat.test <- dftest %>%
  wilcox_test(Freq ~ QFT, p.adjust.method = "BH") 
stat.test 
stat.test <- stat.test %>%
  add_xy_position(x = "QFT", dodge = 0.8)
stat.test$y.position <- -2.1
# Log Transform Response:
dataf <- dftest
dataf$Freq <- log(dataf$Freq+0.00001)
ggplot(dataf, aes(y=Freq, x=QFT)) +
  geom_boxplot(alpha=1, color=c("#660000", "#000066"), size =1, aes(fill = QFT), show.legend = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("Logged % of CD4 cells expressing IL-17") +
  xlab("")+
  ylim(-11.6,-2)+
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  scale_fill_manual(values=c("#FF0000", "blue"))+
  stat_pvalue_manual(stat.test, 
                     label = "p", 
                     tip.length = 0.02,
                     #bracket.nudge.y = 0.01, 
                     size = 10, 
                     bracket.size = 0.8) 


#----------------------------------------------------------------- Peak ESAT6:
dftest <- il17[il17[,"timepnt"]=="peak"& il17[,"Stim"]=="ESAT6",]
# Wilcoxon Test:
stat.test <- dftest %>%
  wilcox_test(Freq ~ QFT, p.adjust.method = "BH") 
stat.test 
stat.test$p <- round(stat.test$p, 3)
stat.test <- stat.test %>%
  add_xy_position(x = "QFT", dodge = 0.8)
stat.test$y.position <- -2.8
# Log Transform Response:
dataf <- dftest
dataf$Freq <- log(dataf$Freq+0.00001)
ggplot(dataf, aes(y=Freq, x=QFT)) +
  geom_boxplot(alpha=1, color=c("#660000", "#000066"), size =1, aes(fill = QFT), show.legend = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 28)) +
  ylab("Logged % of CD4 cells expressing IL-17") +
  xlab("")+
  ylim(-11.6,-2.7)+
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  scale_fill_manual(values=c("#FF0000", "blue"))+
  stat_pvalue_manual(stat.test, 
                     label = "p", 
                     tip.length = 0.02,
                     #bracket.nudge.y = 0.01, 
                     size = 10, 
                     bracket.size = 0.8) 



#----------------------------------------------------------------- Memory Ag85B:
dftest <- il17[il17[,"timepnt"]=="memory" & il17[,"Stim"]=="Ag85B",]
# Wilcoxon Test:
stat.test <- dftest %>%
  wilcox_test(Freq ~ QFT, p.adjust.method = "BH") 
stat.test 
stat.test <- stat.test %>%
  add_xy_position(x = "QFT", dodge = 0.8)
# Log Transform Response:
dataf <- dftest
dataf$Freq <- log(dataf$Freq+0.00001)
stat.test$y.position <- -2.3
ggplot(dataf, aes(y=Freq, x=QFT)) +
  geom_boxplot(alpha=1, color=c("#660000", "#000066"), size =1, aes(fill = QFT), show.legend = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 28)) +
  ylab("Logged % of CD4 cells expressing IL-17") +
  xlab("")+
  ylim(-11.6,-2)+
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  scale_fill_manual(values=c("#FF0000", "blue"))+
  stat_pvalue_manual(stat.test, 
                     label = "p", 
                     tip.length = 0.02,
                     #bracket.nudge.y = 0.02, 
                     size = 10, 
                     bracket.size = 0.8)


#----------------------------------------------------------------- Memory ESAT6:
dftest <- il17[il17[,"timepnt"]=="memory" & il17[,"Stim"]=="ESAT6",]
# Wilcoxon Test:
stat.test <- dftest %>%
  wilcox_test(Freq ~ QFT, p.adjust.method = "BH") 
stat.test 
stat.test$p <- round(stat.test$p, 3)
stat.test <- stat.test %>%
  add_xy_position(x = "QFT", dodge = 0.8)
stat.test$y.position <- -3
# Log Transform Response:
dataf <- dftest
dataf$Freq <- log(dataf$Freq+0.00001)
ggplot(dataf, aes(y=Freq, x=QFT)) +
  geom_boxplot(alpha=1, color=c("#660000", "#000066"), size =1, aes(fill = QFT), show.legend = FALSE) +
  theme_minimal()+
  theme(text = element_text(size = 28)) +
  ylab("Logged % of CD4 cells expressing IL-17") +
  xlab("")+
  ylim(-11.6,-2.7)+
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  scale_fill_manual(values=c("#FF0000", "blue"))+
  stat_pvalue_manual(stat.test, 
                     label = "p", 
                     tip.length = 0.02,
                     #bracket.nudge.y = 0.02, 
                     size = 10, 
                     bracket.size = 0.8) 

