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
library(ggnewscale)


#---------------------
# Import Data
#---------------------

cd4_counts <- read.csv("Filtered_cd4_counts_dataset_1.csv", check.names = F)
str(cd4_counts)

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


#-------------------------------
# Calculate TRF using old data
#-------------------------------

# Only keep counts for cytokine positive cells
cd4_trf <- cd4_counts_old[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd4_counts_old), fixed = TRUE)==T)]
str(cd4_trf)

# Convert to Frequencies
cd4_trf[,11:70] <- (cd4_trf[,11:70]/cd4_trf[,10])*100

# Background Subtract
cd4_bstrf <- data.frame()
for(i in 1:length(unique(cd4_counts_old$PID))){
  pid <- unique(cd4_counts_old$PID)[i]
  dat <- cd4_trf[which(cd4_trf$PID == pid),]
  days <- unique(dat$Day)
  for (j in 1:length(days)) {
    day <- days[j]
    df  <- dat[which(dat$Day == day),]
    df[-which(df$Stim == "UNS"),11:70]  <- df[-which(df$Stim == "UNS"),11:70] - df[rep(which(df$Stim == "UNS"), (length(df$Stim)-1)),11:70]
    cd4_bstrf <- rbind(cd4_bstrf,df)
  }
}
summary(cd4_bstrf)
which(is.na(cd4_bstrf)==T)


# Deal with Negative Frequencies:
for (i in 11:ncol(cd4_bstrf)) {
  # select one column of marker combinations
  column <- cd4_bstrf[,i]
  # extract negative values
  negval <- column[which(column<0)]
  
  # only proceed if there were negative values found
  if(length(negval)>0){
    # find upper and lower bounds of negative value
    ci  <- CI(negval,ci = 0.8)
    # calculate threshold
    cut <- abs(ci["upper"]-ci["lower"])
    # find indices of the values below the cut-off
    ind <- which(column<cut)
    # set all values less than the cut-off to 0 (including negative numbers)
    cd4_bstrf[ind,i] <- 0
  }
  
}
# Check:
sapply(sign(cd4_bstrf[,-(1:10)]), table)

# Sum across rows
cd4_bstrf <- cd4_bstrf %>%
  mutate(Total = select(., 11:70) %>% rowSums(na.rm = TRUE))
str(cd4_bstrf)

# Remove unnecessary rows
cd4_bstrf <- cd4_bstrf[,c(1,2,3,5,7,9,71)] 
head(cd4_bstrf)

# Reshape Data
cd4_bstrf <- reshape(cd4_bstrf, idvar = c("Schedule", "Dose",  "timepnt", "QFT", "PID"), timevar = "Stim", direction = "wide")
head(cd4_bstrf)
length(unique(cd4_bstrf$PID))
length(unique(cd4_counts_old$PID)) 



#--------------------------------------
# Create Responders Indicator Variable
#--------------------------------------

# first need separate datasets for Ag85B and ESAT6:
trf_ag85B <- cd4_bstrf[,-c(6,8,9)]
trf_esat6 <- cd4_bstrf[,-c(6,7,9)]

# Remove NAs:
trf_ag85B <- na.omit(trf_ag85B)
trf_esat6 <- na.omit(trf_esat6)

# list responders:
resp_ag85B <- cd4_counts[cd4_counts[,"Stim"]=="Ag85B",c(3,7)]
resp_esat6 <- cd4_counts[cd4_counts[,"Stim"]=="ESAT6",c(3,7)]

# Add indicator variable:
trf_ag85B$responder <- 0.  # 0 for no, 1 for yes
trf_esat6$responder <- 0

# Populate responder indicator: 
for(i in 1:nrow(resp_ag85B)) {
  id  <- resp_ag85B[i,"PID"]
  day <- resp_ag85B[i,"timepnt"]
  trf_ag85B[trf_ag85B[,"PID"]==id & trf_ag85B[,"timepnt"]==day, "responder"] <- 1
}
sum(trf_ag85B$responder) == nrow(resp_ag85B) # T, thus we have the right number of responders

for(i in 1:nrow(resp_esat6)) {
  id  <- resp_esat6[i,"PID"]
  day <- resp_esat6[i,"timepnt"]
  trf_esat6[trf_esat6[,"PID"]==id & trf_esat6[,"timepnt"]==day, "responder"] <- 1
}
sum(trf_esat6$responder) == nrow(resp_esat6) # T, thus we have the right number of responders




#--------------------------------------
# Plot Ag85B
#--------------------------------------

# Set as factor:
trf_ag85B$responder <- as.factor(trf_ag85B$responder)
trf_esat6$responder <- as.factor(trf_esat6$responder)

# Change "Schedule" to "Administrations":
names(trf_ag85B)[1] <- "Admin."
names(trf_esat6)[1] <- "Admin."

#----------------------------------------------------------- QFT- with 5mg:

#geom_point(aes(color = responder), position=position_jitterdodge(),alpha=0.5)+
#geom_point(aes(color = responder),shape = 21,position = position_jitter(seed = 1, width = .05)) +
#geom_jitter(shape=16, position=position_jitter(0.2), aes(color=responder))+
 

ggplot(trf_ag85B[trf_ag85B[,"QFT"]=="QFT neg" & trf_ag85B[,"Dose"]==5,], aes(x = timepnt, y = Total.Ag85B, fill=Admin.)) +
  geom_boxplot(aes(color = Admin.), fill = NA, outlier.shape = NA, lwd=2) +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  scale_color_manual(values = c("#FF0000","#660000")) +
  new_scale_color() +
  geom_point(aes(color = responder, group=Admin.), position=position_jitterdodge(), size=3, alpha=0.6)+    # position_dodge(0.9)
  scale_color_manual(values = c("black", "#FF00FF"), labels = c("NR", "R"), name="")+
  xlab("")+
  ylim(0,2)+
  ylab("TRF (%)") 



#----------------------------------------------------------- QFT+ with 5mg:

ggplot(trf_ag85B[trf_ag85B[,"QFT"]=="QFT pos" & trf_ag85B[,"Dose"]==5,], aes(x = timepnt, y = Total.Ag85B, fill=Admin.)) +
  geom_boxplot(aes(color = Admin.), fill = NA, outlier.shape = NA, lwd=2) +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  scale_color_manual(values = c("blue", "#000066")) +
  new_scale_color() +
  geom_point(aes(color = responder, group=Admin.), position=position_jitterdodge(), size=3, alpha=0.6)+    # position_dodge(0.9)
  scale_color_manual(values = c("black", "#FF00FF"), labels = c("NR", "R"), name="")+
  xlab("")+
  ylim(0,2)+
  ylab("TRF (%)") 




#----------------------------------------------------------- QFT- with 15mg:

ggplot(trf_ag85B[trf_ag85B[,"QFT"]=="QFT neg" & trf_ag85B[,"Dose"]==15,], aes(x = timepnt, y = Total.Ag85B, fill=Admin.)) +
  geom_boxplot(aes(color = Admin.), fill = NA, outlier.shape = NA, lwd=2) +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  scale_color_manual(values =c("#FF9999","#FF0000","#660000")) +
  new_scale_color() +
  geom_point(aes(color = responder, group=Admin.), position=position_jitterdodge(), size=3, alpha=0.6)+    # position_dodge(0.9)
  scale_color_manual(values = c("black", "#FF00FF"), labels =  c("NR", "R"), name="")+
  xlab("")+
  ylim(0,2)+
  ylab("TRF (%)") 



#----------------------------------------------------------- QFT+ with 15mg:

ggplot(trf_ag85B[trf_ag85B[,"QFT"]=="QFT pos" & trf_ag85B[,"Dose"]==15,], aes(x = timepnt, y = Total.Ag85B, fill=Admin.)) +
  geom_boxplot(aes(color = Admin.), fill = NA, outlier.shape = NA, lwd=2) +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  scale_color_manual(values = c("#99CCFF", "blue", "#000066")) +
  new_scale_color() +
  geom_point(aes(color = responder, group=Admin.), position=position_jitterdodge(), size=3, alpha=0.6)+    # position_dodge(0.9)
  scale_color_manual(values = c("black", "#FF00FF"), labels =  c("NR", "R"), name="")+
  xlab("")+
  ylim(0,2)+
  ylab("TRF (%)") 




#----------------------------------------------------------- QFT- with 50mg:

ggplot(trf_ag85B[trf_ag85B[,"QFT"]=="QFT neg" & trf_ag85B[,"Dose"]==50,], aes(x = timepnt, y = Total.Ag85B, fill=Admin.)) +
  geom_boxplot(aes(color = Admin.), fill = NA, outlier.shape = NA, lwd=2) +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  scale_color_manual(values =c("#FF9999","#FF0000","#660000")) +
  new_scale_color() +
  geom_point(aes(color = responder, group=Admin.), position=position_jitterdodge(), size=3, alpha=0.6)+    # position_dodge(0.9)
  scale_color_manual(values = c("black", "#FF00FF"), labels =  c("NR", "R"), name="")+
  xlab("")+
  ylim(0,2)+
  ylab("TRF (%)") 



#----------------------------------------------------------- QFT+ with 50mg:

ggplot(trf_ag85B[trf_ag85B[,"QFT"]=="QFT pos" & trf_ag85B[,"Dose"]==50,], aes(x = timepnt, y = Total.Ag85B, fill=Admin.)) +
  geom_boxplot(aes(color = Admin.), fill = NA, outlier.shape = NA, lwd=2) +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  scale_color_manual(values = c("#99CCFF", "blue", "#000066")) +
  new_scale_color() +
  geom_point(aes(color = responder, group=Admin.), position=position_jitterdodge(), size=3, alpha=0.6)+    # position_dodge(0.9)
  scale_color_manual(values = c("black", "#FF00FF"), labels =  c("NR", "R"), name="")+
  xlab("")+
  ylim(0,2)+
  ylab("TRF (%)") 






#--------------------------------------
# Plot ESAT6
#--------------------------------------

#----------------------------------------------------------- QFT- with 5mg:

ggplot(trf_esat6[trf_esat6[,"QFT"]=="QFT neg" & trf_esat6[,"Dose"]==5,], aes(x = timepnt, y = Total.ESAT6, fill=Admin.)) +
  geom_boxplot(aes(color = Admin.), fill = NA, outlier.shape = NA, lwd=2) +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  scale_color_manual(values = c("#FF0000","#660000")) +
  new_scale_color() +
  geom_point(aes(color = responder, group=Admin.), position=position_jitterdodge(), size=3, alpha=0.6)+    # position_dodge(0.9)
  scale_color_manual(values = c("black", "#FF00FF"), labels =  c("NR", "R"), name="")+
  xlab("")+
  ylim(0,1.85)+
  ylab("TRF (%)") 



#----------------------------------------------------------- QFT+ with 5mg:

ggplot(trf_esat6[trf_esat6[,"QFT"]=="QFT pos" & trf_esat6[,"Dose"]==5,], aes(x = timepnt, y = Total.ESAT6, fill=Admin.)) +
  geom_boxplot(aes(color = Admin.), fill = NA, outlier.shape = NA, lwd=2) +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  scale_color_manual(values = c("blue", "#000066")) +
  new_scale_color() +
  geom_point(aes(color = responder, group=Admin.), position=position_jitterdodge(), size=3, alpha=0.6)+    # position_dodge(0.9)
  scale_color_manual(values = c("black", "#FF00FF"), labels =  c("NR", "R"), name="")+
  xlab("")+
  ylim(0,1.85)+
  ylab("TRF (%)") 




#----------------------------------------------------------- QFT- with 15mg:

ggplot(trf_esat6[trf_esat6[,"QFT"]=="QFT neg" & trf_esat6[,"Dose"]==15,], aes(x = timepnt, y = Total.ESAT6, fill=Admin.)) +
  geom_boxplot(aes(color = Admin.), fill = NA, outlier.shape = NA, lwd=2) +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  scale_color_manual(values =c("#FF9999","#FF0000","#660000")) +
  new_scale_color() +
  geom_point(aes(color = responder, group=Admin.), position=position_jitterdodge(), size=3, alpha=0.6)+    # position_dodge(0.9)
  scale_color_manual(values = c("black", "#FF00FF"), labels =  c("NR", "R"), name="")+
  xlab("")+
  ylim(0,1.85)+
  ylab("TRF (%)") 



#----------------------------------------------------------- QFT+ with 15mg:

ggplot(trf_esat6[trf_esat6[,"QFT"]=="QFT pos" & trf_esat6[,"Dose"]==15,], aes(x = timepnt, y = Total.ESAT6, fill=Admin.)) +
  geom_boxplot(aes(color = Admin.), fill = NA, outlier.shape = NA, lwd=2) +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  scale_color_manual(values = c("#99CCFF", "blue", "#000066")) +
  new_scale_color() +
  geom_point(aes(color = responder, group=Admin.), position=position_jitterdodge(), size=3, alpha=0.6)+    # position_dodge(0.9)
  scale_color_manual(values = c("black", "#FF00FF"), labels =  c("NR", "R"), name="")+
  xlab("")+
  ylim(0,1.85)+
  ylab("TRF (%)") 




#----------------------------------------------------------- QFT- with 50mg:

ggplot(trf_esat6[trf_esat6[,"QFT"]=="QFT neg" & trf_esat6[,"Dose"]==50,], aes(x = timepnt, y = Total.ESAT6, fill=Admin.)) +
  geom_boxplot(aes(color = Admin.), fill = NA, outlier.shape = NA, lwd=2) +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  scale_color_manual(values =c("#FF9999","#FF0000","#660000")) +
  new_scale_color() +
  geom_point(aes(color = responder, group=Admin.), position=position_jitterdodge(), size=3, alpha=0.6)+    # position_dodge(0.9)
  scale_color_manual(values = c("black", "#FF00FF"), labels =  c("NR", "R"), name="")+
  xlab("")+
  ylim(0,1.85)+
  ylab("TRF (%)") 



#----------------------------------------------------------- QFT+ with 50mg:

ggplot(trf_esat6[trf_esat6[,"QFT"]=="QFT pos" & trf_esat6[,"Dose"]==50,], aes(x = timepnt, y = Total.ESAT6, fill=Admin.)) +
  geom_boxplot(aes(color = Admin.), fill = NA, outlier.shape = NA, lwd=2) +
  theme_minimal()+ 
  theme(text = element_text(size = 30))+
  scale_color_manual(values = c("#99CCFF", "blue", "#000066")) +
  new_scale_color() +
  geom_point(aes(color = responder, group=Admin.), position=position_jitterdodge(), size=3, alpha=0.6)+    # position_dodge(0.9)
  scale_color_manual(values = c("black", "#FF00FF"), labels =  c("NR", "R"), name="")+
  xlab("")+
  ylim(0,1.85)+
  ylab("TRF (%)") 



