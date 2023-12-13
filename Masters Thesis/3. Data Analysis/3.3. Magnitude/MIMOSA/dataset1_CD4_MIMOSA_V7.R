# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

#---------------------
# Load Libraries
#---------------------

#install.packages("devtools") #install devtools package
library(devtools) # load it
#install_github("RGLab/MIMOSA",ref="trunk") # Need R 3.0.0, and a bunch of dependencies. The install will fail with various error messages until you install those dependencies.
library(MIMOSA) 

library(tidyr)
library(dplyr)

#---------------------
# Import Data
#---------------------

cd4_counts <- read.csv("clean_cd4_counts_dataset_1.csv", check.names = F)
str(cd4_counts)

#---------------------
# Prepare Data
#---------------------

# Remove PHA:
counts <- cd4_counts[-which(cd4_counts$Stim=="PHA"),]

# Select columns for cells expressing no cytokines:
(neg_ind <- which(grepl("Gc..L2c..L17c..TNFc", names(counts), fixed = TRUE)==T))

# Calculate no. of cells expressing no cytokines:
(total_neg <- rowSums(counts[,neg_ind]))
summary(total_neg)

# Calculate no. of cells expressing all cytokines:
(total_pos <- counts$CD4count - total_neg)
summary(total_pos)

# Make new dataframe with these variables 
dftot <- counts[,4:6]
dftot$Total_pos <- total_pos
dftot$Total_neg <- total_neg
str(dftot)
head(dftot)



#---------------------
# Create MIMOSA input
#---------------------

E <- ConstructMIMOSAExpressionSet(dftot,
                                  reference=Stim%in%"UNS",
                                  measure.columns=c("Total_neg","Total_pos"),
                                  default.cast.formula=component~PID+Day+Stim,
                                  .variables=.(PID,Day)) 
E


#-----------------------
# Conduct MIMOSA: MCMC
#-----------------------
set.seed(100)
resultmcmc <- MIMOSA(Total_neg + Total_pos ~ PID + Day|Stim,
                     data=E, method="mcmc",burn=1000,iter=5000)

# Look at those with fdr < 0.01:
lapply(resultmcmc,function(x)table(fdr(x)<0.01))
pData(resultmcmc)[fdr(resultmcmc)<0.01,]


#--------------------------------------
# Select Responders (fdr < 0.01)
#--------------------------------------
# Responders are those whose FDR < 0.01
# AND who have at least a 3-fold increase in response compared to background

#------------------------------------------------ Select responders whose FDR < 0.01:
r  <- pData(resultmcmc)[fdr(resultmcmc)<0.01,c(1,2,4)]  # extract (PID, Day) responders


#----------------------------------------------- Select those with a 3-fold increase:

# first convert to a frequency:
df <- dftot
df$freq = dftot$Total_pos/(dftot$Total_pos+dftot$Total_neg)
# create wide dataframe:
dftot_w <- reshape(df[,-c(4,5)], idvar = c("PID", "Day"), timevar = "Stim", direction = "wide")
# select responders:
r3_ag85b <- dftot_w[which(dftot_w$freq.Ag85B >= 3*dftot_w$freq.UNS),c(1,2)]
r3_esat6 <- dftot_w[which(dftot_w$freq.ESAT6 >= 3*dftot_w$freq.UNS),c(1,2)]


#-----------------------------------------------  Intersect MIMOSA responders and 3-fold responders:
r_fin <- data.frame()
# Ag85B:
for (i in 1:nrow(r3_ag85b)) {
  row <- r[r[,"Stim"]=="Ag85B" & r[,"PID"]==r3_ag85b[i,"PID"] & r[,"Day"]==r3_ag85b[i,"Day"],] 
  if(dim(row)[1]!=0){
    r_fin <- rbind(r_fin, row)
  }
}
# ESAT6:
for (i in 1:nrow(r3_esat6)) {
  row <- r[r[,"Stim"]=="ESAT6" & r[,"PID"]==r3_esat6[i,"PID"] & r[,"Day"]==r3_esat6[i,"Day"],] 
  if(dim(row)[1]!=0){
    r_fin <- rbind(r_fin, row)
  }
}


# CHECK:
setdiff(r$PID, r_fin$PID)
r3_ag85b[r3_ag85b[,"PID"]=="1007",]
dftot_w[dftot_w[,"PID"]=="1007",]


# -----------------------------------------------  Create filtered CD4 dataset:

dat_filt <- data.frame()
for (i in 1:nrow(r_fin)) {
  
  dat_filt <- rbind(dat_filt, cd4_counts[cd4_counts[,"Stim"]==r_fin[i,"Stim"] & cd4_counts[,"PID"]==r_fin[i,"PID"] & cd4_counts[,"Day"]==r_fin[i,"Day"],])
  
  # Check if this visit has been added before:
  row <- dat_filt[dat_filt[,"PID"]==r_fin[i,"PID"] & dat_filt[,"Day"]==r_fin[i,"Day"],]
  # If visit has not been added before, add UNS and PHA
  if(dim(row)[1] == 1){ # if only one row, then only the above stim has been added. If more than 1, the there should be 4
    # add in corrosponding UNS value
    dat_filt <- rbind(dat_filt, cd4_counts[cd4_counts[,"Stim"]=="UNS" & cd4_counts[,"PID"]==r_fin[i,"PID"] & cd4_counts[,"Day"]==r_fin[i,"Day"],])
    # add in corrosponding PHA value (not really necessary)
    dat_filt <- rbind(dat_filt, cd4_counts[cd4_counts[,"Stim"]=="PHA" & cd4_counts[,"PID"]==r_fin[i,"PID"] & cd4_counts[,"Day"]==r_fin[i,"Day"],])
    
  }
  
}



# CHECKS: 
dim(cd4_counts)[1] - dim(dat_filt)[1]  # lost 592/4 = 148 visits

length(unique(cd4_counts$PID)) - length(unique(dat_filt$PID)) # lost only 4 pids

# Baseline:
(length(which(cd4_counts$Day == 0)) - length(which(dat_filt$Day == 0)))/length(which(cd4_counts$Day == 0))    # lost 33.79% day 0 obs
# Peak:
(length(which(cd4_counts$Day == 14)) + length(which(cd4_counts$Day == 70)) + length(which(cd4_counts$Day == 126))
  - length(which(dat_filt$Day == 14)) - length(which(dat_filt$Day == 70)) -length(which(dat_filt$Day == 126)))/
  (length(which(cd4_counts$Day == 14)) + length(which(cd4_counts$Day == 70)) + length(which(cd4_counts$Day == 126))) # lost 10.22% peak obs
# Memory:
(length(which(cd4_counts$Day == 224)) + length(which(cd4_counts$Day == 210)) + length(which(cd4_counts$Day == 292))
  - length(which(dat_filt$Day == 224)) - length(which(dat_filt$Day == 210)) -length(which(dat_filt$Day == 292)))/
  (length(which(cd4_counts$Day == 224)) + length(which(cd4_counts$Day == 210)) + length(which(cd4_counts$Day == 292))) # lost 16.74% memory obs

dat_filt$QFT <- as.factor(dat_filt$QFT)
dat_filt$Study <- as.factor(dat_filt$Study)
dat_filt$Day <- as.factor(dat_filt$Day)
dat_filt$Stim <- as.factor(dat_filt$Stim)

cd4_counts$QFT <- as.factor(cd4_counts$QFT)
cd4_counts$Study <- as.factor(cd4_counts$Study)
cd4_counts$Day <- as.factor(cd4_counts$Day)
cd4_counts$Stim <- as.factor(cd4_counts$Stim)

summary(cd4_counts[,c(1,2,5,6)])
summary(dat_filt[,c(1,2,5,6)])

# Ag85B:
nrow(dat_filt[dat_filt[,"Stim"]=="Ag85B",]) # number of responders
nrow(cd4_counts[cd4_counts[,"Stim"]=="Ag85B",]) - nrow(dat_filt[dat_filt[,"Stim"]=="Ag85B",]) # number of non-responers
length(unique(dat_filt[dat_filt[,"Stim"]=="Ag85B","PID"])) # final number of pids


# ESAT6
nrow(dat_filt[dat_filt[,"Stim"]=="ESAT6",])  # number of responders
nrow(cd4_counts[cd4_counts[,"Stim"]=="ESAT6",]) - nrow(dat_filt[dat_filt[,"Stim"]=="ESAT6",]) # number of non-responers
length(unique(dat_filt[dat_filt[,"Stim"]=="ESAT6","PID"])) # final number of pids





#---------------------
# Add Time Category
#---------------------
dat_filt <- data.frame(timepnt = NA, dat_filt)

dat_filt$timepnt[which(dat_filt$Day==0)] <- "baseline"
dat_filt$timepnt[which(dat_filt$Day==14)] <- "peak"
dat_filt$timepnt[which(dat_filt$Day==70)] <- "peak"
dat_filt$timepnt[which(dat_filt$Day==126)] <- "peak"
dat_filt$timepnt[which(dat_filt$Day==224)] <- "memory"
dat_filt$timepnt[which(dat_filt$Day==292)] <- "memory"
dat_filt$timepnt[which(dat_filt$Day==210)] <- "memory"

which(is.na(dat_filt$timepnt) == T)

dat_filt$timepnt <- as.factor(dat_filt$timepnt)
dat_filt$timepnt <- factor(dat_filt$timepnt, levels=c("baseline", "peak", "memory"))
levels(dat_filt$timepnt)


#---------------------
# Add Dose
#---------------------
# Need to add in dose and schedule information
dat_filt <- data.frame(Dose = NA, dat_filt)

dat_filt[intersect(which(dat_filt$Study == "H56-035"),which(dat_filt$Group == 1)),1] <- 5
dat_filt[intersect(which(dat_filt$Study == "H56-035"),which(dat_filt$Group == 2)),1] <- 5
dat_filt[intersect(which(dat_filt$Study == "H56-035"),which(dat_filt$Group == 3)),1] <- 5
dat_filt[intersect(which(dat_filt$Study == "H56-035"),which(dat_filt$Group == 4)),1] <- 5

dat_filt[intersect(which(dat_filt$Study == "H1"),which(dat_filt$Group == 1)),1] <- 15
dat_filt[intersect(which(dat_filt$Study == "H1"),which(dat_filt$Group == 2)),1] <- 50
dat_filt[intersect(which(dat_filt$Study == "H1"),which(dat_filt$Group == 3)),1] <- 15
dat_filt[intersect(which(dat_filt$Study == "H1"),which(dat_filt$Group == 4)),1] <- 50

dat_filt[intersect(which(dat_filt$Study == "H56-032"),which(dat_filt$Group == 1)),1] <- 50
dat_filt[intersect(which(dat_filt$Study == "H56-032"),which(dat_filt$Group == 2)),1] <- 15
dat_filt[intersect(which(dat_filt$Study == "H56-032"),which(dat_filt$Group == 3)),1] <- 50

# check
dat_filt[which(is.na(dat_filt$Dose)==T),]
# convert to factor
dat_filt$Dose <- as.factor(dat_filt$Dose)


#---------------------
# Add Schedule
#---------------------
# Need to add in dose and schedule information
dat_filt <- data.frame(Schedule = NA, dat_filt)

dat_filt[intersect(which(dat_filt$Study == "H56-035"),which(dat_filt$Group == 1)),1] <- 2
dat_filt[intersect(which(dat_filt$Study == "H56-035"),which(dat_filt$Group == 2)),1] <- 3
dat_filt[intersect(which(dat_filt$Study == "H56-035"),which(dat_filt$Group == 3)),1] <- 2
dat_filt[intersect(which(dat_filt$Study == "H56-035"),which(dat_filt$Group == 4)),1] <- 3

dat_filt[intersect(which(dat_filt$Study == "H1"),which(dat_filt$Group == 1)),1] <- 2
dat_filt[intersect(which(dat_filt$Study == "H1"),which(dat_filt$Group == 2)),1] <- 2
dat_filt[intersect(which(dat_filt$Study == "H1"),which(dat_filt$Group == 3)),1] <- 1
dat_filt[intersect(which(dat_filt$Study == "H1"),which(dat_filt$Group == 4)),1] <- 1

dat_filt[intersect(which(dat_filt$Study == "H56-032"),which(dat_filt$Group == 1)),1] <- 3
dat_filt[intersect(which(dat_filt$Study == "H56-032"),which(dat_filt$Group == 2)),1] <- 3
dat_filt[intersect(which(dat_filt$Study == "H56-032"),which(dat_filt$Group == 3)),1] <- 3

# check 
dat_filt[which(is.na(dat_filt$Schedule)==T),]
# convert to factor
dat_filt$Schedule <- as.factor(dat_filt$Schedule)
head(dat_filt)



summary(dat_filt)

#-------------------------- 
# EXPORT DATA
#-------------------------- 


write.csv(dat_filt,"Filtered_cd4_counts_dataset_1.csv", row.names = FALSE) 
