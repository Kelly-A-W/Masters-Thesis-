# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load Libraries
library(dplyr)
library(tidyr)

#---------------------
# Combine Datasets
#---------------------


pids <- c(1002, 1004, 1007, 1013, 1016, 1019, 1025, 1026, 1027, 1033, 1034, 1040, 1045,
          1047, 1048, 2001, 2002, 2006, 2007, 2008, 2009, 2010, 2012, 2014, 2015,
          2016, 3001, 3007, 3008, 3009, 3011, 3019, 3023, 3024, 3026, 3028, 3030, 
          3006, 3013, 3022, 3031, 3032, 3003, 3005, 3014, 3015, 3017, 3027, 3020)

dates <- c(2022016, 4022016, 8122015, 9122015, 18122015, 22122015, 10022016)


#-------------------------- 
# IMPORT DATA AND COMBINE
#-------------------------- 

# CD4 COUNTS
count_cd4 <- data.frame()
for(i in 1:length(pids)){
  x         <- read.csv(paste0("cd4_counts_pid_", pids[i], ".csv"), check.names = F)
  count_cd4 <- rbind(count_cd4,x)
}
for(i in 1:length(dates)){
  x         <- read.csv(paste0("cd4_counts_date_", dates[i], ".csv"), check.names = F)
  count_cd4 <- rbind(count_cd4,x)
}
str(count_cd4)


# CD8 COUNTS
count_cd8 <- data.frame()
for(i in 1:length(pids)){
  x         <- read.csv(paste0("cd8_counts_pid_", pids[i], ".csv"), check.names = F)
  count_cd8 <- rbind(count_cd8,x)
}
for(i in 1:length(dates)){
  x         <- read.csv(paste0("cd8_counts_date_", dates[i], ".csv"), check.names = F)
  count_cd8 <- rbind(count_cd8,x)
}
str(count_cd8)


# Check for NAs
count_cd4[which(is.na(count_cd4$Stim)==T),1:6]    # these NAs are fine, they are just for Rv stim
count_cd8[which(is.na(count_cd8$Stim)==T),1:6]   



# CHECK FOR DUPLICATE PIDS
nrow(count_cd4)-nrow(unique(count_cd4)) # 0, therefore no duplicated PIDs
nrow(count_cd8)-nrow(unique(count_cd8))
# CHECK FOR DUPLICATE ROWS
nrow(count_cd4)-nrow(unique(count_cd4[,-1])) # 0, therefore no duplicated rows
nrow(count_cd8)-nrow(unique(count_cd8[,-1]))


# Remove extra pids from counts_date_
# this is because when there are some unnecessary pids in the counts_dates fcs files that are not included in our study
count_cd4 <- count_cd4 %>% filter(PID %in% pids)
sum(unique(count_cd4$PID)[order(unique(count_cd4$PID))] != pids[order(pids)])  # 0, so pids vector is same as PID in dataframe
count_cd8 <- count_cd8 %>% filter(PID %in% pids)
sum(unique(count_cd8$PID)[order(unique(count_cd8$PID))] != pids[order(pids)]) 



#-------------------------- 
# ADD GROUP DETAILS 
#-------------------------- 

grp1 <- c(1002, 1004, 1007, 1013, 1016, 1019, 1025, 1026, 1027, 1033, 1034, 1040, 1045, 1047, 1048)
grp2 <- c(2001, 2002, 2006, 2007, 2008, 2009, 2010, 2012, 2014, 2015, 2016)
grp3 <- c(3003, 3005, 3006, 3013, 3014, 3015, 3017, 3020, 3022, 3027, 3031, 3032)
grp4 <- c(3001, 3007, 3008, 3009, 3011, 3019, 3023, 3024, 3026, 3028, 3030)


# Check that there are no duplicates within each group
# 0 means no duplicates
which(duplicated(grp1) == T)
which(duplicated(grp2) == T)
which(duplicated(grp3) == T)
which(duplicated(grp4) == T)
# no duplicates !

# Check that no pids have been assigned to multiple groups
# 0 means no pids have been assigned to multiple groups
which(duplicated(c(grp1, grp2, grp3, grp4)) == T)
# 0, therefore fine !

# Check if there are pids in our dataset that are not assigned to a group
setdiff(unique(count_cd4$PID), c(grp1, grp2, grp3, grp4))   #setdiff(a,b) tells you whats in a but not b
# 0, therefore none

# Check if there are any pids in our group assignments that are not in our dataset
setdiff(c(grp1, grp2, grp3, grp4), unique(count_cd4$PID))
# 0, therefore none

#--------------------------------------------------------------------------------Add group labels to dataset

# CD4 Counts:
count_cd4 <-data.frame(Group=NA,
                       count_cd4)
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp1))),1] <- 1
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp2))),1] <- 2
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp3))),1] <- 3
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp4))),1] <- 4
# Any NAs left ?
length(na.omit(count_cd4$Group))==length(count_cd4$Group)
# T, therefore all PIDs have been assigned to a group
str(count_cd4)





# CD8 Counts:
count_cd8 <-data.frame(Group=NA,
                       count_cd8)
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% grp1))),1] <- 1
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% grp2))),1] <- 2
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% grp3))),1] <- 3
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% grp4))),1] <- 4
# Any NAs left ?
length(na.omit(count_cd8$Group))==length(count_cd8$Group)
# T, therefore all PIDs have been assigned to a group





#-------------------------- 
# Add QFT category
#-------------------------- 

# QFT status info according to "SSI samples for R01.xlsx" sheet, Table 2

QFTneg <- c(grp1, grp2)
QFTpos <- c(grp3,grp4)


# CD4 Counts:
count_cd4 <-data.frame(QFT=NA,
                       count_cd4)
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% QFTpos))),1] <- "QFT pos"
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% QFTneg))),1] <- "QFT neg"

# Any NAs left ?
length(na.omit(count_cd4$QFT))==length(count_cd4$QFT)
# T, therefore all PIDs have been assigned to a group
str(count_cd4)


# CD8 Counts:
count_cd8 <-data.frame(QFT=NA,
                       count_cd8)
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% QFTpos))),1] <- "QFT pos"
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% QFTneg))),1] <- "QFT neg"
# Any NAs left ?
length(na.omit(count_cd8$QFT))==length(count_cd8$QFT)
# T, therefore all PIDs have been assigned to a group






#-------------------------------------------
# CHECK ALL PIDS AND OBSERVATIONS ARE HERE
#-------------------------------------------
# how many pids do we have for each dataset ?
length(unique(count_cd4$PID))
length(unique(count_cd8$PID))

# How many are in each group?
length(unique(count_cd4$PID[count_cd4[,"Group"]==1 & count_cd4[,"QFT"]=="QFT neg"]))
length(unique(count_cd4$PID[count_cd4[,"Group"]==1 & count_cd4[,"QFT"]=="QFT pos"]))
length(unique(count_cd4$PID[count_cd4[,"Group"]==2 & count_cd4[,"QFT"]=="QFT neg"]))
length(unique(count_cd4$PID[count_cd4[,"Group"]==2 & count_cd4[,"QFT"]=="QFT pos"]))
length(unique(count_cd4$PID[count_cd4[,"Group"]==3 & count_cd4[,"QFT"]=="QFT neg"]))
length(unique(count_cd4$PID[count_cd4[,"Group"]==3 & count_cd4[,"QFT"]=="QFT pos"]))
length(unique(count_cd4$PID[count_cd4[,"Group"]==4 & count_cd4[,"QFT"]=="QFT neg"]))
length(unique(count_cd4$PID[count_cd4[,"Group"]==4 & count_cd4[,"QFT"]=="QFT pos"]))

length(unique(count_cd8$PID[count_cd8[,"Group"]==1 & count_cd8[,"QFT"]=="QFT neg"]))
length(unique(count_cd8$PID[count_cd8[,"Group"]==1 & count_cd8[,"QFT"]=="QFT pos"]))
length(unique(count_cd8$PID[count_cd8[,"Group"]==2 & count_cd8[,"QFT"]=="QFT neg"]))
length(unique(count_cd8$PID[count_cd8[,"Group"]==2 & count_cd8[,"QFT"]=="QFT pos"]))
length(unique(count_cd8$PID[count_cd8[,"Group"]==3 & count_cd8[,"QFT"]=="QFT neg"]))
length(unique(count_cd8$PID[count_cd8[,"Group"]==3 & count_cd8[,"QFT"]=="QFT pos"]))
length(unique(count_cd8$PID[count_cd8[,"Group"]==4 & count_cd8[,"QFT"]=="QFT neg"]))
length(unique(count_cd8$PID[count_cd8[,"Group"]==4 & count_cd8[,"QFT"]=="QFT pos"]))



# Now check missing days and stimuli
df_rows <- count_cd4 %>% 
  group_by(PID) %>%
  summarise(no_rows = length(PID))
df_rows$PID[which(df_rows$no_rows<12)]
count_cd4[which(count_cd4$PID==3001),1:6] # Missing Ag85B day 0
count_cd4[which(count_cd4$PID==3020),1:6] # Missing day 292
count_cd4[which(count_cd4$PID==3024),1:6] # Missing Ag85B day 0
count_cd4[which(count_cd4$PID==3026),1:6] # Missing UNS day 0



df_rows <- count_cd8 %>% 
  group_by(PID) %>%
  summarise(no_rows = length(PID))
df_rows$PID[which(df_rows$no_rows<12)]













#-------------------------- 
# EXPORT DATA
#-------------------------- 


write.csv(count_cd4, "cd4_counts_H56_035.csv", row.names = FALSE) 

write.csv(count_cd8, "cd8_counts_H56_035.csv", row.names = FALSE)  



