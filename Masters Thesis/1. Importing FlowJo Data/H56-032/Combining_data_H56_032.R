# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()


#---------------------
# Combine Datasets
#---------------------


pids <- c("08406", "07705", "06204", "08707", "03001", "02402", "02803", "03814", "01015", "09708", "10309", 
          "05013", "02916", "09317", "07818", "00210","00411", "01712", "08819", "04921", "07620", "12123",
          "08922", "10524", "11625")


#-------------------------- 
# IMPORT DATA AND COMBINE
#-------------------------- 

# CD4 COUNTS
count_cd4 <- data.frame()
for(i in 1:length(pids)){
  x         <- read.csv(paste0("cd4_counts_pid_", pids[i], ".csv"), check.names = F)
  count_cd4 <- rbind(count_cd4,x)
}
str(count_cd4)


# CD8 COUNTS
count_cd8 <- data.frame()
for(i in 1:length(pids)){
  x         <- read.csv(paste0("cd8_counts_pid_", pids[i], ".csv"), check.names = F)
  count_cd8 <- rbind(count_cd8,x)
}
str(count_cd8)


# Check for NAs
count_cd4[which(is.na(count_cd4$Stim)==T),1:6]   # none
count_cd8[which(is.na(count_cd8$Stim)==T),1:6]   # none

# number of pids
length(unique(count_cd4$PID))

# CHECK FOR DUPLICATE PIDS
nrow(count_cd4)-nrow(unique(count_cd4)) # 0, therefore no duplicated PIDs
nrow(count_cd8)-nrow(unique(count_cd8))




#-------------------------- 
# ADD GROUP DETAILS 
#-------------------------- 

grp1 <- c(08406, 07705, 06204, 08707, 03001, 02402, 02803, 09708)
grp2 <- c(03814, 01015, 10309, 05013, 02916, 00210, 00411, 01712)
grp3 <- c(09317, 07818, 08819, 04921, 07620, 12123, 08922, 10524, 11625)


# Check that there are no duplicates within each group
# 0 means no duplicates
which(duplicated(grp1) == T)
which(duplicated(grp2) == T)
which(duplicated(grp3) == T)
# no duplicates !

# Check that no pids have been assigned to multiple groups
# 0 means no pids have been assigned to multiple groups
which(duplicated(c(grp1, grp2, grp3)) == T)
# 0, therefore fine !

# Check if there are pids in our dataset that are not assigned to a group
setdiff(unique(count_cd4$PID), c(grp1, grp2, grp3))   #setdiff(a,b) tells you whats in a but not b
# 0, therefore none

# Check if there are any pids in our group assignments that are not in our dataset
setdiff(c(grp1, grp2, grp3), unique(count_cd4$PID))
# 0, therefore none


#--------------------------------------------------------------------------------Add group labels to dataset

# CD4 Counts:
count_cd4 <-data.frame(Group=NA,
                       count_cd4)
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp1))),1] <- 1
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp2))),1] <- 2
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp3))),1] <- 3
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
# Any NAs left ?
length(na.omit(count_cd8$Group))==length(count_cd8$Group)
# T, therefore all PIDs have been assigned to a group



#-------------------------- 
# Add QFT category
#-------------------------- 

QFTpos <- c(00210, 00411, 01015, 01712, 02916, 03814, 05013, 10309, 04921, 07620, 07818,
            08819, 08922, 09317, 10524, 11625, 12123)
QFTneg<- c(02402, 02803, 03001, 06204, 07705, 08406, 08707, 09708)


# Check that there are no duplicates within each group
# 0 means no duplicates
which(duplicated(QFTneg) == T)
which(duplicated(QFTpos) == T)
# no duplicates !

# Check that no pids have been assigned to multiple groups
# 0 means no pids have been assigned to multiple groups
which(duplicated(c(QFTneg, QFTpos)) == T)
# 0, therefore fine !

# Check if there are pids in our dataset that are not assigned to a group
setdiff(unique(count_cd4$PID), c(QFTneg, QFTpos))   #setdiff(a,b) tells you whats in a but not b
# 0, therefore none

# Check if there are any pids in our group assignments that are not in our dataset
setdiff(c(QFTneg, QFTpos), unique(count_cd4$PID))
# 0, therefore none


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

length(unique(count_cd8$PID[count_cd8[,"Group"]==1 & count_cd8[,"QFT"]=="QFT neg"]))
length(unique(count_cd8$PID[count_cd8[,"Group"]==1 & count_cd8[,"QFT"]=="QFT pos"]))
length(unique(count_cd8$PID[count_cd8[,"Group"]==2 & count_cd8[,"QFT"]=="QFT neg"]))
length(unique(count_cd8$PID[count_cd8[,"Group"]==2 & count_cd8[,"QFT"]=="QFT pos"]))
length(unique(count_cd8$PID[count_cd8[,"Group"]==3 & count_cd8[,"QFT"]=="QFT neg"]))
length(unique(count_cd8$PID[count_cd8[,"Group"]==3 & count_cd8[,"QFT"]=="QFT pos"]))


# Now check missing days and stimuli
df_rows <- count_cd4 %>% 
  group_by(PID) %>%
  summarise(no_rows = length(PID))
df_rows$PID[which(df_rows$no_rows<12)]
count_cd4[which(count_cd4$PID==411),1:6]  # only has day 0 visit


df_rows <- count_cd8 %>% 
  group_by(PID) %>%
  summarise(no_rows = length(PID))
df_rows$PID[which(df_rows$no_rows<12)]
count_cd8[which(count_cd8$PID==411),1:6]  








#-------------------------- 
# EXPORT DATA
#-------------------------- 


write.csv(count_cd4, "cd4_counts_H56_032.csv", row.names = FALSE) 

write.csv(count_cd8, "cd8_counts_H56_032.csv", row.names = FALSE)  


