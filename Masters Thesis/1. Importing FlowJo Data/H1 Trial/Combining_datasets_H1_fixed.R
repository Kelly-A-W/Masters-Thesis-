# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()


#---------------------
# Load Libraries
#---------------------
library(tidyr)
library(dplyr)


#---------------------
# Combine Datasets
#---------------------


pids <- c("001","002","004","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022",
          "023","025","026","027","028","030","031","032","035","036","041","043","044","046","047","049","053","055",
          "056","061","063","064","065","066","067","068","071","072","073","074","075","076","077","079","080","081",
          "082","083","084","086","087","088","089","090","091","092","094","095","096","097","098","099","100","101",
          "102","103","105","106","107","108","109","112","113","114","115","116","117","118","119","120","123",124, "125",
          "126","127","128","129","130","131","132",133, "134","136","137","139","140","141","142",145,"147",148,149,150,
          151,152,153,154,155,156,164,165,167,168,169,170,171,172,173,174,175,179,180,181,183,184,185,186,190,192,193,194,
          195,196,197,198,199,205,206,207,209,210,212,213,214,215,216,217,218,219,221,222,223,224,225,226,228,229,231,
          232,233,234,235,236,237,238,239,240)


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

# check for NAs
sum(is.na(count_cd4$Stim))
sum(is.na(count_cd8$Stim))
sum(is.na(count_cd4$PID))
sum(is.na(count_cd8$PID))
sum(is.na(count_cd4$Day))
sum(is.na(count_cd8$Day))


# CHECK FOR DUPLICATE PIDS
nrow(count_cd4)-nrow(unique(count_cd4)) # 0, therefore no duplicated PIDs
nrow(count_cd8)-nrow(unique(count_cd8))




#-------------------------- 
# ADD GROUP DETAILS 
#-------------------------- 


# From "H1 GROUPS PBMC FCS.xlsx":
# Only looking at pids with PBMC samples available ("PBMC = OK" in excel sheet)
grp1 <- c(101, 105, 113, 117, 14, 15, 16, 19, 2, 22, 28, 46, 47, 53, 56, 67, 71, 74, 79, 80, 86, 90, 98, 99, 126,
          131, 147, 152, 172, 190, 193, 195, 198, 206, 225, 229, 234)
grp2 <- c(102, 104, 109, 118, 21, 23, 26, 30, 32, 41, 43, 44, 66, 68, 73, 77, 8, 82, 85, 88, 9, 91, 95, 123, 125, 129, 
          130, 132, 136, 139, 140, 142, 148, 149, 151, 153, 156, 168, 174, 183, 194, 212, 213, 214, 215, 218, 221, 231, 
          232, 237)
grp3 <- c(1, 106, 108, 112, 114, 116, 120, 13, 17, 18, 20, 31, 36, 55, 61, 72, 75, 83, 84, 89, 92, 94, 127, 128, 134, 
          137, 150, 154, 165, 167, 169, 171, 175, 181, 186, 192, 197, 199, 205, 207, 209, 210, 222, 226, 228, 235, 238)
grp4 <- c(10, 100, 103, 107, 11, 115, 119, 12, 25, 27, 35, 4, 49, 63, 64, 65, 76, 81, 87, 96, 97, 121, 122, 124, 133, 
          141, 145, 155, 164, 170, 173, 179, 180, 184, 185, 196, 216, 217, 219, 223, 224, 233, 236, 239, 240)

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
# 104  85 121 122
# these pids were all flagged as red on HD, hence why they were excluded


# Now lets check another excel sheet just to make sure !

# From "VisitsDone_EN screening vs enrolment pID.xlsx":
# Only looking at pids with workspaces available ("WSP = Y" in excel sheet)
x1 <- c(002, 014, 015, 016, 019, 022, 028, 046, 047, 053, 056, 067, 071, 074, 079, 080, 086, 090, 098, 099, 101, 105,
        113, 117, 126, 131, 147, 152, 172, 190, 193, 195, 198, 206, 225, 229, 234)
x2 <- c(008, 009, 021, 023, 026, 030, 032, 041, 043, 044, 066, 068, 073, 077, 082, 085, 088, 091, 095, 102, 104, 109,
        118, 123, 125, 129, 130, 132, 136, 139, 140, 142, 148, 149, 151, 153, 156, 168, 174, 183, 194, 212, 213, 214,
        215, 218, 221, 231, 232, 237)
x3 <- c(001, 013, 017, 018, 020, 031, 036, 055, 061, 072, 075, 083, 084, 089, 092, 094, 106, 108, 112, 114, 116, 120,
        127, 128, 134, 137, 150, 154, 165, 167, 169, 171, 175, 181, 186, 192, 197, 199, 205, 207, 209, 210, 222, 226,
        228, 235, 238)
x4 <- c(004, 010, 011, 012, 025, 027, 035, 049, 063, 064, 065, 076, 081, 087, 096, 097, 100, 103, 107, 115, 119, 121, 
        122, 124, 133, 141, 145, 155, 164, 170, 173, 179, 180, 184, 185, 196, 216, 217, 219, 223, 224, 233, 236, 239, 
        240)

# GROUP 1
setdiff(grp1, x1)
setdiff(x1, grp1)
# both zero, so no difference in group 1 assignment between the two excel sheets

# GROUP 2
setdiff(grp2, x2) 
setdiff(x2, grp2)
# both zero, so no difference in group 2 assignment between the two excel sheets

# GROUP 3
setdiff(grp3, x3) # none
setdiff(x3, grp3)
# both zero, so no difference in group 3 assignment between the two excel sheets

# GROUP 4
setdiff(grp4, x4) # none
setdiff(x4, grp4)
# both zero, so no difference in group 4 assignment between the two excel sheets





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


QFTneg <- c(101, 105, 113, 117, 14, 15, 16, 19, 2, 22, 28, 46, 47, 53, 56, 67, 71, 74, 79, 80, 86, 90, 98, 99,102,
            104, 109, 118, 21, 23, 26, 30, 32, 41, 43, 44, 66, 68, 73, 77, 8, 82, 85, 88, 9, 91, 95, 1, 106, 108, 
            112, 114, 116, 120, 13, 17, 18, 20, 31, 36, 55, 61, 72, 75, 83, 84, 89, 92, 94,10, 100, 103, 107, 11, 
            115, 119, 12, 25, 27, 35, 4, 49, 63, 64, 65, 76, 81, 87, 96, 97)


QFTpos <- c(126,131, 147, 152, 172, 190, 193, 195, 198, 206, 225, 229, 234, 123, 125, 129, 130, 132, 136, 139, 140, 
            142, 148, 149, 151, 153, 156, 168, 174, 183, 194, 212, 213, 214, 215, 218, 221, 231, 232, 237, 127, 128, 134, 
            137, 150, 154, 165, 167, 169, 171, 175, 181, 186, 192, 197, 199, 205, 207, 209, 210, 222, 226, 228, 235, 238, 
            121, 122, 124, 133, 141, 145, 155, 164, 170, 173, 179, 180, 184, 185, 196, 216, 217, 219, 223, 224, 233, 
            236, 239, 240)



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
# 104  85 121 122
# these pids were all flagged as red on HD, hence why they were excluded



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









#-------------------------- 
# FIX MISSLABELED DAYS
#-------------------------- 
# For groups three and four, I labeled the peak response day as day 70 when it should actually be day 14!

a <- which(count_cd4$Group==3)
b <- which(count_cd4$Day==70)
count_cd4$Day[intersect(a,b)] <- 14

a <- which(count_cd4$Group==4)
b <- which(count_cd4$Day==70)
count_cd4$Day[intersect(a,b)] <- 14




# CD8: 
a <- which(count_cd8$Group==3)
b <- which(count_cd8$Day==70)
count_cd8$Day[intersect(a,b)] <- 14

a <- which(count_cd8$Group==4)
b <- which(count_cd8$Day==70)
count_cd8$Day[intersect(a,b)] <- 14


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
count_cd4[which(count_cd4$PID==173),1:6]  # missing PHA for day 0 and day 14
count_cd4[which(count_cd4$PID==184),1:6]  # missing PHA for all three days
count_cd4[which(count_cd4$PID==185),1:6]  # missing PHA for day 0
count_cd4[which(count_cd4$PID==198),1:6]  # missing PHA for day 70
count_cd4[which(count_cd4$PID==223),1:6]  # missing PHA for day 0 and 224

df_rows <- count_cd8 %>% 
  group_by(PID) %>%
  summarise(no_rows = length(PID))
df_rows$PID[which(df_rows$no_rows<12)]
count_cd8[which(count_cd8$PID==173),1:6]  # missing PHA for day 0 and day 14
count_cd8[which(count_cd8$PID==184),1:6]  # missing PHA for all three days
count_cd8[which(count_cd8$PID==185),1:6]  # missing PHA for day 0
count_cd8[which(count_cd8$PID==198),1:6]  # missing PHA for day 70
count_cd8[which(count_cd8$PID==223),1:6]  # missing PHA for day 0 and 224




#-------------------------- 
# EXPORT DATA
#-------------------------- 


write.csv(count_cd4, "cd4_counts_H1.csv", row.names = FALSE) 


write.csv(count_cd8, "cd8_counts_H1.csv", row.names = FALSE)  






