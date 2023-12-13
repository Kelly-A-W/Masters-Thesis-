# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load Libraries
library(dplyr)

#---------------------
# Import Data
#---------------------

# CD4
cd4_counts_H1      <- read.csv("cd4_counts_H1.csv", check.names = F)
cd4_counts_H56_032 <- read.csv("cd4_counts_H56_032.csv", check.names = F)
cd4_counts_H56_035 <- read.csv("cd4_counts_H56_035.csv", check.names = F)


# CD8
cd8_counts_H1      <- read.csv("cd8_counts_H1.csv", check.names = F)
cd8_counts_H56_032 <- read.csv("cd8_counts_H56_032.csv", check.names = F)
cd8_counts_H56_035 <- read.csv("cd8_counts_H56_035.csv", check.names = F)



#---------------------
# Merge Data
#---------------------

# First add columns indicating what study the data belongs to
cd4_counts_H1      <- data.frame(Study="H1", cd4_counts_H1)
cd4_counts_H56_032 <- data.frame(Study="H56-032", cd4_counts_H56_032)
cd4_counts_H56_035 <- data.frame(Study="H56-035", cd4_counts_H56_035)
cd8_counts_H1      <- data.frame(Study="H1", cd8_counts_H1)
cd8_counts_H56_032 <- data.frame(Study="H56-032", cd8_counts_H56_032)
cd8_counts_H56_035 <- data.frame(Study="H56-035", cd8_counts_H56_035)

# Check if pids are unique between studies
intersect(unique(cd4_counts_H1$PID), unique(cd4_counts_H56_032$PID))
intersect(unique(cd4_counts_H1$PID), unique(cd4_counts_H56_035$PID))
intersect(unique(cd4_counts_H56_035$PID), unique(cd4_counts_H56_032$PID))
# two pids are the same between studies: 210 between H1 and H56-032 and 3001 for H56-035 and H56-032
# Give these four participants unique pids:
cd4_counts_H1$PID[which(cd4_counts_H1$PID == 210)]            <- "210 H1"
cd4_counts_H56_032$PID[which(cd4_counts_H56_032$PID == 210)]  <- "210 H56-032"
cd4_counts_H56_032$PID[which(cd4_counts_H56_032$PID == 3001)] <- "3001 H56-032"
cd4_counts_H56_035$PID[which(cd4_counts_H56_035$PID == 3001)] <- "3001 H56-035"
# Check:
intersect(unique(cd4_counts_H1$PID), unique(cd4_counts_H56_032$PID))
intersect(unique(cd4_counts_H1$PID), unique(cd4_counts_H56_035$PID))
intersect(unique(cd4_counts_H56_035$PID), unique(cd4_counts_H56_032$PID))
# 0, so sorted !
# Do the same for CD8:
intersect(unique(cd8_counts_H1$PID), unique(cd8_counts_H56_032$PID))
intersect(unique(cd8_counts_H1$PID), unique(cd8_counts_H56_035$PID))
intersect(unique(cd8_counts_H56_035$PID), unique(cd8_counts_H56_032$PID))
# Again, two pids are the same between studies: 210 between H1 and H56-032 and 3001 for H56-035 and H56-032
# Give these four participants unique pids:
cd8_counts_H1$PID[which(cd8_counts_H1$PID == 210)]            <- "210 H1"
cd8_counts_H56_032$PID[which(cd8_counts_H56_032$PID == 210)]  <- "210 H56-032"
cd8_counts_H56_032$PID[which(cd8_counts_H56_032$PID == 3001)] <- "3001 H56-032"
cd8_counts_H56_035$PID[which(cd8_counts_H56_035$PID == 3001)] <- "3001 H56-035"
# Check:
intersect(unique(cd8_counts_H1$PID), unique(cd8_counts_H56_032$PID))
intersect(unique(cd8_counts_H1$PID), unique(cd8_counts_H56_035$PID))
intersect(unique(cd8_counts_H56_035$PID), unique(cd8_counts_H56_032$PID))
# 0, so sorted !

# Check columns line up:
which(names(cd4_counts_H1)     !=  names(cd4_counts_H56_032))
which(names(cd4_counts_H1)     !=  names(cd4_counts_H56_035))
which(names(cd4_counts_H56_035)!=  names(cd4_counts_H56_032))
which(names(cd8_counts_H1)     !=  names(cd8_counts_H56_032))
which(names(cd8_counts_H1)     !=  names(cd8_counts_H56_035))
which(names(cd8_counts_H56_035)!=  names(cd8_counts_H56_032))
# all zero, thus all columns line up
# so all we have to do is put dataframes one under each other 


# Combine datasets 
cd4_counts_dat1 <- rbind(cd4_counts_H1, cd4_counts_H56_032, cd4_counts_H56_035)
cd8_counts_dat1 <- rbind(cd8_counts_H1, cd8_counts_H56_032, cd8_counts_H56_035)


# Check Dimensions
dim(cd4_counts_dat1)[1] == dim(cd4_counts_H1)[1] + dim(cd4_counts_H56_032)[1] + dim(cd4_counts_H56_035)[1]
dim(cd4_counts_dat1)[2] == dim(cd4_counts_H1)[2]
dim(cd8_counts_dat1)[1] == dim(cd8_counts_H1)[1] + dim(cd8_counts_H56_032)[1] + dim(cd8_counts_H56_035)[1]
dim(cd8_counts_dat1)[2] == dim(cd8_counts_H1)[2]
# T, therefore dimensions line up


#-------------------------------------
# Remove Non-consenting Participants
#-------------------------------------
# 3 participants from trial H1 did not consent for long-term storage of their samples
# so need to remove pid 81, 94 and 131

cd4_counts_dat1 <- cd4_counts_dat1[-which(cd4_counts_dat1$PID == 81),]
cd4_counts_dat1 <- cd4_counts_dat1[-which(cd4_counts_dat1$PID == 94),]
cd4_counts_dat1 <- cd4_counts_dat1[-which(cd4_counts_dat1$PID == 131),]


#---------------------
# Export Data
#---------------------


write.csv(cd4_counts_dat1, "cd4_counts_dataset_1.csv", row.names = FALSE) 

write.csv(cd8_counts_dat1, "cd8_counts_dataset_1.csv", row.names = FALSE)  





