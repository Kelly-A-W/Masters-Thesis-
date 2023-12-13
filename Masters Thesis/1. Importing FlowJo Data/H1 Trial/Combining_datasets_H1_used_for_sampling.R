# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()


#---------------------
# Combine Datasets
#---------------------


pids <- c("001","002","004","008","009","010","011","012","013","014","015","016","017","018","019","020","021","022",
          "023","025","026","027","028","030","031","032","035","036","041","043","044","046","047","049","053","055",
          "056","061","063","064","065","066","067","068","071","072","073","074","075","076","077","079","080","081",
          "082","083","084","086","087","088","089","090","091","092","094","095","096","097","098","099","100","101",
          "102","103","105","106","107","108","109","112","113","114","115","116","117","118","119","120","123","125",
          "126","127","128","129","130","131","132","134","136","137","139","140","141","142","147",148,149,150,
          151,152,153,154,155,156,164,165,167,168,169,170,171,172,173,174,179,180,181,183,184,185,186,190,192,193,194,
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


# CHECK FOR DUPLICATE PIDS
nrow(count_cd4)-nrow(distinct(count_cd4)) # 0, therefore no duplicated PIDs
nrow(count_cd8)-nrow(distinct(count_cd8))


# CD4 FREQUENCY
freq_cd4 <- data.frame()
for(i in 1:length(pids)){
  x         <- read.csv(paste0("cd4_frequency_pid_", pids[i], ".csv"), check.names = F)
  freq_cd4 <- rbind(freq_cd4,x)
}
str(freq_cd4)


# CD8 FREQUENCY
freq_cd8 <- data.frame()
for(i in 1:length(pids)){
  x         <- read.csv(paste0("cd8_frequency_pid_", pids[i], ".csv"), check.names = F)
  freq_cd8 <- rbind(freq_cd8,x)
}
str(freq_cd8)


# CHECK FOR DUPLICATE PIDS
nrow(freq_cd4)-nrow(distinct(freq_cd4)) # no duplicated PIDs 
nrow(freq_cd8)-nrow(distinct(freq_cd8)) # no duplicated PIDs 





#-------------------------- 
# ADD GROUP DETAILS 
#-------------------------- 

grp1_QFTneg <- unique(c(101,105,113,117,14,15,16,19,2,22,28,6,47,53,56,67,71,74,79,80,86,90,98,99,46))
grp1_QFTpos <- unique(c(126,131,147,152,172,190,193,195,198,206,225,229,234))
grp2_QFTneg <- unique(c(102,104,109,118,21,23,26,30,32,41,43,44,66,68,73,77,8,82,85,88,9,91,95))
grp2_QFTpos <- unique(c(123,125,129,130,132,136,139,140,142,148,149,151,153,156,168,174,183,194,212,213,214,215,218,
                        221,231,232,237))
grp3_QFTneg <- unique(c(1,106,108,112,114,116,120,13,17,18,20,31,36,55,61,72,75,83,84,89,92,94))
grp3_QFTpos <- unique(c(127,128,134,137,150,154,165,167,169,171,175,181,186,192,197,199,205,207,209,210,222,226,228,235,238))
grp4_QFTneg <- unique(c(10,100,103,107,11,115,119,12,25,27,35,4,49,63,64,65,76,81,87,96,97))
grp4_QFTpos <- unique(c(121,122,124,133,141,145,155,164,170,173,179,180,184,185,196,216,217,219,223,224,233,236,239,240))


# Check PIDs in multiple groups
pids <- c(grp1_QFTneg,grp1_QFTpos,grp2_QFTneg,grp2_QFTpos,grp3_QFTneg,grp3_QFTpos,grp4_QFTneg,grp4_QFTpos)
length(pids) == length(unique(pids))
# T, therefore no PIDs in more than one group



#--------------------------------------------------------------------------------Add group labels to dataset

# CD4 Counts:
count_cd4 <-data.frame(Group=NA,
                        count_cd4)
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp1_QFTneg))),1] <- "Group 1 QFTneg"
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp1_QFTpos))),1] <- "Group 1 QFTpos"
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp2_QFTneg))),1] <- "Group 2 QFTneg"
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp2_QFTpos))),1] <- "Group 2 QFTpos"
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp3_QFTneg))),1] <- "Group 3 QFTneg"
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp3_QFTpos))),1] <- "Group 3 QFTpos"
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp4_QFTneg))),1] <- "Group 4 QFTneg"
count_cd4[which(apply(data.frame(count_cd4$PID), 1, function(r) any(r %in% grp4_QFTpos))),1] <- "Group 4 QFTpos"
# Any NAs left ?
length(na.omit(count_cd4$Group))==length(count_cd4$Group)
# T, therefore all PIDs have been assigned to a group



# CD4 Frequency:
which(count_cd4[,c(2,3,4)] != freq_cd4[,c(1,2,3)])    # checking freq_cd4 rows match count_cd4 rows (they should !)
freq_cd4 <-data.frame(Group=count_cd4$Group, freq_cd4)
# check
str(freq_cd4)   
which(freq_cd4[,c(1,2,3,4)] != freq_cd4[,c(1,2,3,4)])    # 0, therefore looks good!






# CD8 Counts:
count_cd8 <-data.frame(Group=NA,
                       count_cd8)
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% grp1_QFTneg))),1] <- "Group 1 QFTneg"
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% grp1_QFTpos))),1] <- "Group 1 QFTpos"
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% grp2_QFTneg))),1] <- "Group 2 QFTneg"
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% grp2_QFTpos))),1] <- "Group 2 QFTpos"
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% grp3_QFTneg))),1] <- "Group 3 QFTneg"
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% grp3_QFTpos))),1] <- "Group 3 QFTpos"
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% grp4_QFTneg))),1] <- "Group 4 QFTneg"
count_cd8[which(apply(data.frame(count_cd8$PID), 1, function(r) any(r %in% grp4_QFTpos))),1] <- "Group 4 QFTpos"
length(na.omit(count_cd8$Group))==length(count_cd8$Group)
# T, therefore all PIDs have been assigned to a group



# CD8 Frequency:
which(count_cd8[,c(2,3,4)] != freq_cd8[,c(1,2,3)])    # checking freq_cd8 rows match count_cd8 rows (they should !)
freq_cd8 <-data.frame(Group=count_cd8$Group,freq_cd8)
# check
str(freq_cd8)   
which(count_cd8[,c(1,2,3,4)] != freq_cd8[,c(1,2,3,4)])    # looks good!






#-------------------------- 
# FIX MISSLABELED DAYS
#-------------------------- 
# For groups three and four, I labeled the peak response day as day 70 when it should actually be day 14!

a <- which(count_cd4$Group=="Group 3 QFTneg")
b <- which(count_cd4$Day==70)
count_cd4$Day[intersect(a,b)] <- 14
a <- which(count_cd4$Group=="Group 3 QFTpos")
b <- which(count_cd4$Day==70)
count_cd4$Day[intersect(a,b)] <- 14
a <- which(count_cd4$Group=="Group 4 QFTneg")
b <- which(count_cd4$Day==70)
count_cd4$Day[intersect(a,b)] <- 14
a <- which(count_cd4$Group=="Group 4 QFTpos")
b <- which(count_cd4$Day==70)
count_cd4$Day[intersect(a,b)] <- 14

freq_cd4$Day <- count_cd4$Day



# CD8: 
a <- which(count_cd8$Group=="Group 3 QFTneg")
b <- which(count_cd8$Day==70)
count_cd8$Day[intersect(a,b)] <- 14
a <- which(count_cd8$Group=="Group 3 QFTpos")
b <- which(count_cd8$Day==70)
count_cd8$Day[intersect(a,b)] <- 14
a <- which(count_cd8$Group=="Group 4 QFTneg")
b <- which(count_cd8$Day==70)
count_cd8$Day[intersect(a,b)] <- 14
a <- which(count_cd8$Group=="Group 4 QFTpos")
b <- which(count_cd8$Day==70)
count_cd8$Day[intersect(a,b)] <- 14

freq_cd8$Day <- count_cd8$Day



#-------------------------- 
# EXPORT DATA
#-------------------------- 


write.csv(count_cd4, "cd4_counts.csv", row.names = FALSE) 
write.csv(freq_cd4, "cd4_frequency.csv", row.names = FALSE) 

write.csv(count_cd8, "cd8_counts.csv", row.names = FALSE)  
write.csv(freq_cd8, "cd8_frequency.csv", row.names = FALSE) 






