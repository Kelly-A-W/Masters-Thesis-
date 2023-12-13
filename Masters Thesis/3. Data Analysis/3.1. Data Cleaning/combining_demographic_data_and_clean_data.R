# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load Libraries
library(dplyr)
library("readxl")

#---------------------
# Import Data
#---------------------

h1dat <- read_excel("Elisa_H1_032_035Data_11072023.xlsx", sheet =1)
h56032dat  <- read_excel("Elisa_H1_032_035Data_11072023.xlsx", sheet =3)
h56035dat  <- read_excel("Elisa_H1_032_035Data_11072023.xlsx", sheet =2)

cd4_counts <- read.csv("clean_cd4_counts_dataset_1.csv", check.names = F)
str(cd4_counts)

#---------------------
# Add Time Category
#---------------------
cd4_counts <- data.frame(timepnt = NA, cd4_counts)

cd4_counts$timepnt[which(cd4_counts$Day==0)]   <- "baseline"
cd4_counts$timepnt[which(cd4_counts$Day==14)]  <- "peak"
cd4_counts$timepnt[which(cd4_counts$Day==70)]  <- "peak"
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





#---------------------
# Merge Data
#---------------------

# Create demographic columns:
cd4_counts <- data.frame(sex = NA, cd4_counts)
cd4_counts <- data.frame(ethnicity = NA, cd4_counts)
cd4_counts <- data.frame(age = NA, cd4_counts)


# Find the PIDs that had special (non-numeric) labels:
cd4_counts$PID[which(is.na(as.numeric(cd4_counts$PID))==T)]


# Create new data without these special labels (will add back in later):
fulldata <- cd4_counts[-which(is.na(as.numeric(cd4_counts$PID))==T),]
fulldata$PID[which(is.na(as.numeric(fulldata$PID))==T)] # check done, =0 thus all special labels were removed


# H1 Trial:
as.numeric(fulldata[fulldata[,"Study"]=="H1",10])
as.numeric(h1dat$ScreenNo)
# Check that there are no pids in our data that we don't have demographic data for: 
setdiff(unique(as.numeric(fulldata[fulldata[,"Study"]=="H1",10])), as.numeric(h1dat$ScreenNo))   #setdiff(a,b) tells you whats in a but not b.  # 0, therefore none
# Check which pids we have demographic data for but are not in our data:
setdiff(as.numeric(h1dat$ScreenNo), unique(as.numeric(fulldata[fulldata[,"Study"]=="H1",10])))  
sum(1:240 != as.numeric(h1dat$ScreenNo)) # = 0, thus our pids in demo data for H1 is the same as 1:240, thus the pids == position in dataset
# Only keep the demographic data for the PIDs we have:
temp_demo <- h1dat[intersect(as.numeric(fulldata[fulldata[,"Study"]=="H1",10]), as.numeric(h1dat$ScreenNo)),]
# Check:
sum(unique(as.numeric(fulldata[fulldata[,"Study"]=="H1",10])) != as.numeric(temp_demo$ScreenNo)) # 0, so everything matches !
# Merge H1:
for(i in 1:dim(temp_demo)[1]){
  pid <- as.numeric(temp_demo$ScreenNo[i])
  rows <- intersect(which(fulldata$Study=="H1"), which(as.numeric(fulldata$PID)==pid))
  fulldata[rows, 1] <- temp_demo$AgeYrs[i]
  fulldata[rows, 2] <- temp_demo$Ethnicity[i]
  fulldata[rows, 3] <- temp_demo$Sex[i]
}
# Check:
which(is.na(fulldata[fulldata[,"Study"]=="H1",1])==T)
which(is.na(fulldata[fulldata[,"Study"]=="H1",2])==T)
which(is.na(fulldata[fulldata[,"Study"]=="H1",3])==T) # All are zero, thus all H1 rows were populated



# H56:032 Trial:
as.numeric(fulldata[fulldata[,"Study"]=="H56-032",10])
as.numeric(h56032dat$PTID)
# Check that there are no pids in our data that we don't have demographic data for: 
setdiff(unique(as.numeric(fulldata[fulldata[,"Study"]=="H56-032",10])), as.numeric(h56032dat$PTID))   #setdiff(a,b) tells you whats in a but not b.  # 0, therefore none
# Check which pids we have demographic data for but are not in our data:
(diff <- setdiff(as.numeric(h56032dat$PTID), unique(as.numeric(fulldata[fulldata[,"Study"]=="H56-032",10]))))
# Only keep the demographic data for the PIDs we have:
temp_demo <- h56032dat[-c(which(as.numeric(h56032dat$PTID)==diff[1]), which(as.numeric(h56032dat$PTID)==diff[2])),]
# Check:
length(unique(as.numeric(fulldata[fulldata[,"Study"]=="H56-032",10]))) == length(as.numeric(temp_demo$PTID)) # T, so looks good
# Merge H56-032:
for(i in 1:dim(temp_demo)[1]){
  pid  <- as.numeric(temp_demo$PTID[i])
  rows <- intersect(which(fulldata$Study=="H56-032"), which(as.numeric(fulldata$PID)==pid))
  fulldata[rows, 1] <- temp_demo$AgeYrs[i]
  fulldata[rows, 2] <- temp_demo$Ethnicity[i]
  fulldata[rows, 3] <- temp_demo$Sex[i]
}
# Check:
which(is.na(fulldata[fulldata[,"Study"]=="H56-032",1])==T)
which(is.na(fulldata[fulldata[,"Study"]=="H56-032",2])==T)
which(is.na(fulldata[fulldata[,"Study"]=="H56-032",3])==T) # All are zero, thus all H1 rows were populated



# H56:035 Trial:
as.numeric(fulldata[fulldata[,"Study"]=="H56-035",10])
as.numeric(h56035dat$PTID)
# Check that there are no pids in our data that we don't have demographic data for: 
setdiff(unique(as.numeric(fulldata[fulldata[,"Study"]=="H56-035",10])), as.numeric(h56035dat$PTID))   #setdiff(a,b) tells you whats in a but not b.  # 0, therefore none
# Check which pids we have demographic data for but are not in our data:
(diff <- setdiff(as.numeric(h56035dat$PTID), unique(as.numeric(fulldata[fulldata[,"Study"]=="H56-035",10]))))
# Only keep the demographic data for the PIDs we have:
rerows <- c()
for (i in 1:length(diff)) { # find all the rows in demographic data that I dont need
  rerows <- append(rerows, which(as.numeric(h56035dat$PTID)==diff[i]))
}
rerows
temp_demo <- h56035dat[-rerows,]
# Check:
length(unique(as.numeric(fulldata[fulldata[,"Study"]=="H56-035",10]))) == length(as.numeric(temp_demo$PTID)) # T, so looks good
# Merge H56-035:
for(i in 1:dim(temp_demo)[1]){
  pid  <- as.numeric(temp_demo$PTID[i])
  rows <- intersect(which(fulldata$Study=="H56-035"), which(as.numeric(fulldata$PID)==pid))
  fulldata[rows, 1] <- temp_demo$AgeYrs[i]
  fulldata[rows, 2] <- temp_demo$Ethnicity[i]
  fulldata[rows, 3] <- temp_demo$Sex[i]
}
# Check:
which(is.na(fulldata$age)==T)
which(is.na(fulldata$ethnicity)==T)
which(is.na(fulldata$sex)==T) 
# PID 3032 has missing sex data - was checked on spreadsheet, this info is indeed missing



# Create df with the rows with special PID labels:
df1 <- data.frame(age=h1dat$AgeYrs[210], ethnicity=h1dat$Ethnicity[210], sex=h1dat$Sex[210],
                  cd4_counts[which(cd4_counts$PID=="210 H1"),-c(1,2,3)])
df2 <- data.frame(age=h56032dat$AgeYrs[1], ethnicity=h56032dat$Ethnicity[1], sex=h56032dat$Sex[1],
                  cd4_counts[which(cd4_counts$PID=="3001 H56-032"),-c(1,2,3)])
df3 <- data.frame(age=h56032dat$AgeYrs[10], ethnicity=h56032dat$Ethnicity[10], sex=h56032dat$Sex[10],
                  cd4_counts[which(cd4_counts$PID=="210 H56-032"),-c(1,2,3)])
df4 <- data.frame(age=h56035dat$AgeYrs[49], ethnicity=h56035dat$Ethnicity[49], sex=h56035dat$Sex[49],
                  cd4_counts[which(cd4_counts$PID=="3001 H56-035"),-c(1,2,3)])

# Merge special rows to full data:
fulldata <- rbind(fulldata, df1, df2, df3, df4)

# Check:
which(is.na(fulldata$age)==T)
which(is.na(fulldata$ethnicity)==T)
which(is.na(fulldata$sex)==T) 

#---------------------
# Export Data
#---------------------


write.csv(fulldata, "clean_data_with_demographic_info.csv", row.names = FALSE) 



