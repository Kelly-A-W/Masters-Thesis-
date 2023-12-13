# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

#---------------------
# Load Libraries
#---------------------
#install.packages("BiocManager")
#library(BiocManager)
#BiocManager::install("COMPASS")
library(COMPASS)
library(ggplot2)
library(dplyr)
library(rstatix) # wilcox test
library(ggpubr) # add p-values to plots

# Get rid of scientific notation:
options(scipen=999) 

#---------------------
# Import Data
#---------------------

cd4_counts <- read.csv("clean_data_with_demographic_info.csv", check.names = F)
str(cd4_counts)

#---------------------
# Prepare Data
#---------------------

# Remove PHA:
counts <- cd4_counts[-which(cd4_counts$Stim=="PHA"),]

# Sum over memory markers:
df = counts[,14:29] + counts[,30:45] + counts[,46:61] + counts[,62:77]

# Change column names to correct format:
colnames(df) = c("G..L2..L17..TNF",     "G..L2..L17..!TNF",    "G..L2..!L17..TNF",   
                 "G..L2..!L17..!TNF",   "G..!L2..L17..TNF",    "G..!L2..L17..!TNF",  
                 "G..!L2..!L17..TNF",   "G..!L2..!L17..!TNF",  "!G..L2..L17..TNF",   
                 "!G..L2..L17..!TNF",   "!G..L2..!L17..TNF",   "!G..L2..!L17..!TNF", 
                 "!G..!L2..L17..TNF",   "!G..!L2..L17..!TNF",  "!G..!L2..!L17..TNF", 
                 "!G..!L2..!L17..!TNF")
# check no more c's left:
sum(grepl("c", colnames(df), fixed = TRUE)) # 0, so none !
# replace .. with &
colnames(df) = gsub("..", "&", colnames(df), fixed = TRUE)

# Add column with subject PID and time point
df$iid = paste(counts$PID, counts$timepnt)

# Add columns:
df$Stim = counts$Stim
df$Schedule = counts$Schedule
df$Dose = counts$Dose
df$QFT = counts$QFT
df$timepnt = counts$timepnt

# Set factor:
df$timepnt <- factor(df$timepnt, levels=c("baseline", "peak", "memory"))
df$QFT = as.factor(df$QFT)
levels(df$QFT) = c("QFT-", "QFT+")
df$Schedule = as.factor(df$Schedule)
df$Dose = as.factor(df$Dose)

# Subset by Stimulus:
df_esat6 = df[df[,"Stim"]=="ESAT6" | df[,"Stim"]=="UNS", ]
df_esat6 <- df_esat6[order(df_esat6$iid), ] # reorder by iid

df_ag85B = df[df[,"Stim"]=="Ag85B" | df[,"Stim"]=="UNS", ]
df_ag85B <- df_ag85B[order(df_ag85B$iid), ] # reorder by iid


#---------------------
# Perform COMPASS
#---------------------

# Create Count Matrices
ns_esat6 = df_esat6[df_esat6[,"Stim"]=="ESAT6",1:16] 
nu_esat6 = df_esat6[df_esat6[,"Stim"]=="UNS",1:16]
rownames(ns_esat6) = unique(df_esat6$iid) # set row names to iid
rownames(nu_esat6) = rownames(ns_esat6)

# Create meta data frame:
meta_esat6 = df_esat6[, c(17,19,20,21,22)] # Stim removed
meta_esat6 = meta_esat6[!duplicated(meta_esat6), ] # remove duplicated rows (we only want one row per subject)

# Fit COMPASS:
comp_esat6 = SimpleCOMPASS(n_s = ns_esat6, n_u = nu_esat6, meta = meta_esat6, individual_id="iid")




# Get scores:
FS <- FunctionalityScore(comp_esat6)
PFS <- PolyfunctionalityScore(comp_esat6)


# Plot a heatmap of the mean probability of response, to visualize differences 
## in expression for each category
plot(comp_esat6, fontsize = 13)
plot(comp_esat6, fontsize = 13, row_annotation = "QFT")
plot(comp_esat6, fontsize = 13, row_annotation = "timepnt")
plot(comp_esat6, fontsize = 13, row_annotation = "Schedule")
plot(comp_esat6, fontsize = 13, row_annotation = "Dose")
plot(comp_esat6, fontsize = 14, row_annotation = c("timepnt", "QFT"), pheatmap(cluster_cols=F))


df_esat6_compass = comp_esat6$data$meta
# check ordering:
sum(names(FS) != df_esat6_compass$iid)
sum(names(PFS) != df_esat6_compass$iid) # both 0, so both are in the same order !
# Add columns:
df_esat6_compass$FS = FS
df_esat6_compass$PFS = PFS


# Extract Posterior Probabilities:
sum(rownames(comp_esat6$fit$mean_gamma)!= df_esat6_compass$iid) # match ! 
# prepare column names
pp = comp_esat6$fit$mean_gamma
colnames(pp) = c("G..L2..L17..TNF",     "G..L2..L17..TNFc",    "G..L2..L17c..TNF",   
                 "G..L2..L17c..TNFc",   "G..L2c..L17..TNF",    "G..L2c..L17..TNFc",  
                 "G..L2c..L17c..TNF",   "G..L2c..L17c..TNFc",  "Gc..L2..L17..TNF",   
                 "Gc..L2..L17..TNFc",   "Gc..L2..L17c..TNF",   "Gc..L2..L17c..TNFc", 
                 "Gc..L2c..L17..TNF",   "Gc..L2c..L17..TNFc",  "Gc..L2c..L17c..TNF", 
                 "Gc..L2c..L17c..TNFc")
# check no more !'s left:
sum(grepl("!", colnames(pp), fixed = TRUE)) # 0, so none

df_esat6_compass_all = data.frame(df_esat6_compass, pp)
head(df_esat6_compass_all) 

#-----------------------
# Export Data
#-----------------------
write.csv(df_esat6_compass_all,"compass_esat6_dataset_1.csv", row.names = FALSE) 




