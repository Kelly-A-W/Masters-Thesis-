# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

# Load Libraries
library(dplyr)
library(ggplot2)


#---------------------
# Import Data
#---------------------

cd4_counts <- read.csv("cd4_counts_dataset_1.csv", check.names = F)

# Remove rows with NA stim (i.e. stim = RV2660)
cd4_counts <- cd4_counts[-which(is.na(cd4_counts$Stim) == T),]


#---------------------
# Point (a):
#---------------------
# "Unstimulated control is present and interpretable for each set of samples."
# First, remove all visits for which UNS is not available
# Then Calculate Total Response Frequencies (TRF)
# Plot UNS TRF for each visit and identify outliers
# Look at Anele's comments to decide if outliers should be removed 

#---------------------------------------------- Remove Visits with missing UNS

# Create a dataframe with only UNS rows
df_unst <- cd4_counts[which(cd4_counts$Stim =="UNS"),]

# Count how many rows each PID has
df_rows <- df_unst %>% 
  group_by(PID) %>%
  summarise(no_rows = length(PID))

# PIDs with less than 3 unstim rows:
(unsmis <- df_rows$PID[which(df_rows$no_rows<3)]) 
cd4_counts[which(cd4_counts$PID==unsmis[1]),c(4,5,6)] # only 2 days, none missing unstim, so keep
cd4_counts[which(cd4_counts$PID==unsmis[2]),c(4,5,6)] # remove day 0
cd4_counts[which(cd4_counts$PID==unsmis[3]),c(4,5,6)] # only 1 day, unstim not missing, so keep

# Remove day 0 visit for pid 3026 from H56-035 trial:
cd4_counts <- cd4_counts[-which(cd4_counts$PID == unsmis[2] & cd4_counts$Day == 0), ]



#---------------------------------------------- Calculate TRF

# Only keep counts for cells expressing at least one of INFg, IL2,TNF or L17
cd4_trf <- cd4_counts[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd4_counts), fixed = TRUE)==T)]
str(cd4_trf)

# Sum across rows
cd4_trf <- cd4_trf %>%
  mutate(Total = select(., 8:67) %>% rowSums(na.rm = TRUE))

# Convert to frequencies
cd4_trf$Total <- (cd4_trf$Total/cd4_trf$CD4count)*100

# Remove unnecessary rows
cd4_trf <- cd4_trf[,c(1,2,3,4,5,6,68)] 
head(cd4_trf)

# Reshape Data
cd4_trf <- reshape(cd4_trf, idvar = c("Study", "QFT",  "Group", "PID", "Day"), timevar = "Stim", direction = "wide")
head(cd4_trf)


#---------------------------------------------- Plot TRF

ggplot(cd4_trf,aes(y=Total.UNS, x=PID)) +
  geom_point() +
  geom_text(aes(label=ifelse(Total.UNS>0.2,PID,'')),hjust=.7,vjust=-.6, size=7)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        text = element_text(size = 30)) +
  ylab("Unstimulated Total Response Frequencies")
  

ggplot(cd4_trf,aes(y=Total.UNS, x=PID, color=Study)) +
  geom_point() +
  geom_text(aes(label=ifelse(Total.UNS>0.2,PID,'')),hjust=.7,vjust=-.6, size=7)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        text = element_text(size = 30)) +
  ylab("Unstimulated Total Response Frequencies")




df_out <- cd4_trf[which(cd4_trf$Total.UNS>0.2),-c(7,8,9)]
df_out <- df_out[order(df_out$Total.UNS, decreasing = T),]
df_out <- data.frame(df_out, 
                     Comment = c("NC", "NC", "Missing PHA stim for all time points", "High IFNg background", "NC",
                                 "Time Gate D70 on axis", "High INFg background on Day 224", "High IFNg background",
                                 "High background on UNS", "NC", "Background on UNS", "Missing PHA stim for all time points",
                                 "NC", "NC", "No row in excel sheet!", 
                                 "Background on UNS", "Very high background CD8; CCR7 gate on CD4 too high"))
df_out

write.csv(df_out, "outliers.csv", row.names = FALSE) 
# For now keep all flagged visits


#---------------------
# Point (b):
#---------------------
# "Frequencies of PHA-induced total cytokine-expressing CD4 or CD8 T cells are greater than the 
# median + 3MAD (median absolute deviation) of the total cytokine+ CD4 or CD8 T cells
# of the unstimulated controls of the entire cohort."
# PHA is the positive control
# Method:
# 1. Calculate Median TRF for CD4 unstim
# 2. Calculate MAD
# 3. Apply cut-off median + 3MAD


#------------------------------------ 1. Median Unstim NTR
# median unstimulated TRF for cd4 

(med_cd4 <- median(cd4_trf$Total.UNS))


#------------------------------------ 2. Median Absolute Deviation (MAD)
# calculate the difference between the median and the actual TRF each visit
# take the absolute value of this
# then take the median of this

(mad_cd4 <- median(abs(cd4_trf$Total.UNS - med_cd4)))


#------------------------------------ 4. Apply cut-off median + 3MAD

# Calculate cut-off:
(cut_cd4 <- med_cd4 + 3*mad_cd4)



# Apply cut-off:
# Check PHA
(rem_pha <- cd4_trf[which(cd4_trf$Total.PHA <= cut_cd4), ])  # 2 visits don't make cut-off
# Check Ag85B
rem_pha[which(rem_pha$Total.Ag85B <= cut_cd4), ]   # None
# Therefore both pass and we keep both visits



# Check which PHA TRF = NA
(pha_na <- which(is.na(cd4_trf$Total.PHA) == T))
cd4_trf[pha_na,]
# use Ag85B instead
(rem_ag85b <- cd4_trf[pha_na[which(cd4_trf$Total.Ag85B[pha_na] <= cut_cd4)],])  # two don't pass
# check ESAT6
(rem_esat6 <- cd4_trf[pha_na[which(cd4_trf$Total.ESAT6[pha_na] <= cut_cd4)],])  # different two don't pass
# So all pass 



#---------------------
# Point (c):
#---------------------
# "For each sample, the frequency of PHA-induced total cytokine+ CD4 or CD8 T cells are greater than the frequency 
# of the same cell population in its respective unstimulated control."

(rem <- cd4_trf[which(cd4_trf$Total.PHA <= cd4_trf$Total.UNS), ]) # two visits

# Check Ag85B and ESAT6:
rem[which(rem$Total.Ag85B <= rem$Total.UNS), ]  # None
rem[which(rem$Total.ESAT6 <= rem$Total.UNS), ]  # Day 14 for pid 170, H1 trial
# keep both since they both past Ag85B

#cd4_trf    <- cd4_trf[-which(cd4_trf$PID == rem$PID[1] & cd4_trf$Day == rem$Day[1]),]   # Remove day 14 for pid 170, study H1
#cd4_counts <- cd4_counts[-which(cd4_trf$PID == rem$PID[1] & cd4_trf$Day == rem$Day[1]),] 


# Check for PHA TRF = NA
cd4_trf[pha_na,]
# use Ag85B
(rem_ag85b <- cd4_trf[pha_na[which(cd4_trf$Total.Ag85B[pha_na] <= cd4_trf$Total.UNS[pha_na])],])  # none
# check ESAT6
(rem_esat6 <- cd4_trf[pha_na[which(cd4_trf$Total.ESAT6[pha_na] <= cd4_trf$Total.UNS[pha_na])],])  # one
# So pass all





#---------------------
# Export Data
#---------------------


write.csv(cd4_counts, "clean_cd4_counts_dataset_1.csv", row.names = FALSE) 










