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

# Get rid of scientific notation:
options(scipen=999) 

#---------------------
# Import Data
#---------------------

cd4_counts <- read.csv("clean_cd4_counts_dataset_1.csv", check.names = F)
str(cd4_counts)

#-------------------------------------
# Calculate TRF
#-------------------------------------

# Only keep counts for cells expressing at least one of INFg, IL2, TNF or L17
cd4_trf <- cd4_counts[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd4_counts), fixed = TRUE)==T)]
str(cd4_trf)

# Convert to Frequencies
cd4_trf[,8:67] <- (cd4_trf[,8:67]/cd4_trf[,7])*100

# Background Subtract
cd4_bstrf <- data.frame()
for(i in 1:length(unique(cd4_counts$PID))){
  pid <- unique(cd4_counts$PID)[i]
  dat <- cd4_trf[which(cd4_trf$PID == pid),]
  days <- unique(dat$Day)
  for (j in 1:length(days)) {
    day <- days[j]
    df  <- dat[which(dat$Day == day),]
    df[-which(df$Stim == "UNS"),8:67]  <- df[-which(df$Stim == "UNS"),8:67] - df[rep(which(df$Stim == "UNS"), (length(df$Stim)-1)),8:67]
    cd4_bstrf <- rbind(cd4_bstrf,df)
  }
}
summary(cd4_bstrf)
which(is.na(cd4_bstrf)==T)


# Deal with Negative Frequencies:
for (i in 8:ncol(cd4_bstrf)) {
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
sapply(sign(cd4_bstrf[,-(1:7)]), table)

# Sum across rows
cd4_bstrf <- cd4_bstrf %>%
  mutate(Total = select(., 8:67) %>% rowSums(na.rm = TRUE))
str(cd4_bstrf)

# Remove unnecessary rows
cd4_bstrf <- cd4_bstrf[,c(1,2,3,4,5,6,68)] 
head(cd4_bstrf)

# Reshape Data
cd4_bstrf <- reshape(cd4_bstrf, idvar = c("Study", "QFT",  "Group", "PID", "Day"), timevar = "Stim", direction = "wide")
head(cd4_bstrf)
length(unique(cd4_bstrf$PID))
length(unique(cd4_counts$PID)) 

# Set Factors
cd4_bstrf$Day <- as.factor(cd4_bstrf$Day)


#-----------------------------------------------------------
# Extract Ag85B TRF (named TVR because less changes to code)
#-----------------------------------------------------------

TVR <- cd4_bstrf$Total.ESAT6
cd4_tvr <- data.frame(cd4_bstrf, TVR=TVR)
head(cd4_tvr)
cd4_tvr <- cd4_tvr[,-c(6,7,8,9)]
cd4_tvr <- na.omit(cd4_tvr)

#---------------------
# Add Time Category
#---------------------
cd4_tvr <- data.frame(timepnt = NA, cd4_tvr)

cd4_tvr$timepnt[which(cd4_tvr$Day==0)] <- "baseline"
cd4_tvr$timepnt[which(cd4_tvr$Day==14)] <- "peak"
cd4_tvr$timepnt[which(cd4_tvr$Day==70)] <- "peak"
cd4_tvr$timepnt[which(cd4_tvr$Day==126)] <- "peak"
cd4_tvr$timepnt[which(cd4_tvr$Day==224)] <- "memory"
cd4_tvr$timepnt[which(cd4_tvr$Day==292)] <- "memory"
cd4_tvr$timepnt[which(cd4_tvr$Day==210)] <- "memory"

which(is.na(cd4_tvr$timepnt) == T)

cd4_tvr$timepnt <- as.factor(cd4_tvr$timepnt)
cd4_tvr$timepnt <- factor(cd4_tvr$timepnt, levels=c("baseline", "peak", "memory"))
levels(cd4_tvr$timepnt)


#---------------------
# Add Dose
#---------------------
# Need to add in dose and schedule information
cd4_tvr <- data.frame(cd4_tvr, Dose = NA)

cd4_tvr[intersect(which(cd4_tvr$Study == "H56-035"),which(cd4_tvr$Group == 1)),8] <- 5
cd4_tvr[intersect(which(cd4_tvr$Study == "H56-035"),which(cd4_tvr$Group == 2)),8] <- 5
cd4_tvr[intersect(which(cd4_tvr$Study == "H56-035"),which(cd4_tvr$Group == 3)),8] <- 5
cd4_tvr[intersect(which(cd4_tvr$Study == "H56-035"),which(cd4_tvr$Group == 4)),8] <- 5

cd4_tvr[intersect(which(cd4_tvr$Study == "H1"),which(cd4_tvr$Group == 1)),8] <- 15
cd4_tvr[intersect(which(cd4_tvr$Study == "H1"),which(cd4_tvr$Group == 2)),8] <- 50
cd4_tvr[intersect(which(cd4_tvr$Study == "H1"),which(cd4_tvr$Group == 3)),8] <- 15
cd4_tvr[intersect(which(cd4_tvr$Study == "H1"),which(cd4_tvr$Group == 4)),8] <- 50

cd4_tvr[intersect(which(cd4_tvr$Study == "H56-032"),which(cd4_tvr$Group == 1)),8] <- 50
cd4_tvr[intersect(which(cd4_tvr$Study == "H56-032"),which(cd4_tvr$Group == 2)),8] <- 15
cd4_tvr[intersect(which(cd4_tvr$Study == "H56-032"),which(cd4_tvr$Group == 3)),8] <- 50

# check
cd4_tvr[which(is.na(cd4_tvr$Dose)==T),]
# convert to factor
cd4_tvr$Dose <- as.factor(cd4_tvr$Dose)


#---------------------
# Add Schedule
#---------------------
# Need to add in dose and schedule information
cd4_tvr <- data.frame(cd4_tvr, Schedule = NA)

cd4_tvr[intersect(which(cd4_tvr$Study == "H56-035"),which(cd4_tvr$Group == 1)),9] <- 2
cd4_tvr[intersect(which(cd4_tvr$Study == "H56-035"),which(cd4_tvr$Group == 2)),9] <- 3
cd4_tvr[intersect(which(cd4_tvr$Study == "H56-035"),which(cd4_tvr$Group == 3)),9] <- 2
cd4_tvr[intersect(which(cd4_tvr$Study == "H56-035"),which(cd4_tvr$Group == 4)),9] <- 3

cd4_tvr[intersect(which(cd4_tvr$Study == "H1"),which(cd4_tvr$Group == 1)),9] <- 2
cd4_tvr[intersect(which(cd4_tvr$Study == "H1"),which(cd4_tvr$Group == 2)),9] <- 2
cd4_tvr[intersect(which(cd4_tvr$Study == "H1"),which(cd4_tvr$Group == 3)),9] <- 1
cd4_tvr[intersect(which(cd4_tvr$Study == "H1"),which(cd4_tvr$Group == 4)),9] <- 1

cd4_tvr[intersect(which(cd4_tvr$Study == "H56-032"),which(cd4_tvr$Group == 1)),9] <- 3
cd4_tvr[intersect(which(cd4_tvr$Study == "H56-032"),which(cd4_tvr$Group == 2)),9] <- 3
cd4_tvr[intersect(which(cd4_tvr$Study == "H56-032"),which(cd4_tvr$Group == 3)),9] <- 3

# check 
cd4_tvr[which(is.na(cd4_tvr$Schedule)==T),]
# convert to factor
cd4_tvr$Schedule <- as.factor(cd4_tvr$Schedule)
head(cd4_tvr)

#---------------------
# Median
#---------------------

# create function to calculate median
med.fun <- function(data, idx)
{
  df <- data[idx,]
  median(df[,7])
}


meds <- data.frame(Dose = c(rep(5,12), rep(15,15), rep(50, 18)),
                   Schedule = c(rep(2,6),rep(3,6), rep(1,6), rep(2,6), rep(3,3), rep(1,6), rep(2,6), rep(3,6)),
                   QFT = c(rep(c(rep("QFT neg",3), rep("QFT pos",3)), 4), rep("QFT pos",3), rep(c(rep("QFT neg",3), rep("QFT pos",3)), 3)),
                   timepnt = rep(c("baseline", "peak", "memory"),15),
                   Medians = NA, 
                   LowerQ = NA,
                   UpperQ = NA)


for (i in 1:nrow(meds)) {
  
  y <- cd4_tvr[cd4_tvr[,"Dose"]== meds$Dose[i] & cd4_tvr[,"Schedule"]== meds$Schedule[i] & cd4_tvr[,"QFT"]==meds$QFT[i] & cd4_tvr[,"timepnt"]==meds$timepnt[i], "TVR"]
  
  meds[i,"Medians"] <- median(y)
  
  meds[i,"LowerQ"] <- as.numeric(summary(y)[2])
  meds[i,"UpperQ"] <- as.numeric(summary(y)[5])
  
  
}

head(meds)

# convert to factors
meds$timepnt <- as.factor(meds$timepnt)
meds$timepnt <- factor(meds$timepnt, levels=c("baseline", "peak", "memory"))
levels(meds$timepnt)
meds$Dose <- as.factor(meds$Dose)
meds$Schedule <- as.factor(meds$Schedule)



############################### ONLY MEDIANS ##################################################

# Create Grob function for highlighting significant p-values later:
find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}



# Select 5mg and QFT-:
dftest <- cd4_tvr[cd4_tvr[,"Dose"]== 5 &  cd4_tvr[,"QFT"]=="QFT neg",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(TVR ~ Schedule))
# Adjusted P-values:
(padj <- p.adjust(stat.test$p, method = "BH"))
# Create Plot:
p <- ggplot(meds[meds[,"Dose"]== 5 & meds[,"QFT"]=="QFT neg",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("TVR (%)") +
  xlab("")+
  ylim(0,0.9)+
  scale_x_discrete(labels = c("", "", "")) +
  scale_color_manual(values=c("#FF0000","#660000"), name = "Admin.") 
# Make grobtable:
df1         <- data.frame(Test = "2 vs 3", Baseline = padj[1], 
                          Peak = padj[2], Memory = padj[3])
names(df1)  <- c("Test", "Baseline", "Peak", "Memory")
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 25))
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Test), limits = c(rev(df1$Test), ""))+
  theme(axis.text.y = element_text(face="bold", size=25))
# Patchwork magic:
p + p3 + plot_layout(ncol = 1, heights = c(12, 1.5))



# Select 5mg and QFT+:
dftest <- cd4_tvr[cd4_tvr[,"Dose"]== 5 &  cd4_tvr[,"QFT"]=="QFT pos",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(TVR ~ Schedule))
# Adjusted P-values:
(padj <- p.adjust(stat.test$p, method = "BH"))
# Create Plot:
p <- ggplot(meds[meds[,"Dose"]== 5 & meds[,"QFT"]=="QFT pos",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("TVR (%)") +
  xlab("")+
  ylim(0,0.9)+
  scale_x_discrete(labels = c("", "", "")) +
  scale_color_manual(values=c("blue", "#000066"), name = "Admin.")
# Make grobtable:
df1         <- data.frame(Test = "2 vs 3", Baseline = padj[1], 
                          Peak = padj[2], Memory = padj[3])
names(df1)  <- c("Test", "Baseline", "Peak", "Memory")
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 25))
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Test), limits = c(rev(df1$Test), ""))+
  theme(axis.text.y = element_text(face="bold", size=25))
# Patchwork magic
p + p3 + plot_layout(ncol = 1, heights = c(12, 1.5))



# Select 15mg and QFT-:
dftest <- cd4_tvr[cd4_tvr[,"Dose"]== 15 &  cd4_tvr[,"QFT"]=="QFT neg",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(TVR ~ Schedule))
# Adjusted P-values:
(padj <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
p <- ggplot(meds[meds[,"Dose"]== 15 & meds[,"QFT"]=="QFT neg",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("TVR (%)") +
  xlab("")+
  ylim(0,0.9)+
  scale_x_discrete(labels = c("", "", "")) +
  scale_color_manual(values=c("#FF9999","#FF0000"), name = "Admin.")
# Make grobtable:
df1         <- data.frame(Test = "1 vs 2", Baseline = padj[1], 
                          Peak = padj[2], Memory = padj[3])
names(df1)  <- c("Test", "Baseline", "Peak", "Memory")
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 25))
# Set significant p-values to bold:

# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Test), limits = c(rev(df1$Test), ""))+
  theme(axis.text.y = element_text(face="bold", size=25))
# Patchwork magic:
p + p3 + plot_layout(ncol = 1, heights = c(12, 1.5))


# Select 15mg and QFT+:
dftest     <- cd4_tvr[cd4_tvr[,"Dose"]== 15 &  cd4_tvr[,"QFT"]=="QFT pos",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(TVR ~ Schedule))
# Adjusted P-values:
(padj      <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot: 
p <- ggplot(meds[meds[,"Dose"]== 15 & meds[,"QFT"]=="QFT pos",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("TVR (%)") +
  xlab("")+
  ylim(0,0.9)+
  scale_x_discrete(labels = c("", "", "")) +
  scale_color_manual(values=c("#99CCFF", "blue", "#000066"), name = "Admin.")
# Make grobtable:
df1         <- data.frame(Test = c("1 vs 2", "1 vs 3", "2 vs 3"), Baseline = padj[1:3], 
                          Peak = padj[4:6], Memory = padj[7:9]) 
names(df1)  <- c("Test", "Baseline", "Peak", "Memory")
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 25))
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Test), limits = c(rev(df1$Test), ""))+
  theme(axis.text.y = element_text(face="bold", size=25))
# Patchwork magic:
p + p3 + plot_layout(ncol = 1, heights = c(12, 3))



# Select 50mg and QFT-:
dftest     <- cd4_tvr[cd4_tvr[,"Dose"]== 50 &  cd4_tvr[,"QFT"]=="QFT neg",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(TVR ~ Schedule))
# Adjusted P-values:
(padj      <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot: 
p <- ggplot(meds[meds[,"Dose"]== 50 & meds[,"QFT"]=="QFT neg",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,  # Bottom margin
                             l = 10)) +
  ylab("TVR (%)") +
  xlab("")+
  ylim(0,0.9)+
  scale_x_discrete(labels = c("", "", "")) +
  scale_color_manual(values=c("#FF9999","#FF0000","#660000"), name = "Admin.")
# Make grobtable:
df1         <- data.frame(Test = c("1 vs 2", "1 vs 3", "2 vs 3"), Baseline = padj[1:3], 
                          Peak = padj[4:6], Memory = padj[7:9]) 
names(df1)  <- c("Test", "Baseline", "Peak", "Memory")
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 25))
# Set significant p-values to bold:
ind1 <- find_cell(tab, 2, 2, "core-fg")
ind2 <- find_cell(tab, 2, 3, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind2][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Test),limits = c(rev(df1$Test), ""))+
  theme(axis.text.y = element_text(face="bold", size=25))
# Patchwork magic
p + p3 + plot_layout(ncol = 1, heights = c(12, 3))



# Select 50mg and QFT+:
dftest <- cd4_tvr[cd4_tvr[,"Dose"]== 50 &  cd4_tvr[,"QFT"]=="QFT pos",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(TVR ~ Schedule))
# Adjusted P-values:
(padj      <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot: 
p <- ggplot(meds[meds[,"Dose"]== 50 & meds[,"QFT"]=="QFT pos",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("TVR (%)") +
  xlab("")+
  scale_x_discrete(labels = c("", "", "")) +
  ylim(0,0.9)+
  scale_color_manual(values=c("#99CCFF", "blue", "#000066"), name = "Admin.")
# Make grobtable:
df1         <- data.frame(Test = c("1 vs 2", "1 vs 3", "2 vs 3"), Baseline = padj[1:3], 
                          Peak = padj[4:6], Memory = padj[7:9]) 
names(df1)  <- c("Test", "Baseline", "Peak", "Memory")
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 25))
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Test), limits = c(rev(df1$Test), ""))+
  theme(axis.text.y = element_text(face="bold", size=25))
# Patchwork magic
p + p3 + plot_layout(ncol = 1, heights = c(12, 3))






#----------------------------------
# AUC
#----------------------------------

str(cd4_tvr)
dat         <- cd4_tvr
dat$timepnt <- as.numeric(dat$timepnt)
dat$PID     <- as.factor(dat$PID)

dat$AUC <- NA

for(i in 1:length(levels(dat$PID))){
  pid <- levels(dat$PID)[i]
  a   <- AUC(dat$timepnt[which(dat$PID == pid)], dat$TVR[which(dat$PID == pid)])
  dat$AUC[which(dat$PID == pid)] <- a
}

str(dat)



# Remove duplicates
dat <- dat %>% distinct()

# Only want one row for each pid:
# so remove the three things that differ within each pid: timpnt, TVR and Day
dat <- dat[,-c(1,6,7)]
# now remove duplicates 
dat <- dat %>% distinct()
head(dat)
length(unique(dat$PID)) == nrow(dat) # T, so only one row per pid

################# PLOT ############################################


# QFT-:
dftest     <- dat[dat[,"QFT"]=="QFT neg",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(Dose) %>%
    wilcox_test(AUC ~ Schedule))
stat.test <- stat.test %>%
  add_xy_position(x = "Dose", dodge = 0.8)
stat.test$y.position <- c(1.5, 2.2, 2.4, 2.6, 2.2, 2.4, 2.6)
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
ggplot(dftest, aes(y=AUC, x=Dose)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = Schedule)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("TRF AUC") +
  xlab(expression(paste("Concentration (", mu, "g)")))+
  ylim(0,2.8)+
  scale_fill_manual(values=c("#FF9999","#FF0000","#660000"), name="Admin.")+
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     bracket.nudge.y = 0.5, 
                     size = 10, 
                     bracket.size = 0.8, 
                     step.increase = 0.1) 




# QFT+:
dftest <- dat[dat[,"QFT"]=="QFT pos",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(Dose) %>%
    wilcox_test(AUC ~ Schedule))
stat.test <- stat.test %>%
  add_xy_position(x = "Dose", dodge = 0.8)
stat.test$y.position <- c(1.5, 2.2, 2.4, 2.6, 2.2, 2.4, 2.6)
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
ggplot(dftest, aes(y=AUC, x=Dose)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = Schedule)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("TRF AUC") +
  xlab(expression(paste("Concentration (", mu, "g)")))+
  ylim(0,2.8)+
  scale_fill_manual(values=c("#99CCFF", "blue", "#000066"), name="Admin.")+
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     #bracket.nudge.y = -0.5, 
                     size = 10, 
                     bracket.size = 0.8, 
                     step.increase = 0.01) 
# One NA in dataframe because PID 411 only has baseline measurements, thus AUC cannot be calculated



#----------------------------------
# Select Dose 2
#----------------------------------

dat_d2     <- cd4_tvr[which(cd4_tvr$Schedule==2),]
dat_d2$QFT <- as.factor(dat_d2$QFT)


# Plot Peak:
dftest     <- dat_d2[dat_d2[,"timepnt"]=="peak",]
# Calculate p-values:
stat.test  <- dftest %>%
  group_by(QFT) %>%
  wilcox_test(TVR ~ Dose) 
(stat.test <- stat.test %>%
    add_xy_position(x = "QFT", dodge = 0.8))
stat.test$y.position <- c(1.1, 1.25, 1.4, 1.1, 1.25, 1.4)
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
ggplot(dftest, aes(y=TVR, x=QFT)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Dose), color ="#003300") +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("TVR (%)") +
  xlab("")+
  ylim(0,1.45)+
  scale_fill_manual(values=c("#FFFFFF","#00CC66","#006633"), name="", 
                    labels=c(expression(paste("5", mu, "g")), 
                             expression(paste("15", mu, "g")), 
                             expression(paste("50", mu, "g"))))+ 
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8) 



# Plot Memory:
dftest <- dat_d2[dat_d2[,"timepnt"]=="memory",]
# Calculate p-values:
stat.test  <- dftest %>%
  group_by(QFT) %>%
  wilcox_test(TVR ~ Dose) 
(stat.test <- stat.test %>%
    add_xy_position(x = "QFT", dodge = 0.8))
stat.test$y.position <- c(1.1, 1.25, 1.4, 1.1, 1.25, 1.4)
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
ggplot(dftest, aes(y=TVR, x=QFT)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Dose), color ="#330066") +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("TVR (%)") +
  xlab("")+
  ylim(0,1.45)+
  scale_fill_manual(values=c("#FFFFFF","#CC99FF","#660099"), name="", 
                    labels=c(expression(paste("5", mu, "g")), 
                             expression(paste("15", mu, "g")), 
                             expression(paste("50", mu, "g"))))+ 
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8)

