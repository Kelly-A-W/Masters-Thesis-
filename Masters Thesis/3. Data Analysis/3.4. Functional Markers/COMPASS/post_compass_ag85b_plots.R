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
library(grid)
library(gridExtra)
library(patchwork)
library(RColorBrewer)

# Get rid of scientific notation:
options(scipen=999) 

#---------------------
# Import Data
#---------------------

scores <- read.csv("compass_ag85b_dataset_1.csv", check.names = F)
str(scores)

########################
# CALCULATE MEDIANS
########################

meds_FS <- data.frame(Dose = c(rep(5,12), rep(15,15), rep(50, 18)),
                      Schedule = c(rep(2,6),rep(3,6), rep(1,6), rep(2,6), rep(3,3), rep(1,6), rep(2,6), rep(3,6)),
                      QFT = c(rep(c(rep("QFT-",3), rep("QFT+",3)), 4), rep("QFT+",3), rep(c(rep("QFT-",3), rep("QFT+",3)), 3)),
                      timepnt = rep(c("baseline", "peak", "memory"),15),
                      Medians = NA, 
                      LowerQ = NA,
                      UpperQ = NA)
meds_PFS = meds_FS

for (i in 1:nrow(meds_FS)) {
  
  y <- scores[scores[,"Dose"]== meds_PFS$Dose[i] & scores[,"Schedule"]== meds_PFS$Schedule[i] & scores[,"QFT"]==meds_PFS$QFT[i] & scores[,"timepnt"]==meds_PFS$timepnt[i], ]
  
  meds_FS[i,"Medians"] <- median(y$FS)
  meds_FS[i,"LowerQ"] <- as.numeric(summary(y$FS)[2])
  meds_FS[i,"UpperQ"] <- as.numeric(summary(y$FS)[5])
  
  meds_PFS[i,"Medians"] <- median(y$PFS)
  meds_PFS[i,"LowerQ"] <- as.numeric(summary(y$PFS)[2])
  meds_PFS[i,"UpperQ"] <- as.numeric(summary(y$PFS)[5])
}

str(meds_FS)
str(meds_PFS)

# convert to factors
meds_FS$timepnt <- as.factor(meds_FS$timepnt)
meds_FS$timepnt <- factor(meds_FS$timepnt, levels=c("baseline", "peak", "memory"))
levels(meds_FS$timepnt)
meds_FS$Dose <- as.factor(meds_FS$Dose)
meds_FS$Schedule <- as.factor(meds_FS$Schedule)

meds_PFS$timepnt <- as.factor(meds_PFS$timepnt)
meds_PFS$timepnt <- factor(meds_PFS$timepnt, levels=c("baseline", "peak", "memory"))
levels(meds_PFS$timepnt)
meds_PFS$Dose <- as.factor(meds_PFS$Dose)
meds_PFS$Schedule <- as.factor(meds_PFS$Schedule)





########################
# PLOT FS MEDIANS
########################

# Create Grob function for highlighting significant p-values later:
find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}

#-------------------------
# 5mg and QFT-
#-------------------------

# Select 5mg and QFT-:
dftest <- scores[scores[,"Dose"]== 5 &  scores[,"QFT"]=="QFT-",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(FS ~ Schedule))
# Adjusted P-values:
(padj <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
p <- ggplot(meds_FS[meds_FS[,"Dose"]== 5 &  meds_FS[,"QFT"]=="QFT-",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("FS") +
  xlab(NULL)+
  ylim(0,0.85)+
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









#-------------------------
# 5mg and QFT+
#-------------------------
dftest <- scores[scores[,"Dose"]== 5 &  scores[,"QFT"]=="QFT+",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(FS ~ Schedule))
# Adjusted P-values:
(padj <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
p <- ggplot(meds_FS[meds_FS[,"Dose"]== 5 & meds_FS[,"QFT"]=="QFT+",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("FS") +
  xlab(NULL)+
  ylim(0,0.85)+
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


#-------------------------
# 15mg and QFT-
#-------------------------
dftest <- scores[scores[,"Dose"]== 15 &  scores[,"QFT"]=="QFT-",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(FS ~ Schedule))
# Adjusted P-values:
(padj <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
p <- ggplot(meds_FS[meds_FS[,"Dose"]== 15 & meds_FS[,"QFT"]=="QFT-",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("FS") +
  xlab(NULL)+
  ylim(0,0.85)+
  scale_x_discrete(labels = c("", "", "")) +
  scale_color_manual(values=c("#FF9999","#FF0000"), name = "Admin.")
# Make grobtable:
df1         <- data.frame(Test = "1 vs 2", Baseline = padj[1], 
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


#-------------------------
# 15mg and QFT+
#-------------------------
dftest     <- scores[scores[,"Dose"]== 15 &  scores[,"QFT"]=="QFT+",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(FS ~ Schedule))
# Adjusted P-values:
(padj      <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot: 
p <- ggplot(meds_FS[meds_FS[,"Dose"]== 15 & meds_FS[,"QFT"]=="QFT+",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("FS") +
  xlab(NULL)+
  ylim(0,0.85)+
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

#-------------------------
# 50mg and QFT-
#-------------------------
dftest     <- scores[scores[,"Dose"]== 50 &  scores[,"QFT"]=="QFT-",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(FS ~ Schedule))
# Adjusted P-values:
(padj      <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot: 
p <- ggplot(meds_FS[meds_FS[,"Dose"]== 50 & meds_FS[,"QFT"]=="QFT-",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,  # Bottom margin
                             l = 10)) +
  ylab("FS") +
  xlab(NULL)+
  ylim(0,0.85)+
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
ind3 <- find_cell(tab, 3, 1, "core-fg")
ind4 <- find_cell(tab, 3, 2, "core-fg")
ind5 <- find_cell(tab, 3, 3, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind2][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind3][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind4][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind5][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
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


#-------------------------
# 50mg and QFT+
#-------------------------
dftest <- scores[scores[,"Dose"]== 50 &  scores[,"QFT"]=="QFT+",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(FS ~ Schedule))
# Adjusted P-values:
(padj      <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot: 
p <- ggplot(meds_FS[meds_FS[,"Dose"]== 50 & meds_FS[,"QFT"]=="QFT+",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("FS") +
  xlab(NULL)+
  ylim(0,0.85)+
  scale_x_discrete(labels = c("", "", "")) +
  scale_color_manual(values=c("#99CCFF", "blue", "#000066"), name = "Admin.")
# Make grobtable:
df1         <- data.frame(Test = c("1 vs 2", "1 vs 3", "2 vs 3"), Baseline = padj[1:3], 
                          Peak = padj[4:6], Memory = padj[7:9]) 
names(df1)  <- c("Test", "Baseline", "Peak", "Memory")
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 25))
# Set significant p-values to bold:
ind1 <- find_cell(tab, 3, 2, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
ind1 <- find_cell(tab, 3, 3, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
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






########################
# PLOT PFS MEDIANS
########################

#-------------------------
# 5mg and QFT-
#-------------------------

# Select 5mg and QFT-:
dftest <- scores[scores[,"Dose"]== 5 &  scores[,"QFT"]=="QFT-",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(PFS ~ Schedule))
# Adjusted P-values:
(padj <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
p <- ggplot(meds_PFS[meds_PFS[,"Dose"]== 5 &  meds_PFS[,"QFT"]=="QFT-",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("PFS") +
  xlab(NULL)+
  ylim(0,0.85)+
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









#-------------------------
# 5mg and QFT+
#-------------------------
dftest <- scores[scores[,"Dose"]== 5 &  scores[,"QFT"]=="QFT+",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(PFS ~ Schedule))
# Adjusted P-values:
(padj <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
p <- ggplot(meds_PFS[meds_PFS[,"Dose"]== 5 & meds_PFS[,"QFT"]=="QFT+",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("PFS") +
  ylim(0,0.85)+
  xlab(NULL)+
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


#-------------------------
# 15mg and QFT-
#-------------------------
dftest <- scores[scores[,"Dose"]== 15 &  scores[,"QFT"]=="QFT-",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(PFS ~ Schedule))
# Adjusted P-values:
(padj <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
p <- ggplot(meds_PFS[meds_PFS[,"Dose"]== 15 & meds_PFS[,"QFT"]=="QFT-",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("PFS") +
  xlab(NULL)+
  ylim(0,0.85)+
  scale_x_discrete(labels = c("", "", "")) +
  scale_color_manual(values=c("#FF9999","#FF0000"), name = "Admin.")
# Make grobtable:
df1         <- data.frame(Test = "1 vs 2", Baseline = padj[1], 
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


#-------------------------
# 15mg and QFT+
#-------------------------
dftest     <- scores[scores[,"Dose"]== 15 &  scores[,"QFT"]=="QFT+",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(PFS ~ Schedule))
# Adjusted P-values:
(padj      <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot: 
p <- ggplot(meds_PFS[meds_PFS[,"Dose"]== 15 & meds_PFS[,"QFT"]=="QFT+",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("PFS") +
  xlab(NULL)+
  ylim(0,0.85)+
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

#-------------------------
# 50mg and QFT-
#-------------------------
dftest     <- scores[scores[,"Dose"]== 50 &  scores[,"QFT"]=="QFT-",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(PFS ~ Schedule))
# Adjusted P-values:
(padj      <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot: 
p <- ggplot(meds_PFS[meds_PFS[,"Dose"]== 50 & meds_PFS[,"QFT"]=="QFT-",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,  # Bottom margin
                             l = 10)) +
  ylab("PFS") +
  xlab(NULL)+
  ylim(0,0.85)+
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
ind3 <- find_cell(tab, 3, 1, "core-fg")
ind4 <- find_cell(tab, 3, 2, "core-fg")
ind5 <- find_cell(tab, 3, 3, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind2][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind3][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind4][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind5][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
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


#-------------------------
# 50mg and QFT+
#-------------------------
dftest <- scores[scores[,"Dose"]== 50 &  scores[,"QFT"]=="QFT+",]
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(timepnt) %>%
    wilcox_test(PFS ~ Schedule))
# Adjusted P-values:
(padj      <- round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot: 
p <- ggplot(meds_PFS[meds_PFS[,"Dose"]== 50 & meds_PFS[,"QFT"]=="QFT+",], aes(y=Medians, x=timepnt, group = Schedule)) +
  geom_line(aes(color=Schedule),linewidth=2)+
  geom_point(size=5, aes(color = Schedule))+
  geom_errorbar(aes(ymin=LowerQ, ymax=UpperQ,color=Schedule), width=.2, position=position_dodge(0.05), linewidth=0.8) +
  theme_minimal()+
  theme(text = element_text(size = 30), 
        plot.margin = margin(t = 20,  # Top margin
                             r = 10,  # Right margin
                             b = 0,   # Bottom margin
                             l = 10)) +
  ylab("PFS") +
  xlab(NULL)+
  ylim(0,0.85)+
  scale_x_discrete(labels = c("", "", "")) +
  scale_color_manual(values=c("#99CCFF", "blue", "#000066"), name = "Admin.")
# Make grobtable:
df1         <- data.frame(Test = c("1 vs 2", "1 vs 3", "2 vs 3"), Baseline = padj[1:3], 
                          Peak = padj[4:6], Memory = padj[7:9]) 
names(df1)  <- c("Test", "Baseline", "Peak", "Memory")
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 25))
ind1 <- find_cell(tab, 3, 2, "core-fg")
ind2 <- find_cell(tab, 3, 3, "core-fg")
ind3 <- find_cell(tab, 4, 1, "core-fg")
ind4 <- find_cell(tab, 4, 2, "core-fg")
ind5 <- find_cell(tab, 4, 3, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind2][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind3][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind4][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
tab$grobs[ind5][[1]][["gp"]] <- gpar(fontsize=25, fontface="bold")
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






############################
# PLOTS FOR 2 ADMIN
############################
scores2 = scores[scores[,"Schedule"]== 2,]
#-------------------------
# FS Baseline
#-------------------------
dftest     <- scores2[scores2[,"timepnt"]== "baseline",]
dftest$Dose = as.factor(dftest$Dose)
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(Dose) %>%
    wilcox_test(FS ~ QFT))
stat.test <- stat.test %>%
  add_xy_position(x = "Dose", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
ggplot(dftest, aes(y=FS, x=Dose)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("FS") +
  xlab(expression(paste("Concentration (", mu, "g)")))+
  ylim(0,1)+
  scale_fill_manual(values=c("#FF0000","blue"), name="")+
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8)





#-------------------------
# FS Peak
#-------------------------
dftest     <- scores2[scores2[,"timepnt"]== "peak",]
dftest$Dose = as.factor(dftest$Dose)
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(Dose) %>%
    wilcox_test(FS ~ QFT))
stat.test <- stat.test %>%
  add_xy_position(x = "Dose", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
ggplot(dftest, aes(y=FS, x=Dose)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("FS") +
  xlab(expression(paste("Concentration (", mu, "g)")))+
  ylim(0,1)+
  scale_fill_manual(values=c("#FF0000","blue"), name="")+
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8) 


#-------------------------
# FS Memory
#-------------------------
dftest     <- scores2[scores2[,"timepnt"]== "memory",]
dftest$Dose = as.factor(dftest$Dose)
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(Dose) %>%
    wilcox_test(FS ~ QFT))
stat.test <- stat.test %>%
  add_xy_position(x = "Dose", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
ggplot(dftest, aes(y=FS, x=Dose)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("FS") +
  xlab(expression(paste("Concentration (", mu, "g)")))+
  ylim(0,1)+
  scale_fill_manual(values=c("#FF0000","blue"), name="")+
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8) 


#-------------------------
# PFS Baseline
#-------------------------
dftest     <- scores2[scores2[,"timepnt"]== "baseline",]
dftest$Dose = as.factor(dftest$Dose)
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(Dose) %>%
    wilcox_test(PFS ~ QFT))
stat.test <- stat.test %>%
  add_xy_position(x = "Dose", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
ggplot(dftest, aes(y=PFS, x=Dose)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("PFS") +
  xlab(expression(paste("Concentration (", mu, "g)")))+
  ylim(0,1)+
  scale_fill_manual(values=c("#FF0000","blue"), name="")+
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8) 





#-------------------------
# PFS Peak
#-------------------------
dftest     <- scores2[scores2[,"timepnt"]== "peak",]
dftest$Dose = as.factor(dftest$Dose)
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(Dose) %>%
    wilcox_test(PFS ~ QFT))
stat.test <- stat.test %>%
  add_xy_position(x = "Dose", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
ggplot(dftest, aes(y=PFS, x=Dose)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("PFS") +
  xlab(expression(paste("Concentration (", mu, "g)")))+
  ylim(0,1)+
  scale_fill_manual(values=c("#FF0000","blue"), name="")+
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8) 


#-------------------------
# PFS Memory
#-------------------------
dftest     <- scores2[scores2[,"timepnt"]== "memory",]
dftest$Dose = as.factor(dftest$Dose)
# Calculate p-values:
(stat.test <- dftest %>%
    group_by(Dose) %>%
    wilcox_test(PFS ~ QFT))
stat.test <- stat.test %>%
  add_xy_position(x = "Dose", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
# Create Plot:
ggplot(dftest, aes(y=PFS, x=Dose)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("PFS") +
  xlab(expression(paste("Concentration (", mu, "g)")))+
  ylim(0,1)+
  scale_fill_manual(values=c("#FF0000","blue"), name="")+
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8) 






############################
# PLOTS CHOOSING CONC
############################
# again only use 2 admin data
#-------------------------
# FS Peak
#-------------------------
dftest     <- scores2[scores2[,"timepnt"]== "peak",]
dftest$Dose = as.factor(dftest$Dose)
# Calculate p-values:
stat.test  <- dftest %>%
  group_by(QFT) %>%
  wilcox_test(FS ~ Dose) 
(stat.test <- stat.test %>%
    add_xy_position(x = "QFT", dodge = 0.8))
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$y.position = c(0.81, 0.89, 0.83, 0.8, 0.88, 0.82)
# Create Plot:
ggplot(dftest, aes(y=FS, x=QFT)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Dose)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("FS") +
  xlab(NULL)+
  ylim(0,0.91)+
  scale_fill_manual(values = brewer.pal(n = 3, name = "BuGn"),name="Conc.", 
                    labels=c(expression(paste("5", mu, "g")), 
                             expression(paste("15", mu, "g")), 
                             expression(paste("50", mu, "g"))))+ 
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8) 


#-------------------------
# FS Memory
#-------------------------
dftest     <- scores2[scores2[,"timepnt"]== "memory",]
dftest$Dose = as.factor(dftest$Dose)
# Calculate p-values:
stat.test  <- dftest %>%
  group_by(QFT) %>%
  wilcox_test(FS ~ Dose) 
(stat.test <- stat.test %>%
    add_xy_position(x = "QFT", dodge = 0.8))
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$y.position = c(0.71, 0.79, 0.73, 0.74, 0.86, 0.8)
# Create Plot:
ggplot(dftest, aes(y=FS, x=QFT)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Dose)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("FS") +
  xlab(NULL)+
  ylim(0,0.91)+
  scale_fill_manual(values = brewer.pal(n = 3, name = "BuGn"),name="Conc.", 
                    labels=c(expression(paste("5", mu, "g")), 
                             expression(paste("15", mu, "g")), 
                             expression(paste("50", mu, "g"))))+ 
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8) 


#-------------------------
# PFS Peak
#-------------------------
dftest     <- scores2[scores2[,"timepnt"]== "peak",]
dftest$Dose = as.factor(dftest$Dose)
# Calculate p-values:
stat.test  <- dftest %>%
  group_by(QFT) %>%
  wilcox_test(PFS ~ Dose) 
(stat.test <- stat.test %>%
    add_xy_position(x = "QFT", dodge = 0.8))
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$y.position = c(0.91, 1, 0.93, 0.92, 1.01, 0.94)
# Create Plot:
ggplot(dftest, aes(y=PFS, x=QFT)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Dose)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("PFS") +
  xlab(NULL)+
  ylim(0,1.1)+
  scale_fill_manual(values = brewer.pal(n = 3, name = "BuGn"),name="Conc.", 
                    labels=c(expression(paste("5", mu, "g")), 
                             expression(paste("15", mu, "g")), 
                             expression(paste("50", mu, "g"))))+ 
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8) 



#-------------------------
# PFS Memory
#-------------------------
dftest     <- scores2[scores2[,"timepnt"]== "memory",]
dftest$Dose = as.factor(dftest$Dose)
# Calculate p-values:
stat.test  <- dftest %>%
  group_by(QFT) %>%
  wilcox_test(PFS ~ Dose) 
(stat.test <- stat.test %>%
    add_xy_position(x = "QFT", dodge = 0.8))
# Calculate adjusted p-values:
(stat.test$p.adj <-round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$y.position = c(0.81, 0.91, 0.84, 0.87, 0.97, 0.9)
# Create Plot:
ggplot(dftest, aes(y=PFS, x=QFT)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Dose)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("PFS") +
  xlab(NULL)+
  ylim(0,1.1)+
  scale_fill_manual(values = brewer.pal(n = 3, name = "BuGn"),name="Conc.", 
                    labels=c(expression(paste("5", mu, "g")), 
                             expression(paste("15", mu, "g")), 
                             expression(paste("50", mu, "g"))))+ 
  scale_x_discrete(labels = c("QFT-", "QFT+")) +
  stat_pvalue_manual(stat.test, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8) 

