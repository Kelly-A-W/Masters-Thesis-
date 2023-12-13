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
library(patchwork)
library(Rmisc) # to get CI
library(gridExtra)

# Get rid of scientific notation:
options(scipen=999) 

#---------------------
# Import Data
#---------------------

cd4_counts <- read.csv("Filtered_cd4_counts_dataset_1.csv", check.names = F)
str(cd4_counts)


# Set Factors:
cd4_counts$timepnt <- as.factor(cd4_counts$timepnt)
cd4_counts$timepnt <- factor(cd4_counts$timepnt, levels=c("baseline", "peak", "memory"))
levels(cd4_counts$timepnt)
cd4_counts$Dose <- as.factor(cd4_counts$Dose)
cd4_counts$Schedule <- as.factor(cd4_counts$Schedule)


#---------------------------------------
# Frequency of cytockine+ cells
#---------------------------------------

# Keep only ESAT6 and Ag85B:
cd4_counts <- cd4_counts[cd4_counts[,"Stim"] == "ESAT6" | cd4_counts[,"Stim"] == "Ag85B",]

# Only keep counts for cells expressing at least one of INFg, IL2, TNF or L17
cd4_counts <- cd4_counts[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd4_counts), fixed = TRUE)==T)]
str(cd4_counts)

# Sum over cytokines:
cd4_counts$R7..RA..   <- rowSums(cd4_counts[,11:25])
cd4_counts$R7..RAc..  <- rowSums(cd4_counts[,26:40])
cd4_counts$R7c..RA..  <- rowSums(cd4_counts[,41:55])
cd4_counts$R7c..RAc.. <- rowSums(cd4_counts[,56:70])

# Remove unnecessary columns:
cd4_counts = cd4_counts[,-c(4,6,8,10:70)]

# Convert to Frequencies (%):
cd4_counts[,7:10] <- (cd4_counts[,7:10]/rowSums(cd4_counts[,7:10]))*100
# check:
sum(is.na(cd4_counts$R7..RA..))
sum(is.na(cd4_counts$R7..RAc..))
sum(is.na(cd4_counts$R7c..RA..))
sum(is.na(cd4_counts$R7c..RAc..)) # All 0, so no NAs!


#---------------------
# Reshape Data
#---------------------

dat_long <- cd4_counts %>% pivot_longer(cols='R7..RA..':'R7c..RAc..',
                                        names_to='boolean',
                                        values_to='frequency')


# Make boolean labels nice
dat_long$boolean <- gsub("[[:punct:]]", "+", dat_long$boolean)
dat_long$boolean <- gsub("\\++", "+", dat_long$boolean)
dat_long$boolean <- gsub("c\\+", "-", dat_long$boolean)

# Convert to factors
dat_long$boolean <- as.factor(dat_long$boolean)


# Reorder levels of boolean factor 
dat_long$boolean <- factor(dat_long$boolean, 
                           levels=c("R7+RA+", "R7+RA-", "R7-RA-","R7-RA+"))
levels(dat_long$boolean)


###############################
# COMPARING ANTIGENS:
###############################

# only 2 admin of 5mg
df_ca = dat_long[dat_long[,"Schedule"]==2 & dat_long[,"Dose"]==5,]

######################## 
# Baseline QFT-:
######################## 

# Create Grob function for highlighting significant p-values later:
find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}


dftest    <- df_ca[df_ca[,"timepnt"]=="baseline" & df_ca[,"QFT"]=="QFT neg",]

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=Stim),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = Stim, group=Stim), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#669999", "#996699"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 95)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7"), 
                          x1 = c("+", "+"),
                          x2 = c("-", "+"),
                          x3 = c("-", "-"),
                          x4 = c("-", "+")) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 30), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=18))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -1.2, 3))




######################## 
# Peak QFT-:
######################## 

dftest    <- df_ca[df_ca[,"timepnt"]=="peak" & df_ca[,"QFT"]=="QFT neg",]

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=Stim),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = Stim, group=Stim), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#669999", "#996699"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 95)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7"), 
                          x1 = c("+", "+"),
                          x2 = c("-", "+"),
                          x3 = c("-", "-"),
                          x4 = c("-", "+")) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 30), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=18))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -1.2, 3))



######################## 
# Memory QFT-:
######################## 

dftest    <- df_ca[df_ca[,"timepnt"]=="memory" & df_ca[,"QFT"]=="QFT neg",]


# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=Stim),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = Stim, group=Stim), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#669999", "#996699"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 95)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7"), 
                          x1 = c("+", "+"),
                          x2 = c("-", "+"),
                          x3 = c("-", "-"),
                          x4 = c("-", "+")) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 30), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=18))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -1.2, 3))



######################## 
# Baseline QFT+:
######################## 

dftest    <- df_ca[df_ca[,"timepnt"]=="baseline" & df_ca[,"QFT"]=="QFT pos",]

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=Stim),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = Stim, group=Stim), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#669999", "#996699"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 95)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7"), 
                          x1 = c("+", "+"),
                          x2 = c("-", "+"),
                          x3 = c("-", "-"),
                          x4 = c("-", "+")) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 30), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=18))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -1.2, 3))





######################## 
# Peak QFT+:
######################## 
dftest    <- df_ca[df_ca[,"timepnt"]=="peak" & df_ca[,"QFT"]=="QFT pos",]

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=Stim),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = Stim, group=Stim), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#669999", "#996699"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 95)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7"), 
                          x1 = c("+", "+"),
                          x2 = c("-", "+"),
                          x3 = c("-", "-"),
                          x4 = c("-", "+")) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 30), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=18))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -1.2, 3))


######################## 
# Memory QFT+:
######################## 
dftest    <- df_ca[df_ca[,"timepnt"]=="memory" & df_ca[,"QFT"]=="QFT pos",]

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=Stim),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = Stim, group=Stim), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#669999", "#996699"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 95)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7"), 
                          x1 = c("+", "+"),
                          x2 = c("-", "+"),
                          x3 = c("-", "-"),
                          x4 = c("-", "+")) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 30), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=18))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -1.2, 3))







#########################################
# COMPARING QFT:
#########################################

# only 2 admin vaccine:
df_cq = dat_long[dat_long[,"Schedule"]==2,]


######################## 
# Memory Ag85B 5mg
######################## 
# Conduct Mann Whitney Test
dftest    <- df_cq[df_cq[,"timepnt"]=="memory" & df_cq[,"Stim"]=="Ag85B" & df_cq[,"Dose"]==5,]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))


# Label positions:
pos <- aggregate(dftest$frequency, by = list(dftest$boolean), max)$x

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=QFT),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = QFT, group=QFT), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 80)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7", "P-value"), 
                          x1 = c("+", "+", stat.test$padj[1]),
                          x2 = c("-", "+", stat.test$padj[2]),
                          x3 = c("-", "-", stat.test$padj[3]),
                          x4 = c("-", "+", stat.test$padj[4])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 20), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", "CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))








######################## 
# Memory Ag85B 15mg
######################## 
# Conduct Mann Whitney Test
dftest    <- df_cq[df_cq[,"timepnt"]=="memory" & df_cq[,"Stim"]=="Ag85B" & df_cq[,"Dose"]==15,]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))


# Label positions:
pos <- aggregate(dftest$frequency, by = list(dftest$boolean), max)$x

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=QFT),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = QFT, group=QFT), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 80)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7", "P-value"), 
                          x1 = c("+", "+", stat.test$padj[1]),
                          x2 = c("-", "+", stat.test$padj[2]),
                          x3 = c("-", "-", stat.test$padj[3]),
                          x4 = c("-", "+", stat.test$padj[4])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 20), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", "CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))








######################## 
# Memory Ag85B 50mg
######################## 
# Conduct Mann Whitney Test
dftest    <- df_cq[df_cq[,"timepnt"]=="memory" & df_cq[,"Stim"]=="Ag85B" & df_cq[,"Dose"]==50,]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))


# Label positions:
pos <- aggregate(dftest$frequency, by = list(dftest$boolean), max)$x

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=QFT),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = QFT, group=QFT), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 80)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7", "P-value"), 
                          x1 = c("+", "+", stat.test$padj[1]),
                          x2 = c("-", "+", stat.test$padj[2]),
                          x3 = c("-", "-", stat.test$padj[3]),
                          x4 = c("-", "+", stat.test$padj[4])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 20), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", "CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))








######################## 
# Memory ESAT6 5mg
######################## 
# Conduct Mann Whitney Test
dftest    <- df_cq[df_cq[,"timepnt"]=="memory" & df_cq[,"Stim"]=="ESAT6" & df_cq[,"Dose"]==5,]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))


# Label positions:
pos <- aggregate(dftest$frequency, by = list(dftest$boolean), max)$x

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=QFT),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = QFT, group=QFT), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 80)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7", "P-value"), 
                          x1 = c("+", "+", stat.test$padj[1]),
                          x2 = c("-", "+", stat.test$padj[2]),
                          x3 = c("-", "-", stat.test$padj[3]),
                          x4 = c("-", "+", stat.test$padj[4])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 20), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", "CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))








######################## 
# Memory ESAT6 15mg
######################## 
# Conduct Mann Whitney Test
dftest    <- df_cq[df_cq[,"timepnt"]=="memory" & df_cq[,"Stim"]=="ESAT6" & df_cq[,"Dose"]==15,]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))


# Label positions:
pos <- aggregate(dftest$frequency, by = list(dftest$boolean), max)$x

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=QFT),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = QFT, group=QFT), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 80)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7", "P-value"), 
                          x1 = c("+", "+", stat.test$padj[1]),
                          x2 = c("-", "+", stat.test$padj[2]),
                          x3 = c("-", "-", stat.test$padj[3]),
                          x4 = c("-", "+", stat.test$padj[4])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 20), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", "CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))








######################## 
# Memory ESAT6 50mg
######################## 
# Conduct Mann Whitney Test
dftest    <- df_cq[df_cq[,"timepnt"]=="memory" & df_cq[,"Stim"]=="ESAT6" & df_cq[,"Dose"]==50,]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))


# Label positions:
pos <- aggregate(dftest$frequency, by = list(dftest$boolean), max)$x

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=QFT),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = QFT, group=QFT), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 80)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7", "P-value"), 
                          x1 = c("+", "+", stat.test$padj[1]),
                          x2 = c("-", "+", stat.test$padj[2]),
                          x3 = c("-", "-", stat.test$padj[3]),
                          x4 = c("-", "+", stat.test$padj[4])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 20), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", "CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))








############################################
# COMPARING QFT: 2 ADMIN OF 5mg ONLY
############################################

# only 2 admin and 5mg vaccine:
df_cq25 = dat_long[dat_long[,"Schedule"]==2 & dat_long[,"Dose"]==5,]

######################## 
# Peak Ag85B
######################## 
# Conduct Mann Whitney Test
dftest    <- df_cq25[df_cq25[,"timepnt"]=="peak" & df_cq25[,"Stim"]=="Ag85B",]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))


# Label positions:
pos <- aggregate(dftest$frequency, by = list(dftest$boolean), max)$x

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=QFT),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = QFT, group=QFT), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 80)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7", "P-value"), 
                          x1 = c("+", "+", stat.test$padj[1]),
                          x2 = c("-", "+", stat.test$padj[2]),
                          x3 = c("-", "-", stat.test$padj[3]),
                          x4 = c("-", "+", stat.test$padj[4])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 20), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", "CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))









######################## 
# Peak ESAT6
######################## 
# Conduct Mann Whitney Test
df_cq25    <- df_cq25[df_cq25[,"timepnt"]=="peak" & df_cq25[,"Stim"]=="ESAT6",]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))


# Label positions:
pos <- aggregate(dftest$frequency, by = list(dftest$boolean), max)$x

# plot
p <- ggplot(dftest, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, fill = NA, aes(color=QFT),outlier.shape = NA)+
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  geom_point(aes(color = QFT, group=QFT), position=position_jitterdodge(), size=3, alpha=0.6)+    
  scale_color_manual(values = c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")+
  ylab("% of cytokine+ cells")+
  xlab(NULL)+
  ylim(0, 80)+
  scale_x_discrete(labels = rep("", 4))
# Make grobtable:
df1         <- data.frame(Marker = c("RA", "R7", "P-value"), 
                          x1 = c("+", "+", stat.test$padj[1]),
                          x2 = c("-", "+", stat.test$padj[2]),
                          x3 = c("-", "-", stat.test$padj[3]),
                          x4 = c("-", "+", stat.test$padj[4])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 20), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", "CCR7", "CD45RA"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))

