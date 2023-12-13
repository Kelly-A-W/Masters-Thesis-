
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

# Get rid of scientific notation:
options(scipen=999) 

#---------------------
# Import Data
#---------------------

scores <- read.csv("compass_esat6_dataset_1.csv", check.names = F)
str(scores)


############################
# PLOTS ESAT6 VS AG85B
############################
# Import Ag85B:
scores_ag85b <- read.csv("compass_ag85b_dataset_1.csv", check.names = F)
str(scores_ag85b)

# combine datasets:
df1 = data.frame(scores, Stim = "ESAT6")
df2 = data.frame(scores_ag85b, Stim = "Ag85B")
all_scores = rbind(df1,df2)

all_scores$Stim = as.factor(all_scores$Stim)

#------------------------
# FS Peak
#------------------------
dftest     <- all_scores[all_scores[,"timepnt"]== "peak",]
# Calculate p-values:
(stat.test1 <- dftest %>%
    group_by(Stim) %>%
    wilcox_test(FS ~ QFT, paired = F))
stat.test1 <- stat.test1 %>%
  add_xy_position(x = "Stim", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test1$p.adj <-p.adjust(stat.test1$p, method = "BH"))
stat.test1$p.adj = c("0.09", "<0.00001")
# Calculate p-values:
(stat.test2 <- dftest %>%
    group_by(QFT) %>%
    wilcox_test(FS ~ Stim, paired=T))
stat.test2 <- stat.test2 %>%
  add_xy_position(x = "QFT", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test2$p.adj <-p.adjust(stat.test2$p, method = "BH"))
stat.test2$p.adj = c("<0.00001", "<0.00001")
stat.test2$xmin = as.numeric(stat.test1[1,12:13])
stat.test2$xmax = as.numeric(stat.test1[2,12:13])
stat.test2$y.position = c(0.9,0.98)

# Create Plot:
ggplot(dftest, aes(y=FS, x=Stim)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("FS") +
  xlab(NULL)+
  #ylim(0,0.8)+
  scale_fill_manual(values=c("#FF0000","blue"), name="")+
  stat_pvalue_manual(stat.test1, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8, 
                     fontface=c("plain",  "plain", "plain", "bold", "plain", "plain"))+
  stat_pvalue_manual(stat.test2, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8, 
                     fontface="bold") 


#------------------------
# FS Memory
#------------------------
dftest     <- all_scores[all_scores[,"timepnt"]== "memory",]
# Calculate p-values:
(stat.test1 <- dftest %>%
    group_by(Stim) %>%
    wilcox_test(FS ~ QFT, paired = F))
stat.test1 <- stat.test1 %>%
  add_xy_position(x = "Stim", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test1$p.adj <-p.adjust(stat.test1$p, method = "BH"))
stat.test1$p.adj = c("0.351", "<0.00001")
# Calculate p-values:
(stat.test2 <- dftest %>%
    group_by(QFT) %>%
    wilcox_test(FS ~ Stim, paired=T))
stat.test2 <- stat.test2 %>%
  add_xy_position(x = "QFT", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test2$p.adj <-p.adjust(stat.test2$p, method = "BH"))
stat.test2$p.adj = c("<0.00001", "<0.00001")
stat.test2$xmin = as.numeric(stat.test1[1,12:13])
stat.test2$xmax = as.numeric(stat.test1[2,12:13])
stat.test2$y.position = c(0.9,0.98)

# Create Plot:
ggplot(dftest, aes(y=FS, x=Stim)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("FS") +
  xlab(NULL)+
  #ylim(0,0.8)+
  scale_fill_manual(values=c("#FF0000","blue"), name="")+
  stat_pvalue_manual(stat.test1, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8, 
                     fontface=c("plain",  "plain", "plain", "bold", "plain", "plain"))+
  stat_pvalue_manual(stat.test2, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8, 
                     fontface="bold") 





#------------------------
# PFS Peak
#------------------------
dftest     <- all_scores[all_scores[,"timepnt"]== "peak",]
# Calculate p-values:
(stat.test1 <- dftest %>%
    group_by(Stim) %>%
    wilcox_test(PFS ~ QFT, paired = F))
stat.test1 <- stat.test1 %>%
  add_xy_position(x = "Stim", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test1$p.adj <-p.adjust(stat.test1$p, method = "BH"))
stat.test1$p.adj = c("0.18", "<0.00001")
# Calculate p-values:
(stat.test2 <- dftest %>%
    group_by(QFT) %>%
    wilcox_test(PFS ~ Stim, paired=T))
stat.test2 <- stat.test2 %>%
  add_xy_position(x = "QFT", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test2$p.adj <-p.adjust(stat.test2$p, method = "BH"))
stat.test2$p.adj = c("<0.00001", "<0.00001")
stat.test2$xmin = as.numeric(stat.test1[1,12:13])
stat.test2$xmax = as.numeric(stat.test1[2,12:13])
stat.test2$y.position = c(1,1.08)

# Create Plot:
ggplot(dftest, aes(y=PFS, x=Stim)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("PFS") +
  xlab(NULL)+
  #ylim(0,0.8)+
  scale_fill_manual(values=c("#FF0000","blue"), name="")+
  stat_pvalue_manual(stat.test1, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8, 
                     fontface=c("plain",  "plain", "plain", "bold", "plain", "plain"))+
  stat_pvalue_manual(stat.test2, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8, 
                     fontface="bold") 


#------------------------
# PFS Memory
#------------------------
dftest     <- all_scores[all_scores[,"timepnt"]== "memory",]
# Calculate p-values:
(stat.test1 <- dftest %>%
    group_by(Stim) %>%
    wilcox_test(PFS ~ QFT, paired = F))
stat.test1 <- stat.test1 %>%
  add_xy_position(x = "Stim", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test1$p.adj <-p.adjust(stat.test1$p, method = "BH"))
stat.test1$p.adj = c("0.684", "<0.00001")
# Calculate p-values:
(stat.test2 <- dftest %>%
    group_by(QFT) %>%
    wilcox_test(PFS ~ Stim, paired=T))
stat.test2 <- stat.test2 %>%
  add_xy_position(x = "QFT", dodge = 0.8)
# Calculate adjusted p-values:
(stat.test2$p.adj <-p.adjust(stat.test2$p, method = "BH"))
stat.test2$p.adj = c("<0.00001", "<0.00001")
stat.test2$xmin = as.numeric(stat.test1[1,12:13])
stat.test2$xmax = as.numeric(stat.test1[2,12:13])
stat.test2$y.position = c(0.92,1)

# Create Plot:
ggplot(dftest, aes(y=PFS, x=Stim)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("PFS") +
  xlab(NULL)+
  #ylim(0,0.8)+
  scale_fill_manual(values=c("#FF0000","blue"), name="")+
  stat_pvalue_manual(stat.test1, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8, 
                     fontface=c("plain",  "plain", "plain", "bold", "plain", "plain"))+
  stat_pvalue_manual(stat.test2, 
                     label = "p.adj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8, 
                     fontface="bold") 


