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
library(data.table) # for melt function

# Get rid of scientific notation:
options(scipen=999) 

#---------------------
# Import Data
#---------------------

cd4_counts <- read.csv("clean_data_with_demographic_info.csv", check.names = F)
str(cd4_counts)

# Remove demographic variables:
cd4_counts = cd4_counts[,-c(1,2,3)]

# Set Factors:
cd4_counts$timepnt <- as.factor(cd4_counts$timepnt)
cd4_counts$timepnt <- factor(cd4_counts$timepnt, levels=c("baseline", "peak", "memory"))
levels(cd4_counts$timepnt)
cd4_counts$Dose <- as.factor(cd4_counts$Dose)
cd4_counts$Schedule <- as.factor(cd4_counts$Schedule)


#---------------------
# BS Frequency of CD4
#---------------------
# Only keep counts for cells expressing at least one of INFg, IL2, TNF or L17 (i.e. cytokine+)
cd4_counts <- cd4_counts[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd4_counts), fixed = TRUE)==T)]
str(cd4_counts)


# Sum over memory markers:
cyto_counts = data.frame(cd4_counts[,1:10],
                         cd4_counts[,11:25] + cd4_counts[,26:40] + cd4_counts[,41:55] + cd4_counts[,56:70])
names(cyto_counts)[11:25] = gsub("R7..RA..", "", names(cyto_counts)[11:25])


# Convert to Frequencies (%):
cyto_counts[11:25] <- (cyto_counts[11:25]/cyto_counts$CD4count)*100


# Background Subtract
cd4_bs <- data.frame()
for(i in 1:length(unique(cyto_counts$PID))){
  pid <- unique(cyto_counts$PID)[i]
  dat <- cyto_counts[which(cyto_counts$PID == pid),]
  days <- unique(dat$Day)
  for (j in 1:length(days)) {
    day <- days[j]
    df  <- dat[which(dat$Day == day),]
    df[-which(df$Stim == "UNS"),11:25]  <- df[-which(df$Stim == "UNS"),11:25] - df[rep(which(df$Stim == "UNS"), (length(df$Stim)-1)),11:25]
    cd4_bs <- rbind(cd4_bs,df)
  }
}
summary(cd4_bs)
which(is.na(cd4_bs)==T)
str(cd4_bs)

# Deal with Negative Frequencies:
for (i in 11:ncol(cd4_bs)) {
  # select one column of marker combinations
  column <- cd4_bs[,i]
  # extract negative values
  negval <- column[which(column<0)]
  
  # only proceed if there were negative values found
  if(length(negval)>1){ # have to have more than one value to get a CI
    # find upper and lower bounds of negative value
    ci  <- CI(negval,ci = 0.8)
    # calculate threshold
    cut <- abs(ci["upper"]-ci["lower"])
    # find indices of the values below the cut-off
    ind <- which(column<cut)
    # set all values less than the cut-off to 0 (including negative numbers)
    cd4_bs[ind,i] <- 0 }
  
  if(length(negval)==1){ # If only 1 value, abandon CI and set neg value to 0
    # find index of the value below the 0
    ind <- which(column<0)
    # set value to 0 
    cd4_bs[ind,i] <- 0 }
  
}
# Check:
sapply(sign(cd4_bs[,-(1:10)]), table)


# Only keep Ag85B and ESAT6
cd4_bs <- cd4_bs[cd4_bs[,"Stim"] == "ESAT6" | cd4_bs[,"Stim"] == "Ag85B",]


#---------------------
# Reshape Data
#---------------------

# Three-Level Grouping:
long <- melt(setDT(cd4_bs), id.vars = c("Schedule","Dose", "timepnt", "Study", "QFT", "Group", "PID", "Day", "Stim", "CD4count"),
             variable.name = "Marker")
long$Marker <- as.factor(long$Marker)
long = as.data.frame(long)
head(long)
dim(long)



# Make Marker labels nice
long$Marker <- gsub("[[:punct:]]", "+", long$Marker)
long$Marker <- gsub("\\++", "+", long$Marker)
long$Marker <- gsub("c\\+", "-", long$Marker)
long$Marker <- gsub("TNF", "TNF+", long$Marker)
long$Marker <- gsub("\\+c", "-", long$Marker)

# Convert to factors
long$Marker <- as.factor(long$Marker)



long$Marker <- factor(long$Marker,
                      levels=c("G-L2-L17+TNF-", "G-L2+L17-TNF-", "G-L2+L17+TNF-", 
                      "G-L2-L17-TNF+", "G-L2-L17+TNF+", "G-L2+L17-TNF+",
                      "G-L2+L17+TNF+", "G+L2+L17-TNF+", "G+L2+L17+TNF+", 
                      "G+L2+L17-TNF-", "G+L2+L17+TNF-", "G+L2-L17-TNF+", 
                      "G+L2-L17+TNF+",  "G+L2-L17-TNF-", "G+L2-L17+TNF-"))
levels(long$Marker)

long$QFT = as.factor(long$QFT)

###############################
# COMPARING ANTIGENS:
###############################

# only 2 admin of 5mg
df_ca = long[long[,"Schedule"]==2 & long[,"Dose"]==5,]

######################## 
# Baseline QFT-:
######################## 

# Create Grob function for highlighting significant p-values later:
find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}

# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="baseline" & df_ca[,"QFT"]=="QFT neg",]
(stat.test <- dftest %>%
    group_by(Marker) %>%
    wilcox_test(value ~ Stim, paired = T))
stat.test <- stat.test %>%
  add_xy_position(x = "Marker", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))

# Log Transform
df <- dftest
df$value <- log(df$value+0.00001)

# Label positions:
pos <- aggregate(df$value, by = list(df$Marker), max)$x

# plot
p <- ggplot(df, aes(y=value, x=Marker)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Stim)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12.5, -0.9)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values = c("#669999", "#996699"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("IL17", "TNF", "IL2", "G", "P-value"), 
                          x1 = c("+","-",  "-", "-", stat.test$padj[1]),
                          x2 = c( "-", "-","+", "-", stat.test$padj[2]),
                          x3 = c( "+","-", "+", "-", stat.test$padj[3]),
                          x4 = c("-", "+", "-", "-", stat.test$padj[4]),
                          x5 = c("+", "+", "-", "-", stat.test$padj[5]),
                          x6 = c("-", "+", "+", "-", stat.test$padj[6]),
                          x7 = c("+", "+", "+", "-", stat.test$padj[7]),
                          x8 = c("-", "+", "+", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "+", stat.test$padj[9]),
                          x10 = c("-", "-", "+", "+", stat.test$padj[10]),
                          x11 = c("+", "-", "+", "+", stat.test$padj[11]),
                          x12 = c("-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "+", "-", "+", stat.test$padj[13]),
                          x14 = c("-", "-", "-", "+", stat.test$padj[14]),
                          x15 = c("+", "-", "-", "+", stat.test$padj[15])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF", "IL17"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))



######################## 
# Peak QFT-:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="peak" & df_ca[,"QFT"]=="QFT neg",]
(stat.test <- dftest %>%
    group_by(Marker) %>%
    wilcox_test(value ~ Stim, paired = T))
stat.test <- stat.test %>%
  add_xy_position(x = "Marker", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$padj[6] = 0.0005
stat.test$padj[7] = 0.0005

# Log Transform
df <- dftest
df$value <- log(df$value+0.00001)

# Label positions:
pos <- aggregate(df$value, by = list(df$Marker), max)$x

# plot
p <- ggplot(df, aes(y=value, x=Marker)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Stim)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12.5, -0.9)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values = c("#669999", "#996699"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("IL17", "TNF", "IL2", "G", "P-value"), 
                          x1 = c("+","-",  "-", "-", stat.test$padj[1]),
                          x2 = c( "-", "-","+", "-", stat.test$padj[2]),
                          x3 = c( "+","-", "+", "-", stat.test$padj[3]),
                          x4 = c("-", "+", "-", "-", stat.test$padj[4]),
                          x5 = c("+", "+", "-", "-", stat.test$padj[5]),
                          x6 = c("-", "+", "+", "-", stat.test$padj[6]),
                          x7 = c("+", "+", "+", "-", stat.test$padj[7]),
                          x8 = c("-", "+", "+", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "+", stat.test$padj[9]),
                          x10 = c("-", "-", "+", "+", stat.test$padj[10]),
                          x11 = c("+", "-", "+", "+", stat.test$padj[11]),
                          x12 = c("-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "+", "-", "+", stat.test$padj[13]),
                          x14 = c("-", "-", "-", "+", stat.test$padj[14]),
                          x15 = c("+", "-", "-", "+", stat.test$padj[15])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set significant p-values to bold:
ind <- find_cell(tab, 5, 2, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 3, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 4, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 5, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 6, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=14, fontface="bold")
ind <- find_cell(tab, 5, 7, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=14, fontface="bold")
ind <- find_cell(tab, 5, 8, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 9, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 10, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 11, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 13, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF", "IL17"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))




######################## 
# Memory QFT-:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="memory" & df_ca[,"QFT"]=="QFT neg",]
(stat.test <- dftest %>%
    group_by(Marker) %>%
    wilcox_test(value ~ Stim, paired = T))
stat.test <- stat.test %>%
  add_xy_position(x = "Marker", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))

# Log Transform
df <- dftest
df$value <- log(df$value+0.00001)

# Label positions:
pos <- aggregate(df$value, by = list(df$Marker), max)$x

# plot
p <- ggplot(df, aes(y=value, x=Marker)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Stim)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12.5, -0.9)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values = c("#669999", "#996699"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("IL17", "TNF", "IL2", "G", "P-value"), 
                          x1 = c("+","-",  "-", "-", stat.test$padj[1]),
                          x2 = c( "-", "-","+", "-", stat.test$padj[2]),
                          x3 = c( "+","-", "+", "-", stat.test$padj[3]),
                          x4 = c("-", "+", "-", "-", stat.test$padj[4]),
                          x5 = c("+", "+", "-", "-", stat.test$padj[5]),
                          x6 = c("-", "+", "+", "-", stat.test$padj[6]),
                          x7 = c("+", "+", "+", "-", stat.test$padj[7]),
                          x8 = c("-", "+", "+", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "+", stat.test$padj[9]),
                          x10 = c("-", "-", "+", "+", stat.test$padj[10]),
                          x11 = c("+", "-", "+", "+", stat.test$padj[11]),
                          x12 = c("-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "+", "-", "+", stat.test$padj[13]),
                          x14 = c("-", "-", "-", "+", stat.test$padj[14]),
                          x15 = c("+", "-", "-", "+", stat.test$padj[15])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set significant p-values to bold:
ind <- find_cell(tab, 5, 1, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 2, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 6, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 7, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 8, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 10, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF", "IL17"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))




######################## 
# Baseline QFT+:
######################## 

# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="baseline" & df_ca[,"QFT"]=="QFT pos",]
(stat.test <- dftest %>%
    group_by(Marker) %>%
    wilcox_test(value ~ Stim, paired = T))
stat.test <- stat.test %>%
  add_xy_position(x = "Marker", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))

# Log Transform
df <- dftest
df$value <- log(df$value+0.00001)

# Label positions:
pos <- aggregate(df$value, by = list(df$Marker), max)$x

# plot
p <- ggplot(df, aes(y=value, x=Marker)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Stim)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12.5, -0.9)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values = c("#669999", "#996699"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("IL17", "TNF", "IL2", "G", "P-value"), 
                          x1 = c("+","-",  "-", "-", stat.test$padj[1]),
                          x2 = c( "-", "-","+", "-", stat.test$padj[2]),
                          x3 = c( "+","-", "+", "-", stat.test$padj[3]),
                          x4 = c("-", "+", "-", "-", stat.test$padj[4]),
                          x5 = c("+", "+", "-", "-", stat.test$padj[5]),
                          x6 = c("-", "+", "+", "-", stat.test$padj[6]),
                          x7 = c("+", "+", "+", "-", stat.test$padj[7]),
                          x8 = c("-", "+", "+", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "+", stat.test$padj[9]),
                          x10 = c("-", "-", "+", "+", stat.test$padj[10]),
                          x11 = c("+", "-", "+", "+", stat.test$padj[11]),
                          x12 = c("-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "+", "-", "+", stat.test$padj[13]),
                          x14 = c("-", "-", "-", "+", stat.test$padj[14]),
                          x15 = c("+", "-", "-", "+", stat.test$padj[15]))  
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set significant p-values to bold:
ind <- find_cell(tab, 5, 6, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 12, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF", "IL17"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))



######################## 
# Peak QFT+:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="peak" & df_ca[,"QFT"]=="QFT pos",]
(stat.test <- dftest %>%
    group_by(Marker) %>%
    wilcox_test(value ~ Stim, paired = T))
stat.test <- stat.test %>%
  add_xy_position(x = "Marker", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))

# Log Transform
df <- dftest
df$value <- log(df$value+0.00001)

# Label positions:
pos <- aggregate(df$value, by = list(df$Marker), max)$x

# plot
p <- ggplot(df, aes(y=value, x=Marker)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Stim)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12.5, -0.9)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values = c("#669999", "#996699"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("IL17", "TNF", "IL2", "G", "P-value"), 
                          x1 = c("+","-",  "-", "-", stat.test$padj[1]),
                          x2 = c( "-", "-","+", "-", stat.test$padj[2]),
                          x3 = c( "+","-", "+", "-", stat.test$padj[3]),
                          x4 = c("-", "+", "-", "-", stat.test$padj[4]),
                          x5 = c("+", "+", "-", "-", stat.test$padj[5]),
                          x6 = c("-", "+", "+", "-", stat.test$padj[6]),
                          x7 = c("+", "+", "+", "-", stat.test$padj[7]),
                          x8 = c("-", "+", "+", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "+", stat.test$padj[9]),
                          x10 = c("-", "-", "+", "+", stat.test$padj[10]),
                          x11 = c("+", "-", "+", "+", stat.test$padj[11]),
                          x12 = c("-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "+", "-", "+", stat.test$padj[13]),
                          x14 = c("-", "-", "-", "+", stat.test$padj[14]),
                          x15 = c("+", "-", "-", "+", stat.test$padj[15])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set significant p-values to bold:
ind <- find_cell(tab, 5, 2, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 12, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 14, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF", "IL17"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))




######################## 
# Memory QFT+:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="memory" & df_ca[,"QFT"]=="QFT pos",]
(stat.test <- dftest %>%
    group_by(Marker) %>%
    wilcox_test(value ~ Stim, paired = T))
stat.test <- stat.test %>%
  add_xy_position(x = "Marker", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))

# Log Transform
df <- dftest
df$value <- log(df$value+0.00001)

# Label positions:
pos <- aggregate(df$value, by = list(df$Marker), max)$x

# plot
p <- ggplot(df, aes(y=value, x=Marker)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Stim)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12.5, -0.9)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values = c("#669999", "#996699"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("IL17", "TNF", "IL2", "G", "P-value"), 
                          x1 = c("+","-",  "-", "-", stat.test$padj[1]),
                          x2 = c( "-", "-","+", "-", stat.test$padj[2]),
                          x3 = c( "+","-", "+", "-", stat.test$padj[3]),
                          x4 = c("-", "+", "-", "-", stat.test$padj[4]),
                          x5 = c("+", "+", "-", "-", stat.test$padj[5]),
                          x6 = c("-", "+", "+", "-", stat.test$padj[6]),
                          x7 = c("+", "+", "+", "-", stat.test$padj[7]),
                          x8 = c("-", "+", "+", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "+", stat.test$padj[9]),
                          x10 = c("-", "-", "+", "+", stat.test$padj[10]),
                          x11 = c("+", "-", "+", "+", stat.test$padj[11]),
                          x12 = c("-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "+", "-", "+", stat.test$padj[13]),
                          x14 = c("-", "-", "-", "+", stat.test$padj[14]),
                          x15 = c("+", "-", "-", "+", stat.test$padj[15])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set significant p-values to bold:
ind <- find_cell(tab, 5, 10, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 12, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 14, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF", "IL17"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))







###############################
# COMPARING QFT:
###############################

######################## 
# Peak ESAT6:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="peak" & df_ca[,"Stim"]=="ESAT6",]
(stat.test <- dftest %>%
    group_by(Marker) %>%
    wilcox_test(value ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "Marker", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$padj[10] = 0.0002
stat.test$padj[14] = 0.0002

# Log Transform
df <- dftest
df$value <- log(df$value+0.00001)

# Label positions:
pos <- aggregate(df$value, by = list(df$Marker), max)$x

# plot
p <- ggplot(df, aes(y=value, x=Marker)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color=rep(c("#660000", "#000066"), 15)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12.5, -0.9)+
  scale_x_discrete(labels = rep("", 15))+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("IL17", "TNF", "IL2", "G", "P-value"), 
                          x1 = c("+","-",  "-", "-", stat.test$padj[1]),
                          x2 = c( "-", "-","+", "-", stat.test$padj[2]),
                          x3 = c( "+","-", "+", "-", stat.test$padj[3]),
                          x4 = c("-", "+", "-", "-", stat.test$padj[4]),
                          x5 = c("+", "+", "-", "-", stat.test$padj[5]),
                          x6 = c("-", "+", "+", "-", stat.test$padj[6]),
                          x7 = c("+", "+", "+", "-", stat.test$padj[7]),
                          x8 = c("-", "+", "+", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "+", stat.test$padj[9]),
                          x10 = c("-", "-", "+", "+", stat.test$padj[10]),
                          x11 = c("+", "-", "+", "+", stat.test$padj[11]),
                          x12 = c("-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "+", "-", "+", stat.test$padj[13]),
                          x14 = c("-", "-", "-", "+", stat.test$padj[14]),
                          x15 = c("+", "-", "-", "+", stat.test$padj[15])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set significant p-values to bold:
ind <- find_cell(tab, 5, 10, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=14, fontface="bold")
ind <- find_cell(tab, 5, 13, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 14, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=14, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF", "IL17"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))









######################## 
# Memory ESAT6:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="memory" & df_ca[,"Stim"]=="ESAT6",]
p1 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[1],], value ~ QFT, paired = F)$p
p2 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[2],], value ~ QFT, paired = F)$p
p3 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[3],], value ~ QFT, paired = F)$p
p4 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[4],], value ~ QFT, paired = F)$p
p5 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[5],], value ~ QFT, paired = F)$p
p6 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[6],], value ~ QFT, paired = F)$p
p7 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[7],], value ~ QFT, paired = F)$p
p8 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[8],], value ~ QFT, paired = F)$p
p9 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[9],], value ~ QFT, paired = F)$p
p10 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[10],], value ~ QFT, paired = F)$p
#p11 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[11],], value ~ QFT, paired = F)$p. #Can't get this p-value because too zero-inflated
p12 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[12],], value ~ QFT, paired = F)$p
p13 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[13],], value ~ QFT, paired = F)$p
p14 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[14],], value ~ QFT, paired = F)$p
p15 = wilcox_test(dftest[dftest[,"Marker"]==levels(dftest$Marker)[15],], value ~ QFT, paired = F)$p

# Adjust p-values:
(padj <- round(p.adjust(c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p12,p13,p14,p15), method = "BH"),3))
padj[8] = 0.0004

# Log Transform
df <- dftest
df$value <- log(df$value+0.00001)

# plot
p <- ggplot(df, aes(y=value, x=Marker)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color=rep(c("#660000", "#000066"), 15)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12.5, -0.9)+
  scale_x_discrete(labels = rep("", 15))+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("IL17", "TNF", "IL2", "G", "P-value"), 
                          x1 = c("+","-",  "-", "-", padj[1]),
                          x2 = c( "-", "-","+", "-", padj[2]),
                          x3 = c( "+","-", "+", "-", padj[3]),
                          x4 = c("-", "+", "-", "-", padj[4]),
                          x5 = c("+", "+", "-", "-",padj[5]),
                          x6 = c("-", "+", "+", "-", padj[6]),
                          x7 = c("+", "+", "+", "-", padj[7]),
                          x8 = c("-", "+", "+", "+", padj[8]),
                          x9 = c("+", "+", "+", "+", padj[9]),
                          x10 = c("-", "-", "+", "+", padj[10]),
                          x11 = c("+", "-", "+", "+", "N/A"),
                          x12 = c("-", "+", "-", "+", padj[11]),
                          x13 = c("+", "+", "-", "+", padj[12]),
                          x14 = c("-", "-", "-", "+", padj[13]),
                          x15 = c("+", "-", "-", "+", padj[14])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set significant p-values to bold:
ind <- find_cell(tab, 5, 6, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 8, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=16, fontface="bold")
ind <- find_cell(tab, 5, 10, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 12, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 5, 14, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF", "IL17"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))



######################## 
# Peak Ag85B:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="peak" & df_ca[,"Stim"]=="Ag85B",]
(stat.test <- dftest %>%
    group_by(Marker) %>%
    wilcox_test(value ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "Marker", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))

# Log Transform
df <- dftest
df$value <- log(df$value+0.00001)

# Label positions:
pos <- aggregate(df$value, by = list(df$Marker), max)$x

# plot
p <- ggplot(df, aes(y=value, x=Marker)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color=rep(c("#660000", "#000066"), 15)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12.5, -0.9)+
  scale_x_discrete(labels = rep("", 15))+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("IL17", "TNF", "IL2", "G", "P-value"), 
                          x1 = c("+","-",  "-", "-", stat.test$padj[1]),
                          x2 = c( "-", "-","+", "-", stat.test$padj[2]),
                          x3 = c( "+","-", "+", "-", stat.test$padj[3]),
                          x4 = c("-", "+", "-", "-", stat.test$padj[4]),
                          x5 = c("+", "+", "-", "-", stat.test$padj[5]),
                          x6 = c("-", "+", "+", "-", stat.test$padj[6]),
                          x7 = c("+", "+", "+", "-", stat.test$padj[7]),
                          x8 = c("-", "+", "+", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "+", stat.test$padj[9]),
                          x10 = c("-", "-", "+", "+", stat.test$padj[10]),
                          x11 = c("+", "-", "+", "+", stat.test$padj[11]),
                          x12 = c("-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "+", "-", "+", stat.test$padj[13]),
                          x14 = c("-", "-", "-", "+", stat.test$padj[14]),
                          x15 = c("+", "-", "-", "+", stat.test$padj[15])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set significant p-values to bold:
ind <- find_cell(tab, 5, 7, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF", "IL17"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))



######################## 
# Memory Ag85B:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="memory" & df_ca[,"Stim"]=="Ag85B",]
(stat.test <- dftest %>%
    group_by(Marker) %>%
    wilcox_test(value ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "Marker", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))

# Log Transform
df <- dftest
df$value <- log(df$value+0.00001)

# Label positions:
pos <- aggregate(df$value, by = list(df$Marker), max)$x

# plot
p <- ggplot(df, aes(y=value, x=Marker)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color=rep(c("#660000", "#000066"), 15)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12.5, -0.9)+
  scale_x_discrete(labels = rep("", 15))+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("IL17", "TNF", "IL2", "G", "P-value"), 
                          x1 = c("+","-",  "-", "-", stat.test$padj[1]),
                          x2 = c( "-", "-","+", "-", stat.test$padj[2]),
                          x3 = c( "+","-", "+", "-", stat.test$padj[3]),
                          x4 = c("-", "+", "-", "-", stat.test$padj[4]),
                          x5 = c("+", "+", "-", "-", stat.test$padj[5]),
                          x6 = c("-", "+", "+", "-", stat.test$padj[6]),
                          x7 = c("+", "+", "+", "-", stat.test$padj[7]),
                          x8 = c("-", "+", "+", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "+", stat.test$padj[9]),
                          x10 = c("-", "-", "+", "+", stat.test$padj[10]),
                          x11 = c("+", "-", "+", "+", stat.test$padj[11]),
                          x12 = c("-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "+", "-", "+", stat.test$padj[13]),
                          x14 = c("-", "-", "-", "+", stat.test$padj[14]),
                          x15 = c("+", "-", "-", "+", stat.test$padj[15])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF", "IL17"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))




###############################
# COMPARING CONCENTRATIONS:
###############################

# Create IL2+TNF+IFNg- and IL2+TNF+IFNg+ columns
df_cc1 = data.frame(cd4_bs[,c(1,2,3,5,7,9)], 
                    Gc..L2..TNF = cd4_bs$Gc..L2..L17..TNF + cd4_bs$Gc..L2..L17c..TNF)
df_cc2 = data.frame(cd4_bs[,c(1,2,3,5,7,9)], 
                    G..L2..TNF = cd4_bs$G..L2..L17..TNF + cd4_bs$G..L2..L17c..TNF)

# Keep only 2 admin:
df_cc1 = df_cc1[df_cc1[,"Schedule"]==2,]
df_cc2 = df_cc2[df_cc2[,"Schedule"]==2,]

# Only at Memory:
df_cc1 = df_cc1[df_cc1[,"timepnt"]=="memory",]
df_cc2 = df_cc2[df_cc2[,"timepnt"]=="memory",]

# set factors
df_cc1$QFT = as.factor(df_cc1$QFT)
df_cc2$QFT = as.factor(df_cc2$QFT)
df_cc1$Dose = as.factor(df_cc1$Dose)
df_cc2$Dose = as.factor(df_cc2$Dose)

######################## 
# ESAT6 IL2+TNF+IFNg-:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_cc1[df_cc1[,"Stim"]=="ESAT6",]
(stat.test <- dftest %>%
    group_by(Dose) %>%
    wilcox_test(Gc..L2..TNF ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "Dose", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$y.position = c(-3, -2.5,-1.5)

# Log Transform
df <- dftest
df$Gc..L2..TNF <- log(df$Gc..L2..TNF+0.00001)

# plot
ggplot(df, aes(y=Gc..L2..TNF, x=Dose)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("Logged % of CD4 cells")+
  xlab(expression(paste("Concentration (", mu, "g)")))+
  ylim(-12.5, 0)+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")+
  stat_pvalue_manual(stat.test, 
                     label = "padj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8,
                     fontface=c("bold",  "plain", "plain", "plain", "plain", "plain", "plain", "plain", "plain"))














######################## 
# ESAT6 IL2+TNF+IFNg+:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_cc2[df_cc2[,"Stim"]=="ESAT6",]
(stat.test <- dftest %>%
    group_by(Dose) %>%
    wilcox_test(G..L2..TNF ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "Dose", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$padj[1] = 0.00004
stat.test$padj[2] = 0.00004
stat.test$y.position = c(-1, -1,-1)

# Log Transform
df <- dftest
df$G..L2..TNF <- log(df$G..L2..TNF+0.00001)

# plot
ggplot(df, aes(y=G..L2..TNF, x=Dose)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("Logged % of CD4 cells")+
  xlab(expression(paste("Concentration (", mu, "g)")))+
  ylim(-12.5, 0)+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")+
  stat_pvalue_manual(stat.test, 
                     label = "padj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8,
                     fontface=c("bold",  "plain", "plain", "bold", "plain", "plain", "bold", "plain", "plain"))



######################## 
# Ag85B IL2+TNF+IFNg-:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_cc1[df_cc1[,"Stim"]=="Ag85B",]
(stat.test <- dftest %>%
    group_by(Dose) %>%
    wilcox_test(Gc..L2..TNF ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "Dose", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$y.position = c(-1.5, -1,-1.5)

# Log Transform
df <- dftest
df$Gc..L2..TNF <- log(df$Gc..L2..TNF+0.00001)

# plot
ggplot(df, aes(y=Gc..L2..TNF, x=Dose)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("Logged % of CD4 cells")+
  xlab(expression(paste("Concentration (", mu, "g)")))+
  ylim(-12.5, 0)+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")+
  stat_pvalue_manual(stat.test, 
                     label = "padj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8)




######################## 
# Ag85B IL2+TNF+IFNg+:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_cc2[df_cc2[,"Stim"]=="ESAT6",]
(stat.test <- dftest %>%
    group_by(Dose) %>%
    wilcox_test(G..L2..TNF ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "Dose", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$padj[1] = 0.00004
stat.test$padj[2] = 0.00004
stat.test$y.position = c(-1, -1,-1)

# Log Transform
df <- dftest
df$G..L2..TNF <- log(df$G..L2..TNF+0.00001)

# plot
ggplot(df, aes(y=G..L2..TNF, x=Dose)) +
  geom_boxplot(alpha=1, color="black", size =1, aes(fill = QFT)) +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("Logged % of CD4 cells")+
  xlab(expression(paste("Concentration (", mu, "g)")))+
  ylim(-12.5, 0)+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")+
  stat_pvalue_manual(stat.test, 
                     label = "padj", 
                     tip.length = 0.01,
                     size = 10, 
                     bracket.size = 0.8,
                     fontface=c("bold",  "plain", "plain", "bold", "plain", "plain", "bold", "plain", "plain"))


