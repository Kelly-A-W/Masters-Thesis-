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

cd4_counts <- read.csv("clean_cd4_counts_dataset_1.csv", check.names = F)
str(cd4_counts)



#---------------------
# Add Time Category
#---------------------
cd4_counts <- data.frame(timepnt = NA, cd4_counts)

cd4_counts$timepnt[which(cd4_counts$Day==0)] <- "baseline"
cd4_counts$timepnt[which(cd4_counts$Day==14)] <- "peak"
cd4_counts$timepnt[which(cd4_counts$Day==70)] <- "peak"
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
cd4_counts <- data.frame(Dose = NA,cd4_counts)

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
# convert to factor
cd4_counts$Dose <- as.factor(cd4_counts$Dose)


#---------------------
# Add Schedule
#---------------------
# Need to add in dose and schedule information
cd4_counts <- data.frame(Schedule = NA,cd4_counts)

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
# convert to factor
cd4_counts$Schedule <- as.factor(cd4_counts$Schedule)
head(cd4_counts)

#---------------------------------------------------------------
# Extract only pids receiving 2 administrations of 5 micrograms 
#-----------------------------------------------------------------

cd4_counts <- cd4_counts[cd4_counts[,"Schedule"]==2 & cd4_counts[,"Dose"]==5, ]

# check
cd4_counts$Schedule
cd4_counts$Dose


#---------------------
# BS Frequency of CD4
#---------------------

# Convert to Frequencies
freq <- cd4_counts
freq[,11:74] <- (cd4_counts[,11:74]/cd4_counts[,10])*100

# Background Subtract
cd4_bs <- data.frame()
for(i in 1:length(unique(cd4_counts$PID))){
  pid <- unique(cd4_counts$PID)[i]
  dat <- freq[which(freq$PID == pid),]
  days <- unique(dat$Day)
  for (j in 1:length(days)) {
    day <- days[j]
    df  <- dat[which(dat$Day == day),]
    df[-which(df$Stim == "UNS"),11:74]  <- df[-which(df$Stim == "UNS"),11:74] - df[rep(which(df$Stim == "UNS"), (length(df$Stim)-1)),11:74]
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


# Separate by stimulus
bsfreq_ag85b <- cd4_bs[cd4_bs[,"Stim"]=="Ag85B", -c(1,2,4,6,7,8,9,10)]
bsfreq_esat6 <- cd4_bs[cd4_bs[,"Stim"]=="ESAT6",-c(1,2,4,6,7,8,9,10)]





#----------------------------
# Get Rid of Unwanted Combos
#----------------------------

# Remove booleans that express no cytokines (note that we are ignoring il17)
bsfreq_ag85b <- bsfreq_ag85b[,-c(which(names(bsfreq_ag85b) == "R7c..RAc..Gc..L2c..L17c..TNFc"),
                                 which(names(bsfreq_ag85b) == "R7..RAc..Gc..L2c..L17c..TNFc"),
                                 which(names(bsfreq_ag85b) == "R7c..RA..Gc..L2c..L17c..TNFc"),
                                 which(names(bsfreq_ag85b) == "R7..RA..Gc..L2c..L17c..TNFc"),
                                 which(names(bsfreq_ag85b) == "R7c..RAc..Gc..L2c..L17..TNFc"), # here the only cytokine expressed is il17 - but we are ignoring 1l17 so this counts as no cytokines being expressed
                                 which(names(bsfreq_ag85b) == "R7..RAc..Gc..L2c..L17..TNFc"),
                                 which(names(bsfreq_ag85b) == "R7c..RA..Gc..L2c..L17..TNFc"),
                                 which(names(bsfreq_ag85b) == "R7..RA..Gc..L2c..L17..TNFc"))]
bsfreq_esat6 <- bsfreq_esat6[,-c(which(names(bsfreq_esat6) == "R7c..RAc..Gc..L2c..L17c..TNFc"),
                                 which(names(bsfreq_esat6) == "R7..RAc..Gc..L2c..L17c..TNFc"),
                                 which(names(bsfreq_esat6) == "R7c..RA..Gc..L2c..L17c..TNFc"),
                                 which(names(bsfreq_esat6) == "R7..RA..Gc..L2c..L17c..TNFc"),
                                 which(names(bsfreq_esat6) == "R7c..RAc..Gc..L2c..L17..TNFc"),
                                 which(names(bsfreq_esat6) == "R7..RAc..Gc..L2c..L17..TNFc"),
                                 which(names(bsfreq_esat6) == "R7c..RA..Gc..L2c..L17..TNFc"),
                                 which(names(bsfreq_esat6) == "R7..RA..Gc..L2c..L17..TNFc"))]

# Sum over R7-RA+ and R7+RA+
bsfreq_ag85b$R7..RA.. <- rowSums(bsfreq_ag85b[,3:16])
bsfreq_ag85b$R7c..RA.. <- rowSums(bsfreq_ag85b[,31:44])
bsfreq_esat6$R7..RA.. <- rowSums(bsfreq_esat6[,3:16])
bsfreq_esat6$R7c..RA.. <- rowSums(bsfreq_esat6[,31:44])


# Need to ignore IL17
# lets do sum checks to see how we can sum over IL17
x <-cbind(names(bsfreq_ag85b)[which(grepl("..L17c..", names(bsfreq_ag85b), fixed = TRUE)==T)],
          names(bsfreq_ag85b)[which(grepl("..L17..", names(bsfreq_ag85b), fixed = TRUE)==T)])
y <- gsub("L17c", "L17", x)
which(y[,1] != y[,2]) # 0, therefore the only difference between the two columns in any given row is 
# whether IL17 is expressed or not

# Now sum over IL17
ag85b_dat <- matrix(NA, nrow=nrow(bsfreq_ag85b), ncol = nrow(x))
ag85b_dat <- as.data.frame(ag85b_dat)
for (i in 1:nrow(x)) {
  ag85b_dat[,i] <- bsfreq_ag85b[,x[i,1]] + bsfreq_ag85b[,x[i,2]]
  names(ag85b_dat)[i] <- gsub("L17c..", "", x[i,1])
}
ag85b_dat$timepnt <- bsfreq_ag85b$timepnt
ag85b_dat$QFT <- bsfreq_ag85b$QFT
ag85b_dat$R7..RA.. <- bsfreq_ag85b$R7..RA..
ag85b_dat$R7c..RA.. <- bsfreq_ag85b$R7c..RA..
str(ag85b_dat)

esat6_dat <- matrix(NA, nrow=nrow(bsfreq_esat6), ncol = nrow(x))
esat6_dat <- as.data.frame(esat6_dat)
for (i in 1:nrow(x)) {
  esat6_dat[,i] <- bsfreq_esat6[,x[i,1]] + bsfreq_esat6[,x[i,2]]
  names(esat6_dat)[i] <- gsub("L17c..", "", x[i,1])
}
esat6_dat$timepnt <- bsfreq_esat6$timepnt
esat6_dat$QFT <- bsfreq_esat6$QFT
esat6_dat$R7..RA.. <- bsfreq_esat6$R7..RA..
esat6_dat$R7c..RA.. <- bsfreq_esat6$R7c..RA..
str(esat6_dat)


# Keep only columns of interest
ag85b_dat_smol <- ag85b_dat[,c("R7..RA..", "R7..RAc..Gc..L2..TNFc", "R7c..RAc..Gc..L2..TNFc", 
                               "R7..RAc..Gc..L2c..TNF", "R7c..RAc..Gc..L2c..TNF", "R7..RAc..Gc..L2..TNF",  
                               "R7c..RAc..Gc..L2..TNF", "R7..RAc..G..L2..TNF", "R7c..RAc..G..L2..TNF", 
                               "R7..RAc..G..L2..TNFc", "R7c..RAc..G..L2..TNFc", "R7..RAc..G..L2c..TNF",   
                               "R7c..RAc..G..L2c..TNF",  "R7..RAc..G..L2c..TNFc", "R7c..RAc..G..L2c..TNFc", 
                               "R7c..RA..", "timepnt", "QFT")]

esat6_dat_smol <- esat6_dat[,c("R7..RA..", "R7..RAc..Gc..L2..TNFc", "R7c..RAc..Gc..L2..TNFc", 
                               "R7..RAc..Gc..L2c..TNF", "R7c..RAc..Gc..L2c..TNF", "R7..RAc..Gc..L2..TNF",  
                               "R7c..RAc..Gc..L2..TNF", "R7..RAc..G..L2..TNF", "R7c..RAc..G..L2..TNF", 
                               "R7..RAc..G..L2..TNFc", "R7c..RAc..G..L2..TNFc", "R7..RAc..G..L2c..TNF",   
                               "R7c..RAc..G..L2c..TNF",  "R7..RAc..G..L2c..TNFc", "R7c..RAc..G..L2c..TNFc", 
                               "R7c..RA..", "timepnt", "QFT")]







#---------------------
# Reshape Data
#---------------------



#---------------------------------------------- Ag85B

ag85b_dat_long <- ag85b_dat_smol %>% pivot_longer(cols='R7..RA..':'R7c..RA..',
                                                  names_to='boolean',
                                                  values_to='frequency')


# Make boolean labels nice
ag85b_dat_long$boolean <- gsub("[[:punct:]]", "+", ag85b_dat_long$boolean)
ag85b_dat_long$boolean <- gsub("\\++", "+", ag85b_dat_long$boolean)
ag85b_dat_long$boolean <- gsub("c\\+", "-", ag85b_dat_long$boolean)
ag85b_dat_long$boolean <- gsub("TNF", "TNF+", ag85b_dat_long$boolean)
ag85b_dat_long$boolean <- gsub("\\+c", "-", ag85b_dat_long$boolean)

# Convert to factors
ag85b_dat_long$boolean <- as.factor(ag85b_dat_long$boolean)
ag85b_dat_long$QFT <- as.factor(ag85b_dat_long$QFT)

# Reorder levels of boolean factor 
ag85b_dat_long$boolean <- factor(ag85b_dat_long$boolean, 
                                 levels=c("R7+RA+", "R7+RA-G-L2+TNF-", "R7-RA-G-L2+TNF-",
                                          "R7+RA-G-L2-TNF+", "R7-RA-G-L2-TNF+", "R7+RA-G-L2+TNF+",
                                          "R7-RA-G-L2+TNF+", "R7+RA-G+L2+TNF+", "R7-RA-G+L2+TNF+", 
                                          "R7+RA-G+L2+TNF-", "R7-RA-G+L2+TNF-", "R7+RA-G+L2-TNF+", 
                                          "R7-RA-G+L2-TNF+",  "R7+RA-G+L2-TNF-", "R7-RA-G+L2-TNF-", 
                                          "R7-RA+"))
levels(ag85b_dat_long$boolean)

#--------------------------- PEAK:

# Create Grob function for highlighting significant p-values later:
find_cell <- function(table, row, col, name="core-fg"){
  l <- table$layout
  which(l$t==row & l$l==col & l$name==name)
}

# Conduct Mann Whitney Test
dftest    <- ag85b_dat_long[ag85b_dat_long[,"timepnt"]=="peak",]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))

# Log Transform
df <- dftest
df$frequency <- log(df$frequency+0.00001)

# Label positions:
pos <- aggregate(df$frequency, by = list(df$boolean), max)$x

# plot
p <- ggplot(df, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color=rep(c("#660000", "#000066"), 16)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab("")+
  ylim(-12, -0.8)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("TNF", "IL2", "G", "RA", "R7", "P-value"), 
                          x1 = c("", "", "", "+", "+", stat.test$padj[1]),
                          x2 = c("-", "+", "-", "-", "+", stat.test$padj[2]),
                          x3 = c("-", "+", "-", "-", "-", stat.test$padj[3]),
                          x4 = c("+", "-", "-", "-", "+", stat.test$padj[4]),
                          x5 = c("+", "-", "-", "-", "-", stat.test$padj[5]),
                          x6 = c("+", "+", "-", "-", "+", stat.test$padj[6]),
                          x7 = c("+", "+", "-", "-", "-", stat.test$padj[7]),
                          x8 = c("+", "+", "+", "-", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "-", "-", stat.test$padj[9]),
                          x10 = c("-", "+", "+", "-", "+", stat.test$padj[10]),
                          x11 = c("-", "+", "+", "-", "-", stat.test$padj[11]),
                          x12 = c("+", "-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "-", "+", "-", "-", stat.test$padj[13]),
                          x14 = c("-", "-", "+", "-", "+", stat.test$padj[14]),
                          x15 = c("-", "-", "+", "-", "-", stat.test$padj[15]),
                          x16 = c("", "", "", "+", "-", stat.test$padj[16])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set significant p-values to bold:
ind <- find_cell(tab, 6, 6, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", "CCR7", "CD45RA", expression(paste("IFN", gamma)), "IL2", "TNF"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))







#--------------------------- MEMORY:
# Conduct Mann Whitney Test
dftest <- ag85b_dat_long[ag85b_dat_long[,"timepnt"]=="memory",]
stat.test <- dftest %>%
  group_by(boolean) %>%
  wilcox_test(frequency ~ QFT, paired = F)
(stat.test <- stat.test %>%
    add_xy_position(x = "boolean", dodge = 0.8))
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))

# Log Transform
df <- dftest
df$frequency <- log(df$frequency+0.00001)

# Label positions:
pos <- aggregate(df$frequency, by = list(df$boolean), max)$x

# plot
p <- ggplot(df, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color=rep(c("#660000", "#000066"), 16)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25)) +
  ylab("Logged % of CD4 cells") +
  xlab("")+
  ylim(-12, -0.8)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("TNF", "IL2", "G", "RA", "R7", "P-value"), 
                          x1 = c("", "", "", "+", "+", stat.test$padj[1]),
                          x2 = c("-", "+", "-", "-", "+", stat.test$padj[2]),
                          x3 = c("-", "+", "-", "-", "-", stat.test$padj[3]),
                          x4 = c("+", "-", "-", "-", "+", stat.test$padj[4]),
                          x5 = c("+", "-", "-", "-", "-", stat.test$padj[5]),
                          x6 = c("+", "+", "-", "-", "+", stat.test$padj[6]),
                          x7 = c("+", "+", "-", "-", "-", stat.test$padj[7]),
                          x8 = c("+", "+", "+", "-", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "-", "-", stat.test$padj[9]),
                          x10 = c("-", "+", "+", "-", "+", stat.test$padj[10]),
                          x11 = c("-", "+", "+", "-", "-", stat.test$padj[11]),
                          x12 = c("+", "-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "-", "+", "-", "-", stat.test$padj[13]),
                          x14 = c("-", "-", "+", "-", "+", stat.test$padj[14]),
                          x15 = c("-", "-", "+", "-", "-", stat.test$padj[15]),
                          x16 = c("", "", "", "+", "-", stat.test$padj[16])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", "CCR7", "CD45RA", expression(paste("IFN", gamma)), "IL2", "TNF"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))






#---------------------------------------------- ESAT6

esat6_dat_long <- esat6_dat_smol %>% pivot_longer(cols='R7..RA..':'R7c..RA..',
                                                  names_to='boolean',
                                                  values_to='frequency')


# Make boolean labels nice
esat6_dat_long$boolean <- gsub("[[:punct:]]", "+", esat6_dat_long$boolean)
esat6_dat_long$boolean <- gsub("\\++", "+", esat6_dat_long$boolean)
esat6_dat_long$boolean <- gsub("c\\+", "-", esat6_dat_long$boolean)
esat6_dat_long$boolean <- gsub("TNF", "TNF+", esat6_dat_long$boolean)
esat6_dat_long$boolean <- gsub("\\+c", "-", esat6_dat_long$boolean)

# Convert to factors
esat6_dat_long$boolean <- as.factor(esat6_dat_long$boolean)
esat6_dat_long$QFT <- as.factor(esat6_dat_long$QFT)

# Reorder levels of boolean factor 
esat6_dat_long$boolean <- factor(esat6_dat_long$boolean, 
                                 levels=c("R7+RA+", "R7+RA-G-L2+TNF-", "R7-RA-G-L2+TNF-",
                                          "R7+RA-G-L2-TNF+", "R7-RA-G-L2-TNF+", "R7+RA-G-L2+TNF+",
                                          "R7-RA-G-L2+TNF+", "R7+RA-G+L2+TNF+", "R7-RA-G+L2+TNF+", 
                                          "R7+RA-G+L2+TNF-", "R7-RA-G+L2+TNF-", "R7+RA-G+L2-TNF+", 
                                          "R7-RA-G+L2-TNF+",  "R7+RA-G+L2-TNF-", "R7-RA-G+L2-TNF-", 
                                          "R7-RA+"))
levels(esat6_dat_long$boolean)

#--------------------------- PEAK:
# Conduct Mann Whitney Test
dftest <- esat6_dat_long[esat6_dat_long[,"timepnt"]=="peak",]
stat.test <- dftest %>%
  group_by(boolean) %>%
  wilcox_test(frequency ~ QFT,  paired = F) 
(stat.test <- stat.test %>%
    add_xy_position(x = "boolean", dodge = 0.8))
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$padj[11] <- round(p.adjust(stat.test$p, method = "BH"),4)[11]  # don't want zero here
stat.test$padj[15] <- round(p.adjust(stat.test$p, method = "BH"),4)[15] # don't want zero here

# Log Transform
df <- dftest
df$frequency <- log(df$frequency+0.00001)

# Label positions:
pos <- aggregate(df$frequency, by = list(df$boolean), max)$x

# plot
p <- ggplot(df, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color=rep(c("#660000", "#000066"), 16)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25)) +
  ylab("Logged % of CD4 cells") +
  xlab("")+
  ylim(-12, -0.8)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("TNF", "IL2", "G", "RA", "R7", "P-value"), 
                          x1 = c("", "", "", "+", "+", stat.test$padj[1]),
                          x2 = c("-", "+", "-", "-", "+", stat.test$padj[2]),
                          x3 = c("-", "+", "-", "-", "-", stat.test$padj[3]),
                          x4 = c("+", "-", "-", "-", "+", stat.test$padj[4]),
                          x5 = c("+", "-", "-", "-", "-", stat.test$padj[5]),
                          x6 = c("+", "+", "-", "-", "+", stat.test$padj[6]),
                          x7 = c("+", "+", "-", "-", "-", stat.test$padj[7]),
                          x8 = c("+", "+", "+", "-", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "-", "-", stat.test$padj[9]),
                          x10 = c("-", "+", "+", "-", "+", stat.test$padj[10]),
                          x11 = c("-", "+", "+", "-", "-", stat.test$padj[11]),
                          x12 = c("+", "-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "-", "+", "-", "-", stat.test$padj[13]),
                          x14 = c("-", "-", "+", "-", "+", stat.test$padj[14]),
                          x15 = c("-", "-", "+", "-", "-", stat.test$padj[15]),
                          x16 = c("", "", "", "+", "-", stat.test$padj[16])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 16), cols=NULL)
# Set significant p-values to bold:
ind1 <- find_cell(tab, 6, 3, "core-fg")
ind2 <- find_cell(tab, 6, 10, "core-fg")
ind3 <- find_cell(tab, 6, 11, "core-fg")
ind4 <- find_cell(tab, 6, 14, "core-fg")
ind5 <- find_cell(tab, 6, 15, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=16, fontface="bold")
tab$grobs[ind2][[1]][["gp"]] <- gpar(fontsize=16, fontface="bold")
tab$grobs[ind3][[1]][["gp"]] <- gpar(fontsize=16, fontface="bold")
tab$grobs[ind4][[1]][["gp"]] <- gpar(fontsize=16, fontface="bold")
tab$grobs[ind5][[1]][["gp"]] <- gpar(fontsize=16, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", "CCR7", "CD45RA", expression(paste("IFN", gamma)), "IL2", "TNF"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))





#--------------------------- MEMORY:
# Conduct Mann Whitney Test
dftest <- esat6_dat_long[esat6_dat_long[,"timepnt"]=="memory",]
stat.test <- dftest %>%
  group_by(boolean) %>%
  wilcox_test(frequency ~ QFT,  paired = F) 
(stat.test <- stat.test %>%
    add_xy_position(x = "boolean", dodge = 0.8))
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$padj[8] <- round(p.adjust(stat.test$p, method = "BH"),4)[8]  # don't want zero here
stat.test$padj[9] <- round(p.adjust(stat.test$p, method = "BH"),4)[9]  # don't want zero here

# Log Transform
df <- dftest
df$frequency <- log(df$frequency+0.00001)

# Label positions:
pos <- aggregate(df$frequency, by = list(df$boolean), max)$x

# plot
p <- ggplot(df, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color=rep(c("#660000", "#000066"), 16)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25)) +
  ylab("Logged % of CD4 cells") +
  xlab("")+
  ylim(-12, -0.8)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("TNF", "IL2", "G", "RA", "R7", "P-value"), 
                          x1 = c("", "", "", "+", "+", stat.test$padj[1]),
                          x2 = c("-", "+", "-", "-", "+", stat.test$padj[2]),
                          x3 = c("-", "+", "-", "-", "-", stat.test$padj[3]),
                          x4 = c("+", "-", "-", "-", "+", stat.test$padj[4]),
                          x5 = c("+", "-", "-", "-", "-", stat.test$padj[5]),
                          x6 = c("+", "+", "-", "-", "+", stat.test$padj[6]),
                          x7 = c("+", "+", "-", "-", "-", stat.test$padj[7]),
                          x8 = c("+", "+", "+", "-", "+", stat.test$padj[8]),
                          x9 = c("+", "+", "+", "-", "-", stat.test$padj[9]),
                          x10 = c("-", "+", "+", "-", "+", stat.test$padj[10]),
                          x11 = c("-", "+", "+", "-", "-", stat.test$padj[11]),
                          x12 = c("+", "-", "+", "-", "+", stat.test$padj[12]),
                          x13 = c("+", "-", "+", "-", "-", stat.test$padj[13]),
                          x14 = c("-", "-", "+", "-", "+", stat.test$padj[14]),
                          x15 = c("-", "-", "+", "-", "-", stat.test$padj[15]),
                          x16 = c("", "", "", "+", "-", stat.test$padj[16])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 15), cols=NULL)
# Set significant p-values to bold:
ind1 <- find_cell(tab, 6, 1, "core-fg")
ind2 <- find_cell(tab, 6, 2, "core-fg")
ind3 <- find_cell(tab, 6, 4, "core-fg")
ind4 <- find_cell(tab, 6, 6, "core-fg")
ind5 <- find_cell(tab, 6, 8, "core-fg")
ind6 <- find_cell(tab, 6, 9, "core-fg")
ind7 <- find_cell(tab, 6, 10, "core-fg")
ind8 <- find_cell(tab, 6, 11, "core-fg")
ind9 <- find_cell(tab, 6, 12, "core-fg")
ind10 <- find_cell(tab, 6, 13, "core-fg")
ind11 <- find_cell(tab, 6, 14, "core-fg")
ind12 <- find_cell(tab, 6, 15, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
tab$grobs[ind2][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
tab$grobs[ind3][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
tab$grobs[ind4][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
tab$grobs[ind5][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
tab$grobs[ind6][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
tab$grobs[ind7][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
tab$grobs[ind8][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
tab$grobs[ind9][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
tab$grobs[ind10][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
tab$grobs[ind11][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
tab$grobs[ind12][[1]][["gp"]] <- gpar(fontsize=15, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", "CCR7", "CD45RA", expression(paste("IFN", gamma)), "IL2", "TNF"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))




#----------------------------------------------
# Summing Over the Memory Markers
#----------------------------------------------

ag85b_dat_nm <- data.frame( 
  "Gc..L2..TNFc" = ag85b_dat$R7..RAc..Gc..L2..TNFc + ag85b_dat$R7c..RAc..Gc..L2..TNFc +  ag85b_dat$R7..RA..Gc..L2..TNFc + ag85b_dat$R7c..RA..Gc..L2..TNFc,
  "Gc..L2c..TNF" = ag85b_dat$R7..RAc..Gc..L2c..TNF + ag85b_dat$R7c..RAc..Gc..L2c..TNF + ag85b_dat$R7..RA..Gc..L2c..TNF + ag85b_dat$R7c..RA..Gc..L2c..TNF, 
  "Gc..L2..TNF" = ag85b_dat$R7..RAc..Gc..L2..TNF + ag85b_dat$R7c..RAc..Gc..L2..TNF + ag85b_dat$R7..RA..Gc..L2..TNF + ag85b_dat$R7c..RA..Gc..L2..TNF, 
  "G..L2..TNF" = ag85b_dat$R7..RAc..G..L2..TNF + ag85b_dat$R7c..RAc..G..L2..TNF + ag85b_dat$R7..RA..G..L2..TNF + ag85b_dat$R7c..RA..G..L2..TNF, 
  "G..L2..TNFc" = ag85b_dat$R7..RAc..G..L2..TNFc + ag85b_dat$R7c..RAc..G..L2..TNFc + ag85b_dat$R7..RA..G..L2..TNFc + ag85b_dat$R7c..RA..G..L2..TNFc,
  "G..L2c..TNF" = ag85b_dat$R7..RAc..G..L2c..TNF + ag85b_dat$R7c..RAc..G..L2c..TNF + ag85b_dat$R7..RA..G..L2c..TNF + ag85b_dat$R7c..RA..G..L2c..TNF,
  "G..L2c..TNFc" = ag85b_dat$R7..RAc..G..L2c..TNFc + ag85b_dat$R7c..RAc..G..L2c..TNFc + ag85b_dat$R7..RA..G..L2c..TNFc + ag85b_dat$R7c..RA..G..L2c..TNFc, 
  "timepnt" = ag85b_dat$timepnt, 
  "QFT" = ag85b_dat$QFT)

ag85b_dat_nm_long <- ag85b_dat_nm %>% pivot_longer(cols='Gc..L2..TNFc':'G..L2c..TNFc',
                                                   names_to='boolean',
                                                   values_to='frequency')


# Make boolean labels nice
ag85b_dat_nm_long$boolean <- gsub("[[:punct:]]", "+", ag85b_dat_nm_long$boolean)
ag85b_dat_nm_long$boolean <- gsub("\\++", "+", ag85b_dat_nm_long$boolean)
ag85b_dat_nm_long$boolean <- gsub("c\\+", "-", ag85b_dat_nm_long$boolean)
ag85b_dat_nm_long$boolean <- gsub("TNF", "TNF+", ag85b_dat_nm_long$boolean)
ag85b_dat_nm_long$boolean <- gsub("\\+c", "-", ag85b_dat_nm_long$boolean)

# Convert to factors
ag85b_dat_nm_long$boolean <- as.factor(ag85b_dat_nm_long$boolean)
ag85b_dat_nm_long$QFT <- as.factor(ag85b_dat_nm_long$QFT)

# Reorder levels of boolean factor 
ag85b_dat_nm_long$boolean <- factor(ag85b_dat_nm_long$boolean, 
                                    levels=c( 
                                      "G-L2+TNF-", 
                                      "G-L2-TNF+",  
                                      "G-L2+TNF+",
                                      "G+L2+TNF+", 
                                      "G+L2+TNF-", 
                                      "G+L2-TNF+", 
                                      "G+L2-TNF-"))
levels(ag85b_dat_nm_long$boolean)



#--------------------------- PEAK:
# Conduct Mann Whitney Test
dftest <- ag85b_dat_nm_long[ag85b_dat_nm_long[,"timepnt"]=="peak",]
stat.test <- dftest %>%
  group_by(boolean) %>%
  wilcox_test(frequency ~ QFT,  paired = F) 
(stat.test <- stat.test %>%
    add_xy_position(x = "boolean", dodge = 0.8))
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))

# Log Transform
df <- dftest
df$frequency <- log(df$frequency+0.00001)

# Label positions:
pos <- aggregate(df$frequency, by = list(df$boolean), max)$x

# plot
p <- ggplot(df, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color=rep(c("#660000", "#000066"), 7)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
         axis.title = element_text(size=25)) +
  ylab("Logged % of CD4 cells") +
  xlab("")+
  ylim(-12, -0.8)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")

# Make grobtable:
df1         <- data.frame(Marker = c("TNF", "IL2", "G", "P-value"), 
                          x1 = c("-", "+", "-", stat.test$padj[1]),
                          x2 = c("+", "-", "-", stat.test$padj[2]),
                          x3 = c("+", "+", "-", stat.test$padj[3]),
                          x4 = c("+", "+", "+", stat.test$padj[4]),
                          x5 = c("-", "+", "+", stat.test$padj[5]),
                          x6 = c("+", "-", "+", stat.test$padj[6]),
                          x7 = c("-", "-", "+", stat.test$padj[7])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 30), cols=NULL)
# Set significant p-values to bold:
which(stat.test$padj<=0.05)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF"))+
  theme(axis.text.y = element_text(face="bold", size=22))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))




#--------------------------- MEMORY:
# Conduct Mann Whitney Test
dftest <- ag85b_dat_nm_long[ag85b_dat_nm_long[,"timepnt"]=="memory",]
stat.test <- dftest %>%
  group_by(boolean) %>%
  wilcox_test(frequency ~ QFT,  paired = F) 
(stat.test <- stat.test %>%
    add_xy_position(x = "boolean", dodge = 0.8))
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))

# Log Transform
df <- dftest
df$frequency <- log(df$frequency+0.00001)

# Label positions:
pos <- aggregate(df$frequency, by = list(df$boolean), max)$x

# plot
p <- ggplot(df, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color=rep(c("#660000", "#000066"), 7)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25)) +
  ylab("Logged % of CD4 cells") +
  xlab("")+
  ylim(-12, -0.8)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("TNF", "IL2", "G", "P-value"), 
                          x1 = c("-", "+", "-", stat.test$padj[1]),
                          x2 = c("+", "-", "-", stat.test$padj[2]),
                          x3 = c("+", "+", "-", stat.test$padj[3]),
                          x4 = c("+", "+", "+", stat.test$padj[4]),
                          x5 = c("-", "+", "+", stat.test$padj[5]),
                          x6 = c("+", "-", "+", stat.test$padj[6]),
                          x7 = c("-", "-", "+", stat.test$padj[7])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 30), cols=NULL)
# Set significant p-values to bold:
which(stat.test$padj<=0.05)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF"))+
  theme(axis.text.y = element_text(face="bold", size=22))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))





#--------------------------------------------------------------------------- ESAT6

esat6_dat_nm <- data.frame(
  "Gc..L2..TNFc" = esat6_dat$R7..RAc..Gc..L2..TNFc + esat6_dat$R7c..RAc..Gc..L2..TNFc + esat6_dat$R7..RA..Gc..L2..TNFc + esat6_dat$R7c..RA..Gc..L2..TNFc,
  "Gc..L2c..TNF" = esat6_dat$R7..RAc..Gc..L2c..TNF + esat6_dat$R7c..RAc..Gc..L2c..TNF + esat6_dat$R7..RA..Gc..L2c..TNF + esat6_dat$R7c..RA..Gc..L2c..TNF, 
  "Gc..L2..TNF" = esat6_dat$R7..RAc..Gc..L2..TNF + esat6_dat$R7c..RAc..Gc..L2..TNF + esat6_dat$R7..RA..Gc..L2..TNF + esat6_dat$R7c..RA..Gc..L2..TNF, 
  "G..L2..TNF" = esat6_dat$R7..RAc..G..L2..TNF + esat6_dat$R7c..RAc..G..L2..TNF + esat6_dat$R7..RA..G..L2..TNF + esat6_dat$R7c..RA..G..L2..TNF, 
  "G..L2..TNFc" = esat6_dat$R7..RAc..G..L2..TNFc + esat6_dat$R7c..RAc..G..L2..TNFc + esat6_dat$R7..RA..G..L2..TNFc + esat6_dat$R7c..RA..G..L2..TNFc,
  "G..L2c..TNF" = esat6_dat$R7..RAc..G..L2c..TNF + esat6_dat$R7c..RAc..G..L2c..TNF + esat6_dat$R7..RA..G..L2c..TNF + esat6_dat$R7c..RA..G..L2c..TNF,
  "G..L2c..TNFc" = esat6_dat$R7..RAc..G..L2c..TNFc + esat6_dat$R7c..RAc..G..L2c..TNFc + esat6_dat$R7..RA..G..L2c..TNFc + esat6_dat$R7c..RA..G..L2c..TNFc, 
  
  "timepnt" = esat6_dat$timepnt, 
  "QFT" = esat6_dat$QFT)

esat6_dat_nm_long <- esat6_dat_nm %>% pivot_longer(cols='Gc..L2..TNFc':'G..L2c..TNFc',
                                                   names_to='boolean',
                                                   values_to='frequency')


# Make boolean labels nice
esat6_dat_nm_long$boolean <- gsub("[[:punct:]]", "+", esat6_dat_nm_long$boolean)
esat6_dat_nm_long$boolean <- gsub("\\++", "+", esat6_dat_nm_long$boolean)
esat6_dat_nm_long$boolean <- gsub("c\\+", "-", esat6_dat_nm_long$boolean)
esat6_dat_nm_long$boolean <- gsub("TNF", "TNF+", esat6_dat_nm_long$boolean)
esat6_dat_nm_long$boolean <- gsub("\\+c", "-", esat6_dat_nm_long$boolean)

# Convert to factors
esat6_dat_nm_long$boolean <- as.factor(esat6_dat_nm_long$boolean)
esat6_dat_nm_long$QFT <- as.factor(esat6_dat_nm_long$QFT)

# Reorder levels of boolean factor 
esat6_dat_nm_long$boolean <- factor(esat6_dat_nm_long$boolean, 
                                    levels=c(
                                      "G-L2+TNF-", 
                                      "G-L2-TNF+",  
                                      "G-L2+TNF+",
                                      "G+L2+TNF+", 
                                      "G+L2+TNF-", 
                                      "G+L2-TNF+", 
                                      "G+L2-TNF-"))
levels(esat6_dat_nm_long$boolean)



#--------------------------- PEAK:
# Conduct Mann Whitney Test
dftest <- esat6_dat_nm_long[esat6_dat_nm_long[,"timepnt"]=="peak",]
stat.test <- dftest %>%
  group_by(boolean) %>%
  wilcox_test(frequency ~ QFT,  paired = F) 
(stat.test <- stat.test %>%
    add_xy_position(x = "boolean", dodge = 0.8))
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
(stat.test$padj[5] <- round(p.adjust(stat.test$p, method = "BH"),4)[5])  # don't want zero here
(stat.test$padj[7] <- round(p.adjust(stat.test$p, method = "BH"),4)[7])  # don't want zero here

# Log Transform
df <- dftest
df$frequency <- log(df$frequency+0.00001)

# Label positions:
pos <- aggregate(df$frequency, by = list(df$boolean), max)$x

# plot
p <- ggplot(df, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color=rep(c("#660000", "#000066"), 7)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25)) +
  ylab("Logged % of CD4 cells") +
  xlab("")+
  ylim(-12,-0.8)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("TNF", "IL2", "G", "P-value"), 
                          x1 = c("-", "+", "-", stat.test$padj[1]),
                          x2 = c("+", "-", "-", stat.test$padj[2]),
                          x3 = c("+", "+", "-", stat.test$padj[3]),
                          x4 = c("+", "+", "+", stat.test$padj[4]),
                          x5 = c("-", "+", "+", stat.test$padj[5]),
                          x6 = c("+", "-", "+", stat.test$padj[6]),
                          x7 = c("-", "-", "+", stat.test$padj[7])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 30), cols=NULL)
# Set significant p-values to bold:
which(stat.test$padj<=0.05)
ind1 <- find_cell(tab, 4, 5, "core-fg")
ind2 <- find_cell(tab, 4, 7, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=30, fontface="bold")
tab$grobs[ind2][[1]][["gp"]] <- gpar(fontsize=30, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF"))+
  theme(axis.text.y = element_text(face="bold", size=22))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))





#--------------------------- MEMORY:
# Conduct Mann Whitney Test
dftest <- esat6_dat_nm_long[esat6_dat_nm_long[,"timepnt"]=="memory",]
stat.test <- dftest %>%
  group_by(boolean) %>%
  wilcox_test(frequency ~ QFT, paired = F) 
(stat.test <- stat.test %>%
    add_xy_position(x = "boolean", dodge = 0.8))
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
(stat.test$padj[4] <- round(p.adjust(stat.test$p, method = "BH"),4)[4])  # don't want zero here


# Log Transform
df <- dftest
df$frequency <- log(df$frequency+0.00001)

# Label positions:
pos <- aggregate(df$frequency, by = list(df$boolean), max)$x

# plot
p <- ggplot(df, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color=rep(c("#660000", "#000066"), 7)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
         axis.title = element_text(size=25)) +
  ylab("Logged % of CD4 cells") +
  xlab("")+
  ylim(-12,-0.8)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")
# Make grobtable:
df1         <- data.frame(Marker = c("TNF", "IL2", "G", "P-value"), 
                          x1 = c("-", "+", "-", stat.test$padj[1]),
                          x2 = c("+", "-", "-", stat.test$padj[2]),
                          x3 = c("+", "+", "-", stat.test$padj[3]),
                          x4 = c("+", "+", "+", stat.test$padj[4]),
                          x5 = c("-", "+", "+", stat.test$padj[5]),
                          x6 = c("+", "-", "+", stat.test$padj[6]),
                          x7 = c("-", "-", "+", stat.test$padj[7])) 
tab         <- tableGrob(df1[,-1], rows=NULL, theme = ttheme_default(base_size = 30), cols=NULL)
# Set significant p-values to bold:
which(stat.test$padj<=0.05)
ind1 <- find_cell(tab, 4, 2, "core-fg")
ind2 <- find_cell(tab, 4, 3, "core-fg")
ind3 <- find_cell(tab, 4, 4, "core-fg")
ind4 <- find_cell(tab, 4, 5, "core-fg")
ind5 <- find_cell(tab, 4, 6, "core-fg")
ind6 <- find_cell(tab, 4, 7, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=30, fontface="bold")
tab$grobs[ind2][[1]][["gp"]] <- gpar(fontsize=30, fontface="bold")
tab$grobs[ind3][[1]][["gp"]] <- gpar(fontsize=30, fontface="bold")
tab$grobs[ind4][[1]][["gp"]] <- gpar(fontsize=30, fontface="bold")
tab$grobs[ind5][[1]][["gp"]] <- gpar(fontsize=30, fontface="bold")
tab$grobs[ind6][[1]][["gp"]] <- gpar(fontsize=30, fontface="bold")
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(df1$Marker),limits = rev(df1$Marker), 
                   labels = c("P-value", expression(paste("IFN", gamma)), "IL2", "TNF"))+
  theme(axis.text.y = element_text(face="bold", size=22))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))














