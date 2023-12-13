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


# Only keep Ag85B and ESAT6
cd4_bs <- cd4_bs[cd4_bs[,"Stim"] == "ESAT6" | cd4_bs[,"Stim"] == "Ag85B",]

#----------------------------
# Get Rid of Unwanted Combos
#----------------------------

# Remove booleans that express no cytokines (note that we are ignoring il17)
bsfreq <- cd4_bs[,-c(which(names(cd4_bs) == "R7c..RAc..Gc..L2c..L17c..TNFc"),
                     which(names(cd4_bs) == "R7..RAc..Gc..L2c..L17c..TNFc"),
                     which(names(cd4_bs) == "R7c..RA..Gc..L2c..L17c..TNFc"),
                     which(names(cd4_bs) == "R7..RA..Gc..L2c..L17c..TNFc"),
                     which(names(cd4_bs) == "R7c..RAc..Gc..L2c..L17..TNFc"), # here the only cytokine expressed is il17 - but we are ignoring 1l17 so this counts as no cytokines being expressed
                     which(names(cd4_bs) == "R7..RAc..Gc..L2c..L17..TNFc"),
                     which(names(cd4_bs) == "R7c..RA..Gc..L2c..L17..TNFc"),
                     which(names(cd4_bs) == "R7..RA..Gc..L2c..L17..TNFc"))]


# Sum over R7-RA+ and R7+RA+
bsfreq$R7..RA.. <- rowSums(bsfreq[,11:24])
bsfreq$R7c..RA.. <- rowSums(bsfreq[,39:52])


# Need to ignore IL17
# lets do sum checks to see how we can sum over IL17
x <-cbind(names(bsfreq)[which(grepl("..L17c..", names(bsfreq), fixed = TRUE)==T)],
          names(bsfreq)[which(grepl("..L17..", names(bsfreq), fixed = TRUE)==T)])
y <- gsub("L17c", "L17", x)
which(y[,1] != y[,2]) # 0, therefore the only difference between the two columns in any given row is 
# whether IL17 is expressed or not

# Now sum over IL17
dat <- matrix(NA, nrow=nrow(bsfreq), ncol = nrow(x))
dat <- as.data.frame(dat)
for (i in 1:nrow(x)) {
  dat[,i] <- bsfreq[,x[i,1]] + bsfreq[,x[i,2]]
  names(dat)[i] <- gsub("L17c..", "", x[i,1])
}
dat$Schedule = bsfreq$Schedule
dat$Dose = bsfreq$Dose
dat$timepnt <- bsfreq$timepnt
dat$QFT <- bsfreq$QFT
dat$PID = bsfreq$PID
dat$Stim = bsfreq$Stim
dat$R7..RA.. <- bsfreq$R7..RA..
dat$R7c..RA.. <- bsfreq$R7c..RA..
str(dat)




# Keep only columns of interest and reorder
dat_smol <- dat[,c("R7..RA..", "R7..RAc..Gc..L2..TNFc", "R7c..RAc..Gc..L2..TNFc", 
                   "R7..RAc..Gc..L2c..TNF", "R7c..RAc..Gc..L2c..TNF", "R7..RAc..Gc..L2..TNF",  
                   "R7c..RAc..Gc..L2..TNF", "R7..RAc..G..L2..TNF", "R7c..RAc..G..L2..TNF", 
                   "R7..RAc..G..L2..TNFc", "R7c..RAc..G..L2..TNFc", "R7..RAc..G..L2c..TNF",   
                   "R7c..RAc..G..L2c..TNF",  "R7..RAc..G..L2c..TNFc", "R7c..RAc..G..L2c..TNFc", 
                   "R7c..RA..","Schedule", "Dose", "timepnt", "QFT", "PID", "Stim")]





#---------------------
# Reshape Data
#---------------------

dat_long <- dat_smol %>% pivot_longer(cols='R7..RA..':'R7c..RA..',
                                      names_to='boolean',
                                      values_to='frequency')


# Make boolean labels nice
dat_long$boolean <- gsub("[[:punct:]]", "+", dat_long$boolean)
dat_long$boolean <- gsub("\\++", "+", dat_long$boolean)
dat_long$boolean <- gsub("c\\+", "-", dat_long$boolean)
dat_long$boolean <- gsub("TNF", "TNF+", dat_long$boolean)
dat_long$boolean <- gsub("\\+c", "-", dat_long$boolean)

# Convert to factors
dat_long$boolean <- as.factor(dat_long$boolean)


# Reorder levels of boolean factor 
dat_long$boolean <- factor(dat_long$boolean, 
                           levels=c("R7+RA+", "R7+RA-G-L2+TNF-", "R7-RA-G-L2+TNF-",
                                    "R7+RA-G-L2-TNF+", "R7-RA-G-L2-TNF+", "R7+RA-G-L2+TNF+",
                                    "R7-RA-G-L2+TNF+", "R7+RA-G+L2+TNF+", "R7-RA-G+L2+TNF+", 
                                    "R7+RA-G+L2+TNF-", "R7-RA-G+L2+TNF-", "R7+RA-G+L2-TNF+", 
                                    "R7-RA-G+L2-TNF+",  "R7+RA-G+L2-TNF-", "R7-RA-G+L2-TNF-", 
                                    "R7-RA+"))
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

# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="baseline" & df_ca[,"QFT"]=="QFT neg",]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ Stim, paired = T))
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
  geom_boxplot(alpha=1, size =1, aes(fill = Stim)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12, -1)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values = c("#669999", "#996699"), name="")
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


######################## 
# Peak QFT-:
######################## 
# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="peak" & df_ca[,"QFT"]=="QFT neg",]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ Stim, paired = T))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$padj[6] = 0.0005
stat.test$padj[7] = 0.0005

# Log Transform
df <- dftest
df$frequency <- log(df$frequency+0.00001)

# Label positions:
pos <- aggregate(df$frequency, by = list(df$boolean), max)$x

# plot
p <- ggplot(df, aes(y=frequency, x=boolean)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Stim)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12, -1)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values = c("#669999", "#996699"), name="")
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
ind <- find_cell(tab, 6, 2, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 3, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 4, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 5, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 6, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=14, fontface="bold")
ind <- find_cell(tab, 6, 7, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=14, fontface="bold")
ind <- find_cell(tab, 6, 8, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 9, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 10, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 11, "core-fg")
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


######################## 
# Memory QFT-:
######################## 

# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="memory" & df_ca[,"QFT"]=="QFT neg",]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ Stim, paired = T))
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
  geom_boxplot(alpha=1, size =1, aes(fill = Stim)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12, -1)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values = c("#669999", "#996699"), name="")
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
# Set significant p-values to bold:
ind1 <- find_cell(tab, 6, 2, "core-fg")
ind3 <- find_cell(tab, 6, 6, "core-fg")
ind4 <- find_cell(tab, 6, 7, "core-fg")
ind5 <- find_cell(tab, 6, 8, "core-fg")
ind6 <- find_cell(tab, 6, 9, "core-fg")
ind7 <- find_cell(tab, 6, 11, "core-fg")
ind8 <- find_cell(tab, 6, 3, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
tab$grobs[ind3][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
tab$grobs[ind4][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
tab$grobs[ind5][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
tab$grobs[ind6][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
tab$grobs[ind7][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
tab$grobs[ind8][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
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




######################## 
# Baseline QFT+:
######################## 

# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="baseline" & df_ca[,"QFT"]=="QFT pos",]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ Stim, paired = T))
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
  geom_boxplot(alpha=1, size =1, aes(fill = Stim)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12, -1)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values = c("#669999", "#996699"), name="")
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
ind1 <- find_cell(tab, 6, 13, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
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


######################## 
# Peak QFT+:
######################## 

# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="peak" & df_ca[,"QFT"]=="QFT pos",]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ Stim, paired = T))
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
  geom_boxplot(alpha=1, size =1, aes(fill = Stim)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12, -1)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values = c("#669999", "#996699"), name="")
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
ind1 <- find_cell(tab, 6, 13, "core-fg")
tab$grobs[ind1][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind2 <- find_cell(tab, 6, 15, "core-fg")
tab$grobs[ind2][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
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


######################## 
# Memory QFT+:
######################## 

# Conduct Mann Whitney Test
dftest    <- df_ca[df_ca[,"timepnt"]=="memory" & df_ca[,"QFT"]=="QFT pos",]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ Stim, paired = T))
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
  geom_boxplot(alpha=1, size =1, aes(fill = Stim)) +
  theme_minimal()+
  theme(text = element_text(size = 30),axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title = element_text(size=25))+
  ylab("Logged % of CD4 cells")+
  xlab(NULL)+
  ylim(-12, -1)+
  scale_x_discrete(labels = rep("", 16))+
  scale_fill_manual(values = c("#669999", "#996699"), name="")
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
ind <- find_cell(tab, 6, 1, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 8, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 11, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 13, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 15, "core-fg")
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















#########################################
# COMPARING QFT:
#########################################

# only 2 admin vaccine:
df_cq = dat_long[dat_long[,"Schedule"]==2,]

#####################
# Memory Ag85B 5mg
#####################

# Conduct Mann Whitney Test
dftest    <- df_cq[df_cq[,"timepnt"]=="memory" & df_cq[,"Stim"]=="Ag85B" & df_cq[,"Dose"]==5,]
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
  xlab(NULL)+
  ylim(-12, -0.6)+
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
ind <- find_cell(tab, 6, 2, "core-fg")
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



#####################
# Memory Ag85B 15mg
#####################

# Conduct Mann Whitney Test
dftest    <- df_cq[df_cq[,"timepnt"]=="memory" & df_cq[,"Stim"]=="Ag85B" & df_cq[,"Dose"]==15,]
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
  xlab(NULL)+
  ylim(-12, -0.5)+
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


#####################
# Memory Ag85B 50mg
#####################

# Conduct Mann Whitney Test
dftest    <- df_cq[df_cq[,"timepnt"]=="memory" & df_cq[,"Stim"]=="Ag85B" & df_cq[,"Dose"]==50,]
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
  xlab(NULL)+
  ylim(-12, -0.6)+
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



#####################
# Memory ESAT6 5mg
#####################

# Conduct Mann Whitney Test
dftest    <- df_cq[df_cq[,"timepnt"]=="memory" & df_cq[,"Stim"]=="ESAT6" & df_cq[,"Dose"]==5,]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$padj[8] = 0.0004
stat.test$padj[9] = 0.0004

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
  xlab(NULL)+
  ylim(-12, -0.6)+
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
ind <- find_cell(tab, 6, 1, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 2, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 4, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 6, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 8, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=14, fontface="bold")
ind <- find_cell(tab, 6, 9, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=14, fontface="bold")
ind <- find_cell(tab, 6, 11, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 12, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 13, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 15, "core-fg")
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



#####################
# Memory ESAT6 15mg
#####################

# Conduct Mann Whitney Test
dftest    <- df_cq[df_cq[,"timepnt"]=="memory" & df_cq[,"Stim"]=="ESAT6" & df_cq[,"Dose"]==15,]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$padj[9] = 0.0004

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
  xlab(NULL)+
  ylim(-12, -0.6)+
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
ind <- find_cell(tab, 6, 5, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 8, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 9, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=14, fontface="bold")
ind <- find_cell(tab, 6, 10, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 11, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 13, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 14, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 15, "core-fg")
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


#####################
# Memory ESAT6 50mg
#####################

# Conduct Mann Whitney Test
dftest    <- df_cq[df_cq[,"timepnt"]=="memory" & df_cq[,"Stim"]=="ESAT6" & df_cq[,"Dose"]==50,]
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
  xlab(NULL)+
  ylim(-12, -0.6)+
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
ind <- find_cell(tab, 6, 8, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 9, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 13, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 14, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 15, "core-fg")
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




############################################
# COMPARING QFT: 2 ADMIN OF 5mg ONLY
############################################

# only 2 admin vaccine:
df_cq25 = dat_long[dat_long[,"Schedule"]==2 & dat_long[,"Dose"]==5,]

##################
# Ag85B Peak:
##################
# Conduct Mann Whitney Test
dftest    <- df_cq25[df_cq25[,"timepnt"]=="peak" & df_cq25[,"Stim"]=="Ag85B",]
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
  xlab(NULL)+
  ylim(-12, -0.6)+
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


##################
# ESAT6 Peak:
##################
# Conduct Mann Whitney Test
dftest    <- df_cq25[df_cq25[,"timepnt"]=="peak" & df_cq25[,"Stim"]=="ESAT6",]
(stat.test <- dftest %>%
    group_by(boolean) %>%
    wilcox_test(frequency ~ QFT, paired = F))
stat.test <- stat.test %>%
  add_xy_position(x = "boolean", dodge = 0.8)
(stat.test$padj <- round(p.adjust(stat.test$p, method = "BH"),3))
stat.test$padj[11] = 0.0003
stat.test$padj[15] = 0.0003

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
  xlab(NULL)+
  ylim(-12, -0.6)+
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
ind <- find_cell(tab, 6, 3, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 10, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 11, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=14, fontface="bold")
ind <- find_cell(tab, 6, 14, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=17, fontface="bold")
ind <- find_cell(tab, 6, 15, "core-fg")
tab$grobs[ind][[1]][["gp"]] <- gpar(fontsize=14, fontface="bold")
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

