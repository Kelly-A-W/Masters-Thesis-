# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

#---------------------
# Load Libraries
#---------------------

library(devtools) 
library(tidyr)
library(dplyr)
library(Rmisc) # to get CI
library(fitdistrplus)
library(ggcorrplot)
library(data.table) # for melt function
library(tweedie)
library(cplm)
library(gamlss) # defines cdf and pdf of ZIP and ZINB


#---------------------
# Import Data
#---------------------

cd4_counts <- read.csv("clean_cd4_counts_dataset_1.csv", check.names = F)
str(cd4_counts)



#-------------------------------------------
# Calculate Memory Response Frequency (MRF)
#-------------------------------------------

# Only keep counts for cells expressing at least one of INFg, IL2, TNF or L17 (i.e. cytokine+)
cd4_counts <- cd4_counts[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd4_counts), fixed = TRUE)==T)]
str(cd4_counts)


# Sum over cytokines:
cd4_counts$R7..RA..   <- rowSums(cd4_counts[,8:22])
cd4_counts$R7..RAc..  <- rowSums(cd4_counts[,23:37])
cd4_counts$R7c..RA..  <- rowSums(cd4_counts[,38:52])
cd4_counts$R7c..RAc.. <- rowSums(cd4_counts[,53:67])


# Convert to Frequencies (%):
cd4_counts[,8:71] <- (cd4_counts[,8:71]/cd4_counts$CD4count)*100

# Check:
summary(rowSums(cd4_counts[,8:67]))
summary(rowSums(cd4_counts[,68:71]))
# Both are the same

# Remove unwanted columns:
cd4_counts <- cd4_counts[,-c(7:67)]

# Background Subtract
cd4_mrf <- data.frame()
for(i in 1:length(unique(cd4_counts$PID))){
  pid  <- unique(cd4_counts$PID)[i]
  dat  <- cd4_counts[which(cd4_counts$PID == pid),]
  days <- unique(dat$Day)
  for (j in 1:length(days)) {
    day <- days[j]
    df  <- dat[which(dat$Day == day),]
    df[-which(df$Stim == "UNS"),7:10]  <- df[-which(df$Stim == "UNS"),7:10] - df[rep(which(df$Stim == "UNS"), (length(df$Stim)-1)),7:10]
    cd4_mrf <- rbind(cd4_mrf,df)
  }
}
summary(cd4_mrf)
which(is.na(cd4_mrf)==T)


# Deal with Negative Frequencies:
for (i in 7:ncol(cd4_mrf)) {
  # select one column of marker combinations
  column <- cd4_mrf[,i]
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
    cd4_mrf[ind,i] <- 0
  }
  
}
# Check:
sapply(sign(cd4_mrf[,-(1:6)]), table)


# Remove UNS and PHA
cd4_mrf <- cd4_mrf[cd4_mrf[,"Stim"]=="Ag85B" | cd4_mrf[,"Stim"]=="ESAT6",]
cd4_mrf$Stim <- as.factor(cd4_mrf$Stim)
summary(cd4_mrf)


#---------------------
# Add Time Category
#---------------------
cd4_mrf <- data.frame(timepnt = NA, cd4_mrf)

cd4_mrf$timepnt[which(cd4_mrf$Day==0)]   <- "baseline"
cd4_mrf$timepnt[which(cd4_mrf$Day==14)]  <- "peak"
cd4_mrf$timepnt[which(cd4_mrf$Day==70)]  <- "peak"
cd4_mrf$timepnt[which(cd4_mrf$Day==126)] <- "peak"
cd4_mrf$timepnt[which(cd4_mrf$Day==224)] <- "memory"
cd4_mrf$timepnt[which(cd4_mrf$Day==292)] <- "memory"
cd4_mrf$timepnt[which(cd4_mrf$Day==210)] <- "memory"

which(is.na(cd4_mrf$timepnt) == T)

cd4_mrf$timepnt <- as.factor(cd4_mrf$timepnt)
cd4_mrf$timepnt <- factor(cd4_mrf$timepnt, levels=c("baseline", "peak", "memory"))
levels(cd4_mrf$timepnt)


#---------------------
# Add Dose
#---------------------
# Need to add in dose and schedule information
cd4_mrf <- data.frame(Dose = NA, cd4_mrf)

cd4_mrf[intersect(which(cd4_mrf$Study == "H56-035"),which(cd4_mrf$Group == 1)),1] <- 5
cd4_mrf[intersect(which(cd4_mrf$Study == "H56-035"),which(cd4_mrf$Group == 2)),1] <- 5
cd4_mrf[intersect(which(cd4_mrf$Study == "H56-035"),which(cd4_mrf$Group == 3)),1] <- 5
cd4_mrf[intersect(which(cd4_mrf$Study == "H56-035"),which(cd4_mrf$Group == 4)),1] <- 5

cd4_mrf[intersect(which(cd4_mrf$Study == "H1"),which(cd4_mrf$Group == 1)),1] <- 15
cd4_mrf[intersect(which(cd4_mrf$Study == "H1"),which(cd4_mrf$Group == 2)),1] <- 50
cd4_mrf[intersect(which(cd4_mrf$Study == "H1"),which(cd4_mrf$Group == 3)),1] <- 15
cd4_mrf[intersect(which(cd4_mrf$Study == "H1"),which(cd4_mrf$Group == 4)),1] <- 50

cd4_mrf[intersect(which(cd4_mrf$Study == "H56-032"),which(cd4_mrf$Group == 1)),1] <- 50
cd4_mrf[intersect(which(cd4_mrf$Study == "H56-032"),which(cd4_mrf$Group == 2)),1] <- 15
cd4_mrf[intersect(which(cd4_mrf$Study == "H56-032"),which(cd4_mrf$Group == 3)),1] <- 50

# check
cd4_mrf[which(is.na(cd4_mrf$Dose)==T),]



#---------------------
# Add Schedule
#---------------------
# Need to add in dose and schedule information
cd4_mrf <- data.frame(Schedule = NA, cd4_mrf)

cd4_mrf[intersect(which(cd4_mrf$Study == "H56-035"),which(cd4_mrf$Group == 1)),1] <- 2
cd4_mrf[intersect(which(cd4_mrf$Study == "H56-035"),which(cd4_mrf$Group == 2)),1] <- 3
cd4_mrf[intersect(which(cd4_mrf$Study == "H56-035"),which(cd4_mrf$Group == 3)),1] <- 2
cd4_mrf[intersect(which(cd4_mrf$Study == "H56-035"),which(cd4_mrf$Group == 4)),1] <- 3

cd4_mrf[intersect(which(cd4_mrf$Study == "H1"),which(cd4_mrf$Group == 1)),1] <- 2
cd4_mrf[intersect(which(cd4_mrf$Study == "H1"),which(cd4_mrf$Group == 2)),1] <- 2
cd4_mrf[intersect(which(cd4_mrf$Study == "H1"),which(cd4_mrf$Group == 3)),1] <- 1
cd4_mrf[intersect(which(cd4_mrf$Study == "H1"),which(cd4_mrf$Group == 4)),1] <- 1

cd4_mrf[intersect(which(cd4_mrf$Study == "H56-032"),which(cd4_mrf$Group == 1)),1] <- 3
cd4_mrf[intersect(which(cd4_mrf$Study == "H56-032"),which(cd4_mrf$Group == 2)),1] <- 3
cd4_mrf[intersect(which(cd4_mrf$Study == "H56-032"),which(cd4_mrf$Group == 3)),1] <- 3

# check 
cd4_mrf[which(is.na(cd4_mrf$Schedule)==T),]


#----------------------------
# Deal with Missing Values:
#----------------------------

# Number of missing observations:
which(table(cd4_mrf$PID)<6) # 5 pids with missing visits: 3020, 3026 and 411
cd4_mrf[cd4_mrf[,"PID"]=="3001 H56-035",] # Ag85B missing baseline
cd4_mrf[cd4_mrf[,"PID"]==3020,] # Ag85B and ESAT6 missing memory
cd4_mrf[cd4_mrf[,"PID"]==3024,]  # Ag85B missing baseline
cd4_mrf[cd4_mrf[,"PID"]==3026,] # ESAT6 and AG85B missing baseline
cd4_mrf[cd4_mrf[,"PID"]==411,]  # ESAT6 and Ag85B missing peak and memory - think  we should remove

# Remove 411:
cd4_mrf <- cd4_mrf[-which(cd4_mrf$PID == 411),]


#---------------------
# Prepare Data
#---------------------


# Two-Level Grouping:
long <- melt(setDT(cd4_mrf), id.vars = c("Schedule","Dose", "timepnt", "Study", "QFT", "Group", "PID", "Day", "Stim"),
             variable.name = "Marker")
long$Marker <- as.factor(long$Marker)
long = as.data.frame(long)
head(long)
dim(long)
esat6 <- long[long[,"Stim"] == "ESAT6",]
ag85b <- long[long[,"Stim"] == "Ag85B",]
head(esat6)
head(ag85b)



#---------------------------------------------------------------------
# Univariate EDA: Choosing an Appropriate Distribution for Y:
#---------------------------------------------------------------------



#---------------------
# Density Plots:
#---------------------


#----------------------------------------------------- Two-Level Grouping (all markers):
ggplot(esat6, aes(x = value)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("Memory Response Frequency")

ggplot(ag85b, aes(x = value)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("Memory Response Frequency")


#----------------------------------------------------- Two-Level Grouping (2 markers):
esat6_2 = esat6[esat6[,"Marker"]=="R7..RAc.." | esat6[,"Marker"]=="R7c..RAc..",]
ag85b_2 = ag85b[ag85b[,"Marker"]=="R7..RAc.." | ag85b[,"Marker"]=="R7c..RAc..",]

ggplot(esat6_2, aes(x = value)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("Memory Response Frequency")

ggplot(ag85b_2, aes(x = value)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("Memory Response Frequency")


#-------------------------
# Fitting a Distribution
#-------------------------

#-------------------------------------------------------------------------------
#----------------------------------------------------------- Two-Level Grouping: 2 markers only
#-------------------------------------------------------------------------------

#---------------------------------------------------------------------- ESAT6
f1.1 <- fitdist((esat6_2$value+0.001), "lnorm")
f1.3 <- fitdist((esat6_2$value+0.001), "gamma")
f1.4 <- fitdist((esat6_2$value+0.001), "weibull")
model_1 <- cpglm(
  value ~1,
  data = esat6_2,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)



# QQ PLOTS:
ggplot(esat6_2, aes(sample=value)) +
  stat_qq(distribution = qtweedie, dparams = list(phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients)))+
  stat_qq_line(distribution = qtweedie, dparams = list(phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients)), color="red", linewidth=1)+
  theme_minimal(base_size = 27)+
  xlab("Theoretical quantiles")+
  ylab("Empirical quantiles")

qqcomp(f1.1, addlegend=F, plotstyle = "ggplot") + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)

qqcomp(f1.3, addlegend=F, plotstyle = "ggplot") + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)

qqcomp(f1.4, addlegend=F, plotstyle = "ggplot") + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)




#---------------------------------------------------------------------- Ag85B
f1.1 <- fitdist((ag85b_2$value+0.001), "lnorm")
f1.3 <- fitdist((ag85b_2$value+0.001), "gamma")
f1.4 <- fitdist((ag85b_2$value+0.001), "weibull")
model_1 <- cpglm(
  value ~1,
  data = ag85b_2,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)


# QQ PLOTS:
ggplot(ag85b_2, aes(sample=value)) +
  stat_qq(distribution = qtweedie, dparams = list(phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients)))+
  stat_qq_line(distribution = qtweedie, dparams = list(phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients)), color="red", linewidth=1)+
  theme_minimal(base_size = 27)+
  xlab("Theoretical quantiles")+
  ylab("Empirical quantiles")

qqcomp(f1.1, addlegend=F, plotstyle = "ggplot") + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)

qqcomp(f1.3, addlegend=F, plotstyle = "ggplot") + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)

qqcomp(f1.4, addlegend=F, plotstyle = "ggplot") + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)


