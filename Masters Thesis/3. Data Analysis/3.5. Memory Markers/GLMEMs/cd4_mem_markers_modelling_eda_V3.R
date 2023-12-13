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
library(data.table) # fot melt function
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

# Single-Level Grouping:
# subset data by stim (will fit a response to each marker column of each dataset)
esat6_mrf <- cd4_mrf[cd4_mrf[,"Stim"] == "ESAT6",]
ag85b_mrf <- cd4_mrf[cd4_mrf[,"Stim"] == "Ag85B",]
head(esat6_mrf)
head(ag85b_mrf)
# 8 models

# Two-Level Grouping:
# fit a separate model to each column
cd4_mrf <- cd4_mrf[cd4_mrf[,"Stim"]=="Ag85B" | cd4_mrf[,"Stim"]=="ESAT6",]
cd4_mrf$Stim <- as.factor(cd4_mrf$Stim)
summary(cd4_mrf)
dim(cd4_mrf)
# 4 models

# Three-Level Grouping:
# fit only one model
long <- melt(setDT(cd4_mrf), id.vars = c("Schedule","Dose", "timepnt", "Study", "QFT", "Group", "PID", "Day", "Stim"),
             variable.name = "Marker")
long$Marker <- as.factor(long$Marker)
head(long)
dim(long)





#---------------------------------------------------------------------
# Univariate EDA: Choosing an Appropriate Distribution for Y:
#---------------------------------------------------------------------



#---------------------
# Density Plots:
#---------------------

#----------------------------------------------------- Single-Level Grouping:

ggplot(esat6_mrf, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("CCR7+CD45RA+ Frequency")

ggplot(esat6_mrf, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("CCR7+CD45RA- Frequency")

ggplot(esat6_mrf, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("CCR7-CD45RA+ Frequency")

ggplot(esat6_mrf, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("CCR7-CD45RA- Frequency")

ggplot(ag85b_mrf, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("CCR7+CD45RA+ Frequency")

ggplot(ag85b_mrf, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("CCR7+CD45RA- Frequency")

ggplot(ag85b_mrf, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("CCR7-CD45RA+ Frequency")

ggplot(ag85b_mrf, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("CCR7-CD45RA- Frequency")

# all extremely zero-inflated, although R7-RA- less so
# plots are fairly similar for same memory combo but different stim, so maybe could get away with fitting a multi-level MME
# grouped by stimulus




#----------------------------------------------------- Two-Level Grouping:
ggplot(cd4_mrf, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("CCR7+CD45RA+ Frequency")

ggplot(cd4_mrf, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("CCR7+CD45RA- Frequency")

ggplot(cd4_mrf, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("CCR7-CD45RA+ Frequency")

ggplot(cd4_mrf, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("CCR7-CD45RA- Frequency")



#----------------------------------------------------- Three-Level Grouping:

ggplot(long, aes(x = value)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("Memory Response Frequency")


#-------------------------------------
# Classical Descriptive Statistics:
#-------------------------------------

#----------------------------------------------------- Single-Level Grouping:

set.seed(123)
descdist(esat6_mrf$R7..RA.., boot = 1000)
descdist(esat6_mrf$R7..RAc.., boot = 1000)  
descdist(esat6_mrf$R7c..RA.., boot = 1000)  
descdist(esat6_mrf$R7c..RAc.., boot = 1000)
descdist(ag85b_mrf$R7..RA.., boot = 1000)  
descdist(ag85b_mrf$R7..RAc.., boot = 1000)  
descdist(ag85b_mrf$R7c..RA.., boot = 1000)  
descdist(ag85b_mrf$R7c..RAc.., boot = 1000) 


#----------------------------------------------------- Two-Level Grouping:

descdist(cd4_mrf$R7..RA.., boot = 1000)    
descdist(cd4_mrf$R7..RAc.., boot = 1000)   
descdist(cd4_mrf$R7c..RA.., boot = 1000)   
descdist(cd4_mrf$R7c..RAc.., boot = 1000) 


#----------------------------------------------------- Three-Level Grouping:

descdist(long$value, boot = 1000)    



#-------------------------
# Fitting a Distribution
#-------------------------

#----------------------------------------------------------------- Single-Level Grouping:

# ------------------------------------------------- ESAT6 R7..RA..
f1.1 <- fitdist((esat6_mrf$R7..RA..+0.001), "lnorm")
f1.3 <- fitdist((esat6_mrf$R7..RA..+0.001), "gamma")
f1.4 <- fitdist((esat6_mrf$R7..RA..+0.001), "weibull")
f1.5 <- fitdist((10 + esat6_mrf$R7c..RA..), "SICHEL", start = list(mu = 10, sigma = 0.016))
model_1 <- cpglm(
  R7..RA.. ~1,
  data = esat6_mrf,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)



# DENSITY PLOTS
x1 <- seq(0,0.08, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(esat6_mrf, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA+ Frequency") +
  ylim(0,350)


# Need histogram to be Y+0.001 not Y:
df <- esat6_mrf
df$R7..RA.. <- df$R7..RA.. + 0.001

x1 <- seq(0,0.08, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA+ Frequency")+
  ylim(0,350)


x1 <- seq(0,0.08, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA+ Frequency")+
  ylim(0,350)

x1 <- seq(0,0.08, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA+ Frequency")+
  ylim(0,350)

m1 <-histDist(df$R7..RA.., "GG" , density=TRUE)
plot(m1)

m2 <-histDist(df$R7..RA.., "GIG" , density=TRUE)
plot(m2)

m3 <-histDist(df$R7..RA.., "GIG" , density=TRUE)
plot(m3)


# QQ PLOTS:

ggplot(esat6_mrf, aes(sample=R7..RA..)) +
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







#--------------------------------------------------------------- ESAT6 R7..RAc..
f1.1 <- fitdist((esat6_mrf$R7..RAc..+0.001), "lnorm")
f1.3 <- fitdist((esat6_mrf$R7..RAc..+0.001), "gamma")
f1.4 <- fitdist((esat6_mrf$R7..RAc..+0.001), "weibull")
model_1 <- cpglm(
  R7..RAc.. ~1,
  data = esat6_mrf,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)

# DENSITY PLOTS
x1 <- seq(0,0.7, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(esat6_mrf, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA- Frequency")+
  ylim(0,40)

# Need histogram to be Y+0.001 not Y:
df <- esat6_mrf
df$R7..RAc.. <- df$R7..RAc.. + 0.001

x1 <- seq(0,0.7, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA- Frequency")+
  ylim(0,40)


x1 <- seq(0,0.7, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA- Frequency")+
  ylim(0,40)


x1 <- seq(0,0.7, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA- Frequency")+
  ylim(0,40)



# QQ PLOTS:

ggplot(esat6_mrf, aes(sample=R7..RAc..)) +
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




#------------------------------------------------------------- ESAT6 R7c..RA..
f1.1 <- fitdist((esat6_mrf$R7c..RA..+0.001), "lnorm")
f1.3 <- fitdist((esat6_mrf$R7c..RA..+0.001), "gamma")
f1.4 <- fitdist((esat6_mrf$R7c..RA..+0.001), "weibull")
model_1 <- cpglm(
  R7c..RA.. ~1,
  data = esat6_mrf,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)

# DENSITY PLOTS
x1 <- seq(0,0.5, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(esat6_mrf, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,200)

# Need histogram to be Y+0.001 not Y:
df <- esat6_mrf
df$R7c..RA.. <- df$R7c..RA.. + 0.001

x1 <- seq(0,0.5, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,200)


x1 <- seq(0,0.5, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,200)

x1 <- seq(0,0.5, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,200)



# QQ PLOTS:

ggplot(esat6_mrf, aes(sample=R7c..RA..)) +
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


#### CROP OUT OUTLIER
# DENSITY PLOTS
x1 <- seq(0,0.08, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(esat6_mrf, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,450)+
  xlim(0,0.08)

# Need histogram to be Y+0.001 not Y:
df <- esat6_mrf
df$R7c..RA.. <- df$R7c..RA.. + 0.001

x1 <- seq(0,0.08, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,450)+
  xlim(0,0.08)


x1 <- seq(0,0.08, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,450)+
  xlim(0,0.08)

x1 <- seq(0,0.08, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,450)+
  xlim(0,0.08)



# QQ PLOTS:

ggplot(esat6_mrf, aes(sample=R7c..RA..)) +
  stat_qq(distribution = qtweedie, dparams = list(phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients)))+
  stat_qq_line(distribution = qtweedie, dparams = list(phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients)), color="red", linewidth=1)+
  theme_minimal(base_size = 27)+
  xlab("Theoretical quantiles")+
  ylab("Empirical quantiles")+
  ylim(0,0.1)

qqcomp(f1.1, addlegend=F, ylim=c(0,0.1), plotstyle = "ggplot") + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)


qqcomp(f1.3, addlegend=F, ylim=c(0,0.1), plotstyle = "ggplot") + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)

qqcomp(f1.4, addlegend=F, ylim=c(0,0.1), plotstyle = "ggplot") + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)



#----------------------------------------------------------- ESAT6 R7c..RAc..
f1.1 <- fitdist((esat6_mrf$R7c..RAc..+0.001), "lnorm")
f1.3 <- fitdist((esat6_mrf$R7c..RAc..+0.001), "gamma")
f1.4 <- fitdist((esat6_mrf$R7c..RAc..+0.001), "weibull")
model_1 <- cpglm(
  R7c..RAc.. ~1,
  data = esat6_mrf,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)


# DENSITY PLOTS
x1 <- seq(0,1.75, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(esat6_mrf, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA- Frequency")+
  ylim(0,13)

# Need histogram to be Y+0.001 not Y:
df <- esat6_mrf
df$R7c..RAc.. <- df$R7c..RAc.. + 0.001

x1 <- seq(0,1.75, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,13)


x1 <- seq(0,1.75, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,13)

x1 <- seq(0,1.75, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,13)



# QQ PLOTS:

ggplot(esat6_mrf, aes(sample=R7c..RAc..)) +
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







#--------------------------------------------------------------- Ag85B R7..RA..
f1.1 <- fitdist((ag85b_mrf$R7..RA..+0.001), "lnorm")
f1.3 <- fitdist((ag85b_mrf$R7..RA..+0.001), "gamma")
f1.4 <- fitdist((ag85b_mrf$R7..RA..+0.001), "weibull")
summary(f1.1)
summary(f1.3)
summary(f1.4)
model_1 <- cpglm(
  R7..RA.. ~1,
  data = ag85b_mrf,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)

# DENSITY PLOTS
x1 <- seq(0,0.12, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(ag85b_mrf, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA+ Frequency")+
  ylim(0,200)

# Need histogram to be Y+0.001 not Y:
df <- ag85b_mrf
df$R7..RA.. <- df$R7..RA.. + 0.001

x1 <- seq(0,0.12, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA+ Frequency")+
  ylim(0,200)


x1 <- seq(0,0.12, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA+ Frequency")+
  ylim(0,200)

x1 <- seq(0,0.12, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA+ Frequency")+
  ylim(0,200)



# QQ PLOTS:

ggplot(ag85b_mrf, aes(sample=R7..RA..)) +
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







# -------------------------------------------------------------- Ag85B R7..RAc..
f1.1 <- fitdist((ag85b_mrf$R7..RAc..+0.001), "lnorm")
f1.3 <- fitdist((ag85b_mrf$R7..RAc..+0.001), "gamma")
f1.4 <- fitdist((ag85b_mrf$R7..RAc..+0.001), "weibull")
summary(f1.1)
summary(f1.3)
summary(f1.4)
model_1 <- cpglm(
  R7..RAc.. ~1,
  data = ag85b_mrf,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)

# DENSITY PLOTS
x1 <- seq(0,1, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(ag85b_mrf, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA- Frequency")+
  ylim(0,20)

# Need histogram to be Y+0.001 not Y:
df <- ag85b_mrf
df$R7..RAc.. <- df$R7..RAc.. + 0.001

x1 <- seq(0,1, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA- Frequency")+
  ylim(0,20)


x1 <- seq(0,1, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA- Frequency")+
  ylim(0,20)

x1 <- seq(0,1, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA- Frequency")+
  ylim(0,20)



# QQ PLOTS:

ggplot(ag85b_mrf, aes(sample=R7..RAc..)) +
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





# ------------------------------------------------------------- Ag85B R7c..RA..
f1.1 <- fitdist((ag85b_mrf$R7c..RA..+0.001), "lnorm")
f1.3 <- fitdist((ag85b_mrf$R7c..RA..+0.001), "gamma")
f1.4 <- fitdist((ag85b_mrf$R7c..RA..+0.001), "weibull")
summary(f1.1)
summary(f1.3)
summary(f1.4)
model_1 <- cpglm(
  R7c..RA.. ~1,
  data = ag85b_mrf,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)

# DENSITY PLOTS
x1 <- seq(0,0.6, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(ag85b_mrf, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,100)

# Need histogram to be Y+0.001 not Y:
df <- ag85b_mrf
df$R7c..RA.. <- df$R7c..RA.. + 0.001

x1 <- seq(0,0.6, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,100)


x1 <- seq(0,0.6, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,100)

x1 <- seq(0,0.6, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,100)



# QQ PLOTS:

ggplot(ag85b_mrf, aes(sample=R7c..RA..)) +
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


### CROP OUTLIER:

ggplot(ag85b_mrf, aes(sample=R7c..RA..)) +
  stat_qq(distribution = qtweedie, dparams = list(phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients)))+
  stat_qq_line(distribution = qtweedie, dparams = list(phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients)), color="red", linewidth=1)+
  theme_minimal(base_size = 27)+
  xlab("Theoretical quantiles")+
  ylab("Empirical quantiles")+
  ylim(0,0.11)

qqcomp(f1.1, addlegend=F, plotstyle = "ggplot", ylim = c(0,0.11)) + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)


qqcomp(f1.3, addlegend=F, plotstyle = "ggplot", ylim = c(0,0.11)) + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)

qqcomp(f1.4, addlegend=F, plotstyle = "ggplot", ylim = c(0,0.11)) + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)




#------------------------------------------------------------- Ag85B R7c..RAc..
f1.1 <- fitdist((ag85b_mrf$R7c..RAc..+0.001), "lnorm")
f1.3 <- fitdist((ag85b_mrf$R7c..RAc..+0.001), "gamma")
f1.4 <- fitdist((ag85b_mrf$R7c..RAc..+0.001), "weibull")
summary(f1.1)
summary(f1.3)
summary(f1.4)
model_1 <- cpglm(
  R7c..RAc.. ~1,
  data = ag85b_mrf,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)

# DENSITY PLOTS
x1 <- seq(0,1.5, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(ag85b_mrf, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA- Frequency")+
  ylim(0,12)

# Need histogram to be Y+0.001 not Y:
df <- ag85b_mrf
df$R7c..RAc.. <- df$R7c..RAc.. + 0.001

x1 <- seq(0,1.5, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA- Frequency")+
  ylim(0,12)


x1 <- seq(0,1.5, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA- Frequency")+
  ylim(0,12)

x1 <- seq(0,1.5, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA- Frequency")+
  ylim(0,12)



# QQ PLOTS:

ggplot(ag85b_mrf, aes(sample=R7c..RAc..)) +
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




#-------------------------------------------------------------------------------
#----------------------------------------------------------- Two-Level Grouping:
#-------------------------------------------------------------------------------

#---------------------------------------------------------------------- R7..RA..
f1.1 <- fitdist((cd4_mrf$R7..RA..+0.001), "lnorm")
f1.3 <- fitdist((cd4_mrf$R7..RA..+0.001), "gamma")
f1.4 <- fitdist((cd4_mrf$R7..RA..+0.001), "weibull")
model_1 <- cpglm(
  R7..RA.. ~1,
  data = cd4_mrf,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)

# DENSITY PLOTS
x1 <- seq(0,0.12, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(cd4_mrf, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA+ Frequency")+
  ylim(0,225)


# Need histogram to be Y+0.001 not Y:
df <- cd4_mrf
df$R7..RA.. <- df$R7..RA.. + 0.001

x1 <- seq(0,0.12, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA+ Frequency")+
  ylim(0,225)


x1 <- seq(0,0.12, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA+ Frequency")+
  ylim(0,225)

x1 <- seq(0,0.12, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA+ Frequency")+
  ylim(0,225)



# QQ PLOTS:

ggplot(cd4_mrf, aes(sample=R7..RA..)) +
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




#---------------------------------------------------------------------- R7..RAc..
f1.1 <- fitdist((cd4_mrf$R7..RAc..+0.001), "lnorm")
f1.3 <- fitdist((cd4_mrf$R7..RAc..+0.001), "gamma")
f1.4 <- fitdist((cd4_mrf$R7..RAc..+0.001), "weibull")
model_1 <- cpglm(
  R7..RAc.. ~1,
  data = cd4_mrf,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)

# DENSITY PLOTS
x1 <- seq(0,1, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(cd4_mrf, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA- Frequency")+
  ylim(0,30)

# Need histogram to be Y+0.001 not Y:
df <- cd4_mrf
df$R7..RAc.. <- df$R7..RAc.. + 0.001

x1 <- seq(0,1, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA- Frequency")+
  ylim(0,30)

x1 <- seq(0,1, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA- Frequency")+
  ylim(0,30)

x1 <- seq(0,1, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7+CD45RA- Frequency")+
  ylim(0,30)



# QQ PLOTS:

ggplot(cd4_mrf, aes(sample=R7..RAc..)) +
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


#---------------------------------------------------------------------- R7c..RA..
f1.1 <- fitdist((cd4_mrf$R7c..RA..+0.001), "lnorm")
f1.3 <- fitdist((cd4_mrf$R7c..RA..+0.001), "gamma")
f1.4 <- fitdist((cd4_mrf$R7c..RA..+0.001), "weibull")
model_1 <- cpglm(
  R7c..RA.. ~1,
  data = cd4_mrf,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)

# DENSITY PLOTS
x1 <- seq(0,0.6, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(cd4_mrf, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,100)

# Need histogram to be Y+0.001 not Y:
df <- cd4_mrf
df$R7c..RA.. <- df$R7c..RA.. + 0.001

x1 <- seq(0,0.6, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,100)


x1 <- seq(0,0.6, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,100)

x1 <- seq(0,0.6, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7c..RA..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA+ Frequency")+
  ylim(0,100)



# QQ PLOTS:

ggplot(cd4_mrf, aes(sample=R7c..RA..)) +
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


# CLOES UP:
ggplot(cd4_mrf, aes(sample=R7c..RA..)) +
  stat_qq(distribution = qtweedie, dparams = list(phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients)))+
  stat_qq_line(distribution = qtweedie, dparams = list(phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients)), color="red", linewidth=1)+
  theme_minimal(base_size = 27)+
  xlab("Theoretical quantiles")+
  ylab("Empirical quantiles")+
  ylim(0,0.11)
qqcomp(f1.1, addlegend=F, plotstyle = "ggplot", ylim=c(0,0.11)) + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)
qqcomp(f1.3, addlegend=F, plotstyle = "ggplot", ylim=c(0,0.11)) + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)
qqcomp(f1.4, addlegend=F, plotstyle = "ggplot", ylim=c(0,0.11)) + 
  theme_minimal(base_size = 27) + 
  geom_point(shape=19, size=2, fill="black", color="black")+
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)



#---------------------------------------------------------------------- R7c..RAc..
f1.1 <- fitdist((cd4_mrf$R7c..RAc..+0.001), "lnorm")
f1.3 <- fitdist((cd4_mrf$R7c..RAc..+0.001), "gamma")
f1.4 <- fitdist((cd4_mrf$R7c..RAc..+0.001), "weibull")
model_1 <- cpglm(
  R7c..RAc.. ~1,
  data = cd4_mrf,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)


# DENSITY PLOTS
x1 <- seq(0,1.75, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(cd4_mrf, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA- Frequency")+
  ylim(0,12)

# Need histogram to be Y+0.001 not Y:
df <- cd4_mrf
df$R7c..RAc.. <- df$R7c..RAc.. + 0.001

x1 <- seq(0,1.75, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA- Frequency")+
  ylim(0,12)


x1 <- seq(0,1.75, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA- Frequency")+
  ylim(0,12)

x1 <- seq(0,1.75, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = R7c..RAc..)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("CCR7-CD45RA- Frequency")+
  ylim(0,12)



# QQ PLOTS:
ggplot(cd4_mrf, aes(sample=R7c..RAc..)) +
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




#-------------------------------------------------------------------------------
#----------------------------------------------------- Three-Level Grouping:
#-------------------------------------------------------------------------------

f1.1 <- fitdist((long$value+0.001), "lnorm")
f1.3 <- fitdist((long$value+0.001), "gamma")
f1.4 <- fitdist((long$value+0.001), "weibull")
model_1 <- cpglm(
  value ~1,
  data = long,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)




# DENSITY PLOTS
x1 <- seq(0,1.75, length=50000)
y1 <- dtweedie(x1, phi=model_1$phi, power=model_1$p, mu=exp(model_1$coefficients))
ggplot(long, aes(x = value)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("Frequency")+
  ylim(0,30)

# Need histogram to be Y+0.001 not Y:
df <- long
df$value <- df$value + 0.001

x1 <- seq(0,1.75, length=50000)
y1 <- dlnorm(x1, meanlog=f1.1$estimate[1], sdlog=f1.1$estimate[2])
ggplot(df, aes(x = value)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("Frequency")+
  ylim(0,30)


x1 <- seq(0,1.75, length=50000)
y1 <- dgamma(x1, shape=f1.3$estimate[1], rate=f1.3$estimate[2])
ggplot(df, aes(x = value)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("Frequency")+
  ylim(0,30)

x1 <- seq(0,1.75, length=50000)
y1 <- dweibull(x1, shape=f1.4$estimate[1], scale=f1.4$estimate[2])
ggplot(df, aes(x = value)) + 
  theme_minimal(base_size = 23)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_line(data = data.frame(x1,y1), aes(x = x1 , y = y1),lwd = 1, colour = 4)+
  ylab("Density") +
  geom_area(data = data.frame(x1,y1), aes(y = y1, x=x1), color=4, fill= 4, alpha=0.25)+
  xlab("Frequency")+
  ylim(0,30)



# QQ PLOTS:
ggplot(long, aes(sample=value)) +
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














#---------------------------------------------------------------------
# Bivariate EDA: Investigating Response Profiles:
#---------------------------------------------------------------------






#-------------------------------------
# Interaction Plot 
#-------------------------------------


#----------------------------------------------------- Single-Level Grouping:

long$Marker <- as.factor(long$Marker)
ggplot(data=long, aes(x = Day, y = value, group = PID, colour = Marker)) +
  geom_line()+
  ylab("Frequency")+
  xlab("Day") +
  theme_minimal(base_size = 27)
# too busy

# Use Medians:
ggplot() +
  aes(x = long$Day, color = long$Marker, group = long$Marker, y = long$value) +
  stat_summary(fun = median, geom = "point", size=4) +
  stat_summary(fun = median, geom = "line", linewidth=1.5) +
  ylab("Frequency")+
  xlab("Day") +
  theme_minimal(base_size = 27) 
# this looks really confusing because we have combined all no.s of administrations and plotted 
# by day rather than time point - misleading. Either need to separate by no. of admin or plot by timepnt

# Plot by timepnt:
ggplot() +
  aes(x = long$timepnt, color = long$Marker, group = long$Marker, y = long$value) +
  stat_summary(fun = median, geom = "point", size=4) +
  stat_summary(fun = median, geom = "line", linewidth=1.5) +
  ylab("Frequency")+
  xlab("Day") +
  theme_minimal(base_size = 27) 
# medians are all zero, try mean

ggplot() +
  aes(x = long$timepnt, color = long$Marker, group = long$Marker, y = long$value) +
  stat_summary(fun = mean, geom = "point", size=4) +
  stat_summary(fun = mean, geom = "line", linewidth=1.5) +
  ylab("Frequency")+
  xlab("Day") +
  theme_minimal(base_size = 27) 

corr <- cor(esat6_mrf[,c(6,8,10,11,12,13)])
ggcorrplot(corr)








