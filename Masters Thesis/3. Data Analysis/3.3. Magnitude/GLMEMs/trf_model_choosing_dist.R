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
cd4_pos <- data.frame(cd4_bstrf[,1:7], total = rowSums(cd4_bstrf[,8:67]))
  
# remove cd4 counts column
cd4_pos = cd4_pos[,-7]

#---------------------
# Add Time Category
#---------------------
cd4_pos <- data.frame(timepnt = NA, cd4_pos)

cd4_pos$timepnt[which(cd4_pos$Day==0)] <- "baseline"
cd4_pos$timepnt[which(cd4_pos$Day==14)] <- "peak"
cd4_pos$timepnt[which(cd4_pos$Day==70)] <- "peak"
cd4_pos$timepnt[which(cd4_pos$Day==126)] <- "peak"
cd4_pos$timepnt[which(cd4_pos$Day==224)] <- "memory"
cd4_pos$timepnt[which(cd4_pos$Day==292)] <- "memory"
cd4_pos$timepnt[which(cd4_pos$Day==210)] <- "memory"

which(is.na(cd4_pos$timepnt) == T)

cd4_pos$timepnt <- as.factor(cd4_pos$timepnt)
cd4_pos$timepnt <- factor(cd4_pos$timepnt, levels=c("baseline", "peak", "memory"))
levels(cd4_pos$timepnt)


#---------------------
# Add Dose
#---------------------
# Need to add in dose and schedule information
cd4_pos <- data.frame(Dose = NA,cd4_pos)

cd4_pos[intersect(which(cd4_pos$Study == "H56-035"),which(cd4_pos$Group == 1)),1] <- 5
cd4_pos[intersect(which(cd4_pos$Study == "H56-035"),which(cd4_pos$Group == 2)),1] <- 5
cd4_pos[intersect(which(cd4_pos$Study == "H56-035"),which(cd4_pos$Group == 3)),1] <- 5
cd4_pos[intersect(which(cd4_pos$Study == "H56-035"),which(cd4_pos$Group == 4)),1] <- 5

cd4_pos[intersect(which(cd4_pos$Study == "H1"),which(cd4_pos$Group == 1)),1] <- 15
cd4_pos[intersect(which(cd4_pos$Study == "H1"),which(cd4_pos$Group == 2)),1] <- 50
cd4_pos[intersect(which(cd4_pos$Study == "H1"),which(cd4_pos$Group == 3)),1] <- 15
cd4_pos[intersect(which(cd4_pos$Study == "H1"),which(cd4_pos$Group == 4)),1] <- 50

cd4_pos[intersect(which(cd4_pos$Study == "H56-032"),which(cd4_pos$Group == 1)),1] <- 50
cd4_pos[intersect(which(cd4_pos$Study == "H56-032"),which(cd4_pos$Group == 2)),1] <- 15
cd4_pos[intersect(which(cd4_pos$Study == "H56-032"),which(cd4_pos$Group == 3)),1] <- 50

# check
cd4_pos[which(is.na(cd4_pos$Dose)==T),]
# convert to factor
cd4_pos$Dose <- as.factor(cd4_pos$Dose)


#---------------------
# Add Schedule
#---------------------
# Need to add in dose and schedule information
cd4_pos <- data.frame(Schedule = NA,cd4_pos)

cd4_pos[intersect(which(cd4_pos$Study == "H56-035"),which(cd4_pos$Group == 1)),1] <- 2
cd4_pos[intersect(which(cd4_pos$Study == "H56-035"),which(cd4_pos$Group == 2)),1] <- 3
cd4_pos[intersect(which(cd4_pos$Study == "H56-035"),which(cd4_pos$Group == 3)),1] <- 2
cd4_pos[intersect(which(cd4_pos$Study == "H56-035"),which(cd4_pos$Group == 4)),1] <- 3

cd4_pos[intersect(which(cd4_pos$Study == "H1"),which(cd4_pos$Group == 1)),1] <- 2
cd4_pos[intersect(which(cd4_pos$Study == "H1"),which(cd4_pos$Group == 2)),1] <- 2
cd4_pos[intersect(which(cd4_pos$Study == "H1"),which(cd4_pos$Group == 3)),1] <- 1
cd4_pos[intersect(which(cd4_pos$Study == "H1"),which(cd4_pos$Group == 4)),1] <- 1

cd4_pos[intersect(which(cd4_pos$Study == "H56-032"),which(cd4_pos$Group == 1)),1] <- 3
cd4_pos[intersect(which(cd4_pos$Study == "H56-032"),which(cd4_pos$Group == 2)),1] <- 3
cd4_pos[intersect(which(cd4_pos$Study == "H56-032"),which(cd4_pos$Group == 3)),1] <- 3

# check 
cd4_pos[which(is.na(cd4_pos$Schedule)==T),]
# convert to factor
cd4_pos$Schedule <- as.factor(cd4_pos$Schedule)
head(cd4_pos)


#-----------------------------------------------------------
# Prepare Data
#-----------------------------------------------------------


# 2 levels:
tvr <- cd4_pos[cd4_pos[,"Stim"] == "ESAT6" | cd4_pos[,"Stim"] == "Ag85B",]


# 1 level:
trf_esat6 = tvr[tvr[,"Stim"]== "ESAT6", ]
trf_ag85b = tvr[tvr[,"Stim"]== "Ag85B", ]


#---------------------------------------------------------------------
# Univariate EDA: Choosing an Appropriate Distribution for Y:
#---------------------------------------------------------------------



#---------------------
# Density Plots:
#---------------------


#----------------------------------------------------- Two-Level Grouping:
ggplot(tvr, aes(x = total)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("Total Response Frequency")


#----------------------------------------------------- One-Level Grouping:
ggplot(trf_esat6, aes(x = total)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("Total Response Frequency")

ggplot(trf_ag85b, aes(x = total)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("Total Response Frequency")



#-------------------------
# Fitting a Distribution
#-------------------------


#-------------------------------------------------------------------------------
#----------------------------------------------------------- Two-Level Grouping: 
#-------------------------------------------------------------------------------
f1.1 <- fitdist((tvr$total+0.001), "lnorm")
f1.3 <- fitdist((tvr$total+0.001), "gamma")
f1.4 <- fitdist((tvr$total+0.001), "weibull")
model_1 <- cpglm(
  total ~1,
  data = tvr,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)


# QQ PLOTS:
ggplot(tvr, aes(sample=total)) +
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
#----------------------------------------------------------- One-Level Grouping: 
#-------------------------------------------------------------------------------

#---------------------------------------------------------------------- ESAT6
f1.1 <- fitdist((trf_esat6$total+0.001), "lnorm")
f1.3 <- fitdist((trf_esat6$total+0.001), "gamma")
f1.4 <- fitdist((trf_esat6$total+0.001), "weibull")
model_1 <- cpglm(
  total ~1,
  data = trf_esat6,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)



# QQ PLOTS:
ggplot(trf_esat6, aes(sample=total)) +
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
f1.1 <- fitdist((trf_ag85b$total+0.001), "lnorm")
f1.3 <- fitdist((trf_ag85b$total+0.001), "gamma")
f1.4 <- fitdist((trf_ag85b$total+0.001), "weibull")
model_1 <- cpglm(
  total ~1,
  data = trf_ag85b,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)


# QQ PLOTS:
ggplot(trf_ag85b, aes(sample=total)) +
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




