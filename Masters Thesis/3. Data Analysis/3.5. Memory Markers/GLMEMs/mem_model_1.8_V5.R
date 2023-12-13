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
library(lme4) # glmer function
library(optimx)  # More optimisers to use when fitting a glmer
library(cplm)
library(lmerTest) # to get p-values!
library(glmmTMB)
library(gamlss)
library(DHARMa) # get quantile residuals for glmmTMB
library(RColorBrewer)
library(ppcc) #Filliben correlation coefficient

#---------------------
# Import Data
#---------------------

cd4_counts <- read.csv("clean_data_with_demographic_info.csv", check.names = F)
str(cd4_counts)



#-------------------------------------------
# Calculate Memory Response Frequency (MRF)
#-------------------------------------------

# Only keep counts for cells expressing at least one of INFg, IL2, TNF or L17 (i.e. cytokine+)
cd4_counts <- cd4_counts[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd4_counts), fixed = TRUE)==T)]
str(cd4_counts)


# Sum over cytokines:
cd4_counts$R7..RA..   <- rowSums(cd4_counts[,14:28])
cd4_counts$R7..RAc..  <- rowSums(cd4_counts[,29:43])
cd4_counts$R7c..RA..  <- rowSums(cd4_counts[,44:58])
cd4_counts$R7c..RAc.. <- rowSums(cd4_counts[,59:73])


# Convert to Frequencies (%):
cd4_counts[,14:77] <- (cd4_counts[,14:77]/cd4_counts$CD4count)*100

# Check:
summary(rowSums(cd4_counts[,14:73]))
summary(rowSums(cd4_counts[,74:77]))
# Both are the same

# Remove unwanted columns:
cd4_counts <- cd4_counts[,-c(13:73)]

# Background Subtract
cd4_mrf <- data.frame()
for(i in 1:length(unique(cd4_counts$PID))){
  pid  <- unique(cd4_counts$PID)[i]
  dat  <- cd4_counts[which(cd4_counts$PID == pid),]
  days <- unique(dat$Day)
  for (j in 1:length(days)) {
    day <- days[j]
    df  <- dat[which(dat$Day == day),]
    df[-which(df$Stim == "UNS"),13:16]  <- df[-which(df$Stim == "UNS"),13:16] - df[rep(which(df$Stim == "UNS"), (length(df$Stim)-1)),13:16]
    cd4_mrf <- rbind(cd4_mrf,df)
  }
}
summary(cd4_mrf)
which(is.na(cd4_mrf)==T) # only th 12 observations (for the same pid) where sex is missing


# Deal with Negative Frequencies:
for (i in 13:ncol(cd4_mrf)) {
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
sapply(sign(cd4_mrf[,-(1:12)]), table)


# Remove UNS and PHA
cd4_mrf <- cd4_mrf[cd4_mrf[,"Stim"]=="Ag85B" | cd4_mrf[,"Stim"]=="ESAT6",]
cd4_mrf$Stim <- as.factor(cd4_mrf$Stim)
summary(cd4_mrf)



#---------------------
# Remove PID 411
#---------------------


# Remove 411 because only baseline available
cd4_mrf <- cd4_mrf[-which(cd4_mrf$PID == 411),]


# Set Factors:
cd4_mrf$ethnicity <- as.factor(cd4_mrf$ethnicity)
cd4_mrf$sex       <- as.factor(cd4_mrf$sex)
cd4_mrf$timepnt   <- as.factor(cd4_mrf$timepnt)
cd4_mrf$Study     <- as.factor(cd4_mrf$Study)
cd4_mrf$QFT       <- as.factor(cd4_mrf$QFT)
cd4_mrf$PID       <- as.factor(cd4_mrf$PID)
cd4_mrf$Stim      <- as.factor(cd4_mrf$Stim)


cd4_mrf$timepnt <- as.factor(cd4_mrf$timepnt)
cd4_mrf$timepnt <- factor(cd4_mrf$timepnt, levels=c("baseline", "peak", "memory"))
levels(cd4_mrf$timepnt)



#---------------------
# Prepare Data
#---------------------

# Subset data by stim:
esat6_mrf <- cd4_mrf[cd4_mrf[,"Stim"] == "ESAT6",]
ag85b_mrf <- cd4_mrf[cd4_mrf[,"Stim"] == "Ag85B",]
head(esat6_mrf)
head(ag85b_mrf)



#---------------------------------------------------------
# Choose Fixed Effects: Step-wise Selection
#---------------------------------------------------------
df <- ag85b_mrf
df$R7c..RAc.. <- df$R7c..RAc..+0.001

# Keep Only complete cases:
df <- df[,-c(13,14,15)]
df1<-df[complete.cases(df),]

# Base Model:
base_model <-gamlss(R7c..RAc.. ~ 
                      QFT + 
                      Schedule + 
                      Dose +
                      timepnt +
                      re(random=~1|PID, method = "REML"),  
                    data=df1, 
                    family=WEI(mu.link=log, sigma.link = identity)) # NB identity link because lmm
summary(base_model)



# Step-wise Selection:
mod1 <- stepGAIC(base_model, 
                 scope = list(lower= ~ QFT + Schedule + Dose + timepnt + re(random=~1|PID, method = "REML"),
                              upper = ~ QFT + Schedule + Dose + timepnt + age + ethnicity + sex + Study + Day + re(random=~1|PID, method = "REML")),
                 k=round(log(length(df1$R7c..RAc..)), 2))
summary(mod1)



# Since Sex was not included in our model, refit model using observation with missing sex data included 
# (note that the model will not be comparebale in terms of AIC and BIC to the previous model)
df <- ag85b_mrf
df$R7c..RAc.. <- df$R7c..RAc..+0.001

# Keep Only complete cases:
df <- df[,-c(3,13,14,15)]
df1<-df[complete.cases(df),]


# Base Model:
mod1 <-gamlss(R7c..RAc.. ~ 
                QFT + 
                Schedule + 
                Dose +
                timepnt +
                re(random=~1|PID, method = "REML"), 
              data=df1, 
              family=WEI(mu.link=log, sigma.link = identity)) # NB idenity link because lmm


#---------------------------------------------------------
# Choose Interactions: Step-wise Selection
#---------------------------------------------------------

# Step-wise Selection:
mod2 <- stepGAIC(mod1, 
                 scope = list(lower= ~ QFT + Schedule + Dose + timepnt + re(random=~1|PID, method = "REML"),
                              upper = ~ Schedule*Dose + timepnt*Schedule + QFT*Dose + QFT*Schedule + QFT*timepnt + Dose*timepnt + re(random=~1|PID, method = "REML")),
                 k=round(log(length(df1$R7c..RAc..)), 2))
summary(mod2) # Choose model2

#---------------------------------------------------------
# Choose Sigma Model: Step-wise Selection
#---------------------------------------------------------

# Step-wise Selection for Sigma: 
mod3.0 <-gamlss(R7c..RAc.. ~ 
                  QFT*timepnt + 
                  Schedule + 
                  Dose +
                  re(random=~1|PID, method = "REML"),
                data=df1, 
                sigma.formula = ~ 1, 
                family=WEI(mu.link=log, sigma.link = log)) # log link for sigma because adding formula

mod3 <- stepGAIC(mod3.0, 
                 scope = list(lower= ~ 1,
                              upper = ~ QFT + Schedule + Dose), # For some reason timepnt results in algorithm not converging, so i'm just going to try other variables
                 k=round(log(length(df1$R7c..RAc..)), 2), 
                 what = "sigma")
summary(mod3) # No variables selected




# FINAL MODEL:
summary(mod2)
intervals(getSmo(mod2))
GAIC(mod2, k=round(log(length(df1$R7c..RAc..)), 2))


#-----------------------------------
# Residuals by category:
#-----------------------------------
# Residuals vs Timepnt:
df2 <- data.frame(residuals = resid(mod2),  time = df1$timepnt)
ggplot(df2, aes(y=residuals, x=time)) +
  geom_boxplot(alpha=1, size =1, aes(fill = time), color ="black") +
  theme_minimal()+
  theme(text = element_text(size = 30), legend.position = "none") +
  ylab("Normalized quantile residuals") +
  xlab(NULL)+
  scale_fill_manual(values = brewer.pal(n = 4, name = "BuGn"))
# Residuals vs QFT:
df2 <- data.frame(residuals = resid(mod2),  QFT = df1$QFT)
ggplot(df2, aes(y=residuals, x=QFT)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color ="black") +
  theme_minimal()+
  theme(text = element_text(size = 30), legend.position = "none") +
  ylab("Normalized quantile residuals") +
  xlab(NULL)+
  scale_fill_manual(values=c("#FF0000", "blue"))+ 
  scale_x_discrete(labels=c("QFT neg" = "QFT-", "QFT pos" = "QFT+"))
# Residuals vs admin:
df2 <- data.frame(residuals = resid(mod2), Schedule = df1$Schedule)
df2$Schedule <- as.factor(df2$Schedule)
ggplot(df2, aes(y=residuals, x=Schedule)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Schedule), color ="black") +
  theme_minimal()+
  theme(text = element_text(size = 30), legend.position = "none") +
  ylab("Normalized quantile residuals") +
  xlab("Number of Administrations")+
  scale_fill_manual(values = brewer.pal(n = 3, name = "BuGn"))
# Residuals vs concentration:
df2 <- data.frame(residuals = resid(mod2),  Dose = df1$Dose)
df2$Dose <- as.factor(df2$Dose)
ggplot(df2, aes(y=residuals, x=Dose)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Dose), color ="black") +
  theme_minimal()+
  theme(text = element_text(size = 30), legend.position = "none") +
  ylab("Normalized quantile residuals") +
  xlab(expression(paste("Concentration (", mu, "g)")))+
  scale_fill_manual(values = brewer.pal(n = 3, name = "BuGn"))
# QFT*timepnt:
df2 <- data.frame(residuals = resid(mod2),  time = df1$timepnt, QFT=df1$QFT)
ggplot(df2, aes(y=residuals, x=time)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color ="black") +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("Normalized quantile residuals") +
  xlab(NULL)+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")



#-----------------------------------------
# Check Normality of Residuals 
#-----------------------------------------

df2 <- data.frame(residuals = resid(mod2))
ggplot(df2, aes(sample=residuals))+
  stat_qq()+
  theme_minimal(base_size = 27) + 
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)+
  ylab("Sample Quantiles")+
  xlab("Normal Quantiles")


#---------------------------------------------
# Check Marginal Normality of Random Errors
#---------------------------------------------

df2 <- data.frame(randeff = ranef(getSmo(mod2)))
ggplot(df2, aes(sample=X.Intercept.))+
  stat_qq()+
  theme_minimal(base_size = 27) + 
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)+
  ylab("Sample Quantiles")+
  xlab("Normal Quantiles")


#---------------------------------------------
# Check Overall Model fit:
#---------------------------------------------

# Residuals vs fitted:
df2 <- data.frame(residuals = resid(mod2), fitted = fitted(mod2))
ggplot(data = df2, aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Normalized quantile residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)
