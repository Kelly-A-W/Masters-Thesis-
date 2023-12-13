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



#-------------------------------------
# Calculate TRF
#-------------------------------------

# Only keep counts for cells expressing at least one of INFg, IL2, TNF or L17
cd4_trf <- cd4_counts[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd4_counts), fixed = TRUE)==T)]
str(cd4_trf)

# Convert to Frequencies
cd4_trf[,14:73] <- (cd4_trf[,14:73]/cd4_trf$CD4count)*100

# Background Subtract
cd4_bstrf <- data.frame()
for(i in 1:length(unique(cd4_trf$PID))){
  pid <- unique(cd4_trf$PID)[i]
  dat <- cd4_trf[which(cd4_trf$PID == pid),]
  days <- unique(dat$Day)
  for (j in 1:length(days)) {
    day <- days[j]
    df  <- dat[which(dat$Day == day),]
    df[-which(df$Stim == "UNS"),14:73]  <- df[-which(df$Stim == "UNS"),14:73] - df[rep(which(df$Stim == "UNS"), (length(df$Stim)-1)),14:73]
    cd4_bstrf <- rbind(cd4_bstrf,df)
  }
}
summary(cd4_bstrf)
which(is.na(cd4_bstrf)==T) # only the 12 observations (for the same pid) where sex is missing
which(is.na(cd4_bstrf$sex)==T)


# Deal with Negative Frequencies:
for (i in 14:ncol(cd4_bstrf)) {
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
sapply(sign(cd4_bstrf[,-(1:13)]), table)

# Sum across rows
cd4_pos <- data.frame(cd4_bstrf[,c(1,2,3,4,5,6,7,8,10,12)], total = rowSums(cd4_bstrf[,14:73]))
str(cd4_pos)

#---------------------
# Remove PID 411
#---------------------


# Remove 411 because only baseline available
cd4_pos <- cd4_pos[-which(cd4_pos$PID == 411),]


# Set Factors:
cd4_pos$ethnicity <- as.factor(cd4_pos$ethnicity)
cd4_pos$sex       <- as.factor(cd4_pos$sex)
cd4_pos$timepnt   <- as.factor(cd4_pos$timepnt)
cd4_pos$Study     <- as.factor(cd4_pos$Study)
cd4_pos$QFT       <- as.factor(cd4_pos$QFT)
cd4_pos$PID       <- as.factor(cd4_pos$PID)
cd4_pos$Stim      <- as.factor(cd4_pos$Stim)


cd4_pos$timepnt <- as.factor(cd4_pos$timepnt)
cd4_pos$timepnt <- factor(cd4_pos$timepnt, levels=c("baseline", "peak", "memory"))
levels(cd4_pos$timepnt)



#-----------------------------------------------------------
# Prepare Data
#-----------------------------------------------------------


# 2 levels:
tvr <- cd4_pos[cd4_pos[,"Stim"] == "ESAT6" | cd4_pos[,"Stim"] == "Ag85B",]


# 1 level:
trf_esat6 = tvr[tvr[,"Stim"]== "ESAT6", ]
trf_ag85b = tvr[tvr[,"Stim"]== "Ag85B", ]




#---------------------------------------------------------
# Choose Fixed Effects: Step-wise Selection
#---------------------------------------------------------

# Convert data to log scale:
df <- trf_esat6
df$total <- df$total + 0.001 

# Keep Only complete cases:
df1<-df[complete.cases(df),]


# Base Model:
base_model <-gamlss(total ~ 
                      QFT + 
                      Schedule + 
                      Dose +
                      timepnt +
                      re(random=~1|PID, method = "REML"),  
                    data=df1, 
                    family=WEI(mu.link=log, sigma.link = identity)) # NB identity link because lmm



# Step-wise Selection:
mod1 <- stepGAIC(base_model, 
                 scope = list(lower= ~ QFT + Schedule + Dose + timepnt + re(random=~1|PID, method = "REML"),
                              upper = ~ QFT + Schedule + Dose + timepnt + age + ethnicity + sex + Study + Day + re(random=~1|PID, method = "REML")),
                 k=round(log(length(df1$total)), 2))
summary(mod1)  
# no fixed effects added


# Since Sex was not included in our model, refit model using observation with missing sex data included 
# (note that the model will not be comparebale in terms of AIC and BIC to the previous model)
# Convert data to log scale:
df <- trf_esat6
df$total <- df$total + 0.001  # fitting a lognormal distr

# Keep Only complete cases:
df <- df[,-3]
df1<-df[complete.cases(df),]


# Base Model:
mod1 <-gamlss(total ~ 
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
                 k=round(log(length(df1$total)), 2))
summary(mod2) # interactions added between schedule and timepont and QFT and timepnt



#----------------------------
# 3-way Interactions
#----------------------------

mod3 <-gamlss(total ~ 
                  QFT*Schedule*timepnt + 
                  Dose +
                  timepnt +
                  re(random=~1|PID, method = "REML"),
                data=df1, 
                family=WEI(mu.link=log, sigma.link = identity))


GAIC(mod2, k=round(log(length(df1$total)), 2)) # choose mod2
GAIC(mod3, k=round(log(length(df1$total)), 2)) 



#-------------------------
# Model Sigma:
#-------------------------

# Step-wise Selection for Sigma: 
mod4.0 <-gamlss(total ~ 
                  QFT*Schedule + 
                  Schedule*timepnt + 
                  Dose +
                  timepnt +
                  re(random=~1|PID, method = "REML"),
                data=df1, 
                sigma.formula = ~ 1, 
                family=WEI(mu.link=log, sigma.link = log)) # log link for sigma because adding formula

mod4 <- stepGAIC(mod4.0, 
                 scope = list(lower= ~ 1,
                              upper = ~ Schedule + Dose + QFT),
                 k=round(log(length(df1$total)), 2), 
                 what = "sigma")
summary(mod4) # schedule added


# Compare to old model:
GAIC(mod4.0, k=round(log(length(df1$total)), 2))
GAIC(mod4, k=round(log(length(df1$total)), 2)) # choose mod4




# FINAL MODEL:
summary(mod4)
getSmo(mod4) 
GAIC(mod4, k=round(log(length(df1$total)), 2))
intervals(getSmo(mod4))








#----------------------------
# Interaction Plots 
#---------------------------
qft_pos_1 = c(0, 0.9921135 + 0.0013404*1, -0.1351799 + 0.3145422*1) + 2.6294794 - 0.0093157*1 -0.5589063*1
qft_neg_1 = c(0, 0.9921135 + 0.0013404*1, -0.1351799 + 0.3145422*1) - 0.0093157*1
qft_pos_2 = c(0, 0.9921135 + 0.0013404*2, -0.1351799 + 0.3145422*2) + 2.6294794 - 0.0093157*2 -0.5589063*2
qft_neg_2 = c(0, 0.9921135 + 0.0013404*2, -0.1351799 + 0.3145422*2) - 0.0093157*2
qft_pos_3 = c(0, 0.9921135 + 0.0013404*3, -0.1351799 + 0.3145422*3) + 2.6294794 - 0.0093157*3 -0.5589063*3
qft_neg_3 = c(0, 0.9921135 + 0.0013404*3, -0.1351799 + 0.3145422*3) - 0.0093157*3
time = rep(c("baseline", "peak", "memory"),6)

dfi = data.frame(value = c(qft_pos_1, qft_neg_1, qft_pos_2, qft_neg_2, qft_pos_3, qft_neg_3),
                 time = time,
                 cat = c(rep("QFT+ 1 admin.", 3), rep("QFT- 1 admin.", 3), rep("QFT+ 2 admin.", 3), rep("QFT- 2 admin.", 3), rep("QFT+ 3 admin.", 3), rep("QFT- 3 admin.", 3)))


dfi$time <- factor(dfi$time, levels=c("baseline", "peak", "memory"))
levels(dfi$time)
dfi$cat = as.factor(dfi$cat)

ggplot(dfi, aes(y=value, x=time, group = cat)) +
  geom_line(aes(color=cat),linewidth=2)+
  geom_point(size=5, aes(color = cat))+
  theme_minimal()+
  xlab(NULL)+
  ylab(expression(paste("log(", mu[ij], ")")))+
  theme(text = element_text(size = 28)) +
  scale_color_manual(values=c("#FF9999","#FF0000","#660000","#99CCFF","blue", "#000066"), name = NULL)


#-------------------------
# Model Checks:
#-------------------------

# Residuals vs Timepnt:
df2 <- data.frame(residuals = resid(mod4),  time = df1$timepnt)
ggplot(df2, aes(y=residuals, x=time)) +
  geom_boxplot(alpha=1, size =1, aes(fill = time), color ="black") +
  theme_minimal()+
  theme(text = element_text(size = 30), legend.position = "none") +
  ylab("Normalized quantile residuals") +
  xlab(NULL)+
  scale_fill_manual(values = brewer.pal(n = 4, name = "BuGn"))
# Residuals vs QFT:
df2 <- data.frame(residuals = resid(mod4),  QFT = df1$QFT)
ggplot(df2, aes(y=residuals, x=QFT)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color ="black") +
  theme_minimal()+
  theme(text = element_text(size = 30), legend.position = "none") +
  ylab("Normalized quantile residuals") +
  xlab(NULL)+
  scale_fill_manual(values=c("#FF0000", "blue"))+ 
  scale_x_discrete(labels=c("QFT neg" = "QFT-", "QFT pos" = "QFT+"))
# Residuals vs admin:
df2 <- data.frame(residuals = resid(mod4), Schedule = df1$Schedule)
df2$Schedule <- as.factor(df2$Schedule)
ggplot(df2, aes(y=residuals, x=Schedule)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Schedule), color ="black") +
  theme_minimal()+
  theme(text = element_text(size = 30), legend.position = "none") +
  ylab("Normalized quantile residuals") +
  xlab("Number of Administrations")+
  scale_fill_manual(values = brewer.pal(n = 3, name = "BuGn"))
# Residuals vs concentration:
df2 <- data.frame(residuals = resid(mod4),  Dose = df1$Dose)
df2$Dose <- as.factor(df2$Dose)
ggplot(df2, aes(y=residuals, x=Dose)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Dose), color ="black") +
  theme_minimal()+
  theme(text = element_text(size = 30), legend.position = "none") +
  ylab("Normalized quantile residuals") +
  xlab(expression(paste("Concentration (", mu, "g)")))+
  scale_fill_manual(values = brewer.pal(n = 3, name = "BuGn"))

# admin*timepnt:
df2 <- data.frame(residuals = resid(mod4),  time = df1$timepnt, Schedule=df1$Schedule)
df2$Schedule <- as.factor(df2$Schedule)
ggplot(df2, aes(y=residuals, x=time)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Schedule), color ="black") +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("Normalized quantile residuals") +
  xlab(NULL)+
  scale_fill_manual(values = brewer.pal(n = 3, name = "BuGn"), name="No. Admin.")

# admin*QFT
df2 <- data.frame(residuals = resid(mod4),  Schedule = df1$Schedule, QFT=df1$QFT)
df2$Schedule <- as.factor(df2$Schedule)
ggplot(df2, aes(y=residuals, x=Schedule)) +
  geom_boxplot(alpha=1, size =1, aes(fill = QFT), color ="black") +
  theme_minimal()+
  theme(text = element_text(size = 30)) +
  ylab("Normalized quantile residuals") +
  xlab("Number of Administrations")+
  scale_fill_manual(values=c("#FF0000", "blue"), labels = c("QFT-", "QFT+"), name="")



#-----------------------------------------
# Check Normality of Residuals 
#-----------------------------------------

df2 <- data.frame(residuals = resid(mod4))
ggplot(df2, aes(sample=residuals))+
  stat_qq()+
  theme_minimal(base_size = 27) + 
  geom_abline(linewidth=1, color="red")+ 
  theme(legend.position = "none")+
  labs(title = NULL, subtitle = NULL)+
  ylab("Sample Quantiles")+
  xlab("Normal Quantiles")




#---------------------------------------------
# Check Marginal Normality of Random Effects
#---------------------------------------------

df2 <- data.frame(randeff = ranef(getSmo(mod4)))
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
df2 <- data.frame(residuals = resid(mod4), fitted = fitted(mod4))
ggplot(data = df2, aes(x = fitted, y = residuals)) +
  geom_point(size=2) + 
  ylab("Normalized quantile residuals")+
  xlab("Fitted values")+
  theme_minimal(base_size = 22)



