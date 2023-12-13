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
library(RColorBrewer)


#---------------------
# Import Data
#---------------------

ag85b_df <- read.csv("Ag85b_compass_filtered_data_for_mods.csv", check.names = F)
str(ag85b_df)

#---------------------
# Set Factors
#---------------------
ag85b_df$ethnicity = as.factor(ag85b_df$ethnicity)
ag85b_df$sex = as.factor(ag85b_df$sex)
ag85b_df$timepnt = as.factor(ag85b_df$timepnt)
ag85b_df$QFT = as.factor(ag85b_df$QFT)
ag85b_df$PID = as.factor(ag85b_df$PID)
ag85b_df$Marker = as.factor(ag85b_df$Marker)

ag85b_df$timepnt <- factor(ag85b_df$timepnt, levels=c("baseline", "peak", "memory"))
levels(ag85b_df$timepnt)

#---------------------------------------------------------
# Choose Random Effect
#---------------------------------------------------------


#--------------------------------------------------- RE on Intercept Only

# Convert data:
df <- ag85b_df
df$value <- log(df$value+0.001)

# Keep Only complete cases:
df1<-df[complete.cases(df),]


# Intercept Model:
int_model <-gamlss(value ~ 
                     QFT + 
                     Schedule + 
                     Dose +
                     timepnt +
                     Marker +
                     re(random=~1|PID, method = "REML"),  
                   data=df1, 
                   family=NO(mu.link=identity, sigma.link = identity)) 
summary(int_model)
getSmo(int_model)  # NB getSmo convert my model to a lme object 
intervals(getSmo(int_model) )



#--------------------------------------------------- RE on Intercept and Slope

mod_slo <- gamlss(value ~ 
                    QFT + 
                    Schedule + 
                    Dose +
                    timepnt +
                    Marker +
                    re(random=~timepnt|PID, method = "REML"),  
                  data=df1, 
                  family=NO(mu.link=identity, sigma.link = identity))


# re on marker did not converge



# GAIC:
GAIC(int_model, k=round(log(length(df1$value)), 2))    # choose int_model
GAIC(mod_slo, k=round(log(length(df1$value)), 2)) 

BIC(int_model)
BIC(mod_slo) 


#---------------------------------------------------------
# Choose Fixed Effects: Step-wise Selection
#---------------------------------------------------------

# Step-wise Selection:
mod1 <- stepGAIC(int_model, 
                 scope = list(lower= ~ QFT + Schedule + Dose + timepnt + Marker + re(random=~1|PID, method = "REML"),
                              upper = ~ QFT + Schedule + Dose + timepnt + Marker + age + ethnicity + sex + Study + Day + re(random=~1|PID, method = "REML")),
                 k=round(log(length(df1$value)), 2))
summary(mod1)




# Since Sex was not included in our model, refit model using observation with missing sex data included 
# (note that the model will not be comparebale in terms of AIC and BIC to the previous model)
df <- ag85b_df
df$value <- log(df$value+0.001)

# Keep Only complete cases:
df <- df[,-c(1,2,3)]
df1<-df[complete.cases(df),]


# Base Model:
mod1 <-gamlss(value ~ 
                QFT + 
                Schedule + 
                Dose +
                timepnt +
                Marker +
                re(random=~1|PID, method = "REML"),  
              data=df1, 
              family=NO(mu.link=identity, sigma.link = identity)) # NB idenity link because lmm



#---------------------------------------------------------
# Choose Interactions: Step-wise Selection
#---------------------------------------------------------

# Step-wise Selection:
mod2 <- stepGAIC(mod1, 
                 scope = list(lower= ~ QFT + Schedule + Dose + timepnt + Marker + re(random=~1|PID, method = "REML"),
                              upper = ~ Schedule*Dose + timepnt*Schedule + QFT*Dose + QFT*Schedule + QFT*timepnt +
                                Dose*timepnt + Marker*QFT + Marker*Schedule + Marker*Dose + Marker*timepnt + re(random=~1|PID, method = "REML")),
                 k=round(log(length(df1$value)), 2))
summary(mod2) # choose mod2


#----------------------------
# 3-way Interactions
#----------------------------

mod3.1 <-update(mod2, ~.+QFT*Marker*timepnt)
mod3.2 <-update(mod2, ~.+Dose*Marker*timepnt)
mod3.3 <-update(mod2, ~.+Schedule*Marker*timepnt)

GAIC(mod2, k=round(log(length(df1$value)), 2))  # choose mod2
GAIC(mod3.1, k=round(log(length(df1$value)), 2)) 
GAIC(mod3.2, k=round(log(length(df1$value)), 2)) 
GAIC(mod3.3, k=round(log(length(df1$value)), 2)) 

# don't add any 3-way interactions

#---------------------------------------------------------
# Choose Sigma Model: Step-wise Selection
#---------------------------------------------------------
# Step-wise Selection for Sigma: 
mod4.0 <-gamlss(value ~ 
                  Schedule*Marker + 
                  Dose*Marker +
                  Dose*timepnt+
                  timepnt*Schedule +
                  Marker*timepnt + 
                  QFT*Marker +
                  QFT*timepnt+
                  re(random=~1|PID, method = "REML"),  
                data=df1, 
                sigma.formula = ~ 1, 
                family=NO(mu.link=identity, sigma.link = log)) # log link for sigma because adding formula

mod4 <- stepGAIC(mod4.0, 
                 scope = list(lower= ~ 1,
                              upper = ~ QFT + Marker + Schedule + Dose+timepnt), 
                 k=round(log(length(df1$value)), 2), 
                 what = "sigma")
summary(mod4) # Marker, QFT, dose and schedule selected



# FINAL MODEL:
summary(mod4)
intervals(getSmo(mod4))
GAIC(mod4, k=round(log(length(df1$value)), 2))


#-----------------------------------
# Residuals by category:
#-----------------------------------
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
# Residuals vs Marker:
df2 <- data.frame(residuals = resid(mod4),  Marker = df1$Marker)
df2$Marker <- as.factor(df2$Marker)
p=ggplot(df2, aes(y=residuals, x=Marker)) +
  geom_boxplot(alpha=1, size =1, aes(fill = Marker), color ="black") +
  theme_minimal()+
  theme(text = element_text(size = 30), legend.position = "none") +
  ylab("Normalized quantile residuals") +
  xlab(NULL)+
  scale_x_discrete(labels = rep("", 9))+
  scale_fill_manual(values = brewer.pal(n = 9, name = "BuGn"))
# Make grobtable:
grob         <- data.frame(Marker = c("IL17", "TNF", "IL2", "G"), 
                           x1 = c("+","+",  "+", "+"),
                           x2 = c("-","+",  "+", "+"),
                           x3 = c("-","-",  "+", "+"),
                           x4 = c("-","+",  "-", "+"),
                           x5 = c("-","-",  "-", "+"),
                           x6 = c("+","+",  "+", "-"),
                           x7 = c("-","+",  "+", "-"),
                           x8 = c("-","-",  "+", "-"),
                           x9 = c("-","+",  "-", "-"))
tab         <- tableGrob(grob[,-1], rows=NULL, theme = ttheme_default(base_size = 17), cols=NULL)
# Set widths/heights to 'fill whatever space I have':
tab$widths  <- unit(rep(1, ncol(tab)), "null")
tab$heights <- unit(rep(1, nrow(tab)), "null")
# Format table as plot:
p3 <- ggplot() +
  annotation_custom(tab)+
  scale_y_discrete(breaks = rev(grob$Marker),limits = rev(grob$Marker), 
                   labels = c(expression(paste("IFN", gamma)), "IL2", "TNF", "IL17"))+
  theme(axis.text.y = element_text(face="bold", size=15))
# Patchwork magic
p + plot_spacer() + p3 + plot_layout(ncol = 1, heights = c(12, -2, 3))






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
# Check Marginal Normality of Random Errors
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

