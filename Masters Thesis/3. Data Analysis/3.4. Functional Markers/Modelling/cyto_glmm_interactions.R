# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

#---------------------
# Load Libraries
#---------------------
library(ggplot2)

#---------------------
# ESAT-6
#---------------------
a1c5= c(0, 1.945 - 0.378*1 - 0.035*5 + 0.014*1*5, 0.630 + 0.031*1 - 0.020*5 + 0.007*1*5)-0.056*1 + 0.002*5 - 0.001*1*5
a1c15 = c(0, 1.945 - 0.378*1 - 0.035*15 + 0.014*1*15, 0.630 + 0.031*1 - 0.020*15 + 0.007*1*15)-0.056*1 + 0.002*15 - 0.001*1*15
a1c50 = c(0, 1.945 - 0.378*1 - 0.035*50 + 0.014*1*50, 0.630 + 0.031*1 - 0.020*50 + 0.007*1*50)-0.056*1 + 0.002*50 - 0.001*1*50

a2c5= c(0, 1.945 - 0.378*2 - 0.035*5 + 0.014*2*5, 0.630 + 0.031*2 - 0.020*5 + 0.007*2*5)-0.056*2 + 0.002*5 - 0.001*2*5
a2c15 = c(0, 1.945 - 0.378*2 - 0.035*15 + 0.014*2*15, 0.630 + 0.031*2 - 0.020*15 + 0.007*2*15)-0.056*2 + 0.002*15 - 0.001*2*15
a2c50 = c(0, 1.945 - 0.378*2 - 0.035*50 + 0.014*2*50, 0.630 + 0.031*2 - 0.020*50 + 0.007*2*50)-0.056*2 + 0.002*50 - 0.001*2*50

a3c5= c(0, 1.945 - 0.378*3 - 0.035*5 + 0.014*3*5, 0.630 + 0.031*3 - 0.020*5 + 0.007*3*5)-0.056*3 + 0.002*5 - 0.001*3*5
a3c15 = c(0, 1.945 - 0.378*3 - 0.035*15 + 0.014*3*15, 0.630 + 0.031*3 - 0.020*15 + 0.007*3*15)-0.056*3 + 0.002*15 - 0.001*3*15
a3c50 = c(0, 1.945 - 0.378*3 - 0.035*50 + 0.014*3*50, 0.630 + 0.031*3 - 0.020*50 + 0.007*3*50)-0.056*3 + 0.002*50 - 0.001*3*50



time = rep(c("baseline", "peak", "memory"),9)

df = data.frame(value = c(a1c5, a1c15, a1c50, a2c5, a2c15, a2c50, a3c5, a3c15, a3c50),
                time = time,
                cat = c(rep("1 admin, 5 conc", 3), rep("1 admin, 15 conc", 3), rep("1 admin, 50 conc", 3), 
                        rep("2 admin, 5 conc", 3), rep("2 admin, 15 conc", 3), rep("2 admin, 50 conc", 3),
                        rep("3 admin, 5 conc", 3), rep("3 admin, 15 conc", 3), rep("3 admin, 50 conc", 3)))


df$time <- factor(df$time, levels=c("baseline", "peak", "memory"))
levels(df$time)
df$cat = factor(df$cat, levels=c("1 admin, 5 conc", "1 admin, 15 conc", "1 admin, 50 conc", 
                                 "2 admin, 5 conc", "2 admin, 15 conc", "2 admin, 50 conc",
                                 "3 admin, 5 conc", "3 admin, 15 conc", "3 admin, 50 conc"))

ggplot(df, aes(y=value, x=time, group = cat)) +
  geom_line(aes(color=cat),linewidth=2)+
  geom_point(size=5, aes(color = cat))+
  theme_minimal()+
  xlab(NULL)+
  ylab(expression(paste("log(", mu[ij], ")")))+
  theme(text = element_text(size = 28)) +
  scale_color_manual(values=c("#ADDD8E","#41AB5D","#006837",
                              "#9ECAE1", "#4292C6", "#08519C",
                              "#FA9FB5", "#DD3497", "#7A0177"), name = NULL)

#---------------------
# Ag85B
#---------------------
a1c5= c(0, 0.508 + 0.198*1 -0.005*5, 0.018 + 0.224*1 + 0.0004*5) + 0.069*1 -0.003*5
a1c15 = c(0, 0.508 + 0.198*1 -0.005*15, 0.018 + 0.224*1 + 0.0004*15) + 0.069*1 -0.003*15
a1c50 = c(0, 0.508 + 0.198*1 -0.005*50, 0.018 + 0.224*1 + 0.0004*50) + 0.069*1 -0.003*50

a2c5= c(0, 0.508 + 0.198*2 -0.005*5, 0.018 + 0.224*2 + 0.0004*5) + 0.069*2 -0.003*5
a2c15 = c(0, 0.508 + 0.198*2 -0.005*15, 0.018 + 0.224*2 + 0.0004*15) + 0.069*2 -0.003*15
a2c50 = c(0, 0.508 + 0.198*2 -0.005*50, 0.018 + 0.224*2 + 0.0004*50) + 0.069*2 -0.003*50

a3c5= c(0, 0.508 + 0.198*3 -0.005*5, 0.018 + 0.224*3 + 0.0004*5) + 0.069*3 -0.003*5
a3c15 = c(0, 0.508 + 0.198*3 -0.005*15, 0.018 + 0.224*3 + 0.0004*15) + 0.069*3 -0.003*15
a3c50 = c(0, 0.508 + 0.198*3 -0.005*50, 0.018 + 0.224*3 + 0.0004*50) + 0.069*3 -0.003*50



time = rep(c("baseline", "peak", "memory"),9)

df = data.frame(value = c(a1c5, a1c15, a1c50, a2c5, a2c15, a2c50, a3c5, a3c15, a3c50),
                time = time,
                cat = c(rep("1 admin, 5 conc", 3), rep("1 admin, 15 conc", 3), rep("1 admin, 50 conc", 3), 
                        rep("2 admin, 5 conc", 3), rep("2 admin, 15 conc", 3), rep("2 admin, 50 conc", 3),
                        rep("3 admin, 5 conc", 3), rep("3 admin, 15 conc", 3), rep("3 admin, 50 conc", 3)))


df$time <- factor(df$time, levels=c("baseline", "peak", "memory"))
levels(df$time)
df$cat = factor(df$cat, levels=c("1 admin, 5 conc", "1 admin, 15 conc", "1 admin, 50 conc", 
                                 "2 admin, 5 conc", "2 admin, 15 conc", "2 admin, 50 conc",
                                 "3 admin, 5 conc", "3 admin, 15 conc", "3 admin, 50 conc"))

ggplot(df, aes(y=value, x=time, group = cat)) +
  geom_line(aes(color=cat),linewidth=2)+
  geom_point(size=5, aes(color = cat))+
  theme_minimal()+
  xlab(NULL)+
  ylab(expression(paste("log(", mu[ij], ")")))+
  theme(text = element_text(size = 28)) +
  scale_color_manual(values=c("#ADDD8E","#41AB5D","#006837",
                              "#9ECAE1", "#4292C6", "#08519C",
                              "#FA9FB5", "#DD3497", "#7A0177"), name = NULL)

