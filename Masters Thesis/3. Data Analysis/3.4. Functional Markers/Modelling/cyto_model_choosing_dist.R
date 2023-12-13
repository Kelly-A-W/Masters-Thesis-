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

cd4_counts <- read.csv("clean_data_with_demographic_info.csv", check.names = F)
str(cd4_counts)


pp_esat6 = read.csv("compass_esat6_dataset_1.csv", check.names = F)
pp_ag85b = read.csv("compass_ag85b_dataset_1.csv", check.names = F)
str(pp_esat6)
str(pp_ag85b)

#---------------------
# BS Frequency of CD4
#---------------------
# Only keep counts for cells expressing at least one of INFg, IL2, TNF or L17 (i.e. cytokine+)
cd4_counts <- cd4_counts[,-which(grepl("Gc..L2c..L17c..TNFc", names(cd4_counts), fixed = TRUE)==T)]
str(cd4_counts)


# Sum over memory markers:
cyto_counts = data.frame(cd4_counts[,1:13],
                         cd4_counts[,14:28] + cd4_counts[,29:43] + cd4_counts[,44:58] + cd4_counts[,59:73])
names(cyto_counts)[14:28] = gsub("R7..RA..", "", names(cyto_counts)[14:28])


# Convert to Frequencies (%):
cyto_counts[14:28] <- (cyto_counts[14:28]/cyto_counts$CD4count)*100


# Background Subtract
cd4_bs <- data.frame()
for(i in 1:length(unique(cyto_counts$PID))){
  pid <- unique(cyto_counts$PID)[i]
  dat <- cyto_counts[which(cyto_counts$PID == pid),]
  days <- unique(dat$Day)
  for (j in 1:length(days)) {
    day <- days[j]
    df  <- dat[which(dat$Day == day),]
    df[-which(df$Stim == "UNS"),14:28]  <- df[-which(df$Stim == "UNS"),14:28] - df[rep(which(df$Stim == "UNS"), (length(df$Stim)-1)),14:28]
    cd4_bs <- rbind(cd4_bs,df)
  }
}
summary(cd4_bs)
which(is.na(cd4_bs)==T)
str(cd4_bs)

# Deal with Negative Frequencies:
for (i in 14:ncol(cd4_bs)) {
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
sapply(sign(cd4_bs[,-(1:14)]), table)


# Only keep Ag85B and ESAT6
cd4_bs <- cd4_bs[cd4_bs[,"Stim"] == "ESAT6" | cd4_bs[,"Stim"] == "Ag85B",]


#----------------------------
# Deal with Missing Values:
#----------------------------

# Number of missing observations:
which(table(cd4_bs$PID)<6) # 5 pids with missing visits: 3020, 3026 and 411
cd4_bs[cd4_bs[,"PID"]=="3001 H56-035",] # Ag85B missing baseline
cd4_bs[cd4_bs[,"PID"]==3020,] # Ag85B and ESAT6 missing memory
cd4_bs[cd4_bs[,"PID"]==3024,]  # Ag85B missing baseline
cd4_bs[cd4_bs[,"PID"]==3026,] # ESAT6 and AG85B missing baseline
cd4_bs[cd4_bs[,"PID"]==411,]  # ESAT6 and Ag85B missing peak and memory - think  we should remove

# Remove 411:
cd4_bs <- cd4_bs[-which(cd4_bs$PID == 411),]



#---------------------
# Filter Data
#---------------------

# create dataframes to filter
cd4_bs = data.frame(cd4_bs)
esat6 <- cd4_bs[cd4_bs[,"Stim"] == "ESAT6",]
ag85b <- cd4_bs[cd4_bs[,"Stim"] == "Ag85B",]

# identify columns to remove
(x1 = which(colSums(pp_esat6[,8:22]>0.1)/nrow(pp_esat6)<0.3)) # ignore last column because that is Gc..L2c..L17c..TNFc
(x2 = which(colSums(pp_ag85b[,8:22]>0.1)/nrow(pp_ag85b)<0.3))

# check that names are in the same order:
sum(names(pp_esat6[,8:22]) != names(esat6[,14:28])) # 0, so same order
sum(names(pp_ag85b[,8:22]) != names(ag85b[,14:28])) # 0, so same order

# create filtered datasets:
esat6_fil = esat6[,-(13+as.numeric(x1))] # removing all marker combos identified as non-repsonsive
ag85b_fil = ag85b[,-(13+as.numeric(x2))] 


#---------------------
# Prepare Data
#---------------------

# Two-Level Grouping:
long_esat6 <- melt(setDT(esat6_fil), id.vars = c("age", "ethnicity", "sex", "Schedule","Dose", "timepnt", "Study", "QFT", "Group", "PID", "Day", "Stim", "CD4count"),
             variable.name = "Marker")
long_esat6$Marker <- as.factor(long_esat6$Marker)
long_esat6 = as.data.frame(long_esat6)

long_ag85b <- melt(setDT(ag85b_fil), id.vars = c("age", "ethnicity", "sex","Schedule","Dose", "timepnt", "Study", "QFT", "Group", "PID", "Day", "Stim", "CD4count"),
                   variable.name = "Marker")
long_ag85b$Marker <- as.factor(long_ag85b$Marker)
long_ag85b = as.data.frame(long_ag85b)




#-------------------------
# Export Filtered Data
#-------------------------


write.csv(long_esat6,"Esat6_compass_filtered_data_for_mods.csv", row.names = FALSE) 
write.csv(long_ag85b,"Ag85b_compass_filtered_data_for_mods.csv", row.names = FALSE) 



#---------------------------------------------------------------------
# Univariate EDA: Choosing an Appropriate Distribution for Y:
#---------------------------------------------------------------------



#---------------------
# Density Plots:
#---------------------


#----------------------------------------------------- Two-Level Grouping:
ggplot(long_esat6, aes(x = value)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("Frequency")

ggplot(long_ag85b, aes(x = value)) + 
  theme_minimal(base_size = 27)+
  geom_histogram(aes(y =after_stat(density)), colour = 1, fill = "white", bins=40, boundary = 0) +
  geom_density(lwd = 1, colour = 4,
               fill = 4, alpha = 0.25)+
  ylab("Density") +
  xlab("Frequency")


#-------------------------
# Fitting a Distribution
#-------------------------

#-------------------------------------------------------------------------------
#----------------------------------------------------------- Two-Level Grouping: 
#-------------------------------------------------------------------------------

#---------------------------------------------------------------------- ESAT6
f1.1 <- fitdist((long_esat6$value+0.001), "lnorm")
f1.3 <- fitdist((long_esat6$value+0.001), "gamma")
f1.4 <- fitdist((long_esat6$value+0.001), "weibull")
model_1 <- cpglm(
  value ~1,
  data = long_esat6,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)



# QQ PLOTS:
ggplot(long_esat6, aes(sample=value)) +
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
f1.1 <- fitdist((long_ag85b$value+0.001), "lnorm")
f1.3 <- fitdist((long_ag85b$value+0.001), "gamma")
f1.4 <- fitdist((long_ag85b$value+0.001), "weibull")
model_1 <- cpglm(
  value ~1,
  data = long_ag85b,
  control = list(max.iter = 5000, trace = 1),
  optimizer = "nlminb",
  doFit = TRUE,
  nAGQ = 1L)
summary(model_1)
summary(f1.1)
summary(f1.3)
summary(f1.4)


# QQ PLOTS:
ggplot(long_ag85b, aes(sample=value)) +
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




