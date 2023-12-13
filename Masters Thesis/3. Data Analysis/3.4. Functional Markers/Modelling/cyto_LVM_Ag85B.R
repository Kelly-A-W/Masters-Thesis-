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
library(cplm)
library(RColorBrewer)
library(boral)
library(fastDummies)
library(corrplot)

#---------------------
# Import Data
#---------------------

ag85b_df <- read.csv("Ag85b_compass_filtered_data_for_mods.csv", check.names = F)
str(ag85b_df)

#---------------------
# Prepare Model Data
#---------------------
# log transform data:
ag85b_df$value <- log(ag85b_df$value+0.001)

# convert to wide format:
df_long = ag85b_df[,c(4,5,6,8,10,14,15)]
df_wide = reshape(df_long, idvar = c("Schedule", "Dose", "timepnt", "QFT", "PID"), timevar = "Marker", direction = "wide")
str(df_wide)

# change cytokine combo names for later:
names(df_wide)[6:14] = c("IFNg+IL2+TNF+IL17+ ", 
                         "IFNg+IL2+TNF+IL17-", 
                         "IFNg+IL2+TNF-IL17-", 
                         "IFNg+IL2-TNF+IL17-", 
                         "IFNg+IL2-TNF-IL17-", 
                         "IFNg-IL2+TNF+IL17+", 
                         "IFNg-IL2+TNF+IL17-",
                         "IFNg-IL2+TNF-IL17-",
                         "IFNg-IL2-TNF+IL17-")

# create dummy variables:
d1 = dummy_cols(df_wide[,3])[,-c(1,2)]
d2 = dummy_cols(df_wide[,4])[,-c(1,2)]

# create covariate df X
X = cbind(df_wide[,c(1,2)], d1, d2)
names(X) = c("Schedule", "Dose", "memory", "peak", "QFTpos")

# create covariate matrix Xmat
Xmat <-model.matrix(~Schedule + Dose + memory + peak + QFTpos-1,
                          data = X)

# create response matrix Y
Y = df_wide[,6:14]

# Set Row IDs
row.ids <- matrix(c(as.integer(as.factor(df_wide$PID))),
                  nrow = length(c(as.integer(as.factor(df_wide$PID)))))

#---------------------
# Base Model
#---------------------

mod <- boral(Y,
             X = Xmat,
             family = "normal",
             lv.control = list(num.lv = 2),
             row.eff = "random",
             save.model = T,
             row.ids = row.ids)
summary(mod)

# 95% highest posterior density (credible) intervals
mod$hpdintervals

# perform residual analysis 
plot(mod,cex.lab=1.5, cex.axis=1.5)
grid()

# ordination plot
par(cex.lab=1.5, cex.axis=1.5)
lvsplot(mod)
grid()

# correlations plots
envcors <- get.enviro.cor(mod)
colnames(envcors$sig.cor) = c("IFNg+IL2+TNF+IL17+ ", 
                              "IFNg+IL2+TNF+IL17-", 
                              "IFNg+IL2+TNF-IL17-", 
                              "IFNg+IL2-TNF+IL17-", 
                              "IFNg+IL2-TNF-IL17-", 
                              "IFNg-IL2+TNF+IL17+", 
                              "IFNg-IL2+TNF+IL17-",
                              "IFNg-IL2+TNF-IL17-",
                              "IFNg-IL2-TNF+IL17-")
rownames(envcors$sig.cor) = colnames(envcors$sig.cor)
rescors <- get.residual.cor(mod)
colnames(rescors$sig.cor) = colnames(envcors$sig.cor)
rownames(rescors$sig.cor) = colnames(envcors$sig.cor)                                                       
# Correlations due to covariates: (Only the significant correlations, as based on the 95% credible intervals excluding zero, have been plotted)
corrplot(envcors$sig.cor, type = "lower", diag = FALSE, tl.srt = 45, tl.cex = 1.2,mar = c(4, 0, 0, 0))
# Residual correlations: (Only the significant correlations, as based on the 95% credible intervals excluding zero, have been plotted)
corrplot(rescors$sig.cor, type = "lower", diag = FALSE, tl.srt = 45, tl.cex = 1.2,mar = c(4, 0, 0, 0))

# quantify how much of the cytokine combination co-occurrence is explained by covariates:
rescors$trace


#--------------------------
# Model with interactions
#--------------------------

Xmat_int <-model.matrix(~Schedule*memory + Schedule*peak +
                          Dose*memory + Dose*peak +
                          QFTpos*memory + QFTpos*peak - 1,
                    data = X)
mod2 <- boral(Y,
             X = Xmat_int,
             family = "normal",
             lv.control = list(num.lv = 2),
             row.eff = "random",
             save.model = T,
             row.ids = row.ids)

summary(mod2)
# order nicely for tables:
summary(mod2)$X.coefficients[,c(3,2,5,1,4,7,6)]
cbind(summary(mod2)$X.coefficients[,c(9,8,11,10)],summary(mod2)$coefficients)

# 95% highest posterior density (credible) intervals
mod2$hpdintervals
round(mod2$hpdintervals$X.coefs[,c(3,2,5,1,4,7,6),1],3) # lower bound
round(mod2$hpdintervals$X.coefs[,c(3,2,5,1,4,7,6),2],3) # upper bound
round(cbind(mod2$hpdintervals$X.coefs[,c(9,8,11,10),1], mod2$hpdintervals$lv.coefs[,,1]),3) # lower bound
round(cbind(mod2$hpdintervals$X.coefs[,c(9,8,11,10),2], mod2$hpdintervals$lv.coefs[,,2]),3) # upper bound

# perform residual analysis 
plot(mod2,cex.lab=1.5, cex.axis=1.5)
grid()

# ordination plot
par(cex.lab=1.5, cex.axis=1.5)
lvsplot(mod2, cex=1.2)
grid()

#coef2 = data.frame(x=summary(mod2)$coefficients[,2], 
                   #y=summary(mod2)$coefficients[,3], 
                   #combo = c("IFNg+IL2+TNF+IL17+ ", 
                             #"IFNg+IL2+TNF+IL17-", 
                             #"IFNg+IL2+TNF-IL17-", 
                             #"IFNg+IL2-TNF+IL17-", 
                             #"IFNg+IL2-TNF-IL17-", 
                             #"IFNg-IL2+TNF+IL17+", 
                             #"IFNg-IL2+TNF+IL17-",
                             #"IFNg-IL2+TNF-IL17-",
                             #"IFNg-IL2-TNF+IL17-"))
#scores2 = as.data.frame(summary(mod2)$lvs)  
#plot <- ggplot(scores2, aes(x=lv1, y=lv2))+geom_point()+xlim(-4,4)+ylim(-4,4)
#plot <- plot + coord_equal() + geom_text(data=coef2, aes(x=x, y=y, label=combo), size = 5, vjust=1, color="red")
#plot <- plot + geom_segment(data=coef2, aes(x=0, y=0, xend=x, yend=y), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75, color="red")
#plot




# correlations plots
envcors <- get.enviro.cor(mod2)
colnames(envcors$sig.cor) = c("IFNg+IL2+TNF+IL17+ ", 
                              "IFNg+IL2+TNF+IL17-", 
                              "IFNg+IL2+TNF-IL17-", 
                              "IFNg+IL2-TNF+IL17-", 
                              "IFNg+IL2-TNF-IL17-", 
                              "IFNg-IL2+TNF+IL17+", 
                              "IFNg-IL2+TNF+IL17-",
                              "IFNg-IL2+TNF-IL17-",
                              "IFNg-IL2-TNF+IL17-")
rownames(envcors$sig.cor) = colnames(envcors$sig.cor)
rescors <- get.residual.cor(mod2)
colnames(rescors$sig.cor) = colnames(envcors$sig.cor)
rownames(rescors$sig.cor) = colnames(envcors$sig.cor)                                                       
# Correlations due to covariates: (Only the significant correlations, as based on the 95% credible intervals excluding zero, have been plotted)
corrplot(envcors$sig.cor, type = "lower", diag = FALSE, tl.srt = 45, tl.cex = 1.2,mar = c(4, 0, 0, 0))
# Residual correlations: (Only the significant correlations, as based on the 95% credible intervals excluding zero, have been plotted)
corrplot(rescors$sig.cor, type = "lower", diag = FALSE, tl.srt = 45, tl.cex = 1.2,mar = c(4, 0, 0, 0))

# quantify how much of the cytokine combination co-occurrence is explained by covariates:
rescors$trace


#--------------------------
# Model with traits
#--------------------------

traits <- matrix(
  c(1,1,1,1,  # IFNg IL2 IL17 TNF
    1,1,0,1,
    1,1,0,0,
    1,0,0,1,
    1,0,0,0,
    0,1,1,1,
    0,1,0,1,
    0,1,0,0,
    0,0,0,1),
  nrow = 9, 
  byrow = T)
which_traits_1 <- vector("list",ncol(Xmat_int)+1)
for(i in 1:length(which_traits_1))
  which_traits_1[[i]] <- 1:ncol(traits)



mod3 <- boral(Y,
              X = Xmat_int,
              family = "normal",
              traits = traits,
              which.traits = which_traits_1,
              lv.control = list(num.lv = 2),
              row.eff = "random",
              save.model = T,
              row.ids = row.ids)
  
summary(mod3)
# order nicely for tables:
summary(mod3)$traits.coefficients[c(1,4,3,6,2,5,8,7,10,9,12,11), c(1,2,3,5,4,6)] # original kappa order: IFNg IL2 IL17 TNF; swap to IFNg IL2 TNF IL17
summary(mod3)$coefficients[,-1] #ignore beta0

# 95% highest posterior density (credible) intervals
mod3$hpdintervals
round(mod3$hpdintervals$traits.coefs[c(1,4,3,6,2,5,8,7,10,9,12,11),c(1,2,3,5,4,6),1],3) # lower bound
round(mod3$hpdintervals$traits.coefs[c(1,4,3,6,2,5,8,7,10,9,12,11),c(1,2,3,5,4,6),2],3) # upper bound
round(mod3$hpdintervals$lv.coefs[,-1,1],3) # lower bound
round(mod3$hpdintervals$lv.coefs[,-1,2],3) # upper bound

# perform residual analysis 
plot(mod3,cex.lab=1.5, cex.axis=1.5)
grid()

# ordination plot
par(cex.lab=1.5, cex.axis=1.5)
lvsplot(mod3, cex=1.2)
grid()

# correlations plots
envcors <- get.enviro.cor(mod3)
colnames(envcors$sig.cor) = c("IFNg+IL2+TNF+IL17+ ", 
                              "IFNg+IL2+TNF+IL17-", 
                              "IFNg+IL2+TNF-IL17-", 
                              "IFNg+IL2-TNF+IL17-", 
                              "IFNg+IL2-TNF-IL17-", 
                              "IFNg-IL2+TNF+IL17+", 
                              "IFNg-IL2+TNF+IL17-",
                              "IFNg-IL2+TNF-IL17-",
                              "IFNg-IL2-TNF+IL17-")
rownames(envcors$sig.cor) = colnames(envcors$sig.cor)
rescors <- get.residual.cor(mod3)
colnames(rescors$sig.cor) = colnames(envcors$sig.cor)
rownames(rescors$sig.cor) = colnames(envcors$sig.cor)                                                       
# Correlations due to covariates: (Only the significant correlations, as based on the 95% credible intervals excluding zero, have been plotted)
corrplot(envcors$sig.cor, type = "lower", diag = FALSE, tl.srt = 45, tl.cex = 1.2,mar = c(4, 0, 0, 0))
# Residual correlations: (Only the significant correlations, as based on the 95% credible intervals excluding zero, have been plotted)
corrplot(rescors$sig.cor, type = "lower", diag = FALSE, tl.srt = 45, tl.cex = 1.2,mar = c(4, 0, 0, 0))

# quantify how much of the cytokine combination co-occurrence is explained by covariates:
rescors$trace
