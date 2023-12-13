# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()


#------------------------------------
# Load Libraries
#------------------------------------

library("readxl")
library(DescTools)

#------------------------------------
# Load Data
#------------------------------------


# PID 063
my_cd4_counts_063 <- read.csv("cd4_counts_pid_063.csv", check.names = F) #ignore columns 1-16
my_cd8_counts_063 <- read.csv("cd8_counts_pid_063.csv", check.names = F)
true_cd4_counts_063 <- read.csv("Flowjo_CD4_Counts_063.csv", check.names = F)  #ignore last column
true_cd4_counts_063 <- true_cd4_counts_063[,-66]
true_cd8_counts_063 <- read.csv("Flowjo_CD8_Counts_063.csv", check.names = F)  
true_cd8_counts_063 <- true_cd8_counts_063[,-66]


# PID 065
my_cd4_counts_065 <- read.csv("cd4_counts_pid_065.csv", check.names = F) #ignore columns 1-16
my_cd8_counts_065 <- read.csv("cd8_counts_pid_065.csv", check.names = F)
true_cd4_counts_065 <- read_excel("Flowjo_CD4_counts_065.xls")  #ignore last column
true_cd4_counts_065 <- as.data.frame(true_cd4_counts_065)
true_cd8_counts_065 <- read.csv("Flowjo_CD8_counts_065.csv", check.names = F)  
true_cd8_counts_065 <- true_cd8_counts_065[,-66]


# PID 131
my_cd4_freq_131 <- read.csv("cd4_frequency_pid_131.csv", check.names = F) #ignore columns 1-16
true_cd4_freq_131 <- read.csv("CD4 Boolean 131.csv", check.names = F)  #ignore last column
true_cd4_freq_131 <- true_cd4_freq_131[,-66]


# PID 185
my_cd4_freq_185 <- read.csv("cd4_frequency_pid_185.csv", check.names = F) #ignore columns 1-16
true_cd4_freq_185 <- read.csv("CD4 Boolean 185.csv", check.names = F)  #ignore last column
true_cd4_freq_185 <- true_cd4_freq_185[,-66]




#------------------------------------
# Match Columns
#------------------------------------

# PID 063 CD4
names(my_cd4_counts_063)
names(true_cd4_counts_063)
names(true_cd4_counts_063) <- gsub("CCR7, Time+ & CD3, Time+ & CD8, Time+ & IFN-g, Time+/Singlets/Lymphocytes/keeper gate/CD3+/CD8-/CD4+/",
                                "", names(true_cd4_counts_063), fixed = TRUE)
names(true_cd4_counts_063) <- gsub(" | Count",
                                "", names(true_cd4_counts_063), fixed = TRUE)
names(true_cd4_counts_063) <- gsub("-",
                                "c ", names(true_cd4_counts_063), fixed = TRUE)
names(true_cd4_counts_063) <- gsub("+",
                                " ", names(true_cd4_counts_063), fixed = TRUE)
names(true_cd4_counts_063) <- gsub("2",
                                "L2", names(true_cd4_counts_063), fixed = TRUE)
names(true_cd4_counts_063) <- gsub("17",
                                "L17", names(true_cd4_counts_063), fixed = TRUE)
names(true_cd4_counts_063) <- gsub("T ",
                                "TNF", names(true_cd4_counts_063), fixed = TRUE)
names(true_cd4_counts_063) <- gsub("Tc ",
                                "TNFc", names(true_cd4_counts_063), fixed = TRUE)
names(true_cd4_counts_063)[1] <- "Day and Stim"
names(true_cd4_counts_063)
names(my_cd4_counts_063) <- gsub(",",
                              "", names(my_cd4_counts_063), fixed = TRUE)
names(my_cd4_counts_063)
which(names(true_cd4_counts_063[,-1]) != names(my_cd4_counts_063[,-c(1,2,3,4)]))


# PID 063 CD8
names(my_cd8_counts_063)
names(true_cd8_counts_063)
names(true_cd8_counts_063) <- gsub("CCR7, Time+ & CD3, Time+ & CD8, Time+ & IFN-g, Time+/Singlets/Lymphocytes/keeper gate/CD3+/CD4-/CD8+/",
                                "", names(true_cd8_counts_063), fixed = TRUE)
names(true_cd8_counts_063) <- gsub(" | Count",
                                "", names(true_cd8_counts_063), fixed = TRUE)
names(true_cd8_counts_063) <- gsub("-",
                                "c ", names(true_cd8_counts_063), fixed = TRUE)
names(true_cd8_counts_063) <- gsub("+",
                                " ", names(true_cd8_counts_063), fixed = TRUE)
names(true_cd8_counts_063) <- gsub("2",
                                "L2", names(true_cd8_counts_063), fixed = TRUE)
names(true_cd8_counts_063) <- gsub("17",
                                "L17", names(true_cd8_counts_063), fixed = TRUE)
names(true_cd8_counts_063) <- gsub("T ",
                                "TNF", names(true_cd8_counts_063), fixed = TRUE)
names(true_cd8_counts_063) <- gsub("Tc ",
                                "TNFc", names(true_cd8_counts_063), fixed = TRUE)
names(true_cd8_counts_063)[1] <- "Day and Stim"
names(true_cd8_counts_063)
names(my_cd8_counts_063) <- gsub(",",
                              "", names(my_cd8_counts_063), fixed = TRUE)
names(my_cd8_counts_063)
which(names(true_cd8_counts_063[,-1]) != names(my_cd8_counts_063[,-c(1,2,3,4)]))


# PID 065 CD4
names(my_cd4_counts_065)
names(true_cd4_counts_065)
names(true_cd4_counts_065) <- gsub("CCR7, Time+ & CD3, Time+ & CD8, Time+ & IFN-g, Time+/Singlets/Lymphocytes/keeper gate/CD3+/CD8-/CD4+/",
                                "", names(true_cd4_counts_065), fixed = TRUE)
names(true_cd4_counts_065) <- gsub(" | Count",
                                "", names(true_cd4_counts_065), fixed = TRUE)
names(true_cd4_counts_065) <- gsub("-",
                                "c ", names(true_cd4_counts_065), fixed = TRUE)
names(true_cd4_counts_065) <- gsub("+",
                                " ", names(true_cd4_counts_065), fixed = TRUE)
names(true_cd4_counts_065) <- gsub("2",
                                "L2", names(true_cd4_counts_065), fixed = TRUE)
names(true_cd4_counts_065) <- gsub("17",
                                "L17", names(true_cd4_counts_065), fixed = TRUE)
names(true_cd4_counts_065) <- gsub("T ",
                                "TNF", names(true_cd4_counts_065), fixed = TRUE)
names(true_cd4_counts_065) <- gsub("Tc ",
                                "TNFc", names(true_cd4_counts_065), fixed = TRUE)
names(true_cd4_counts_065)[1] <- "Day and Stim"
names(true_cd4_counts_065)
names(my_cd4_counts_065) <- gsub(",",
                              "", names(my_cd4_counts_065), fixed = TRUE)
names(my_cd4_counts_065)
which(names(true_cd4_counts_065[,-1]) != names(my_cd4_counts_065[,-c(1,2,3,4)]))



# PID 065 CD8
names(my_cd8_counts_065)
names(true_cd8_counts_065)
names(true_cd8_counts_065) <- gsub("CCR7, Time+ & CD3, Time+ & CD8, Time+ & IFN-g, Time+/Singlets/Lymphocytes/keeper gate/CD3+/CD4-/CD8+/",
                                "", names(true_cd8_counts_065), fixed = TRUE)
names(true_cd8_counts_065) <- gsub(" | Count",
                                "", names(true_cd8_counts_065), fixed = TRUE)
names(true_cd8_counts_065) <- gsub("-",
                                "c ", names(true_cd8_counts_065), fixed = TRUE)
names(true_cd8_counts_065) <- gsub("+",
                                " ", names(true_cd8_counts_065), fixed = TRUE)
names(true_cd8_counts_065) <- gsub("2",
                                "L2", names(true_cd8_counts_065), fixed = TRUE)
names(true_cd8_counts_065) <- gsub("17",
                                "L17", names(true_cd8_counts_065), fixed = TRUE)
names(true_cd8_counts_065) <- gsub("T ",
                                "TNF", names(true_cd8_counts_065), fixed = TRUE)
names(true_cd8_counts_065) <- gsub("Tc ",
                                "TNFc", names(true_cd8_counts_065), fixed = TRUE)
names(true_cd8_counts_065)[1] <- "Day and Stim"
names(true_cd8_counts_065)
names(my_cd8_counts_065) <- gsub(",",
                              "", names(my_cd8_counts_065), fixed = TRUE)
names(my_cd8_counts_065)
which(names(true_cd8_counts_065[,-1]) != names(my_cd8_counts_065[,-c(1,2,3,4)]))



# PID 131 CD4
names(my_cd4_freq_131)
names(true_cd4_freq_131)
names(true_cd4_freq_131) <- gsub("CCR7, Time+ & CD3, Time+ & CD8, Time+ & IFN-g, Time+/Singlets/Lymphocytes/keeper gate/CD3+/CD8-/CD4+/",
                                "", names(true_cd4_freq_131), fixed = TRUE)
names(true_cd4_freq_131) <- gsub(" | Freq. of Parent (%)",
                                "", names(true_cd4_freq_131), fixed = TRUE)
names(true_cd4_freq_131) <- gsub("-",
                                "c ", names(true_cd4_freq_131), fixed = TRUE)
names(true_cd4_freq_131) <- gsub("+",
                                " ", names(true_cd4_freq_131), fixed = TRUE)
names(true_cd4_freq_131) <- gsub("2",
                                "L2", names(true_cd4_freq_131), fixed = TRUE)
names(true_cd4_freq_131) <- gsub("17",
                                "L17", names(true_cd4_freq_131), fixed = TRUE)
names(true_cd4_freq_131) <- gsub("T ",
                                "TNF", names(true_cd4_freq_131), fixed = TRUE)
names(true_cd4_freq_131) <- gsub("Tc ",
                                "TNFc", names(true_cd4_freq_131), fixed = TRUE)
names(true_cd4_freq_131) <- gsub("7",
                                   "R7", names(true_cd4_freq_131), fixed = TRUE)
names(true_cd4_freq_131) <- gsub("45",
                                   "RA", names(true_cd4_freq_131), fixed = TRUE)
names(true_cd4_freq_131) <- gsub("L1R7",
                                 "L17", names(true_cd4_freq_131), fixed = TRUE)
names(true_cd4_freq_131)[1] <- "Day and Stim"
names(true_cd4_freq_131)
names(my_cd4_freq_131) <- gsub(",",
                              "", names(my_cd4_freq_131), fixed = TRUE)
names(my_cd4_freq_131)
which(names(true_cd4_freq_131[,-1]) != names(my_cd4_freq_131[,-c(1,2,3,4)]))



# PID 185 CD4
names(my_cd4_freq_185)
names(true_cd4_freq_185)
names(true_cd4_freq_185) <- gsub("CCR7, Time+ & CD3, Time+ & CD8, Time+ & IFN-g, Time+/Singlets/Lymphocytes/keeper gate/CD3+/CD8-/CD4+/",
                                "", names(true_cd4_freq_185), fixed = TRUE)
names(true_cd4_freq_185) <- gsub(" | Freq. of Parent (%)",
                                "", names(true_cd4_freq_185), fixed = TRUE)
names(true_cd4_freq_185) <- gsub("-",
                                "c ", names(true_cd4_freq_185), fixed = TRUE)
names(true_cd4_freq_185) <- gsub("+",
                                " ", names(true_cd4_freq_185), fixed = TRUE)
names(true_cd4_freq_185) <- gsub("2",
                                "L2", names(true_cd4_freq_185), fixed = TRUE)
names(true_cd4_freq_185) <- gsub("17",
                                "L17", names(true_cd4_freq_185), fixed = TRUE)
names(true_cd4_freq_185) <- gsub("T ",
                                "TNF", names(true_cd4_freq_185), fixed = TRUE)
names(true_cd4_freq_185) <- gsub("Tc ",
                                "TNFc", names(true_cd4_freq_185), fixed = TRUE)
names(true_cd4_freq_185) <- gsub("7",
                                 "R7", names(true_cd4_freq_185), fixed = TRUE)
names(true_cd4_freq_185) <- gsub("45",
                                 "RA", names(true_cd4_freq_185), fixed = TRUE)
names(true_cd4_freq_185) <- gsub("L1R7",
                                 "L17", names(true_cd4_freq_185), fixed = TRUE)
names(true_cd4_freq_185)[1] <- "Day and Stim"
names(true_cd4_freq_185)
names(my_cd4_freq_185) <- gsub(",",
                              "", names(my_cd4_freq_185), fixed = TRUE)
names(my_cd4_freq_185)
which(names(true_cd4_freq_185[,-1]) != names(my_cd4_freq_185[,-c(1,2,3,4)]))








#-------------------------------------------
# Set Values as Numerics
#-------------------------------------------

# PID 063 CD4
true_cd4_counts_063[1,-1] <- as.numeric(sub(",", ".", true_cd4_counts_063[1,-1], fixed = T))
true_cd4_counts_063[2,-1] <- as.numeric(sub(",", ".", true_cd4_counts_063[2,-1], fixed = T))
true_cd4_counts_063[3,-1] <- as.numeric(sub(",", ".", true_cd4_counts_063[3,-1], fixed = T))
true_cd4_counts_063[4,-1] <- as.numeric(sub(",", ".", true_cd4_counts_063[4,-1], fixed = T))

true_cd4_counts_063[5,-1] <- as.numeric(sub(",", ".", true_cd4_counts_063[5,-1], fixed = T))
true_cd4_counts_063[6,-1] <- as.numeric(sub(",", ".", true_cd4_counts_063[6,-1], fixed = T))
true_cd4_counts_063[7,-1] <- as.numeric(sub(",", ".", true_cd4_counts_063[7,-1], fixed = T))
true_cd4_counts_063[8,-1] <- as.numeric(sub(",", ".", true_cd4_counts_063[8,-1], fixed = T))

true_cd4_counts_063[9,-1] <- as.numeric(sub(",", ".", true_cd4_counts_063[9,-1], fixed = T))
true_cd4_counts_063[10,-1] <- as.numeric(sub(",", ".", true_cd4_counts_063[10,-1], fixed = T))
true_cd4_counts_063[11,-1] <- as.numeric(sub(",", ".", true_cd4_counts_063[11,-1], fixed = T))
true_cd4_counts_063[12,-1] <- as.numeric(sub(",", ".", true_cd4_counts_063[12,-1], fixed = T))

# ignore last 2 rows:
true_cd4_counts_063 <- true_cd4_counts_063[-c(13,14),]
true_cd4_counts_063[2:65] <- sapply(true_cd4_counts_063[2:65],as.numeric)
str(true_cd4_counts_063)


# PID 063 CD8
true_cd8_counts_063[1,-1] <- as.numeric(sub(",", ".", true_cd8_counts_063[1,-1], fixed = T))
true_cd8_counts_063[2,-1] <- as.numeric(sub(",", ".", true_cd8_counts_063[2,-1], fixed = T))
true_cd8_counts_063[3,-1] <- as.numeric(sub(",", ".", true_cd8_counts_063[3,-1], fixed = T))
true_cd8_counts_063[4,-1] <- as.numeric(sub(",", ".", true_cd8_counts_063[4,-1], fixed = T))

true_cd8_counts_063[5,-1] <- as.numeric(sub(",", ".", true_cd8_counts_063[5,-1], fixed = T))
true_cd8_counts_063[6,-1] <- as.numeric(sub(",", ".", true_cd8_counts_063[6,-1], fixed = T))
true_cd8_counts_063[7,-1] <- as.numeric(sub(",", ".", true_cd8_counts_063[7,-1], fixed = T))
true_cd8_counts_063[8,-1] <- as.numeric(sub(",", ".", true_cd8_counts_063[8,-1], fixed = T))

true_cd8_counts_063[9,-1] <- as.numeric(sub(",", ".", true_cd8_counts_063[9,-1], fixed = T))
true_cd8_counts_063[10,-1] <- as.numeric(sub(",", ".", true_cd8_counts_063[10,-1], fixed = T))
true_cd8_counts_063[11,-1] <- as.numeric(sub(",", ".", true_cd8_counts_063[11,-1], fixed = T))
true_cd8_counts_063[12,-1] <- as.numeric(sub(",", ".", true_cd8_counts_063[12,-1], fixed = T))

# ignore last 2 rows:
true_cd8_counts_063 <- true_cd8_counts_063[-c(13,14),]
true_cd8_counts_063[2:65] <- sapply(true_cd8_counts_063[2:65],as.numeric)
str(true_cd8_counts_063)




# PID 065 CD4
# already numerics
# ignore last 2 rows:
true_cd4_counts_065 <- true_cd4_counts_065[-c(13,14),]
true_cd4_counts_065[2:65] <- sapply(true_cd4_counts_065[2:65],as.numeric)
str(true_cd4_counts_065)


# PID 065 CD8
true_cd8_counts_065[1,-1] <- as.numeric(sub(",", ".", true_cd8_counts_065[1,-1], fixed = T))
true_cd8_counts_065[2,-1] <- as.numeric(sub(",", ".", true_cd8_counts_065[2,-1], fixed = T))
true_cd8_counts_065[3,-1] <- as.numeric(sub(",", ".", true_cd8_counts_065[3,-1], fixed = T))
true_cd8_counts_065[4,-1] <- as.numeric(sub(",", ".", true_cd8_counts_065[4,-1], fixed = T))

true_cd8_counts_065[5,-1] <- as.numeric(sub(",", ".", true_cd8_counts_065[5,-1], fixed = T))
true_cd8_counts_065[6,-1] <- as.numeric(sub(",", ".", true_cd8_counts_065[6,-1], fixed = T))
true_cd8_counts_065[7,-1] <- as.numeric(sub(",", ".", true_cd8_counts_065[7,-1], fixed = T))
true_cd8_counts_065[8,-1] <- as.numeric(sub(",", ".", true_cd8_counts_065[8,-1], fixed = T))

true_cd8_counts_065[9,-1] <- as.numeric(sub(",", ".", true_cd8_counts_065[9,-1], fixed = T))
true_cd8_counts_065[10,-1] <- as.numeric(sub(",", ".", true_cd8_counts_065[10,-1], fixed = T))
true_cd8_counts_065[11,-1] <- as.numeric(sub(",", ".", true_cd8_counts_065[11,-1], fixed = T))
true_cd8_counts_065[12,-1] <- as.numeric(sub(",", ".", true_cd8_counts_065[12,-1], fixed = T))

# ignore last 2 rows:
true_cd8_counts_065 <- true_cd8_counts_065[-c(13,14),]
true_cd8_counts_065[2:65] <- sapply(true_cd8_counts_065[2:65],as.numeric)
str(true_cd8_counts_065)



# PID 131 CD4
str(true_cd4_freq_131)
true_cd4_freq_131[1,-1] <- as.numeric(sub(",", ".", true_cd4_freq_131[1,-1], fixed = T))
true_cd4_freq_131[2,-1] <- as.numeric(sub(",", ".", true_cd4_freq_131[2,-1], fixed = T))
true_cd4_freq_131[3,-1] <- as.numeric(sub(",", ".", true_cd4_freq_131[3,-1], fixed = T))
true_cd4_freq_131[4,-1] <- as.numeric(sub(",", ".", true_cd4_freq_131[4,-1], fixed = T))

true_cd4_freq_131[5,-1] <- as.numeric(sub(",", ".", true_cd4_freq_131[5,-1], fixed = T))
true_cd4_freq_131[6,-1] <- as.numeric(sub(",", ".", true_cd4_freq_131[6,-1], fixed = T))
true_cd4_freq_131[7,-1] <- as.numeric(sub(",", ".", true_cd4_freq_131[7,-1], fixed = T))
true_cd4_freq_131[8,-1] <- as.numeric(sub(",", ".", true_cd4_freq_131[8,-1], fixed = T))

true_cd4_freq_131[9,-1] <- as.numeric(sub(",", ".", true_cd4_freq_131[9,-1], fixed = T))
true_cd4_freq_131[10,-1] <- as.numeric(sub(",", ".", true_cd4_freq_131[10,-1], fixed = T))
true_cd4_freq_131[11,-1] <- as.numeric(sub(",", ".", true_cd4_freq_131[11,-1], fixed = T))
true_cd4_freq_131[12,-1] <- as.numeric(sub(",", ".", true_cd4_freq_131[12,-1], fixed = T))

# ignore last 2 rows:
true_cd4_freq_131 <- true_cd4_freq_131[-c(13,14),]
true_cd4_freq_131[2:65] <- sapply(true_cd4_freq_131[2:65],as.numeric)
str(true_cd4_freq_131)


# PID 185 CD4
str(true_cd4_freq_185)
true_cd4_freq_185[1,-1] <- as.numeric(sub(",", ".", true_cd4_freq_185[1,-1], fixed = T))
true_cd4_freq_185[2,-1] <- as.numeric(sub(",", ".", true_cd4_freq_185[2,-1], fixed = T))
true_cd4_freq_185[3,-1] <- as.numeric(sub(",", ".", true_cd4_freq_185[3,-1], fixed = T))
true_cd4_freq_185[4,-1] <- as.numeric(sub(",", ".", true_cd4_freq_185[4,-1], fixed = T))

true_cd4_freq_185[5,-1] <- as.numeric(sub(",", ".", true_cd4_freq_185[5,-1], fixed = T))
true_cd4_freq_185[6,-1] <- as.numeric(sub(",", ".", true_cd4_freq_185[6,-1], fixed = T))
true_cd4_freq_185[7,-1] <- as.numeric(sub(",", ".", true_cd4_freq_185[7,-1], fixed = T))
true_cd4_freq_185[8,-1] <- as.numeric(sub(",", ".", true_cd4_freq_185[8,-1], fixed = T))

true_cd4_freq_185[9,-1] <- as.numeric(sub(",", ".", true_cd4_freq_185[9,-1], fixed = T))
true_cd4_freq_185[10,-1] <- as.numeric(sub(",", ".", true_cd4_freq_185[10,-1], fixed = T))
true_cd4_freq_185[11,-1] <- as.numeric(sub(",", ".", true_cd4_freq_185[11,-1], fixed = T))
true_cd4_freq_185[12,-1] <- as.numeric(sub(",", ".", true_cd4_freq_185[12,-1], fixed = T))

# ignore last 2 rows:
true_cd4_freq_185 <- true_cd4_freq_185[-c(13,14),]
true_cd4_freq_185[2:65] <- sapply(true_cd4_freq_185[2:65],as.numeric)
str(true_cd4_freq_185)


#------------------------------------
#  PID 063 CD4
#------------------------------------

my_comb_dat <- c(as.numeric(my_cd4_counts_063[1,-c(1,2,3,4)]),as.numeric(my_cd4_counts_063[2,-c(1,2,3,4)]),
                 as.numeric(my_cd4_counts_063[3,-c(1,2,3,4)]), as.numeric(my_cd4_counts_063[4,-c(1,2,3,4)]),
                 as.numeric(my_cd4_counts_063[5,-c(1,2,3,4)]), as.numeric(my_cd4_counts_063[6,-c(1,2,3,4)]),
                 as.numeric(my_cd4_counts_063[7,-c(1,2,3,4)]),as.numeric(my_cd4_counts_063[8,-c(1,2,3,4)]),
                 as.numeric(my_cd4_counts_063[9,-c(1,2,3,4)]),as.numeric(my_cd4_counts_063[10,-c(1,2,3,4)]),
                 as.numeric(my_cd4_counts_063[11,-c(1,2,3,4)]),as.numeric(my_cd4_counts_063[12,-c(1,2,3,4)]))

true_comb_dat <- c(as.numeric(true_cd4_counts_063[2,-1]),as.numeric(true_cd4_counts_063[5,-1]), 
                   as.numeric(true_cd4_counts_063[8,-1]),as.numeric(true_cd4_counts_063[12,-1]),
                   as.numeric(true_cd4_counts_063[3,-1]),as.numeric(true_cd4_counts_063[6,-1]),
                   as.numeric(true_cd4_counts_063[9,-1]),as.numeric(true_cd4_counts_063[11,-1]),
                   as.numeric(true_cd4_counts_063[1,-1]), as.numeric(true_cd4_counts_063[4,-1]),
                   as.numeric(true_cd4_counts_063[7,-1]),as.numeric(true_cd4_counts_063[10,-1]))


ccc <- CCC(my_comb_dat, true_comb_dat)$rho.c

# Plot Values
options(scipen=999)
par(mar=c(5.5, 7, 4.5, 2.5))
plot(x=my_comb_dat, y=true_comb_dat,
     xlab="", ylab= "", main = "Pid 063", pch = 19, cex.main=5, cex.lab=4, cex.axis=3)
title(ylab="FlowJo CD4 Counts", line=4, cex.lab=4)
title(xlab="R CD4 Counts", line=4.4, cex.lab=4)
grid()
abline(a=0,b=1, col="red",lty=2)
mtext(paste0("CCC = ", round(ccc[1],6)), side=3, cex = 4, line=-3.5)


# Plot Errors
x <- my_comb_dat - true_comb_dat
options(scipen=999)
par(mar=c(5.5, 7, 4.5, 2.5))
plot(x,
     xlab="", ylab= "", main = "Pid 063", pch = 19, cex.main=5, cex.lab=4, cex.axis=3)
title(ylab="CD4 Count Errors (R-FlowJo)", line=4, cex.lab=4)
grid()
abline(a=0,b=0, col="red",lty=2)


# Plot Error Percentages 
x <- abs(my_comb_dat - true_comb_dat)/true_comb_dat
x <- na.omit(x) # ignore when true count  == 0
plot(x,
     xlab="", ylab= "", main = "Pid 063", pch = 19, cex.main=5, cex.lab=4, cex.axis=3)
title(ylab="CD4 Count Errors (%)", line=4, cex.lab=4)
grid()
abline(a=0,b=0, col="red",lty=2)

y <- which(true_comb_dat==0)
z <- which(my_comb_dat[y]!=0)
z
my_comb_dat[y][z]







#------------------------------------
#  PID 063 CD8
#------------------------------------

my_comb_dat <- c(as.numeric(my_cd8_counts_063[1,-c(1,2,3,4)]),as.numeric(my_cd8_counts_063[2,-c(1,2,3,4)]),
                 as.numeric(my_cd8_counts_063[3,-c(1,2,3,4)]), as.numeric(my_cd8_counts_063[4,-c(1,2,3,4)]),
                 as.numeric(my_cd8_counts_063[5,-c(1,2,3,4)]), as.numeric(my_cd8_counts_063[6,-c(1,2,3,4)]),
                 as.numeric(my_cd8_counts_063[7,-c(1,2,3,4)]),as.numeric(my_cd8_counts_063[8,-c(1,2,3,4)]),
                 as.numeric(my_cd8_counts_063[9,-c(1,2,3,4)]),as.numeric(my_cd8_counts_063[10,-c(1,2,3,4)]),
                 as.numeric(my_cd8_counts_063[11,-c(1,2,3,4)]),as.numeric(my_cd8_counts_063[12,-c(1,2,3,4)]))

true_comb_dat <- c(as.numeric(true_cd8_counts_063[2,-1]),as.numeric(true_cd8_counts_063[5,-1]), 
                   as.numeric(true_cd8_counts_063[8,-1]),as.numeric(true_cd8_counts_063[12,-1]),
                   as.numeric(true_cd8_counts_063[3,-1]),as.numeric(true_cd8_counts_063[6,-1]),
                   as.numeric(true_cd8_counts_063[9,-1]),as.numeric(true_cd8_counts_063[11,-1]),
                   as.numeric(true_cd8_counts_063[1,-1]), as.numeric(true_cd8_counts_063[4,-1]),
                   as.numeric(true_cd8_counts_063[7,-1]),as.numeric(true_cd8_counts_063[10,-1]))


ccc <- CCC(my_comb_dat, true_comb_dat)$rho.c

# Plot Values
plot(x=my_comb_dat, y=true_comb_dat,
     xlab= "R CD8 Counts", ylab= "FlowJo CD8 Counts", pch = 19, main = "Pid 063")
grid()
abline(a=0,b=1, col="red",lty=2)
mtext(paste0("CCC = ", round(ccc[1],8)), side=3)


# Plot Errors
x <- my_comb_dat - true_comb_dat
plot(x, pch = 19, main = "Pid 063", ylab="CD8 Count Errors (R-FlowJo)")
grid()
abline(a=0,b=0, col="red",lty=2)


# Plot Error Percentages 
x <- abs(my_comb_dat - true_comb_dat)/true_comb_dat
x <- na.omit(x) # ignore when true count  == 0
plot(x, pch = 19, main = "Pid 063", ylab="CD8 Count Errors (%)")
grid()
abline(a=0,b=0, col="red",lty=2)

y <- which(true_comb_dat==0)
z <- which(my_comb_dat[y]!=0)
z
my_comb_dat[y][z]






#------------------------------------
#  PID 065 CD4
#------------------------------------

my_comb_dat <- c(as.numeric(my_cd4_counts_065[1,-c(1,2,3,4)]),as.numeric(my_cd4_counts_065[2,-c(1,2,3,4)]),
                 as.numeric(my_cd4_counts_065[3,-c(1,2,3,4)]), as.numeric(my_cd4_counts_065[4,-c(1,2,3,4)]),
                 as.numeric(my_cd4_counts_065[5,-c(1,2,3,4)]), as.numeric(my_cd4_counts_065[6,-c(1,2,3,4)]),
                 as.numeric(my_cd4_counts_065[7,-c(1,2,3,4)]),as.numeric(my_cd4_counts_065[8,-c(1,2,3,4)]),
                 as.numeric(my_cd4_counts_065[9,-c(1,2,3,4)]),as.numeric(my_cd4_counts_065[10,-c(1,2,3,4)]),
                 as.numeric(my_cd4_counts_065[11,-c(1,2,3,4)]),as.numeric(my_cd4_counts_065[12,-c(1,2,3,4)]))

true_comb_dat <- c(as.numeric(true_cd4_counts_065[3,-1]),as.numeric(true_cd4_counts_065[6,-1]), 
                   as.numeric(true_cd4_counts_065[8,-1]),as.numeric(true_cd4_counts_065[12,-1]),
                   as.numeric(true_cd4_counts_065[1,-1]),as.numeric(true_cd4_counts_065[6,-1]),
                   as.numeric(true_cd4_counts_065[7,-1]),as.numeric(true_cd4_counts_065[10,-1]),
                   as.numeric(true_cd4_counts_065[2,-1]), as.numeric(true_cd4_counts_065[5,-1]),
                   as.numeric(true_cd4_counts_065[9,-1]),as.numeric(true_cd4_counts_065[11,-1]))


ccc <- CCC(my_comb_dat, true_comb_dat)$rho.c

# Plot Values
plot(x=my_comb_dat, y=true_comb_dat,
     xlab= "R CD4 Counts", ylab= "FlowJo CD4 Counts", pch = 19, main = "Pid 065")
grid()
abline(a=0,b=1, col="red",lty=2)
mtext(paste0("CCC = ", round(ccc[1],8)), side=3)


# Plot Errors
x <- my_comb_dat - true_comb_dat
plot(x, pch = 19, main = "Pid 065", ylab="CD4 Count Errors (R-FlowJo)")
grid()
abline(a=0,b=0, col="red",lty=2)


# Plot Error Percentages 
x <- abs(my_comb_dat - true_comb_dat)/true_comb_dat
x <- na.omit(x) # ignore when true count  == 0
plot(x, pch = 19, main = "Pid 065", ylab="CD4 Count Errors (%)")
grid()
abline(a=0,b=0, col="red",lty=2)

y <- which(true_comb_dat==0)
z <- which(my_comb_dat[y]!=0)
z
my_comb_dat[y][z]



#------------------------------------
#  PID 065 CD8
#------------------------------------

my_comb_dat <- c(as.numeric(my_cd8_counts_065[1,-c(1,2,3,4)]),as.numeric(my_cd8_counts_065[2,-c(1,2,3,4)]),
                 as.numeric(my_cd8_counts_065[3,-c(1,2,3,4)]), as.numeric(my_cd8_counts_065[4,-c(1,2,3,4)]),
                 as.numeric(my_cd8_counts_065[5,-c(1,2,3,4)]), as.numeric(my_cd8_counts_065[6,-c(1,2,3,4)]),
                 as.numeric(my_cd8_counts_065[7,-c(1,2,3,4)]),as.numeric(my_cd8_counts_065[8,-c(1,2,3,4)]),
                 as.numeric(my_cd8_counts_065[9,-c(1,2,3,4)]),as.numeric(my_cd8_counts_065[10,-c(1,2,3,4)]),
                 as.numeric(my_cd8_counts_065[11,-c(1,2,3,4)]),as.numeric(my_cd8_counts_065[12,-c(1,2,3,4)]))

true_comb_dat <- c(as.numeric(true_cd8_counts_065[3,-1]),as.numeric(true_cd8_counts_065[6,-1]), 
                   as.numeric(true_cd8_counts_065[8,-1]),as.numeric(true_cd8_counts_065[12,-1]),
                   as.numeric(true_cd8_counts_065[1,-1]),as.numeric(true_cd8_counts_065[6,-1]),
                   as.numeric(true_cd8_counts_065[7,-1]),as.numeric(true_cd8_counts_065[10,-1]),
                   as.numeric(true_cd8_counts_065[2,-1]), as.numeric(true_cd8_counts_065[5,-1]),
                   as.numeric(true_cd8_counts_065[9,-1]),as.numeric(true_cd8_counts_065[11,-1]))


ccc <- CCC(my_comb_dat, true_comb_dat)$rho.c

# Plot Values
plot(x=my_comb_dat, y=true_comb_dat,
     xlab= "R CD8 Counts", ylab= "FlowJo CD8 Counts", pch = 19, main = "Pid 065")
grid()
abline(a=0,b=1, col="red",lty=2)
mtext(paste0("CCC = ", round(ccc[1],8)), side=3)


# Plot Errors
x <- my_comb_dat - true_comb_dat
plot(x, pch = 19, main = "Pid 065", ylab="CD8 Count Errors (R-FlowJo)")
grid()
abline(a=0,b=0, col="red",lty=2)


# Plot Error Percentages 
x <- abs(my_comb_dat - true_comb_dat)/true_comb_dat
x <- na.omit(x) # ignore when true count  == 0
plot(x, pch = 19, main = "Pid 065", ylab="CD8 Count Errors (%)")
grid()
abline(a=0,b=0, col="red",lty=2)

y <- which(true_comb_dat==0)
z <- which(my_comb_dat[y]!=0)
z
my_comb_dat[y][z]







#------------------------------------
#  PID 131 CD4
#------------------------------------

my_comb_dat <- c(as.numeric(my_cd4_freq_131[1,-c(1,2,3,4)]),as.numeric(my_cd4_freq_131[2,-c(1,2,3,4)]),
                 as.numeric(my_cd4_freq_131[3,-c(1,2,3,4)]), as.numeric(my_cd4_freq_131[4,-c(1,2,3,4)]),
                 as.numeric(my_cd4_freq_131[5,-c(1,2,3,4)]), as.numeric(my_cd4_freq_131[6,-c(1,2,3,4)]),
                 as.numeric(my_cd4_freq_131[7,-c(1,2,3,4)]),as.numeric(my_cd4_freq_131[8,-c(1,2,3,4)]),
                 as.numeric(my_cd4_freq_131[9,-c(1,2,3,4)]),as.numeric(my_cd4_freq_131[10,-c(1,2,3,4)]),
                 as.numeric(my_cd4_freq_131[11,-c(1,2,3,4)]),as.numeric(my_cd4_freq_131[12,-c(1,2,3,4)]))

true_comb_dat <- c(as.numeric(true_cd4_freq_131[6,-1]),as.numeric(true_cd4_freq_131[9,-1]), 
                   as.numeric(true_cd4_freq_131[12,-1]),as.numeric(true_cd4_freq_131[3,-1]),
                   as.numeric(true_cd4_freq_131[4,-1]),as.numeric(true_cd4_freq_131[7,-1]),
                   as.numeric(true_cd4_freq_131[10,-1]),as.numeric(true_cd4_freq_131[1,-1]),
                   as.numeric(true_cd4_freq_131[5,-1]), as.numeric(true_cd4_freq_131[8,-1]),
                   as.numeric(true_cd4_freq_131[11,-1]),as.numeric(true_cd4_freq_131[2,-1]))


ccc <- CCC(my_comb_dat, true_comb_dat)$rho.c



options(scipen=999)
par(mar=c(5.5, 7, 4.5, 2.5))
plot(x=my_comb_dat, y=true_comb_dat,
     xlab="", ylab= "", main = "Pid 131", pch = 19, cex.main=5, cex.lab=4, cex.axis=3)
title(ylab="FlowJo CD4 Frequencies", line=4, cex.lab=4)
title(xlab="R CD4 Frequencies", line=4.4, cex.lab=4)
grid()
abline(a=0,b=1, col="red",lty=2)
mtext(paste0("CCC = ", round(ccc[1],6)), side=3, cex = 4, line=-3.5)


# Plot Errors
x <- my_comb_dat - true_comb_dat
options(scipen=999)
par(mar=c(5.5, 7, 4.5, 2.5))
plot(x,
     xlab="", ylab= "", main = "Pid 131", pch = 19, cex.main=5, cex.lab=4, cex.axis=3)
title(ylab="CD4 Frequency Errors (R-FlowJo)", line=4, cex.lab=4)
grid()
abline(a=0,b=0, col="red",lty=2)


# Plot Error Percentages 
x <- abs(my_comb_dat - true_comb_dat)/true_comb_dat
x <- na.omit(x) # ignore when true count  == 0
plot(x,
     xlab="", ylab= "", main = "Pid 131", pch = 19, cex.main=5, cex.lab=4, cex.axis=3)
title(ylab="CD4 Frequency Errors (%)", line=4, cex.lab=4)
grid()
abline(a=0,b=0, col="red",lty=2)



y <- which(true_comb_dat==0)
z <- which(my_comb_dat[y]!=0)
z
my_comb_dat[y][z]






#------------------------------------
#  PID 185 CD4
#------------------------------------

my_comb_dat <- c(as.numeric(my_cd4_freq_185[1,-c(1,2,3,4)]),as.numeric(my_cd4_freq_185[2,-c(1,2,3,4)]),
                 as.numeric(my_cd4_freq_185[3,-c(1,2,3,4)]), 
                 as.numeric(my_cd4_freq_185[4,-c(1,2,3,4)]), as.numeric(my_cd4_freq_185[5,-c(1,2,3,4)]),
                 as.numeric(my_cd4_freq_185[6,-c(1,2,3,4)]),as.numeric(my_cd4_freq_185[7,-c(1,2,3,4)]),
                 as.numeric(my_cd4_freq_185[8,-c(1,2,3,4)]),as.numeric(my_cd4_freq_185[9,-c(1,2,3,4)]),
                 as.numeric(my_cd4_freq_185[10,-c(1,2,3,4)]),as.numeric(my_cd4_freq_185[11,-c(1,2,3,4)]))

true_comb_dat <- c(as.numeric(true_cd4_freq_185[3,-1]),as.numeric(true_cd4_freq_185[4,-1]), 
                   as.numeric(true_cd4_freq_185[7,-1]),
                   as.numeric(true_cd4_freq_185[1,-1]),as.numeric(true_cd4_freq_185[5,-1]),
                   as.numeric(true_cd4_freq_185[8,-1]),as.numeric(true_cd4_freq_185[11,-1]),
                   as.numeric(true_cd4_freq_185[2,-1]), as.numeric(true_cd4_freq_185[6,-1]),
                   as.numeric(true_cd4_freq_185[9,-1]),as.numeric(true_cd4_freq_185[10,-1]))


ccc <- CCC(my_comb_dat, true_comb_dat)$rho.c

options(scipen=999)
par(mar=c(5.5, 7, 4.5, 2.5))
plot(x=my_comb_dat, y=true_comb_dat,
     xlab="", ylab= "", main = "Pid 185", pch = 19, cex.main=5, cex.lab=4, cex.axis=3)
title(ylab="FlowJo CD4 Frequencies", line=4, cex.lab=4)
title(xlab="R CD4 Frequencies", line=4.4, cex.lab=4)
grid()
abline(a=0,b=1, col="red",lty=2)
mtext(paste0("CCC = ", round(ccc[1],6)), side=3, cex = 4, line=-3.5)


# Plot Errors
x <- my_comb_dat - true_comb_dat
options(scipen=999)
par(mar=c(5.5, 7, 4.5, 2.5))
plot(x,
     xlab="", ylab= "", main = "Pid 185", pch = 19, cex.main=5, cex.lab=4, cex.axis=3)
title(ylab="CD4 Frequency Errors (R-FlowJo)", line=4, cex.lab=4)
grid()
abline(a=0,b=0, col="red",lty=2)


# Plot Error Percentages 
x <- abs(my_comb_dat - true_comb_dat)/true_comb_dat
x <- na.omit(x) # ignore when true count  == 0
plot(x,
     xlab="", ylab= "", main = "Pid 185", pch = 19, cex.main=5, cex.lab=4, cex.axis=3)
title(ylab="CD4 Frequency Errors (%)", line=4, cex.lab=4)
grid()
abline(a=0,b=0, col="red",lty=2)

y <- which(true_comb_dat==0)
z <- which(my_comb_dat[y]!=0)
z









