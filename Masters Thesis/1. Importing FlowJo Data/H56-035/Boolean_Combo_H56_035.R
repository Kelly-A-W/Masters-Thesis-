# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()


#---------------------
# Load packages
#---------------------
# for the tricky packages, ALWAYS choose 'a' or 'Yes' when applicable
# ignore 'force = TRUE' 

library(devtools) 
# devtools::install_github("RGLab/flowCore")
if (!requireNamespace("flowCore", quietly = TRUE)) {
  BiocManager::install("flowCore")   # need for fcs data
}
library(flowCore)
if (!requireNamespace("flowWorkspaceData", quietly = TRUE)) {
  BiocManager::install("flowWorkspaceData")  # need for importing workspace
}
library(flowWorkspaceData)
if (!requireNamespace("CytoML", quietly = TRUE)) {
  BiocManager::install("CytoML") # need for importing workspace
}
library(CytoML)
library(flowWorkspace)                    # need to get gates and compensation matrix
library(sp)
library(ggcyto)
library(tibble)
library(dplyr)



#------------------------------------------
# DEFINE SUBSETS / COMBINATIONS OF MARKERS
#-----------------------------------------

marker <- c("R7", "R7c", "RA", "RAc", "G", "Gc", "L2", "L2c", "L17", "L17c", "TNF", "TNFc")
marker_combo <- combn(marker, 6)
marker_combo <- t(marker_combo)
op1 <- which(apply(marker_combo, 1, function(r) any(r %in% c("R7", "R7c")))) # if neither R7 or R7c is in a combo, that means there must be at least one other cytokine and its complement in the combo
op2 <- which(apply(marker_combo, 1, function(r) any(r %in% c("RA", "RAc"))))
op3 <- which(apply(marker_combo, 1, function(r) any(r %in% c("G", "Gc"))))
op4 <- which(apply(marker_combo, 1, function(r) any(r %in% c("L2", "L2c"))))
op5 <- which(apply(marker_combo, 1, function(r) any(r %in% c("L17", "L17c"))))
op6 <- which(apply(marker_combo, 1, function(r) any(r %in% c("TNF", "TNFc"))))
int <- Reduce(intersect, list(op1,op2,op3,op4,op5,op6))
comb_final <- marker_combo[int,]




#-------------------------
# BOOLEAN GATING FUNCTION
#-------------------------


combo_cd4 <- function(ws, pid, days, comb=comb_final){
  out <- data.frame()
  
  for(i in 1:length(days)){
    gates <- flowjo_to_gatingset(ws, name = (i+2))
    day   <- days[i]
    stim  <- stringr::str_extract(sampleNames(gates), "UNS|Ag85B|Esat6|PHA")
    
    for(j in 1:length(stim)){
      
      if("try-error" %in% class(try(gates[[j]], silent = T))){
        x  <- matrix(NA, nrow = 1, ncol = nrow(comb))
        x  <- as.data.frame(x)
        for(k in 1:nrow(comb)){
          colnames(x)[k] <- toString(comb[k,])
        }
        y  <- cbind(PID=pid, 
                    Day=day,
                    Stim=stim[j],
                    CD4count=NA,
                    x)
      }else{
        
        gh     <- gates[[j]]
        cd4    <- gh_pop_get_indices(gh, "CD4+/")
        cd4_in <- which(cd4 == T)
        R7     <- gh_pop_get_indices(gh, "CD4+/CCR7+")[cd4_in]
        RA     <- gh_pop_get_indices(gh, "CD4+/CD45RA+")[cd4_in]
        G      <- gh_pop_get_indices(gh, "CD4+/IFN-g+")[cd4_in]
        L2     <- gh_pop_get_indices(gh, "CD4+/IL-2+")[cd4_in]
        L17    <- gh_pop_get_indices(gh, "CD4+/IL17+")[cd4_in]
        TNF    <- gh_pop_get_indices(gh, "CD4+/TNF+")[cd4_in]
        df     <- cbind(R7, RA, G, L2, L17, TNF)
        df     <- as.data.frame(df)
        
        df_big      <- df
        df_big$R7c  <- !df$R7
        df_big$RAc  <- !df$RA
        df_big$Gc   <- !df$G
        df_big$L2c  <- !df$L2
        df_big$L17c <- !df$L17
        df_big$TNFc <- !df$TNF
        
        x <- matrix(NA, nrow = 1, ncol = nrow(comb))
        x <- as.data.frame(x)
        
        for(k in 1:nrow(comb)){
          count          <- length(which(rowSums(df_big[,comb[k,]]) == 6))
          colnames(x)[k] <- toString(comb[k,])
          x[1,k]         <- count
        }
        
        y      <- cbind(PID=pid, 
                        Day=day,
                        Stim=stim[j],
                        CD4count=gh_pop_get_count(gh, "CD4+"),
                        x)
      }
      
      out    <- rbind(out, y)
      
    }
    
  }
  out
}













combo_cd8 <- function(ws, pid, days, comb=comb_final){
  out <- data.frame()
  
  for(i in 1:length(days)){
    gates <- flowjo_to_gatingset(ws, name = (i+2))
    day   <- days[i]
    stim  <- stringr::str_extract(sampleNames(gates), "UNS|Ag85B|Esat6|PHA")
    
    
    for(j in 1:length(stim)){
      
      if("try-error" %in% class(try(gates[[j]], silent = T))){
        x  <- matrix(NA, nrow = 1, ncol = nrow(comb))
        x  <- as.data.frame(x)
        for(k in 1:nrow(comb)){
          colnames(x)[k] <- toString(comb[k,])
        }
        y  <- cbind(PID=pid, 
                    Day=day,
                    Stim=stim[j],
                    CD8count=NA,
                    x)
      }else{
        
        gh     <- gates[[j]]
        cd8    <- gh_pop_get_indices(gh, "CD8+/")
        cd8_in <- which(cd8 == T)
        R7     <- gh_pop_get_indices(gh, "CD8+/CCR7+")[cd8_in]
        RA     <- gh_pop_get_indices(gh, "CD8+/CD45RA+")[cd8_in]
        G      <- gh_pop_get_indices(gh, "CD8+/IFN-g+")[cd8_in]
        L2     <- gh_pop_get_indices(gh, "CD8+/IL-2+")[cd8_in]
        L17    <- gh_pop_get_indices(gh, "CD8+/IL17+")[cd8_in]
        TNF    <- gh_pop_get_indices(gh, "CD8+/TNF+")[cd8_in]
        df     <- cbind(R7, RA, G, L2, L17, TNF)
        df     <- as.data.frame(df)
        
        df_big      <- df
        df_big$R7c  <- !df$R7
        df_big$RAc  <- !df$RA
        df_big$Gc   <- !df$G
        df_big$L2c  <- !df$L2
        df_big$L17c <- !df$L17
        df_big$TNFc <- !df$TNF
        
        x <- matrix(NA, nrow = 1, ncol = nrow(comb))
        x <- as.data.frame(x)
        
        for(k in 1:nrow(comb)){
          count          <- length(which(rowSums(df_big[,comb[k,]]) == 6))
          colnames(x)[k] <- toString(comb[k,])
          x[1,k]         <- count
        }
        
        y      <- cbind(PID=pid, 
                        Day=day,
                        Stim=stim[j],
                        CD8count=gh_pop_get_count(gh, "CD8+"),
                        x)
      }
      
      out    <- rbind(out, y)
      
    }
    
  }
  out
}




#------------------------------------
# GET MATRICES FOR CD4+ AND CD8+
#------------------------------------

my_days <- c(0,70)


setwd("/Volumes/Elements/WBA Workspace & Fcs Final/H56-035/G1")

# UNS|Ag85B|ESAT6|PHA
my_pid  <- "1034"                               ####### CHANGE #######
#my_ws   <- open_flowjo_xml(paste0(my_pid, "_AG.wsp")) 
my_ws   <- open_flowjo_xml("035-456 01034 04082014_AG_Final.wsp")
#my_ws   <- open_flowjo_xml("3003_AG.wsp")
my_ws 



counts_cd4 <- combo_cd4(ws=my_ws, pid=my_pid, days=my_days)
freq_cd4   <- counts_cd4
freq_cd4   <- freq_cd4 %>%
  mutate(across(5:68, ~ ./counts_cd4$CD4count))
freq_cd4[,5:68] <- freq_cd4[,5:68]*100  


counts_cd8 <- combo_cd8(ws=my_ws, pid=my_pid, days=my_days)
freq_cd8   <- counts_cd8
freq_cd8   <- freq_cd8 %>%
  mutate(across(5:68, ~ ./counts_cd8$CD8count))
freq_cd8[,5:68] <- freq_cd8[,5:68]*100  



# Relabel Esat6 to ESAT6 for consistency 
counts_cd4$Stim <- as.factor(counts_cd4$Stim)
levels(counts_cd4$Stim)[2] <- "ESAT6"
freq_cd4$Stim <- as.factor(freq_cd4$Stim)
levels(freq_cd4$Stim)[2] <- "ESAT6"
counts_cd8$Stim <- as.factor(counts_cd8$Stim)
levels(counts_cd8$Stim)[2] <- "ESAT6"
freq_cd8$Stim <- as.factor(freq_cd8$Stim)
levels(freq_cd8$Stim)[2] <- "ESAT6"


#------------------------------------
# EXPORT TO CSV FILE
#------------------------------------


setwd("/Users/kellywilliams/Documents/Masters Thesis/FlowJo Booleans Dataset 1/H56-035/Data")

write.csv(counts_cd4, paste0("cd4_counts_pid_",my_pid, ".csv"), row.names = FALSE) 
write.csv(freq_cd4, paste0("cd4_frequency_pid_",my_pid, ".csv"), row.names = FALSE) 

write.csv(counts_cd8, paste0("cd8_counts_pid_",my_pid, ".csv"), row.names = FALSE)  
write.csv(freq_cd8, paste0("cd8_frequency_pid_",my_pid, ".csv"), row.names = FALSE) 




