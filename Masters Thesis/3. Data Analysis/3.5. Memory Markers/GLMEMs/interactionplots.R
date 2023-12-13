# Clear environment
rm(list=ls())
# Clear all plots
if(!is.null(dev.list())) dev.off()

#---------------------
# Load Libraries
#---------------------
library(ggplot2)



############################################
# Interaction Plot for CCR7-CD45RA- ESAT-6:
############################################
qft_pos = (0.303-0.611)*c(1,2,3)+2.984
qft_neg = (0.303)*c(1,2,3)
admin = rep(c(1,2,3),2)
df = data.frame(value = c(qft_pos, qft_neg),
                admin = admin,
                qft = c(rep("QFT+",3), rep("QFT-",3)))
df$qft = as.factor(df$qft)
df$admin = as.factor(df$admin)

ggplot(df, aes(y=value, x=admin, group = qft)) +
  geom_line(aes(color=qft),linewidth=2)+
  geom_point(size=5, aes(color = qft))+
  theme_minimal()+
  theme(text = element_text(size = 28)) +
  xlab("Number of administrations")+
  ylab(expression(paste("log(", mu[ij], ")")))+
  scale_color_manual(values=c("#FF0000", "blue"), name = NULL)


############################################
# Interaction Plot for CCR7+CD45RA- Ag85B:
############################################
qft_pos = c(0, 2.015-0.788, 1.625-0.875) + 1.260
qft_neg = c(0, 2.015, 1.625)
time = rep(c("baseline", "peak", "memory"),2)
df = data.frame(value = c(qft_pos, qft_neg),
                time = time,
                qft = c(rep("QFT+",3), rep("QFT-",3)))
df$qft = as.factor(df$qft)
df$time <- factor(df$time, levels=c("baseline", "peak", "memory"))
levels(df$time)

ggplot(df, aes(y=value, x=time, group = qft)) +
  geom_line(aes(color=qft),linewidth=2)+
  geom_point(size=5, aes(color = qft))+
  xlab(NULL)+
  ylab(expression(paste("log(", mu[ij], ")")))+
  theme_minimal()+
  theme(text = element_text(size = 28))+
  scale_color_manual(values=c("#FF0000", "blue"), name = NULL)




############################################
# Interaction Plot for CCR7-CD45RA- Ag85B:
############################################
qft_pos = c(0, 2.098-0.801, 1.293-0.691) + 1.185
qft_neg = c(0, 2.098, 1.293)
time = rep(c("baseline", "peak", "memory"),2)
df = data.frame(value = c(qft_pos, qft_neg),
                time = time,
                qft = c(rep("QFT+",3), rep("QFT-",3)))
df$qft = as.factor(df$qft)
df$time <- factor(df$time, levels=c("baseline", "peak", "memory"))
levels(df$time)

ggplot(df, aes(y=value, x=time, group = qft)) +
  geom_line(aes(color=qft),linewidth=2)+
  geom_point(size=5, aes(color = qft))+
  theme_minimal()+
  theme(text = element_text(size = 28)) +
  xlab(NULL)+
  ylab(expression(paste("log(", mu[ij], ")")))+
  scale_color_manual(values=c("#FF0000", "blue"), name = NULL)




###############################################
# Interaction Plot for CCR7+CD45RA- BOTH Stim:
###############################################
esat6_qft_pos = c(0, 2.029-0.993-0.877, 1.556-0.930-0.730 )+1.365-1.142+1.370
esat6_qft_neg = c(0,2.029-0.877, 1.556-0.730)-1.142
ag85b_qft_pos = c(0, 2.029-0.993 , 1.556-0.930 )+1.365
ag85b_qft_neg =c(0, 2.029, 1.556)
time = rep(c("baseline", "peak", "memory"),4)

df = data.frame(value = c(esat6_qft_pos, esat6_qft_neg, ag85b_qft_pos, ag85b_qft_neg),
                time = time,
                cat = c(rep("ESAT6 QFT+", 3), rep("ESAT6 QFT-", 3), rep("Ag85B QFT+", 3), rep("Ag85B QFT-", 3)))


df$time <- factor(df$time, levels=c("baseline", "peak", "memory"))
levels(df$time)
df$cat = as.factor(df$cat)

ggplot(df, aes(y=value, x=time, group = cat)) +
  geom_line(aes(color=cat),linewidth=2)+
  geom_point(size=5, aes(color = cat))+
  theme_minimal()+
  xlab(NULL)+
  ylab(expression(paste("log(", mu[ij], ")")))+
  theme(text = element_text(size = 28)) +
  scale_color_manual(values=c("#FF9999","#99CCFF","#FF0000","blue"), name = NULL)


########################################
# Interaction Plot for CCR7+CD45RA-:
########################################
esat6_qft_pos = c(0,1.752-0.620 ,1.103-0.587)+1.114-1.048+1.174
esat6_qft_neg = c(0,1.752,1.103)-1.048
ag85b_qft_pos = c(0,1.752-0.620 ,1.103-0.587)+1.114
ag85b_qft_neg =c(0,1.752,1.103)
time = rep(c("baseline", "peak", "memory"),4)

df2 = data.frame(value = c(esat6_qft_pos, esat6_qft_neg, ag85b_qft_pos, ag85b_qft_neg),
                time = time,
                cat =  c(rep("ESAT6 QFT+", 3), rep("ESAT6 QFT-", 3), rep("Ag85B QFT+", 3), rep("Ag85B QFT-", 3)))


df2$time <- factor(df2$time, levels=c("baseline", "peak", "memory"))
levels(df2$time)
df2$cat = as.factor(df2$cat)

ggplot(df2, aes(y=value, x=time, group = cat)) +
  geom_line(aes(color=cat),linewidth=2)+
  geom_point(size=5, aes(color = cat))+
  theme_minimal()+
  xlab(NULL)+
  ylab(expression(paste("log(", mu[ij], ")")))+
  theme(text = element_text(size = 28)) +
  scale_color_manual(values=c("#FF9999","#99CCFF","#FF0000","blue"), name = NULL)




###############################################
# Interaction Plot for ESAT6 BOTH Markers:
###############################################

mm_qft_pos = c(0, 0.959+0.774-0.699 , 0.577+0.569-0.671)+1.182+2.367
mm_qft_neg = c(0, 0.959+0.774, 0.577+0.569)+1.182
mp_qft_pos = c(0, 0.959-0.699 , 0.577-0.671)+2.367
mp_qft_neg =c(0, 0.959, 0.577)
time = rep(c("baseline", "peak", "memory"),4)

df = data.frame(value = c(mm_qft_pos, mm_qft_neg, mp_qft_pos, mp_qft_neg),
                time = time,
                cat = c(rep("CCR7-CD45RA- QFT+", 3), rep("CCR7-CD45RA- QFT-", 3), rep("CCR7+CD45RA- QFT+", 3), rep("CCR7+CD45RA- QFT-", 3)))


df$time <- factor(df$time, levels=c("baseline", "peak", "memory"))
levels(df$time)
df$cat = as.factor(df$cat)

ggplot(df, aes(y=value, x=time, group = cat)) +
  geom_line(aes(color=cat),linewidth=2)+
  geom_point(size=5, aes(color = cat))+
  theme_minimal()+
  xlab(NULL)+
  ylab(expression(paste("log(", mu[ij], ")")))+
  theme(text = element_text(size = 28)) +
  scale_color_manual(values=c("#FF9999","#99CCFF","#FF0000","blue"), name = NULL)





###############################################
# Interaction Plot for Ag85B BOTH Markers:
###############################################

qft_pos = c(0, 2.064-0.783, 1.388-0.702) + 1.207
qft_neg = c(0, 2.064, 1.388) 
time = rep(c("baseline", "peak", "memory"),2)
df = data.frame(value = c(qft_pos, qft_neg),
                time = time,
                qft = c(rep("QFT+",3), rep("QFT-",3)))
df$qft = as.factor(df$qft)
df$time <- factor(df$time, levels=c("baseline", "peak", "memory"))
levels(df$time)

ggplot(df, aes(y=value, x=time, group = qft)) +
  geom_line(aes(color=qft),linewidth=2)+
  geom_point(size=5, aes(color = qft))+
  xlab(NULL)+
  ylab(expression(paste("log(", mu[ij], ")")))+
  theme_minimal()+
  theme(text = element_text(size = 28))+
  scale_color_manual(values=c("#FF0000", "blue"), name = NULL)



