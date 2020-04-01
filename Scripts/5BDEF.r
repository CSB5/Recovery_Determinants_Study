####
# Tannistha Nandi
# Gut Manuscript
# 16th September 2019
####

library(dplyr)
library(Rmisc) 
library(ggplot2)
library(vegan)
library(gridExtra)
library(ggExtra)
library(ggpubr)
require(cowplot)
library(extrafont)

# function to compute median summary
summarySE_2 = function (data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                        conf.interval = 0.95, .drop = TRUE)
{
  length2 <- function(x, na.rm = FALSE) {
    if (na.rm)
      sum(!is.na(x))
    else length(x)
  }
  datac <- ddply(data, groupvars, .drop = .drop, .fun = function(xx, col, na.rm) {
    c(N = length2(xx[, col], na.rm = na.rm),
      mean = mean(xx[, col], na.rm = na.rm),
      median = median(xx[, col], na.rm = na.rm),
      sd = sd(xx[, col], na.rm = na.rm),
      mad=mad(xx[, col], na.rm = na.rm))
  }, measurevar, na.rm)
  datac <- rename(datac, c(mean = measurevar))
  datac$se <- datac$sd/sqrt(datac$N)
  ciMult <- qt(conf.interval/2 + 0.5, datac$N - 1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

####################################################################
# Figure 5B Total Biomass normalised by plant
# read data
data=read.csv("../Data/totalBiomass_plant_inp.txt", header=T, sep ="\t")
pd <- position_dodge(0.1)
tgc <- summarySE_2(data, measurevar="Total_Biomass_plant2", .drop= T,groupvars=c("Treatment", "Days"))
tgc=filter(tgc,Treatment != "BC" & Treatment != "BCBA") #filter off BC treated mice

p <- ggplot(tgc, aes(x=Days, y=median, colour=Treatment, linetype=Treatment)  ) +  scale_y_log10(breaks=c(10,100,1000)) +
  geom_errorbar(aes(ymin=median - mad, ymax=median + mad),  width=.5, position=pd) +
  geom_line(position=pd, size=1.5)+  
  geom_point(position=pd, size=4, shape=16) + 
  xlab("Days since start of experiment") +
  ylab("Bacterial Biomass") +    
  scale_linetype_manual("", values=c("twodash", "twodash", "twodash", "twodash"), guide="none")+
  scale_color_manual(values=c("darkkhaki", "cyan", "darkblue", "grey"),  labels=c("Ba", "Bt", "Bt + Ba",  "Vehicle"))+
  scale_x_continuous(limits=c(-1,23), breaks=c(0,3,6,7,10,13,16,19,22)) +       
  theme_bw() +
  theme(legend.justification="center",  legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6), text=element_text(family="Arial"),
        legend.text=element_text(family="Arial", size=15),legend.title=element_text(family="Arial", size=20), 
        axis.text.x = element_text(family="Arial", colour="grey20",size=15,angle=0,hjust=0.5,vjust=0.5,face="plain"),
        axis.text.y = element_text(family="Arial", colour="grey20",size=15,angle=0,hjust=1,vjust=0.5,face="plain"),  
        axis.title.x = element_text(family="Arial", colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(family="Arial", colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))              # Position legend in bottom right

p1=p+removeGrid(x = TRUE, y = TRUE)
p1

y=data[which(data$Treatment=="BA" & data$Days==10 | data$Treatment=="PBS" & data$Days==10 ),]
x=data[which(data$Treatment=="BTBA" & data$Days==10 & data$Cage <3 | data$Treatment=="BT" & data$Days==10 & data$Cage <3) ,]
wilcox.test(x$Total_Biomass_plant2, y$Total_Biomass_plant2, alternative = 'g')
y=data[which(data$Treatment=="BA" & data$Days==13 | data$Treatment=="PBS" & data$Days==13 ),]
x=data[which(data$Treatment=="BTBA" & data$Days==13 & data$Cage <3 | data$Treatment=="BT" & data$Days==13 & data$Cage <3) ,]
wilcox.test(x$Total_Biomass_plant2, y$Total_Biomass_plant2, alternative = 'g')
y=data[which(data$Treatment=="BA" & data$Days==16 | data$Treatment=="PBS" & data$Days==16 ),]
x=data[which(data$Treatment=="BTBA" & data$Days==16 & data$Cage <3 | data$Treatment=="BT" & data$Days==16 & data$Cage <3) ,]
wilcox.test(x$Total_Biomass_plant2, y$Total_Biomass_plant2, alternative = 'g')
y=data[which(data$Treatment=="BA" & data$Days==19 | data$Treatment=="PBS" & data$Days==19 ),]
x=data[which(data$Treatment=="BTBA" & data$Days==19 & data$Cage <3 | data$Treatment=="BT" & data$Days==19 & data$Cage <3) ,]
wilcox.test(x$Total_Biomass_plant2, y$Total_Biomass_plant2, alternative = 'g')
y=data[which(data$Treatment=="BA" & data$Days==22 | data$Treatment=="PBS" & data$Days==22 ),]
x=data[which(data$Treatment=="BTBA" & data$Days==22 & data$Cage <3 | data$Treatment=="BT" & data$Days==22 & data$Cage <3) ,]
wilcox.test(x$Total_Biomass_plant2, y$Total_Biomass_plant2, alternative = 'g')

####################################################################################
# Fig 5D
data=read.csv("../Data/CaZyme_inp.txt", header=T, sep ="\t")
pd <- position_dodge(0.1) 

tgc1 <- summarySE_2(data, measurevar="cell_wall_N", .drop= T,groupvars=c("Treatment", "Days"))
tgc1=filter(tgc1,Treatment != "BC" & Treatment != "BCBA") #filter off BC treated mice

p <- ggplot(tgc1, aes(x=Days, y=median, colour=Treatment, linetype=Treatment)  ) + scale_y_log10(breaks=c(20,30,50,60, 80, 100, 120))+ 
  geom_errorbar(aes(ymin=median - mad, ymax=median + mad),  width=.5, position=pd) +
  geom_line(position=pd, size=1.5)+  
  geom_point(position=pd, size=4, shape=16) + 
  xlab("Days since start of experiment") +
  ylab("Plant/Animal \n Cell Wall Degradation") +    
  scale_linetype_manual("", values=c("twodash", "twodash", "twodash", "twodash"), guide="none")+
  scale_color_manual(values=c("darkkhaki", "cyan", "darkblue", "grey"),  labels=c("BA", "BT", "BT + BA",  "Vehicle"))+
  ggtitle("") +
  scale_x_continuous(limits=c(-1,23), breaks=c(0,3,6,7,10,13,16,19,22, 25, 28)) +        
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position="right", 
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=0.5,vjust=0.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))              # Position legend in bottom right

p3=p+removeGrid(x = TRUE, y = TRUE) 
p3

#################################################################################
# Fig 5E

tgc2 <- summarySE_2(data, measurevar="mucin_N", .drop= T,groupvars=c("Treatment", "Days"))
tgc2=filter(tgc2,Treatment != "BC" & Treatment != "BCBA")#filter off BC treated mice

p <- ggplot(tgc2, aes(x=Days, y=median, colour=Treatment, linetype=Treatment)  ) + scale_y_log10(breaks=c(1,5,10,15,30,300))+
  geom_errorbar(aes(ymin=median - mad, ymax=median + mad),  width=.5, position=pd) +
  geom_line(position=pd, size=1.5)+  
  geom_point(position=pd, size=4, shape=16) +
  xlab("Days since start of experiment") +
  ylab("Mucin Degradation") +    
  scale_linetype_manual("", values=c("twodash", "twodash", "twodash", "twodash"), guide="none")+
  scale_color_manual(values=c("darkkhaki", "cyan", "darkblue", "grey"),  labels=c("BA", "BT", "BT + BA",  "Vehicle"))+
  ggtitle("") +
  scale_x_continuous(limits=c(-1,23), breaks=c(0,3,6,7,10,13,16,19,22, 25, 28)) +     
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position="right", 
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=0.5,vjust=0.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))              # Position legend in bottom right

p4=p+removeGrid(x = TRUE, y = TRUE) 
p4

####################################################################################################
# Fig 5F
tgc3 <- summarySE_2(data, measurevar="peptidoglycan_N", .drop= T,groupvars=c("Treatment", "Days"))
tgc3=filter(tgc3,Treatment != "BC" & Treatment != "BCBA") #filter off BC treated mice

p <- ggplot(tgc3, aes(x=Days, y=median, colour=Treatment, linetype=Treatment)  ) + scale_y_log10()+
  geom_errorbar(aes(ymin=median - mad, ymax=median + mad),  width=.5, position=pd) +
  geom_line(position=pd, size=1.5)+  
  geom_point(position=pd, size=4, shape=16) +
  xlab("Days since start of experiment") +
  ylab("Peptidoglycan \n Degradation") +    
  scale_linetype_manual("", values=c("twodash", "twodash", "twodash", "twodash"), guide="none")+
  scale_color_manual(values=c("darkkhaki", "cyan", "darkblue", "grey"),  labels=c("BA", "BT", "BT + BA",  "Vehicle"))+
  ggtitle("") +
  scale_x_continuous(limits=c(-1,23), breaks=c(0,3,6,7,10,13,16,19,22, 25, 28)) +        
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position="right", 
        axis.text.x = element_text(colour="grey20",size=15,angle=0,hjust=0.5,vjust=0.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=15,angle=0,hjust=1,vjust=0.5,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))              # Position legend in bottom right

p5=p+removeGrid(x = TRUE, y = TRUE) + scale_y_continuous(breaks=c(1,2,5,10,15)) 
p5

###############################################################################
# wilcoxon test for Fig 5D,E,F

y=data[which(data$Treatment=="BA" & data$Days==10 | data$Treatment=="PBS" & data$Days==10 ),]
x=data[which(data$Treatment=="BTBA" & data$Days==10 & data$Cage <3 | data$Treatment=="BT" & data$Days==10 & data$Cage <3) ,]
wilcox.test(x$cell_wall_N, y$cell_wall_N, alternative = 'g')
wilcox.test(x$mucin_N, y$mucin_N, alternative = 'g')
wilcox.test(x$peptidoglycan_N, y$peptidoglycan_N, alternative = 'g')
y=data[which(data$Treatment=="BA" & data$Days==13 | data$Treatment=="PBS" & data$Days==13 ),]
x=data[which(data$Treatment=="BTBA" & data$Days==13 & data$Cage <3 | data$Treatment=="BT" & data$Days==13 & data$Cage <3) ,]
wilcox.test(x$cell_wall_N, y$cell_wall_N, alternative = 'g')
wilcox.test(x$mucin_N, y$mucin_N, alternative = 'g')
wilcox.test(x$peptidoglycan_N, y$peptidoglycan_N, alternative = 'g')
y=data[which(data$Treatment=="BA" & data$Days==16 | data$Treatment=="PBS" & data$Days==16 ),]
x=data[which(data$Treatment=="BTBA" & data$Days==16 & data$Cage <3 | data$Treatment=="BT" & data$Days==16 & data$Cage <3) ,]
wilcox.test(x$cell_wall_N, y$cell_wall_N, alternative = 'g')
wilcox.test(x$mucin_N, y$mucin_N, alternative = 'g')
wilcox.test(x$peptidoglycan_N, y$peptidoglycan_N, alternative = 'g')
y=data[which(data$Treatment=="BA" & data$Days==19 | data$Treatment=="PBS" & data$Days==19 ),]
x=data[which(data$Treatment=="BTBA" & data$Days==19 & data$Cage <3 | data$Treatment=="BT" & data$Days==19 & data$Cage <3) ,]
wilcox.test(x$cell_wall_N, y$cell_wall_N, alternative = 'g')
wilcox.test(x$mucin_N, y$mucin_N, alternative = 'g')
wilcox.test(x$peptidoglycan_N, y$peptidoglycan_N, alternative = 'g')
y=data[which(data$Treatment=="BA" & data$Days==22 | data$Treatment=="PBS" & data$Days==22 ),]
x=data[which(data$Treatment=="BTBA" & data$Days==22 & data$Cage <3 | data$Treatment=="BT" & data$Days==22 & data$Cage <3) ,]
wilcox.test(x$cell_wall_N, y$cell_wall_N, alternative = 'g')
wilcox.test(x$mucin_N, y$mucin_N, alternative = 'g')
wilcox.test(x$peptidoglycan_N, y$peptidoglycan_N, alternative = 'g')

###############################################################################################3
