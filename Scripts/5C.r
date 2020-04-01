####
# Tannistha Nandi
# 19th July2018
###Estimated diversity in the mouse gut and found differences between the various groups. 
####
#Figure 5C

library(dplyr)
library(Rcpp)
library(tidyr)
library(ggplot2)
library(ggsignif)
library(vegan)
library(gridExtra)
library(plyr)
library(ggpubr)

data=read.csv("../Data/Gavage_adjust_diversityV2.txt", header=T, stringsAsFactors = FALSE, sep="\t", row.names = 1)
metadata=data[,1:3]
data_s=data[,4:887]
data_s=t(data_s)
data_percent= apply(data_s, 2,function(x) x/sum(x))
data_s=t(data_percent)

species = merge(data_s, metadata, by=0)
rownames(species) = species$Row.names
species = species[,-1]

species$Diversity = diversity(species[,1:884], index = "simpson")

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

pd <- position_dodge(0.1)
tgc <- summarySE_2(species, measurevar="Diversity", .drop= T,groupvars=c("Treatment", "Days"))
tgc=filter(tgc,Treatment!="BC"&Treatment!="BCBA") #Filter off BC treated mice

p <- ggplot(tgc, aes(x=Days, y=median, colour=Treatment, linetype=Treatment)  ) +  
  geom_errorbar(aes(ymin=(median - mad), ymax=(median + mad)),  width=.5, position=pd) +
  geom_line(position=pd,size=1.5)+  
  geom_point(position=pd,size=4,shape=16) + 
  xlab("Days since start of experiment") +
  ylab("Simpson Diversity") +    
  scale_linetype_manual("", values=c("twodash", "twodash", "twodash", "twodash"), guide="none")+
  scale_color_manual(values=c("darkkhaki", "cyan", "darkblue", "grey"))+
  scale_x_continuous(limits=c(-1,23), breaks=c(0,3,6,7,10, 13, 16,19, 22)) +       
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.justification="center",  legend.box.background = element_rect(), legend.box.margin = margin(6, 6, 6, 6), text=element_text(family="Arial"),
        legend.text=element_text(family="Arial", size=15),legend.title=element_text(family="Arial", size=20), 
        axis.text.x = element_text(family="Arial", colour="grey20",size=15,angle=0,hjust=0.5,vjust=0.5,face="plain"),
        axis.text.y = element_text(family="Arial", colour="grey20",size=15,angle=0,hjust=1,vjust=0.5,face="plain"),  
        axis.title.x = element_text(family="Arial", colour="grey20",size=15,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(family="Arial", colour="grey20",size=15,angle=90,hjust=.5,vjust=.5,face="bold"))              # Position legend in bottom right
p

######### Wilcoxon Test
x=species[which(species$Treatment=="BTBA" & species$Days==16 & species$Cage <3),]
y=species[which(species$Treatment=="BA" & species$Days==16 | species$Treatment=="BT" & species$Days==16 & species$Cage<3 | species$Treatment=="PBS" & species$Days==16) ,]
wilcox.test(x$Diversity, y$Diversity,alternative="g")

x=species[which(species$Treatment=="BTBA" & species$Days==19 & species$Cage <3),]
y=species[which(species$Treatment=="BA" & species$Days==19 | species$Treatment=="BT" & species$Days==19 & species$Cage<3 | species$Treatment=="PBS" & species$Days==19) ,]
wilcox.test(x$Diversity, y$Diversity,alternative="g")

x=species[which(species$Treatment=="BTBA" & species$Days==22 & species$Cage <3),]
y=species[which(species$Treatment=="BA" & species$Days==22 | species$Treatment=="BT" & species$Days==22 & species$Cage<3 | species$Treatment=="PBS" & species$Days==22) ,]
wilcox.test(x$Diversity, y$Diversity,alternative="g")
