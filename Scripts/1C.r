## 14022020 CKR
#Figure 1C
library(ggplot2)
library(vegan)
library(ggsci)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggpubr)
library(stringr)
library(matrixStats)

Canada=read.table("../Data/Suppl1_CA_Final.txt",sep="\t",head=T,row.names=1)#read in abundances
Canada=t(Canada)#transpose
label = read.table("../Data/TariniLabels_edit.txt",sep="\t",h=T) #read labels
label=filter(label,Country=="CA") %>% filter(str_detect(ID,"E90")|Status=="C") #filter for CA subjects, controls and post-treatment cases 
Canada=select(data.frame(Canada),as.character(label$ID)) #filter for abundances of subject of interest

#calculate distance
dist=as.matrix(vegdist(t(Canada),method="bray"))
dist=data.frame(dist)
C=as.character(filter(label,Status=="C")$ID)
R=as.character(filter(label,Status=="R")$ID)
N=as.character(filter(label,Status=="NR")$ID)
distC=select(dist,C)

a=(as.matrix(distC[N,]))
distN=data.frame(t(as.matrix(distC[N,])))%>%summarise_all(median) #median of non-recoverer distance to controls
distR=data.frame(t(as.matrix(distC[R,])))%>%summarise_all(median) #median of recoverer distance to controls
distN=data.frame(gather(distN),label="N")
distR=data.frame(gather(distR),label="R")
Dist=data.frame(rbind(distR,distN))

#plot  
Dist$label2<-factor(Dist$label,levels=levels(Dist$label)[c(2,1)]) #change order
ggplot(data=Dist, aes(x=label2,y=value,fill=label,grp=label)) +
  geom_boxplot(width=0.5)  +
  theme_classic() +
  scale_fill_manual(values=c("green","red")) + 
  border()
wilcox.test(filter(Dist,label=="R")$value,filter(Dist,label=="N")$value,alternative="l")
