## 13022020 CKR
#Figure 1B
library(vegan)
library(ggsci)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggsci)

Canada=read.table("../Data/Suppl1_CA_Final.txt",sep="\t",head=T,row.names=1) #read abundances
Canada=t(Canada)
label = read.table("../Data/TariniLabels_edit.txt",sep="\t",h=T) #read labels
label=filter(label,Country=="CA") %>% filter(str_detect(ID,"E90")|Status=="C") #filter for CA subjects, controls and post-treatment cases 
Canada=select(data.frame(Canada),as.character(label$ID)) #filter for abundances of subject of interest
#calculate distance and PCoA 
dist.mat=vegdist(t(Canada),method="euclidean")
cmds <- cmdscale(dist.mat, k=3, eig=TRUE)
eigen <- cmds$eig / sum(cmds$eig) * 100
data=data.frame(cmds$points,labels=rownames(data.frame(cmds$points)))
#annotate points on PCoA plot
C=as.character(filter(label,Status=="C")$ID)
R=as.character(filter(label,Status=="R")$ID)
N=as.character(filter(label,Status=="NR")$ID)
dat.merged=mutate(data,status=ifelse(labels%in%C,"C",ifelse(labels%in%N,"N","R")))

#plot 
ggplot(dat.merged%>%arrange(status), aes(x=-X1, y=X2,col=status), lwd=2) +
    geom_point(size=5,shape=18)+
    labs(x=paste0('PC1 (',round(eigen[1], 1),'%)'),
    y=paste0('PC2 (',round(eigen[2], 1),'%)')) +
    scale_color_manual(values=c("lightblue","red","green"))+
    theme_classic()+border()