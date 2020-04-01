## 17022020 CKR
#Figure 2E
library(vegan)
library(ggsci)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggpubr)
library(stringr)
library(matrixStats)

#Canada
Caz_Canada=read.table("../Data/Suppl4_CA.txt",sep="\t",head=T,row.names=1)
Caz_Canada=select(Caz_Canada,contains("Number"))%>%
  mutate(label=rownames(Caz_Canada))
colnames(Caz_Canada)[1]="Number"
#Canada PTR
CanadaPTR=read.table("../Data/Suppl5_CA.txt",sep="\t",head=T,row.names=1)#read in PTR
CanadaPTR=CanadaPTR[1:(nrow(CanadaPTR)-1),] #remove last row
CanadaPTR[is.na(CanadaPTR)]=0 #change NA to 0
CanadaPTR=CanadaPTR[rowSums(CanadaPTR)!=0,]
#filter for common species
CanadaPTR=CanadaPTR[rowSums(CanadaPTR==0)<(ncol(CanadaPTR)/2),] #less than 50% of subjects do not have PTR for the species
CanadaPTR[CanadaPTR<1]=1 #set PTR < 1 to 1 for non-removed species
CanadaPTR=summarise_all(CanadaPTR,median)
CanadaPTR=melt(CanadaPTR)
Canada=merge(CanadaPTR,Caz_Canada,by.x="variable",by.y="label")

#plot  
ggplot(Canada, aes(x=Number,y=value)) +
  geom_point()  +
  geom_smooth(method = "lm")+
  theme_classic()+border()
cor.test(Canada$value,Canada$Number,method="spearman")

#Singapore
Caz_Singapore=read.table("../Data/Suppl4_SG_Final.txt",sep="\t",head=T,check.names=F,row.names=1)
Caz_Singapore=select(Caz_Singapore,contains("Number"))%>%
  mutate(label=rownames(Caz_Singapore))
colnames(Caz_Singapore)[1]="Number"

SingPTR=read.table("../Data/Suppl5_SG_Final.txt",sep="\t",head=T,row.names=1,check.names=F) #read in PTR
SingPTR=SingPTR[1:(nrow(SingPTR)-1),] #remove last row
SingPTR[is.na(SingPTR)]=0 #change NA to 0
SingPTR=SingPTR[rowSums(SingPTR)!=0,] #remove row with no abundances 

#less than 50% of subjects do not have PTR for the species
SingPTR=SingPTR[(rowSums(SingPTR!=0)>(ncol(SingPTR)/2)),]
SingPTR[SingPTR<1]=1
SingPTR=summarise_all(SingPTR,median)
SingPTR=melt(SingPTR)
Singapore=merge(SingPTR,Caz_Singapore,by.x="variable",by.y="label")

#plot  
ggplot(Singapore, aes(x=Number,y=value)) +
  geom_point()  +
  geom_smooth(method = "lm")+
  theme_classic()+border()
cor.test(Singapore$value,Singapore$Number,method="spearman")
