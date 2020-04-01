## 17022020 CKR
#Figure 2C
library(vegan)
library(ggsci)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggpubr)
library(stringr)
library(matrixStats)
library(beanplot)

#Canada
CanadaPTR=read.table("../Data/Suppl5_CA.txt",sep="\t",head=T,row.names=1)#read in PTR
CanadaPTR=CanadaPTR[1:(nrow(CanadaPTR)-1),] #remove last row
CanadaPTR[is.na(CanadaPTR)]=0 #change NA to 0
CanadaPTR=CanadaPTR[rowSums(CanadaPTR)!=0,]

#filter for common species
CanadaPTR=CanadaPTR[rowSums(CanadaPTR==0)<(ncol(CanadaPTR)/2),]
CanadaPTR[CanadaPTR<1]=1 #set PTR < 1 to 1 for non-removed species
CanadaPTR=summarise_all(CanadaPTR,median)
CanadaPTR=melt(CanadaPTR)

label = read.table("../Data/TariniLabels_edit.txt",sep="\t",h=T) #read labels
label=filter(label,Country=="CA")
data_C=merge(CanadaPTR,label,by.x="variable",by.y="ID")

SingPTR=read.table("../Data/Suppl5_SG_Final.txt",sep="\t",head=T,row.names=1,check.names=F) #read in PTR
SingPTR=SingPTR[1:(nrow(SingPTR)-1),] #remove last row
SingPTR[is.na(SingPTR)]=0 #change NA to 0
SingPTR=SingPTR[rowSums(SingPTR)!=0,] #remove row with no abundances 

#less than 50% of subjects do not have PTR for the species
SingPTR=SingPTR[(rowSums(SingPTR!=0)>(ncol(SingPTR)/2)),]
SingPTR[SingPTR<1]=1
SingPTR=summarise_all(SingPTR,median)
SingPTR=melt(SingPTR)

label = read.table("../Data/TariniLabels_edit.txt",sep="\t",h=T) #read labels
label=filter(label,Country=="SG")
data_S=merge(SingPTR,label,by.x="variable",by.y="ID")

final=rbind(data_C,data_S)%>%select(-variable)
final$Country=as.character(final$Country)
final$Status=as.character(final$Status)

#beanplot
beanplot(value~Status*Country,data=final,side="both",col=list("red","green"),border="white",ll=0, ylim=c(0.9,1.4)
         ,beanlines="median",what = c(FALSE, TRUE, TRUE, FALSE), beanlinewd =0.5)
legend("topleft",fill = c("red", "green"), legend = c("Non-Recoverers", "Recoverers"), cex=0.65)

wilcox.test(filter(final,Status=="NR" & Country=="SG")$value,filter(final,Status=="R" & Country=="SG")$value,alternative="l")
wilcox.test(filter(final,Status=="NR" & Country=="CA")$value,filter(final,Status=="R"& Country=="CA")$value,alternative="l")