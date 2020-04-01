## 17022020 CKR
#Figure 2D
library(vegan)
library(ggsci)
library(tidyr)
library(reshape2)
library(dplyr)
library(ggpubr)
library(stringr)
library(beanplot)

#Canada
Canada=read.table("../Data/Suppl1_CA_Final.txt",sep="\t",head=T,row.names=1)
Canada=t(Canada)
Canada_predur=select(data.frame(Canada),matches("E0|E7")) # pre and during timepoints
Canada_predur=Canada_predur[rowSums(Canada_predur)!=0,]

#species present in at least 5 samples
Canada_predur=Canada_predur[rowSums(Canada_predur!=0)>4,]

#read in rabs
Rabs=read.table("../Data/rabs_final.txt",head=T,sep="\t")
Rabs=paste0("s__",as.character(as.matrix(Rabs))) #fix species name to match
Canada_pd_rab=Canada_predur[rownames(Canada_predur)%in%Rabs,] #filter for rab aundances across samples
Canada_pd_rab=data.frame(Canada_pd_rab,species=rownames(Canada_pd_rab))
Canada_pd_rab=melt(Canada_pd_rab,id="species")#change to right format
Canada_pd_rab$variable=gsub("E[0-9]","E",Canada_pd_rab$variable)#change all samples of same subject to same label

Canada_pd_Nonrab=Canada_predur[!(rownames(Canada_predur)%in%Rabs),]#filter for non rab abundance
Canada_pd_Nonrab=data.frame(Canada_pd_Nonrab,species=rownames(Canada_pd_Nonrab))
Canada_pd_Nonrab=melt(Canada_pd_Nonrab,id="species")#change to right format
Canada_pd_Nonrab$variable=gsub("E[0-9]","E",Canada_pd_Nonrab$variable) #change all samples of same subject to same label

CanadaPTR=read.table("../Data/Suppl5_CA.txt",sep="\t",head=T,row.names=1)#read in PTR
CanadaPTR=CanadaPTR[1:(nrow(CanadaPTR)-1),] #remove last row
CanadaPTR[is.na(CanadaPTR)]=0 #change NA to 0
CanadaPTR=CanadaPTR[rowSums(CanadaPTR)!=0,]

#less than 50% of subjects do not have PTR for the species
CanadaPTR=CanadaPTR[rowSums(CanadaPTR==0)<(ncol(CanadaPTR)/2),]
CanadaPTR[CanadaPTR<1]=1 #set PTR < 1 to 1 for non-removed species
CanadaPTR=summarise_all(CanadaPTR,median)
CanadaPTR=select(CanadaPTR,contains("E90")) #post-treatment timepoint
CanadaPTR=melt(CanadaPTR)#change to right format
CanadaPTR$variable=gsub("E90","E",CanadaPTR$variable)

Canada_rabPTR=merge(Canada_pd_rab,CanadaPTR,by="variable")
C_rabPTR=group_by(Canada_rabPTR,species)%>%summarise(spearman=cor(value.x,value.y,method="spearman")) %>%
    mutate(status="Rabs")
Canada_NonrabPTR=merge(Canada_pd_Nonrab,CanadaPTR,by="variable")
C_NonrabPTR=group_by(Canada_NonrabPTR,species)%>%summarise(spearman=cor(value.x,value.y,method="spearman")) %>%
    mutate(status="NonRabs")
C=data.frame(rbind(C_rabPTR,C_NonrabPTR),country="Canada")

#Singapore
Singapore=read.table("../Data/Suppl1_SG.txt",sep="\t",head=T,row.names=1) #read in abundances
colnames(Singapore)=paste0("s__",colnames(Singapore))
Singapore=data.frame(t(Singapore),check.names = F)
Sing_predur=select(data.frame(Singapore,check.names = F),matches("INITIAL|CLOSING")) #select pre and during samples
Sing_predur=Sing_predur[rowSums(Sing_predur)!=0,]

#species present in at least 5 samples
Sing_predur=Sing_predur[rowSums(Sing_predur!=0)>4,]

Rabs=read.table("../Data/rabs_final.txt",head=T,sep="\t")
Rabs=paste0("s__",as.character(as.matrix(Rabs))) #fix species name to match
Sing_rab=Sing_predur[rownames(Sing_predur)%in%Rabs,]
Sing_rab=data.frame(Sing_rab,species=rownames(Sing_rab),check.names = F)
Sing_rab=melt(Sing_rab,id="species")%>%mutate(Variable=str_sub(variable,1,5))%>%select(-variable)

Sing_Nonrab=Sing_predur[!(rownames(Sing_predur)%in%Rabs),]
Sing_Nonrab=data.frame(Sing_Nonrab,species=rownames(Sing_Nonrab),check.names = F)
Sing_Nonrab=melt(Sing_Nonrab,id="species")%>%mutate(Variable=str_sub(variable,1,5))%>%select(-variable)

SingPTR=read.table("../Data/Suppl5_SG_Final.txt",sep="\t",head=T,row.names=1,check.names=F) #read in PTR
SingPTR=SingPTR[1:(nrow(SingPTR)-1),] #remove last row
SingPTR[is.na(SingPTR)]=0 #change NA to 0
SingPTR=SingPTR[rowSums(SingPTR)!=0,] #remove row with no abundances 

#less than 50% of subjects do not have PTR for the species
SingPTR=SingPTR[(rowSums(SingPTR!=0)>(ncol(SingPTR)/2)),]
SingPTR[SingPTR<1]=1
SingPTR=summarise_all(SingPTR,median)
SingPTR=select(SingPTR,contains("POST"))#select for post
SingPTR=melt(SingPTR)%>%mutate(Variable=str_sub(variable,1,5))%>%select(-variable)
SingPTR=group_by(SingPTR,Variable)%>%summarise(MeanPTR=mean(value)) #mean of community PTR of all post samples for same individual

Sing_rabPTR=merge(Sing_rab,SingPTR,by="Variable")
S_rabPTR=group_by(Sing_rabPTR,species)%>%summarise(spearman=cor(value,MeanPTR,method="spearman")) %>%
  mutate(status="Rabs")
Sing_NonrabPTR=merge(Sing_Nonrab,SingPTR,by="Variable")
S_NonrabPTR=group_by(Sing_NonrabPTR,species)%>%summarise(spearman=cor(value,MeanPTR,method="spearman")) %>%
  mutate(status="NonRabs")
S=data.frame(rbind(S_rabPTR,S_NonrabPTR),country="Singapore")

final=rbind(S,C)

#canadaBeanPlot
beanplot(spearman~status,data=filter(final,country=="Canada"),side="both",col=list("red","green"),border="white",ll=0
         ,beanlines="median",what = c(F, T, T, F),beanlinewd=1)

#test for Canada
wilcox.test(filter(final,status=="Rabs"&country=="Canada")$spearman,filter(final,status=="NonRabs"&country=="Canada")$spearman,alternative="g")
#test for Singapore
wilcox.test(filter(final,status=="Rabs"&country=="Singapore")$spearman,filter(final,status=="NonRabs"&country=="Singapore")$spearman, alternative="g")
#test for merged data
wilcox.test(filter(final,status=="Rabs")$spearman,filter(final,status=="NonRabs")$spearman,alternative="g")
