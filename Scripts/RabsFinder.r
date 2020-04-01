## 07022020 CKR
#Rab identification and validation
library(vegan)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(stringr)
library(BisRNA)
#Stage 1
enriched<-function(country,suppl,side){
Label=read.table("../Data/TariniLabels_edit.txt", sep='\t', head = T) #read in labels
LabelCA=filter(Label,Country==country) #filter for country
CA_R=filter(LabelCA,Status=="R")%>%select(ID) #filter for recoverers
CA_NR=filter(LabelCA,Status=="NR")%>%select(ID) #filter for non-recoverers
CA_abun=read.table(suppl,sep='\t', head=T, row.names=1,check.names = F) #read in abundance
CA_abun=as.matrix(t(CA_abun)) #transpose
if(country=="SW"){
  rownames(CA_abun)=str_replace_all(rownames(CA_abun)," ","_")
} #change species name to right format
if(country=="CA"){
  rownames(CA_abun)=str_replace_all(rownames(CA_abun),"s__","")
} #change species name to right format
CA_abun=cbind(CA_abun[,as.matrix(CA_R)],CA_abun[,as.matrix(CA_NR)]) #abundance profile for recoverers and non-recoverers only
CA_abun=CA_abun[rowSums(CA_abun)>0, ] #remove species not present in recoverers and non-recoverers
a=CA_abun[,as.matrix(CA_R)]
b=CA_abun[,as.matrix(CA_NR)]
p=NULL
for(i in 1:nrow(a)){
  p=rbind(p,c(rownames(a)[i],wilcox.test(as.matrix(a[i,]),as.matrix(b[i,]),alternative=side)$p.value))
} #Wilcoxon test
colnames(p)=c("Species","pval")
return(p)
}

#get all one sided p-values
pCan=enriched("CA","../Data/Suppl1_CA_Final.txt","greater")
pSin=enriched("SG","../Data/Suppl1_SG.txt","greater")
pEN=enriched("EN","../Data/Suppl1_EN.txt","greater")
pSW=enriched("SW","../Data/Suppl1_SW.txt","greater")

pCanL=enriched("CA","../Data/Suppl1_CA_Final.txt","less")
pSinL=enriched("SG","../Data/Suppl1_SG.txt","less")
pENL=enriched("EN","../Data/Suppl1_EN.txt","less")
pSWL=enriched("SW","../Data/Suppl1_SW.txt","less")

#merge all pvalues and combine
mergeP<-function(pCan,pSin,pEN,pSW){
pAll=merge(pCan,pSin,by="Species")
colnames(pAll)=c("Species","pCA","pSG")
pAll=merge(pAll,pEN,by="Species")
colnames(pAll)=c("Species","pCA","pSG","pEN")
pAll=merge(pAll,pSW,by="Species")
colnames(pAll)=c("Species","pCA","pSG","pEN","pSW")
pAll$pCA=as.numeric(as.character(pAll$pCA))
pAll$pSG=as.numeric(as.character(pAll$pSG))
pAll$pEN=as.numeric(as.character(pAll$pEN))
pAll$pSW=as.numeric(as.character(pAll$pSW))
return(pAll)
}

pAll=mergeP(pCan,pSin,pEN,pSW) #p-value for greater
pAll_lower=mergeP(pCanL,pSinL,pENL,pSWL) #p-value for lower

#filter for BH corrected fisher combined pval<0.01
pFish=mutate(pAll,pfish=p.adjust(apply(select(pAll,-"Species"),1,fisher.method),method="BH"))%>%
  filter(pfish<0.01)
pFish_lower=mutate(pAll_lower,pfish=p.adjust(apply(select(pAll_lower,-"Species"),1,fisher.method),method="BH"))%>%
  filter(pfish<0.01) #none of the species satisfy the criteria of being significantly lower in non-recoverers
selected=as.matrix(pFish$Species)

#Stage 2
cohortfdr<-function(country,suppl,selected){
Label=read.table("../Data/TariniLabels_edit.txt", sep='\t', head = T) #read in labels
LabelCA=filter(Label,Country==country) #filter for country
CA_R=filter(LabelCA,Status=="R")%>%select(ID) #filter for recoverers
CA_NR=filter(LabelCA,Status=="NR")%>%select(ID) #filter for non-recoverers
CA_abun=read.table(suppl,sep='\t', head=T, row.names=1,check.names=F) #read in abundance
CA_abun=t(CA_abun) #transpose
if(country=="SW"){
  rownames(CA_abun)=str_replace_all(rownames(CA_abun)," ","_")
}#change species name to right format
if(country=="CA"){
  rownames(CA_abun)=str_replace_all(rownames(CA_abun),"s__","")
}#change species name to right format
CA_abun=CA_abun[selected,] #filter for species that passed stage 1
CA_abun=cbind(CA_abun[,as.matrix(CA_R)],CA_abun[,as.matrix(CA_NR)])#abundance profile for recoverers and non-recoverers only
CA_abun=CA_abun[rowSums(CA_abun)>0, ] #remove species not present in the target group
a=CA_abun[,as.matrix(CA_R)]
b=CA_abun[,as.matrix(CA_NR)]
p=NULL
for(i in 1:nrow(a)){
  p=rbind(p,c(rownames(a)[i],wilcox.test(as.matrix(a[i,]),as.matrix(b[i,]),alternative="greater")$p.value))
}#Wilcoxon test

padj=p.adjust(p[,2],method="BH",n=nrow(p))#BH correction
p=cbind(p[,1],padj)
colnames(p)=c("Species","padj")
return(p)
}

CohpCan=cohortfdr("CA","../Data/Suppl1_CA_Final.txt",selected)
CohpSG=cohortfdr("SG","../Data/Suppl1_SG.txt",selected)
CohpEN=cohortfdr("EN","../Data/Suppl1_EN.txt",selected)
CohpSW=cohortfdr("SW","../Data/Suppl1_SW.txt",selected)

#merge all p-values from different cohort
CohpAll=merge(CohpCan,CohpSG,by="Species")
colnames(CohpAll)=c("Species","pCA","pSG")
CohpAll=merge(CohpAll,CohpEN,by="Species")
colnames(CohpAll)=c("Species","pCA","pSG","pEN")
CohpAll=merge(CohpAll,CohpSW,by="Species")
colnames(CohpAll)=c("Species","pCA","pSG","pEN","pSW")
CohpAll$pCA=as.numeric(as.character(CohpAll$pCA))
CohpAll$pSG=as.numeric(as.character(CohpAll$pSG))
CohpAll$pEN=as.numeric(as.character(CohpAll$pEN))
CohpAll$pSW=as.numeric(as.character(CohpAll$pSW))

#filter for species with at least 2 cohorts fdr p-value < 0.05
filtered=CohpAll[rowSums(select(CohpAll,-"Species")<0.05)>=2,]

#include NUH cohort
Data=read.table("../Data/NUH_StoolSamples_MetaPhlAn2.txt", sep='\t', head = T, row.names=1) #read abundances
Data=filter(Data,Group=="Case") %>%
  select(PatientID,TimePoint,contains("s__"),-contains("t__")) #select for cases who took antibiotics and species profile
colnames(Data)=gsub(".*s__","s__",colnames(Data)) #changing species name to right format

##identify recoverers and non-recoverers
DataAnn=filter(Data,TimePoint=="POST")%>%
  select(-TimePoint)
DataAnn=data.frame(DataAnn,row.names=1) #set patient ID to rownames
Diversity=data.frame(Simpson=diversity((DataAnn),index='simpson'),ID=rownames(DataAnn))
Diversity=data.frame(Diversity,row.names=2)
Threshold = Diversity %>%
  summarise(IQR=IQR(Simpson),Median=median(Simpson)) %>%
  summarise(Upper=(0.1*IQR)+Median,Lower=Median-(0.1*IQR))#Define threshold
Diversity=mutate(Diversity,Lib=rownames(Diversity))
Recoverer=filter(Diversity,Simpson>as.numeric(Threshold["Upper"]))%>%
  select(-Simpson) #Define recoverer
NonRec=filter(Diversity,Simpson<as.numeric(Threshold["Lower"])) %>%
  select(-Simpson) #Define non-recoverer

#filter for abundance for R and NR
Metarec=filter(Data,PatientID%in%as.matrix(Recoverer))%>%
  select(-TimePoint,-PatientID)
MetaNon=filter(Data,PatientID%in%as.matrix(NonRec))%>%
  select(-TimePoint,-PatientID)
RecData=t(Metarec) #transpose
NonData=t(MetaNon) #transpose
p=NULL
for(i in 1:nrow(RecData)){
  p=rbind(p,c(rownames(RecData)[i],wilcox.test(as.matrix(RecData[i,]),as.matrix(NonData[i,]),alternative="greater")$p.value))
} #wilcoxon test
p[,1]=str_replace(p[,1],"s__","")#change Species name to right format
colnames(p)=c("Species","pvalue")
p=filter(data.frame(p),Species%in%selected)

#if NUH cohort included in stage 2
c_fdr=mutate(p,padj=p.adjust(as.numeric(as.character(p$pvalue)),method="BH"))%>%
  filter(padj<0.05)
c_fdr[,1]%in%filtered[,1] #add 2 more rab

#independent validation by NUH cohort on rabs
p=filter(p,as.numeric(as.character(p$pvalue))<0.1) #no corrrection
p[,1]%in%filtered[,1] #12 out of 21 validated

