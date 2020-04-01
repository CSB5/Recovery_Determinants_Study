## 07022020 CKR
#Figure 1D
library(vegan)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(stringr)
library(BisRNA)
library(ggpubr)

#select for species
SpecificSpecies<-function(country,suppl,species){
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

a=data.frame(t(CA_abun[,as.matrix(CA_R)]),status="R",subject=rownames(t(CA_abun[,as.matrix(CA_R)])), country=country)
a=melt(a,id=c("status","subject","country")) #convert recoverer data to right format
b=data.frame(t(CA_abun[,as.matrix(CA_NR)]),status="NR",subject=rownames(t(CA_abun[,as.matrix(CA_NR)])),country=country)
b=melt(b,id=c("status","subject","country")) #convert non-recoverer data to right format
final=rbind(a,b)
final=filter(final,variable%in%species) #filter for species to be plotted
return(final)
}


species=c("Bacteroides_uniformis","Alistipes_putredinis","Alistipes_shahii","Parabacteroides_distasonis",
          "Bacteroides_thetaiotaomicron","Bifidobacterium_adolescentis") #define species to be plotted
finalCA=SpecificSpecies("CA","../Data/Suppl1_CA_Final.txt",species)
finalSG=SpecificSpecies("SG","../Data/Suppl1_SG.txt",species)
finalSW=SpecificSpecies("SW","../Data/Suppl1_SW.txt",species)
finalEN=SpecificSpecies("EN","../Data/Suppl1_EN.txt",species)
final=rbind(finalCA,finalSG,finalSW,finalEN)

final$status2<-factor(final$status,levels=levels(final$status)[c(2,1)]) #changing order to plot

#Alistipes_putredinis
final1=filter(final,variable=="Alistipes_putredinis") %>% filter(country=="SG")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,30))

final1=filter(final,variable=="Alistipes_putredinis") %>% filter(country=="CA")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,40))

final1=filter(final,variable=="Alistipes_putredinis") %>% filter(country=="EN")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,4))

final1=filter(final,variable=="Alistipes_putredinis") %>% filter(country=="SW")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,6))

#Alistipes_shahii
final1=filter(final,variable=="Alistipes_shahii") %>% filter(country=="SG")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,5.5))

final1=filter(final,variable=="Alistipes_shahii") %>% filter(country=="CA")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,8))

final1=filter(final,variable=="Alistipes_shahii") %>% filter(country=="EN")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,1.2))

final1=filter(final,variable=="Alistipes_shahii") %>% filter(country=="SW")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,1))

#Bacteroides_uniformis
final1=filter(final,variable=="Bacteroides_uniformis") %>% filter(country=="SG")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,2.5))

final1=filter(final,variable=="Bacteroides_uniformis") %>% filter(country=="CA")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,3))

final1=filter(final,variable=="Bacteroides_uniformis") %>% filter(country=="EN")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,0.8))

final1=filter(final,variable=="Bacteroides_uniformis") %>% filter(country=="SW")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,7))

#Bacteroides_thetaiotaomicron
final1=filter(final,variable=="Bacteroides_thetaiotaomicron") %>% filter(country=="SG")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,1.2))

final1=filter(final,variable=="Bacteroides_thetaiotaomicron") %>% filter(country=="CA")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,0.5))

final1=filter(final,variable=="Bacteroides_thetaiotaomicron") %>% filter(country=="EN")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,3))

final1=filter(final,variable=="Bacteroides_thetaiotaomicron") %>% filter(country=="SW")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,9))

#Bifidobacterium_adolescentis
final1=filter(final,variable=="Bifidobacterium_adolescentis") %>% filter(country=="SG")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,0.12))

final1=filter(final,variable=="Bifidobacterium_adolescentis") %>% filter(country=="CA")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,1.2))

final1=filter(final,variable=="Bifidobacterium_adolescentis") %>% filter(country=="EN")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,0.15))

final1=filter(final,variable=="Bifidobacterium_adolescentis") %>% filter(country=="SW")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,0.11))

#Parabacteroides_distasonis
final1=filter(final,variable=="Parabacteroides_distasonis") %>% filter(country=="SG")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,0.1))

final1=filter(final,variable=="Parabacteroides_distasonis") %>% filter(country=="CA")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,0.5))

final1=filter(final,variable=="Parabacteroides_distasonis") %>% filter(country=="EN")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,2))

final1=filter(final,variable=="Parabacteroides_distasonis") %>% filter(country=="SW")
ggplot(data = final1, aes(x=status2,y=value,fill=status)) +
  geom_boxplot(width=2,outlier.shape=NA) +
  theme(axis.text.x = element_text(angle=40, hjust=1, size=10)) +
  scale_fill_manual(values=c("green", "red"))+
  theme_classic()+
  border()+
  coord_cartesian(ylim=c(0,1.25))
