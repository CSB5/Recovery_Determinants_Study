## 17022020 CKR
#Figure 2B
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
Caz_Canada=read.table("../Data/Suppl4_CA.txt",sep="\t",head=T,row.names=1)
Caz_Canada=select(Caz_Canada,contains("Number"))%>%
  mutate(label=rownames(Caz_Canada))
colnames(Caz_Canada)[1]="Number"

label = read.table("../Data/TariniLabels_edit.txt",sep="\t",h=T) #read labels
label=filter(label,Country=="CA")
data_C=merge(Caz_Canada,label,by.x="label",by.y="ID")

#Singapore
Caz_Singapore=read.table("../Data/Suppl4_SG_Final.txt",sep="\t",head=T,check.names=F,row.names=1)
Caz_Singapore=select(Caz_Singapore,contains("Number"))%>%
  mutate(label=rownames(Caz_Singapore))
colnames(Caz_Singapore)[1]="Number"

label = read.table("../Data/TariniLabels_edit.txt",sep="\t",h=T) #read labels
label=filter(label,Country=="SG")
data_S=merge(Caz_Singapore,label,by.x="label",by.y="ID")

final=rbind(data_C,data_S)%>%select(-label)
final$Country=as.character(final$Country)
final$Status=as.character(final$Status)

#beanplot
beanplot(Number~Status*Country,data=final,side="both",col=list("red","green"),border="white",ll=0
         ,beanlines="median",what = c(FALSE, TRUE, TRUE, FALSE), beanlinewd =0.5)
legend("topleft",fill = c("red", "green"), legend = c("Non-Recoverers", "Recoverers"), cex=0.65)

wilcox.test(filter(final,Status=="NR" & Country=="SG")$Number,filter(final,Status=="R" & Country=="SG")$Number,alternative="l")
wilcox.test(filter(final,Status=="NR" & Country=="CA")$Number,filter(final,Status=="R"& Country=="CA")$Number,alternative="l")
