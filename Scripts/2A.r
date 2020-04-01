#Tarini 18032020
#Figure 2A
library(gplots)
library(dendextend)
library(RColorBrewer)
library(ggplot2)

PlantCellWall <- c("GH1","GH2","GH3","GH4","GH5","GH8","GH9","GH11","GH12","GH15","GH16","GH17","GH26","GH27","GH28","GH29","GH36","GH39","GH43","GH44","GH48","GH51","GH53","GH55","GH67","GH74","GH78","GH93","GH94","GH95","GH115","GH117","GH121","PL1","PL2","PL6","PL7","PL9","PL11","PL15","PL22")
Animal <- c("GH1","GH2","GH3","GH4","GH18","GH19","GH20","GH29","GH33","GH38","GH58","GH79","GH84","GH85","GH88","GH89","GH92","GH95","GH98","GH99","GH101","GH105","GH109","GH110","GH113","PL6","PL8","PL12","PL13","PL21")
Peptidoglycan <- c("GH23","GH24","GH25","GH73","GH102","GH103","GH104","GH108")
Starch_Glycogen <- c("GH13","GH15","GH57","GH77")
Sucrose_Fructans <- c("GH32","GH68","GH70","GH91")
Fungal <- c("GH5","GH8","GH16","GH18","GH19","GH20","GH55","GH64","GH71","GH81")
Dextran <- c("GH66","GH70","GH87")
HumanDigestive <- c("GH13","GH1","GH37","GH31","GH18","GH35","GH9")
Mucin<-c("GH101","GH129","GH2","GH20","GH29","GH33","GH42","GH84","GH89","GH95","GH98")
IndividualGroupsCumulated <- read.table("../Data/IndividualGroupsCumulated.txt",sep="\t",row.names=1,header=TRUE)
RABs <- read.table("../Data/rabs_final.txt",header=TRUE,sep="\t")
RABs_List <- as.character(RABs[,1])

IndividualGroupsCumulated$SpeciesType <- ifelse(rownames(IndividualGroupsCumulated) %in% RABs_List,"RAB","nonRAB")

#For Figure 3A: Comparing total counts of GH and PL families for RABs and non RABs
wilcox.test(IndividualGroupsCumulated$Total~IndividualGroupsCumulated$SpeciesType)

#plot
ggplot(IndividualGroupsCumulated,aes(Total,color=SpeciesType)) + geom_density(mapping=NULL,size=2) + labs(y="Fraction",x="Total CAZyme Families (GH and PL)") + theme_bw()

