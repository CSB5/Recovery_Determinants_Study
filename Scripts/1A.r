#Tarini 10032020
#Figure 1A
library(tidyverse)
library(vegan)

#Getting SG Labels and Median Diversity for the three stages
SG_Species <- read.table("../Data/Suppl1_SG.txt",sep="\t",row.names=1,header=TRUE)
SG_Sample_Details <- data.frame(t(rbind(as.data.frame(strsplit(rownames(SG_Species),"_"))[1,],as.data.frame(strsplit(rownames(SG_Species),"_"))[4,])),row.names=rownames(SG_Species))
colnames(SG_Sample_Details) <- c("SampleID","Stage")
SG_Sample_Details$Stage_Individual <- paste0(SG_Sample_Details$SampleID,"_",SG_Sample_Details$Stage)
SG_MedianDiversity_Stages <- as.data.frame(aggregate(as.data.frame(diversity(SG_Species,index="simpson")),by=list(SG_Sample_Details$Stage_Individual),FUN=median)[,-1])
rownames(SG_MedianDiversity_Stages) <- aggregate(as.data.frame(diversity(SG_Species)),by=list(SG_Sample_Details$Stage_Individual),FUN=median)[,1]
colnames(SG_MedianDiversity_Stages) <- "MedianDiversity"
Labels <- read.table("../Data/TariniLabels_edit.txt",sep="\t",row.names=1,header=TRUE)
SG_Labels <- Labels[Labels$Country == "SG",]
temp_sample_tags <- cbind(SG_Sample_Details[rownames(SG_Labels),],SG_Labels$Status)[,c(2:4)] %>% distinct
rownames(temp_sample_tags) <- temp_sample_tags[,2]
SG_MedianDiversity_Stages <- as.data.frame((cbind(temp_sample_tags,SG_MedianDiversity_Stages[rownames(temp_sample_tags),])))
colnames(SG_MedianDiversity_Stages) <- c("Stage","Stage_Individuals","Status","MedianDiversity")
SG_MedianDiversity_Stages$Stage <- ifelse(SG_MedianDiversity_Stages$Stage == "INITIAL","PRE",ifelse(SG_MedianDiversity_Stages$Stage == "CLOSING","DURING","POST"))

#Getting CA Labels and Median Diversity for the three stages
CA_Species <- read.table("../Data/Suppl1_CA_Final.txt",sep="\t",row.names=1,header=TRUE)
colnames(CA_Species) <- sub("s__","",colnames(CA_Species))
CA_Labels <- Labels[Labels$Country == "CA",]
CA_MedianDiversity_Stages <- as.data.frame(cbind(ifelse(rownames(CA_Labels) %in% grep("E0",rownames(CA_Labels),value=TRUE),"PRE",ifelse(rownames(CA_Labels) %in% grep("E7",rownames(CA_Labels),value=TRUE),"DURING","POST")),rownames(CA_Labels),as.character(CA_Labels$Status),as.numeric(diversity(CA_Species[rownames(CA_Labels),],index="simpson"))))
colnames(CA_MedianDiversity_Stages) <- c("Stage","Stage_Individuals","Status","MedianDiversity")
CA_MedianDiversity_Stages <- CA_MedianDiversity_Stages[CA_MedianDiversity_Stages$Status != "C",]
CA_MedianDiversity_Stages$MedianDiversity <- as.numeric(as.vector(CA_MedianDiversity_Stages$MedianDiversity))

#Rank Median Diversities within Cohorts
#Median Simpson Diversities were range normalized within Cohorts to remove Cohort-specific variations in Simpson Diversities
range_scale <- function(x)
{
        y <- (x-min(x))/(max(x)-min(x));
        return(y);
}
CA_MedianDiversity_Stages$MedianDiversity <- range_scale(CA_MedianDiversity_Stages$MedianDiversity)
SG_MedianDiversity_Stages$MedianDiversity <- range_scale(SG_MedianDiversity_Stages$MedianDiversity)

#Getting Names
SG_MedianDiversity_Stages$Stage_Individuals <- t(as.data.frame(strsplit(as.character(SG_MedianDiversity_Stages$Stage_Individuals),"_"))[1,])[,1]
CA_MedianDiversity_Stages$Stage_Individuals <- t(as.data.frame(strsplit(as.character(CA_MedianDiversity_Stages$Stage_Individuals),"E"))[1,])[,1]

#Combine and retain only Individuals with all the time_points
R_MedianDiversity_Stages <- as.data.frame(rbind(CA_MedianDiversity_Stages[CA_MedianDiversity_Stages$Status=="R",],SG_MedianDiversity_Stages[SG_MedianDiversity_Stages$Status=="R",]))
R_1 <- R_MedianDiversity_Stages
R_remove_rows <- c("218RZ_INITIAL","218RZ_POST","808ZJ_POST","808ZJ_INITIAL")
R_MedianDiversity_Stages <- R_MedianDiversity_Stages[setdiff(rownames(R_MedianDiversity_Stages),R_remove_rows),]
R_MedianDiversity_Stages$Stage <- factor(R_MedianDiversity_Stages$Stage,levels=c("PRE","DURING","POST"))

NR_MedianDiversity_Stages <- as.data.frame(rbind(CA_MedianDiversity_Stages[CA_MedianDiversity_Stages$Status=="NR",],SG_MedianDiversity_Stages[SG_MedianDiversity_Stages$Status=="NR",]))
NR_remove_rows <- c("412AD_INITIAL","412AD_POST","664SB_POST","664SB_CLOSING","987AB_INITIAL","987AB_POST")
NR_1 <- NR_MedianDiversity_Stages
NR_MedianDiversity_Stages <- NR_MedianDiversity_Stages[setdiff(rownames(NR_MedianDiversity_Stages),NR_remove_rows),]
NR_MedianDiversity_Stages$Stage <- factor(NR_MedianDiversity_Stages$Stage,levels=c("PRE","DURING","POST"))

R_summary <- data.frame(Stage=unique(factor(R_MedianDiversity_Stages$Stage,levels=c("PRE","POST","DURING"))), n=tapply(R_MedianDiversity_Stages$MedianDiversity, factor(R_MedianDiversity_Stages$Stage,levels=c("PRE","DURING","POST")), length), mean=tapply(R_MedianDiversity_Stages$MedianDiversity, factor(R_MedianDiversity_Stages$Stage,levels=c("PRE","DURING","POST")), median))
R_summary$sd <- tapply(R_MedianDiversity_Stages$MedianDiversity, factor(R_MedianDiversity_Stages$Stage,levels=c("PRE","DURING","POST")), sd)
R_summary$sem <- R_summary$sd/sqrt(R_summary$n-1)
R_summary$lower_25 <- R_summary$mean + qt((1-0.25)/2, df=R_summary$n-1)*R_summary$sem
R_summary$higher_25 <- R_summary$mean + qt((1+0.25)/2, df=R_summary$n-1)*R_summary$sem
R_summary$lower_50 <- R_summary$mean + qt((1-0.50)/2, df=R_summary$n-1)*R_summary$sem
R_summary$higher_50 <- R_summary$mean + qt((1+0.50)/2, df=R_summary$n-1)*R_summary$sem
R_summary$lower_75 <- R_summary$mean + qt((1-0.75)/2, df=R_summary$n-1)*R_summary$sem
R_summary$higher_75 <- R_summary$mean + qt((1+0.75)/2, df=R_summary$n-1)*R_summary$sem
R_summary$lower_80 <- R_summary$mean + qt((1-0.80)/2, df=R_summary$n-1)*R_summary$sem
R_summary$higher_80 <- R_summary$mean + qt((1+0.80)/2, df=R_summary$n-1)*R_summary$sem
rownames(R_summary) <- c("PRE","DURING","POST")
R_summary$Stage <- c(1,2,3)

NR_summary <- data.frame(Stage=unique(factor(NR_MedianDiversity_Stages$Stage,levels=c("PRE","POST","DURING"))), n=tapply(NR_MedianDiversity_Stages$MedianDiversity, factor(NR_MedianDiversity_Stages$Stage,levels=c("PRE","DURING","POST")), length), mean=tapply(NR_MedianDiversity_Stages$MedianDiversity, factor(NR_MedianDiversity_Stages$Stage,levels=c("PRE","DURING","POST")), median))
NR_summary$sd <- tapply(NR_MedianDiversity_Stages$MedianDiversity, factor(NR_MedianDiversity_Stages$Stage,levels=c("PRE","DURING","POST")), sd)
NR_summary$sem <- NR_summary$sd/sqrt(NR_summary$n-1)
NR_summary$lower_25 <- NR_summary$mean + qt((1-0.25)/2, df=NR_summary$n-1)*NR_summary$sem
NR_summary$higher_25 <- NR_summary$mean + qt((1+0.25)/2, df=NR_summary$n-1)*NR_summary$sem
NR_summary$lower_50 <- NR_summary$mean + qt((1-0.50)/2, df=NR_summary$n-1)*NR_summary$sem
NR_summary$higher_50 <- NR_summary$mean + qt((1+0.50)/2, df=NR_summary$n-1)*NR_summary$sem
NR_summary$lower_75 <- NR_summary$mean + qt((1-0.75)/2, df=NR_summary$n-1)*NR_summary$sem
NR_summary$higher_75 <- NR_summary$mean + qt((1+0.75)/2, df=NR_summary$n-1)*NR_summary$sem
NR_summary$lower_80 <- NR_summary$mean + qt((1-0.80)/2, df=NR_summary$n-1)*NR_summary$sem
NR_summary$higher_80 <- NR_summary$mean + qt((1+0.80)/2, df=NR_summary$n-1)*NR_summary$sem
rownames(NR_summary) <- c("PRE","DURING","POST")
NR_summary$Stage <- c(1,2,3)

#Plot Ribbon R
g1 <- ggplot(R_summary) + geom_line(aes(x=Stage,y=mean),stat="smooth",method="auto",color="blue4",size=2) + geom_line(aes(x=Stage,y=higher_25),stat="smooth",method="auto") + geom_line(aes(x=Stage,y=lower_25),stat="smooth",method="auto")+ geom_line(aes(x=Stage,y=higher_50),stat="smooth",method="auto") + geom_line(aes(x=Stage,y=lower_50),stat="smooth",method="auto")+ geom_line(aes(x=Stage,y=higher_75),stat="smooth",method="auto") + geom_line(aes(x=Stage,y=lower_75),stat="smooth",method="auto")
gg1 <- ggplot_build(g1)
g1_ribbon <- data.frame(x = gg1$data[[1]]$x, ymax_25 = gg1$data[[2]]$y, ymin_25 = gg1$data[[3]]$y, ymax_50 = gg1$data[[4]]$y, ymin_50 = gg1$data[[5]]$y, ymax_75 = gg1$data[[6]]$y, ymin_75 = gg1$data[[7]]$y) 
g1 +  geom_ribbon(data = g1_ribbon, aes(x = x, ymin = ymin_75, ymax = ymax_75), fill = "yellow", alpha = 0.4) +  geom_ribbon(data = g1_ribbon, aes(x = x, ymin = ymin_75, ymax = ymax_75), fill = "darkgoldenrod1", alpha = 0.4) + geom_ribbon(data = g1_ribbon, aes(x = x, ymin = ymin_50, ymax = ymax_50), fill = "orange", alpha = 0.4) + geom_ribbon(data = g1_ribbon, aes(x = x, ymin = ymin_25, ymax = ymax_25), fill = "firebrick4", alpha = 0.4) + theme_bw() + theme(axis.text=element_text(size=20)) +ylim(0.6,1)

g2 <- ggplot(NR_summary) + geom_line(aes(x=Stage,y=mean),stat="smooth",method="auto",color="blue4",size=2) + geom_line(aes(x=Stage,y=higher_25),stat="smooth",method="auto") + geom_line(aes(x=Stage,y=lower_25),stat="smooth",method="auto")+ geom_line(aes(x=Stage,y=higher_50),stat="smooth",method="auto") + geom_line(aes(x=Stage,y=lower_50),stat="smooth",method="auto")+ geom_line(aes(x=Stage,y=higher_75),stat="smooth",method="auto") + geom_line(aes(x=Stage,y=lower_75),stat="smooth",method="auto")
gg2 <- ggplot_build(g2)
g2_ribbon <- data.frame(x = gg2$data[[1]]$x, ymax_25 = gg2$data[[2]]$y, ymin_25 = gg2$data[[3]]$y, ymax_50 = gg2$data[[4]]$y, ymin_50 = gg2$data[[5]]$y, ymax_75 = gg2$data[[6]]$y, ymin_75 = gg2$data[[7]]$y) 
g2 +  geom_ribbon(data = g2_ribbon, aes(x = x, ymin = ymin_75, ymax = ymax_75), fill = "yellow", alpha = 0.4) +  geom_ribbon(data = g2_ribbon, aes(x = x, ymin = ymin_75, ymax = ymax_75), fill = "darkgoldenrod1", alpha = 0.4) + geom_ribbon(data = g2_ribbon, aes(x = x, ymin = ymin_50, ymax = ymax_50), fill = "orange", alpha = 0.4) + geom_ribbon(data = g2_ribbon, aes(x = x, ymin = ymin_25, ymax = ymax_25), fill = "firebrick4", alpha = 0.4) + theme_bw() + theme(axis.text=element_text(size=20)) + ylim(0.4,0.9)




