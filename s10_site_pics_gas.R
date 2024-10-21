#### CRC  barplot for sites
### 2024-9-6
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

## 1 Preparation
##rm(list = ls())

#library("survival")
#library("survminer")
#library(gridExtra)
#library("grid")
#library(reshape2)	  
#library("RColorBrewer")
##library("plyr")

library("ggplot2")
library(reshape2)	
library("RColorBrewer")
library(plyr)

#install.packages("patchwork")
library("patchwork")
source("functions.R")
##install.packages('ggbeeswarm')
library("ggbeeswarm")

library(ggthemes)


###################### 1.2 Inputs
info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

vsv_info<-read.table("01.cleanData/SV_info/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

all_basic<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv",header=T)

load("01.cleanData/SV_all/dsgv.RData")
load("01.cleanData/SV_all/vsgv.RData")

vsgv<-vsgv_sub
dsgv<-dsgv_sub


#######################################################################################
### clinical data and dsgv/vsgv
all_basic$id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"[Eubacterium] rectale DSM 17629:2268_2269",
"[Eubacterium] rectale DSM 17629:1110_1113",
"[Eubacterium] rectale DSM 17629:1122_1124;1124_1125",
"Collinsella sp. 4_8_47FAA:1991_1996",
"Lachnospiraceae bacterium 3_1_46FAA:2740_2741;3118_3119",
"Bilophila wadsworthia ATCC 49260:7_9;2906_2908",
"Ruminococcus sp. SR1/5:1250_1252")]

########################################################################################
##  
dsgv_pic$match_id<-row.names(dsgv_pic)
pic_dsv<-merge(all_basic,dsgv_pic,by.x="match_id",by.y="match_id")


#########################################################################
#####  [Eubacterium] rectale DSM 17629:1122_1124;1124_1125

pic_dsv$temp<-pic_dsv$"[Eubacterium] rectale DSM 17629:1122_1124;1124_1125"
#lg_res_model <- glm(response_code~cohort+temp+age+read_count,family=binomial(logit), data=pic_dsv)
#lg_res <- summary(lg_res_model)

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(CRC)==F,drop=T)
pic_sub$disease<-pic_sub$CRC
table(pic_sub$study_name)

LiuNN_2022_pic<-subset(pic_sub,study_name=="LiuNN_2022",drop=T)
ThomasAM_2018_pic<-subset(pic_sub,study_name=="ThomasAM_2018",drop=T)
VogtmannE_2016_pic<-subset(pic_sub,study_name=="VogtmannE_2016",drop=T)
WirbelJ_2018_pic<-subset(pic_sub,study_name=="WirbelJ_2018",drop=T)
YachidaS_2019_pic<-subset(pic_sub,study_name=="YachidaS_2019",drop=T)
YuJ_2015_pic<-subset(pic_sub,study_name=="YuJ_2015",drop=T)
ZellerG_2014_pic<-subset(pic_sub,study_name=="ZellerG_2014",drop=T)

LiuNN_2022_ratio<-dsv_disease(LiuNN_2022_pic,"LiuNN_2022")
ThomasAM_2018_ratio<-dsv_disease(ThomasAM_2018_pic,"ThomasAM_2018")
VogtmannE_2016_ratio<-dsv_disease(VogtmannE_2016_pic,"VogtmannE_2016")
WirbelJ_2018_ratio<-dsv_disease(WirbelJ_2018_pic,"WirbelJ_2018")
YachidaS_2019_ratio<-dsv_disease(YachidaS_2019_pic,"YachidaS_2019")
YuJ_2015_ratio<-dsv_disease(YuJ_2015_pic,"YuJ_2015")
ZellerG_2014_ratio<-dsv_disease(ZellerG_2014_pic,"ZellerG_2014")


table(LiuNN_2022_pic$CRC)
p1<-ggplot(LiuNN_2022_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("LiuNN\n2022")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = 'Ratio')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("health","disease"),
   labels=c("Health (n=45)","CRC (n=34)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1,color="black"),text=element_text(family="sans"),legend.position="none",
   plot.title=element_text(size=11,face="plain")
   )

table(ThomasAM_2018_pic$CRC)
p2<-ggplot(ThomasAM_2018_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("ThomasAM\n2018")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("health","disease"),
   labels=c("Health (n=25)","CRC (n=30)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain")
   )

table(VogtmannE_2016_pic$CRC)
p3<-ggplot(VogtmannE_2016_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("VogtmannE\n2016")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("health","disease"),
   labels=c("Health (n=42)","CRC (n=41)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain")
   )

table(WirbelJ_2018_pic$CRC)
p4<-ggplot(WirbelJ_2018_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("WirbelJ\n2018")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("health","disease"),
   labels=c("Health (n=51)","CRC (n=10)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain")
   )

table(YachidaS_2019_pic$CRC)
p5<-ggplot(YachidaS_2019_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("YachidaS\n2019")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("health","disease"),
   labels=c("Health (n=122)","CRC (n=133)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain")
   )

table(YuJ_2015_pic$CRC)
p6<-ggplot(YuJ_2015_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("YuJ\n2015")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("health","disease"),
   labels=c("Health (n=24)","CRC (n=20)"),,expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain"),plot.margin = unit(c(0, 0, 0, 0), "pt")
   )


table(ZellerG_2014_pic$CRC)
p7<-ggplot(ZellerG_2014_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("ZellerG\n2014")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("health","disease"),
   labels=c("Health (n=45)","CRC (n=26)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="right",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain"),legend.title=element_blank()
   )

##p_whole<-p1+p2+p3+p4+p5+p6+p7+plot_layout(ncol=7) + plot_annotation(title = "dsv Eubacterium rectale DSM 17629:1122 ~ 1124;1124 ~1125 kbp\nMeta OR=0.46;Meta p=1.25e-4")

p_whole<-p1+p2+p3+p4+p5+p6+p7+plot_layout(ncol=7) 

tiff(file = "pics/sites/F5_dsv_CRC_E.rectale1124_1125_barplot.tiff", width =2500, height =1400, res =300) 
p_whole
dev.off()



##################################################################################
### vsgv picture 

all_basic$id<-row.names(all_basic)

vsgv_pic<-vsgv[,c(
"Lachnospiraceae bacterium 3_1_46FAA:2765_2768",
"Parabacteroides distasonis ATCC 8503:4735_4737"
)]


vsgv_pic$id<-row.names(vsgv_pic)
### mel 
pic_vsv<-merge(all_basic,vsgv_pic,by.x="id",by.y="id")

pic_vsv$temp<-pic_vsv$"Lachnospiraceae bacterium 3_1_46FAA:2765_2768"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(CRC)==F,drop=T)

pic_sub$CRC<-as.factor(pic_sub$CRC)
flevels <- levels(pic_sub$CRC)

LiuNN_2022_pic<-subset(pic_sub,study_name=="LiuNN_2022",drop=T)
ThomasAM_2018_pic<-subset(pic_sub,study_name=="ThomasAM_2018",drop=T)
VogtmannE_2016_pic<-subset(pic_sub,study_name=="VogtmannE_2016",drop=T)
WirbelJ_2018_pic<-subset(pic_sub,study_name=="WirbelJ_2018",drop=T)
YachidaS_2019_pic<-subset(pic_sub,study_name=="YachidaS_2019",drop=T)
YuJ_2015_pic<-subset(pic_sub,study_name=="YuJ_2015",drop=T)
ZellerG_2014_pic<-subset(pic_sub,study_name=="ZellerG_2014",drop=T)


table(LiuNN_2022_pic$CRC)

p1<-ggplot(data=LiuNN_2022_pic,aes(x=CRC,y=temp,color=CRC))+
  geom_quasirandom(aes(color=CRC),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("LiuNN\n2022")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#2983BB","#ECD452") )+
  scale_x_discrete(limits=flevels,labels=(c("health(n=39)","CRC(n=38)")))+
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1))+
  geom_hline(yintercept=c(-3,-2,-1,0,1,2,3),linetype="dashed")+
  theme_cowplot() +
  ylab("Normalized region coverage")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(ThomasAM_2018_pic$CRC)
p2<-ggplot(data=ThomasAM_2018_pic,aes(x=CRC,y=temp,color=CRC))+
  geom_quasirandom(aes(color=CRC),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("ThomasAM\n2018")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#2983BB","#ECD452") )+
  scale_x_discrete(limits=flevels,labels=(c("health(n=9)","CRC(n=16)")))+
  geom_hline(yintercept=c(-3,-2,-1,0,1,2,3),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(VogtmannE_2016_pic$CRC)
p3<-ggplot(data=VogtmannE_2016_pic,aes(x=CRC,y=temp,color=CRC))+
  geom_quasirandom(aes(color=CRC),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("VogtmannE\n2016")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#2983BB","#ECD452") )+
  scale_x_discrete(limits=flevels,labels=(c("health(n=21)","CRC(n=29)")))+
  geom_hline(yintercept=c(-3,-2,-1,0,1,2,3),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(WirbelJ_2018_pic$CRC)
p4<-ggplot(data=WirbelJ_2018_pic,aes(x=CRC,y=temp,color=CRC))+
  geom_quasirandom(aes(color=CRC),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("WirbelJ\n2018")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#2983BB","#ECD452") )+
  scale_x_discrete(limits=flevels,labels=(c("health(n=8)","CRC(n=5)")))+
  geom_hline(yintercept=c(-3,-2,-1,0,1,2,3),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(YachidaS_2019_pic$CRC)

p5<-ggplot(data=YachidaS_2019_pic,aes(x=CRC,y=temp,color=CRC))+
  geom_quasirandom(aes(color=CRC),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("YachidaS\n2019")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#2983BB","#ECD452") )+ggtitle("YachidaS\n2019")+
  scale_x_discrete(limits=flevels,labels=(c("health(n=74)","CRC(n=113)")))+
  geom_hline(yintercept=c(-3,-2,-1,0,1,2,3),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(YuJ_2015_pic$CRC)

p6<-ggplot(data=YuJ_2015_pic,aes(x=CRC,y=temp,color=CRC))+
  geom_quasirandom(aes(color=CRC),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#2983BB","#ECD452") )+ggtitle("YuJ\n2015")+
  scale_x_discrete(limits=flevels,labels=(c("health(n=11)","CRC(n=20)")))+
  geom_hline(yintercept=c(-3,-2,-1,0,1,2,3),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(ZellerG_2014_pic$CRC)
p7<-ggplot(data=ZellerG_2014_pic,aes(x=CRC,y=temp,color=CRC))+
  geom_quasirandom(aes(color=CRC),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#2983BB","#ECD452") )+ggtitle("ZellerG\n2014")+
  scale_x_discrete(limits=flevels,labels=(c("health(n=32)","CRC(n=30)")))+
  geom_hline(yintercept=c(-3,-2,-1,0,1,2,3),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=14,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))



#p_whole<-p1+p2+p3+p4+p5+p6+p7+plot_layout(ncol=7) + plot_annotation(title = "vSV Lachnospiraceae bacterium 3_1_46FAA:2765 ~ 2768 kbp\nMeta OR=0.57;Meta p=1.27e-6")

p_whole<-p1+p2+p3+p4+p5+p6+p7+plot_layout(ncol=7) 

tiff(file = "pics/sites/F5_vsv_CRC_Lachnospiraceaebacterium_barplot.tiff", width =2400, height =1500, res =300) 
p_whole
dev.off()




###################################################################################
## IBD


###################### 1.2 Inputs
info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

vsv_info<-read.table("01.cleanData/SV_info/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

all_basic<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv",header=T)

load("01.cleanData/SV_all/dsgv.RData")
load("01.cleanData/SV_all/vsgv.RData")

vsgv<-vsgv_sub
dsgv<-dsgv_sub


#######################################################################################
### clinical data and dsgv/vsgv
all_basic$id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"[Eubacterium] rectale DSM 17629:1122_1124;1124_1125",
"Ruminococcus sp. SR1/5:1250_1252")]

########################################################################################
##  
dsgv_pic$match_id<-row.names(dsgv_pic)
pic_dsv<-merge(all_basic,dsgv_pic,by.x="match_id",by.y="match_id")

########################################################################################
##  
pic_dsv$temp<-pic_dsv$"[Eubacterium] rectale DSM 17629:1122_1124;1124_1125"

#lg_res_model <- glm(response_code~cohort+temp+age+read_count,family=binomial(logit), data=pic_dsv)
#lg_res <- summary(lg_res_model)

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(IBD)==F,drop=T)
pic_sub$disease<-pic_sub$IBD
table(pic_sub$study_name)

FranzosaEA_2018_pic<-subset(pic_sub,study_name=="FranzosaEA_2018",drop=T)
HallAB_2017_pic<-subset(pic_sub,study_name=="HallAB_2017",drop=T)
HMP_2019_ibdmdb_pic<-subset(pic_sub,study_name=="HMP_2019_ibdmdb",drop=T)
#LiJ_2014_pic<-subset(pic_sub,study_name=="LiJ_2014",drop=T)

FranzosaEA_2018_ratio<-dsv_disease(FranzosaEA_2018_pic,"FranzosaEA_2018")
HallAB_2017_ratio<-dsv_disease(HallAB_2017_pic,"HallAB_2017")
HMP_2019_ibdmdb_ratio<-dsv_disease(HMP_2019_ibdmdb_pic,"HMP_2019_ibdmdb")


table(FranzosaEA_2018_pic$IBD)
p1<-ggplot(FranzosaEA_2018_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("FranzosaEA\n2018")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = 'Ratio')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("health","disease"),
   labels=c("Health (n=48)","IBD (n=75)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1,color="black"),text=element_text(family="sans"),legend.position="none",
   plot.title=element_text(size=12,face="plain")
   )

table(HallAB_2017_pic$IBD)
p2<-ggplot(HallAB_2017_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("HallAB\n2017")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("health","disease"),
   labels=c("Health (n=20)","IBD (n=104)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=12,face="plain")
   )

table(HMP_2019_ibdmdb_pic$IBD)
p3<-ggplot(HMP_2019_ibdmdb_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("HMP_2019_ibdmdb\n")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("health","disease"),
   labels=c("Health (n=16)","IBD (n=42)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1),text=element_text(family="sans"),legend.position="right",
   axis.text.y = element_blank(),plot.title=element_text(size=12,face="plain"),legend.title=element_blank()
   )


#p_whole<-p1+p2+p3+plot_layout(ncol=3) + plot_annotation(title = "dSV Eubacterium rectale DSM 17629:1122 ~ 1124;1124 ~ 1125 kbp\nMeta OR=0.24;Meta p=2.18e-5")

p_whole<-p1+p2+p3+plot_layout(ncol=3) 

tiff(file = "pics/sites/F5_dsv_IBD_Eubacteriumrectale_1124_barplot.tiff", width =1600, height =1400, res =300) 
p_whole
dev.off()



##################################################################################
### vsgv picture (different datasets) 

all_basic$id<-row.names(all_basic)
vsgv_pic<-vsgv[,c(
"Faecalibacterium cf. prausnitzii KLE1255:1372_1379",
"Faecalibacterium cf. prausnitzii KLE1255:866_872;872_874",
"Faecalibacterium cf. prausnitzii KLE1255:1525_1528",
"Faecalibacterium cf. prausnitzii KLE1255:2399_2400"
)]


vsgv_pic$id<-row.names(vsgv_pic)
pic_vsv<-merge(all_basic,vsgv_pic,by.x="id",by.y="id")

pic_vsv$temp<-pic_vsv$"Faecalibacterium cf. prausnitzii KLE1255:1372_1379"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(IBD)==F,drop=T)

pic_sub$IBD<-as.factor(pic_sub$IBD)
flevels <- levels(pic_sub$IBD)

table(pic_sub$study_name)

FranzosaEA_2018_pic<-subset(pic_sub,study_name=="FranzosaEA_2018",drop=T)
HallAB_2017_pic<-subset(pic_sub,study_name=="HallAB_2017",drop=T)
HMP_2019_ibdmdb_pic<-subset(pic_sub,study_name=="HMP_2019_ibdmdb",drop=T)


table(FranzosaEA_2018_pic$IBD)

p1<-ggplot(data=FranzosaEA_2018_pic,aes(x=IBD,y=temp,color=IBD))+
  geom_quasirandom(aes(color=IBD),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("FranzosaEA\n2018")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#207F4C","#FEBA07") )+
  scale_x_discrete(limits=flevels,labels=(c("health(n=45)","IBD(n=69)")))+
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1))+
  geom_hline(yintercept=c(-3,-2,-1,0,1,2,3),linetype="dashed")+
  theme_cowplot() +
  ylab("Normalized region coverage")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(HallAB_2017_pic$IBD)
p2<-ggplot(data=HallAB_2017_pic,aes(x=IBD,y=temp,color=IBD))+
  geom_quasirandom(aes(color=IBD),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("HallAB\n2017")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#207F4C","#FEBA07") )+
  scale_x_discrete(limits=flevels,labels=(c("health(n=18)","IBD(n=86)")))+
  geom_hline(yintercept=c(-3,-2,-1,0,1,2,3),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(HMP_2019_ibdmdb_pic$IBD)
p3<-ggplot(data=HMP_2019_ibdmdb_pic,aes(x=IBD,y=temp,color=IBD))+
  geom_quasirandom(aes(color=IBD),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("HMP_2019\nibdmdb")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#207F4C","#FEBA07") )+
  scale_x_discrete(limits=flevels,labels=(c("health(n=14)","IBD(n=59)")))+
  geom_hline(yintercept=c(-3,-2,-1,0,1,2,3),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))

#p_whole<-p1+p2+p3+plot_layout(ncol=3) + plot_annotation(title = "Faecalibacterium cf. prausnitzii KLE1255:1372 ~ 1379 kbp\nMeta OR=0.47;Meta p=2.2e-6")

p_whole<-p1+p2+p3+plot_layout(ncol=3) 

tiff(file = "pics/sites/F5_vsv_IBD_prausnitzii_barplot.tiff", width =1200, height =1400, res =300) 
p_whole
dev.off()









