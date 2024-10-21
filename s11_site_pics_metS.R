#### Met  barplot for sites
### 2024-9-28
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
### HDL

## clinical data and dsgv/vsgv
all_basic$id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"Escherichia coli:4263_4266",
"Bifidobacterium pseudocatenulatum DSM 20438 = JCM 1200 = LMG 10505:1518_1520 and 2 segments")]

########################################################################################
##  
dsgv_pic$match_id<-row.names(dsgv_pic)
pic_dsv<-merge(all_basic,dsgv_pic,by.x="match_id",by.y="match_id")

#########################################################################
#####  Escherichia coli:4263_4266

pic_dsv$temp<-pic_dsv$"Escherichia coli:4263_4266"
pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(hdl)==F,drop=T)

table(pic_sub$study_name)

pic_sub$temp<-as.factor(pic_sub$temp)
flevels <- levels(pic_sub$temp)

JieZ_2017_pic<-subset(pic_sub,study_name=="JieZ_2017",drop=T)
#KarlssonFH_2013_pic<-subset(pic_sub,study_name=="KarlssonFH_2013",drop=T)
QinJ_2012_pic<-subset(pic_sub,study_name=="QinJ_2012",drop=T)
YuJ_2015_pic<-subset(pic_sub,study_name=="YuJ_2015",drop=T)


table(JieZ_2017_pic$temp)

p1<-ggplot(data=JieZ_2017_pic,aes(x=as.factor(temp),y=hdl,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("JieZ_2017")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#7876B1FF","#207F4C") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=17)","Non-delection(n=130)")))+
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5))+
  geom_hline(yintercept=c(0.5,1,1.5,2,2.5),linetype="dashed")+
  theme_cowplot() +
  ylab("HDL")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(QinJ_2012_pic$temp)
p2<-ggplot(data=QinJ_2012_pic,aes(x=as.factor(temp),y=hdl,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("QinJ_2012")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#7876B1FF","#207F4C") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=16)","Non-delection(n=72)")))+
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5))+
  geom_hline(yintercept=c(0.5,1,1.5,2,2.5),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(YuJ_2015_pic$temp)
p3<-ggplot(data=YuJ_2015_pic,aes(x=as.factor(temp),y=hdl,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("YuJ_2015")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#7876B1FF","#207F4C") )+
  scale_x_discrete(limits=flevels,labels=(c("Delection(n=16)","Non-delection(n=39)")))+
  scale_y_continuous(limits = c(0, 2.5), breaks = seq(0, 2.5, 0.5))+
  geom_hline(yintercept=c(0.5,1,1.5,2,2.5),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))

#p_whole<-p1+p2+p3+plot_layout(ncol=3) + plot_annotation(title = "dSV Escherichia coli:4263 ~ 4266 kbp\nMeta beta = -0.23;Meta p=9.41e-6")

p_whole<-p1+p2+p3+plot_layout(ncol=3) 

tiff(file = "pics/sites/F5_dSV_hdl_Escherichia_colibarplot.tiff", width =1400, height =1300, res =300) 
p_whole
dev.off()




################################################################################
## cholesterol 

## clinical data and dsgv/vsgv
all_basic$id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"Prevotella copri DSM 18205:699_701",
"Roseburia intestinalis XB6B4:182_184")]

########################################################################################
##  
dsgv_pic$match_id<-row.names(dsgv_pic)
pic_dsv<-merge(all_basic,dsgv_pic,by.x="match_id",by.y="match_id")

#########################################################################
#####  Prevotella copri DSM 18205:699_701

pic_dsv$temp<-pic_dsv$"Prevotella copri DSM 18205:699_701"
pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(cholesterol)==F,drop=T)

table(pic_sub$study_name)

pic_sub$temp<-as.factor(pic_sub$temp)
flevels <- levels(pic_sub$temp)

JieZ_2017_pic<-subset(pic_sub,study_name=="JieZ_2017",drop=T)
QinJ_2012_pic<-subset(pic_sub,study_name=="QinJ_2012",drop=T)


table(JieZ_2017_pic$temp)

p1<-ggplot(data=JieZ_2017_pic,aes(x=as.factor(temp),y=cholesterol,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("JieZ_2017")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#7876B1FF","#207F4C") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=25)","Non-delection(n=24)")))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2))+
  geom_hline(yintercept=c(2,4,6),linetype="dashed")+
  theme_cowplot() +
  ylab("Total cholesterol")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(QinJ_2012_pic$temp)
p2<-ggplot(data=QinJ_2012_pic,aes(x=as.factor(temp),y=cholesterol,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("QinJ_2012")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#7876B1FF","#207F4C") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=10)","Non-delection(n=30)")))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2))+
  geom_hline(yintercept=c(2,4,6),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))

#p_whole<-p1+p2+plot_layout(ncol=2) + plot_annotation(title = "dSV Prevotella copri DSM 18205:699~701 kbp\nMeta beta=0.73;Meta p=6.04e-4")

p_whole<-p1+p2+plot_layout(ncol=2) 

tiff(file = "pics/sites/F5_dSV_cholesterol_barplot.tiff", width =1000, height =1500, res =300) 
p_whole
dev.off()



################################################################################
## HDL  

# Bacteroides fragilis NCTC 9343:2487_2488 and 2 segments

## clinical data and dsgv/vsgv
all_basic$id<-row.names(all_basic)

vsgv_pic<-vsgv[,c(
"Bacteroides fragilis NCTC 9343:2487_2488 and 2 segments",
"Prevotella copri DSM 18205:197_198")]

########################################################################################
##  
vsgv_pic$match_id<-row.names(vsgv_pic)
pic_vsv<-merge(all_basic,vsgv_pic,by.x="match_id",by.y="match_id")

#########################################################################
#####  Bacteroides fragilis NCTC 9343:2487_2488 and 2 segments

pic_vsv$temp<-pic_vsv$"Bacteroides fragilis NCTC 9343:2487_2488 and 2 segments"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(hdl)==F & study_name!="KarlssonFH_2013",drop=T)

table(pic_sub$study_name)

# Save the scatter plot in a variable
p <- ggplot(pic_sub, aes(x = temp, y = hdl,color=study_name,group=study_name)) +
  geom_point()+
  geom_smooth(method = lm,aes(color=study_name,fill=study_name))+
  scale_color_manual(values = c("#BC3523","#2983BB","#ECD452")) +
  scale_fill_manual(values = c("#BC3523","#2983BB","#ECD452")) +
  theme_classic() +
  theme(legend.position="bottom") + 
  theme(text = element_text(size=16)) + 
  #ylim(0.0,1.3)+
  xlab("Normalized region coverage")+
  ylab("HDL")+
  #ggtitle("Bacteroides fragilis NCTC 9343:2487_2488 and 2 segments")+
  theme(panel.grid = element_blank(), 
  panel.background = element_rect(color = 'white', fill = 'transparent'), 
  strip.text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(1,1,1,1), "cm"),
  axis.text.x = element_text(angle = 0, hjust = 1))


tiff(file = "pics/sites/F5_vSV_hdl_barplot_Bacteroides_fragilis.tiff", width =1500, height =1300, res =300) 
p
dev.off()

aes(color=study_name)


##################################################################################
### LDL  vsgv picture 

all_basic$id<-row.names(all_basic)

vsgv_pic<-vsgv[,c(
"Bifidobacterium pseudocatenulatum DSM 20438 = JCM 1200 = LMG 10505:1806_1807;1807_1808",
"Parabacteroides distasonis ATCC 8503:4735_4737"
)]

vsgv_pic$id<-row.names(vsgv_pic)
pic_vsv<-merge(all_basic,vsgv_pic,by.x="id",by.y="id")

#################################################################################
##Bifidobacterium pseudocatenulatum DSM 20438 = JCM 1200 = LMG 10505:1806_1807;1807_1808

pic_vsv$temp<-pic_vsv$"Bifidobacterium pseudocatenulatum DSM 20438 = JCM 1200 = LMG 10505:1806_1807;1807_1808"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(hdl)==F,drop=T)

# Save the scatter plot in a variable
p <- ggplot(pic_sub, aes(x = temp, y = hdl,color=study_name)) +
  geom_point(aes(color=study_name))+
  scale_color_manual(values = c("#BC3523","#2983BB","#ECD452","#EA517F")) +
  scale_fill_manual(values = c("#BC3523","#2983BB","#ECD452","#EA517F")) +
  geom_smooth(method = lm,aes(color=study_name,fill=study_name))+
  theme_classic() +
  theme(legend.position.inside=c(0,65)) + 
  theme(text = element_text(size=16)) + 
  #ylim(0.0,1.3)+
  xlab("Normalized region coverage")+
  ylab("HDL")+
  ggtitle("Bifidobacterium pseudocatenulatum DSM 20438 = JCM 1200 = LMG 10505\n1806_1807;1807_1808")+
  theme(panel.grid = element_blank(), 
  panel.background = element_rect(color = 'black', fill = 'transparent'), 
  strip.text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(1,1,1,1), "cm"),
  axis.text.x = element_text(angle = 0, hjust = 1))

tiff(file = "pics/sites/Fig6_HDL_Bif_vsv_2166_2173.tiff", width =1800, height =1800, res =300) 
p
dev.off()


###################################################################################
## LDL

##################################################################################
### vsgv picture 

all_basic$id<-row.names(all_basic)

vsgv_pic<-vsgv[,c(
"Prevotella copri DSM 18205:197_198",
"Parabacteroides distasonis ATCC 8503:4735_4737"
)]

vsgv_pic$id<-row.names(vsgv_pic)
pic_vsv<-merge(all_basic,vsgv_pic,by.x="id",by.y="id")

#################################################################################
## Prevotella copri DSM 18205:197_198

pic_vsv$temp<-pic_vsv$"Prevotella copri DSM 18205:197_198"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(ldl)==F,drop=T)

# Save the scatter plot in a variable
p <- ggplot(pic_sub, aes(x = temp, y = ldl,color=study_name)) +
  geom_point(aes(color=study_name))+
  scale_color_manual(values = c("#ECD452","#207F4C","#EA8958","#2983BB","#BC3523")) +
  scale_fill_manual(values = c("#ECD452","#207F4C","#EA8958","#2983BB","#BC3523")) +
  geom_smooth(method = lm,aes(color=study_name,fill=study_name))+
  theme_classic() +
  guides(col = guide_legend(nrow =2, byrow = TRUE))+
  theme(legend.position="bottom") + 
  theme(text = element_text(size=16)) + 
  #ylim(0.0,1.3)+
  xlab("Normalized region coverage")+
  ylab("LDL")+
  #ggtitle("Prevotella copri DSM 18205:197_198")+
  theme(panel.grid = element_blank(), 
  panel.background = element_rect(color = 'white', fill = 'transparent'), 
  strip.text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(1,1,1,1), "cm"),
  axis.text.x = element_text(angle = 0, hjust = 1))

tiff(file = "pics/sites/Fig6_LDL_P_copri_vsv_197_198.tiff", width =1800, height =1800, res =300) 
p
dev.off()




################################################################################
## triglycerides

# Ruminococcus gnavus ATCC 29149:786_793;793_800

## clinical data and dsgv/vsgv
all_basic$id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"Ruminococcus gnavus ATCC 29149:786_793;793_800",
"Roseburia intestinalis XB6B4:182_184")]

########################################################################################
##  
dsgv_pic$match_id<-row.names(dsgv_pic)
pic_dsv<-merge(all_basic,dsgv_pic,by.x="match_id",by.y="match_id")

#########################################################################
#####  Ruminococcus gnavus ATCC 29149:786_793;793_800

pic_dsv$temp<-pic_dsv$"Ruminococcus gnavus ATCC 29149:786_793;793_800"
pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(hdl)==F,drop=T)

table(pic_sub$study_name)

pic_sub$temp<-as.factor(pic_sub$temp)
flevels <- levels(pic_sub$temp)

JieZ_2017_pic<-subset(pic_sub,study_name=="JieZ_2017",drop=T)
QinJ_2012_pic<-subset(pic_sub,study_name=="QinJ_2012",drop=T)
YuJ_2015_pic<-subset(pic_sub,study_name=="YuJ_2015",drop=T)


table(JieZ_2017_pic$temp)

p1<-ggplot(data=JieZ_2017_pic,aes(x=as.factor(temp),y=triglycerides,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("JieZ_2017")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#7876B1FF","#207F4C") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=71)","Non-delection(n=55)")))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2))+
  geom_hline(yintercept=c(2,4,6),linetype="dashed")+
  theme_cowplot() +
  ylab("Triglycerides")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(QinJ_2012_pic$temp)
p2<-ggplot(data=QinJ_2012_pic,aes(x=as.factor(temp),y=triglycerides,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("QinJ_2012")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#7876B1FF","#207F4C") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=19)","Non-delection(n=26)")))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2))+
  geom_hline(yintercept=c(2,4,6),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(YuJ_2015_pic$temp)
p3<-ggplot(data=YuJ_2015_pic,aes(x=as.factor(temp),y=triglycerides,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("YuJ_2015")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#7876B1FF","#207F4C") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=22)","Non-delection(n=17)")))+
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2))+
  geom_hline(yintercept=c(2,4,6),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


#p_whole<-p1+p2+p3+plot_layout(ncol=3) + plot_annotation(title = "dSV Ruminococcus gnavus ATCC 29149:786~793;793~800 kbp\nMeta beta=-0.13;Meta p=2.87e-4")

p_whole<-p1+p2+p3+plot_layout(ncol=3) 

tiff(file = "pics/sites/F5_dSV_triglycerides_Ruminococcus_gnavus_barplot.tiff", width =1300, height =1300, res =300) 
p_whole
dev.off()


#######################################################################################
##### cirrhosis
### clinical data and dsgv/vsgv
all_basic$id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"Bacteroides xylanisolvens XB1A:337_339;339_340",
"Bacteroides uniformis ATCC 8492:2444_2456")]

########################################################################################
##  
dsgv_pic$match_id<-row.names(dsgv_pic)
pic_dsv<-merge(all_basic,dsgv_pic,by.x="match_id",by.y="match_id")

#########################################################################
#####  Bacteroides xylanisolvens XB1A:337_339;339_340

pic_dsv$temp<-pic_dsv$"Bacteroides xylanisolvens XB1A:337_339;339_340"

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(cirrhosis)==F,drop=T)
pic_sub$disease<-pic_sub$cirrhosis
QinN_2014_ratio<-dsv_disease(pic_sub,"QinN_2014")

table(pic_sub$cirrhosis)

p1<-ggplot(QinN_2014_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("QinN_2014")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","No delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = 'Ratio')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("health","disease"),
   labels=c("Health (n=60)","Cirrhosis (n=69)"))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1,color="black"),text=element_text(family="sans"),legend.position="right",
   plot.title=element_text(size=11,face="plain"),legend.title=element_blank()
   )

#p_whole<-p1+plot_layout(ncol=1) + plot_annotation(title = "dSV Bacteroides xylanisolvens XB1A:337 ~ 339;339 ~ 340 kbp\n OR=0.14;p=2.10e-4")

p_whole<-p1+plot_layout(ncol=1) 

tiff(file = "pics/sites/F5_dsv_cirrhosis_Bacteroides_xylanisolvens_barplot.tiff", width =900, height =1000, res =300) 
p_whole
dev.off()



##################################################################################
### vsgv picture 

all_basic$id<-row.names(all_basic)

vsgv_pic<-vsgv[,c(
"Bacteroides uniformis ATCC 8492:3885_3890",
"Parabacteroides distasonis ATCC 8503:4735_4737"
)]


vsgv_pic$id<-row.names(vsgv_pic)
pic_vsv<-merge(all_basic,vsgv_pic,by.x="id",by.y="id")


#################################################################################
## Bacteroides uniformis ATCC 8492:3885_3890

pic_vsv$temp<-pic_vsv$"Bacteroides uniformis ATCC 8492:3885_3890"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(cirrhosis)==F,drop=T)

pic_sub$cirrhosis<-as.factor(pic_sub$cirrhosis)
flevels <- levels(pic_sub$cirrhosis)
table(pic_sub$cirrhosis)

#pic_sub<-subset(pic_sub,temp<10,drop=T)

pic1<-ggplot(data=pic_sub,aes(x=cirrhosis,
                         y=temp,
                         color=cirrhosis))+
  #geom_jitter(alpha=0.8,
  #            position=position_jitterdodge(jitter.width = 0.65, 
  #                                          jitter.height = 0, 
  #                                          dodge.width = 0.8))+
  geom_quasirandom(aes(color=cirrhosis),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),
              size=0.75)+
  scale_color_manual(values =  c("#97C24E","#98369E") )+
  scale_x_discrete(limits=flevels,labels=(c("health(n=103)","cirrhosis(n=77)")))+
  theme_classic() +
  theme(legend.position="none") + 
  theme(text = element_text(size=12)) + #ylim(0.0,1.3)+
  ylab("Normalized region coverage")+xlab("")+
  ggtitle("Bacteroides uniformis ATCC 8492:3885 ~ 3890")+
  theme(panel.grid = element_blank(), 
  panel.background = element_rect(color = 'black', fill = 'transparent'), 
  strip.text = element_text(size = 12)) +
  theme(axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(1,1,1,1), "cm"),
  axis.text.x = element_text(angle = 45, hjust = 1))

tiff(file = "pics/sites/cirrhosis_vsv_Bacteroides_uniformis.tiff", width =1500, height =1500, res =300) 
pic1
dev.off()




###################################################################################
## creatinine

#######################################################################################
### clinical data and dsgv/vsgv
all_basic$id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"Bacteroides fragilis NCTC 9343:1872_1875",
"Bacteroides xylanisolvens XB1A:510_512 and 14 segments")]

########################################################################################
##  
dsgv_pic$match_id<-row.names(dsgv_pic)
pic_dsv<-merge(all_basic,dsgv_pic,by.x="match_id",by.y="match_id")

#########################################################################
#####  Bacteroides fragilis NCTC 9343:1872_1875
pic_dsv$temp<-pic_dsv$"Bacteroides fragilis NCTC 9343:1872_1875"

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(creatinine)==F,drop=T)
table(pic_sub$study_name)

pic_sub$temp<-as.factor(pic_sub$temp)
flevels <- levels(pic_sub$temp)

JieZ_2017_pic<-subset(pic_sub,study_name=="JieZ_2017",drop=T)
QinN_2014_pic<-subset(pic_sub,study_name=="QinN_2014",drop=T)
YuJ_2015_pic<-subset(pic_sub,study_name=="YuJ_2015",drop=T)


table(JieZ_2017_pic$temp)

p1<-ggplot(data=JieZ_2017_pic,aes(x=as.factor(temp),y=creatinine,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("JieZ_2017")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#207F4C","#EA517F") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=25)","Non-delection(n=43)")))+
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40))+
  geom_hline(yintercept=c(40,80,120,160),linetype="dashed")+
  theme_cowplot() +
  ylab("creatinine")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=11,face="plain"),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(QinN_2014_pic$temp)
p2<-ggplot(data=QinN_2014_pic,aes(x=as.factor(temp),y=creatinine,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("QinN_2014")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#207F4C","#EA517F") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=32)","Non-delection(n=40)")))+
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40))+
  geom_hline(yintercept=c(40,80,120,160),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=11,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(YuJ_2015_pic$temp)
p3<-ggplot(data=YuJ_2015_pic,aes(x=as.factor(temp),y=creatinine,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("YuJ_2015")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#207F4C","#EA517F") )+
  scale_x_discrete(limits=flevels,labels=(c("Delection(n=17)","Non-delection(n=15)")))+
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 40))+
  geom_hline(yintercept=c(40,80,120,160),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=11,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


p_whole<-p1+p2+p3+plot_layout(ncol=3) + plot_annotation(title = "Bacteroides fragilis NCTC 9343:1872 ~ 1875 kbp\nMeta OR=8.15;Meta p=6.23e-4")

tiff(file = "pics/sites/F5_dSV_creatinine_Bacteroides_fragilis_barplot.tiff", width =1500, height =1500, res =300) 
p_whole
dev.off()













