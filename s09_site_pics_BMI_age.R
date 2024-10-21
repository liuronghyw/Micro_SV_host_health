### barplot for sites (BMI and age)
### 2024-9-6
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

## 1 Preparation
##rm(list = ls())

library("survival")
library("survminer")
library(gridExtra)
library("grid")
#library(reshape2)	  
#library("RColorBrewer")
##library("plyr")

#install.packages("patchwork")
library("patchwork")
source("functions.R")


###################### 1.2 Inputs
info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

vsv_info<-read.table("01.cleanData/SV_info/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

all_basic<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv",header=T)

### just includ healthy patients 
all_basic<-subset(all_basic,disease=="healthy",drop=T)

load("01.cleanData/SV_all/dsgv.RData")
load("01.cleanData/SV_all/vsgv.RData")

vsgv<-vsgv_sub
dsgv<-dsgv_sub


#######################################################################################
### clinical data and dsgv/vsgv
all_basic$match_id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"[Eubacterium] hallii DSM 3353:285_287 and 5 segments",
"Bifidobacterium pseudocatenulatum DSM 20438 = JCM 1200 = LMG 10505:1518_1520 and 2 segments")]

### draw pictures
##install.packages('ggbeeswarm')
library("ggbeeswarm")

########################################################################################
##  
dsgv_pic$match_id<-row.names(dsgv_pic)
pic_dsv<-merge(all_basic,dsgv_pic,by.x="match_id",by.y="match_id")


#########################################################################
##### Bifidobacterium pseudocatenulatum DSM 20438 = JCM 1200 = LMG 10505:1518_1520 and 2 segments

pic_dsv$temp<-pic_dsv$"Bifidobacterium pseudocatenulatum DSM 20438 = JCM 1200 = LMG 10505:1518_1520 and 2 segments"

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(BMI)==F,drop=T)
table(pic_sub$study_name)

pic_sub$temp<-as.factor(pic_sub$temp)
flevels <- levels(pic_sub$temp)

AsnicarF_2021_pic<-subset(pic_sub,study_name=="AsnicarF_2021",drop=T)
CosteaPI_2017_pic<-subset(pic_sub,study_name=="CosteaPI_2017",drop=T)
SchirmerM_2016_pic<-subset(pic_sub,study_name=="SchirmerM_2016",drop=T)
YachidaS_2019_pic<-subset(pic_sub,study_name=="YachidaS_2019",drop=T)
ZeeviD_2015_pic<-subset(pic_sub,study_name=="ZeeviD_2015",drop=T)

table(AsnicarF_2021_pic$temp)

p1<-ggplot(data=AsnicarF_2021_pic,aes(x=as.factor(temp),y=BMI,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("AsnicarF_2021")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#207F4C","#EA517F") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=65)","Non-delection(n=48)")))+
  scale_y_continuous(limits = c(15, 45), breaks = seq(15, 45, 15))+
  geom_hline(yintercept=c(15,30,45),linetype="dashed")+
  theme_cowplot() +
  ylab("BMI")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(CosteaPI_2017_pic$temp)
p2<-ggplot(data=CosteaPI_2017_pic,aes(x=as.factor(temp),y=BMI,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("CosteaPI_2017")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#207F4C","#EA517F") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=20)","Non-delection(n=12)")))+
  scale_y_continuous(limits = c(15, 45), breaks = seq(15, 45, 15))+
  geom_hline(yintercept=c(15,30,45),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(SchirmerM_2016_pic$temp)
p3<-ggplot(data=SchirmerM_2016_pic,aes(x=as.factor(temp),y=BMI,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("SchirmerM_2016")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#207F4C","#EA517F") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=29)","Non-delection(n=22)")))+
  scale_y_continuous(limits = c(15, 45), breaks = seq(15, 45, 15))+
  geom_hline(yintercept=c(15,30,45),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))

table(YachidaS_2019_pic$temp)
p4<-ggplot(data=YachidaS_2019_pic,aes(x=as.factor(temp),y=BMI,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("YachidaS_2019")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#207F4C","#EA517F") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=41)","Non-delection(n=27)")))+
  scale_y_continuous(limits = c(15, 45), breaks = seq(15, 45, 15))+
  geom_hline(yintercept=c(15,30,45),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


table(ZeeviD_2015_pic$temp)
p5<-ggplot(data=ZeeviD_2015_pic,aes(x=as.factor(temp),y=BMI,color=as.factor(temp)))+
  geom_quasirandom(aes(color=as.factor(temp)),width = 0.2,size=1,alpha=0.8)+
  geom_boxplot(alpha=0.2,width=0.25,color="black",
               position=position_dodge(width=0.8),
               size=0.75,outlier.colour = NA)+ggtitle("ZeeviD_2015")+
  geom_violin(alpha=0.05,width=0.9,
              position=position_dodge(width=0.8),size=0.75)+
  scale_color_manual(values =  c("#207F4C","#EA517F") )+
  scale_x_discrete(limits=flevels,labels=(c("Detection(n=111)","Non-delection(n=55)")))+
  scale_y_continuous(limits = c(15, 45), breaks = seq(15, 45, 15))+
  geom_hline(yintercept=c(15,30,45),linetype="dashed")+
  theme_cowplot() +
  ylab("")+xlab("")+
  theme(panel.grid = element_blank(),legend.position="none",text = element_text(size=12),
  strip.text = element_text(size = 12),axis.text = element_text(size = 12),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(0,0,0,0), "pt"), plot.title=element_text(size=13,face="plain"),
  axis.text.y = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1))


#p_whole<-p1+p2+p3+p4+p5+plot_layout(ncol=5) + plot_annotation(title = "dSV Bifidobacterium pseudocatenulatum DSM 20438 = JCM 1200 = LMG 10505:\n1518 ~ 1520 and 2 segments kbp\nMeta OR=1.49;Meta p=1.25e-5")

p_whole<-p1+p2+p3+p4+p5+plot_layout(ncol=5) 

tiff(file = "pics/sites/F5_dSV_BMI_Bifidobacterium_barplot.tiff", width =2200, height =1500, res =300) 
p_whole
dev.off()



##################################################################################
### vsgv picture
all_basic$match_id<-row.names(all_basic)

vsgv_pic<-vsgv[,c(
"Faecalibacterium cf. prausnitzii KLE1255:1372_1379",
"Ruminococcus torques L2-14:514_516"
)]

vsgv_pic$match_id<-row.names(vsgv_pic)
### mel 
pic_vsv<-merge(all_basic,vsgv_pic,by.x="match_id",by.y="match_id")

#################################################################################
## Faecalibacterium cf. prausnitzii KLE1255:1372_1379

pic_vsv$temp<-pic_vsv$"Faecalibacterium cf. prausnitzii KLE1255:1372_1379"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(BMI)==F,drop=T)

table(pic_sub$study_name)

pic_sub<-subset(pic_vsv,study_name=="AsnicarF_2021" | study_name=="HansenLBS_2018"
                        | study_name=="MetaCardis_2020_a" | study_name=="ZeeviD_2015",drop=T)


# Save the scatter plot in a variable
p <- ggplot(pic_sub, aes(x = temp, y = BMI,color=study_name)) +
  geom_point(aes(color=study_name))+
  scale_color_manual(values = c("#ECD452","#BB97C5","#EA8958","#2983BB")) +
  scale_fill_manual(values = c("#ECD452","#BB97C5","#EA8958","#2983BB")) +
  geom_smooth(method = lm,aes(color=study_name,fill=study_name))+
  guides(col = guide_legend(nrow =2, byrow = TRUE))+
  theme_classic() +
  #theme(legend.position.inside=c(0,65)) + 
  theme(text = element_text(size=16),legend.position="bottom") + 
  #ylim(0.0,1.3)+
  xlab("Normalized region coverage")+
  ylab("BMI")+
  #ggtitle("Faecalibacterium cf. prausnitzii KLE1255:\n1372 ~ 1379 kbp")+
  theme(panel.grid = element_blank(), 
  panel.background = element_rect(color = 'white', fill = 'transparent'), 
  strip.text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(1,1,1,1), "cm"),
  axis.text.x = element_text(angle = 0, hjust = 1))

tiff(file = "pics/sites/Fig6_BMI_Fprausnitzii_vsv_1379.tiff", width =1800, height =1800, res =300) 
p
dev.off()




##############################################################################
### age
## dsv
## Bifidobacterium adolescentis:1660_1661
## Bacteroides caccae ATCC 43185:4267_4269 and 3 segments

### define age categray

all_basic$age_cat<-NA
nr<-nrow(all_basic)

for(i in 1:nr){
     if(is.na(all_basic$age[i])==T){all_basic$age_cat[i]=NA}
     if(is.na(all_basic$age[i])==F){
      if (all_basic$age[i]>=15 & all_basic$age[i]<30){all_basic$age_cat[i]="15 <= Age <30"}
      else if (all_basic$age[i]>=30 & all_basic$age[i]<45){all_basic$age_cat[i]="30 <= Age <45"}
      else if (all_basic$age[i]>=45 & all_basic$age[i]<60){all_basic$age_cat[i]="45 <= Age <60"}
      else if (all_basic$age[i]>=60 & all_basic$age[i]<75){all_basic$age_cat[i]="60 <= Age <75"}
      else if (all_basic$age[i]>=75 & all_basic$age[i]<90){all_basic$age_cat[i]="75 <= Age <90"}
    }
    else {all_basic$age_cat[i]=NA}
 }



#######################################################################################
### clinical data and dsgv/vsgv
all_basic$match_id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"Bifidobacterium adolescentis:1660_1661",
"Bacteroides caccae ATCC 43185:4267_4269 and 3 segments")]

### draw pictures

########################################################################################
##  age
dsgv_pic$match_id<-row.names(dsgv_pic)
pic_dsv<-merge(all_basic,dsgv_pic,by.x="match_id",by.y="match_id")

#########################################################################
##### Bifidobacterium adolescentis:1660_1661

pic_dsv$temp<-pic_dsv$"Bifidobacterium adolescentis:1660_1661"

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(age_cat)==F,drop=T)
table(pic_sub$study_name)

pic_sub$temp<-as.factor(pic_sub$temp)
flevels <- levels(pic_sub$temp)

AsnicarF_2021_pic<-subset(pic_sub,study_name=="AsnicarF_2021",drop=T)
CosteaPI_2017_pic<-subset(pic_sub,study_name=="CosteaPI_2017",drop=T)
SchirmerM_2016_pic<-subset(pic_sub,study_name=="SchirmerM_2016",drop=T)

table(AsnicarF_2021_pic$temp)

AsnicarF_2021_ratio<-dsv_age(AsnicarF_2021_pic,"AsnicarF_2021")
CosteaPI_2017_ratio<-dsv_age(CosteaPI_2017_pic,"CosteaPI_2017")
SchirmerM_2016_ratio<-dsv_age(SchirmerM_2016_pic,"SchirmerM_2016")


table(AsnicarF_2021_pic$age_cat)

p1<-ggplot(AsnicarF_2021_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("AsnicarF_2021")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = 'Ratio')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("15 <= Age <30","30 <= Age <45","45 <= Age <60","60 <= Age <75"),
   labels=c("15 <= Age <30 (n=116)","30 <= Age <45 (n=207)","45 <= Age <60 (n=274)","60 <= Age <75 (n=66)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black"),text=element_text(family="sans"),legend.position="none",
   plot.title=element_text(size=11,face="plain")
   )

table(CosteaPI_2017_pic$age_cat)
p2<-ggplot(CosteaPI_2017_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("CosteaPI_2017")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("15 <= Age <30","30 <= Age <45","45 <= Age <60","60 <= Age <75"),
   labels=c("15 <= Age <30 (n=4)","30 <= Age <45 (n=39)","45 <= Age <60 (n=32)","60 <= Age <75 (n=9)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black"),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain")
   )


table(SchirmerM_2016_pic$age_cat)
p3<-ggplot(SchirmerM_2016_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("SchirmerM_2016")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("15 <= Age <30","30 <= Age <45","45 <= Age <60","60 <= Age <75"),
   labels=c("15 <= Age <30 (n=298)","30 <= Age <45 (n=17)","45 <= Age <60 (n=23)","60 <= Age <75 (n=20)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black"),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain")
   )

#p_whole<-p1+p2+p3+plot_layout(ncol=3) + plot_annotation(title = "Bifidobacterium adolescentis:1660 ~ 1661 kbp\nMeta beta=-2.77;Meta p=5.61e-5")

p_whole<-p1+p2+p3+plot_layout(ncol=3) 

#p_whole<-plot_grid(p1,p2,p3,p4,
#              rel_widths = c(0.55,0.45,0.45,0.7),align = 'hv',
#              ncol = 4,label_size= 10,vjust = 0)

tiff(file = "pics/sites/F5_dSV_Bifidobacterium_adolescentis_age_barplot.tiff", width =1700, height =1300, res =300) 
p_whole
dev.off()






#########################################################################
##### Bacteroides caccae ATCC 43185:4267_4269 and 3 segments

pic_dsv$temp<-pic_dsv$"Bacteroides caccae ATCC 43185:4267_4269 and 3 segments"

pic_sub<-subset(pic_dsv,is.na(temp)==F & is.na(age_cat)==F,drop=T)
table(pic_sub$study_name)

pic_sub$temp<-as.factor(pic_sub$temp)
flevels <- levels(pic_sub$temp)

AsnicarF_2021_pic<-subset(pic_sub,study_name=="AsnicarF_2021",drop=T)
HansenLBS_2018_pic<-subset(pic_sub,study_name=="HansenLBS_2018",drop=T)
YachidaS_2019_pic<-subset(pic_sub,study_name=="YachidaS_2019",drop=T)


AsnicarF_2021_ratio<-dsv_age(AsnicarF_2021_pic,"AsnicarF_2021")
HansenLBS_2018_ratio<-dsv_age(HansenLBS_2018_pic,"HansenLBS_2018")
YachidaS_2019_ratio<-dsv_age(YachidaS_2019_pic,"YachidaS_2019")

table(AsnicarF_2021_pic$age_cat)

p1<-ggplot(AsnicarF_2021_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("AsnicarF_2021")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = 'Ratio')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("15 <= Age <30","30 <= Age <45","45 <= Age <60","60 <= Age <75"),
   labels=c("15 <= Age <30 (n=113)","30 <= Age <45 (n=192)","45 <= Age <60 (n=280)","60 <= Age <75 (n=83)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1,color="black"),text=element_text(family="sans"),legend.position="none",
   plot.title=element_text(size=11,face="plain")
   )

table(HansenLBS_2018_pic$age_cat)
p2<-ggplot(HansenLBS_2018_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("HansenLBS_2018")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("15 <= Age <30","30 <= Age <45","45 <= Age <60","60 <= Age <75"),
   labels=c("15 <= Age <30 (n=13)","30 <= Age <45 (n=15)","45 <= Age <60 (n=47)","60 <= Age <75 (n=11)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1,color="black"),text=element_text(family="sans"),legend.position="none",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain")
   )


table(YachidaS_2019_pic$age_cat)
p3<-ggplot(YachidaS_2019_ratio, aes(chara, as.numeric(value), fill = variable)) +
   geom_col(position = 'stack', width =0.9) +ggtitle("YachidaS_2019")+
   scale_fill_manual(breaks=c("Delection","Non.delection"),labels=c("Delection","Non.delection"),values =c("#7876B1FF","#207F4C")) +
   labs(x = '', y = '')+
   geom_text(aes(y=new_col, label=value), size=3.5,colour="White")+
   scale_x_discrete(limits=c("15 <= Age <30","30 <= Age <45","45 <= Age <60","60 <= Age <75"),
   labels=c("15 <= Age <30 (n=8)","30 <= Age <45 (n=24)","45 <= Age <60 (n=38)","60 <= Age <75 (n=12)"),expand=c(0,0))+
   theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
   strip.text = element_text(size = 10)) +
   theme_pander()+
   theme(axis.text.x = element_text(angle = 50, hjust = 1,color="black"),text=element_text(family="sans"),legend.position="right",
   axis.text.y = element_blank(),plot.title=element_text(size=11,face="plain")
   )

p_whole<-p1+p2+p3+plot_layout(ncol=3) + plot_annotation(title = "Bacteroides caccae ATCC 43185:4267 ~ 4269 and 3 segments kbp\nMeta beta=-2.36;Meta p=2.91e-4")

#p_whole<-plot_grid(p1,p2,p3,p4,
#              rel_widths = c(0.55,0.45,0.45,0.7),align = 'hv',
#              ncol = 4,label_size= 10,vjust = 0)

tiff(file = "pics/sites/F5_dSV_Bacteroides_caccae_age_barplot.tiff", width =1500, height =1500, res =300) 
p_whole
dev.off()




###########################################################################################################
### vSV
## Faecalibacterium cf. prausnitzii KLE1255:117_118 and 14 segments
## Bifidobacterium longum:119_122;129_130
## Bifidobacterium adolescentis:422_423 and 4 segments
## "Faecalibacterium cf. prausnitzii KLE1255:1528_1534",
## "Faecalibacterium cf. prausnitzii KLE1255:1350_1353;1353_1362"


##################################################################################
### vsgv picture
all_basic$match_id<-row.names(all_basic)

vsgv_pic<-vsgv[,c(
"Faecalibacterium cf. prausnitzii KLE1255:117_118 and 14 segments",
"Bifidobacterium longum:119_122;129_130",
"Bifidobacterium adolescentis:422_423 and 4 segments",
"Faecalibacterium cf. prausnitzii KLE1255:1528_1534",
"Faecalibacterium cf. prausnitzii KLE1255:1350_1353;1353_1362"
)]

vsgv_pic$match_id<-row.names(vsgv_pic)
pic_vsv<-merge(all_basic,vsgv_pic,by.x="match_id",by.y="match_id")

#################################################################################
## Faecalibacterium cf. prausnitzii KLE1255:117_118 and 14 segments

pic_vsv$temp<-pic_vsv$"Faecalibacterium cf. prausnitzii KLE1255:117_118 and 14 segments"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(age)==F,drop=T)

table(pic_sub$study_name)

pic_sub<-subset(pic_vsv,study_name=="AsnicarF_2021" | study_name=="KeohaneDM_2020"
                        | study_name=="MetaCardis_2020_a" | study_name=="ZhuF_2020",drop=T)


# Save the scatter plot in a variable
p <- ggplot(pic_sub, aes(x = temp, y = age,color=study_name)) +
  geom_point(aes(color=study_name))+
  scale_color_manual(values = c("#ECD452","#BB97C5","#EA8958","#2983BB")) +
  scale_fill_manual(values = c("#ECD452","#BB97C5","#EA8958","#2983BB")) +
  geom_smooth(method = lm,aes(color=study_name,fill=study_name))+
  guides(col = guide_legend(nrow =2, byrow = TRUE))+
  theme_classic() +
  #theme(legend.position.inside=c(0,65)) + 
  theme(text = element_text(size=16),legend.position="bottom") + 
  #ylim(0.0,1.3)+
  xlab("Normalized region coverage")+
  ylab("Age (years)")+
  ggtitle("Faecalibacterium cf. prausnitzii KLE1255:117 ~ 118 and 14 segments")+
  theme(panel.grid = element_blank(), 
  panel.background = element_rect(color = 'black', fill = 'transparent'), 
  strip.text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(1,1,1,1), "cm"),
  axis.text.x = element_text(angle = 0, hjust = 1))

tiff(file = "pics/sites/Fig6_age_Fprausnitzii_vsv_1379.tiff", width =1800, height =1800, res =300) 
p
dev.off()



######################################################################################
## "Faecalibacterium cf. prausnitzii KLE1255:1528_1534",

pic_vsv$temp<-pic_vsv$"Faecalibacterium cf. prausnitzii KLE1255:1528_1534"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(age)==F,drop=T)

table(pic_sub$study_name)

pic_sub<-subset(pic_sub,study_name=="AsnicarF_2021" | study_name=="HansenLBS_2018"
                        | study_name=="CosteaPI_2017" | study_name=="ZhuF_2020"
                        | study_name=="HMP_2012" | study_name=="MetaCardis_2020_a"
                        | study_name=="YachidaS_2019",drop=T)


# Save the scatter plot in a variable
p <- ggplot(pic_sub, aes(x = temp, y = age,color=study_name)) +
  geom_point(aes(color=study_name))+
  scale_color_manual(values = c("#ECD452","#BB97C5","#EA8958","#2983BB","#BC3523","#207F4C","#EE4863")) +
  geom_smooth(method = lm)+
  guides(col = guide_legend(nrow =3, byrow = TRUE))+
  theme_classic() +
  #theme(legend.position.inside=c(0,65)) + 
  theme(text = element_text(size=16),legend.position="bottom") + 
  #ylim(0.0,1.3)+
  xlab("Normalized region coverage")+
  ylab("Age (years)")+
  ggtitle("Faecalibacterium cf. prausnitzii KLE1255:\n1528 ~ 1534 kbp")+
  theme(panel.grid = element_blank(), 
  panel.background = element_rect(color = 'black', fill = 'transparent'), 
  strip.text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(1,1,1,1), "cm"),
  axis.text.x = element_text(angle = 0, hjust = 1))

tiff(file = "pics/sites/Fig6_age_Fprausnitzii_vsv_1379.tiff", width =1800, height =1800, res =300) 
p
dev.off()


##############################################################################################
## "Faecalibacterium cf. prausnitzii KLE1255:1350_1353;1353_1362"

pic_vsv$temp<-pic_vsv$"Faecalibacterium cf. prausnitzii KLE1255:1350_1353;1353_1362"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(age)==F,drop=T)

table(pic_sub$study_name)

pic_sub<-subset(pic_sub,study_name=="AsnicarF_2021" | study_name=="CosteaPI_2017" 
                        | study_name=="ZhuF_2020"| study_name=="MetaCardis_2020_a" 
                        | study_name=="SchirmerM_2016" | study_name=="ZhuF_2020",drop=T)


# Save the scatter plot in a variable
p <- ggplot(pic_sub, aes(x = temp, y = age,color=study_name)) +
  geom_point(aes(color=study_name))+
  scale_color_manual(values = c("#ECD452","#BB97C5","#EA8958","#2983BB","#BC3523","#207F4C")) +
  geom_smooth(method = lm)+
  guides(col = guide_legend(nrow =3, byrow = TRUE))+
  theme_classic() +
  #theme(legend.position.inside=c(0,65)) + 
  theme(text = element_text(size=16),legend.position="bottom") + 
  #ylim(0.0,1.3)+
  xlab("Normalized region coverage")+
  ylab("Age (years)")+
  ggtitle("Faecalibacterium cf. prausnitzii KLE1255:\n1353 ~ 1362 kbp")+
  theme(panel.grid = element_blank(), 
  panel.background = element_rect(color = 'black', fill = 'transparent'), 
  strip.text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(1,1,1,1), "cm"),
  axis.text.x = element_text(angle = 0, hjust = 1))

tiff(file = "pics/sites/Fig6_age_Fprausnitzii_vsv_1353.tiff", width =1800, height =1800, res =300) 
p
dev.off()



#################################################################################
## Bifidobacterium longum:119_122;129_130

pic_vsv$temp<-pic_vsv$"Bifidobacterium longum:119_122;129_130"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(age)==F,drop=T)

table(pic_sub$study_name)

pic_sub<-subset(pic_sub,study_name=="AsnicarF_2021" | study_name=="YachidaS_2019",drop=T)


# Save the scatter plot in a variable
p <- ggplot(pic_sub, aes(x = temp, y = age,color=study_name)) +
  geom_point(aes(color=study_name))+
  scale_color_manual(values = c("#ECD452","#BB97C5")) +
  geom_smooth(method = lm)+
  guides(col = guide_legend(nrow =1, byrow = TRUE))+
  theme_classic() +
  #theme(legend.position.inside=c(0,65)) + 
  theme(text = element_text(size=16),legend.position="bottom") + 
  #ylim(0.0,1.3)+
  xlab("Normalized region coverage")+
  ylab("Age (years)")+
  ggtitle("Bifidobacterium longum:119_122;129_130")+
  theme(panel.grid = element_blank(), 
  panel.background = element_rect(color = 'black', fill = 'transparent'), 
  strip.text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(1,1,1,1), "cm"),
  axis.text.x = element_text(angle = 0, hjust = 1))

tiff(file = "pics/sites/Fig6_age_Bifidobacterium_vsv_130.tiff", width =1800, height =1800, res =300) 
p
dev.off()




#################################################################################
## Bifidobacterium adolescentis:422_423 and 4 segments

pic_vsv$temp<-pic_vsv$"Bifidobacterium adolescentis:422_423 and 4 segments"
pic_sub<-subset(pic_vsv,is.na(temp)==F & is.na(age)==F,drop=T)

table(pic_sub$study_name)

pic_sub<-subset(pic_sub,study_name=="AsnicarF_2021" | study_name=="MetaCardis_2020_a" |
                        study_name=="SchirmerM_2016" | study_name=="YachidaS_2019",drop=T)


# Save the scatter plot in a variable
p <- ggplot(pic_sub, aes(x = temp, y = age,color=study_name)) +
  geom_point(aes(color=study_name))+
  scale_color_manual(values = c("#ECD452","#BB97C5","#EA8958","#2983BB")) +
  scale_fill_manual(values = c("#ECD452","#BB97C5","#EA8958","#2983BB")) +
  geom_smooth(method = lm,aes(color=study_name,fill=study_name))+
  guides(col = guide_legend(nrow =2, byrow = TRUE))+
  theme_classic() +
  #theme(legend.position.inside=c(0,65)) + 
  theme(text = element_text(size=16),legend.position="bottom") + 
  #ylim(0.0,1.3)+
  xlab("Normalized region coverage")+
  ylab("Age (years)")+
  #ggtitle("Bifidobacterium adolescentis:\n422 ~ 423 and 4 segments")+
  theme(panel.grid = element_blank(), 
  panel.background = element_rect(color = 'white', fill = 'transparent'), 
  strip.text = element_text(size = 15)) +
  theme(axis.text = element_text(size = 15),  
  axis.title = element_text(size = 12), legend.title = element_blank(), 
  plot.margin = unit(c(1,1,1,1), "cm"),
  axis.text.x = element_text(angle = 0, hjust = 1))

tiff(file = "pics/sites/Fig6_age_Bifidobacterium_adolescentis_vsv_130.tiff", width =1800, height =1800, res =300) 
p
dev.off()






