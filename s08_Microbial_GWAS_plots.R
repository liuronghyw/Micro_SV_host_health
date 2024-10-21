### Microbial GWA 
### 2024-9-4
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

## 1 Preparation
### 1.1 Import
#install.packages("survminer")
rm(list = ls())
#options(scipen = 200)

source("functions.R")
knitr::opts_chunk$set(echo = TRUE)
library(ggthemes)

#library("R.utils")
#library(showtext)
#showtext_auto()

#library("survival")
#library("ggplot2")
#library("survminer")
#library(gridExtra)
#library("grid")
library(reshape2)	  
#library("RColorBrewer")
#library("plyr")

###################### 1.2 Inputs
info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

vsv_info<-read.table("01.cleanData/SV_info/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")

all_basic<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")


########################################################
#####  3.3 Visualization 
####   3.3.1 Combine vSV and dSV associations

###############################################################################################
## all phenotype
## CRC  IBD  T2D  IGT ACVD severe_obesity MetS  CAD
## age BMI fasting_insulin hba1c fasting_glucose
## triglycerides hdl ldl  cholesterol dyastolic_p systolic_p
## hscrp albumine creatinine ast alt total_bilirubin

CRC_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/CRC_vsv.sig.anno.csv",header=T)
CRC_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/CRC_dsv.sig.anno.csv",header=T)
CRC_vsv_meta_pvalue.sig$Pheno<-"CRC"
CRC_dsv_meta_pvalue.sig$Pheno<-"CRC"
CRC_dsv_sub<-CRC_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
CRC_vsv_sub<-CRC_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

#adenoma_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/adenoma_vsv.sig.anno.csv",header=T)
adenoma_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/adenoma_dsv.sig.anno.csv",header=T)
#adenoma_vsv_meta_pvalue.sig$Pheno<-"adenoma"
adenoma_dsv_meta_pvalue.sig$Pheno<-"adenoma"
adenoma_dsv_sub<-adenoma_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
#adenoma_vsv_sub<-adenoma_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

IBD_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/IBD_vsv.sig.anno.csv",header=T)
IBD_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/IBD_dsv.sig.anno.csv",header=T)
IBD_vsv_meta_pvalue.sig$Pheno<-"IBD"
IBD_dsv_meta_pvalue.sig$Pheno<-"IBD"
IBD_dsv_sub<-IBD_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
IBD_vsv_sub<-IBD_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

#T2D_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/T2D_vsv.sig.anno.csv",header=T)
#T2D_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/T2D_dsv.sig.anno.csv",header=T)
#T2D_vsv_meta_pvalue.sig$Pheno<-"T2D"
#T2D_dsv_meta_pvalue.sig$Pheno<-"T2D"
#T2D_dsv_sub<-T2D_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
#T2D_vsv_sub<-T2D_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

ACVD_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/ACVD_vsv.sig.anno.csv",header=T)
ACVD_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/ACVD_dsv.sig.anno.csv",header=T)
ACVD_vsv_meta_pvalue.sig$Pheno<-"ACVD"
ACVD_dsv_meta_pvalue.sig$Pheno<-"ACVD"
ACVD_dsv_sub<-ACVD_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
ACVD_vsv_sub<-ACVD_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

#severe_obesity_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/severe_obesity_vsv.sig.anno.csv",header=T)
#severe_obesity_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/severe_obesity_dsv.sig.anno.csv",header=T)
#severe_obesity_vsv_meta_pvalue.sig$Pheno<-"severe_obesity"
#severe_obesity_dsv_meta_pvalue.sig$Pheno<-"severe_obesity"
#severe_obesity_dsv_sub<-severe_obesity_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
#severe_obesity_vsv_sub<-severe_obesity_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

#CAD_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/CAD_vsv.sig.anno.csv",header=T)
#CAD_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/CAD_dsv.sig.anno.csv",header=T)
#CAD_vsv_meta_pvalue.sig$Pheno<-"CAD"
#CAD_dsv_meta_pvalue.sig$Pheno<-"CAD"
#CAD_dsv_sub<-CAD_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
#CAD_vsv_sub<-CAD_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

#MetS_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/MetS_vsv.sig.anno.csv",header=T)
#MetS_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/MetS_dsv.sig.anno.csv",header=T)
#MetS_vsv_meta_pvalue.sig$Pheno<-"MetS"
#MetS_dsv_meta_pvalue.sig$Pheno<-"MetS"
#MetS_dsv_sub<-MetS_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
#MetS_vsv_sub<-MetS_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

cirrhosis_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/cirrhosis_vsv.sig.anno.csv",header=T)
cirrhosis_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/cirrhosis_dsv.sig.anno.csv",header=T)
cirrhosis_vsv_meta_pvalue.sig$Pheno<-"cirrhosis"
cirrhosis_dsv_meta_pvalue.sig$Pheno<-"cirrhosis"
cirrhosis_dsv_sub<-cirrhosis_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
cirrhosis_vsv_sub<-cirrhosis_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

age_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/age_healthy_vsv.sig.anno.csv",header=T)
age_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/age_healthy_dsv.sig.anno.csv",header=T)
age_vsv_meta_pvalue.sig$Pheno<-"age"
age_dsv_meta_pvalue.sig$Pheno<-"age"
age_dsv_sub<-age_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
age_vsv_sub<-age_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

BMI_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/BMI_healthy_vsv.sig.anno.csv",header=T)
BMI_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/BMI_healthy_dsv.sig.anno.csv",header=T)
BMI_vsv_meta_pvalue.sig$Pheno<-"BMI"
BMI_dsv_meta_pvalue.sig$Pheno<-"BMI"
BMI_dsv_sub<-BMI_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
BMI_vsv_sub<-BMI_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

#fasting_insulin_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/fasting_insulin_vsv.sig.anno.csv",header=T)
#fasting_insulin_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/fasting_insulin_dsv.sig.anno.csv",header=T)
#fasting_insulin_vsv_meta_pvalue.sig$Pheno<-"fasting_insulin"
#fasting_insulin_dsv_meta_pvalue.sig$Pheno<-"fasting_insulin"
#fasting_insulin_dsv_sub<-fasting_insulin_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
#fasting_insulin_vsv_sub<-fasting_insulin_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

hba1c_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/hba1c_vsv.sig.anno.csv",header=T)
hba1c_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/hba1c_dsv.sig.anno.csv",header=T)
hba1c_vsv_meta_pvalue.sig$Pheno<-"hba1c"
hba1c_dsv_meta_pvalue.sig$Pheno<-"hba1c"
hba1c_dsv_sub<-hba1c_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
hba1c_vsv_sub<-hba1c_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

#fasting_glucose_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/fasting_glucose_vsv.sig.anno.csv",header=T)
fasting_glucose_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/fasting_glucose_dsv.sig.anno.csv",header=T)
#fasting_glucose_vsv_meta_pvalue.sig$Pheno<-"fasting_glucose"
fasting_glucose_dsv_meta_pvalue.sig$Pheno<-"fasting_glucose"
fasting_glucose_dsv_sub<-fasting_glucose_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
#fasting_glucose_vsv_sub<-fasting_glucose_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

triglycerides_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/triglycerides_vsv.sig.anno.csv",header=T)
triglycerides_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/triglycerides_dsv.sig.anno.csv",header=T)
triglycerides_vsv_meta_pvalue.sig$Pheno<-"triglycerides"
triglycerides_dsv_meta_pvalue.sig$Pheno<-"triglycerides"
triglycerides_dsv_sub<-triglycerides_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
triglycerides_vsv_sub<-triglycerides_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

hdl_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/hdl_vsv.sig.anno.csv",header=T)
hdl_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/hdl_dsv.sig.anno.csv",header=T)
hdl_vsv_meta_pvalue.sig$Pheno<-"hdl"
hdl_dsv_meta_pvalue.sig$Pheno<-"hdl"
hdl_dsv_sub<-hdl_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
hdl_vsv_sub<-hdl_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

ldl_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/ldl_vsv.sig.anno.csv",header=T)
ldl_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/ldl_dsv.sig.anno.csv",header=T)
ldl_vsv_meta_pvalue.sig$Pheno<-"ldl"
ldl_dsv_meta_pvalue.sig$Pheno<-"ldl"
ldl_dsv_sub<-ldl_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
ldl_vsv_sub<-ldl_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

cholesterol_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/cholesterol_vsv.sig.anno.csv",header=T)
cholesterol_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/cholesterol_dsv.sig.anno.csv",header=T)
cholesterol_vsv_meta_pvalue.sig$Pheno<-"cholesterol"
cholesterol_dsv_meta_pvalue.sig$Pheno<-"cholesterol"
cholesterol_dsv_sub<-cholesterol_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
cholesterol_vsv_sub<-cholesterol_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

mean_press_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/mean_press_vsv.sig.anno.csv",header=T)
mean_press_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/mean_press_dsv.sig.anno.csv",header=T)
mean_press_vsv_meta_pvalue.sig$Pheno<-"mean_press"
mean_press_dsv_meta_pvalue.sig$Pheno<-"mean_press"
mean_press_dsv_sub<-mean_press_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
mean_press_vsv_sub<-mean_press_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

creatinine_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/creatinine_vsv.sig.anno.csv",header=T)
creatinine_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/creatinine_dsv.sig.anno.csv",header=T)
creatinine_vsv_meta_pvalue.sig$Pheno<-"creatinine"
creatinine_dsv_meta_pvalue.sig$Pheno<-"creatinine"
creatinine_dsv_sub<-creatinine_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
creatinine_vsv_sub<-creatinine_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

#ast_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/ast_vsv.sig.anno.csv",header=T)
#ast_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/ast_dsv.sig.anno.csv",header=T)
#ast_vsv_meta_pvalue.sig$Pheno<-"ast"
#ast_dsv_meta_pvalue.sig$Pheno<-"ast"
#ast_dsv_sub<-ast_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
#ast_vsv_sub<-ast_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

alt_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/alt_vsv.sig.anno.csv",header=T)
#alt_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/alt_dsv.sig.anno.csv",header=T)
alt_vsv_meta_pvalue.sig$Pheno<-"alt"
#alt_dsv_meta_pvalue.sig$Pheno<-"alt"
#alt_dsv_sub<-alt_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
alt_vsv_sub<-alt_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]

#total_bilirubin_vsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/total_bilirubin_vsv.sig.anno.csv",header=T)
#total_bilirubin_dsv_meta_pvalue.sig<-read.csv("07.Microbial_GWAS/sig_noadjAbun/total_bilirubin_dsv.sig.anno.csv",header=T)
#total_bilirubin_vsv_meta_pvalue.sig$Pheno<-"total_bilirubin"
#total_bilirubin_dsv_meta_pvalue.sig$Pheno<-"total_bilirubin"
#total_bilirubin_dsv_sub<-total_bilirubin_dsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]
#total_bilirubin_vsv_sub<-total_bilirubin_vsv_meta_pvalue.sig[,c("Pheno","X","pvalue_meta","het_p")]


################################################################################
###  basic characteristic

sv_prog_meta_pvalue.sig.edge<-rbind(
age_dsv_sub,age_vsv_sub,
BMI_dsv_sub,BMI_vsv_sub)

## replace pheno value 
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "age"] = "Age"

colnames(sv_prog_meta_pvalue.sig.edge)[2]<-"SV"

sv_prog_adjAbun.sig.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno %>%
  as.character(.) %>%
  .[!duplicated(.)]

sv_prog_adjAbun.sig.sv <- sv_prog_meta_pvalue.sig.edge$SV %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

## circos plot
sv_prog_count<-NULL

for (phe in sv_prog_adjAbun.sig.pheno) {
  #pheno <-"C4"
  sv_prog_df<-sv_prog_meta_pvalue.sig.edge[sv_prog_meta_pvalue.sig.edge$Pheno==phe,]$SV %>%
    str_replace_all(.,"\\:\\d+_\\d+.*", "") %>%
    table %>%
    as.data.frame
  colnames(sv_prog_df)<-c("Species", "Count")
  sv_prog_df<-data.frame(Pheno = rep(phe,nrow(sv_prog_df)), sv_prog_df)
  sv_prog_count<-rbind(sv_prog_count,sv_prog_df)
}

sv_prog_count$Species<-info$Short_name[match(sv_prog_count$Species, info$organism)]
sv_prog_count <- sv_prog_count[order(sv_prog_count$Count),]

sv_prog_count_species_order<-sv_prog_count %>% group_by(Species) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]
sv_prog_count_pheno_order<-sv_prog_count %>% group_by(Pheno) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]

sv_prog_count<-sv_prog_count[order(match(sv_prog_count$Pheno, sv_prog_count_pheno_order$Pheno),decreasing = F),]

sv_prog_count_species_order_str<-as.character(sv_prog_count_species_order$Species)
sv_prog_count_pheno_order_str<-as.character(sv_prog_count_pheno_order$Pheno)


#pdf("pics/F6_SV_pheno_adjAbun.circos.pdf", width = 26, height = 26)
tiff(file = "pics/F6_SV_index_age_BMI_circos.tiff", width =7000, height =5800, res =300) 

circos.clear()
circos.par(start.degree = 180, "clock.wise" = T) # ,xaxis.clock.wise = F

chordDiagram(sv_prog_count,annotationTrack = "grid",
             #grid.col =  c(wes_palette("Darjeeling1", length(sv_prog_count_pheno_order_str), type = "continuous"),
             #              rep('grey',length(sv_prog_count_species_order_str))),
             grid.col =  c(mycolor2_green_blue,
                           rep('grey',length(sv_prog_count_species_order_str))),
             order = c(rev(sv_prog_count_pheno_order_str),
                       sv_prog_count_species_order_str),
             big.gap = 6,
             preAllocateTracks = list(track.margin = c(0, uh(45, "mm")), 
             #                         #track.height = max(strwidth(unlist(dimnames(sv_prog_count)))))
                                       track.height=0.1)
             )

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,cex = 2.5,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),font=3)
}, bg.border = NA) # here set bg.border to NA is important

dev.off()


#mycolor6
#hist(mtcars$mpg,
#     breaks = 15,
#     col ="#9fa39a",
#)


### "#9fa39a","#621624"
### "#ECD452","#BB97C5","#EA8958","#2983BB","#BC3523","#207F4C","#EE4863","#97C24E","#98369E","#63BBD0","#FEBA07","#EA517F"
### "#E18727FF"  ,"#EA517F","#BC3C29FF","#0072B5FF"

write.csv(sv_prog_count,"07.Microbial_GWAS/age_BMI_sv_prog_count.csv",row.names=F)
unique(sv_prog_count$Species)
unique(sv_prog_meta_pvalue.sig.edge$SV)



################################################################################
###  gastrointestinal disorders

## CRC  IBD  T2D  ACVD severe obesity MetS  CAD
## age BMI fasting_insulin hba1c fasting_glucose
## triglycerides hdl ldl  cholesterol mean_press
## hscrp creatine albumine creatinine ast alt total_bilirubin

sv_prog_meta_pvalue.sig.edge<-rbind(
CRC_dsv_sub,CRC_vsv_sub,
adenoma_dsv_sub,#adenoma_vsv_sub,
IBD_dsv_sub,IBD_vsv_sub)


## replace pheno value 
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "adenoma"] = "Adenoma"

colnames(sv_prog_meta_pvalue.sig.edge)[2]<-"SV"
sv_prog_adjAbun.sig.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno %>%
  as.character(.) %>%
  .[!duplicated(.)]

sv_prog_adjAbun.sig.sv <- sv_prog_meta_pvalue.sig.edge$SV %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

## circos plot
sv_prog_count<-NULL

for (phe in sv_prog_adjAbun.sig.pheno) {
  #pheno <-"C4"
  sv_prog_df<-sv_prog_meta_pvalue.sig.edge[sv_prog_meta_pvalue.sig.edge$Pheno==phe,]$SV %>%
    str_replace_all(.,"\\:\\d+_\\d+.*", "") %>%
    table %>%
    as.data.frame
  colnames(sv_prog_df)<-c("Species", "Count")
  sv_prog_df<-data.frame(Pheno = rep(phe,nrow(sv_prog_df)), sv_prog_df)
  sv_prog_count<-rbind(sv_prog_count,sv_prog_df)
}

sv_prog_count$Species<-info$Short_name[match(sv_prog_count$Species, info$organism)]
sv_prog_count <- sv_prog_count[order(sv_prog_count$Count),]

sv_prog_count_species_order<-sv_prog_count %>% group_by(Species) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]
sv_prog_count_pheno_order<-sv_prog_count %>% group_by(Pheno) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]

sv_prog_count<-sv_prog_count[order(match(sv_prog_count$Pheno, sv_prog_count_pheno_order$Pheno),decreasing = F),]

sv_prog_count_species_order_str<-as.character(sv_prog_count_species_order$Species)
sv_prog_count_pheno_order_str<-as.character(sv_prog_count_pheno_order$Pheno)


#pdf("pics/F6_SV_pheno_adjAbun.circos.pdf", width = 26, height = 26)
tiff(file = "pics/F6_SV_index_gas_circos_font3.tiff", width =7000, height =5800, res =300) 

circos.clear()
circos.par(start.degree = 180, "clock.wise" = T) # ,xaxis.clock.wise = F

chordDiagram(sv_prog_count,annotationTrack = "grid",
             #grid.col =  c(wes_palette("Darjeeling1", length(sv_prog_count_pheno_order_str), type = "continuous"),
             #              rep('grey',length(sv_prog_count_species_order_str))),
             grid.col =  c(c("#207F4C","#EE4863","#2983BB"),
                           rep('grey',length(sv_prog_count_species_order_str))),
             order = c(rev(sv_prog_count_pheno_order_str),
                       sv_prog_count_species_order_str),
             big.gap = 6,
             preAllocateTracks = list(track.margin = c(0, uh(45, "mm")), 
             #                         #track.height = max(strwidth(unlist(dimnames(sv_prog_count)))))
                                       track.height=0.1)
             )

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,cex = 2.5,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),font=3)
}, bg.border = NA) # here set bg.border to NA is important

dev.off()

#mycolor6
#hist(mtcars$mpg,
#     breaks = 15,
#     col = "#EE4C97FF",
#)


### "#9fa39a","#621624"
### "#ECD452","#BB97C5","#EA8958","#2983BB","#BC3523","#207F4C","#EE4863","#97C24E","#98369E","#63BBD0","#FEBA07","#EA517F"

write.csv(sv_prog_count,"07.Microbial_GWAS/gas_disorder_sv_prog_count.csv",row.names=F)
unique(sv_prog_count$Species)
unique(sv_prog_meta_pvalue.sig.edge$SV)



################################################################################
###  Metabolic diseases

## CRC  IBD  T2D  ACVD severe obesity MetS  CAD
## fasting_insulin hba1c fasting_glucose
## triglycerides hdl ldl  cholesterol dyastolic_p systolic_p
## hscrp creatine albumine creatinine ast alt total_bilirubin

sv_prog_meta_pvalue.sig.edge<-rbind(
#T2D_vsv_sub, T2D_dsv_sub,
#MetS_dsv_sub,MetS_vsv_sub,
#severe_obesity_dsv_sub,severe_obesity_vsv_sub,
#CAD_dsv_sub,CAD_vsv_sub,
ACVD_dsv_sub,ACVD_vsv_sub,
hba1c_dsv_sub,hba1c_vsv_sub,
fasting_glucose_dsv_sub, #fasting_glucose_vsv_sub,
#fasting_insulin_dsv_sub,fasting_insulin_vsv_sub,
triglycerides_dsv_sub,triglycerides_vsv_sub,
hdl_dsv_sub,hdl_vsv_sub,
ldl_dsv_sub,ldl_vsv_sub,
cholesterol_dsv_sub,cholesterol_vsv_sub,
mean_press_dsv_sub,mean_press_vsv_sub,
cirrhosis_dsv_sub,cirrhosis_vsv_sub,
creatinine_dsv_sub,creatinine_vsv_sub,
#ast_dsv_sub,ast_vsv_sub,alt_dsv_sub,
alt_vsv_sub
#total_bilirubin_vsv_sub,total_bilirubin_dsv_sub
)

## replace pheno value 
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "severe_obesity"] = "Severe obesity"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "hba1c"] = "HbA1c"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "fasting_glucose"] = "Fasting glucose"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "fasting_insulin"] = "Fasting insulin"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "triglycerides"] = "Triglycerides"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "hdl"] = "HDL cholesterol"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "ldl"] = "LDL cholesterol"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "cholesterol"] = "Total cholesterol"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "mean_press"] = "Mean arterial press"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "cirrhosis"] = "Cirrhosis"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "creatinine"] = "Creatinine"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "ast"] = "Ast"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "alt"] = "Alt"
sv_prog_meta_pvalue.sig.edge[,1][sv_prog_meta_pvalue.sig.edge[,1] == "total_bilirubin"] = "Total bilirubin"


colnames(sv_prog_meta_pvalue.sig.edge)[2]<-"SV"

sv_prog_adjAbun.sig.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno %>%
  as.character(.) %>%
  .[!duplicated(.)]

sv_prog_adjAbun.sig.sv <- sv_prog_meta_pvalue.sig.edge$SV %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

## circos plot
sv_prog_count<-NULL

for (phe in sv_prog_adjAbun.sig.pheno) {
  #pheno <-"C4"
  sv_prog_df<-sv_prog_meta_pvalue.sig.edge[sv_prog_meta_pvalue.sig.edge$Pheno==phe,]$SV %>%
    str_replace_all(.,"\\:\\d+_\\d+.*", "") %>%
    table %>%
    as.data.frame
  colnames(sv_prog_df)<-c("Species", "Count")
  sv_prog_df<-data.frame(Pheno = rep(phe,nrow(sv_prog_df)), sv_prog_df)
  sv_prog_count<-rbind(sv_prog_count,sv_prog_df)
}

sv_prog_count$Species<-info$Short_name[match(sv_prog_count$Species, info$organism)]
sv_prog_count <- sv_prog_count[order(sv_prog_count$Count),]

sv_prog_count_species_order<-sv_prog_count %>% group_by(Species) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]
sv_prog_count_pheno_order<-sv_prog_count %>% group_by(Pheno) %>% summarise(sum(Count)) %>% .[order(.$`sum(Count)`, decreasing = T),]

sv_prog_count<-sv_prog_count[order(match(sv_prog_count$Pheno, sv_prog_count_pheno_order$Pheno),decreasing = F),]

sv_prog_count_species_order_str<-as.character(sv_prog_count_species_order$Species)
sv_prog_count_pheno_order_str<-as.character(sv_prog_count_pheno_order$Pheno)


#pdf("pics/F6_SV_pheno_adjAbun.circos.pdf", width = 26, height = 26)
tiff(file = "pics/F6_SV_index_met_circos_font3.tiff", width =7000, height =5800, res =300) 

circos.clear()
circos.par(start.degree = 180, "clock.wise" = T) # ,xaxis.clock.wise = F

chordDiagram(sv_prog_count,annotationTrack = "grid",
             #grid.col =  c(wes_palette("Darjeeling1", length(sv_prog_count_pheno_order_str), type = "continuous"),
             #              rep('grey',length(sv_prog_count_species_order_str))),
             grid.col =  c(c("#E18727FF","#EA517F","#9FA39A","#2983BB","#BC3523","#207F4C","#EE4863",
                            "#97C24E", "#63BBD0","#FEBA07", "#98369E"),
                           rep('grey',length(sv_prog_count_species_order_str))),
             order = c(rev(sv_prog_count_pheno_order_str),
                       sv_prog_count_species_order_str),
             big.gap = 6,
             preAllocateTracks = list(track.margin = c(0, uh(45, "mm")), 
             #                         #track.height = max(strwidth(unlist(dimnames(sv_prog_count)))))
                                       track.height=0.1)
             )

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,cex = 2.5,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),font=3)
}, bg.border = NA) # here set bg.border to NA is important

dev.off()


#mycolor6
#hist(mtcars$mpg,
#     breaks = 15,
#     col ="#9fa39a",
#)

#c("#E18727FF","#EA517F","#9FA39A","#2983BB","#BC3523","#207F4C","#EE4863",
#                            "#97C24E", "#63BBD0","#FEBA07", "#98369E","#BC3C29FF","#0072B5FF","#9FA39A","#621624")


write.csv(sv_prog_count,"07.Microbial_GWAS/Met_sv_prog_count.csv",row.names=F)
unique(sv_prog_count$Species)
unique(sv_prog_meta_pvalue.sig.edge$SV)



################################################################################################################################
##################### 3.4 Examples
##################### 3.4.1 Heatmap of certain species

##############################################################################
## BMI and age

#############################################################################
### B.wexlerae

sig_p_cutoff<-0.001
het_p_cutoff<-0.05
p_count_cutoff<-2

BMI_vsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/BMI_healthy_vsv.csv",header=T)
BMI_vsv_meta_pvalue.anno<-left_join(BMI_vsv_meta_pvalue, vsv_info, by = c("X" = "SV_Name"))

age_vsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/age_healthy_vsv.csv",header=T)
age_vsv_meta_pvalue.anno<-left_join(age_vsv_meta_pvalue, vsv_info, by = c("X" = "SV_Name"))

BMI_vsv_sig_spe<-subset(BMI_vsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff,drop=T)$X
age_vsv_sig_spe<-subset(age_vsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff,drop=T)$X

vsv_select<-data.frame(c(BMI_vsv_sig_spe,age_vsv_sig_spe))
colnames(vsv_select)<-"X"

#######################################
### dsv

BMI_dsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/BMI_healthy_dsv.csv",header=T)
BMI_dsv_meta_pvalue.anno<-left_join(BMI_dsv_meta_pvalue, dsv_info, by = c("X" = "SV_Name"))

age_dsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/age_healthy_dsv.csv",header=T)
age_dsv_meta_pvalue.anno<-left_join(age_dsv_meta_pvalue, dsv_info, by = c("X" = "SV_Name"))

BMI_dsv_sig_spe<-subset(BMI_dsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff,drop=T)$X
age_dsv_sig_spe<-subset(age_dsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff,drop=T)$X

dsv_select<-data.frame(c(BMI_dsv_sig_spe,age_dsv_sig_spe))
colnames(dsv_select)<-"X"

BMI_vsv_select<-merge(vsv_select,BMI_vsv_meta_pvalue.anno,by.x="X",by.y="X")
age_vsv_select<-merge(vsv_select,age_vsv_meta_pvalue.anno,by.x="X",by.y="X")

BMI_dsv_select<-merge(dsv_select,BMI_dsv_meta_pvalue.anno,by.x="X",by.y="X")
age_dsv_select<-merge(dsv_select,age_dsv_meta_pvalue.anno,by.x="X",by.y="X")

BMI_vsv_select$Pheno<-"BMI"
age_vsv_select$Pheno<-"age"

BMI_dsv_select$Pheno<-"BMI"
age_dsv_select$Pheno<-"age"

BMI_vsv_sub<-BMI_vsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]
age_vsv_sub<-age_vsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]

BMI_dsv_sub<-BMI_dsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]
age_dsv_sub<-age_dsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]

BMI_vsv_sub$SV<-c(paste("vsv:",BMI_vsv_sub$X,sep = ""))
age_vsv_sub$SV<-c(paste("vsv:",age_vsv_sub$X,sep = ""))

BMI_dsv_sub$SV<-c(paste("dsv:",BMI_dsv_sub$X,sep = ""))
age_dsv_sub$SV<-c(paste("dsv:",age_dsv_sub$X,sep = ""))

sv_prog_meta_pvalue.sig.edge<-rbind(BMI_vsv_sub,age_vsv_sub,BMI_dsv_sub,age_dsv_sub)

#sv_prog_meta_pvalue.sig.edge$beta<-log(sv_prog_meta_pvalue.sig.edge$effect_meta)
sv_prog_meta_pvalue.sig.edge$beta<-sv_prog_meta_pvalue.sig.edge$effect_meta
sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)

sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)

sv_prog_meta_pvalue.sig.edge$Sign<-NA

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.05 & sv_prog_meta_pvalue.sig.edge$p_count>=p_count_cutoff &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "*" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.01 & sv_prog_meta_pvalue.sig.edge$p_count>=p_count_cutoff &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "¡î" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.001 & sv_prog_meta_pvalue.sig.edge$p_count>=p_count_cutoff &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "¡ï" 

sv_prog_meta_pvalue.sig.edge<-unique(sv_prog_meta_pvalue.sig.edge)

sv_prog_adjAbun.r<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "beta")] %>% spread("Pheno", "beta") %>% data.frame(row.names = "SV")
sv_prog_adjAbun.sign<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "Sign")] %>% spread("Pheno", "Sign") %>% data.frame(row.names = "SV")

sv_prog_meta_pvalue.sig.edge<-subset(sv_prog_meta_pvalue.sig.edge,is.na(pvalue_meta)==F,drop=T)
# write.csv(sv_prog_meta_pvalue.sig.edge,"test.csv")


##########################################################################################################
### B.wexlerae    Blautia wexlerae DSM 19850

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("wexlerae", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("wexlerae", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("wexlerae", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("wexlerae", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Blautia wexlerae DSM 19850:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-3, 0, length.out=ceiling(i/2) + 1), 
              seq(3/i, 3, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

tiff(file = "pics/F6B_BMI_age_B.wexlerae_sv_heatmap.tiff", width =3000, height =5000, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("BMI","Age"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =2.2,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
         )

dev.off()



#############################################################################
### E.rectale   [Eubacterium] rectale DSM 17629  

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("rectale", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("rectale", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("rectale", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("rectale", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot) <- gsub("\\[.*\\]", "", rownames(sv_adjAbun.r.plot))
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "rectale DSM 17629:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-2, 0, length.out=ceiling(i/2) + 1), 
              seq(2/i, 2, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

tiff(file = "pics/F6B_BMI_age_E.rectale_sv_heatmap.tiff", width =3000, height =2000, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("BMI","Age"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =2.2,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
         )

dev.off()


##########################################################################################################
### R.torques   Ruminococcus torques L2-14

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("torques", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("torques", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("torques", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("wexlerae", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Ruminococcus torques L2-14:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-2.5, 0, length.out=ceiling(i/2) + 1), 
              seq(2.5/i, 2.5, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

tiff(file = "pics/F6B_BMI_age_R.torques_sv_heatmap.tiff", width =3000, height =3000, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("BMI","age","Triglycerides","Hdl cholesterol","Total cholesterol","ACVD"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =2.2,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
         )

dev.off()


##########################################################################################################
### C.comes  Coprococcus comes ATCC 27758

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("comes", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("comes", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("comes", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("comes", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Coprococcus comes ATCC 27758:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-3.5, 0, length.out=ceiling(i/2) + 1), 
              seq(3.5/i,3.5, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

tiff(file = "pics/F6B_age_BMI_C.comes_sv_heatmap.tiff", width =3000, height =4200, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("BMI","Age","HbA1c","Ldl cholesterol","ACVD"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =2.2,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
         )

dev.off()




##########################################################################################################
### F.prausnitzii  Faecalibacterium cf. prausnitzii KLE1255

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("prausnitzii", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("prausnitzii", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("prausnitzii", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("prausnitzii", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Faecalibacterium cf. prausnitzii KLE1255", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-3.5, 0, length.out=ceiling(i/2) + 1), 
              seq(3.5/i,3.5, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

tiff(file = "pics/F6B_age_BMI_F.prausnitzii_sv_heatmap.tiff", width =3800, height =2300, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("BMI","Age"),
          cexCol = 2.4, srtCol = 45, cexRow = 2.4,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =3,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
         )

dev.off()




#############################################################################
### E.rectale  (gas disorder)

sig_p_cutoff<-0.001
het_p_cutoff<-0.05
p_count_cutoff<-2

IBD_vsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/IBD_vsv.csv",header=T)
IBD_vsv_meta_pvalue.anno<-left_join(IBD_vsv_meta_pvalue, vsv_info, by = c("X" = "SV_Name"))

CRC_vsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/CRC_vsv.csv",header=T)
CRC_vsv_meta_pvalue.anno<-left_join(CRC_vsv_meta_pvalue, vsv_info, by = c("X" = "SV_Name"))

adenoma_vsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/adenoma_vsv.csv",header=T)
adenoma_vsv_meta_pvalue.anno<-left_join(adenoma_vsv_meta_pvalue, vsv_info, by = c("X" = "SV_Name"))

IBD_vsv_sig_spe<-subset(IBD_vsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff,drop=T)$X
CRC_vsv_sig_spe<-subset(CRC_vsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff,drop=T)$X
adenoma_vsv_sig_spe<-subset(adenoma_vsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff,drop=T)$X

vsv_select<-data.frame(c(IBD_vsv_sig_spe,CRC_vsv_sig_spe,adenoma_vsv_sig_spe))
colnames(vsv_select)<-"X"


#######################################
### dsv

IBD_dsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/IBD_dsv.csv",header=T)
IBD_dsv_meta_pvalue.anno<-left_join(IBD_dsv_meta_pvalue, dsv_info, by = c("X" = "SV_Name"))

CRC_dsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/CRC_dsv.csv",header=T)
CRC_dsv_meta_pvalue.anno<-left_join(CRC_dsv_meta_pvalue, dsv_info, by = c("X" = "SV_Name"))

adenoma_dsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/adenoma_dsv.csv",header=T)
adenoma_dsv_meta_pvalue.anno<-left_join(adenoma_dsv_meta_pvalue, dsv_info, by = c("X" = "SV_Name"))

IBD_dsv_sig_spe<-subset(IBD_dsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff,drop=T)$X
CRC_dsv_sig_spe<-subset(CRC_dsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff,drop=T)$X
adenoma_dsv_sig_spe<-subset(adenoma_dsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff,drop=T)$X

dsv_select<-data.frame(c(IBD_dsv_sig_spe,CRC_dsv_sig_spe,adenoma_dsv_sig_spe))
colnames(dsv_select)<-"X"

IBD_vsv_select<-merge(vsv_select,IBD_vsv_meta_pvalue.anno,by.x="X",by.y="X")
CRC_vsv_select<-merge(vsv_select,CRC_vsv_meta_pvalue.anno,by.x="X",by.y="X")
adenoma_vsv_select<-merge(vsv_select,adenoma_vsv_meta_pvalue.anno,by.x="X",by.y="X")

IBD_dsv_select<-merge(dsv_select,IBD_dsv_meta_pvalue.anno,by.x="X",by.y="X")
CRC_dsv_select<-merge(dsv_select,CRC_dsv_meta_pvalue.anno,by.x="X",by.y="X")
adenoma_dsv_select<-merge(dsv_select,adenoma_dsv_meta_pvalue.anno,by.x="X",by.y="X")

IBD_vsv_select$Pheno<-"IBD"
CRC_vsv_select$Pheno<-"CRC"
adenoma_vsv_select$Pheno<-"Adenoma"

IBD_dsv_select$Pheno<-"IBD"
CRC_dsv_select$Pheno<-"CRC"
adenoma_dsv_select$Pheno<-"Adenoma"

IBD_vsv_select$effect_meta<-log(IBD_vsv_select$OR_meta)
IBD_dsv_select$effect_meta<-log(IBD_dsv_select$OR_meta)

CRC_vsv_select$effect_meta<-log(CRC_vsv_select$OR_meta)
CRC_dsv_select$effect_meta<-log(CRC_dsv_select$OR_meta)

adenoma_vsv_select$effect_meta<-log(adenoma_vsv_select$OR_meta)
adenoma_dsv_select$effect_meta<-log(adenoma_dsv_select$OR_meta)

IBD_vsv_sub<-IBD_vsv_select[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","p_count")]
CRC_vsv_sub<-CRC_vsv_select[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","p_count")]
adenoma_vsv_sub<-adenoma_vsv_select[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","p_count")]

IBD_dsv_sub<-IBD_dsv_select[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","p_count")]
CRC_dsv_sub<-CRC_dsv_select[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","p_count")]
adenoma_dsv_sub<-adenoma_dsv_select[,c("Pheno","X","pvalue_meta","effect_meta","se_meta","p_count")]

IBD_vsv_sub$SV<-c(paste("vsv:",IBD_vsv_sub$X,sep = ""))
CRC_vsv_sub$SV<-c(paste("vsv:",CRC_vsv_sub$X,sep = ""))
adenoma_vsv_sub$SV<-c(paste("vsv:",adenoma_vsv_sub$X,sep = ""))

IBD_dsv_sub$SV<-c(paste("dsv:",IBD_dsv_sub$X,sep = ""))
CRC_dsv_sub$SV<-c(paste("dsv:",CRC_dsv_sub$X,sep = ""))
adenoma_dsv_sub$SV<-c(paste("vsv:",adenoma_dsv_sub$X,sep = ""))


sv_prog_meta_pvalue.sig.edge<-rbind(IBD_vsv_sub,CRC_vsv_sub,IBD_dsv_sub,CRC_dsv_sub)

#sv_prog_meta_pvalue.sig.edge$beta<-log(sv_prog_meta_pvalue.sig.edge$effect_meta)
sv_prog_meta_pvalue.sig.edge$beta<-sv_prog_meta_pvalue.sig.edge$effect_meta
sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)

sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)

sv_prog_meta_pvalue.sig.edge$Sign<-NA

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.05 & sv_prog_meta_pvalue.sig.edge$p_count>=p_count_cutoff &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "*" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.01 & sv_prog_meta_pvalue.sig.edge$p_count>=p_count_cutoff &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "¡î" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.001 & sv_prog_meta_pvalue.sig.edge$p_count>=p_count_cutoff &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "¡ï" 

sv_prog_meta_pvalue.sig.edge<-unique(sv_prog_meta_pvalue.sig.edge)

sv_prog_adjAbun.r<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "beta")] %>% spread("Pheno", "beta") %>% data.frame(row.names = "SV")
sv_prog_adjAbun.sign<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "Sign")] %>% spread("Pheno", "Sign") %>% data.frame(row.names = "SV")

sv_prog_meta_pvalue.sig.edge<-subset(sv_prog_meta_pvalue.sig.edge,is.na(pvalue_meta)==F,drop=T)
# write.csv(sv_prog_meta_pvalue.sig.edge,"test.csv")


##########################################################################################################
###   E.rectale [Eubacterium] rectale DSM 17629

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("rectale", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("rectale", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("rectale", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("rectale", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot) <- gsub("\\[.*\\]", "", rownames(sv_adjAbun.r.plot))
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "rectale DSM 17629:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-1.5, 0, length.out=ceiling(i/2) + 1), 
              seq(1.5/i, 1.5, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

tiff(file = "pics/F6B_E.rectale_gas_disorder.heatmap.tiff", width =3000, height =2000, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("IBD","CRC"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =3.5,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
         )

dev.off()




#########################################################################################
### for metablism disorder 


############################################################################
###  R.torques   Ruminococcus torques L2-14
## HDL cholesterol / ACVD  / Total cholesterol  Alt

sig_p_cutoff<-0.001
het_p_cutoff<-0.05
p_count_cutoff<-2


hdl_vsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/hdl_vsv.csv",header=T)
hdl_vsv_meta_pvalue.anno<-left_join(hdl_vsv_meta_pvalue, vsv_info, by = c("X" = "SV_Name"))

cholesterol_vsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/cholesterol_vsv.csv",header=T)
cholesterol_vsv_meta_pvalue.anno<-left_join(cholesterol_vsv_meta_pvalue, vsv_info, by = c("X" = "SV_Name"))

ACVD_vsv_meta_pvalue<-read.csv("07.Microbial_GWAS/studies_noadjAbun/ACVD_vsv.csv",header=T)
ACVD_vsv_meta_pvalue.anno<-left_join(ACVD_vsv_meta_pvalue, vsv_info, by = c("X" = "SV_Name"))

alt_vsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/alt_vsv.csv",header=T)
alt_vsv_meta_pvalue.anno<-left_join(alt_vsv_meta_pvalue, vsv_info, by = c("X" = "SV_Name"))

hdl_vsv_sig_spe<-subset(hdl_vsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count >= p_count_cutoff,drop=T)$X
cholesterol_vsv_sig_spe<-subset(cholesterol_vsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count >= p_count_cutoff,drop=T)$X
ACVD_vsv_sig_spe<-subset(ACVD_vsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff,drop=T)$X
alt_vsv_sig_spe<-subset(alt_vsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count >= p_count_cutoff,drop=T)$X

vsv_select<-data.frame(c(hdl_vsv_sig_spe,cholesterol_vsv_sig_spe,ACVD_vsv_sig_spe,alt_vsv_sig_spe))
colnames(vsv_select)<-"X"


#######################################
### dsv

hdl_dsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/hdl_dsv.csv",header=T)
hdl_dsv_meta_pvalue.anno<-left_join(hdl_dsv_meta_pvalue, dsv_info, by = c("X" = "SV_Name"))

cholesterol_dsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/cholesterol_dsv.csv",header=T)
cholesterol_dsv_meta_pvalue.anno<-left_join(cholesterol_dsv_meta_pvalue, dsv_info, by = c("X" = "SV_Name"))

ACVD_dsv_meta_pvalue<-read.csv("07.Microbial_GWAS/studies_noadjAbun/ACVD_dsv.csv",header=T)
ACVD_dsv_meta_pvalue.anno<-left_join(ACVD_dsv_meta_pvalue, dsv_info, by = c("X" = "SV_Name"))

alt_dsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/alt_dsv.csv",header=T)
alt_dsv_meta_pvalue.anno<-left_join(alt_dsv_meta_pvalue, dsv_info, by = c("X" = "SV_Name"))

hdl_dsv_sig_spe<-subset(hdl_dsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count >= p_count_cutoff,drop=T)$X
cholesterol_dsv_sig_spe<-subset(cholesterol_dsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count >= p_count_cutoff,drop=T)$X
ACVD_dsv_sig_spe<-subset(ACVD_dsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff,drop=T)$X
alt_dsv_sig_spe<-subset(alt_dsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count >= p_count_cutoff,drop=T)$X

dsv_select<-data.frame(c(hdl_dsv_sig_spe,cholesterol_dsv_sig_spe,ACVD_dsv_sig_spe,alt_dsv_sig_spe))
colnames(dsv_select)<-"X"

alt_vsv_select<-merge(vsv_select,alt_vsv_meta_pvalue.anno,by.x="X",by.y="X")
hdl_vsv_select<-merge(vsv_select,hdl_vsv_meta_pvalue.anno,by.x="X",by.y="X")
cholesterol_vsv_select<-merge(vsv_select,cholesterol_vsv_meta_pvalue.anno,by.x="X",by.y="X")
ACVD_vsv_select<-merge(vsv_select,ACVD_vsv_meta_pvalue.anno,by.x="X",by.y="X")

alt_dsv_select<-merge(dsv_select,alt_dsv_meta_pvalue.anno,by.x="X",by.y="X")
hdl_dsv_select<-merge(dsv_select,hdl_dsv_meta_pvalue.anno,by.x="X",by.y="X")
cholesterol_dsv_select<-merge(dsv_select,cholesterol_dsv_meta_pvalue.anno,by.x="X",by.y="X")
ACVD_dsv_select<-merge(dsv_select,ACVD_dsv_meta_pvalue.anno,by.x="X",by.y="X")

cholesterol_vsv_select$Pheno<-"cholesterol"
hdl_vsv_select$Pheno<-"hdl"
ACVD_vsv_select$Pheno<-"ACVD"
alt_vsv_select$Pheno<-"alt"

cholesterol_dsv_select$Pheno<-"cholesterol"
hdl_dsv_select$Pheno<-"hdl"
ACVD_dsv_select$Pheno<-"ACVD"
alt_dsv_select$Pheno<-"alt"

ACVD_vsv_select$effect_meta<-log(ACVD_vsv_select$JieZ_2017.OR)
ACVD_dsv_select$effect_meta<-log(ACVD_dsv_select$JieZ_2017.OR)

hdl_vsv_sub<-hdl_vsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]
cholesterol_vsv_sub<-cholesterol_vsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]
ACVD_vsv_sub<-ACVD_vsv_select[,c("Pheno","X","pvalue_meta","effect_meta")]
ACVD_vsv_sub$p_count<-2
alt_vsv_sub<-alt_vsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]

hdl_dsv_sub<-hdl_dsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]
cholesterol_dsv_sub<-cholesterol_dsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]
ACVD_dsv_sub<-ACVD_dsv_select[,c("Pheno","X","pvalue_meta","effect_meta")]
ACVD_dsv_sub$p_count<-2
alt_dsv_sub<-alt_dsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]

hdl_vsv_sub$SV<-c(paste("vsv:",hdl_vsv_sub$X,sep = ""))
cholesterol_vsv_sub$SV<-c(paste("vsv:",cholesterol_vsv_sub$X,sep = ""))
ACVD_vsv_sub$SV<-c(paste("vsv:",ACVD_vsv_sub$X,sep = ""))
alt_vsv_sub$SV<-c(paste("vsv:",alt_vsv_sub$X,sep = ""))

hdl_dsv_sub$SV<-c(paste("dsv:",hdl_dsv_sub$X,sep = ""))
cholesterol_dsv_sub$SV<-c(paste("dsv:",cholesterol_dsv_sub$X,sep = ""))
ACVD_dsv_sub$SV<-c(paste("dsv:",ACVD_dsv_sub$X,sep = ""))
alt_dsv_sub$SV<-c(paste("dsv:",alt_dsv_sub$X,sep = ""))

sv_prog_meta_pvalue.sig.edge<-rbind(hdl_vsv_sub,
cholesterol_vsv_sub,ACVD_vsv_sub,alt_vsv_sub,
hdl_dsv_sub,cholesterol_dsv_sub,ACVD_dsv_sub,alt_dsv_sub)

#sv_prog_meta_pvalue.sig.edge$beta<-log(sv_prog_meta_pvalue.sig.edge$effect_meta)
sv_prog_meta_pvalue.sig.edge$beta<-sv_prog_meta_pvalue.sig.edge$effect_meta
sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)

sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)

sv_prog_meta_pvalue.sig.edge$Sign<-NA

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.05 & sv_prog_meta_pvalue.sig.edge$p_count>=p_count_cutoff &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "*" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.01 & sv_prog_meta_pvalue.sig.edge$p_count>=p_count_cutoff &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "¡î" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.001 & sv_prog_meta_pvalue.sig.edge$p_count>=p_count_cutoff &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "¡ï" 

sv_prog_meta_pvalue.sig.edge<-unique(sv_prog_meta_pvalue.sig.edge)

sv_prog_adjAbun.r<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "beta")] %>% spread("Pheno", "beta") %>% data.frame(row.names = "SV")
sv_prog_adjAbun.sign<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "Sign")] %>% spread("Pheno", "Sign") %>% data.frame(row.names = "SV")

sv_prog_meta_pvalue.sig.edge<-subset(sv_prog_meta_pvalue.sig.edge,is.na(pvalue_meta)==F,drop=T)
# write.csv(sv_prog_meta_pvalue.sig.edge,"test.csv")


##########################################################################################################
### R.torques   Ruminococcus torques L2-14

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("torques", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("torques", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("torques", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("torques", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Ruminococcus torques L2-14:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-3, 0, length.out=ceiling(i/2) + 1), 
              seq(3/i, 3, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

tiff(file = "pics/F6B_metdisorder_R.torques_sv_heatmap.tiff", width =2800, height =2800, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("Hdl cholesterol","Total cholesterol","ACVD","ALT"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =3,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
         )

dev.off()




##########################################################################################
###  P.copri

# ldl
# cholesterol

sig_p_cutoff<-0.001
het_p_cutoff<-0.05
p_count_cutoff<-2


ldl_vsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/ldl_vsv.csv",header=T)
ldl_vsv_meta_pvalue.anno<-left_join(ldl_vsv_meta_pvalue, vsv_info, by = c("X" = "SV_Name"))

cholesterol_vsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/cholesterol_vsv.csv",header=T)
cholesterol_vsv_meta_pvalue.anno<-left_join(cholesterol_vsv_meta_pvalue, vsv_info, by = c("X" = "SV_Name"))

ldl_vsv_sig_spe<-subset(ldl_vsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff ,drop=T)$X
cholesterol_vsv_sig_spe<-subset(cholesterol_vsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff,drop=T)$X

vsv_select<-data.frame(c(ldl_vsv_sig_spe,cholesterol_vsv_sig_spe))
colnames(vsv_select)<-"X"


#######################################
### dsv

ldl_dsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/ldl_dsv.csv",header=T)
ldl_dsv_meta_pvalue.anno<-left_join(ldl_dsv_meta_pvalue, dsv_info, by = c("X" = "SV_Name"))

cholesterol_dsv_meta_pvalue<-read.csv("07.Microbial_GWAS/meta_noadjAbun/cholesterol_dsv.csv",header=T)
cholesterol_dsv_meta_pvalue.anno<-left_join(cholesterol_dsv_meta_pvalue, dsv_info, by = c("X" = "SV_Name"))

ldl_dsv_sig_spe<-subset(ldl_dsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff ,drop=T)$X
cholesterol_dsv_sig_spe<-subset(cholesterol_dsv_meta_pvalue,as.numeric(pvalue_meta) <= sig_p_cutoff & het_p > het_p_cutoff & p_count>=p_count_cutoff,drop=T)$X

dsv_select<-data.frame(c(ldl_dsv_sig_spe,cholesterol_dsv_sig_spe))
colnames(dsv_select)<-"X"

cholesterol_vsv_select<-merge(vsv_select, cholesterol_vsv_meta_pvalue.anno,by.x="X",by.y="X")
ldl_vsv_select<-merge(vsv_select,ldl_vsv_meta_pvalue.anno,by.x="X",by.y="X")

cholesterol_dsv_select<-merge(dsv_select, cholesterol_dsv_meta_pvalue.anno,by.x="X",by.y="X")
ldl_dsv_select<-merge(dsv_select,ldl_dsv_meta_pvalue.anno,by.x="X",by.y="X")

ldl_vsv_select$Pheno<-"ldl"
cholesterol_vsv_select$Pheno<-"cholesterol"

ldl_dsv_select$Pheno<-"ldl"
cholesterol_dsv_select$Pheno<-"cholesterol"

##,"se_meta","cilb","ciub"

ldl_vsv_sub<- ldl_vsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]
cholesterol_vsv_sub<-cholesterol_vsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]

ldl_dsv_sub<-ldl_dsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]
cholesterol_dsv_sub<-cholesterol_dsv_select[,c("Pheno","X","pvalue_meta","effect_meta","p_count")]

ldl_vsv_sub$SV<-c(paste("vsv:", ldl_vsv_sub$X,sep = ""))
cholesterol_vsv_sub$SV<-c(paste("vsv:",cholesterol_vsv_sub$X,sep = ""))

ldl_dsv_sub$SV<-c(paste("dsv:", ldl_dsv_sub$X,sep = ""))
cholesterol_dsv_sub$SV<-c(paste("dsv:",cholesterol_dsv_sub$X,sep = ""))

sv_prog_meta_pvalue.sig.edge<-rbind(cholesterol_vsv_sub,ldl_vsv_sub, cholesterol_dsv_sub,
ldl_dsv_sub)

sv_prog_meta_pvalue.sig.edge$beta<-sv_prog_meta_pvalue.sig.edge$effect_meta
sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)

sv_prog_meta_pvalue.sig.edge<-as.data.frame(sv_prog_meta_pvalue.sig.edge)
sv_prog_meta_pvalue.sig.edge$Sign<-NA

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.05 & sv_prog_meta_pvalue.sig.edge$p_count>=p_count_cutoff &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "*" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.01 & sv_prog_meta_pvalue.sig.edge$p_count>=p_count_cutoff &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "¡î" 

sv_prog_meta_pvalue.sig.edge$Sign[sv_prog_meta_pvalue.sig.edge$pvalue_meta <= 0.001 & sv_prog_meta_pvalue.sig.edge$p_count>=p_count_cutoff &
  !is.na(sv_prog_meta_pvalue.sig.edge$pvalue_meta)] <- "¡ï" 

sv_prog_meta_pvalue.sig.edge<-unique(sv_prog_meta_pvalue.sig.edge)

sv_prog_adjAbun.r<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "beta")] %>% spread("Pheno", "beta") %>% data.frame(row.names = "SV")
sv_prog_adjAbun.sign<- sv_prog_meta_pvalue.sig.edge[,c("Pheno", "SV", "Sign")] %>% spread("Pheno", "Sign") %>% data.frame(row.names = "SV")

sv_prog_meta_pvalue.sig.edge<-subset(sv_prog_meta_pvalue.sig.edge,is.na(pvalue_meta)==F,drop=T)
# write.csv(sv_prog_meta_pvalue.sig.edge,"test.csv")


##########################################################################################################
### Prevotella copri DSM 18205:699_701

plot.pheno  <- sv_prog_meta_pvalue.sig.edge$Pheno[grep("copri", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)]

plot.sv <- sv_prog_meta_pvalue.sig.edge$SV[grep("copri", sv_prog_meta_pvalue.sig.edge$SV)] %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

spe.sv_adjAbun.r<-sv_prog_adjAbun.r[grep("copri", rownames(sv_prog_adjAbun.r)),]
#rownames(spe.sv_adjAbun.r)<-str_replace_all(rownames(spe.sv_adjAbun.r), "^.*sv_", "")

spe.sv_adjAbun.sign<-sv_prog_adjAbun.sign[grep("copri", rownames(sv_prog_adjAbun.sign)),]
#spe.sv_adjAbun.sign<-spe.sv_adjAbun.sign[-15,]
rownames(spe.sv_adjAbun.sign)<-str_replace_all(rownames(spe.sv_adjAbun.sign), "^.*sv_", "")

sv_adjAbun.r.plot    <- spe.sv_adjAbun.r %>%
  .[match(plot.sv,rownames(.)),match(plot.pheno,colnames(.))]

sv_adjAbun.r.plot[is.na(sv_adjAbun.r.plot)]<-0

sv_adjAbun.sign.plot <- spe.sv_adjAbun.sign %>%
  .[match(plot.sv,rownames(spe.sv_adjAbun.r)),
    match(plot.pheno,colnames(spe.sv_adjAbun.r))]

rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "Prevotella copri DSM 18205:", "")
rownames(sv_adjAbun.r.plot)<-str_replace_all(rownames(sv_adjAbun.r.plot), "^.*sv_", "")

i<-256
#myBreaks <- c(seq(min(sv_adjAbun.r.plot), 0, length.out=ceiling(i/2) + 1), 
#              seq(max(sv_adjAbun.r.plot)/i, max(sv_adjAbun.r.plot), length.out=floor(i/2)))

myBreaks <- c(seq(-1, 0, length.out=ceiling(i/2) + 1), 
              seq(1/i, 1, length.out=floor(i/2)))

#sv_adjAbun.r.plot<-as.matrix(sv_adjAbun.r.plot[-15,])

b<-apply(sv_adjAbun.r.plot,2,as.numeric)
rownames(b)<-rownames(sv_adjAbun.r.plot)

tiff(file = "pics/F6B_metdisorder_P.copri_sv_heatmap.tiff", width =3000, height =2300, res =300) 

heatmap.2(b,
          col=colorRampPalette(c("#DE3E69","#F5F5F5","#708090"))(i),
          breaks = myBreaks,
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",labCol = c("cholesterol","Fasting glucose"),
          cexCol = 2, srtCol = 45, cexRow = 2,
          cellnote = sv_adjAbun.sign.plot, notecol = "black",notecex =2.2,
          key.xlab = "Beta coefficient",
          key.par=list(mar=c(4,6,4,4), cex.axis = 2, cex.lab = 2), 
          lm=rbind( c(0, 4, 3), c(2, 1, 0 )), lhei = c(0.5, 1.5),lwid=c(0.5, 4, 1 ),key.title = NA,
          margins=c(12,12) # ("margin.Y", "margin.X")
         )

dev.off()





