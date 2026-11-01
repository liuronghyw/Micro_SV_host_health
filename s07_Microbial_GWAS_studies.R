### Microbial GWA 
### 2026-10-12
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

## 1 Preparation
### 1.1 Import
#install.packages("survminer")
rm(list = ls())
#options(scipen = 200)

source("functions.R", encoding = "UTF-8")
knitr::opts_chunk$set(echo = TRUE)
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
library("metafor")

library(foreach)
library(doParallel)
library(dplyr)
library(magrittr)
library(car)
library("scales")


###################### 1.2 Inputs
info<- read.table("01.cleanData/SV_info_QC/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

all_vsv_info<-read.table("01.cleanData/SV_info_QC/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
all_dsv_info<-read.table("01.cleanData/SV_info_QC/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

if (!all(c("SV_ID","SV_Name","SV_Name_long") %in% names(all_dsv_info)) ||
    !all(c("SV_ID","SV_Name","SV_Name_long") %in% names(all_vsv_info))) {
  stop("QC annotations must provide canonical SV_ID and short/long names")
}

all_basic<-read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv")
if(!"read_count"%in%colnames(all_basic)) stop("read_count is missing from all_basic")
all_basic$read_count<-as.numeric(all_basic$read_count)
#all_basic$hscrp <- as.numeric(replace(all_basic$hscrp, all_basic$hscrp =="<5", NA))

load("01.cleanData/SV_all_QC/dsgv.RData")
load("01.cleanData/SV_all_QC/vsgv.RData")

vsgv<-vsgv_sub
dsgv<-dsgv_sub


# Only genuine QC SV columns are tested; no artificial duplicate is added.
all_abun<-read.table("01.cleanData/mbio_all_QC/SV_species_abun.tsv", check.names = F) 

pc_list<-load_study_pc_list(
  "01.cleanData/SV_all_QC/study_pc",
  unique(all_basic$study_name)
)

## 2 Associations between SVs and prognosis
### 2.1 Preparation
if (!dir.exists("07.Microbial_GWAS")) {dir.create("07.Microbial_GWAS")}
if (!dir.exists("07.Microbial_GWAS/studies_adjAbun")) {dir.create("07.Microbial_GWAS/studies_adjAbun")}
if (!dir.exists("07.Microbial_GWAS/meta_adjAbun")) {dir.create("07.Microbial_GWAS/meta_adjAbun")}
if (!dir.exists("07.Microbial_GWAS/multi_adjAbun")) {dir.create("07.Microbial_GWAS/multi_adjAbun")}
if (!dir.exists("07.Microbial_GWAS/study_wide_sig/adjAbun")) {dir.create("07.Microbial_GWAS/study_wide_sig/adjAbun",recursive=T)}
if (!dir.exists("07.Microbial_GWAS/analysis_wide_sig/adjAbun")) {dir.create("07.Microbial_GWAS/analysis_wide_sig/adjAbun",recursive=T)}
if (!dir.exists("07.Microbial_GWAS/studies_adjAbunPCs")) {dir.create("07.Microbial_GWAS/studies_adjAbunPCs")}
if (!dir.exists("07.Microbial_GWAS/meta_adjAbunPCs")) {dir.create("07.Microbial_GWAS/meta_adjAbunPCs")}
if (!dir.exists("07.Microbial_GWAS/multi_adjAbunPCs")) {dir.create("07.Microbial_GWAS/multi_adjAbunPCs")}
if (!dir.exists("07.Microbial_GWAS/study_wide_sig/adjAbunPCs")) {dir.create("07.Microbial_GWAS/study_wide_sig/adjAbunPCs",recursive=T)}
if (!dir.exists("07.Microbial_GWAS/analysis_wide_sig/adjAbunPCs")) {dir.create("07.Microbial_GWAS/analysis_wide_sig/adjAbunPCs",recursive=T)}



# all
all_abun$Others<-1-rowSums(all_abun)

#all_abun_clr<-abundances(x=as.data.frame(na.omit(all_abun+0.02)), transform="clr") %>%as.data.frame

all_abun_clr<-abundances(x=as.data.frame(na.omit(all_abun)), transform="clr") %>%as.data.frame
all_abun_clr <- all_abun_clr[match(rownames(all_abun), rownames(all_abun_clr)),]
rownames(all_abun_clr) <- rownames(all_abun)
#all_covar<-cbind(all_basic,all_abun_clr)

all_covar <- cbind(all_basic, all_abun_clr)
#all_covar<-cbind(all_basic)

# Use the CLR-transformed abundance for all abundance-adjusted models.


############################################################################################################
## associations between disease and SVs

disease_asso<-function(basic,dis_name,adjust_cov,out_vsv_adjust,out_dsv_adjust,out_vsv_adjust_pcs,out_dsv_adjust_pcs){
   #disease_asso<-function(basic,dis_name,adjust_cov,out_vsv_adjust,out_dsv_adjust,out_vsv_noadjust,out_dsv_noadjust){
   #basic<-basic_CRC
   #dis_name<-"CRC"
   #adjust_cov<-c("age","gender_code","BMI","read_count")
   #prog<-basic[,c(dis_name,dis_name)]
   adjust_cov<-unique(c(adjust_cov,"read_count"))

   prog<-data.frame(basic[,dis_name])
   row.names(prog)<-rownames(basic)
   colnames(prog)<-"Y"

   abun_clr<-all_abun_clr[rownames(basic),]

   dis_vsgv<-vsgv[rownames(basic),]
   dis_dsgv<-dsgv[rownames(basic),]

   vsv_adjAbun<-logistic_adjAbun_vsv(prog,dis_vsgv,basic,adjust_cov,abun_clr,info,9,dis_name)
   dsv_adjAbun<-logistic_adjAbun_dsv(prog,dis_dsgv,basic,adjust_cov,abun_clr,info,9,dis_name)

   vsv_adjAbunPCs<-logistic_adjAbunPCs_vsv(prog,dis_vsgv,basic,adjust_cov,abun_clr,info,9,dis_name,pc_list)
   dsv_adjAbunPCs<-logistic_adjAbunPCs_dsv(prog,dis_dsgv,basic,adjust_cov,abun_clr,info,9,dis_name,pc_list)

   #vsv_noadjAbun<-logistic_noadjAbun_vsv(prog,dis_vsgv,basic,adjust_cov,info,dis_name)
   #dsv_noadjAbun<-logistic_noadjAbun_dsv(prog,dis_dsgv,basic,adjust_cov,info,dis_name)

   write.csv(vsv_adjAbun, file = out_vsv_adjust,row.names=F)
   write.csv(dsv_adjAbun, file = out_dsv_adjust,row.names=F)

   write.csv(vsv_adjAbunPCs, file = out_vsv_adjust_pcs,row.names=F)
   write.csv(dsv_adjAbunPCs, file = out_dsv_adjust_pcs,row.names=F)

   #write.csv(vsv_noadjAbun, file = out_vsv_noadjust,row.names=F)
   #write.csv(dsv_noadjAbun, file = out_dsv_noadjust,row.names=F)
}


### diseaase
### IBD T2D CRC CAD
### ACVD adenoma cirrhosis  MetS PD

basic_CRC<-subset(all_basic,is.na(CRC)==F,drop=T)
a<-table(basic_CRC$study_name,basic_CRC$CRC)
basic_CRC<-subset(all_basic,is.na(CRC)==F & study_name!="FengQ_2015",drop=T)

disease_asso(basic_CRC,"CRC",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/CRC_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/CRC_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/CRC_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/CRC_dsv.csv")

basic_T2D<-subset(all_basic,is.na(T2D)==F,drop=T)
a<-table(basic_T2D$study_name,basic_T2D$T2D)
basic_T2D<-subset(all_basic,is.na(T2D)==F & study_name!="HMP_2019_t2d" & study_name!="LiJ_2014",drop=T)
disease_asso(basic_T2D,"T2D",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/T2D_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/T2D_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/T2D_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/T2D_dsv.csv")

basic_ACVD<-subset(all_basic,is.na(ACVD)==F,drop=T)
table(basic_ACVD$study_name,basic_ACVD$ACVD)
disease_asso(basic_ACVD,"ACVD",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/ACVD_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/ACVD_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/ACVD_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/ACVD_dsv.csv")


basic_CAD<-subset(all_basic,is.na(CAD)==F,drop=T)
table(basic_CAD$study_name,basic_CAD$CAD)
disease_asso(basic_CAD,"CAD",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/CAD_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/CAD_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/CAD_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/CAD_dsv.csv")


basic_cirrhosis<-subset(all_basic,is.na(cirrhosis)==F,drop=T)
table(basic_cirrhosis$study_name,basic_cirrhosis$cirrhosis)
disease_asso(basic_cirrhosis,"cirrhosis",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/cirrhosis_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/cirrhosis_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/cirrhosis_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/cirrhosis_dsv.csv")

basic_MetS<-subset(all_basic,is.na(MetS)==F,drop=T)
table(basic_MetS$study_name,basic_MetS$MetS)
disease_asso(basic_MetS,"MetS",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/MetS_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/MetS_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/MetS_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/MetS_dsv.csv")

basic_adenoma<-subset(all_basic,is.na(adenoma)==F,drop=T)
table(basic_adenoma$study_name,basic_adenoma$adenoma)
basic_adenoma<-subset(all_basic,is.na(adenoma)==F & study_name!="FengQ_2015",drop=T)
disease_asso(basic_adenoma,"adenoma",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/adenoma_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/adenoma_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/adenoma_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/adenoma_dsv.csv")

basic_PD<-subset(all_basic,is.na(PD)==F,drop=T)
table(basic_PD$study_name,basic_PD$PD)
disease_asso(basic_PD,"PD",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/PD_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/PD_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/PD_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/PD_dsv.csv")

basic_IBD<-subset(all_basic,is.na(IBD)==F,drop=T)
table(basic_IBD$study_name,basic_IBD$IBD)
basic_IBD<-subset(all_basic,is.na(IBD)==F & study_name!="LiJ_2014",drop=T)

disease_asso(basic_IBD,"IBD",c("read_count"),
"07.Microbial_GWAS/studies_adjAbun/IBD_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/IBD_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/IBD_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/IBD_dsv.csv")


### 2.3 linear regression model 1
### Linear model with covariates: age, gender, BMI, read count.
### 1 indicate the organism name column of the info file

############################################################################################
#######  index 
## age  BMI fasting_insulin  hba1c  fasting_glucose 
## triglycerides hdl  ldl  cholesterol 
## dyastolic_p systolic_p prothrombin_time inr
## albumine creatinine  alt  ast  total_bilirubin

############# age all
age_basic<-subset(all_basic,is.na(age)==F ,drop=T)
table(age_basic$study_name)

age_basic<-subset(age_basic,study_name!="LiJ_2014",drop=T)

table(age_basic$study_name)

age_prog<-data.frame(age_basic[,"age_norm"])
row.names(age_prog)<-rownames(age_basic)
colnames(age_prog)<-"Y"

age_abun_clr<-all_abun_clr[rownames(age_basic),]

age_vsgv<-vsgv[rownames(age_basic),]
age_dsgv<-dsgv[rownames(age_basic),]

vsv_age_adjAbun<-lm_mats_adjAbun_vsv(age_prog,age_vsgv,age_basic,c("gender_code","read_count"),age_abun_clr,info,9,"age")
dsv_age_adjAbun<-lm_mats_adjAbun_dsv(age_prog,age_dsgv,age_basic,c("gender_code","read_count"),age_abun_clr,info,9,"age")
vsv_age_adjAbunPCs<-lm_mats_adjAbun_vsv(age_prog,age_vsgv,age_basic,c("gender_code","read_count"),age_abun_clr,info,9,"age",pc_list)
dsv_age_adjAbunPCs<-lm_mats_adjAbun_dsv(age_prog,age_dsgv,age_basic,c("gender_code","read_count"),age_abun_clr,info,9,"age",pc_list)

write.csv(vsv_age_adjAbun, file = "07.Microbial_GWAS/studies_adjAbun/age_vsv.csv",row.names=F)
write.csv(dsv_age_adjAbun, file = "07.Microbial_GWAS/studies_adjAbun/age_dsv.csv",row.names=F)
write.csv(vsv_age_adjAbunPCs, file = "07.Microbial_GWAS/studies_adjAbunPCs/age_vsv.csv",row.names=F)
write.csv(dsv_age_adjAbunPCs, file = "07.Microbial_GWAS/studies_adjAbunPCs/age_dsv.csv",row.names=F)



############# BMI all
BMI_basic<-subset(all_basic,is.na(BMI)==F,header=T)
table(BMI_basic$study_name)

BMI_basic<-subset(BMI_basic,study_name!="PiccinnoG_2025_cohort6",drop=T)

BMI_prog<-data.frame(BMI_basic[,"BMI_norm"])
row.names(BMI_prog)<-rownames(BMI_basic)
colnames(BMI_prog)<-"Y"

BMI_abun_clr<-all_abun_clr[rownames(BMI_basic),]

BMI_vsgv<-vsgv[rownames(BMI_basic),]
BMI_dsgv<-dsgv[rownames(BMI_basic),]

vsv_BMI_adjAbun<-lm_mats_adjAbun_vsv(BMI_prog,BMI_vsgv,BMI_basic,c("gender_code","read_count"),BMI_abun_clr,info,9,"BMI")
dsv_BMI_adjAbun<-lm_mats_adjAbun_dsv(BMI_prog,BMI_dsgv,BMI_basic,c("gender_code","read_count"),BMI_abun_clr,info,9,"BMI")
vsv_BMI_adjAbunPCs<-lm_mats_adjAbun_vsv(BMI_prog,BMI_vsgv,BMI_basic,c("gender_code","read_count"),BMI_abun_clr,info,9,"BMI",pc_list)
dsv_BMI_adjAbunPCs<-lm_mats_adjAbun_dsv(BMI_prog,BMI_dsgv,BMI_basic,c("gender_code","read_count"),BMI_abun_clr,info,9,"BMI",pc_list)

write.csv(vsv_BMI_adjAbun, file = "07.Microbial_GWAS/studies_adjAbun/BMI_vsv.csv",row.names=F)
write.csv(dsv_BMI_adjAbun, file = "07.Microbial_GWAS/studies_adjAbun/BMI_dsv.csv",row.names=F)
write.csv(vsv_BMI_adjAbunPCs, file = "07.Microbial_GWAS/studies_adjAbunPCs/BMI_vsv.csv",row.names=F)
write.csv(dsv_BMI_adjAbunPCs, file = "07.Microbial_GWAS/studies_adjAbunPCs/BMI_dsv.csv",row.names=F)


index_asso<-function(clin,index,adjust_covar,out_vsv_adjust,out_dsv_adjust,out_vsv_adjust_pcs,out_dsv_adjust_pcs){
   #basic<-basic_CRC
   #dis_name<-"CRC"
   #adjust_cov<-c("age","gender_code","BMI")
   #prog<-basic[,c(dis_name,dis_name)]
   adjust_covar<-unique(c(adjust_covar,"read_count"))

   prog<-data.frame(clin[,index])
   row.names(prog)<-rownames(clin)
   colnames(prog)<-"Y"

   abun_clr<-all_abun_clr[rownames(clin),]

   index_vsgv<-vsgv[rownames(clin),]
   index_dsgv<-dsgv[rownames(clin),]

   vsv_adjAbun<-lm_mats_adjAbun_vsv(prog,index_vsgv,clin,adjust_covar,abun_clr,info,9,index)
   dsv_adjAbun<-lm_mats_adjAbun_dsv(prog,index_dsgv,clin,adjust_covar,abun_clr,info,9,index)

   vsv_adjAbunPCs<-lm_mats_adjAbun_vsv(prog,index_vsgv,clin,adjust_covar,abun_clr,info,9,index,pc_list)
   dsv_adjAbunPCs<-lm_mats_adjAbun_dsv(prog,index_dsgv,clin,adjust_covar,abun_clr,info,9,index,pc_list)

   #vsv_noadjAbun<-lm_mats_noadjAbun_vsv(prog,index_vsgv,clin,adjust_covar,info,index)
   #dsv_noadjAbun<-lm_mats_noadjAbun_dsv(prog,index_dsgv,clin,adjust_covar,info,index)

   write.csv(vsv_adjAbun, file = out_vsv_adjust,row.names=F)
   write.csv(dsv_adjAbun, file = out_dsv_adjust,row.names=F)
   write.csv(vsv_adjAbunPCs, file = out_vsv_adjust_pcs,row.names=F)
   write.csv(dsv_adjAbunPCs, file = out_dsv_adjust_pcs,row.names=F)

   #write.csv(vsv_noadjAbun, file = out_vsv_noadjust,row.names=F)
   #write.csv(dsv_noadjAbun, file = out_dsv_noadjust,row.names=F)
}


basic_hba1c<-subset(all_basic,is.na(hba1c)==F,drop=T)
index_asso(basic_hba1c,"hba1c_norm",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/hba1c_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/hba1c_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/hba1c_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/hba1c_dsv.csv")

basic_triglycerides<-subset(all_basic,is.na(triglycerides)==F,drop=T)
table(basic_triglycerides$study_name)
index_asso(basic_triglycerides,"triglycerides_norm",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/triglycerides_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/triglycerides_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/triglycerides_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/triglycerides_dsv.csv")


basic_fasting_insulin<-subset(all_basic,is.na(fasting_insulin)==F,drop=T)
table(basic_fasting_insulin$study_name)
index_asso(basic_fasting_insulin,"fasting_insulin_norm",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/fasting_insulin_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/fasting_insulin_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/fasting_insulin_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/fasting_insulin_dsv.csv")

basic_hdl<-subset(all_basic,is.na(hdl)==F,drop=T)
table(basic_hdl$study_name)
index_asso(basic_hdl,"hdl_norm",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/hdl_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/hdl_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/hdl_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/hdl_dsv.csv")

basic_ldl<-subset(all_basic,is.na(ldl)==F,drop=T)
table(basic_ldl$study_name)
index_asso(basic_ldl,"ldl_norm",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/ldl_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/ldl_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/ldl_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/ldl_dsv.csv")

basic_cholesterol<-subset(all_basic,is.na(cholesterol)==F,drop=T)
table(basic_cholesterol$study_name)
index_asso(basic_cholesterol,"cholesterol_norm",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/cholesterol_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/cholesterol_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/cholesterol_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/cholesterol_dsv.csv")


basic_ast<-subset(all_basic,is.na(ast)==F,drop=T)
index_asso(basic_ast,"ast_norm",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/ast_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/ast_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/ast_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/ast_dsv.csv")


basic_mean_press<-subset(all_basic,is.na(mean_press)==F,drop=T)
index_asso(basic_mean_press,"mean_press_norm",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/mean_press_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/mean_press_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/mean_press_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/mean_press_dsv.csv")


basic_alt<-subset(all_basic,is.na(alt)==F,drop=T)
index_asso(basic_alt,"alt_norm",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/alt_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/alt_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/alt_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/alt_dsv.csv")


basic_total_bilirubin<-subset(all_basic,is.na(total_bilirubin)==F,drop=T)
index_asso(basic_total_bilirubin,"total_bilirubin_norm",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/total_bilirubin_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/total_bilirubin_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/total_bilirubin_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/total_bilirubin_dsv.csv")


basic_creatinine<-subset(all_basic,is.na(creatinine)==F,drop=T)
index_asso(basic_creatinine,"creatinine_norm",c("age","gender_code","BMI","read_count"),
"07.Microbial_GWAS/studies_adjAbun/creatinine_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/creatinine_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/creatinine_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/creatinine_dsv.csv")

basic_fasting_glucose<-subset(all_basic,is.na(fasting_glucose)==F,drop=T)
index_asso(basic_fasting_glucose,"fasting_glucose_norm",c("read_count"),
"07.Microbial_GWAS/studies_adjAbun/fasting_glucose_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/fasting_glucose_dsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/fasting_glucose_vsv.csv",
"07.Microbial_GWAS/studies_adjAbunPCs/fasting_glucose_dsv.csv")




#################################################################################################################################
### Model 1 meta-analysis: abundance-adjusted

model1_binary_multi<-c("CRC","IBD","T2D","adenoma","PD")

model1_continuous<-c(
  "age","BMI","hba1c","fasting_glucose",
  "fasting_insulin","triglycerides","hdl","ldl","cholesterol","mean_press",
  "creatinine","ast","alt","total_bilirubin"
)

model1_binary_single<-c("ACVD","CAD","cirrhosis","MetS")

for(outcome in model1_binary_multi){
  for(sv_type in c("vsv","dsv")){
    study_file<-paste0("07.Microbial_GWAS/studies_adjAbun/",outcome,"_",sv_type,".csv")
    meta_file<-paste0("07.Microbial_GWAS/meta_adjAbun/",outcome,"_",sv_type,".csv")
    if(file.exists(study_file)){
      combine_meta_p_binary(read.csv(study_file,header=T),6,meta_file,single_study_I2=1)
    }
  }
  meta_vsv<-paste0("07.Microbial_GWAS/meta_adjAbun/",outcome,"_vsv.csv")
  meta_dsv<-paste0("07.Microbial_GWAS/meta_adjAbun/",outcome,"_dsv.csv")
  if(file.exists(meta_vsv) && file.exists(meta_dsv)){
    multi_adjust_sv(
      meta_vsv,
      meta_dsv,
      paste0("07.Microbial_GWAS/multi_adjAbun/",outcome,"_sv.csv")
    )
  }
}


for(outcome in model1_continuous){
  for(sv_type in c("vsv","dsv")){
    study_file<-paste0("07.Microbial_GWAS/studies_adjAbun/",outcome,"_",sv_type,".csv")
    meta_file<-paste0("07.Microbial_GWAS/meta_adjAbun/",outcome,"_",sv_type,".csv")
    if(file.exists(study_file)){
      combine_meta_p_lr(read.csv(study_file,header=T),4,meta_file,single_study_I2=1)
    }
  }
  meta_vsv<-paste0("07.Microbial_GWAS/meta_adjAbun/",outcome,"_vsv.csv")
  meta_dsv<-paste0("07.Microbial_GWAS/meta_adjAbun/",outcome,"_dsv.csv")
  if(file.exists(meta_vsv) && file.exists(meta_dsv)){
    multi_adjust_sv(
      meta_vsv,
      meta_dsv,
      paste0("07.Microbial_GWAS/multi_adjAbun/",outcome,"_sv.csv")
    )
  }
}

for(outcome in model1_binary_single){
  for(sv_type in c("vsv","dsv")){
    study_file<-paste0("07.Microbial_GWAS/studies_adjAbun/",outcome,"_",sv_type,".csv")
    meta_file<-paste0("07.Microbial_GWAS/meta_adjAbun/",outcome,"_",sv_type,".csv")
    if(file.exists(study_file)){
      meta_data<-read.csv(study_file,header=T)
      meta_data$pvalue_meta<-meta_data[[4]]
      meta_data$I2<-rep(1,nrow(meta_data))
      meta_data$het_p<-rep(1,nrow(meta_data))
      write.csv(meta_data,meta_file,row.names=F)
    }
  }
  meta_vsv<-paste0("07.Microbial_GWAS/meta_adjAbun/",outcome,"_vsv.csv")
  meta_dsv<-paste0("07.Microbial_GWAS/meta_adjAbun/",outcome,"_dsv.csv")
  if(file.exists(meta_vsv) && file.exists(meta_dsv)){
    multi_adjust_sv(
      meta_vsv,
      meta_dsv,
      paste0("07.Microbial_GWAS/multi_adjAbun/",outcome,"_sv.csv")
    )
  }
}

###################################################################################################
### Model 2 meta-analysis: abundance and genetic-structure PCs

model2_binary_multi<-c("CRC","IBD","T2D","adenoma","PD")

model2_continuous<-c(
  "age","BMI","hba1c","fasting_glucose",
  "fasting_insulin","triglycerides","hdl","ldl","cholesterol","mean_press",
  "creatinine","ast","alt","total_bilirubin"
)

model2_binary_single<-c("ACVD","CAD","cirrhosis","MetS")

for(outcome in model2_binary_multi){
  for(sv_type in c("vsv","dsv")){
    study_file<-paste0("07.Microbial_GWAS/studies_adjAbunPCs/",outcome,"_",sv_type,".csv")
    meta_file<-paste0("07.Microbial_GWAS/meta_adjAbunPCs/",outcome,"_",sv_type,".csv")
    if(file.exists(study_file)){
      combine_meta_p_binary(read.csv(study_file,header=T),6,meta_file,single_study_I2=1)
    }
  }
  meta_vsv<-paste0("07.Microbial_GWAS/meta_adjAbunPCs/",outcome,"_vsv.csv")
  meta_dsv<-paste0("07.Microbial_GWAS/meta_adjAbunPCs/",outcome,"_dsv.csv")
  if(file.exists(meta_vsv) && file.exists(meta_dsv)){
    multi_adjust_sv(
      meta_vsv,
      meta_dsv,
      paste0("07.Microbial_GWAS/multi_adjAbunPCs/",outcome,"_sv.csv")
    )
  }
}

for(outcome in model2_continuous){
  for(sv_type in c("vsv","dsv")){
    study_file<-paste0("07.Microbial_GWAS/studies_adjAbunPCs/",outcome,"_",sv_type,".csv")
    meta_file<-paste0("07.Microbial_GWAS/meta_adjAbunPCs/",outcome,"_",sv_type,".csv")
    if(file.exists(study_file)){
      combine_meta_p_lr(read.csv(study_file,header=T),4,meta_file,single_study_I2=1)
    }
  }
  meta_vsv<-paste0("07.Microbial_GWAS/meta_adjAbunPCs/",outcome,"_vsv.csv")
  meta_dsv<-paste0("07.Microbial_GWAS/meta_adjAbunPCs/",outcome,"_dsv.csv")
  if(file.exists(meta_vsv) && file.exists(meta_dsv)){
    multi_adjust_sv(
      meta_vsv,
      meta_dsv,
      paste0("07.Microbial_GWAS/multi_adjAbunPCs/",outcome,"_sv.csv")
    )
  }
}

for(outcome in model2_binary_single){
  for(sv_type in c("vsv","dsv")){
    study_file<-paste0("07.Microbial_GWAS/studies_adjAbunPCs/",outcome,"_",sv_type,".csv")
    meta_file<-paste0("07.Microbial_GWAS/meta_adjAbunPCs/",outcome,"_",sv_type,".csv")
    if(file.exists(study_file)){
      meta_data<-read.csv(study_file,header=T)
      meta_data$pvalue_meta<-meta_data[[4]]
      meta_data$I2<-rep(1,nrow(meta_data))
      meta_data$het_p<-rep(1,nrow(meta_data))
      write.csv(meta_data,meta_file,row.names=F)
    }
  }
  meta_vsv<-paste0("07.Microbial_GWAS/meta_adjAbunPCs/",outcome,"_vsv.csv")
  meta_dsv<-paste0("07.Microbial_GWAS/meta_adjAbunPCs/",outcome,"_dsv.csv")
  if(file.exists(meta_vsv) && file.exists(meta_dsv)){
    multi_adjust_sv(
      meta_vsv,
      meta_dsv,
      paste0("07.Microbial_GWAS/multi_adjAbunPCs/",outcome,"_sv.csv")
    )
  }
}



######################################################################################################
###  significant p values

###################### Inputs
info<- read.table("01.cleanData/SV_info_QC/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)
vsv_info<-read.table("01.cleanData/SV_info_QC/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info_QC/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

vsv_info$type <- "vsv"
dsv_info$type <- "dsv"
sv_info<-rbind(vsv_info,dsv_info)

## Study-wide correction across 23 primary phenotypes and the current QC SVs

# Retain the prespecified 23 primary phenotypes; update the SV denominator after QC.
n_sv_tests <- ncol(dsgv) + ncol(vsgv)
study_wide_cutoff <- 0.05 / (23L * n_sv_tests)
write.csv(data.frame(phenotypes=23L, SVs=n_sv_tests, raw_P_cutoff=study_wide_cutoff),
          "07.Microbial_GWAS/study_wide_threshold.csv", row.names=FALSE)
analysis_wide_cutoff<-0.05
I2_cutoff<-75

select_meta_sig<-function(datain,info,dataout,p_col,cutoff){
    meta_pvalue<-read.csv(datain,header=T)
    meta_pvalue.sig<-meta_pvalue[!is.na(meta_pvalue[[p_col]]) & meta_pvalue[[p_col]] < cutoff & !is.na(meta_pvalue$I2) & meta_pvalue$I2<I2_cutoff,,drop=F]
    join_by <- if ("type" %in% colnames(meta_pvalue.sig) && "type" %in% colnames(info)) {
        c("X" = "SV_Name", "type" = "type")
    } else {
        info <- info[!duplicated(info$SV_Name), , drop = FALSE]
        c("X" = "SV_Name")
    }
    meta_pvalue.sig<-left_join(meta_pvalue.sig, info, by = join_by)
    write.csv(meta_pvalue.sig,dataout,row.names=F)
}


select_model_meta_sig<-function(meta_dir,model_name){
    meta_files<-list.files(meta_dir,pattern="_sv\\.csv$",full.names=T)
    for(datain in meta_files){
        base_name<-basename(datain)
        select_meta_sig(
            datain,sv_info,
            file.path("07.Microbial_GWAS/study_wide_sig",model_name,base_name),
            "pvalue_meta",study_wide_cutoff
        )
        select_meta_sig(
            datain,sv_info,
            file.path("07.Microbial_GWAS/analysis_wide_sig",model_name,base_name),
            "bonferroni",analysis_wide_cutoff
        )
    }
}

## Select significant associations for both models from the combined FDR tables.
select_model_meta_sig("07.Microbial_GWAS/multi_adjAbun", "adjAbun")
select_model_meta_sig("07.Microbial_GWAS/multi_adjAbunPCs", "adjAbunPCs")


########################################################################################
## annotation p values

### diseaase
### IBD T2D CRC severe obesity  CAD
### ACVD adenoma cirrhosis PD
### index 
## age  BMI fasting_insulin  hba1c  fasting_glucose 
## triglycerides hdl  ldl  cholesterol
## dyastolic_p systolic_p 
## creatinine  alt  ast  total_bilirubin 



stopifnot("SV_Name" %in% names(sv_info))

for (suffix in c("", "PCs")) {
  in_dir  <- paste0("07.Microbial_GWAS/multi_adjAbun", suffix)
  out_dir <- paste0("07.Microbial_GWAS/anno_adjAbun", suffix)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(in_dir, "\\.csv$", full.names = TRUE, ignore.case = TRUE)

  for (f in files) {
    dat  <- read.csv(f)
    info <- sv_info
    by   <- c("X" = "SV_Name")

    if ("type" %in% names(dat) && "type" %in% names(info)) {
      by <- c(by, "type" = "type")
    } else {
      info <- info[!duplicated(info$SV_Name), , drop = FALSE]
    }

    write.csv(unique(dplyr::left_join(dat, info, by = by)),
              file.path(out_dir, basename(f)), row.names = FALSE)
  }
}




#test<-read.csv("07.Microbial_GWAS/fdr_noadjAbun/CRC_sv.csv")
#freq<-table(test$X)
#write.csv(freq,"CRC_SV_freq.csv")

#freq<-table(sv_info$SV_Name)
#write.csv(freq,"info_SV_freq.csv")





###################################################################
### select the significant SVs and genes covered 

library(RCurl)
library(stringr)
library(XML)

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")


#################################################################################
## prepare the significant sites 
## functions

convert<-function(datain,dataout=NULL){
  # Retain an SV without gene annotation as one row with an empty gene.
  if(nrow(datain)==0){
    return(data.frame(genes=character(),SV=character(),stringsAsFactors=FALSE))
  }
  gene_lists<-lapply(as.character(datain$Genes),function(genes){
    if(is.na(genes) || genes=="") return("")
    strsplit(genes,";",fixed=TRUE)[[1]]
  })
  data.frame(genes=unlist(gene_lists,use.names=FALSE),
             SV=rep(as.character(datain$SV_Name_long),lengths(gene_lists)),
             stringsAsFactors=FALSE)
}

sv_organism_metadata<-function(organisms,infor){
  organism<-data.frame(organism=unique(organisms),stringsAsFactors=FALSE)
  organism<-merge(organism,infor,by="organism",all.x=TRUE)
  organism$genome_name<-paste(organism$NCBI_taxonomy_id,
                              organism$NCBI_bioproject_accession,"fasta",sep=".")
  organism$genome_name[is.na(organism$NCBI_taxonomy_id) |
                         is.na(organism$NCBI_bioproject_accession)]<-NA_character_
  organism[,c("organism","genome_name")]
}

infor<- read.table("01.cleanData/SV_info_QC/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)


##############################################################################################
###  analysis_wide: model 1

age_sv<-read.csv("07.Microbial_GWAS/analysis_wide_sig/adjAbun/age_sv.csv",header=T)
age_sv_sub<-age_sv[,c("SV_Name_long","Genes")]
age_sv_sites<-convert(age_sv_sub,age_sv_sites)
age_sv_sites$phenotype<-rep("age_sv",nrow(age_sv_sites))
age_sv_beta<-age_sv[,c("SV_Name_long","effect_meta","pvalue_meta","X","type")]
age_sv_sites<-merge(age_sv_sites,age_sv_beta,by.x="SV",by.y="SV_Name_long")

BMI_sv<-read.csv("07.Microbial_GWAS/analysis_wide_sig/adjAbun/BMI_sv.csv",header=T)
BMI_sv_sub<-BMI_sv[,c("SV_Name_long","Genes")]
BMI_sv_sites<-convert(BMI_sv_sub,BMI_sv_sites)
BMI_sv_sites$phenotype<-rep("BMI_sv",nrow(BMI_sv_sites))
BMI_sv_beta<-BMI_sv[,c("SV_Name_long","effect_meta","pvalue_meta","X","type")]
BMI_sv_sites<-merge(BMI_sv_sites,BMI_sv_beta,by.x="SV",by.y="SV_Name_long")

IBD_sv<-read.csv("07.Microbial_GWAS/analysis_wide_sig/adjAbun/IBD_sv.csv",header=T)
IBD_sv_sub<-IBD_sv[,c("SV_Name_long","Genes")]
IBD_sv_sites<-convert(IBD_sv_sub,IBD_sv_sites)
IBD_sv_sites$phenotype<-rep("IBD_sv",nrow(IBD_sv_sites))
IBD_sv_beta<-IBD_sv[,c("SV_Name_long","OR_meta","pvalue_meta","X","type")]
colnames(IBD_sv_beta)[2]<-"effect_meta"
IBD_sv_sites<-merge(IBD_sv_sites,IBD_sv_beta,by.x="SV",by.y="SV_Name_long")

cirrhosis_sv<-read.csv("07.Microbial_GWAS/analysis_wide_sig/adjAbun/cirrhosis_sv.csv",header=T)
cirrhosis_sv_sub<-cirrhosis_sv[,c("SV_Name_long","Genes")]
cirrhosis_sv$OR_meta<-cirrhosis_sv$QinN_2014.OR
cirrhosis_sv_sites<-convert(cirrhosis_sv_sub,cirrhosis_sv_sites)
cirrhosis_sv_sites$phenotype<-rep("cirrhosis_sv",nrow(cirrhosis_sv_sites))
cirrhosis_sv_beta<-cirrhosis_sv[,c("SV_Name_long","OR_meta","pvalue_meta","X","type")]
colnames(cirrhosis_sv_beta)[2]<-"effect_meta"
cirrhosis_sv_sites<-merge(cirrhosis_sv_sites,cirrhosis_sv_beta,by.x="SV",by.y="SV_Name_long")

triglycerides_sv<-read.csv("07.Microbial_GWAS/analysis_wide_sig/adjAbun/triglycerides_sv.csv",header=T)
triglycerides_sv_sub<-triglycerides_sv[,c("SV_Name_long","Genes")]
triglycerides_sv_sites<-convert(triglycerides_sv_sub,triglycerides_sv_sites)
triglycerides_sv_sites$phenotype<-rep("triglycerides_sv",nrow(triglycerides_sv_sites))
triglycerides_sv_beta<-triglycerides_sv[,c("SV_Name_long","effect_meta","pvalue_meta","X","type")]
triglycerides_sv_sites<-merge(triglycerides_sv_sites,triglycerides_sv_beta,by.x="SV",by.y="SV_Name_long")

whole<-rbind(age_sv_sites,BMI_sv_sites,IBD_sv_sites,triglycerides_sv_sites,cirrhosis_sv_sites)

whole$organism<-str_replace_all(whole$SV,"\\:\\d+_\\d+.*","")
organism<-sv_organism_metadata(whole$organism,infor)

write.csv(organism,"07.Microbial_GWAS/models_sig_genes/SV_organism_analysis_wide_model1.csv",row.names=F)

infor_sub<-infor[,c("organism","Short_name")]
whole<-merge(whole,infor_sub,by.x="organism",by.y="organism",all.x=TRUE)
whole<-merge(whole,organism,by.x="organism",by.y="organism",all.x=TRUE)
write.csv(whole,"07.Microbial_GWAS/models_sig_genes/SV_sig_genes_analysis_wide_model1.csv",row.names=F)

unique_gene<-whole[,c("Short_name","organism","genes")]
nrow(unique_gene)

unique_gene<-unique(unique_gene)
nrow(unique_gene)

write.csv(unique_gene,"07.Microbial_GWAS/models_sig_genes/SV_sig_unique_genes_analysis_wide_model1.csv",row.names=F)


##############################################################################################
###  analysis_wide: model 2

age_sv<-read.csv("07.Microbial_GWAS/analysis_wide_sig/adjAbunPCs/age_sv.csv",header=T)
age_sv_sub<-age_sv[,c("SV_Name_long","Genes")]
age_sv_sites<-convert(age_sv_sub,age_sv_sites)
age_sv_sites$phenotype<-rep("age_sv",nrow(age_sv_sites))
age_sv_beta<-age_sv[,c("SV_Name_long","effect_meta","pvalue_meta","X","type")]
age_sv_sites<-merge(age_sv_sites,age_sv_beta,by.x="SV",by.y="SV_Name_long")

BMI_sv<-read.csv("07.Microbial_GWAS/analysis_wide_sig/adjAbunPCs/BMI_sv.csv",header=T)
BMI_sv_sub<-BMI_sv[,c("SV_Name_long","Genes")]
BMI_sv_sites<-convert(BMI_sv_sub,BMI_sv_sites)
BMI_sv_sites$phenotype<-rep("BMI_sv",nrow(BMI_sv_sites))
BMI_sv_beta<-BMI_sv[,c("SV_Name_long","effect_meta","pvalue_meta","X","type")]
BMI_sv_sites<-merge(BMI_sv_sites,BMI_sv_beta,by.x="SV",by.y="SV_Name_long")

triglycerides_sv<-read.csv("07.Microbial_GWAS/analysis_wide_sig/adjAbunPCs/triglycerides_sv.csv",header=T)
triglycerides_sv_sub<-triglycerides_sv[,c("SV_Name_long","Genes")]
triglycerides_sv_sites<-convert(triglycerides_sv_sub,triglycerides_sv_sites)
triglycerides_sv_sites$phenotype<-rep("triglycerides_sv",nrow(triglycerides_sv_sites))
triglycerides_sv_beta<-triglycerides_sv[,c("SV_Name_long","effect_meta","pvalue_meta","X","type")]
triglycerides_sv_sites<-merge(triglycerides_sv_sites,triglycerides_sv_beta,by.x="SV",by.y="SV_Name_long")

whole<-rbind(age_sv_sites,BMI_sv_sites,triglycerides_sv_sites)
whole$organism<-str_replace_all(whole$SV,"\\:\\d+_\\d+.*","")

whole<-merge(whole,infor_sub,by.x="organism",by.y="organism",all.x=TRUE)
whole<-merge(whole,organism,by.x="organism",by.y="organism",all.x=TRUE)
write.csv(whole,"07.Microbial_GWAS/models_sig_genes/SV_sig_genes_analysis_wide_model2.csv",row.names=F)

unique_gene<-whole[,c("Short_name","organism","genes")]
nrow(unique_gene)

unique_gene<-unique(unique_gene)
nrow(unique_gene)

write.csv(unique_gene,"07.Microbial_GWAS/models_sig_genes/SV_sig_unique_genes_analysis_wide_model2.csv",row.names=F)



##############################################################################################
###  study_wide: model 1

age_sv<-read.csv("07.Microbial_GWAS/study_wide_sig/adjAbun/age_sv.csv",header=T)
age_sv_sub<-age_sv[,c("SV_Name_long","Genes")]
age_sv_sites<-convert(age_sv_sub,age_sv_sites)
age_sv_sites$phenotype<-rep("age_sv",nrow(age_sv_sites))
age_sv_beta<-age_sv[,c("SV_Name_long","effect_meta","pvalue_meta","X","type")]
age_sv_sites<-merge(age_sv_sites,age_sv_beta,by.x="SV",by.y="SV_Name_long")

BMI_sv<-read.csv("07.Microbial_GWAS/study_wide_sig/adjAbun/BMI_sv.csv",header=T)
BMI_sv_sub<-BMI_sv[,c("SV_Name_long","Genes")]
BMI_sv_sites<-convert(BMI_sv_sub,BMI_sv_sites)
BMI_sv_sites$phenotype<-rep("BMI_sv",nrow(BMI_sv_sites))
BMI_sv_beta<-BMI_sv[,c("SV_Name_long","effect_meta","pvalue_meta","X","type")]
BMI_sv_sites<-merge(BMI_sv_sites,BMI_sv_beta,by.x="SV",by.y="SV_Name_long")

whole<-rbind(age_sv_sites,BMI_sv_sites)

whole$organism<-str_replace_all(whole$SV,"\\:\\d+_\\d+.*","")
organism<-data.frame(unique(whole$organism))
colnames(organism)<-"organism"

infor_sub<-infor[,c("organism","Short_name")]
whole<-merge(whole,infor_sub,by.x="organism",by.y="organism",all.x=TRUE)
whole<-merge(whole,organism,by.x="organism",by.y="organism",all.x=TRUE)
write.csv(whole,"07.Microbial_GWAS/models_sig_genes/SV_sig_genes_study_wide_model1.csv",row.names=F)

unique_gene<-whole[,c("Short_name","organism","genes")]
nrow(unique_gene)

unique_gene<-unique(unique_gene)
nrow(unique_gene)

write.csv(unique_gene,"07.Microbial_GWAS/models_sig_genes/SV_sig_unique_genes_study_wide_model1.csv",row.names=F)


##############################################################################################
### study_wide: model 2

age_sv<-read.csv("07.Microbial_GWAS/study_wide_sig/adjAbunPCs/age_sv.csv",header=T)
age_sv_sub<-age_sv[,c("SV_Name_long","Genes")]
age_sv_sites<-convert(age_sv_sub,age_sv_sites)
age_sv_sites$phenotype<-rep("age_sv",nrow(age_sv_sites))
age_sv_beta<-age_sv[,c("SV_Name_long","effect_meta","pvalue_meta","X","type")]
age_sv_sites<-merge(age_sv_sites,age_sv_beta,by.x="SV",by.y="SV_Name_long")

BMI_sv<-read.csv("07.Microbial_GWAS/study_wide_sig/adjAbunPCs/BMI_sv.csv",header=T)
BMI_sv_sub<-BMI_sv[,c("SV_Name_long","Genes")]
BMI_sv_sites<-convert(BMI_sv_sub,BMI_sv_sites)
BMI_sv_sites$phenotype<-rep("BMI_sv",nrow(BMI_sv_sites))
BMI_sv_beta<-BMI_sv[,c("SV_Name_long","effect_meta","pvalue_meta","X","type")]
BMI_sv_sites<-merge(BMI_sv_sites,BMI_sv_beta,by.x="SV",by.y="SV_Name_long")

whole<-rbind(age_sv_sites,BMI_sv_sites)

whole$organism<-str_replace_all(whole$SV,"\\:\\d+_\\d+.*","")
organism<-data.frame(unique(whole$organism))
colnames(organism)<-"organism"

whole<-merge(whole,infor_sub,by.x="organism",by.y="organism",all.x=TRUE)
whole<-merge(whole,organism,by.x="organism",by.y="organism",all.x=TRUE)
write.csv(whole,"07.Microbial_GWAS/models_sig_genes/SV_sig_genes_study_wide_model2.csv",row.names=F)

unique_gene<-whole[,c("Short_name","organism","genes")]
nrow(unique_gene)

unique_gene<-unique(unique_gene)
nrow(unique_gene)

write.csv(unique_gene,"07.Microbial_GWAS/models_sig_genes/SV_sig_unique_genes_study_wide_model2.csv",row.names=F)




#############################################################################################################
### overlap

# Retain no-annotation SVs in the main tables; overlap requires a known gene.
overlap_labels<-c("study_wide_model1","study_wide_model2","analysis_wide_model1","analysis_wide_model2")
unknown_gene_ids<-c("","na","n/a","nan","unknown","unannotated","unknown / unannotated")
overlap<-NULL
for(label in overlap_labels){
  gene_set<-read.csv(paste0("07.Microbial_GWAS/models_sig_genes/SV_sig_unique_genes_",label,".csv"),
                     header=T,stringsAsFactors=F)
  gene_set$organism<-trimws(as.character(gene_set$organism))
  gene_set$genes<-trimws(as.character(gene_set$genes))
  known_gene<-!is.na(gene_set$organism) & nzchar(gene_set$organism) &
              !is.na(gene_set$genes) & !tolower(gene_set$genes) %in% unknown_gene_ids
  gene_set<-unique(gene_set[known_gene,c("organism","genes"),drop=F])
  gene_set[[label]]<-rep("T",nrow(gene_set))
  if(is.null(overlap)){
    overlap<-gene_set
  }else{
    overlap<-merge(overlap,gene_set,by=c("organism","genes"),all=T)
  }
}
write.csv(overlap,"07.Microbial_GWAS/models_sig_genes/overlap_genes.csv",row.names=F)




