### Microbial GWA 
### 2024-8-13
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

## 1 Preparation
### 1.1 Import
#install.packages("survminer")
rm(list = ls())
#options(scipen = 200)

source("functions.R")
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
info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

all_vsv_info<-read.table("01.cleanData/SV_info/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
all_dsv_info<-read.table("01.cleanData/SV_info/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

all_basic<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")
#all_basic$hscrp <- as.numeric(replace(all_basic$hscrp, all_basic$hscrp =="<5", NA))

load("01.cleanData/SV_all/dsgv.RData")
load("01.cleanData/SV_all/vsgv.RData")

vsgv<-vsgv_sub
dsgv<-dsgv_sub

vsgv$"Phascolarctobacterium sp. CAG:207:1597_1599;1599_1600_2"<-vsgv$"Phascolarctobacterium sp. CAG:207:1597_1599;1599_1600"
all_abun<-read.table("01.cleanData/mbio_all/SV_species_abun.tsv", check.names = F) 

## 2 Associations between SVs and prognosis
### 2.1 Preparation
if (!dir.exists("07.Microbial_GWAS")) {dir.create("07.Microbial_GWAS")}
if (!dir.exists("07.Microbial_GWAS/studies_adjAbun")) {dir.create("07.Microbial_GWAS/studies_adjAbun")}
if (!dir.exists("07.Microbial_GWAS/studies_noadjAbun")) {dir.create("07.Microbial_GWAS/studies_noadjAbun")}
##if (!dir.exists("07.Microbial_GWAS/meta_adjAbun")) {dir.create("07.Microbial_GWAS/meta_adjAbun")}
if (!dir.exists("07.Microbial_GWAS/meta_noadjAbun")) {dir.create("07.Microbial_GWAS/meta_noadjAbun")}
##if (!dir.exists("07.Microbial_GWAS/fdr_adjAbun")) {dir.create("07.Microbial_GWAS/fdr_adjAbun")}
if (!dir.exists("07.Microbial_GWAS/fdr_noadjAbun")) {dir.create("07.Microbial_GWAS/fdr_noadjAbun")}
##if (!dir.exists("07.Microbial_GWAS/sig_adjAbun")) {dir.create("07.Microbial_GWAS/sig_adjAbun")}
if (!dir.exists("07.Microbial_GWAS/sig_noadjAbun")) {dir.create("07.Microbial_GWAS/sig_noadjAbun")}
##if (!dir.exists("07.Microbial_GWAS/anno_adjAbun")) {dir.create("07.Microbial_GWAS/anno_adjAbun")}
if (!dir.exists("07.Microbial_GWAS/anno_noadjAbun")) {dir.create("07.Microbial_GWAS/anno_noadjAbun")}


# all
all_abun$Others<-1-rowSums(all_abun)
all_abun_clr<-abundances(x=as.data.frame(na.omit(all_abun)), transform="clr") %>%as.data.frame
all_abun_clr <- all_abun_clr[match(rownames(all_abun), rownames(all_abun_clr)),]
rownames(all_abun_clr) <- rownames(all_abun)
#all_covar<-cbind(all_basic,all_abun_clr)

all_covar<-cbind(all_basic,all_abun)
all_abun_clr<-all_abun


############################################################################################################
## associations between disease and SVs

disease_asso<-function(basic,dis_name,adjust_cov,out_vsv_adjust,out_dsv_adjust,out_vsv_noadjust,out_dsv_noadjust){
   #basic<-basic_CRC
   #dis_name<-"CRC"
   #adjust_cov<-c("age","gender_code","BMI")
   #prog<-basic[,c(dis_name,dis_name)]

   prog<-data.frame(basic[,dis_name])
   row.names(prog)<-rownames(basic)
   colnames(prog)<-"Y"

   abun_clr<-all_abun_clr[rownames(basic),]

   dis_vsgv<-vsgv[rownames(basic),]
   dis_dsgv<-dsgv[rownames(basic),]

   vsv_adjAbun<-logistic_adjAbun_vsv(prog,dis_vsgv,basic,adjust_cov,abun_clr,info,1,dis_name)
   dsv_adjAbun<-logistic_adjAbun_dsv(prog,dis_dsgv,basic,adjust_cov,abun_clr,info,1,dis_name)

   vsv_noadjAbun<-logistic_noadjAbun_vsv(prog,dis_vsgv,basic,adjust_cov,info,dis_name)
   dsv_noadjAbun<-logistic_noadjAbun_dsv(prog,dis_dsgv,basic,adjust_cov,info,dis_name)

   write.csv(vsv_adjAbun, file = out_vsv_adjust,row.names=F)
   write.csv(dsv_adjAbun, file = out_dsv_adjust,row.names=F)

   write.csv(vsv_noadjAbun, file = out_vsv_noadjust,row.names=F)
   write.csv(dsv_noadjAbun, file = out_dsv_noadjust,row.names=F)
}


### diseaase
### IBD T2D CRC severe obesity CAD
### ACVD adenoma cirrhosis  MetS

basic_CRC<-subset(all_basic,is.na(CRC)==F,drop=T)
a<-table(basic_CRC$study_name,basic_CRC$CRC)
basic_CRC<-subset(all_basic,is.na(CRC)==F & study_name!="FengQ_2015",drop=T)

disease_asso(basic_CRC,"CRC",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/CRC_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/CRC_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/CRC_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/CRC_dsv.csv")

basic_T2D<-subset(all_basic,is.na(T2D)==F,drop=T)
a<-table(basic_T2D$study_name,basic_T2D$T2D)
basic_T2D<-subset(all_basic,is.na(T2D)==F & study_name!="HMP_2019_t2d" & study_name!="LiJ_2014",drop=T)
disease_asso(basic_T2D,"T2D",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/T2D_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/T2D_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/T2D_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/T2D_dsv.csv")

basic_adenoma<-subset(all_basic,is.na(adenoma)==F,drop=T)
table(basic_adenoma$study_name,basic_adenoma$adenoma)
basic_adenoma<-subset(all_basic,is.na(adenoma)==F & study_name!="FengQ_2015",drop=T)
disease_asso(basic_adenoma,"adenoma",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/adenoma_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/adenoma_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/adenoma_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/adenoma_dsv.csv")

basic_ACVD<-subset(all_basic,is.na(ACVD)==F,drop=T)
table(basic_ACVD$study_name,basic_ACVD$ACVD)
disease_asso(basic_ACVD,"ACVD",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/ACVD_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/ACVD_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/ACVD_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/ACVD_dsv.csv")

basic_severe_obesity<-subset(all_basic,is.na(severe_obesity)==F,drop=T)
table(basic_severe_obesity$study_name,basic_severe_obesity$severe_obesity)
disease_asso(basic_severe_obesity,"severe_obesity",c("age","gender_code"),
"07.Microbial_GWAS/studies_adjAbun/severe_obesity_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/severe_obesity_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/severe_obesity_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/severe_obesity_dsv.csv")

basic_CAD<-subset(all_basic,is.na(CAD)==F,drop=T)
table(basic_CAD$study_name,basic_CAD$CAD)
disease_asso(basic_CAD,"CAD",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/CAD_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/CAD_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/CAD_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/CAD_dsv.csv")

basic_cirrhosis<-subset(all_basic,is.na(cirrhosis)==F,drop=T)
table(basic_cirrhosis$study_name,basic_cirrhosis$cirrhosis)
disease_asso(basic_cirrhosis,"cirrhosis",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/cirrhosis_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/cirrhosis_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/cirrhosis_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/cirrhosis_dsv.csv")

basic_MetS<-subset(all_basic,is.na(MetS)==F,drop=T)
table(basic_MetS$study_name,basic_MetS$MetS)
disease_asso(basic_MetS,"MetS",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/MetS_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/MetS_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/MetS_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/MetS_dsv.csv")


disease_asso_nocovar<-function(basic,dis_name,out_vsv_adjust,out_dsv_adjust,out_vsv_noadjust,out_dsv_noadjust){

   #basic<-basic_IBD
   #dis_name<-"IBD"
   #prog<-basic[,c(dis_name,dis_name)]

   prog<-data.frame(basic[,dis_name])
   row.names(prog)<-rownames(basic)
   colnames(prog)<-"Y"

   abun_clr<-all_abun_clr[rownames(basic),]

   dis_vsgv<-vsgv[rownames(basic),]
   dis_dsgv<-dsgv[rownames(basic),]

   vsv_adjAbun<-logistic_adjAbun_nocovar_vsv(prog,dis_vsgv,basic,abun_clr,info,1,dis_name)
   dsv_adjAbun<-logistic_adjAbun_nocovar_dsv(prog,dis_dsgv,basic,abun_clr,info,1,dis_name)

   vsv_noadjAbun<-logistic_noadjAbun_nocovar_vsv(prog,dis_vsgv,basic,info,dis_name)
   dsv_noadjAbun<-logistic_noadjAbun_nocovar_dsv(prog,dis_dsgv,basic,info,dis_name)

   write.csv(vsv_adjAbun, file = out_vsv_adjust,row.names=F)
   write.csv(dsv_adjAbun, file = out_dsv_adjust,row.names=F)

   write.csv(vsv_noadjAbun, file = out_vsv_noadjust,row.names=F)
   write.csv(dsv_noadjAbun, file = out_dsv_noadjust,row.names=F)
}

basic_IBD<-subset(all_basic,is.na(IBD)==F,drop=T)
table(basic_IBD$study_name,basic_IBD$IBD)
basic_IBD<-subset(all_basic,is.na(IBD)==F & study_name!="LiJ_2014",drop=T)

disease_asso_nocovar(basic_IBD,"IBD",
"07.Microbial_GWAS/studies_adjAbun/IBD_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/IBD_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/IBD_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/IBD_dsv.csv")


### 2.3 linear regression model 1
### Linear model with covariates: age, gender, BMI, read count.
### 1 indicate the organism name column of the info file


############################################################################################
#######  index 
## age  BMI fasting_insulin  hba1c  fasting_glucose 
## triglycerides hdl  ldl  cholesterol 
## dyastolic_p systolic_p prothrombin_time inr
## albumine creatinine  alt  ast  total_bilirubin


############# BMI healthy 
BMI_basic<-subset(all_basic,is.na(BMI)==F & disease=="healthy",drop=T)

BMI_prog<-data.frame(BMI_basic[,"BMI"])
row.names(BMI_prog)<-rownames(BMI_basic)
colnames(BMI_prog)<-"Y"

BMI_abun_clr<-all_abun_clr[rownames(BMI_basic),]

BMI_vsgv<-vsgv[rownames(BMI_basic),]
BMI_dsgv<-dsgv[rownames(BMI_basic),]

vsv_BMI_adjAbun<-lm_mats_adjAbun_vsv(BMI_prog,BMI_vsgv,BMI_basic,c("gender_code"),BMI_abun_clr,info,1,"BMI")
dsv_BMI_adjAbun<-lm_mats_adjAbun_dsv(BMI_prog,BMI_dsgv,BMI_basic,c("gender_code"),BMI_abun_clr,info,1,"BMI")

write.csv(vsv_BMI_adjAbun, file = "07.Microbial_GWAS/studies_adjAbun/BMI_healthy_vsv.csv",row.names=F)
write.csv(dsv_BMI_adjAbun, file = "07.Microbial_GWAS/studies_adjAbun/BMI_healthy_dsv.csv",row.names=F)

vsv_BMI_noadjAbun<-lm_mats_noadjAbun_vsv(BMI_prog,BMI_vsgv,BMI_basic,c("gender_code"),info,"BMI")
dsv_BMI_noadjAbun<-lm_mats_noadjAbun_dsv(BMI_prog,BMI_dsgv,BMI_basic,c("gender_code"),info,"BMI")

write.csv(vsv_BMI_noadjAbun, file = "07.Microbial_GWAS/studies_noadjAbun/BMI_healthy_vsv.csv",row.names=F)
write.csv(dsv_BMI_noadjAbun, file = "07.Microbial_GWAS/studies_noadjAbun/BMI_healthy_dsv.csv",row.names=F)


############# age healthy 
age_basic<-subset(all_basic,is.na(age)==F & disease=="healthy" & study_name!="LiJ_2014" & 
study_name!="FranzosaEA_2018" & study_name!="HMP_2019_t2d" & study_name!="HMP_2019_ibdmdb",drop=T)

table(age_basic$study_name)

age_prog<-data.frame(age_basic[,"age"])
row.names(age_prog)<-rownames(age_basic)
colnames(age_prog)<-"Y"

age_abun_clr<-all_abun_clr[rownames(age_basic),]

age_vsgv<-vsgv[rownames(age_basic),]
age_dsgv<-dsgv[rownames(age_basic),]

vsv_age_adjAbun<-lm_mats_adjAbun_vsv(age_prog,age_vsgv,age_basic,c("gender_code"),age_abun_clr,info,1,"age")
dsv_age_adjAbun<-lm_mats_adjAbun_dsv(age_prog,age_dsgv,age_basic,c("gender_code"),age_abun_clr,info,1,"age")

write.csv(vsv_age_adjAbun, file = "07.Microbial_GWAS/studies_adjAbun/age_healthy_vsv.csv",row.names=F)
write.csv(dsv_age_adjAbun, file = "07.Microbial_GWAS/studies_adjAbun/age_healthy_dsv.csv",row.names=F)

vsv_age_noadjAbun<-lm_mats_noadjAbun_vsv(age_prog,age_vsgv,age_basic,c("gender_code"),info,"age")
dsv_age_noadjAbun<-lm_mats_noadjAbun_dsv(age_prog,age_dsgv,age_basic,c("gender_code"),info,"age")

write.csv(vsv_age_noadjAbun, file = "07.Microbial_GWAS/studies_noadjAbun/age_healthy_vsv.csv",row.names=F)
write.csv(dsv_age_noadjAbun, file = "07.Microbial_GWAS/studies_noadjAbun/age_healthy_dsv.csv",row.names=F)



index_asso<-function(clin,index,adjust_covar,out_vsv_adjust,out_dsv_adjust,out_vsv_noadjust,out_dsv_noadjust){
   #basic<-basic_CRC
   #dis_name<-"CRC"
   #adjust_cov<-c("age","gender_code","BMI")
   #prog<-basic[,c(dis_name,dis_name)]

   prog<-data.frame(clin[,index])
   row.names(prog)<-rownames(clin)
   colnames(prog)<-"Y"

   abun_clr<-all_abun_clr[rownames(clin),]

   index_vsgv<-vsgv[rownames(clin),]
   index_dsgv<-dsgv[rownames(clin),]

   vsv_adjAbun<-lm_mats_adjAbun_vsv(prog,index_vsgv,clin,adjust_covar,abun_clr,info,1,index)
   dsv_adjAbun<-lm_mats_adjAbun_dsv(prog,index_dsgv,clin,adjust_covar,abun_clr,info,1,index)

   vsv_noadjAbun<-lm_mats_noadjAbun_vsv(prog,index_vsgv,clin,adjust_covar,info,index)
   dsv_noadjAbun<-lm_mats_noadjAbun_dsv(prog,index_dsgv,clin,adjust_covar,info,index)

   write.csv(vsv_adjAbun, file = out_vsv_adjust,row.names=F)
   write.csv(dsv_adjAbun, file = out_dsv_adjust,row.names=F)

   write.csv(vsv_noadjAbun, file = out_vsv_noadjust,row.names=F)
   write.csv(dsv_noadjAbun, file = out_dsv_noadjust,row.names=F)
}


basic_hba1c<-subset(all_basic,is.na(hba1c)==F,drop=T)
index_asso(basic_hba1c,"hba1c",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/hba1c_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/hba1c_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/hba1c_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/hba1c_dsv.csv")

basic_triglycerides<-subset(all_basic,is.na(triglycerides)==F,drop=T)
table(basic_triglycerides$study_name)
index_asso(basic_triglycerides,"triglycerides",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/triglycerides_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/triglycerides_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/triglycerides_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/triglycerides_dsv.csv")

basic_fasting_insulin<-subset(all_basic,is.na(fasting_insulin)==F,drop=T)
table(basic_fasting_insulin$study_name)
index_asso(basic_fasting_insulin,"fasting_insulin",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/fasting_insulin_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/fasting_insulin_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/fasting_insulin_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/fasting_insulin_dsv.csv")

basic_hdl<-subset(all_basic,is.na(hdl)==F,drop=T)
table(basic_hdl$study_name)
index_asso(basic_hdl,"hdl",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/hdl_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/hdl_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/hdl_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/hdl_dsv.csv")

basic_ldl<-subset(all_basic,is.na(ldl)==F,drop=T)
table(basic_ldl$study_name)
index_asso(basic_ldl,"ldl",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/ldl_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/ldl_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/ldl_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/ldl_dsv.csv")

basic_cholesterol<-subset(all_basic,is.na(cholesterol)==F,drop=T)
table(basic_cholesterol$study_name)
index_asso(basic_cholesterol,"cholesterol",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/cholesterol_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/cholesterol_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/cholesterol_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/cholesterol_dsv.csv")

basic_mean_press<-subset(all_basic,is.na(mean_press)==F,drop=T)
index_asso(basic_mean_press,"mean_press",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/mean_press_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/mean_press_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/mean_press_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/mean_press_dsv.csv")

basic_ast<-subset(all_basic,is.na(ast)==F,drop=T)
index_asso(basic_ast,"ast",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/ast_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/ast_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/ast_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/ast_dsv.csv")

basic_alt<-subset(all_basic,is.na(alt)==F,drop=T)
index_asso(basic_alt,"alt",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/alt_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/alt_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/alt_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/alt_dsv.csv")

basic_total_bilirubin<-subset(all_basic,is.na(total_bilirubin)==F,drop=T)
index_asso(basic_total_bilirubin,"total_bilirubin",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/total_bilirubin_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/total_bilirubin_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/total_bilirubin_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/total_bilirubin_dsv.csv")

basic_creatinine<-subset(all_basic,is.na(creatinine)==F,drop=T)
index_asso(basic_creatinine,"creatinine",c("age","gender_code","BMI"),
"07.Microbial_GWAS/studies_adjAbun/creatinine_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/creatinine_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/creatinine_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/creatinine_dsv.csv")



index_asso_nocovar<-function(clin,index,out_vsv_adjust,out_dsv_adjust,out_vsv_noadjust,out_dsv_noadjust){
   #basic<-basic_CRC
   #dis_name<-"CRC"
   #adjust_cov<-c("age","gender_code","BMI")
   #prog<-basic[,c(dis_name,dis_name)]

   prog<-data.frame(clin[,index])
   row.names(prog)<-rownames(clin)
   colnames(prog)<-"Y"

   abun_clr<-all_abun_clr[rownames(clin),]

   index_vsgv<-vsgv[rownames(clin),]
   index_dsgv<-dsgv[rownames(clin),]

   vsv_adjAbun<-lm_mats_adjAbun_nocovar_vsv(prog,index_vsgv,clin,abun_clr,info,1,index)
   dsv_adjAbun<-lm_mats_adjAbun_nocovar_dsv(prog,index_dsgv,clin,abun_clr,info,1,index)

   vsv_noadjAbun<-lm_mats_noadjAbun_nocovar_vsv(prog,index_vsgv,clin,info,index)
   dsv_noadjAbun<-lm_mats_noadjAbun_nocovar_dsv(prog,index_dsgv,clin,info,index)

   write.csv(vsv_adjAbun, file = out_vsv_adjust,row.names=F)
   write.csv(dsv_adjAbun, file = out_dsv_adjust,row.names=F)

   write.csv(vsv_noadjAbun, file = out_vsv_noadjust,row.names=F)
   write.csv(dsv_noadjAbun, file = out_dsv_noadjust,row.names=F)
}

## no covariates 
basic_fasting_glucose<-subset(all_basic,is.na(fasting_glucose)==F,drop=T)
index_asso_nocovar(basic_fasting_glucose,"fasting_glucose",
"07.Microbial_GWAS/studies_adjAbun/fasting_glucose_vsv.csv",
"07.Microbial_GWAS/studies_adjAbun/fasting_glucose_dsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/fasting_glucose_vsv.csv",
"07.Microbial_GWAS/studies_noadjAbun/fasting_glucose_dsv.csv")


#################################################################################################################################
###  meta-analysis  (no adjust abundance) more than one datasets 
################### disease 

#### CRC
CRC_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/CRC_dsv.csv",header=T)
combine_meta_p_binary(CRC_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/CRC_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/CRC_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/CRC_dsv.csv")

CRC_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/CRC_vsv.csv",header=T)
combine_meta_p_binary(CRC_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/CRC_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/CRC_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/CRC_vsv.csv")


#### IBD
IBD_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/IBD_dsv.csv",header=T)
combine_meta_p_binary(IBD_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/IBD_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/IBD_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/IBD_dsv.csv")

IBD_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/IBD_vsv.csv",header=T)
combine_meta_p_binary(IBD_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/IBD_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/IBD_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/IBD_vsv.csv")


#### T2D
T2D_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/T2D_dsv.csv",header=T)
combine_meta_p_binary(T2D_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/T2D_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/T2D_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/T2D_dsv.csv")

T2D_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/T2D_vsv.csv",header=T)
combine_meta_p_binary(T2D_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/T2D_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/T2D_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/T2D_vsv.csv")


#### adenoma
adenoma_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/adenoma_dsv.csv",header=T)
combine_meta_p_binary(adenoma_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/adenoma_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/adenoma_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/adenoma_dsv.csv")

adenoma_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/adenoma_vsv.csv",header=T)
combine_meta_p_binary(adenoma_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/adenoma_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/adenoma_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/adenoma_vsv.csv")

###################################################################################################
###  index

####### age
age_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/age_healthy_dsv.csv",header=T)
combine_meta_p_lr(age_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/age_healthy_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/age_healthy_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/age_healthy_dsv.csv")

age_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/age_healthy_vsv.csv",header=T)
combine_meta_p_lr(age_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/age_healthy_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/age_healthy_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/age_healthy_vsv.csv")

####### BMI
BMI_healthy_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/BMI_healthy_dsv.csv",header=T)
combine_meta_p_lr(BMI_healthy_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/BMI_healthy_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/BMI_healthy_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/BMI_healthy_dsv.csv")

BMI_healthy_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/BMI_healthy_vsv.csv",header=T)
combine_meta_p_lr(BMI_healthy_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/BMI_healthy_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/BMI_healthy_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/BMI_healthy_vsv.csv")

#### hba1c
hba1c_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/hba1c_dsv.csv",header=T)
combine_meta_p_lr(hba1c_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/hba1c_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/hba1c_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/hba1c_dsv.csv")

hba1c_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/hba1c_vsv.csv",header=T)
combine_meta_p_lr(hba1c_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/hba1c_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/hba1c_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/hba1c_vsv.csv")


#### fasting_glucose
fasting_glucose_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/fasting_glucose_dsv.csv",header=T)
combine_meta_p_lr(fasting_glucose_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/fasting_glucose_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/fasting_glucose_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/fasting_glucose_dsv.csv")

fasting_glucose_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/fasting_glucose_vsv.csv",header=T)
combine_meta_p_lr(fasting_glucose_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/fasting_glucose_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/fasting_glucose_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/fasting_glucose_vsv.csv")


#### fasting_insulin
fasting_insulin_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/fasting_insulin_dsv.csv",header=T)
combine_meta_p_lr(fasting_insulin_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/fasting_insulin_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/fasting_insulin_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/fasting_insulin_dsv.csv")

fasting_insulin_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/fasting_insulin_vsv.csv",header=T)
combine_meta_p_lr(fasting_insulin_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/fasting_insulin_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/fasting_insulin_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/fasting_insulin_vsv.csv")


#### triglycerides
triglycerides_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/triglycerides_dsv.csv",header=T)
combine_meta_p_lr(triglycerides_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/triglycerides_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/triglycerides_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/triglycerides_dsv.csv")

triglycerides_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/triglycerides_vsv.csv",header=T)
combine_meta_p_lr(triglycerides_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/triglycerides_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/triglycerides_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/triglycerides_vsv.csv")


#### hdl
hdl_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/hdl_dsv.csv",header=T)
combine_meta_p_lr(hdl_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/hdl_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/hdl_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/hdl_dsv.csv")

hdl_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/hdl_vsv.csv",header=T)
combine_meta_p_lr(hdl_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/hdl_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/hdl_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/hdl_vsv.csv")


#### ldl
ldl_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/ldl_dsv.csv",header=T)
combine_meta_p_lr(ldl_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/ldl_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/ldl_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/ldl_dsv.csv")

ldl_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/ldl_vsv.csv",header=T)
combine_meta_p_lr(ldl_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/ldl_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/ldl_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/ldl_vsv.csv")


#### cholesterol
cholesterol_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/cholesterol_dsv.csv",header=T)
combine_meta_p_lr(cholesterol_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/cholesterol_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/cholesterol_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/cholesterol_dsv.csv")

cholesterol_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/cholesterol_vsv.csv",header=T)
combine_meta_p_lr(cholesterol_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/cholesterol_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/cholesterol_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/cholesterol_vsv.csv")

#### mean_press
mean_press_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/mean_press_dsv.csv",header=T)
combine_meta_p_lr(mean_press_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/mean_press_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/mean_press_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/mean_press_dsv.csv")

mean_press_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/mean_press_vsv.csv",header=T)
combine_meta_p_lr(mean_press_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/mean_press_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/mean_press_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/mean_press_vsv.csv")

#### creatinine
creatinine_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/creatinine_dsv.csv",header=T)
combine_meta_p_lr(creatinine_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/creatinine_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/creatinine_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/creatinine_dsv.csv")

creatinine_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/creatinine_vsv.csv",header=T)
combine_meta_p_lr(creatinine_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/creatinine_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/creatinine_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/creatinine_vsv.csv")

#### ast
ast_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/ast_dsv.csv",header=T)
combine_meta_p_lr(ast_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/ast_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/ast_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/ast_dsv.csv")

ast_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/ast_vsv.csv",header=T)
combine_meta_p_lr(ast_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/ast_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/ast_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/ast_vsv.csv")

#### alt
alt_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/alt_dsv.csv",header=T)
combine_meta_p_lr(alt_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/alt_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/alt_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/alt_dsv.csv")

alt_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/alt_vsv.csv",header=T)
combine_meta_p_lr(alt_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/alt_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/alt_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/alt_vsv.csv")

#### total_bilirubin
total_bilirubin_dsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/total_bilirubin_dsv.csv",header=T)
combine_meta_p_lr(total_bilirubin_dsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/total_bilirubin_dsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/total_bilirubin_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/total_bilirubin_dsv.csv")

total_bilirubin_vsv_noadjAbun<-read.csv("07.Microbial_GWAS/studies_noadjAbun/total_bilirubin_vsv.csv",header=T)
combine_meta_p_lr(total_bilirubin_vsv_noadjAbun,"07.Microbial_GWAS/meta_noadjAbun/total_bilirubin_vsv.csv")
multi_adjust("07.Microbial_GWAS/meta_noadjAbun/total_bilirubin_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/total_bilirubin_vsv.csv")


######################################################################################################
###  significant p values

###################### Inputs
info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)
vsv_info<-read.table("01.cleanData/SV_info/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")

sig_p_cutoff<-0.001
het_p_cutoff<-0.05
p_count_cutoff<-2

select_sig<-function(datain,info,dataout){
    meta_pvalue<-read.csv(datain,header=T)
    meta_pvalue.sig<-subset(meta_pvalue,pvalue_meta <=sig_p_cutoff & het_p >= het_p_cutoff & p_count >= p_count_cutoff ,drop=T)
    meta_pvalue.sig<-left_join(meta_pvalue.sig, info, by = c("X" = "SV_Name"))
    write.csv(meta_pvalue.sig,dataout,row.names=F)
    }

select_sig("07.Microbial_GWAS/fdr_noadjAbun/CRC_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/CRC_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/CRC_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/CRC_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/adenoma_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/adenoma_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/adenoma_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/adenoma_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/IBD_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/IBD_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/IBD_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/IBD_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/T2D_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/T2D_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/T2D_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/T2D_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/age_healthy_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/age_healthy_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/age_healthy_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/age_healthy_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/BMI_healthy_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/BMI_healthy_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/BMI_healthy_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/BMI_healthy_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/fasting_insulin_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/fasting_insulin_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/fasting_insulin_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/fasting_insulin_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/hba1c_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/hba1c_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/hba1c_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/hba1c_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/fasting_glucose_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/fasting_glucose_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/fasting_glucose_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/fasting_glucose_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/triglycerides_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/triglycerides_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/triglycerides_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/triglycerides_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/hdl_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/hdl_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/hdl_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/hdl_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/ldl_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/ldl_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/ldl_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/ldl_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/cholesterol_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/cholesterol_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/cholesterol_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/cholesterol_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/mean_press_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/mean_press_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/mean_press_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/mean_press_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/creatinine_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/creatinine_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/creatinine_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/creatinine_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/alt_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/alt_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/alt_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/alt_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/ast_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/ast_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/ast_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/ast_dsv.sig.anno.csv")

select_sig("07.Microbial_GWAS/fdr_noadjAbun/total_bilirubin_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/total_bilirubin_vsv.sig.anno.csv")
select_sig("07.Microbial_GWAS/fdr_noadjAbun/total_bilirubin_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/total_bilirubin_dsv.sig.anno.csv")


########################################################################################################
#####  just one dataset 

select_sig_one<-function(datain,info,dataout){
    meta_pvalue<-read.csv(datain,header=T)
    meta_pvalue.sig<-subset(meta_pvalue,pvalue_meta <=sig_p_cutoff & het_p >= het_p_cutoff ,drop=T)
    meta_pvalue.sig<-left_join(meta_pvalue.sig, info, by = c("X" = "SV_Name"))
    write.csv(meta_pvalue.sig,dataout,row.names=F)
    }


ACVD_vsv_pvalue<-read.csv("07.Microbial_GWAS/studies_noadjAbun/ACVD_vsv.csv",header=T)
colnames(ACVD_vsv_pvalue)[4]<-"pvalue_meta"
ACVD_vsv_pvalue$het_p<-1
write.csv(ACVD_vsv_pvalue,"07.Microbial_GWAS/studies_noadjAbun/ACVD_vsv.csv",row.names=F)
multi_adjust("07.Microbial_GWAS/studies_noadjAbun/ACVD_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/ACVD_vsv.csv")
select_sig_one("07.Microbial_GWAS/fdr_noadjAbun/ACVD_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/ACVD_vsv.sig.anno.csv")

ACVD_dsv_pvalue<-read.csv("07.Microbial_GWAS/studies_noadjAbun/ACVD_dsv.csv",header=T)
colnames(ACVD_dsv_pvalue)[4]<-"pvalue_meta"
ACVD_dsv_pvalue$het_p<-1
write.csv(ACVD_dsv_pvalue,"07.Microbial_GWAS/studies_noadjAbun/ACVD_dsv.csv",row.names=F)
multi_adjust("07.Microbial_GWAS/studies_noadjAbun/ACVD_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/ACVD_dsv.csv")
select_sig_one("07.Microbial_GWAS/fdr_noadjAbun/ACVD_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/ACVD_dsv.sig.anno.csv")

severe_obesity_vsv_pvalue<-read.csv("07.Microbial_GWAS/studies_noadjAbun/severe_obesity_vsv.csv",header=T)
colnames(severe_obesity_vsv_pvalue)[4]<-"pvalue_meta"
severe_obesity_vsv_pvalue$het_p<-1
write.csv(severe_obesity_vsv_pvalue,"07.Microbial_GWAS/studies_noadjAbun/severe_obesity_vsv.csv",row.names=F)
multi_adjust("07.Microbial_GWAS/studies_noadjAbun/severe_obesity_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/severe_obesity_vsv.csv")
select_sig_one("07.Microbial_GWAS/fdr_noadjAbun/severe_obesity_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/severe_obesity_vsv.sig.anno.csv")

severe_obesity_dsv_pvalue<-read.csv("07.Microbial_GWAS/studies_noadjAbun/severe_obesity_dsv.csv",header=T)
colnames(severe_obesity_dsv_pvalue)[4]<-"pvalue_meta"
severe_obesity_dsv_pvalue$het_p<-1
write.csv(severe_obesity_dsv_pvalue,"07.Microbial_GWAS/studies_noadjAbun/severe_obesity_dsv.csv",row.names=F)
multi_adjust("07.Microbial_GWAS/studies_noadjAbun/severe_obesity_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/severe_obesity_dsv.csv")
select_sig_one("07.Microbial_GWAS/fdr_noadjAbun/severe_obesity_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/severe_obesity_dsv.sig.anno.csv")

CAD_vsv_pvalue<-read.csv("07.Microbial_GWAS/studies_noadjAbun/CAD_vsv.csv",header=T)
colnames(CAD_vsv_pvalue)[4]<-"pvalue_meta"
CAD_vsv_pvalue$het_p<-1
write.csv(CAD_vsv_pvalue,"07.Microbial_GWAS/studies_noadjAbun/CAD_vsv.csv",row.names=F)
multi_adjust("07.Microbial_GWAS/studies_noadjAbun/CAD_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/CAD_vsv.csv")
select_sig_one("07.Microbial_GWAS/fdr_noadjAbun/CAD_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/CAD_vsv.sig.anno.csv")

CAD_dsv_pvalue<-read.csv("07.Microbial_GWAS/studies_noadjAbun/CAD_dsv.csv",header=T)
colnames(CAD_dsv_pvalue)[4]<-"pvalue_meta"
CAD_dsv_pvalue$het_p<-1
write.csv(CAD_dsv_pvalue,"07.Microbial_GWAS/studies_noadjAbun/CAD_dsv.csv",row.names=F)
multi_adjust("07.Microbial_GWAS/studies_noadjAbun/CAD_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/CAD_dsv.csv")
select_sig_one("07.Microbial_GWAS/fdr_noadjAbun/CAD_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/CAD_dsv.sig.anno.csv")

cirrhosis_dsv_pvalue<-read.csv("07.Microbial_GWAS/studies_noadjAbun/cirrhosis_dsv.csv",header=T)
colnames(cirrhosis_dsv_pvalue)[4]<-"pvalue_meta"
cirrhosis_dsv_pvalue$het_p<-1
write.csv(cirrhosis_dsv_pvalue,"07.Microbial_GWAS/studies_noadjAbun/cirrhosis_dsv.csv",row.names=F)
multi_adjust("07.Microbial_GWAS/studies_noadjAbun/cirrhosis_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/cirrhosis_dsv.csv")
select_sig_one("07.Microbial_GWAS/fdr_noadjAbun/cirrhosis_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/cirrhosis_dsv.sig.anno.csv")

cirrhosis_vsv_pvalue<-read.csv("07.Microbial_GWAS/studies_noadjAbun/cirrhosis_vsv.csv",header=T)
colnames(cirrhosis_vsv_pvalue)[4]<-"pvalue_meta"
cirrhosis_vsv_pvalue$het_p<-1
write.csv(cirrhosis_vsv_pvalue,"07.Microbial_GWAS/studies_noadjAbun/cirrhosis_vsv.csv",row.names=F)
multi_adjust("07.Microbial_GWAS/studies_noadjAbun/cirrhosis_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/cirrhosis_vsv.csv")
select_sig_one("07.Microbial_GWAS/fdr_noadjAbun/cirrhosis_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/cirrhosis_vsv.sig.anno.csv")

MetS_dsv_pvalue<-read.csv("07.Microbial_GWAS/studies_noadjAbun/MetS_dsv.csv",header=T)
colnames(MetS_dsv_pvalue)[4]<-"pvalue_meta"
MetS_dsv_pvalue$het_p<-1
write.csv(MetS_dsv_pvalue,"07.Microbial_GWAS/studies_noadjAbun/MetS_dsv.csv",row.names=F)
multi_adjust("07.Microbial_GWAS/studies_noadjAbun/MetS_dsv.csv","07.Microbial_GWAS/fdr_noadjAbun/MetS_dsv.csv")
select_sig_one("07.Microbial_GWAS/fdr_noadjAbun/MetS_dsv.csv",dsv_info,"07.Microbial_GWAS/sig_noadjAbun/MetS_dsv.sig.anno.csv")

MetS_vsv_pvalue<-read.csv("07.Microbial_GWAS/studies_noadjAbun/MetS_vsv.csv",header=T)
colnames(MetS_vsv_pvalue)[4]<-"pvalue_meta"
MetS_vsv_pvalue$het_p<-1
write.csv(MetS_vsv_pvalue,"07.Microbial_GWAS/studies_noadjAbun/MetS_vsv.csv",row.names=F)
multi_adjust("07.Microbial_GWAS/studies_noadjAbun/MetS_vsv.csv","07.Microbial_GWAS/fdr_noadjAbun/MetS_vsv.csv")
select_sig_one("07.Microbial_GWAS/fdr_noadjAbun/MetS_vsv.csv",vsv_info,"07.Microbial_GWAS/sig_noadjAbun/MetS_vsv.sig.anno.csv")


########################################################################################
## annotation p values

### diseaase
### IBD T2D CRC severe obesity MetS CAD
### ACVD adenoma cirrhosis 
### index 
## age  BMI fasting_insulin  hba1c  fasting_glucose 
## triglycerides hdl  ldl  cholesterol
## dyastolic_p systolic_p 
## creatinine  alt  ast  total_bilirubin 

annotation<-function(datain,info,dataout){
    meta_pvalue<-read.csv(datain,header=T)
    meta_pvalue.sig<-left_join(meta_pvalue, info, by = c("X" = "SV_Name"))
    write.csv(meta_pvalue.sig,dataout,row.names=F)
    }

annotation("07.Microbial_GWAS/fdr_noadjAbun/CRC_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/CRC_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/CRC_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/CRC_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/IBD_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/IBD_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/IBD_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/IBD_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/T2D_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/T2D_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/T2D_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/T2D_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/ACVD_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/ACVD_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/ACVD_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/ACVD_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/severe_obesity_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/severe_obesity_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/severe_obesity_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/severe_obesity_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/MetS_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/MetS_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/MetS_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/MetS_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/CAD_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/CAD_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/CAD_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/CAD_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/adenoma_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/adenoma_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/adenoma_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/adenoma_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/cirrhosis_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/cirrhosis_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/cirrhosis_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/cirrhosis_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/age_healthy_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/age_healthy_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/age_healthy_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/age_healthy_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/BMI_healthy_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/BMI_healthy_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/BMI_healthy_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/BMI_healthy_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/fasting_insulin_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/fasting_insulin_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/fasting_insulin_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/fasting_insulin_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/hba1c_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/hba1c_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/hba1c_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/hba1c_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/fasting_glucose_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/fasting_glucose_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/fasting_glucose_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/fasting_glucose_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/triglycerides_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/triglycerides_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/triglycerides_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/triglycerides_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/hdl_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/hdl_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/hdl_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/hdl_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/ldl_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/ldl_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/ldl_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/ldl_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/cholesterol_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/cholesterol_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/cholesterol_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/cholesterol_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/mean_press_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/mean_press_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/mean_press_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/mean_press_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/creatinine_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/creatinine_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/creatinine_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/creatinine_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/ast_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/ast_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/ast_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/ast_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/alt_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/alt_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/alt_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/alt_vsv.csv")

annotation("07.Microbial_GWAS/fdr_noadjAbun/total_bilirubin_dsv.csv",dsv_info,"07.Microbial_GWAS/anno_noadjAbun/total_bilirubin_dsv.csv")
annotation("07.Microbial_GWAS/fdr_noadjAbun/total_bilirubin_vsv.csv",vsv_info,"07.Microbial_GWAS/anno_noadjAbun/total_bilirubin_vsv.csv")





