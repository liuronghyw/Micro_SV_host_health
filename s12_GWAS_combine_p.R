### Microbial GWA 
### 2024-9-11
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

## 1 Preparation
### 1.1 Import
#install.packages("survminer")
rm(list = ls())
options(scipen = 200)

source("functions.R")
knitr::opts_chunk$set(echo = TRUE)
library("R.utils")

#library("survival")
#library("ggplot2")
#library("survminer")
#library(gridExtra)
#library("grid")
library(reshape2)	  
#library("RColorBrewer")
#library("plyr")


###################################################################################################
### 3 Results
### 3.1 Clean results
#### 3.1.1 vSV associations

## merge result tables

###################### 1.2 Inputs
info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

vsv_info<-read.table("01.cleanData/SV_info/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F,quote = "")


##### different variables

disease_select_sig<-function(meta_p,sig_p,sv_type,dataout){
  temp<-read.csv(meta_p,header=T)
  #temp<-read.csv("07.Microbial_GWAS/meta_noadjAbun/CRC_vsv.csv",header=T)
  #sv_type<-vsv_info
  temp<-temp[,-which(grepl("SE",colnames(temp))=="TRUE")]
  temp<-temp[,-which(grepl("se",colnames(temp))=="TRUE")]
  temp<-temp[,-which(grepl("I2",colnames(temp))=="TRUE")]

  temp<-rename(temp,c("meta.OR"=OR_meta,"meta.OR_left"=cilb,"meta.OR_right"=ciub,"meta.het_p"=het_p,
  "Number_of_dataset_with_p_less_than_0.05"=p_count,"meta.p"=pvalue_meta))

  temp_sig<-subset(temp,meta.p<0.05,drop=T)
  temp_sig_final<-left_join(temp_sig, sv_type, by = c("X" = "SV_Name"))

  temp_select<-read.csv(sig_p,header=T)
  #  temp_select<-read.csv("07.Microbial_GWAS/sig_noadjAbun/CRC_vsv.sig.anno.csv",header=T)
  temp_select$Replicate_Sig<-"Yes"
  temp_select<-temp_select[,c("X","Replicate_Sig")]
  temp_final<-merge(temp_sig_final,temp_select,by.x="X",by.y="X",all.x=T)

  temp_final<-subset(temp_final,is.na(SV_Name_long)==F,drop=T)
  row.names(temp_final)<- temp_final$SV_Name_long
  temp_final <- temp_final[,-which(names(temp_final) %in% c("X","Taxonomy_Name","Taxonomy_ID","SV_ID","Adjust_Abundance","SV_Name_long"))]

  temp_final$Replicate_Sig[is.na(temp_final$Replicate_Sig) == T] <- "No"

  nc<-ncol(temp_final)
  temp_final<-temp_final[,c(seq((nc-10),nc),seq(1,nc-11))]

  write.csv(temp_final,dataout)
  }


### no significant SVs
disease_select<-function(meta_p,sv_type,dataout){
  temp<-read.csv(meta_p,header=T)
  #temp<-read.csv("07.Microbial_GWAS/meta_noadjAbun/adenoma_dsv.csv",header=T)
  temp<-temp[,-which(grepl("SE",colnames(temp))=="TRUE")]
  temp<-temp[,-which(grepl("se",colnames(temp))=="TRUE")]
  temp<-temp[,-which(grepl("I2",colnames(temp))=="TRUE")]

  temp<-rename(temp,c("meta.OR"=OR_meta,"meta.OR_left"=cilb,"meta.OR_right"=ciub,"meta.het_p"=het_p,
  "Number_of_dataset_with_p_less_than_0.05"=p_count,"meta.p"=pvalue_meta))

  temp_sig<-subset(temp,meta.p<0.05,drop=T)
  temp_sig_final<-left_join(temp_sig, sv_type, by = c("X" = "SV_Name"))

  temp_sig_final$Replicate_Sig<-"No"
  temp_final<-temp_sig_final
  temp_final<-subset(temp_final,is.na(SV_Name_long)==F,drop=T)
  row.names(temp_final)<- temp_final$SV_Name_long
  temp_final <- temp_final[,-which(names(temp_final) %in% c("X","Taxonomy_Name","Taxonomy_ID","SV_ID","Adjust_Abundance","SV_Name_long"))]

  nc<-ncol(temp_final)
  temp_final<-temp_final[,c(seq((nc-10),nc),seq(1,nc-11))]

  write.csv(temp_final,dataout)
  }


###  adenoma
disease_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/adenoma_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/adenoma_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/adenoma_dsv.csv")

disease_select(
"07.Microbial_GWAS/meta_noadjAbun/adenoma_vsv.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/adenoma_vsv.csv")


## CRC
disease_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/CRC_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/CRC_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/CRC_dsv.csv")

disease_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/CRC_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/CRC_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/CRC_vsv.csv")


## IBD
disease_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/IBD_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/IBD_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/IBD_dsv.csv")

disease_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/IBD_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/IBD_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/IBD_vsv.csv")


## T2D
disease_select(
"07.Microbial_GWAS/meta_noadjAbun/T2D_dsv.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/T2D_dsv.csv")

disease_select(
"07.Microbial_GWAS/meta_noadjAbun/T2D_vsv.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/T2D_vsv.csv")



#### just one dataset

disease_select_one_sig<-function(meta_p,sig_p,sv_type,dataout){
  temp<-read.csv(meta_p,header=T)
  #temp<-read.csv("07.Microbial_GWAS/studies_noadjAbun/ACVD_dsv.csv",header=T)
  temp<-temp[,-which(grepl("SE",colnames(temp))=="TRUE")]

  temp_sig<-subset(temp,pvalue_meta<0.05,drop=T)
  temp_sig_final<-left_join(temp_sig, sv_type, by = c("X" = "SV_Name"))

  temp_select<-read.csv(sig_p,header=T)
  #  temp_select<-read.csv("07.Microbial_GWAS/sig_noadjAbun/CRC_vsv.sig.anno.csv",header=T)
  temp_select$Replicate_Sig<-"Yes"
  temp_select<-temp_select[,c("X","Replicate_Sig")]
  temp_final<-merge(temp_sig_final,temp_select,by.x="X",by.y="X",all.x=T)

  temp_final<-subset(temp_final,is.na(SV_Name_long)==F,drop=T)
  row.names(temp_final)<- temp_final$SV_Name_long
  temp_final <- temp_final[,-which(names(temp_final) %in% c("X","Taxonomy_Name","Taxonomy_ID","SV_ID","Adjust_Abundance","SV_Name_long"))]

  temp_final$Replicate_Sig[is.na(temp_final$Replicate_Sig) == T] <- "No"

  nc<-ncol(temp_final)
  temp_final<-temp_final[,c(seq((nc-5),nc),seq(1,nc-6))]

  write.csv(temp_final,dataout)
  }


### no significant SVs
disease_select_one<-function(meta_p,sv_type,dataout){
  temp<-read.csv(meta_p,header=T)
  #temp<-read.csv("07.Microbial_GWAS/meta_noadjAbun/adenoma_dsv.csv",header=T)
  temp<-temp[,-which(grepl("SE",colnames(temp))=="TRUE")]

  temp_sig<-subset(temp,pvalue_meta<0.05,drop=T)
  temp_sig_final<-left_join(temp_sig, sv_type, by = c("X" = "SV_Name"))

  temp_sig_final$Replicate_Sig<-"No"
  temp_final<-temp_sig_final
  temp_final<-subset(temp_final,is.na(SV_Name_long)==F,drop=T)
  row.names(temp_final)<- temp_final$SV_Name_long
  temp_final <- temp_final[,-which(names(temp_final) %in% c("X","Taxonomy_Name","Taxonomy_ID","SV_ID","Adjust_Abundance","SV_Name_long"))]

  nc<-ncol(temp_final)
  temp_final<-temp_final[,c(seq((nc-5),nc),seq(1,nc-6))]

  write.csv(temp_final,dataout)
  }



## ACVD
disease_select_one_sig(
"07.Microbial_GWAS/studies_noadjAbun/ACVD_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/ACVD_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/ACVD_dsv.csv")

disease_select_one_sig(
"07.Microbial_GWAS/studies_noadjAbun/ACVD_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/ACVD_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/ACVD_vsv.csv")

## CAD
disease_select_one(
"07.Microbial_GWAS/studies_noadjAbun/CAD_dsv.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/CAD_dsv.csv")

disease_select_one(
"07.Microbial_GWAS/studies_noadjAbun/CAD_vsv.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/CAD_vsv.csv")

## cirrhosis
disease_select_one_sig(
"07.Microbial_GWAS/studies_noadjAbun/cirrhosis_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/cirrhosis_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/cirrhosis_dsv.csv")

disease_select_one_sig(
"07.Microbial_GWAS/studies_noadjAbun/cirrhosis_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/cirrhosis_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/cirrhosis_vsv.csv")


## MetS
disease_select_one(
"07.Microbial_GWAS/studies_noadjAbun/MetS_dsv.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/MetS_dsv.csv")

disease_select_one(
"07.Microbial_GWAS/studies_noadjAbun/MetS_vsv.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/MetS_vsv.csv")



####################################################################################
### indexes

index_select_sig<-function(meta_p,sig_p,sv_type,dataout){
  temp<-read.csv(meta_p,header=T)
  #temp<-read.csv("07.Microbial_GWAS/meta_noadjAbun/fasting_glucose_dsv.csv",header=T)
  temp<-temp[,-which(grepl("SE",colnames(temp))=="TRUE")]
  temp<-temp[,-which(grepl("se",colnames(temp))=="TRUE")]
  temp<-temp[,-which(grepl("I2",colnames(temp))=="TRUE")]

  temp<-rename(temp,c("meta.beta"=effect_meta,"meta.het_p"=het_p,
  "Number_of_dataset_with_p_less_than_0.05"=p_count,"meta.p"=pvalue_meta))

  temp_sig<-subset(temp,meta.p<0.05,drop=T)
  temp_sig_final<-left_join(temp_sig, sv_type, by = c("X" = "SV_Name"))

  temp_select<-read.csv(sig_p,header=T)
  #  temp_select<-read.csv("07.Microbial_GWAS/sig_noadjAbun/fasting_glucose_dsv.sig.anno.csv",header=T)
  temp_select$Replicate_Sig<-"Yes"
  temp_select<-temp_select[,c("X","Replicate_Sig")]
  temp_final<-merge(temp_sig_final,temp_select,by.x="X",by.y="X",all.x=T)
  temp_final<-subset(temp_final,is.na(SV_Name_long)==F,drop=T)

  row.names(temp_final)<- temp_final$SV_Name_long
  temp_final <- temp_final[,-which(names(temp_final) %in% c("X","Taxonomy_Name","Taxonomy_ID","SV_ID","Adjust_Abundance","SV_Name_long"))]

  temp_final$Replicate_Sig[is.na(temp_final$Replicate_Sig) == T] <- "No"

  nc<-ncol(temp_final)
  temp_final<-temp_final[,c(seq((nc-8),nc),seq(1,nc-9))]

  write.csv(temp_final,dataout)
  }



### no significant SVs
index_select<-function(meta_p,sv_type,dataout){
  temp<-read.csv(meta_p,header=T)
  #temp<-read.csv("07.Microbial_GWAS/meta_noadjAbun/adenoma_dsv.csv",header=T)
  temp<-temp[,-which(grepl("SE",colnames(temp))=="TRUE")]
  temp<-temp[,-which(grepl("se",colnames(temp))=="TRUE")]
  temp<-temp[,-which(grepl("I2",colnames(temp))=="TRUE")]

  temp<-rename(temp,c("meta.beta"=effect_meta,"meta.het_p"=het_p,
  "Number_of_dataset_with_p_less_than_0.05"=p_count,"meta.p"=pvalue_meta))

  temp_sig<-subset(temp,meta.p<0.05,drop=T)
  temp_sig_final<-left_join(temp_sig, sv_type, by = c("X" = "SV_Name"))

  temp_sig_final$Replicate_Sig<-"No"
  temp_final<-temp_sig_final
  temp_final<-subset(temp_final,is.na(SV_Name_long)==F,drop=T)
  row.names(temp_final)<- temp_final$SV_Name_long
  temp_final <- temp_final[,-which(names(temp_final) %in% c("X","Taxonomy_Name","Taxonomy_ID","SV_ID","Adjust_Abundance","SV_Name_long"))]

  nc<-ncol(temp_final)
  temp_final<-temp_final[,c(seq((nc-8),nc),seq(1,nc-9))]

  write.csv(temp_final,dataout)
  }



index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/age_healthy_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/age_healthy_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/age_healthy_dsv.csv")


index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/age_healthy_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/age_healthy_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/age_healthy_vsv.csv")


index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/BMI_healthy_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/BMI_healthy_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/BMI_healthy_dsv.csv")


index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/BMI_healthy_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/BMI_healthy_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/BMI_healthy_vsv.csv")


index_select(
"07.Microbial_GWAS/meta_noadjAbun/fasting_insulin_dsv.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/fasting_insulin_dsv.csv")


index_select(
"07.Microbial_GWAS/meta_noadjAbun/fasting_insulin_vsv.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/fasting_insulin_vsv.csv")


###  fasting_glucose
index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/fasting_glucose_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/fasting_glucose_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/fasting_glucose_dsv.csv")

index_select(
"07.Microbial_GWAS/meta_noadjAbun/fasting_glucose_vsv.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/fasting_glucose_vsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/hba1c_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/hba1c_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/hba1c_dsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/hba1c_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/hba1c_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/hba1c_vsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/triglycerides_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/triglycerides_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/triglycerides_dsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/triglycerides_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/triglycerides_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/triglycerides_vsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/hdl_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/hdl_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/hdl_dsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/hdl_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/hdl_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/hdl_vsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/ldl_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/ldl_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/ldl_dsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/ldl_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/ldl_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/ldl_vsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/cholesterol_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/cholesterol_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/cholesterol_dsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/cholesterol_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/cholesterol_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/cholesterol_vsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/mean_press_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/mean_press_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/mean_press_dsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/mean_press_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/mean_press_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/mean_press_vsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/creatinine_dsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/creatinine_dsv.sig.anno.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/creatinine_dsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/creatinine_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/creatinine_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/creatinine_vsv.csv")

index_select(
"07.Microbial_GWAS/meta_noadjAbun/alt_dsv.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/alt_dsv.csv")

index_select_sig(
"07.Microbial_GWAS/meta_noadjAbun/alt_vsv.csv",
"07.Microbial_GWAS/sig_noadjAbun/alt_vsv.sig.anno.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/alt_vsv.csv")

index_select(
"07.Microbial_GWAS/meta_noadjAbun/ast_dsv.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/ast_dsv.csv")

index_select(
"07.Microbial_GWAS/meta_noadjAbun/ast_vsv.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/ast_vsv.csv")

index_select(
"07.Microbial_GWAS/meta_noadjAbun/total_bilirubin_dsv.csv",
dsv_info,
"07.Microbial_GWAS/combine_noadjAbun/total_bilirubin_dsv.csv")

index_select(
"07.Microbial_GWAS/meta_noadjAbun/total_bilirubin_vsv.csv",
vsv_info,
"07.Microbial_GWAS/combine_noadjAbun/total_bilirubin_vsv.csv")










#### the species count

dsv_os<-read.csv("07.Microbial_GWAS/combine/melanoma_dsv_os.csv",header=T)$Taxonomy_Name
vsv_os<-read.csv("07.Microbial_GWAS/combine/melanoma_vsv_os.csv",header=T)$Taxonomy_Name
os<-c(dsv_os,vsv_os)
length(unique(os))




vsv_resp<-read.csv("07.Microbial_GWAS/combine/melanoma_vsv_resp.csv",header=T)$Taxonomy_Name
dsv_resp<-read.csv("07.Microbial_GWAS/combine/melanoma_dsv_resp.csv",header=T)$Taxonomy_Name
resp<-c(vsv_resp,dsv_resp)
length(unique(resp))



