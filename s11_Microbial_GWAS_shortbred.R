### Microbial GWA 
### 2026-6-16
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
library(data.table)


###################### 1.2 Inputs
all_basic<-read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv")

sig_SVs_shortbred<-fread("00.rawData/shortBRED/sig_SVs/matched_rows.csv",header=T)
sig_SVs_shortbred_t<-t(sig_SVs_shortbred[,-1])
colnames(sig_SVs_shortbred_t)<-sig_SVs_shortbred[,1]$Family

MAGs_genes_shortbred<-fread("00.rawData/shortBRED/MAGs/unmatched_filtered_zeros_lt80pct.csv",header=T)
MAGs_genes_shortbred_t<-t(MAGs_genes_shortbred[,-1])
colnames(MAGs_genes_shortbred_t)<-MAGs_genes_shortbred[,1]$Family

#sub_sample<-data.frame(setdiff(all_basic$match_id, colnames(sig_SVs_shortbred)[-1]))
#colnames(sub_sample)<-"match_id"
#study_sample<-all_basic[,c("match_id","study_name")]
#study_sample_sub<-merge(study_sample,sub_sample,by.x="match_id",by.y="match_id")
#write.csv(study_sample_sub,"shortbred_sub.csv")


## 2 Associations between genes and prognosis
### 2.1 Preparation
if (!dir.exists("08.Microbial_shortbred")) {dir.create("08.Microbial_shortbred")}
if (!dir.exists("08.Microbial_shortbred/sig_SVs")) {dir.create("08.Microbial_shortbred/sig_SVs")}
if (!dir.exists("08.Microbial_shortbred/sig_SVs/studies")) {dir.create("08.Microbial_shortbred/sig_SVs/studies")}
if (!dir.exists("08.Microbial_shortbred/sig_SVs/meta")) {dir.create("08.Microbial_shortbred/sig_SVs/meta")}
if (!dir.exists("08.Microbial_shortbred/sig_SVs/multi")) {dir.create("08.Microbial_shortbred/sig_SVs/multi")}
if (!dir.exists("08.Microbial_shortbred/sig_SVs/sig")) {dir.create("08.Microbial_shortbred/sig_SVs/sig")}
if (!dir.exists("08.Microbial_shortbred/sig_SVs/anno")) {dir.create("08.Microbial_shortbred/sig_SVs/anno")}

if (!dir.exists("08.Microbial_shortbred")) {dir.create("08.Microbial_shortbred")}
if (!dir.exists("08.Microbial_shortbred/MAGs")) {dir.create("08.Microbial_shortbred/MAGs")}
if (!dir.exists("08.Microbial_shortbred/MAGs/studies")) {dir.create("08.Microbial_shortbred/MAGs/studies")}
if (!dir.exists("08.Microbial_shortbred/MAGs/meta")) {dir.create("08.Microbial_shortbred/MAGs/meta")}
if (!dir.exists("08.Microbial_shortbred/MAGs/multi")) {dir.create("08.Microbial_shortbred/MAGs/multi")}
if (!dir.exists("08.Microbial_shortbred/MAGs/sig")) {dir.create("08.Microbial_shortbred/MAGs/sig")}
if (!dir.exists("08.Microbial_shortbred/MAGs/anno")) {dir.create("08.Microbial_shortbred/MAGs/anno")}



### disease
### IBD T2D CRC CAD
### ACVD adenoma cirrhosis MetS PD

disease_shortbred_asso <- function(basic, index_name, adjust_cov, out_SVs_shortbred,out_MAGs_shortbred) {
  # basic <- basic_hba1c
  # index_name <- "hba1c_norm"
  # adjust_cov <- c("age", "gender_code", "BMI", "read_count")

  prog <- data.frame(basic[, index_name])
  row.names(prog) <- rownames(basic)
  colnames(prog) <- "Y"

  common_rows <- intersect(rownames(basic), rownames(sig_SVs_shortbred_t))
  #SVs_shortbred <- sig_SVs_shortbred_t[rownames(basic), ]
  SVs_shortbred <- sig_SVs_shortbred_t[common_rows, ]

  SVs_shortbred <- logistic_shortbred(prog, SVs_shortbred, basic, adjust_cov, index_name)
  write.csv(SVs_shortbred, file = out_SVs_shortbred, row.names = FALSE)

  common_rows <- intersect(rownames(basic), rownames(MAGs_genes_shortbred_t))
  #MAGs_shortbred <- MAGs_genes_shortbred_t[rownames(basic), ]
  MAGs_shortbred <- MAGs_genes_shortbred_t[common_rows, ]

  MAGs_shortbred <- logistic_shortbred(prog, MAGs_shortbred, basic, adjust_cov, index_name)
  write.csv(MAGs_shortbred, file = out_MAGs_shortbred, row.names = FALSE)
}


basic_CRC<-subset(all_basic,is.na(CRC)==F,drop=T)
a<-table(basic_CRC$study_name,basic_CRC$CRC)
basic_CRC<-subset(all_basic,is.na(CRC)==F & study_name!="FengQ_2015",drop=T)
disease_shortbred_asso(basic_CRC,"CRC",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/CRC.csv",
"08.Microbial_shortbred/MAGs/studies/CRC.csv")

basic_T2D<-subset(all_basic,is.na(T2D)==F,drop=T)
a<-table(basic_T2D$study_name,basic_T2D$T2D)
basic_T2D<-subset(all_basic,is.na(T2D)==F & study_name!="HMP_2019_t2d" & study_name!="LiJ_2014",drop=T)
disease_shortbred_asso(basic_T2D,"T2D",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/T2D.csv",
"08.Microbial_shortbred/MAGs/studies/T2D.csv")

basic_ACVD<-subset(all_basic,is.na(ACVD)==F,drop=T)
table(basic_ACVD$study_name,basic_ACVD$ACVD)
disease_shortbred_asso(basic_ACVD,"ACVD",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/ACVD.csv",
"08.Microbial_shortbred/MAGs/studies/ACVD.csv")

basic_CAD<-subset(all_basic,is.na(CAD)==F,drop=T)
table(basic_CAD$study_name,basic_CAD$CAD)
disease_shortbred_asso(basic_CAD,"CAD",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/CVD.csv",
"08.Microbial_shortbred/MAGs/studies/CVD.csv")

basic_cirrhosis<-subset(all_basic,is.na(cirrhosis)==F,drop=T)
table(basic_cirrhosis$study_name,basic_cirrhosis$cirrhosis)
disease_shortbred_asso(basic_cirrhosis,"cirrhosis",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/cirrhosis.csv",
"08.Microbial_shortbred/MAGs/studies/cirrhosis.csv")

basic_MetS<-subset(all_basic,is.na(MetS)==F,drop=T)
table(basic_MetS$study_name,basic_MetS$MetS)
disease_shortbred_asso(basic_MetS,"MetS",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/MetS.csv",
"08.Microbial_shortbred/MAGs/studies/MetS.csv")

basic_adenoma<-subset(all_basic,is.na(adenoma)==F,drop=T)
table(basic_adenoma$study_name,basic_adenoma$adenoma)
basic_adenoma<-subset(all_basic,is.na(adenoma)==F & study_name!="FengQ_2015",drop=T)
disease_shortbred_asso(basic_adenoma,"adenoma",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/adenoma.csv",
"08.Microbial_shortbred/MAGs/studies/adenoma.csv")

basic_PD<-subset(all_basic,is.na(PD)==F,drop=T)
table(basic_PD$study_name,basic_PD$PD)
disease_shortbred_asso(basic_PD,"PD",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/PD.csv",
"08.Microbial_shortbred/MAGs/studies/PD.csv")

basic_IBD<-subset(all_basic,is.na(IBD)==F,drop=T)
table(basic_IBD$study_name,basic_IBD$IBD)
basic_IBD<-subset(all_basic,is.na(IBD)==F & study_name!="LiJ_2014",drop=T)
disease_shortbred_asso(basic_IBD,"IBD",c(),
"08.Microbial_shortbred/sig_SVs/studies/IBD.csv",
"08.Microbial_shortbred/MAGs/studies/IBD.csv")



### disease
### IBD T2D CRC CAD
### ACVD adenoma cirrhosis MetS PD


############################################################################################
#######  index 
#######  2.3 linear regression model 1
## age  BMI fasting_insulin  hba1c  fasting_glucose 
## triglycerides hdl  ldl  cholesterol 
## dyastolic_p systolic_p prothrombin_time inr
## albumine creatinine  alt  ast  total_bilirubin

## 区分dsv和vSV

SV_sig_sites_genes_m1<-read.csv("07.Microbial_GWAS/models_sig_genes/SV_sig_genes_analysis_wide_model1.csv",header=T)
SV_sig_sites_genes_m2<-read.csv("07.Microbial_GWAS/models_sig_genes/SV_sig_genes_analysis_wide_model2.csv",header=T)

SV_sig_sites_genes_type<-unique(rbind(SV_sig_sites_genes_m1[,c("genes","type")],SV_sig_sites_genes_m2[,c("genes","type")]))

SV_sig_sites_genes_dsv<-subset(SV_sig_sites_genes_type,type=="dsv",drop=T)
SV_sig_sites_genes_vsv<-subset(SV_sig_sites_genes_type,type=="vsv",drop=T)

dsv_genes_to_keep <- unique(SV_sig_sites_genes_dsv$genes)
vsv_genes_to_keep <- unique(SV_sig_sites_genes_vsv$genes)


############# age all
age_basic<-subset(all_basic,is.na(age)==F ,drop=T)
table(age_basic$study_name)

age_basic<-subset(age_basic,study_name!="LiJ_2014",drop=T)
table(age_basic$study_name)

age_prog<-data.frame(age_basic[,"age_norm"])
row.names(age_prog)<-rownames(age_basic)
colnames(age_prog)<-"Y"

common_rows <- intersect(rownames(age_basic), rownames(sig_SVs_shortbred_t))
#index_SVs_shortbred <- sig_SVs_shortbred_t[rownames(age_basic), ]
index_SVs_shortbred <- sig_SVs_shortbred_t[common_rows, ]

dsv_cols_keep <- grep(paste(dsv_genes_to_keep, collapse = "|"), 
                  colnames(index_SVs_shortbred), 
                  value = TRUE)

vsv_cols_keep <- grep(paste(vsv_genes_to_keep, collapse = "|"), 
                  colnames(index_SVs_shortbred), 
                  value = TRUE)

index_dSVs_shortbred <- index_SVs_shortbred[, dsv_cols_keep]
index_vSVs_shortbred <- index_SVs_shortbred[, vsv_cols_keep]

SVs_shortbred_vsv <- linear_shortbred_vsv(age_prog, index_vSVs_shortbred, age_basic,c("gender_code"), "age_norm")
write.csv(SVs_shortbred_vsv, file = "08.Microbial_shortbred/sig_SVs/studies/age_vsv.csv", row.names = FALSE)

SVs_shortbred_dsv <- linear_shortbred_dsv(age_prog, index_dSVs_shortbred, age_basic,c("gender_code"), "age_norm")
write.csv(SVs_shortbred_dsv, file = "08.Microbial_shortbred/sig_SVs/studies/age_dsv.csv", row.names = FALSE)


common_rows <- intersect(rownames(age_basic), rownames(MAGs_genes_shortbred_t))
#index_MAGs_shortbred <- MAGs_genes_shortbred_t[rownames(age_basic), ]
index_MAGs_shortbred <- MAGs_genes_shortbred_t[common_rows, ]

MAGs_shortbred <- linear_shortbred(age_prog, index_MAGs_shortbred, age_basic, c("gender_code"),"age_norm")
write.csv(MAGs_shortbred, file = "08.Microbial_shortbred/MAGs/studies/age.csv", row.names = FALSE)



####################################################################################################
############### BMI all
BMI_basic<-subset(all_basic,is.na(BMI)==F,header=T)
table(BMI_basic$study_name)

BMI_basic<-subset(BMI_basic,study_name!="PiccinnoG_2025_cohort6",drop=T)

BMI_prog<-data.frame(BMI_basic[,"BMI_norm"])
row.names(BMI_prog)<-rownames(BMI_basic)
colnames(BMI_prog)<-"Y"

common_rows <- intersect(rownames(BMI_basic), rownames(sig_SVs_shortbred_t))
#index_SVs_shortbred <- sig_SVs_shortbred_t[rownames(BMI_basic), ]
index_SVs_shortbred <- sig_SVs_shortbred_t[common_rows, ]

dsv_cols_keep <- grep(paste(dsv_genes_to_keep, collapse = "|"), 
                  colnames(index_SVs_shortbred), 
                  value = TRUE)

vsv_cols_keep <- grep(paste(vsv_genes_to_keep, collapse = "|"), 
                  colnames(index_SVs_shortbred), 
                  value = TRUE)

index_dSVs_shortbred <- index_SVs_shortbred[, dsv_cols_keep]
index_vSVs_shortbred <- index_SVs_shortbred[, vsv_cols_keep]

dSVs_shortbred <- linear_shortbred_dsv(BMI_prog, index_dSVs_shortbred, BMI_basic,c("gender_code"), "BMI_norm")
write.csv(dSVs_shortbred, file = "08.Microbial_shortbred/sig_SVs/studies/BMI_dsv.csv", row.names = FALSE)

vSVs_shortbred <- linear_shortbred_vsv(BMI_prog, index_vSVs_shortbred, BMI_basic,c("gender_code"), "BMI_norm")
write.csv(vSVs_shortbred, file = "08.Microbial_shortbred/sig_SVs/studies/BMI_vsv.csv", row.names = FALSE)

common_rows <- intersect(rownames(BMI_basic), rownames(MAGs_genes_shortbred_t))
#index_MAGs_shortbred <- MAGs_genes_shortbred_t[rownames(BMI_basic), ]
index_MAGs_shortbred <- MAGs_genes_shortbred_t[common_rows, ]
MAGs_shortbred <- linear_shortbred(BMI_prog, index_MAGs_shortbred, BMI_basic, c("gender_code"), "BMI_norm")
write.csv(MAGs_shortbred, file = "08.Microbial_shortbred/MAGs/studies/BMI.csv", row.names = FALSE)




############################################################################################################
## associations between index and genes

index_shortbred_asso <- function(basic, index_name, adjust_cov, out_dSVs_shortbred_adjust,out_vSVs_shortbred_adjust,out_MAGs_shortbred_adjust) {
  # basic <- basic_hba1c
  # index_name <- "hba1c_norm"
  # adjust_cov <- c("age", "gender_code", "BMI", "read_count")

  prog <- data.frame(basic[, index_name])
  row.names(prog) <- rownames(basic)
  colnames(prog) <- "Y"

  # 找出 basic 的行名中哪些不在 sig_SVs_shortbred_t 中
  #missing_rows <- rownames(basic)[!rownames(basic) %in% rownames(sig_SVs_shortbred_t)]
  #print(missing_rows)
  #print(length(missing_rows)

  common_rows <- intersect(rownames(basic), rownames(sig_SVs_shortbred_t))
  #index_SVs_shortbred <- sig_SVs_shortbred_t[rownames(basic), ]
  index_SVs_shortbred <- sig_SVs_shortbred_t[common_rows, ]

  index_dSVs_shortbred <- index_SVs_shortbred[, dsv_cols_keep]
  index_vSVs_shortbred <- index_SVs_shortbred[, vsv_cols_keep]

   dSVs_shortbred <- linear_shortbred_dsv(prog, index_SVs_shortbred, basic, adjust_cov, index_name)
   write.csv(dSVs_shortbred, file = out_dSVs_shortbred_adjust, row.names = FALSE)

   vSVs_shortbred <- linear_shortbred_vsv(prog, index_SVs_shortbred, basic, adjust_cov, index_name)
   write.csv(vSVs_shortbred, file =  out_vSVs_shortbred_adjust, row.names = FALSE)

   common_rows <- intersect(rownames(basic), rownames(MAGs_genes_shortbred_t))
   #index_MAGs_shortbred <- MAGs_genes_shortbred_t[rownames(basic), ]
   index_MAGs_shortbred <- MAGs_genes_shortbred_t[common_rows, ]
   MAGs_shortbred <- linear_shortbred(prog, index_MAGs_shortbred, basic, adjust_cov, index_name)
   write.csv(MAGs_shortbred, file = out_MAGs_shortbred_adjust, row.names = FALSE)
}


####################################################################################
#### index

basic_triglycerides<-subset(all_basic,is.na(triglycerides)==F,drop=T)
table(basic_triglycerides$study_name)
index_shortbred_asso(basic_triglycerides,"triglycerides_norm",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/triglycerides_dsv.csv",
"08.Microbial_shortbred/sig_SVs/studies/triglycerides_vsv.csv",
"08.Microbial_shortbred/MAGs/studies/triglycerides.csv")

basic_hdl<-subset(all_basic,is.na(hdl)==F,drop=T)
table(basic_hdl$study_name)
index_shortbred_asso(basic_hdl,"hdl_norm",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/hdl_dsv.csv",
"08.Microbial_shortbred/sig_SVs/studies/hdl_vsv.csv",
"08.Microbial_shortbred/MAGs/studies/hdl.csv")

basic_hba1c<-subset(all_basic,is.na(hba1c)==F,drop=T)
table(basic_hba1c$study_name)
index_shortbred_asso(basic_hba1c,"hba1c_norm",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/hba1c_dsv.csv",
"08.Microbial_shortbred/sig_SVs/studies/hba1c_vsv.csv",
"08.Microbial_shortbred/MAGs/studies/hba1c.csv")

basic_ldl<-subset(all_basic,is.na(ldl)==F,drop=T)
table(basic_ldl$study_name)
index_shortbred_asso(basic_ldl,"ldl_norm",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/ldl_dsv.csv",
"08.Microbial_shortbred/sig_SVs/studies/ldl_vsv.csv",
"08.Microbial_shortbred/MAGs/studies/ldl.csv")

basic_cholesterol<-subset(all_basic,is.na(cholesterol)==F,drop=T)
table(basic_cholesterol$study_name)
index_shortbred_asso(basic_cholesterol,"cholesterol_norm",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/cholesterol_dsv.csv",
"08.Microbial_shortbred/sig_SVs/studies/cholesterol_vsv.csv",
"08.Microbial_shortbred/MAGs/studies/cholesterol.csv")

basic_mean_press<-subset(all_basic,is.na(mean_press)==F,drop=T)
table(basic_mean_press$study_name)
index_shortbred_asso(basic_mean_press,"mean_press_norm",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/mean_press_dsv.csv",
"08.Microbial_shortbred/sig_SVs/studies/mean_press_vsv.csv",
"08.Microbial_shortbred/MAGs/studies/mean_press.csv")

basic_creatinine<-subset(all_basic,is.na(creatinine)==F,drop=T)
table(basic_creatinine$study_name)
index_shortbred_asso(basic_creatinine,"creatinine_norm",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/creatinine_dsv.csv",
"08.Microbial_shortbred/sig_SVs/studies/creatinine_vsv.csv",
"08.Microbial_shortbred/MAGs/studies/creatinine.csv")

basic_ast<-subset(all_basic,is.na(ast)==F,drop=T)
table(basic_ast$study_name)
index_shortbred_asso(basic_ast,"ast_norm",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/ast_dsv.csv",
"08.Microbial_shortbred/sig_SVs/studies/ast_vsv.csv",
"08.Microbial_shortbred/MAGs/studies/ast.csv")

basic_alt<-subset(all_basic,is.na(alt)==F,drop=T)
table(basic_alt$study_name)
index_shortbred_asso(basic_alt,"alt_norm",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/alt_dsv.csv",
"08.Microbial_shortbred/sig_SVs/studies/alt_vsv.csv",
"08.Microbial_shortbred/MAGs/studies/alt.csv")

basic_total_bilirubin<-subset(all_basic,is.na(total_bilirubin)==F,drop=T)
table(basic_total_bilirubin$study_name)
index_shortbred_asso(basic_total_bilirubin,"total_bilirubin_norm",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/total_bilirubin_dsv.csv",
"08.Microbial_shortbred/sig_SVs/studies/total_bilirubin_vsv.csv",
"08.Microbial_shortbred/MAGs/studies/total_bilirubin.csv")

basic_fasting_insulin<-subset(all_basic,is.na(fasting_insulin)==F,drop=T)
table(basic_fasting_insulin$study_name)
index_shortbred_asso(basic_fasting_insulin,"fasting_insulin_norm",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/fasting_insulin_dsv.csv",
"08.Microbial_shortbred/sig_SVs/studies/fasting_insulin_vsv.csv",
"08.Microbial_shortbred/MAGs/studies/fasting_insulin.csv")

basic_fasting_glucose <-subset(all_basic,is.na(fasting_glucose)==F,drop=T)
table(basic_fasting_glucose$study_name)
index_shortbred_asso(basic_fasting_glucose ,"fasting_glucose_norm",c("age","gender_code","BMI"),
"08.Microbial_shortbred/sig_SVs/studies/fasting_glucose_dsv.csv",
"08.Microbial_shortbred/sig_SVs/studies/fasting_glucose_vsv.csv",
"08.Microbial_shortbred/MAGs/studies/fasting_glucose.csv")



#################################################################################################################################
###  meta-analysis  (no adjust abundance) more than one datasets 
################### disease 

binary_multi<-c("CRC","IBD","T2D","adenoma","PD")
continuous<-c(
  "age","BMI","hba1c","fasting_glucose",
  "fasting_insulin","triglycerides","hdl","ldl","cholesterol","mean_press",
  "creatinine","ast","alt","total_bilirubin"
)
binary_single<-c("ACVD","CAD","cirrhosis","MetS")

################### Preparation
meta_results<-data.frame(gene_set=character(),outcome=character(),meta_file=character(),stringsAsFactors=F)
for(gene_set in c("sig_SVs","MAGs")){
  for(out_dir in c("meta","multi","anno","study_wide_sig","analysis_wide_sig")){
    dir.create(paste0("08.Microbial_shortbred/",gene_set,"/",out_dir),recursive=T,showWarnings=F)
  }
}

################### Binary phenotypes: multiple datasets
for(gene_set in c("sig_SVs","MAGs")){
  for(outcome in binary_multi){
    study_file<-paste0("08.Microbial_shortbred/",gene_set,"/studies/",outcome,".csv")
    meta_file<-paste0("08.Microbial_shortbred/",gene_set,"/meta/",outcome,".csv")
    if(!file.exists(study_file)) stop("Missing study input: ",study_file)
    combine_meta_p_binary(read.csv(study_file,header=T),6,meta_file)
    meta_data<-read.csv(meta_file,header=T)
    meta_data$gene_model<-rep("abundance",nrow(meta_data))
    write.csv(meta_data,meta_file,row.names=F)
    meta_results<-rbind(meta_results,data.frame(gene_set=gene_set,outcome=outcome,meta_file=meta_file))
  }
}

###################################################################################################
### Continuous phenotypes
for(gene_set in c("sig_SVs","MAGs")){
  for(outcome in continuous){
    meta_file<-paste0("08.Microbial_shortbred/",gene_set,"/meta/",outcome,".csv")
    if(gene_set=="sig_SVs"){
      study_dsv<-paste0("08.Microbial_shortbred/sig_SVs/studies/",outcome,"_dsv.csv")
      study_vsv<-paste0("08.Microbial_shortbred/sig_SVs/studies/",outcome,"_vsv.csv")
      if(!file.exists(study_dsv)) stop("Missing study input: ",study_dsv)
      if(!file.exists(study_vsv)) stop("Missing study input: ",study_vsv)
      dsv_data<-read.csv(study_dsv,header=T)
      vsv_data<-read.csv(study_vsv,header=T)
      if(!identical(colnames(dsv_data),colnames(vsv_data))) stop("Study columns differ: ",outcome)
      study_data<-rbind(dsv_data,vsv_data)
      # These are different predictor models, even when their gene IDs overlap.
      gene_model<-c(rep("presence",nrow(dsv_data)),rep("abundance",nrow(vsv_data)))
    }else{
      study_file<-paste0("08.Microbial_shortbred/MAGs/studies/",outcome,".csv")
      if(!file.exists(study_file)) stop("Missing study input: ",study_file)
      study_data<-read.csv(study_file,header=T)
      gene_model<-rep("abundance",nrow(study_data))
    }
    # ShortBRED linear regression has six columns per cohort, including CI and N.
    combine_meta_p_lr(study_data,6,meta_file)
    meta_data<-read.csv(meta_file,header=T)
    if(!identical(as.character(meta_data$X),as.character(study_data$X))) stop("Meta gene order changed: ",outcome)
    meta_data$gene_model<-gene_model
    write.csv(meta_data,meta_file,row.names=F)
    meta_results<-rbind(meta_results,data.frame(gene_set=gene_set,outcome=outcome,meta_file=meta_file))
  }
}


########################################################################################################
##### Just one dataset
for(gene_set in c("sig_SVs","MAGs")){
  for(outcome in binary_single){
    file_outcome<-outcome
    if(outcome=="CAD") file_outcome<-"CVD"
    study_file<-paste0("08.Microbial_shortbred/",gene_set,"/studies/",file_outcome,".csv")
    meta_file<-paste0("08.Microbial_shortbred/",gene_set,"/meta/",outcome,".csv")
    if(!file.exists(study_file)) stop("Missing study input: ",study_file)
    meta_data<-read.csv(study_file,header=T)
    if(ncol(meta_data)!=7) stop("Expected one six-column study block: ",study_file)
    valid<-is.finite(meta_data[,2]) & meta_data[,2]>0 & is.finite(meta_data[,3]) & meta_data[,3]>0
    meta_data$OR_meta<-ifelse(valid,meta_data[,2],NA_real_)
    meta_data$se_meta<-ifelse(valid,meta_data[,3],NA_real_)
    meta_data$pvalue_meta<-ifelse(valid,meta_data[,4],NA_real_)
    meta_data$cilb<-ifelse(valid,meta_data[,5],NA_real_)
    meta_data$ciub<-ifelse(valid,meta_data[,6],NA_real_)
    meta_data$I2<-rep(NA_real_,nrow(meta_data))
    meta_data$het_p<-ifelse(valid,1,NA_real_)
    meta_data$p_count<-ifelse(valid,as.integer(meta_data$pvalue_meta<=0.05),NA_integer_)
    meta_data$gene_model<-rep("abundance",nrow(meta_data))
    write.csv(meta_data,meta_file,row.names=F)
    meta_results<-rbind(meta_results,data.frame(gene_set=gene_set,outcome=outcome,meta_file=meta_file))
  }
}


######################################################################################################
### Multiple testing correction
for(i in seq_len(nrow(meta_results))){
  gene_set<-meta_results$gene_set[i]
  outcome<-meta_results$outcome[i]
  meta_file<-meta_results$meta_file[i]
  meta_data<-read.csv(meta_file,header=T)
  if(anyDuplicated(meta_data[,c("X","gene_model"),drop=F])) stop("Duplicate gene-model results: ",meta_file)
  meta_data$study_count<-rep(0L,nrow(meta_data))
  effect_cols<-grep("\\.(OR|beta)$",colnames(meta_data),value=T)
  for(effect_col in effect_cols){
    se_col<-sub("\\.(OR|beta)$",".SE",effect_col)
    valid<-is.finite(meta_data[[effect_col]]) & is.finite(meta_data[[se_col]]) & meta_data[[se_col]]>0
    if(outcome %in% c(binary_multi,binary_single)) valid<-valid & meta_data[[effect_col]]>0
    meta_data$study_count<-meta_data$study_count+as.integer(valid)
  }
  if(any(!is.na(meta_data$pvalue_meta) & (!is.finite(meta_data$pvalue_meta) | meta_data$pvalue_meta<0 | meta_data$pvalue_meta>1))){
    stop("Invalid Meta P values: ",meta_file)
  }
  tested<-which(!is.na(meta_data$pvalue_meta))
  meta_data$bonferroni<-rep(NA_real_,nrow(meta_data))
  meta_data$fdr<-rep(NA_real_,nrow(meta_data))
  meta_data$bonferroni[tested]<-p.adjust(meta_data$pvalue_meta[tested],method="bonferroni")
  meta_data$fdr[tested]<-p.adjust(meta_data$pvalue_meta[tested],method="fdr")
  write.csv(meta_data,paste0("08.Microbial_shortbred/",gene_set,"/multi/",outcome,".csv"),row.names=F)
}


######################################################################################################
### Significant p values

###################### Inputs
gene_info<-read.csv("00.rawData/shortBRED/MAGs/HRGM2_select_gene_info_0.15_add_ko.csv",header=T,stringsAsFactors=F)
id_mapping<-read.table("00.rawData/shortBRED/SV_sig_HRGM2_0.15_mapping.txt",header=T,stringsAsFactors=F)
colnames(id_mapping)<-c("New_Gene_Name","Adjusted_Gene_Name")
gene_info<-merge(gene_info,id_mapping,by.x="New_Gene_Name",by.y="New_Gene_Name")

# Combine annotations for duplicate IDs without multiplying association rows.
gene_info<-unique(gene_info[!is.na(gene_info$Adjusted_Gene_Name) & nzchar(gene_info$Adjusted_Gene_Name),,drop=F])
duplicate_ids<-unique(gene_info$Adjusted_Gene_Name[duplicated(gene_info$Adjusted_Gene_Name)])
for(gene_id in duplicate_ids){
  rows<-which(gene_info$Adjusted_Gene_Name==gene_id)
  for(column in setdiff(colnames(gene_info),"Adjusted_Gene_Name")){
    values<-unique(as.character(gene_info[rows,column]))
    values<-values[!is.na(values) & nzchar(values)]
    gene_info[rows[1],column]<-if(length(values)>0) paste(values,collapse="; ") else NA_character_
  }
}
gene_info<-gene_info[!duplicated(gene_info$Adjusted_Gene_Name),,drop=F]

study_wide_cutoff<-4.42e-5
analysis_wide_cutoff<-0.05
I2_cutoff<-75

for(i in seq_len(nrow(meta_results))){
  gene_set<-meta_results$gene_set[i]
  outcome<-meta_results$outcome[i]
  meta_pvalue<-read.csv(paste0("08.Microbial_shortbred/",gene_set,"/multi/",outcome,".csv"),header=T)
  heterogeneity_ok<-(!is.na(meta_pvalue$I2) & meta_pvalue$I2<I2_cutoff) | meta_pvalue$study_count==1

  meta_pvalue.sig<-meta_pvalue[which(!is.na(meta_pvalue$pvalue_meta) & meta_pvalue$pvalue_meta<study_wide_cutoff & heterogeneity_ok),,drop=F]
  meta_pvalue.sig$annotation_matched<-meta_pvalue.sig$X %in% gene_info$Adjusted_Gene_Name
  meta_pvalue.sig<-left_join(meta_pvalue.sig,gene_info,by=c("X"="Adjusted_Gene_Name"))
  write.csv(meta_pvalue.sig,paste0("08.Microbial_shortbred/",gene_set,"/study_wide_sig/",outcome,".sig.anno.csv"),row.names=F)

  meta_pvalue.sig<-meta_pvalue[which(!is.na(meta_pvalue$bonferroni) & meta_pvalue$bonferroni<analysis_wide_cutoff & heterogeneity_ok),,drop=F]
  meta_pvalue.sig$annotation_matched<-meta_pvalue.sig$X %in% gene_info$Adjusted_Gene_Name
  meta_pvalue.sig<-left_join(meta_pvalue.sig,gene_info,by=c("X"="Adjusted_Gene_Name"))
  write.csv(meta_pvalue.sig,paste0("08.Microbial_shortbred/",gene_set,"/analysis_wide_sig/",outcome,".sig.anno.csv"),row.names=F)
}



########################################################################################
## Annotation of all association results
for(i in seq_len(nrow(meta_results))){
  gene_set<-meta_results$gene_set[i]
  outcome<-meta_results$outcome[i]
  gene_pvalue<-read.csv(paste0("08.Microbial_shortbred/",gene_set,"/multi/",outcome,".csv"),header=T)
  gene_pvalue$annotation_matched<-gene_pvalue$X %in% gene_info$Adjusted_Gene_Name
  gene_pvalue<-left_join(gene_pvalue,gene_info,by=c("X"="Adjusted_Gene_Name"))
  write.csv(gene_pvalue,paste0("08.Microbial_shortbred/",gene_set,"/anno/",outcome,".csv"),row.names=F)
  if(nrow(gene_pvalue)>0 && !any(gene_pvalue$annotation_matched)){
    warning("No gene_info matches for ",gene_set,"/",outcome,"; results retained without annotations")
  }
}



########################################################################################
### BMI: MAG gene association volcano plot
library(ggplot2)
BMI_genes<-read.csv("08.Microbial_shortbred/MAGs/anno/BMI.csv",check.names=FALSE)
stopifnot(all(c("X","effect_meta","pvalue_meta","bonferroni","I2","Orthology","Annotation") %in% names(BMI_genes)),
          !anyDuplicated(BMI_genes$X))
BMI_genes$significant<-is.finite(BMI_genes$bonferroni) & BMI_genes$bonferroni<0.05 &
  is.finite(BMI_genes$I2) & BMI_genes$I2<75
BMI_genes$direction<-ifelse(BMI_genes$effect_meta>0,"Upregulated",
                           ifelse(BMI_genes$effect_meta<0,"Downregulated","Zero effect"))
BMI_genes$species<-sub("_gene.*$","",BMI_genes$X)
BMI_genes$species<-gsub("_",".",BMI_genes$species)
BMI_genes$species[BMI_genes$species=="Oscillibacter.sp"]<-"Oscillibacter sp. ER4"
BMI_genes$category<-ifelse(BMI_genes$significant,BMI_genes$direction,"Not significant")
BMI_genes$category<-factor(BMI_genes$category,levels=c("Upregulated","Downregulated","Not significant","Zero effect"))
BMI_plot<-subset(BMI_genes,is.finite(effect_meta) & is.finite(bonferroni))
stopifnot(nrow(BMI_plot)>0L,all(BMI_plot$bonferroni>=0 & BMI_plot$bonferroni<=1))
BMI_plot<-BMI_plot[order(BMI_plot$significant),]
BMI_sig<-subset(BMI_genes,significant)
stopifnot(nrow(BMI_sig)>0L,all(is.finite(BMI_sig$effect_meta)))
BMI_colors<-setNames(c(mycolor3_red_blue_yellow[1:2],mycolor2_pink_grey[2],"grey40"),levels(BMI_genes$category))
dir.create("pics",showWarnings=FALSE)
BMI_out<-"08.Microbial_shortbred/pic_files/BMI_MAG_Bonferroni"
dir.create(BMI_out,recursive=TRUE,showWarnings=FALSE)
write.csv(BMI_sig,file.path(BMI_out,"significant_genes.csv"),row.names=FALSE)

BMI_volcano<-ggplot(BMI_plot,aes(effect_meta,-log10(pmax(bonferroni,.Machine$double.xmin))))+
  geom_hline(yintercept=-log10(0.05),linetype="dashed",colour="grey50",linewidth=0.4)+
  geom_vline(xintercept=0,linetype="dotted",colour="grey50",linewidth=0.4)+
  geom_point(aes(colour=category),size=1.3,alpha=0.75)+
  scale_colour_manual(values=BMI_colors)+
  labs(x="Meta-analysis effect (beta)",y=expression(-log[10](P[Bonferroni])),
       title="BMI: MAG gene associations",
       subtitle=paste0(nrow(BMI_sig)," significant genes\nBonferroni P < 0.05 and I\u00b2 < 75%"),colour=NULL)+
  guides(colour=guide_legend(nrow=2,byrow=TRUE))+
  theme_classic(base_size=25)+
  theme(axis.text=element_text(size=22),axis.title=element_text(size=25),
    legend.text=element_text(size=22),legend.title=element_text(size=25),
    legend.position="bottom",plot.title=element_text(size=25),
    plot.subtitle=element_text(size=21))
# Label the requested gene using its functional annotation.
BMI_label<-BMI_plot[BMI_plot$X=="Oscillibacter_sp_gene204",,drop=FALSE]
stopifnot(nrow(BMI_label)==1L,!is.na(BMI_label$Annotation),nzchar(BMI_label$Annotation))
BMI_label$label<-paste(strwrap(BMI_label$Annotation,width=23),collapse="\n")
BMI_label$label_x<-BMI_label$effect_meta-0.03*diff(range(BMI_plot$effect_meta))
BMI_label$label_y<--log10(pmax(BMI_label$bonferroni,.Machine$double.xmin))+
  0.12*max(-log10(pmax(BMI_plot$bonferroni,.Machine$double.xmin)))
BMI_volcano<-BMI_volcano+
  geom_point(data=BMI_label,aes(colour=category),size=2.2,show.legend=FALSE)+
  geom_segment(data=BMI_label,
    aes(x=effect_meta,y=-log10(pmax(bonferroni,.Machine$double.xmin)),
        xend=label_x,yend=label_y,colour=category),
    linewidth=0.35,inherit.aes=FALSE,show.legend=FALSE)+
  geom_label(data=BMI_label,aes(label_x,label_y,label=label,colour=category),
    hjust=1,size=3.5+10/ggplot2::.pt,fill="white",label.size=0.25,
    label.padding=grid::unit(0.15,"lines"),lineheight=1.1,
    inherit.aes=FALSE,show.legend=FALSE)
write.csv(BMI_label,file.path(BMI_out,"volcano_gene_label.csv"),row.names=FALSE)
saveRDS(BMI_volcano,file.path(BMI_out,"volcano.plot.rds"))
ragg::agg_tiff("pics/FS10A_BMI_MAG_gene_volcano.tiff",
               width=11,height=12,units="in",res=300,compression="lzw",background="white")
print(BMI_volcano);dev.off()

## Show gene-level significance by species; order by the smallest adjusted P.
BMI_species<-as.data.frame(table(BMI_sig$species),stringsAsFactors=FALSE)
names(BMI_species)<-c("species","genes")
BMI_species$min_bonferroni<-vapply(BMI_species$species,
  function(x) min(BMI_sig$bonferroni[BMI_sig$species==x]),numeric(1))
BMI_species<-BMI_species[order(BMI_species$min_bonferroni,BMI_species$species),]
stopifnot(sum(BMI_species$genes)==nrow(BMI_sig),!anyDuplicated(BMI_sig$X))
write.csv(BMI_species,file.path(BMI_out,"species_counts.csv"),row.names=FALSE)
BMI_species_direction<-as.data.frame(table(species=BMI_sig$species,direction=BMI_sig$direction))
BMI_species_direction$species<-factor(BMI_species_direction$species,levels=rev(BMI_species$species))
write.csv(BMI_species_direction,file.path(BMI_out,"species_direction_counts.csv"),row.names=FALSE)

BMI_species_points<-BMI_sig
BMI_species_points$species<-factor(BMI_species_points$species,levels=rev(BMI_species$species))
BMI_species_points$neg_log10_bonferroni<--log10(pmax(BMI_species_points$bonferroni,.Machine$double.xmin))
BMI_species_labels<-setNames(sprintf("italic('%s')~~'(n = %d)'",
  BMI_species$species,BMI_species$genes),BMI_species$species)
stopifnot(all(BMI_species_points$direction %in% c("Upregulated","Downregulated")))
BMI_species_plot<-ggplot(BMI_species_points,aes(neg_log10_bonferroni,species))+
  geom_vline(xintercept=-log10(0.05),linetype="dashed",colour="grey60",linewidth=0.35)+
  geom_point(aes(colour=direction),size=1.6,alpha=0.7,
    position=position_jitter(width=0,height=0.18,seed=20261027))+
  scale_colour_manual(values=c(Upregulated=mycolor3_red_blue_yellow[1],
    Downregulated=mycolor3_red_blue_yellow[2]),
    breaks=c("Upregulated","Downregulated"),labels=c("Positive","Negative"),name=NULL)+
  scale_y_discrete(labels=parse(text=BMI_species_labels[levels(BMI_species_points$species)]))+
  scale_x_continuous(limits=c(0,NA),expand=expansion(mult=c(0,0.04)))+
  labs(x=expression(-log[10](P[Bonferroni])),y=NULL,
    title="BMI-associated MAG genes by species",
    subtitle=paste0(nrow(BMI_sig)," significant genes; Bonferroni P < 0.05 and I\u00b2 < 75%"))+
  guides(colour=guide_legend(override.aes=list(size=3,alpha=1)))+
  theme_classic(base_size=20,base_family="Arial")+
  theme(axis.text=element_text(size=14*0.8+6),
    axis.title=element_text(size=20),legend.text=element_text(size=16.6),
    plot.title=element_text(size=14*1.2+6),plot.subtitle=element_text(size=16),
    panel.grid.major.y=element_line(colour="grey93",linewidth=0.3),
    plot.title.position="plot",legend.position="bottom")
write.csv(BMI_species_points,file.path(BMI_out,"species_significance_points.csv"),row.names=FALSE)
saveRDS(BMI_species_plot,file.path(BMI_out,"species_significance.plot.rds"))
ragg::agg_tiff("pics/FS10B_BMI_MAG_gene_species.tiff",
  width=9,height=12,units="in",res=300,compression="lzw",background="white")
print(BMI_species_plot);dev.off()

## Assign significant genes to microbial KEGG BRITE L1 categories.
kegg_class<-read.csv("00.rawData/KEGG/KEGG_classification.csv",stringsAsFactors=FALSE)
kegg_class<-unique(kegg_class[grepl("^091(00|20|30|40) ",kegg_class$L1),c("KO","L1")])
kegg_class$category<-sub("^[0-9]+ ","",kegg_class$L1)
gene_KO<-data.frame(X=character(),KO=character())
for(i in seq_len(nrow(BMI_sig))){
  if(is.na(BMI_sig$Orthology[i])) next
  ko<-unique(regmatches(BMI_sig$Orthology[i],gregexpr("K[0-9]{5}",BMI_sig$Orthology[i]))[[1]])
  if(length(ko)) gene_KO<-rbind(gene_KO,data.frame(X=BMI_sig$X[i],KO=ko))
}
gene_KEGG<-unique(merge(gene_KO,kegg_class[,c("KO","category")],by="KO"))
BMI_sig$functional_class<-"Unknown / unannotated"
for(i in seq_len(nrow(BMI_sig))){
  terms<-unique(gene_KEGG$category[gene_KEGG$X==BMI_sig$X[i]])
  if(BMI_sig$X[i] %in% gene_KO$X) BMI_sig$functional_class[i]<-"KO annotated, no functional class"
  if(length(terms)==1L) BMI_sig$functional_class[i]<-terms
  if(length(terms)>1L) BMI_sig$functional_class[i]<-"Multiple KEGG categories"
}
write.csv(gene_KEGG,file.path(BMI_out,"gene_KEGG_mapping.csv"),row.names=FALSE)
write.csv(BMI_sig,file.path(BMI_out,"gene_functional_classes.csv"),row.names=FALSE)

## Ring on the left, functional classes and counts on the right.
KEGG_panels<-data.frame(prefix=c("FS10C","FS10D","FS10E","FS10F"),
  scope=c("ER4","ER4","All","All"),
  direction=c("Upregulated","Downregulated","Upregulated","Downregulated"))
KEGG_categories<-c("Metabolism","Genetic Information Processing",
  "Environmental Information Processing","Cellular Processes","Multiple KEGG categories",
  "KO annotated, no functional class","Unknown / unannotated")
KEGG_colors<-setNames(c("#207F4C","#ECD452","#2983BB","#C72F1E",
  "#EA8958","#C4C4C4","#696969"),KEGG_categories)
KEGG_labels<-c("Metabolism","Genetic information\nprocessing",
  "Environmental information\nprocessing","Cellular processes","Multiple KEGG categories",
  "KO annotated, no\nfunctional class","Unknown / unannotated")
stopifnot(all(BMI_sig$functional_class %in% KEGG_categories),!anyDuplicated(BMI_sig$X))
for(i in seq_len(nrow(KEGG_panels))){
  panel<-KEGG_panels[i,]
  genes<-subset(BMI_sig,direction==panel$direction)
  if(panel$scope=="ER4") genes<-subset(genes,species=="Oscillibacter sp. ER4")
  counts<-data.frame(category=factor(KEGG_categories,levels=KEGG_categories),
    genes=as.integer(table(factor(genes$functional_class,levels=KEGG_categories))))
  counts$percentage<-100*counts$genes/max(1,nrow(genes))
  counts$ymax<-cumsum(counts$genes);counts$ymin<-counts$ymax-counts$genes
  counts$ymid<-(counts$ymin+counts$ymax)/2
  counts$text_colour<-ifelse(KEGG_categories %in% c("Genetic Information Processing",
    "Multiple KEGG categories","KO annotated, no functional class"),"black","white")
  stopifnot(sum(counts$genes)==nrow(genes))
  write.csv(genes,file.path(BMI_out,paste0(panel$prefix,"_genes.csv")),row.names=FALSE)
  write.csv(counts[,c("category","genes","percentage")],
    file.path(BMI_out,paste0(panel$prefix,"_KEGG_composition.csv")),row.names=FALSE)
  legend_labels<-paste0(KEGG_labels,"\n",counts$genes," (",sprintf("%.1f",counts$percentage),"%)")
  positive<-panel$direction=="Upregulated"
  ring<-ggplot(counts)+
    geom_rect(aes(xmin=1.1,xmax=2.9,ymin=ymin,ymax=ymax,fill=category),
      colour="white",linewidth=0.35,key_glyph="point")+
    geom_text(data=counts[counts$percentage>=5,],
      aes(x=2.1,y=ymid,label=sprintf("%.1f%%",percentage),colour=text_colour),size=4.5)+
    annotate("text",x=0,y=0,label=paste0(nrow(genes),"\nMAG genes"),
      colour=if(positive) "#BC3523" else "#2983BB",fontface="bold",size=5,lineheight=1.3)+
    scale_x_continuous(limits=c(0,2.9),expand=c(0,0))+
    scale_y_continuous(limits=c(0,max(1,nrow(genes))),expand=c(0,0))+
    scale_fill_manual(values=KEGG_colors,breaks=KEGG_categories,labels=legend_labels,drop=FALSE)+
    scale_colour_identity()+coord_polar(theta="y",start=0)+
    guides(fill=guide_legend(override.aes=list(shape=22,size=5,colour=NA,stroke=0)))+
    labs(title=paste0(if(panel$scope=="ER4") "Oscillibacter sp. ER4" else "BMI",
      ": ",if(positive) "positive" else "negative"," associations"),
      subtitle=paste0("Meta beta ",if(positive) ">" else "<",
        " 0; Bonferroni P < 0.05; I\u00b2 < 75%"),fill=NULL)+
    theme_void(base_size=13,base_family="Arial")+
    theme(plot.title=element_text(size=18,face="bold",hjust=0.5,margin=margin(b=5)),
      plot.subtitle=element_text(size=14,hjust=0.5,margin=margin(b=16)),
      legend.position="right",legend.text=element_text(size=12,margin=margin(t=3,b=3)),
      legend.key.size=grid::unit(0.42,"cm"),legend.spacing.y=grid::unit(0.12,"cm"),
      plot.margin=margin(8,6,8,6),legend.box.spacing=grid::unit(3,"pt"),
      legend.margin=margin(0,0,0,0))
  prefix<-paste0(panel$prefix,"_BMI_MAG_",panel$scope,"_",tolower(panel$direction),"_KEGG")
  ragg::agg_tiff(file.path("pics",paste0(prefix,".tiff")),
    width=8.8,height=7,units="in",res=300,compression="lzw")
  print(ring);dev.off()
  saveRDS(ring,file.path(BMI_out,paste0(prefix,".plot.rds")))
  cat(panel$prefix,panel$scope,panel$direction,nrow(genes),"genes; KEGG BRITE L1\n")
}

#######################################################################################
## A. hadrus gene-presence and abundance associations with BMI 



SV_sig_genes_analysis_wide_model1<-read.csv("07.Microbial_GWAS/models_sig_genes/SV_sig_genes_analysis_wide_model1.csv",header=T)

AH_1577<-subset(SV_sig_genes_analysis_wide_model1,X=="Anaerostipes hadrus SSC/2:1577_1579 and 9 segments" & phenotype=="BMI_sv",drop=FALSE)

AH_1577$gene_species<-paste(AH_1577$genes,AH_1577$Short_name,sep="_")
AH_1577$gene_species<-gsub("\\.", "_", AH_1577$gene_species)

for (model in c("presence","abundance")) {
BMI_shortbred<-read.csv("08.Microbial_shortbred/sig_SVs/meta/BMI.csv",header=T)
# Plot presence and abundance effects separately.
stopifnot("gene_model" %in% names(BMI_shortbred))
BMI_shortbred<-subset(BMI_shortbred,gene_model==model)
colnames(BMI_shortbred)[1]<-"gene_species"
stopifnot(nrow(AH_1577)>0L,!anyDuplicated(AH_1577$gene_species),
          !anyDuplicated(BMI_shortbred$gene_species),
          all(AH_1577$gene_species %in% BMI_shortbred$gene_species))

AH_1577_pic_file<-merge(AH_1577,BMI_shortbred,by.x="gene_species",by.y="gene_species")[,c("genes","effect_meta.y","pvalue_meta.y","se_meta")]

# Match locus IDs to gene symbols and functional names; retain the original association results.
AH_1577_annotation<-read.csv("00.rawData/A.H/dSV_gene.csv",stringsAsFactors=FALSE)
stopifnot(all(c("order","locas","gene","product") %in% names(AH_1577_annotation)))
AH_1577_annotation$locas<-trimws(AH_1577_annotation$locas)
stopifnot(!anyNA(AH_1577_annotation$locas),all(nzchar(AH_1577_annotation$locas)),
          !anyDuplicated(AH_1577_annotation$locas),all(is.finite(AH_1577_annotation$order)))
AH_1577_annotation_index<-match(AH_1577_pic_file$genes,AH_1577_annotation$locas)
stopifnot(!anyNA(AH_1577_annotation_index))
AH_1577_pic_file$gene_symbol<-trimws(ifelse(is.na(AH_1577_annotation$gene[AH_1577_annotation_index]),"",AH_1577_annotation$gene[AH_1577_annotation_index]))
AH_1577_pic_file$product<-AH_1577_annotation$product[AH_1577_annotation_index]
AH_1577_pic_file$gene_order<-AH_1577_annotation$order[AH_1577_annotation_index]
stopifnot(!anyNA(AH_1577_pic_file$product),all(nzchar(trimws(AH_1577_pic_file$product))))
AH_1577_pic_file<-AH_1577_pic_file[order(AH_1577_pic_file$gene_order),,drop=FALSE]

###################### AH_1577: BMI meta-analysis forest plot
library(ggplot2)
stopifnot(nrow(AH_1577_pic_file)>0,
          all(is.finite(AH_1577_pic_file$effect_meta.y)),
          all(is.finite(AH_1577_pic_file$se_meta)),
          all(AH_1577_pic_file$se_meta>=0),
          all(is.finite(AH_1577_pic_file$pvalue_meta.y)),
          all(AH_1577_pic_file$pvalue_meta.y>=0 & AH_1577_pic_file$pvalue_meta.y<=1))
AH_1577_forest<-AH_1577_pic_file
AH_1577_forest$gene_label<-ifelse(nzchar(AH_1577_forest$gene_symbol),paste(AH_1577_forest$gene_symbol,AH_1577_forest$product,sep=" | "),AH_1577_forest$product)
AH_1577_forest$ci_low<-AH_1577_forest$effect_meta.y-1.96*AH_1577_forest$se_meta
AH_1577_forest$ci_high<-AH_1577_forest$effect_meta.y+1.96*AH_1577_forest$se_meta
AH_1577_forest$meta_p<-ifelse(AH_1577_forest$pvalue_meta.y<0.001,
                            formatC(AH_1577_forest$pvalue_meta.y,format="e",digits=2),
                            formatC(AH_1577_forest$pvalue_meta.y,format="f",digits=3))
AH_1577_forest$gene_row<-factor(seq_len(nrow(AH_1577_forest)),
                              levels=rev(seq_len(nrow(AH_1577_forest))))
AH_1577_xlim<-range(c(AH_1577_forest$ci_low,AH_1577_forest$ci_high,0))
AH_1577_xlim<-AH_1577_xlim+c(-1,1)*diff(AH_1577_xlim)*0.05
AH_1577_p_x<-AH_1577_xlim[2]+diff(AH_1577_xlim)*0.04
AH_1577_forest_plot<-ggplot(AH_1577_forest,aes(y=gene_row))+
  geom_vline(xintercept=0,linetype="dashed",colour="grey60",linewidth=0.35)+
  geom_segment(aes(x=ci_low,xend=ci_high,yend=gene_row),linewidth=0.4,colour="#2983BB")+
  geom_point(aes(x=effect_meta.y),size=1.6,colour="#2983BB")+
  geom_text(aes(label=meta_p),x=AH_1577_p_x,hjust=0,size=8/ggplot2::.pt)+
  annotate("text",x=AH_1577_p_x,y=nrow(AH_1577_forest)+1,label="Meta P",
           hjust=0,size=8/ggplot2::.pt,fontface="bold")+
  scale_y_discrete(labels=rev(AH_1577_forest$gene_label),expand=expansion(add=c(0.6,1.5)))+
  scale_x_continuous(expand=expansion(mult=0))+
  coord_cartesian(xlim=AH_1577_xlim,clip="off")+
  labs(title=expression(bolditalic("Anaerostipes hadrus")~bold("SSC/2:1577_1579 and 9 segments")),
       x="Meta-analysis effect (beta)",y=NULL)+
  theme_classic(base_size=9,base_family="sans")+
  theme(axis.text=element_text(size=7.4),
        axis.text.y=element_text(size=8,colour="black"),
        plot.title=element_text(size=10.6,face="bold"),plot.title.position="plot",
        plot.margin=margin(5.5,65,5.5,5.5))
dir.create("pics",showWarnings=FALSE)
forest_file<-if (model=="presence") "pics/FS8_AH.forest.tiff" else "pics/FS8_AH.abundance.forest.tiff"
ggsave(forest_file,AH_1577_forest_plot,
       device="tiff",width=10,height=7,units="in",dpi=300,compression="lzw")
}



#######################################################################################
## A. hadrus carA: BMI by zero versus positive ShortBRED abundance.
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(patchwork)

carA_annotation <- read.csv("00.rawData/A.H/dSV_gene.csv", stringsAsFactors = FALSE)
carA_id <- paste0(carA_annotation$locas[carA_annotation$gene == "carA" & !is.na(carA_annotation$gene)], "_A_hadrus")
stopifnot(length(carA_id) == 1L)
carA_raw <- fread("00.rawData/shortBRED/sig_SVs/matched_rows.csv", check.names = FALSE)
carA_raw <- carA_raw[Family == carA_id]
stopifnot(nrow(carA_raw) == 1L)
carA_abundance <- data.frame(match_id = names(carA_raw)[-1], abundance = as.numeric(carA_raw[1,-1]))
carA_clinical <- read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv", header = TRUE)
carA_clinical$match_id <- row.names(carA_clinical)
carA_data <- merge(carA_clinical, carA_abundance, by = "match_id")
# Reproduce the existing BMI presence model: all participants, adjusting for sex.
carA_data <- subset(carA_data, study_name != "PiccinnoG_2025_cohort6" &
                     is.finite(BMI) & is.finite(BMI_norm) & is.finite(gender_code) & is.finite(abundance))
stopifnot(all(carA_data$abundance >= 0))
carA_data$temp <- as.integer(carA_data$abundance > 0)
carA_data$study_name <- normalize_study_name(carA_data$study_name)
carA_results <- read.csv("08.Microbial_shortbred/sig_SVs/studies/BMI_dsv.csv", check.names = FALSE)
names(carA_results) <- normalize_study_name(names(carA_results))
result_row <- which(carA_results$X == carA_id)
stopifnot(length(result_row) == 1L)
# Retain all cohorts with more than 100 model samples, ordered by sample size.
sample_columns <- grep("[.]N$", names(carA_results), value = TRUE)
study_n <- setNames(as.numeric(carA_results[result_row, sample_columns]),
                    sub("[.]N$", "", sample_columns))
study_n <- sort(study_n[is.finite(study_n) & study_n > 100L], decreasing = TRUE)
study_n <- study_n[!names(study_n) %in% c("LiJ_2014", "QinJ_2012")]
study_order <- names(study_n)
if (!length(study_order)) stop("No cohorts have more than 100 association samples.")
carA_data <- subset(carA_data, study_name %in% study_order, drop = TRUE)
stopifnot(all(as.numeric(table(carA_data$study_name)[study_order]) == study_n))
study_labels <- setNames(gsub("_", " ", study_order), study_order)
study_labels["PiccinnoG_2025_cohort5"] <- "PiccinnoG 2025\ncohort 5"
study_labels["Fern\u00e1ndez_Calvet_2024"] <- "Fern\u00e1ndez Calvet\n2024"
study_palette <- c("0" = "#7876B1FF", "1" = "#207F4C")

plot_list <- vector("list", length(study_order))
names(plot_list) <- study_order
delta_rows <- vector("list", length(study_order))
names(delta_rows) <- study_order
bmi_limits <- c(min(15, min(carA_data$BMI, na.rm = TRUE)), max(carA_data$BMI, na.rm = TRUE))
bmi_breaks <- seq(15, max(45, floor(max(bmi_limits) / 15) * 15), by = 15)

for (study in study_order) {
  plot_data <- subset(carA_data, study_name == study, drop = TRUE)
  n_del <- sum(plot_data$temp == 0)
  n_no_del <- sum(plot_data$temp == 1)
  median_del <- median(plot_data$BMI[plot_data$temp == 0], na.rm = TRUE)
  median_no_del <- median(plot_data$BMI[plot_data$temp == 1], na.rm = TRUE)
  delta_median <- median_del - median_no_del
  delta_rows[[study]] <- data.frame(
    study_name = study,
    n_deletion = n_del,
    n_no_deletion = n_no_del,
    median_BMI_deletion = median_del,
    median_BMI_no_deletion = median_no_del,
    delta_median_BMI_deletion_minus_no_deletion = delta_median,
    stringsAsFactors = FALSE
  )

  plot_data$SV_group <- factor(
    as.character(plot_data$temp),
    levels = c("0", "1"),
    labels = c(
      paste0("Abundance = 0\n(n=", n_del, ")"),
      paste0("Abundance > 0\n(n=", n_no_del, ")")
    )
  )
  plot_list[[study]] <- ggplot(
    plot_data,
    aes(x = SV_group, y = BMI, color = as.character(temp))
  ) +
    geom_hline(yintercept = bmi_breaks, linetype = "dashed", color = "grey65", linewidth = 0.8) +
    geom_violin(alpha = 0.05, width = 0.9, linewidth = 1.1, trim = TRUE) +
    ggbeeswarm::geom_quasirandom(width = 0.18, size = 0.7, alpha = 0.55) +
    geom_boxplot(
      alpha = 0.18, width = 0.24, color = "black", linewidth = 0.9,
      outlier.shape = NA
    ) +
    scale_color_manual(values = study_palette, guide = "none") +
    scale_y_continuous(limits = bmi_limits, breaks = bmi_breaks) +
    labs(
      title = paste0(unname(study_labels[study]), "\nDelta = ", ifelse(is.finite(delta_median), sprintf("%+.2f", delta_median), "not estimable")),
      x = NULL, y = if (study %in% study_order[c(1,8,15,22)]) "BMI" else NULL
    ) +
    theme_cowplot() +
    theme(
      panel.grid = element_blank(), legend.position = "none",
      axis.line = element_line(linewidth = 0.9), axis.ticks = element_line(linewidth = 0.8),
      text = element_text(size = 20), axis.text = element_text(size = 19),
      axis.title = element_text(size = 20),
      axis.title.y = if (study %in% study_order[c(1,8,15,22)]) element_text(size = 20) else element_blank(),
      plot.title = element_text(size = 21, face = "plain", hjust = 0.5, vjust = 0,
                                margin = margin(t = 12, b = 6)),
      axis.text.y = if (study %in% study_order[c(1,8,15,22)]) element_text(size = 19) else element_blank(),
      axis.ticks.y = if (study %in% study_order[c(1,8,15,22)]) element_line() else element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = unit(c(2, 2, 2, 2), "pt")
    )
}

delta_summary <- do.call(rbind, delta_rows)
delta_summary$study_label <- unname(study_labels[delta_summary$study_name])

beta_rows <- lapply(study_order, function(study) {
  beta_col <- paste0(study, ".beta")
  se_col <- paste0(study, ".SE")
  if (!all(c(beta_col, se_col) %in% names(carA_results))) {
    stop("Missing beta or SE columns for study: ", study)
  }
  beta <- as.numeric(carA_results[result_row, beta_col])
  se <- as.numeric(carA_results[result_row, se_col])
  data.frame(
    study_name = study,
    beta = beta,
    SE = se,
    CI_low = beta - 1.96 * se,
    CI_high = beta + 1.96 * se,
    stringsAsFactors = FALSE
  )
})
beta_summary <- do.call(rbind, beta_rows)
beta_summary$study_label <- gsub("\n", " ", unname(study_labels[beta_summary$study_name]))
beta_summary$study_label <- factor(
  beta_summary$study_label,
  levels = rev(gsub("\n", " ", unname(study_labels[study_order])))
)
beta_summary$beta_CI_label <- sprintf(
  "%+.2f (%.2f, %.2f)",
  beta_summary$beta, beta_summary$CI_low, beta_summary$CI_high
)

beta_summary$beta_CI_label[!is.finite(beta_summary$beta) | !is.finite(beta_summary$SE)] <- "Not estimated"

forest_xlim <- range(c(beta_summary$CI_low, beta_summary$CI_high), na.rm = TRUE)
forest_span <- diff(forest_xlim)
if (!is.finite(forest_span) || forest_span == 0) forest_span <- 1
forest_xlim <- forest_xlim + c(-0.08, 0.08) * forest_span
label_x <- forest_xlim[2] + 0.12 * diff(forest_xlim)
forest_right <- label_x + 2 * diff(forest_xlim)
forest_ticks <- pretty(forest_xlim, n = 2)
forest_ticks <- forest_ticks[forest_ticks >= forest_xlim[1] & forest_ticks <= forest_xlim[2]]

p_forest <- ggplot(beta_summary, aes(y = study_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.8) +
  geom_segment(
    aes(x = CI_low, xend = CI_high, yend = study_label),
    linewidth = 1.1, color = "#3A3A3A"
  ) +
  geom_point(aes(x = beta), size = 3, color = "#2166AC") +
  geom_text(
    aes(label = beta_CI_label), x = label_x,
    hjust = 0, size = 4 + 3 / ggplot2::.pt, lineheight = 0.95, color = "black"
  ) +
  scale_x_continuous(limits = c(forest_xlim[1], forest_right),
                     breaks = forest_ticks) +
  labs(x = "Rank-normalized BMI (beta)", y = NULL) +
  theme_cowplot() +
  theme(
    panel.grid = element_blank(), axis.text = element_text(size = 15),
    axis.line = element_line(linewidth = 0.9), axis.ticks = element_line(linewidth = 0.8),
    axis.title = element_text(size = 23),
    plot.margin = unit(c(5, 12, 5, 5), "pt")
  )

# Arrange cohort plots in four rows and seven columns, with the forest on the right.
p_cohorts <- patchwork::wrap_plots(plot_list, nrow = 4, ncol = 7)
p_whole <- patchwork::wrap_plots(p_cohorts, p_forest, nrow = 1, widths = c(7, 2.2))
f6d_width <- 30
if (!dir.exists("pics/sites")) dir.create("pics/sites", recursive = TRUE)

write.csv(
  merge(delta_summary, beta_summary[, c("study_name", "beta", "SE", "CI_low", "CI_high")],
        by = "study_name", sort = FALSE),
  file = "08.Microbial_shortbred/pic_files/carA_BMI_effect_summary.csv",
  row.names = FALSE
)

tiff(
  file = "pics/FS8_carA_BMI_presence_effects.tiff",
  width = f6d_width, height = 24, units = "in", res = 300, compression = "lzw"
)
print(p_whole)
dev.off()

