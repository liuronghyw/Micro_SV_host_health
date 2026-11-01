### Microbial GWA 
### 2026-10-20
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
library("ggplot2")
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


# Default: reuse completed results and regenerate annotations and FS9 figures.
# Set either switch to TRUE only when those analyses need to be rerun.
run_associations <- FALSE
run_meta_analysis <- FALSE
if (run_associations) {
###################### 1.2 Inputs
all_basic<-read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv")

SVs_humann3<-fread("00.rawData/humann3/SV_no_species_merged_genefamilies.tsv",header=T)

colnames(SVs_humann3) <- gsub("_kneaddata-RPKs", "", colnames(SVs_humann3))

SVs_humann3_t<-t(SVs_humann3[,-1])
colnames(SVs_humann3_t)<-SVs_humann3[,1]$"Gene Family"


## 2 Associations between genes and prognosis
### 2.1 Preparation
if (!dir.exists("10.Microbial_humann3")) {dir.create("10.Microbial_humann3")}
if (!dir.exists("10.Microbial_humann3/studies")) {dir.create("10.Microbial_humann3/studies")}
if (!dir.exists("10.Microbial_humann3/meta")) {dir.create("10.Microbial_humann3/meta")}


### disease
### IBD T2D CRC CAD
### ACVD adenoma cirrhosis MetS PD

disease_humann3_asso <- function(basic, index_name, adjust_cov, out_humann3) {
  # basic <- basic_MetS
  # index_name <- "MetS"
  # adjust_cov <- c("age", "gender_code", "BMI", "read_count")

  prog <- data.frame(basic[, index_name])
  row.names(prog) <- rownames(basic)
  colnames(prog) <- "Y"

  common_rows <- intersect(rownames(basic), rownames(SVs_humann3_t))
  #SVs_humann3 <- SVs_humann3_t[rownames(basic), ]
  SVs_humann3 <- SVs_humann3_t[common_rows, ]

  SVs_humann3_out <- logistic_shortbred(prog, SVs_humann3, basic, adjust_cov, index_name)
  write.csv(SVs_humann3_out, file = out_humann3, row.names = FALSE)
}


basic_CRC<-subset(all_basic,is.na(CRC)==F,drop=T)
a<-table(basic_CRC$study_name,basic_CRC$CRC)
basic_CRC<-subset(all_basic,is.na(CRC)==F & study_name!="FengQ_2015",drop=T)
disease_humann3_asso(basic_CRC,"CRC",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/CRC2.csv")

basic_adenoma<-subset(all_basic,is.na(adenoma)==F,drop=T)
table(basic_adenoma$study_name,basic_adenoma$adenoma)
basic_adenoma<-subset(all_basic,is.na(adenoma)==F & study_name!="FengQ_2015",drop=T)
disease_humann3_asso(basic_adenoma,"adenoma",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/adenoma.csv")

basic_T2D<-subset(all_basic,is.na(T2D)==F,drop=T)
a<-table(basic_T2D$study_name,basic_T2D$T2D)
basic_T2D<-subset(all_basic,is.na(T2D)==F & study_name!="HMP_2019_t2d" & study_name!="LiJ_2014",drop=T)
disease_humann3_asso(basic_T2D,"T2D",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/T2D.csv")

basic_ACVD<-subset(all_basic,is.na(ACVD)==F,drop=T)
table(basic_ACVD$study_name,basic_ACVD$ACVD)
disease_humann3_asso(basic_ACVD,"ACVD",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/ACVD.csv")

basic_CAD<-subset(all_basic,is.na(CAD)==F,drop=T)
table(basic_CAD$study_name,basic_CAD$CAD)
disease_humann3_asso(basic_CAD,"CAD",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/CAD.csv")

basic_cirrhosis<-subset(all_basic,is.na(cirrhosis)==F,drop=T)
table(basic_cirrhosis$study_name,basic_cirrhosis$cirrhosis)
disease_humann3_asso(basic_cirrhosis,"cirrhosis",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/cirrhosis.csv")

basic_MetS<-subset(all_basic,is.na(MetS)==F,drop=T)
table(basic_MetS$study_name,basic_MetS$MetS)
disease_humann3_asso(basic_MetS,"MetS",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/MetS.csv")

basic_PD<-subset(all_basic,is.na(PD)==F,drop=T)
table(basic_PD$study_name,basic_PD$PD)
disease_humann3_asso(basic_PD,"PD",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/PD.csv")

basic_IBD<-subset(all_basic,is.na(IBD)==F,drop=T)
table(basic_IBD$study_name,basic_IBD$IBD)
basic_IBD<-subset(all_basic,is.na(IBD)==F & study_name!="LiJ_2014",drop=T)
disease_humann3_asso(basic_IBD,"IBD",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/IBD.csv")


############################################################################################
#######  index 
#######  2.3 linear regression model 1
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

common_rows <- intersect(rownames(age_basic), rownames(SVs_humann3_t))
#index_SVs_shortbred <- SVs_humann3_t[rownames(basic), ]
age_SVs_humann3 <- SVs_humann3_t[common_rows, ]
SVs_humann3 <- linear_shortbred(age_prog, age_SVs_humann3, age_basic,c("gender_code","read_count"), "age_norm")
write.csv(SVs_humann3, file = "10.Microbial_humann3/studies/age.csv", row.names = FALSE)


############# BMI all
BMI_basic<-subset(all_basic,is.na(BMI)==F,header=T)
table(BMI_basic$study_name)

BMI_basic<-subset(BMI_basic,study_name!="PiccinnoG_2025_cohort6",drop=T)

BMI_prog<-data.frame(BMI_basic[,"BMI_norm"])
row.names(BMI_prog)<-rownames(BMI_basic)
colnames(BMI_prog)<-"Y"

common_rows <- intersect(rownames(BMI_basic), rownames(SVs_humann3_t))
#SVs_humann3 <- SVs_humann3_t[rownames(BMI_basic), ]
BMI_SVs_humann3 <- SVs_humann3_t[common_rows, ]
SVs_humann3 <- linear_shortbred(BMI_prog, BMI_SVs_humann3, BMI_basic,c("gender_code","read_count"), "BMI_norm")
write.csv(SVs_humann3, file = "10.Microbial_humann3/studies/BMI.csv", row.names = FALSE)


############################################################################################################
## associations between index and genes

index_humann3_asso <- function(basic, index_name, adjust_cov, out_SVs_humann3_adjust) {
  # basic <- basic_hba1c
  # index_name <- "hba1c_norm"
  # adjust_cov <- c("age", "gender_code", "BMI", "read_count")

  prog <- data.frame(basic[, index_name])
  row.names(prog) <- rownames(basic)
  colnames(prog) <- "Y"

  common_rows <- intersect(rownames(basic), rownames(SVs_humann3_t))
  #index_SVs_humann3 <- SVs_humann3_t[rownames(basic), ]
  index_SVs_humann3 <- SVs_humann3_t[common_rows, ]
  SVs_humann3 <- linear_shortbred(prog, index_SVs_humann3, basic, adjust_cov, index_name)
  write.csv(SVs_humann3, file = out_SVs_humann3_adjust, row.names = FALSE)
}

### index
basic_hba1c<-subset(all_basic,is.na(hba1c)==F,drop=T)
table(basic_hba1c$study_name)
index_humann3_asso(basic_hba1c,"hba1c_norm",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/hba1c.csv")

basic_triglycerides<-subset(all_basic,is.na(triglycerides)==F,drop=T)
table(basic_triglycerides$study_name)
index_humann3_asso(basic_triglycerides,"triglycerides_norm",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/triglycerides.csv")

basic_fasting_insulin<-subset(all_basic,is.na(fasting_insulin)==F,drop=T)
table(basic_fasting_insulin$study_name)
index_humann3_asso(basic_fasting_insulin,"fasting_insulin_norm",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/fasting_insulin.csv")

basic_hdl<-subset(all_basic,is.na(hdl)==F,drop=T)
table(basic_hdl$study_name)
index_humann3_asso(basic_hdl,"hdl_norm",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/hdl.csv")

basic_ldl<-subset(all_basic,is.na(ldl)==F,drop=T)
table(basic_ldl$study_name)
index_humann3_asso(basic_ldl,"ldl_norm",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/ldl.csv")

basic_cholesterol<-subset(all_basic,is.na(cholesterol)==F,drop=T)
table(basic_cholesterol$study_name)
index_humann3_asso(basic_cholesterol,"cholesterol_norm",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/cholesterol.csv")

basic_ast<-subset(all_basic,is.na(ast)==F,drop=T)
index_humann3_asso(basic_ast,"ast_norm",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/ast.csv")

basic_mean_press<-subset(all_basic,is.na(mean_press)==F,drop=T)
index_humann3_asso(basic_mean_press,"mean_press_norm",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/mean_press.csv")

basic_alt<-subset(all_basic,is.na(alt)==F,drop=T)
index_humann3_asso(basic_alt,"alt_norm",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/alt.csv")

basic_total_bilirubin<-subset(all_basic,is.na(total_bilirubin)==F,drop=T)
index_humann3_asso(basic_total_bilirubin,"total_bilirubin_norm",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/total_bilirubin.csv")

basic_creatinine<-subset(all_basic,is.na(creatinine)==F,drop=T)
index_humann3_asso(basic_creatinine,"creatinine_norm",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/creatinine.csv")

basic_fasting_glucose<-subset(all_basic,is.na(fasting_glucose)==F,drop=T)
index_humann3_asso(basic_fasting_glucose,"fasting_glucose_norm",c("age","gender_code","BMI","read_count"),
"10.Microbial_humann3/studies/fasting_glucose.csv")





######################################################################################################
} # End optional cohort association analysis
if (run_meta_analysis) {
### Meta-analysis without abundance adjustment
library(data.table)
binary_multi<-c("CRC","IBD","T2D","adenoma","PD")
continuous<-c("age","BMI","hba1c","fasting_glucose","fasting_insulin",
              "triglycerides","hdl","ldl","cholesterol","mean_press",
              "creatinine","ast","alt","total_bilirubin")
binary_single<-c("ACVD","CAD","cirrhosis","MetS")
outcomes<-c(binary_multi,continuous,binary_single)

## 0.05/410000/23
study_wide_cutoff<-5.30e-9
analysis_wide_cutoff<-0.05
I2_cutoff<-75

reuse_existing_meta<-identical(tolower(Sys.getenv("SV_REUSE_META","false")),"true")
for(folder in c("meta","multi","anno","study_wide_sig","analysis_wide_sig"))
  dir.create(paste0("10.Microbial_humann3/",folder),recursive=T,showWarnings=F)

### Reuse 17 workers; reserve five of the 22 logical threads.
meta_workers<-17L
meta_pending<-c(binary_multi,continuous)
meta_pending<-meta_pending[file.exists(paste0("10.Microbial_humann3/studies/",meta_pending,".csv"))]
if(reuse_existing_meta)
  meta_pending<-meta_pending[!file.exists(paste0("10.Microbial_humann3/meta/",meta_pending,".csv"))]
meta_cl<-if(length(meta_pending)) parallel::makePSOCKcluster(meta_workers) else NULL
tryCatch({
  if(!is.null(meta_cl)) parallel::clusterEvalQ(meta_cl,{
    Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,VECLIB_MAXIMUM_THREADS=1)
    library(metafor)
    NULL
  })
  ### Binary outcomes: multiple cohorts
  for(outcome in binary_multi){
    study_file<-paste0("10.Microbial_humann3/studies/",outcome,".csv")
    meta_file<-paste0("10.Microbial_humann3/meta/",outcome,".csv")
    if(!file.exists(study_file) || reuse_existing_meta && file.exists(meta_file)) next
    study_data<-as.data.frame(fread(study_file,check.names=F))
    stopifnot("X" %in% names(study_data),!anyNA(study_data$X),!anyDuplicated(study_data$X),
              (ncol(study_data)-1)%%6==0)
    combine_meta_p_binary(study_data,6,meta_file,cluster=meta_cl)
    rm(study_data);gc()
  }

  ### Continuous outcomes
  for(outcome in continuous){
    study_file<-paste0("10.Microbial_humann3/studies/",outcome,".csv")
    meta_file<-paste0("10.Microbial_humann3/meta/",outcome,".csv")
    if(!file.exists(study_file) || reuse_existing_meta && file.exists(meta_file)) next
    study_data<-as.data.frame(fread(study_file,check.names=F))
    stopifnot("X" %in% names(study_data),!anyNA(study_data$X),!anyDuplicated(study_data$X),
              (ncol(study_data)-1)%%6==0)
    combine_meta_p_lr(study_data,6,meta_file,cluster=meta_cl)
    rm(study_data);gc()
  }
},finally={
  foreach::registerDoSEQ()
  if(!is.null(meta_cl)) parallel::stopCluster(meta_cl)
})

### Binary outcomes: one cohort
for(outcome in binary_single){
  study_file<-paste0("10.Microbial_humann3/studies/",outcome,".csv")
  meta_file<-paste0("10.Microbial_humann3/meta/",outcome,".csv")
  if(!file.exists(study_file) || reuse_existing_meta && file.exists(meta_file)) next
  study_data<-as.data.frame(fread(study_file,check.names=F))
  stopifnot("X" %in% names(study_data),!anyNA(study_data$X),!anyDuplicated(study_data$X),ncol(study_data)==7)
  valid<-is.finite(study_data[,2]) & study_data[,2]>0 & is.finite(study_data[,3]) & study_data[,3]>0
  study_data$OR_meta<-ifelse(valid,study_data[,2],NA_real_)
  study_data$se_meta<-ifelse(valid,study_data[,3],NA_real_)
  study_data$pvalue_meta<-ifelse(valid,study_data[,4],NA_real_)
  study_data$cilb<-ifelse(valid,study_data[,5],NA_real_)
  study_data$ciub<-ifelse(valid,study_data[,6],NA_real_)
  study_data$I2<-rep(NA_real_,nrow(study_data))
  study_data$het_p<-ifelse(valid,1,NA_real_)
  study_data$p_count<-ifelse(valid,as.integer(study_data$pvalue_meta<=0.05),NA_integer_)
  write.csv(study_data,meta_file,row.names=F)
  rm(study_data);gc()
}

######################################################################################################
### Multiple testing correction and significant results
meta_summary<-data.frame(outcome=character(),status=character(),features=integer(),
                         testable=integer(),study_wide_n=integer(),analysis_wide_n=integer())
missing_inputs<-data.frame(outcome=character(),study_exists=logical(),meta_exists=logical(),reason=character())
for(outcome in outcomes){
  study_file<-paste0("10.Microbial_humann3/studies/",outcome,".csv")
  meta_file<-paste0("10.Microbial_humann3/meta/",outcome,".csv")
  if(!file.exists(study_file) || !file.exists(meta_file)){
    missing_inputs<-rbind(missing_inputs,data.frame(outcome=outcome,
      study_exists=file.exists(study_file),meta_exists=file.exists(meta_file),
      reason="Missing input; no fitting or substitute output performed"))
    warning("Missing study or meta input: ",outcome)
    next
  }
  meta_pvalue<-as.data.frame(fread(meta_file,check.names=F))
  study_header<-names(fread(study_file,nrows=0,check.names=F))
  study_ids<-fread(study_file,select="UniRef90_id")$UniRef90_id
  stopifnot(all(c("X","pvalue_meta","I2") %in% names(meta_pvalue)),
            !anyNA(meta_pvalue$X),!anyDuplicated(meta_pvalue$X),!anyDuplicated(study_ids),
            identical(as.character(meta_pvalue$X),as.character(study_ids)),
            identical(names(meta_pvalue)[seq_along(study_header)],study_header),(length(study_header)-1)%%6==0)
  effect_cols<-grep("\\.(OR|beta)$",study_header,value=T)
  stopifnot(length(effect_cols)==(length(study_header)-1)/6)
  meta_pvalue$study_count<-rep(0L,nrow(meta_pvalue))
  for(effect_col in effect_cols){
    se_col<-sub("\\.(OR|beta)$",".SE",effect_col)
    valid<-is.finite(meta_pvalue[[effect_col]]) & is.finite(meta_pvalue[[se_col]]) & meta_pvalue[[se_col]]>0
    if(outcome %in% c(binary_multi,binary_single)) valid<-valid & meta_pvalue[[effect_col]]>0
    meta_pvalue$study_count<-meta_pvalue$study_count+as.integer(valid)
  }
  stopifnot(all(is.na(meta_pvalue$pvalue_meta) | is.finite(meta_pvalue$pvalue_meta) &
                meta_pvalue$pvalue_meta>=0 & meta_pvalue$pvalue_meta<=1))
  tested<-which(!is.na(meta_pvalue$pvalue_meta))
  meta_pvalue$bonferroni<-rep(NA_real_,nrow(meta_pvalue))
  meta_pvalue$fdr<-rep(NA_real_,nrow(meta_pvalue))
  meta_pvalue$bonferroni[tested]<-p.adjust(meta_pvalue$pvalue_meta[tested],"bonferroni")
  meta_pvalue$fdr[tested]<-p.adjust(meta_pvalue$pvalue_meta[tested],"BH")
  # Exactly one contributing cohort remains eligible when I2 is unavailable.
  heterogeneity_ok<-(!is.na(meta_pvalue$I2) & meta_pvalue$I2<I2_cutoff) | meta_pvalue$study_count==1
  study_index<-which(meta_pvalue$pvalue_meta<study_wide_cutoff & heterogeneity_ok)
  analysis_index<-which(meta_pvalue$bonferroni<analysis_wide_cutoff & heterogeneity_ok)
  fwrite(meta_pvalue,paste0("10.Microbial_humann3/multi/",outcome,".csv"),na="NA",quote="auto")
  fwrite(meta_pvalue[study_index,,drop=F],paste0("10.Microbial_humann3/study_wide_sig/",outcome,".csv"),na="NA",quote="auto")
  fwrite(meta_pvalue[analysis_index,,drop=F],paste0("10.Microbial_humann3/analysis_wide_sig/",outcome,".csv"),na="NA",quote="auto")
  meta_summary<-rbind(meta_summary,data.frame(outcome=outcome,status="validated_meta_postprocessed",
    features=nrow(meta_pvalue),testable=length(tested),study_wide_n=length(study_index),analysis_wide_n=length(analysis_index)))
  rm(meta_pvalue,study_ids);gc()
}
write.csv(meta_summary,"10.Microbial_humann3/meta_summary.csv",row.names=F)
write.csv(missing_inputs,"10.Microbial_humann3/missing_meta_inputs.csv",row.names=F)
print(meta_summary)




} # End optional meta-analysis and full-outcome correction
## Prepare FS9 inputs from completed meta-analysis; do not refit associations.
fs9_annotation <- fread("00.rawData/humann3/match_file/SGV_Finder_genes_vs_uniref90_annotation.tsv")
fs9_annotation <- fs9_annotation[, lapply(.SD, function(x) paste(unique(x[!is.na(x) & nzchar(x)]), collapse=";")),
  by=UniRef90_id, .SDcols=c("gene","ko","name","SV_location","dest")]
dir.create("10.Microbial_humann3/anno",showWarnings=FALSE)
dir.create("10.Microbial_humann3/multi",showWarnings=FALSE)
for (trait in c("BMI","age")) {
  completed <- fread(paste0("10.Microbial_humann3/meta/",trait,".csv"))
  if (!"UniRef90_id" %in% names(completed) && "X" %in% names(completed)) setnames(completed,"X","UniRef90_id")
  stopifnot("UniRef90_id" %in% names(completed), !anyDuplicated(completed$UniRef90_id))
  stopifnot(all(is.na(completed$pvalue_meta) | (is.finite(completed$pvalue_meta) & completed$pvalue_meta >= 0 & completed$pvalue_meta <= 1)))
  completed$bonferroni <- NA_real_; completed$fdr <- NA_real_
  tested <- which(!is.na(completed$pvalue_meta))
  completed$bonferroni[tested] <- p.adjust(completed$pvalue_meta[tested], "bonferroni")
  completed$fdr[tested] <- p.adjust(completed$pvalue_meta[tested], "BH")
  fwrite(completed,paste0("10.Microbial_humann3/multi/",trait,".csv"),na="NA")
  completed <- merge(completed,fs9_annotation,by="UniRef90_id",all.x=TRUE,sort=FALSE)
  fwrite(completed,paste0("10.Microbial_humann3/anno/",trait,".csv"),na="NA")
}
## 3.1 Reference-gene location, summarized once per UniRef90 family
SV_gene_location<-fread("00.rawData/humann3/match_file/SGV_Finder_genes_vs_uniref90_annotation.tsv",
  select=c("UniRef90_id","SV_location"),na.strings=c("","NA"))
SV_gene_location$SV_location<-trimws(SV_gene_location$SV_location)
stopifnot(!anyNA(SV_gene_location$UniRef90_id),!anyNA(SV_gene_location$SV_location),
          all(SV_gene_location$SV_location %in% c("Yes","No")))
# A HUMAnN3 family can map to both locations. Any Yes is displayed as a triangle;
# this denotes an SV-region reference-gene mapping, not SV-specific abundance.
SV_family_location<-SV_gene_location[,.(
  SV_location=paste(sort(unique(SV_location)),collapse=";"),
  SV_region=if(any(SV_location=="Yes")) "Includes SV genes" else "Non-SV genes only",
  SV_location_mixed=uniqueN(SV_location)>1L),by=UniRef90_id]
stopifnot(!anyDuplicated(SV_family_location$UniRef90_id))
SV_region_levels<-c("Non-SV genes only","Includes SV genes")
SV_region_shapes<-c("Non-SV genes only"=16,"Includes SV genes"=17)
SV_region_labels<-c("Non-SV genes only"="genes not covered by SV",
                    "Includes SV genes"="genes covered by SV")
dir.create("10.Microbial_humann3/pic_files",showWarnings=FALSE)
fwrite(SV_family_location,"10.Microbial_humann3/pic_files/UniRef90_SV_location.csv",na="NA")

## 4 Shared FS9 figure settings and microbial KEGG BRITE L1 classes
KEGG_L1<-unique(fread("00.rawData/KEGG/KEGG_classification.csv",
  select=c("KO","L1"),data.table=FALSE))
KEGG_L1<-KEGG_L1[grepl("^091(00|20|30|40) ",KEGG_L1$L1),]
KEGG_L1$L1<-sub("^[0-9]+ ","",KEGG_L1$L1)
KEGG_L1<-unique(KEGG_L1)
category_order<-c("Metabolism","Genetic Information Processing",
  "Environmental Information Processing","Cellular Processes","Multiple KEGG categories",
  "KO annotated, no functional class","Unknown / unannotated")
category_colours<-setNames(c("#207F4C","#ECD452","#2983BB","#C72F1E",
                            "#EA8958","#C4C4C4","#696969"),category_order)
category_labels<-c("Metabolism","Genetic information\nprocessing",
  "Environmental information\nprocessing","Cellular processes","Multiple KEGG categories",
  "KO annotated, no\nfunctional class","Unknown / unannotated")
direction_colours<-c("Not significant"="#696969","Negative"="#2983BB","Positive"="#BC3523")
dir.create("pics",showWarnings=FALSE)
dir.create("10.Microbial_humann3/pic_files",showWarnings=FALSE)

## 5 BMI and age, in the same order and with the same classification rules
for(trait in c("BMI","age")) {
  trait_label<-if(trait=="BMI") "BMI" else "Age"
  associations<-as.data.frame(fread(paste0("10.Microbial_humann3/anno/",trait,".csv"),
    select=c("UniRef90_id","gene","ko","name","effect_meta","pvalue_meta","bonferroni","fdr","I2")))
  stopifnot(!anyNA(associations$UniRef90_id),!anyDuplicated(associations$UniRef90_id),
    all(is.na(associations$bonferroni) | associations$bonferroni>=0 & associations$bonferroni<=1))
  location_index<-match(associations$UniRef90_id,SV_family_location$UniRef90_id)
  stopifnot(!anyNA(location_index))
  associations$SV_location<-SV_family_location$SV_location[location_index]
  associations$SV_region<-factor(SV_family_location$SV_region[location_index],levels=SV_region_levels)
  associations$SV_location_mixed<-SV_family_location$SV_location_mixed[location_index]
  tested_ids<-fread(paste0("10.Microbial_humann3/multi/",trait,".csv"),select="UniRef90_id")$UniRef90_id
  stopifnot(setequal(associations$UniRef90_id,tested_ids))
  valid<-is.finite(associations$effect_meta) & is.finite(associations$bonferroni)
  associations$significant<-valid & associations$bonferroni<0.05
  stopifnot(any(valid),!any(associations$significant & associations$effect_meta==0))
  associations$Direction<-"Not significant"
  associations$Direction[associations$significant & associations$effect_meta<0]<-"Negative"
  associations$Direction[associations$significant & associations$effect_meta>0]<-"Positive"
  associations$Direction<-factor(associations$Direction,levels=names(direction_colours))
  associations$neg_log10_bonf<--log10(pmax(associations$bonferroni,.Machine$double.xmin))
  significant<-associations[associations$significant,]
  write.csv(significant,paste0("10.Microbial_humann3/pic_files/FS9_",trait,"_significant_families.csv"),row.names=FALSE)

  ### 5.1 FS9A: the existing volcano plot, with SV-region shape encoding
  volcano_data<-associations[valid & abs(associations$effect_meta)<=0.01,]
  volcano_data<-volcano_data[order(volcano_data$significant,volcano_data$SV_region),]
  volcano<-ggplot(volcano_data,aes(effect_meta,neg_log10_bonf))+
    geom_hline(yintercept=-log10(.05),linetype="dashed",colour="grey50",linewidth=.3)+
    geom_vline(xintercept=0,linetype="dotted",colour="grey50",linewidth=.3)+
    geom_point(aes(colour=Direction,shape=SV_region),size=1.1,alpha=.6)+
    scale_colour_manual(values=direction_colours,drop=FALSE)+
    scale_shape_manual(values=SV_region_shapes,labels=SV_region_labels,drop=FALSE,name=NULL)+
    labs(title=paste0(trait_label,": HUMAnN3 gene-family associations"),
      subtitle=sprintf("%s significant families: %s positive, %s negative",
        format(nrow(significant),big.mark=","),format(sum(significant$effect_meta>0),big.mark=","),
        format(sum(significant$effect_meta<0),big.mark=",")),
      x=expression("Meta-analysis "*beta),y=expression(-log[10](P[Bonferroni])),colour=NULL)+
    guides(colour=guide_legend(order=1,nrow=1,override.aes=list(size=3,alpha=1,shape=16)),
      shape=guide_legend(order=2,nrow=1,override.aes=list(size=3,alpha=1,colour="black")))+
    theme_classic(base_size=13,base_family="Arial")+
    theme(legend.position="bottom",legend.box="vertical",plot.title=element_text(face="bold"),
          plot.margin=margin(10,12,10,10))
  # Keep the previous labels and add metabolic functions among the top 30 associations.
  if(trait=="BMI") {
    volcano<-volcano+theme(text=element_text(size=16),
      axis.text=element_text(size=13.4),axis.title=element_text(size=16),
      plot.title=element_text(size=18.6),plot.subtitle=element_text(size=16),
      legend.text=element_text(size=13.4),legend.title=element_text(size=16))
    BMI_label_names<-c("Fructose-1,6-bisphosphate aldolase class II",
      "Pyruvate, phosphate dikinase","Ion-translocating oxidoreductase complex subunit A",
      "Anaerobic ribonucleoside-triphosphate reductase","Glucose-1-phosphate thymidylyltransferase",
      "Glutamate dehydrogenase","Phosphoserine aminotransferase",
      "4-deoxy-L-threo-5-hexosulose-uronate ketol-isomerase",
      "Uracil phosphoribosyltransferase","Malic enzyme, NAD binding domain protein",
      "Hypoxanthine phosphoribosyltransferase","Carbamoyl-phosphate synthase (glutamine-hydrolyzing)")
    BMI_top30<-head(significant[order(significant$bonferroni,significant$pvalue_meta,
                                    significant$UniRef90_id),],30)
    BMI_top30$significance_rank<-seq_len(nrow(BMI_top30))
    # Ignore punctuation when matching; use the requested names on the figure.
    label_index<-match(gsub("[^a-z0-9]","",tolower(BMI_top30$name)),
                       gsub("[^a-z0-9]","",tolower(BMI_label_names)))
    BMI_top30$labelled<-!is.na(label_index)
    write.csv(BMI_top30,"10.Microbial_humann3/pic_files/FS9_BMI_top30_annotation_review.csv",row.names=FALSE)
    label_data<-BMI_top30[!is.na(label_index),]
    label_data$label_name<-BMI_label_names[label_index[!is.na(label_index)]]
    stopifnot(nrow(label_data)==12L,setequal(label_data$label_name,BMI_label_names),
              all(label_data$UniRef90_id %in% volcano_data$UniRef90_id))
    label_data$label<-vapply(label_data$label_name,
      function(x) paste(strwrap(x,width=36),collapse="\n"),character(1))
    volcano<-volcano+
      geom_point(data=label_data,aes(colour=Direction,shape=SV_region),size=1.8,show.legend=FALSE)+
      scale_y_continuous(expand=c(0,0))+
      coord_cartesian(xlim=c(-.011,.011),
        ylim=c(-.3,max(volcano_data$neg_log10_bonf)*1.75),expand=FALSE)
    label_data$label_x<-label_data$label_y<-NA_real_
    for(side in c("Negative","Positive")) {
      label_rows<-which(label_data$Direction==side)
      label_rows<-label_rows[order(-label_data$neg_log10_bonf[label_rows])]
      label_data$label_x[label_rows]<-if(side=="Negative") -.0036 else .0046
      label_heights<-1.1*lengths(strsplit(label_data$label[label_rows],"\n",fixed=TRUE))+.8
      label_centres<-cumsum(label_heights)-label_heights/2
      top_y<-if(side=="Negative") max(label_data$neg_log10_bonf[label_rows]) else max(volcano_data$neg_log10_bonf)*1.67
      label_data$label_y[label_rows]<-top_y-(label_centres-label_centres[1])/
        diff(range(label_centres))*(top_y-9.5)
      volcano<-volcano+
        geom_segment(data=label_data[label_rows,],
          aes(x=effect_meta,y=neg_log10_bonf,xend=label_x,yend=label_y,colour=Direction),
          linewidth=.3,alpha=.8,inherit.aes=FALSE,show.legend=FALSE)+
        geom_label(data=label_data[label_rows,],aes(x=label_x,y=label_y,label=label,colour=Direction),
          hjust=if(side=="Negative") 1 else 0,size=2.6+3/ggplot2::.pt,fill="white",label.size=.25,
          label.padding=grid::unit(.15,"lines"),lineheight=1.1,inherit.aes=FALSE,show.legend=FALSE)
    }
    write.csv(label_data,"10.Microbial_humann3/pic_files/FS9_BMI_volcano_gene_labels.csv",row.names=FALSE)
  }
  ragg::agg_tiff(paste0("pics/FS9A_",trait,"_HUMAnN3_volcano.tiff"),
    width=if(trait=="BMI") 11.25 else 7.5,height=if(trait=="BMI") 8 else 5.5,units="in",res=300,compression="lzw")
  print(volcano);dev.off()
  saveRDS(volcano,paste0("10.Microbial_humann3/pic_files/FS9_",trait,"_volcano.plot.rds"))
  write.csv(as.data.frame(with(associations,table(SV_region,SV_location_mixed,Direction))),
    paste0("10.Microbial_humann3/pic_files/FS9_",trait,"_SV_location_summary.csv"),row.names=FALSE)

  ### 5.2 Assign each significant family to exactly one L1 category
  ko_text<-ifelse(is.na(significant$ko),"",significant$ko)
  ko_list<-regmatches(ko_text,gregexpr("K[0-9]{5}",ko_text))
  family_KO<-unique(data.frame(UniRef90_id=rep(significant$UniRef90_id,lengths(ko_list)),
    KO=as.character(unlist(ko_list,use.names=FALSE))))
  membership<-unique(merge(family_KO,KEGG_L1,by="KO")[,c("UniRef90_id","L1")])
  family_class<-significant[,c("UniRef90_id","Direction","SV_location","SV_region","SV_location_mixed")]
  family_class$Category<-"Unknown / unannotated"
  family_class$Category[family_class$UniRef90_id %in% family_KO$UniRef90_id]<-"KO annotated, no functional class"
  single<-membership[!duplicated(membership$UniRef90_id) & !duplicated(membership$UniRef90_id,fromLast=TRUE),]
  family_class$Category[match(single$UniRef90_id,family_class$UniRef90_id)]<-single$L1
  multiple<-names(which(table(membership$UniRef90_id)>1))
  family_class$Category[family_class$UniRef90_id %in% multiple]<-"Multiple KEGG categories"
  stopifnot(!anyDuplicated(family_class$UniRef90_id),all(family_class$Category %in% category_order))
  write.csv(membership,paste0("10.Microbial_humann3/pic_files/FS9_",trait,"_family_KEGG_L1_mapping.csv"),row.names=FALSE)
  write.csv(family_class,paste0("10.Microbial_humann3/pic_files/FS9_",trait,"_family_L1_classes.csv"),row.names=FALSE)

  ### 5.3 FS9B/C: positive and negative associations, ring left and legend right
  composition<-data.frame()
  for(direction in c("Positive","Negative")) {
    selected<-family_class[family_class$Direction==direction,]
    counts<-as.integer(table(factor(selected$Category,levels=category_order)))
    stopifnot(sum(counts)==nrow(selected))
    ring_data<-data.frame(Category=factor(category_order,levels=category_order),Families=counts)
    ring_data$Percent<-100*counts/max(1,sum(counts))
    ring_data$ymax<-cumsum(counts);ring_data$ymin<-ring_data$ymax-counts
    ring_data$ymid<-(ring_data$ymin+ring_data$ymax)/2
    ring_data$text_colour<-ifelse(category_order %in% c("Genetic Information Processing",
      "Multiple KEGG categories","KO annotated, no functional class"),"black","white")
    legend_labels<-paste0(category_labels,"\n",counts," (",sprintf("%.1f",ring_data$Percent),"%)")
    ring<-ggplot(ring_data)+
      geom_rect(aes(xmin=1.1,xmax=2.9,ymin=ymin,ymax=ymax,fill=Category),
        colour="white",linewidth=.35,key_glyph="point")+
      geom_text(data=ring_data[ring_data$Percent>=5,],
        aes(x=2.1,y=ymid,label=sprintf("%.1f%%",Percent),colour=text_colour),size=4.5)+
      annotate("text",x=0,y=0,label=paste0(sum(counts),"\nUniRef90\nfamilies"),
        colour=direction_colours[direction],fontface="bold",size=5,lineheight=1.3)+
      scale_x_continuous(limits=c(0,2.9),expand=c(0,0))+
      scale_y_continuous(limits=c(0,max(1,sum(counts))),expand=c(0,0))+
      scale_fill_manual(values=category_colours,breaks=category_order,labels=legend_labels,drop=FALSE)+
      scale_colour_identity()+coord_polar(theta="y",start=0)+
      guides(fill=guide_legend(override.aes=list(shape=22,size=5,colour=NA,stroke=0)))+
      labs(title=paste0(trait_label,": ",tolower(direction)," associations"),
        subtitle=paste0("Meta beta ",if(direction=="Positive") ">" else "<"," 0; Bonferroni P < 0.05"),
        fill=NULL)+
      theme_void(base_size=13,base_family="Arial")+
      theme(plot.title=element_text(size=18,face="bold",hjust=.5,margin=margin(b=5)),
        plot.subtitle=element_text(size=14,hjust=.5,margin=margin(b=16)),
        legend.position="right",legend.text=element_text(size=12,margin=margin(t=3,b=3)),
        legend.key.size=grid::unit(.42,"cm"),legend.spacing.y=grid::unit(.12,"cm"),
        plot.margin=margin(14,18,14,14))
    # Tighten BMI donut exports without changing the ring or legend contents.
    if(trait=="BMI") ring<-ring+theme(plot.margin=margin(8,6,8,6),
      legend.box.spacing=grid::unit(3,"pt"),legend.margin=margin(0,0,0,0))
    prefix<-paste0("FS9",if(direction=="Positive") "B" else "C","_",trait,
                   "_HUMAnN3_",if(direction=="Positive") "up" else "down","_KEGG")
    ragg::agg_tiff(file.path("pics",paste0(prefix,".tiff")),
      width=if(trait=="BMI") 8.8 else 10.8,height=7,units="in",res=300,compression="lzw")
    print(ring);dev.off()
    saveRDS(ring,file.path("10.Microbial_humann3/pic_files",paste0(prefix,".plot.rds")))
    composition<-rbind(composition,data.frame(Direction=direction,
      ring_data[,c("Category","Families","Percent")]))
  }
  write.csv(composition,paste0("10.Microbial_humann3/pic_files/FS9_",trait,"_KEGG_L1_composition.csv"),row.names=FALSE)
  cat("FS9",trait,":",nrow(associations),"families;",nrow(significant),"significant;",
      sum(significant$effect_meta>0),"positive;",sum(significant$effect_meta<0),"negative\n")
}
