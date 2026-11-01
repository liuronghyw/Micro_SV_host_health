### clinical information
### Liurong
### 2026-9-21
# Donor selection below is sequential; missing IDs are retained and chronology is audited.

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/00.rawData")


#source("functions.R")
##install.packages("ggmosaic")
##library("plyr")
##library("dplyr")

library(foreach)
library(doParallel)
library(dplyr)
library(magrittr)
library(car)
library("scales")

##install.packages("readr")
library("readr")

###################################################################################
clinical<-unique(read.csv("clinical/metadata_all_datasets.csv",header=T))
normalize_study_name <- function(x) sub("Fern[^_]*ndez_Calvet_2024", "Fern\u00e1ndez_Calvet_2024", as.character(x))
clinical$study_name <- normalize_study_name(clinical$study_name)
unique(clinical$study_name)

## select stool sample 
table(clinical$body_site)
clinical <- subset(clinical,body_site=="stool",drop=T)
unique(clinical$study_name)

##select samples age adult or senior
table(clinical$age_category)

clinical <- subset(clinical,as.numeric(age)>=18 | age_category=="adult" | age_category=="senior",drop=T)
unique(clinical$study_name)

##select no pregnant 
table(clinical$pregnant)
clinical<- subset(clinical,pregnant!="yes" | is.na(pregnant)==T,drop=T)
unique(clinical$study_name)

nr<-nrow(clinical)
  for(i in 1:nr){
    if(clinical$study_name[i]=="ThomasAM_2018a"){clinical$study_name[i]="ThomasAM_2018"}
    else if (clinical$study_name[i]=="ThomasAM_2018b"){clinical$study_name[i]="ThomasAM_2018"}
    }

SV_study<-read.csv("dataset_select/SV_study_names.csv",header=T)
SV_study$study_name <- normalize_study_name(SV_study$study_name)
SV_study<-data.frame(unique(SV_study$study_name))
colnames(SV_study)<-"study_name"

clinical_SV<-unique(merge(SV_study,clinical,by.x="study_name",by.y="study_name"))
test_study<-unique(clinical_SV$study_name)

#test<-subset(clinical,disease=="PD",drop=T)
#write.csv(clinical_SV,"clinical/clinical_SV.csv")
#write.csv(test_study,"clinical/test_study.csv")


##################################################
### determine the phenotype

disease<-table(clinical_SV$disease)
write.csv(disease,"dataset_select/disease.csv")

clin_merge<-clinical_SV
study_contry<-table(clin_merge$study_name,clin_merge$country)
table(clin_merge$country)

nrow(subset(clin_merge,is.na(age)==F,drop=T))
nrow(subset(clin_merge,is.na(BMI)==F,drop=T))
nrow(subset(clin_merge,is.na(ast)==F,drop=T))
nrow(subset(clin_merge,is.na(alt)==F,drop=T))
nrow(subset(clin_merge,is.na(triglycerides)==F,drop=T))
nrow(subset(clin_merge,is.na(hdl)==F,drop=T))
nrow(subset(clin_merge,is.na(ldl)==F,drop=T))
nrow(subset(clin_merge,is.na(hba1c)==F,drop=T))
nrow(subset(clin_merge,is.na(cholesterol)==F,drop=T))
nrow(subset(clin_merge,is.na(creatinine)==F,drop=T))
nrow(subset(clin_merge,is.na(total_bilirubin)==F,drop=T))
nrow(subset(clin_merge,is.na(prothrombin_time)==F,drop=T))
nrow(subset(clin_merge,is.na(dyastolic_p)==F,drop=T))
nrow(subset(clin_merge,is.na(systolic_p)==F,drop=T))
nrow(subset(clin_merge,is.na(fasting_insulin)==F,drop=T))
nrow(subset(clin_merge,is.na(fasting_glucose)==F,drop=T))


### sample size
nrow(subset(clin_merge,is.na(eGFR)==F,drop=T))
nrow(subset(clin_merge,is.na(leptin)==F,drop=T))
nrow(subset(clin_merge,is.na(adiponectin)==F,drop=T))
nrow(subset(clin_merge,is.na(calprotectin)==F,drop=T))
nrow(subset(clin_merge,is.na(direct_bilirubin)==F,drop=T))
nrow(subset(clin_merge,is.na(inr)==F,drop=T))
nrow(subset(clin_merge,is.na(bristol_score)==F,drop=T))
nrow(subset(clin_merge,is.na(glutamate_decarboxylase_2_antibody)==F,drop=T))
nrow(subset(clin_merge,is.na(ESR)==F,drop=T))
nrow(subset(clin_merge,is.na(BASDAI)==F,drop=T))
nrow(subset(clin_merge,is.na(BASFI)==F,drop=T))
nrow(subset(clin_merge,is.na(c_peptide)==F,drop=T))
nrow(subset(clin_merge,is.na(insulin)==F,drop=T))
nrow(subset(clin_merge,is.na(glp_1)==F,drop=T))
nrow(subset(clin_merge,is.na(cd163)==F,drop=T))
nrow(subset(clin_merge,is.na(il_1)==F,drop=T))
nrow(subset(clin_merge,is.na(fgf_19)==F,drop=T))
nrow(subset(clin_merge,is.na(glutamate_decarboxylase_2_antibody)==F,drop=T))
nrow(subset(clin_merge,is.na(autoantibody_positive)==F,drop=T))
nrow(subset(clin_merge,is.na(shigatoxin_2_elisa)==F,drop=T))
nrow(subset(clin_merge,is.na(globulin)==F,drop=T)) 
nrow(subset(clin_merge,is.na(urea_nitrogen)==F,drop=T))
nrow(subset(clin_merge,is.na(albumine)==F,drop=T))
nrow(subset(clin_merge,is.na(hscrp)==F,drop=T))



#################################################################################
###��select the samples (disease or index)

### diseaase
### IBD T2D CRC severe obesity MetS CAD
### ACVD adenoma cirrhosis IGT PD

### index 
## age  BMI fasting_insulin  hba1c  fasting_glucose  adiponectin
## triglycerides hdl  ldl  cholesterol hscrp
## dyastolic_p systolic_p
## albumine  globulin creatinine  alt  ast  total_bilirubin


disease_pmid<-table(subset(clin_merge,disease=="IBD" | disease=="T2D" | disease=="CRC" |
   disease=="severe obesity" | disease=="chronic CAD" | disease=="acute CAD" |
   disease=="ACVD" | disease=="chronic CAD;HF" | disease=="adenoma"|
   disease=="cirrhosis" | disease=="PD",drop=T)$study_name)

clin_select<-subset(clin_merge,
 study_name=="FengQ_2015" | study_name=="Fern\u00e1ndez_Calvet_2024" |study_name=="FranzosaEA_2018" | study_name=="HallAB_2017" | 
 study_name=="HMP_2019_ibdmdb"  | study_name=="JieZ_2017" | 
 study_name=="JoSY_2022" | study_name=="KarlssonFH_2013" | study_name=="LiJ_2014" | 
 study_name=="LiuNN_2022" |  study_name=="MaoLM_2021" | study_name=="MetaCardis_2020_a" |
 study_name=="PiccinnoG_2025_cohort1" |  study_name=="PiccinnoG_2025_cohort2" | study_name=="PiccinnoG_2025_cohort3" |
 study_name=="PiccinnoG_2025_cohort4" |  study_name=="PiccinnoG_2025_cohort5" | 
 study_name=="QinJ_2012" | study_name=="QinN_2014" |
 study_name=="ThomasAM_2018"  | study_name=="VogtmannE_2016" | study_name=="WirbelJ_2018" |  
 study_name=="YachidaS_2019" | study_name=="YangJ_2020" |  study_name=="YangY_2021" | 
  study_name=="YuJ_2015" | study_name=="ZellerG_2014" |
 is.na(age)==F | is.na(BMI)==F | is.na(fasting_insulin)==F | 
 is.na(hba1c)==F |  is.na(fasting_glucose)==F |  
 is.na(triglycerides)==F |  is.na(hdl)==F | is.na(ldl)==F |  is.na(cholesterol)==F | 
 is.na(dyastolic_p)==F |  is.na(systolic_p)==F | 
 is.na(creatinine)==F | is.na(ast)==F |   is.na(alt)==F |
 is.na(total_bilirubin)==F ,drop=T)

study_names<-data.frame(table(clin_select$study_name))
write.csv(study_names,"dataset_select/study_name_select.csv")

colnames(study_names)<-c("study_name","freq")


study_name_SRA<-read.csv("dataset_select/study_name_SRA_id_match.csv",header=T)

studies<-unique(merge(study_names,study_name_SRA,by.x="study_name",by.y="study_name",all.x=T))
write.csv(studies,"dataset_select/studies_SV_final.csv",row.names=F)

## all the studies were selected
## information for table S1

pmid<-unique(clin_select[,c("study_name","PMID")])
table_S1<-unique(merge(studies,pmid,by.x="study_name",by.y="study_name"))
write.csv(table_S1,"dataset_select/table_S1.csv",row.names=F)



##################################################################################################
###  get the matched id 

## AsnicarF_2021	PRJEB39223

PRJEB39223<-subset(clinical,study_name=="AsnicarF_2021",drop=T)
PRJEB39223$match_id<-PRJEB39223$NCBI_accession
nrow(PRJEB39223)

PRJEB39223_Count<-unique(read.table("read_number/PRJEB39223_read_number.txt",header=F))
colnames(PRJEB39223_Count)<-c("match_id","read_count")
PRJEB39223_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB39223_Count$match_id, fixed = TRUE)
PRJEB39223_whole<-merge(PRJEB39223,PRJEB39223_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB39223_whole)

#### CosteaPI_2017	PRJEB17632

PRJEB17632<-subset(clinical,study_name=="CosteaPI_2017",drop=T)
PRJEB17632$match_id<-PRJEB17632$sample_id
nrow(PRJEB17632)

PRJEB17632_Count<-unique(read.table("read_number/PRJEB17632_read_number.txt",header=F))
colnames(PRJEB17632_Count)<-c("match_id","read_count")
PRJEB17632_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB17632_Count$match_id, fixed = TRUE)
PRJEB17632_Count$match_id<-gsub("_kneaddata.fastq", "", PRJEB17632_Count$match_id, fixed = TRUE)

PRJEB17632_whole<-merge(PRJEB17632,PRJEB17632_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB17632_whole)


### FengQ_2015	PRJEB7774
PRJEB7774<-subset(clinical,study_name=="FengQ_2015",drop=T)
PRJEB7774$match_id<-PRJEB7774$sample_id
nrow(PRJEB7774)

PRJEB7774_Count<-unique(read.table("read_number/PRJEB7774_read_number.txt",header=F))
colnames(PRJEB7774_Count)<-c("match_id","read_count")
PRJEB7774_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB7774_Count$match_id, fixed = TRUE)
PRJEB7774_whole<-merge(PRJEB7774,PRJEB7774_Count,by.x="match_id",by.y="match_id")
nrow(PRJEB7774_whole)

### FranzosaEA_2018	PRJNA400072
PRJNA400072<-subset(clinical,study_name=="FranzosaEA_2018",drop=T)
PRJNA400072$match_id<-PRJNA400072$NCBI_accession
nrow(PRJNA400072)

PRJNA400072_Count<-unique(read.table("read_number/PRJNA400072_read_number.txt",header=F))
colnames(PRJNA400072_Count)<-c("match_id","read_count")
PRJNA400072_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA400072_Count$match_id, fixed = TRUE)
PRJNA400072_whole<-merge(PRJNA400072,PRJNA400072_Count,by.x="match_id",by.y="match_id")
nrow(PRJNA400072_whole)

### HallAB_2017	PRJNA385949

PRJNA385949<-subset(clinical,study_name=="HallAB_2017",drop=T)
PRJNA385949$match_id<-PRJNA385949$NCBI_accession
nrow(PRJNA385949)

PRJNA385949_Count<-unique(read.table("read_number/PRJNA385949_read_number.txt",header=F))
colnames(PRJNA385949_Count)<-c("match_id","read_count")
PRJNA385949_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA385949_Count$match_id, fixed = TRUE)
PRJNA385949_whole<-merge(PRJNA385949,PRJNA385949_Count,by.x="match_id",by.y="match_id")	
nrow(PRJNA385949_whole)


### HansenLBS_2018	PRJNA491335
PRJNA491335<-subset(clinical,study_name=="HansenLBS_2018",drop=T)
PRJNA491335$match_id<-PRJNA491335$NCBI_accession
nrow(PRJNA491335)

PRJNA491335_Count<-unique(read.table("read_number/PRJNA491335_read_number.txt",header=F))
colnames(PRJNA491335_Count)<-c("match_id","read_count")
PRJNA491335_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA491335_Count$match_id, fixed = TRUE)
PRJNA491335_whole<-merge(PRJNA491335,PRJNA491335_Count,by.x="match_id",by.y="match_id")
nrow(PRJNA491335_whole)


### HMP_2012

HMP_2012_clin<-subset(clinical,study_name=="HMP_2012",drop=T)
HMP_2012_clin$match_id<-HMP_2012_clin$sample_id
nrow(HMP_2012_clin)

HMP_2012_Count<-unique(read.table("read_number/HMP_2012_read_number.txt",header=F))
colnames(HMP_2012_Count)<-c("match_id","read_count")
HMP_2012_Count$match_id<-gsub(".denovo_duplicates_marked.trimmed.1_kneaddata_paired_1.fastq", "", HMP_2012_Count$match_id, fixed = TRUE)
HMP_2012_whole<-merge(HMP_2012_clin,HMP_2012_Count,by.x="match_id",by.y="match_id")
nrow(HMP_2012_whole)



### HMP_2019_ibdmdb	
PRJNA389280_id<-read.csv("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/2_pipe_qc_snp_sv_abun/2_submit_index/samples/unique/HMP_2019_ibdmdb_PRJNA389280.csv",header=T)
PRJNA389280_id<-PRJNA389280_id[,c(1,4)]
colnames(PRJNA389280_id)[2]<-"match_id"

PRJNA389280_clin<-subset(clinical,study_name=="HMP_2019_ibdmdb",drop=T)
PRJNA389280<-merge(PRJNA389280_clin,PRJNA389280_id,by.x="sample_id",by.y="sample_id")
nrow(PRJNA389280)

PRJNA389280_Count<-unique(read.table("read_number/PRJNA389280_read_number.txt",header=F))
colnames(PRJNA389280_Count)<-c("match_id","read_count")
PRJNA389280_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA389280_Count$match_id, fixed = TRUE)
PRJNA389280_whole<-merge(PRJNA389280,PRJNA389280_Count,by.x="match_id",by.y="match_id")
nrow(PRJNA389280_whole)


### HMP_2019_ibdmdb	PRJNA398089 
PRJNA398089_id<-read.csv("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/2_pipe_qc_snp_sv_abun/2_submit_index/samples/unique/HMP_2019_ibdmdb_PRJNA398089.csv",header=T)

PRJNA398089_id<-PRJNA398089_id[,c(1,4)]
colnames(PRJNA398089_id)[2]<-"match_id"

PRJNA398089_clin<-subset(clinical,study_name=="HMP_2019_ibdmdb",drop=T)
PRJNA398089<-merge(PRJNA398089_clin,PRJNA398089_id,by.x="sample_id",by.y="sample_id")
nrow(PRJNA398089)

PRJNA398089_Count<-unique(read.table("read_number/PRJNA398089_read_number.txt",header=F))
colnames(PRJNA398089_Count)<-c("match_id","read_count")
PRJNA398089_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA398089_Count$match_id, fixed = TRUE)
PRJNA398089_whole<-merge(PRJNA398089,PRJNA398089_Count,by.x="match_id",by.y="match_id")
nrow(PRJNA398089_whole)

## dataset have replicated sample id 
PRJNA398089_whole<-subset(PRJNA398089_whole,sample_id!="CSM67U9H_P" & sample_id!="MSM9VZMA_TR" 
                           & sample_id !="CSM67UEW_TR" & sample_id!="CSM79HLA_TR" & sample_id!="CSM67UEW_P",drop=T)


### JieZ_2017	PRJEB21528
PRJEB21528<-subset(clinical,study_name=="JieZ_2017",drop=T)
PRJEB21528$match_id<-PRJEB21528$NCBI_accession
nrow(PRJEB21528)

PRJEB21528_Count<-unique(read.table("read_number/PRJEB21528_read_number.txt",header=F))
colnames(PRJEB21528_Count)<-c("match_id","read_count")
PRJEB21528_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB21528_Count$match_id, fixed = TRUE)
PRJEB21528_whole<-merge(PRJEB21528,PRJEB21528_Count,by.x="match_id",by.y="match_id")
nrow(PRJEB21528_whole)

##KarlssonFH_2013	PRJEB1786
PRJEB1786<-subset(clinical,study_name=="KarlssonFH_2013",drop=T)
PRJEB1786$match_id<-PRJEB1786$NCBI_accession
nrow(PRJEB1786)

PRJEB1786_Count<-unique(read.table("read_number/PRJEB1786_read_number.txt",header=F))
colnames(PRJEB1786_Count)<-c("match_id","read_count")
PRJEB1786_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB1786_Count$match_id, fixed = TRUE)
PRJEB1786_whole<-merge(PRJEB1786,PRJEB1786_Count,by.x="match_id",by.y="match_id")
nrow(PRJEB1786_whole)


### KeohaneDM_2020	PRJEB36820

PRJEB36820<-subset(clinical,study_name=="KeohaneDM_2020",drop=T)
PRJEB36820$match_id<-PRJEB36820$NCBI_accession
nrow(PRJEB36820)

PRJEB36820_Count<-unique(read.table("read_number/PRJEB36820_read_number.txt",header=F))
colnames(PRJEB36820_Count)<-c("match_id","read_count")
PRJEB36820_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB36820_Count$match_id, fixed = TRUE)
PRJEB36820_whole<-merge(PRJEB36820,PRJEB36820_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB36820_whole)


### LiJ_2014	PRJEB5224 
PRJEB5224<-subset(clinical,study_name=="LiJ_2014",drop=T)
PRJEB5224$match_id<-PRJEB5224$sample_id
nrow(PRJEB5224)

PRJEB5224_Count<-unique(read.table("read_number/PRJEB5224_read_number.txt",header=F))
colnames(PRJEB5224_Count)<-c("match_id","read_count")
PRJEB5224_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB5224_Count$match_id, fixed = TRUE)
PRJEB5224_Count$match_id<-gsub("_kneaddata.fastq", "", PRJEB5224_Count$match_id, fixed = TRUE)
PRJEB5224_whole<-merge(PRJEB5224,PRJEB5224_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB5224_whole)


### LiuNN_2022	PRJNA731589 
PRJNA731589<-subset(clinical,study_name=="LiuNN_2022",drop=T)

PRJNA731589<-subset(PRJNA731589,sample_id!="HS_32",drop=T)
PRJNA731589<-subset(PRJNA731589,sample_id!="HS_1",drop=T)
PRJNA731589<-subset(PRJNA731589,sample_id!="CRC_5",drop=T)

PRJNA731589$match_id<-PRJNA731589$NCBI_accession
nrow(PRJNA731589)

PRJNA731589_Count<-unique(read.table("read_number/PRJNA731589_read_number.txt",header=F))
colnames(PRJNA731589_Count)<-c("match_id","read_count")
PRJNA731589_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA731589_Count$match_id, fixed = TRUE)
PRJNA731589_whole<-unique(merge(PRJNA731589,PRJNA731589_Count,by.x="match_id",by.y="match_id"))
nrow(PRJNA731589_whole)

### MetaCardis_2020_a	PRJEB37249

PRJEB37249<-subset(clinical,study_name=="MetaCardis_2020_a",drop=T)
PRJEB37249$match_id<-PRJEB37249$sample_id
nrow(PRJEB37249)

PRJEB37249_Count<-unique(read.table("read_number/PRJEB37249_read_number.txt",header=F))
colnames(PRJEB37249_Count)<-c("match_id","read_count")
PRJEB37249_Count$match_id<-gsub("_kneaddata.fastq", "", PRJEB37249_Count$match_id, fixed = TRUE)
PRJEB37249_whole<-merge(PRJEB37249,PRJEB37249_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB37249_whole)

### MetaCardis_2020_a	PRJEB38742

PRJEB38742<-subset(clinical,study_name=="MetaCardis_2020_a",drop=T)
PRJEB38742$match_id<-PRJEB38742$sample_id
nrow(PRJEB38742)

PRJEB38742_Count<-unique(read.table("read_number/PRJEB38742_read_number.txt",header=F))
colnames(PRJEB38742_Count)<-c("match_id","read_count")
PRJEB38742_Count$match_id<-gsub("_kneaddata.fastq", "", PRJEB38742_Count$match_id, fixed = TRUE)
PRJEB38742_whole<-merge(PRJEB38742,PRJEB38742_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB38742_whole)

## PasolliE_2019	PRJNA485056

PRJNA485056<-subset(clinical,study_name=="PasolliE_2019",drop=T)
PRJNA485056$match_id<-PRJNA485056$NCBI_accession
nrow(PRJNA485056)

PRJNA485056_Count<-unique(read.table("read_number/PRJNA485056_read_number.txt",header=F))
colnames(PRJNA485056_Count)<-c("match_id","read_count")
PRJNA485056_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA485056_Count$match_id, fixed = TRUE)
PRJNA485056_whole<-merge(PRJNA485056,PRJNA485056_Count,by.x="match_id",by.y="match_id")	
nrow(PRJNA485056_whole)


### QinJ_2012	PRJNA422434
PRJNA422434_id<-read.csv("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/2_pipe_qc_snp_sv_abun/2_submit_index/samples/unique/QinJ_2012_PRJNA422434.csv",header=T)
PRJNA422434_id<-PRJNA422434_id[,c(1,4)]
colnames(PRJNA422434_id)[1]<-"match_id"

PRJNA422434_clin<-subset(clinical,study_name=="QinJ_2012",drop=T)
PRJNA422434<-merge(PRJNA422434_clin,PRJNA422434_id,by.x="sample_id",by.y="sample_id")
nrow(PRJNA422434)

PRJNA422434_Count<-unique(read.table("read_number/PRJNA422434_read_number.txt",header=F))
colnames(PRJNA422434_Count)<-c("match_id","read_count")
PRJNA422434_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA422434_Count$match_id, fixed = TRUE)
PRJNA422434_whole<-merge(PRJNA422434,PRJNA422434_Count,by.x="match_id",by.y="match_id")
nrow(PRJNA422434_whole)


###QinN_2014	PRJEB6337

PRJEB6337<-subset(clinical,study_name=="QinN_2014",drop=T)
PRJEB6337$match_id<-PRJEB6337$sample_id
nrow(PRJEB6337)

PRJEB6337_Count<-unique(read.table("read_number/PRJEB6337_read_number.txt",header=F))
colnames(PRJEB6337_Count)<-c("match_id","read_count")
PRJEB6337_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB6337_Count$match_id, fixed = TRUE)
PRJEB6337_whole<-merge(PRJEB6337,PRJEB6337_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB6337_whole)


### RubelMA_2020	PRJNA547591

PRJNA547591<-subset(clinical,study_name=="RubelMA_2020",drop=T)
PRJNA547591$match_id<-PRJNA547591$NCBI_accession
nrow(PRJNA547591)

PRJNA547591_Count<-unique(read.table("read_number/PRJNA547591_read_number.txt",header=F))
colnames(PRJNA547591_Count)<-c("match_id","read_count")
PRJNA547591_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA547591_Count$match_id, fixed = TRUE)
PRJNA547591_whole<-merge(PRJNA547591,PRJNA547591_Count,by.x="match_id",by.y="match_id")	
nrow(PRJNA547591_whole)

### SchirmerM_2016	PRJNA319574

PRJNA319574<-subset(clinical,study_name=="SchirmerM_2016",drop=T)
PRJNA319574$match_id<-PRJNA319574$sample_id
nrow(PRJNA319574)

PRJNA319574_Count<-unique(read.table("read_number/PRJNA319574_read_number.txt",header=F))
colnames(PRJNA319574_Count)<-c("match_id","read_count")
PRJNA319574_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA319574_Count$match_id, fixed = TRUE)
PRJNA319574_whole<-merge(PRJNA319574,PRJNA319574_Count,by.x="match_id",by.y="match_id")	
nrow(PRJNA319574_whole)

### ThomasAM_2018a ThomasAM_2018b	PRJNA447983
PRJNA447983_id<-read.csv("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/2_pipe_qc_snp_sv_abun/2_submit_index/samples/unique/ThomasAM_2018ab_PRJNA447983.csv",header=T)
PRJNA447983_id<-PRJNA447983_id[,c(2,4)]
colnames(PRJNA447983_id)[1]<-"match_id"

PRJNA447983_clin<-subset(clinical,study_name=="ThomasAM_2018",drop=T)
PRJNA447983<-merge(PRJNA447983_clin,PRJNA447983_id,by.x="sample_id",by.y="sample_id")
nrow(PRJNA447983)

PRJNA447983_Count<-unique(read.table("read_number/PRJNA447983_read_number.txt",header=F))
colnames(PRJNA447983_Count)<-c("match_id","read_count")
PRJNA447983_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA447983_Count$match_id, fixed = TRUE)
PRJNA447983_whole<-merge(PRJNA447983,PRJNA447983_Count,by.x="match_id",by.y="match_id")
nrow(PRJNA447983_whole)


## VincentC_2016	PRJNA297252 (delete the replicated samples)
PRJNA297252<-subset(clinical,study_name=="VincentC_2016",drop=T)
PRJNA297252_select<-read.csv("clinical/PRJNA297252_select.csv",header=T)

PRJNA297252<-merge(PRJNA297252_select,PRJNA297252,by.x="sample_id",by.y="sample_id")
PRJNA297252$match_id<-PRJNA297252$sample_id
nrow(PRJNA297252)

##write.csv(PRJNA297252,"PRJNA297252.csv")

PRJNA297252_Count<-unique(read.table("read_number/PRJNA297252_read_number.txt",header=F))
colnames(PRJNA297252_Count)<-c("match_id","read_count")
PRJNA297252_Count$match_id<-gsub("_kneaddata.fastq", "", PRJNA297252_Count$match_id, fixed = TRUE)
PRJNA297252_whole<-merge(PRJNA297252,PRJNA297252_Count,by.x="match_id",by.y="match_id")	
nrow(PRJNA297252_whole)


### VogtmannE_2016	PRJEB12449
PRJEB12449<-subset(clinical,study_name=="VogtmannE_2016",drop=T)
PRJEB12449$match_id<-PRJEB12449$sample_id
nrow(PRJEB12449)

PRJEB12449_Count<-unique(read.table("read_number/PRJEB12449_read_number.txt",header=F))
colnames(PRJEB12449_Count)<-c("match_id","read_count")
PRJEB12449_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB12449_Count$match_id, fixed = TRUE)
PRJEB12449_whole<-merge(PRJEB12449,PRJEB12449_Count,by.x="match_id",by.y="match_id")
nrow(PRJEB12449_whole)

### WirbelJ_2018	PRJEB27928 
PRJEB27928<-subset(clinical,study_name=="WirbelJ_2018",drop=T)
PRJEB27928<-PRJEB27928[,-21]

PRJEB27928_id<-read.csv("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/2_pipe_qc_snp_sv_abun/2_submit_index/samples/unique/WirbelJ_2018_PRJEB27928.csv",header=T)
PRJEB27928_id<-PRJEB27928_id[,c(1,4)]

PRJEB27928<-merge(PRJEB27928,PRJEB27928_id,by.x="sample_id",by.y="sample_id")
nrow(PRJEB27928)

PRJEB27928$match_id<-PRJEB27928$NCBI_accession

PRJEB27928_Count<-unique(read.table("read_number/PRJEB27928_read_number.txt",header=F))
colnames(PRJEB27928_Count)<-c("match_id","read_count")
PRJEB27928_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB27928_Count$match_id, fixed = TRUE)

PRJEB27928_whole<-merge(PRJEB27928,PRJEB27928_Count,by.x="match_id",by.y="match_id")
nrow(PRJEB27928_whole)


## XieH_2016	PRJEB9576

PRJEB9576<-subset(clinical,study_name=="XieH_2016",drop=T)
PRJEB9576$match_id<-PRJEB9576$NCBI_accession
nrow(PRJEB9576)

PRJEB9576_Count<-unique(read.table("read_number/PRJEB9576_read_number.txt",header=F))
colnames(PRJEB9576_Count)<-c("match_id","read_count")
PRJEB9576_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB9576_Count$match_id, fixed = TRUE)
PRJEB9576_whole<-merge(PRJEB9576,PRJEB9576_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB9576_whole)

### YachidaS_2019 PRJDB4176
PRJDB4176_id<-read.csv("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/2_pipe_qc_snp_sv_abun/2_submit_index/samples/unique/YachidaS_2019_PRJDB4176.csv",header=T)
PRJDB4176_id<-PRJDB4176_id[,c(2,4)]
colnames(PRJDB4176_id)[1]<-"match_id"

PRJDB4176_clin<-subset(clinical,study_name=="YachidaS_2019",drop=T)
PRJDB4176<-merge(PRJDB4176_clin,PRJDB4176_id,by.x="subject_id",by.y="subject_id")
nrow(PRJDB4176)

PRJDB4176_Count<-unique(read.table("read_number/PRJDB4176_read_number.txt",header=F))
colnames(PRJDB4176_Count)<-c("match_id","read_count")
PRJDB4176_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJDB4176_Count$match_id, fixed = TRUE)
PRJDB4176_whole<-merge(PRJDB4176,PRJDB4176_Count,by.x="match_id",by.y="match_id")
nrow(PRJDB4176_whole)

### YuJ_2015	PRJEB10878
PRJEB10878<-subset(clinical,study_name=="YuJ_2015",drop=T)
PRJEB10878$match_id<-PRJEB10878$NCBI_accession
nrow(PRJEB10878)

PRJEB10878_Count<-unique(read.table("read_number/PRJEB10878_read_number.txt",header=F))
colnames(PRJEB10878_Count)<-c("match_id","read_count")
PRJEB10878_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB10878_Count$match_id, fixed = TRUE)
PRJEB10878_whole<-merge(PRJEB10878,PRJEB10878_Count,by.x="match_id",by.y="match_id")
nrow(PRJEB10878_whole)

### ZeeviD_2015	PRJEB11532

PRJEB11532<-subset(clinical,study_name=="ZeeviD_2015",drop=T)
PRJEB11532$match_id<-PRJEB11532$sample_id
nrow(PRJEB11532)

PRJEB11532_Count<-unique(read.table("read_number/PRJEB11532_read_number.txt",header=F))
colnames(PRJEB11532_Count)<-c("match_id","read_count")
PRJEB11532_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB11532_Count$match_id, fixed = TRUE)
PRJEB11532_Count$match_id<-gsub("_kneaddata.fastq", "", PRJEB11532_Count$match_id, fixed = TRUE)
PRJEB11532_whole<-merge(PRJEB11532,PRJEB11532_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB11532_whole)

### ZellerG_2014	PRJEB6070
PRJEB6070<-subset(clinical,study_name=="ZellerG_2014",drop=T)
PRJEB6070$match_id<-PRJEB6070$sample_id
nrow(PRJEB6070)

PRJEB6070_Count<-unique(read.table("read_number/PRJEB6070_read_number.txt",header=F))
colnames(PRJEB6070_Count)<-c("match_id","read_count")
PRJEB6070_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB6070_Count$match_id, fixed = TRUE)
PRJEB6070_whole<-merge(PRJEB6070,PRJEB6070_Count,by.x="match_id",by.y="match_id")
nrow(PRJEB6070_whole)

## ZhuF_2020	PRJEB29127

PRJEB29127<-subset(clinical,study_name=="ZhuF_2020",drop=T)
PRJEB29127$match_id<-PRJEB29127$NCBI_accession
nrow(PRJEB29127)

PRJEB29127_Count<-unique(read.table("read_number/PRJEB29127_read_number.txt",header=F))
colnames(PRJEB29127_Count)<-c("match_id","read_count")
PRJEB29127_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB29127_Count$match_id, fixed = TRUE)
PRJEB29127_whole<-merge(PRJEB29127,PRJEB29127_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB29127_whole)


### PehrssonE_2016	PRJNA300541

PRJNA300541<-subset(clinical,study_name=="PehrssonE_2016",drop=T)
PRJNA300541$match_id<-PRJNA300541$sample_id
nrow(PRJNA300541)

PRJNA300541_Count<-unique(read.table("read_number/PRJNA300541_read_number.txt",header=F))
colnames(PRJNA300541_Count)<-c("match_id","read_count")
PRJNA300541_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA300541_Count$match_id, fixed = TRUE)
PRJNA300541_whole<-merge(PRJNA300541,PRJNA300541_Count,by.x="match_id",by.y="match_id")	
nrow(PRJNA300541_whole)



### Fernandez_Calvet_2024  PRJNA834801

PRJNA834801<-subset(clinical,study_name=="Fern\u00e1ndez_Calvet_2024",drop=T)
PRJNA834801$match_id<-PRJNA834801$NCBI_accession
nrow(PRJNA834801)

PRJNA834801_Count<-unique(read.table("read_number/PRJNA834801_read_number.txt",header=F))
colnames(PRJNA834801_Count)<-c("match_id","read_count")
PRJNA834801_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA834801_Count$match_id, fixed = TRUE)
PRJNA834801_whole<-merge(PRJNA834801,PRJNA834801_Count,by.x="match_id",by.y="match_id")	
nrow(PRJNA834801_whole)


### PiccinnoG_2025_cohort1  PRJEB72523

PRJEB72523<-subset(clinical,study_name=="PiccinnoG_2025_cohort1",drop=T)
PRJEB72523$match_id<-PRJEB72523$NCBI_accession
nrow(PRJEB72523)

PRJEB72523_Count<-unique(read.table("read_number/PRJEB72523_read_number.txt",header=F))
colnames(PRJEB72523_Count)<-c("match_id","read_count")
PRJEB72523_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB72523_Count$match_id, fixed = TRUE)
PRJEB72523_whole<-merge(PRJEB72523,PRJEB72523_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB72523_whole)


### PiccinnoG_2025_cohort2  PRJEB72524

PRJEB72524<-subset(clinical,study_name=="PiccinnoG_2025_cohort2",drop=T)
PRJEB72524$match_id<-PRJEB72524$NCBI_accession
nrow(PRJEB72524)

PRJEB72524_Count<-unique(read.table("read_number/PRJEB72524_read_number.txt",header=F))
colnames(PRJEB72524_Count)<-c("match_id","read_count")
PRJEB72524_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB72524_Count$match_id, fixed = TRUE)
PRJEB72524_whole<-merge(PRJEB72524,PRJEB72524_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB72524_whole)



### PiccinnoG_2025_cohort3  PRJEB72525

PRJEB72525<-subset(clinical,study_name=="PiccinnoG_2025_cohort3",drop=T)
PRJEB72525$match_id<-PRJEB72525$NCBI_accession
nrow(PRJEB72525)

PRJEB72525_Count<-unique(read.table("read_number/PRJEB72525_read_number.txt",header=F))
colnames(PRJEB72525_Count)<-c("match_id","read_count")
PRJEB72525_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB72525_Count$match_id, fixed = TRUE)
PRJEB72525_whole<-merge(PRJEB72525,PRJEB72525_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB72525_whole)


### PiccinnoG_2025_cohort4  PRJEB72526

PRJEB72526<-subset(clinical,study_name=="PiccinnoG_2025_cohort4",drop=T)
PRJEB72526$match_id<-PRJEB72526$NCBI_accession
nrow(PRJEB72526)

PRJEB72526_Count<-unique(read.table("read_number/PRJEB72526_read_number.txt",header=F))
colnames(PRJEB72526_Count)<-c("match_id","read_count")
PRJEB72526_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB72526_Count$match_id, fixed = TRUE)
PRJEB72526_whole<-merge(PRJEB72526,PRJEB72526_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB72526_whole)


### PiccinnoG_2025_cohort5  PRJNA1237248

PRJNA1237248<-subset(clinical,study_name=="PiccinnoG_2025_cohort5",drop=T)
PRJNA1237248$match_id<-PRJNA1237248$NCBI_accession
nrow(PRJNA1237248)

PRJNA1237248_Count<-unique(read.table("read_number/PRJNA1237248_read_number.txt",header=F))
colnames(PRJNA1237248_Count)<-c("match_id","read_count")
PRJNA1237248_Count$match_id<-gsub("_kneaddata.fastq", "", PRJNA1237248_Count$match_id, fixed = TRUE)
PRJNA1237248_whole<-merge(PRJNA1237248,PRJNA1237248_Count,by.x="match_id",by.y="match_id")	
nrow(PRJNA1237248_whole)


### YangJ_2020  PRJNA429097

PRJNA429097<-subset(clinical,study_name=="YangJ_2020",drop=T)
PRJNA429097$match_id<-PRJNA429097$NCBI_accession
nrow(PRJNA429097)

PRJNA429097_Count<-unique(read.table("read_number/PRJNA429097_read_number.txt",header=F))
colnames(PRJNA429097_Count)<-c("match_id","read_count")
PRJNA429097_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA429097_Count$match_id, fixed = TRUE)
PRJNA429097_whole<-merge(PRJNA429097,PRJNA429097_Count,by.x="match_id",by.y="match_id")	
nrow(PRJNA429097_whole)


### YangY_2021  PRJNA763023

PRJNA763023<-subset(clinical,study_name=="YangY_2021",drop=T)
PRJNA763023$match_id<-PRJNA763023$NCBI_accession
nrow(PRJNA763023)

PRJNA763023_Count<-unique(read.table("read_number/PRJNA763023_read_number.txt",header=F))
colnames(PRJNA763023_Count)<-c("match_id","read_count")
PRJNA763023_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA763023_Count$match_id, fixed = TRUE)
PRJNA763023_whole<-merge(PRJNA763023,PRJNA763023_Count,by.x="match_id",by.y="match_id")	
nrow(PRJNA763023_whole)


### ClasenF_2025  PRJEB79944

PRJEB79944<-subset(clinical,study_name=="ClasenF_2025",drop=T)
PRJEB79944$match_id<-PRJEB79944$NCBI_accession
nrow(PRJEB79944)

PRJEB79944_Count<-unique(read.table("read_number/PRJEB79944_read_number.txt",header=F))
colnames(PRJEB79944_Count)<-c("match_id","read_count")
PRJEB79944_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJEB79944_Count$match_id, fixed = TRUE)
PRJEB79944_whole<-merge(PRJEB79944,PRJEB79944_Count,by.x="match_id",by.y="match_id")	
nrow(PRJEB79944_whole)


### JoSY_2022  PRJNA743718

PRJNA743718<-subset(clinical,study_name=="JoSY_2022",drop=T)
PRJNA743718$match_id<-PRJNA743718$NCBI_accession
nrow(PRJNA743718)

PRJNA743718_Count<-unique(read.table("read_number/PRJNA743718_read_number.txt",header=F))
colnames(PRJNA743718_Count)<-c("match_id","read_count")
PRJNA743718_Count$match_id<-gsub("_1_kneaddata_paired_1.fastq", "", PRJNA743718_Count$match_id, fixed = TRUE)
PRJNA743718_whole<-merge(PRJNA743718,PRJNA743718_Count,by.x="match_id",by.y="match_id")	
nrow(PRJNA743718_whole)



## merge the clinical informations 

clin_merge<-rbind(PRJEB39223_whole,PRJEB17632_whole,PRJEB7774_whole,PRJNA400072_whole,
PRJNA385949_whole,PRJNA491335_whole,HMP_2012_whole,PRJNA389280_whole,PRJNA398089_whole,
PRJEB21528_whole,PRJEB1786_whole,PRJEB36820_whole,PRJEB5224_whole,
PRJNA731589_whole,PRJEB37249_whole,PRJEB38742_whole,PRJNA485056_whole,PRJNA422434_whole,
PRJEB6337_whole,PRJNA547591_whole,PRJNA319574_whole,PRJNA447983_whole,
PRJEB12449_whole,PRJEB27928_whole,PRJEB9576_whole,PRJDB4176_whole,PRJEB10878_whole,
PRJEB11532_whole,PRJEB6070_whole,PRJEB29127_whole,PRJNA300541_whole,PRJNA834801_whole,
PRJEB72523_whole,PRJEB72524_whole,PRJEB72525_whole,PRJEB72526_whole,PRJNA1237248_whole,
PRJNA429097_whole,PRJNA763023_whole,PRJEB79944_whole,
PRJNA743718_whole)

##clin_merge$hscrp <- as.numeric(replace(clin_merge$hscrp, clin_merge$hscrp =="<5", NA))
##ThomasAM_2018<-subset(clin_merge,study_name=="ThomasAM_2018",drop=T)


## select stool sample 
table(clin_merge$body_site)
clin_merge <- subset(clin_merge,body_site=="stool",drop=T)
unique(clin_merge$study_name)

##select samples age adult or senior
table(clin_merge$age_category)

clin_merge <- subset(clin_merge,as.numeric(age)>=18 | age_category=="adult" | age_category=="senior",drop=T)
unique(clin_merge$study_name)

##select no pregnant 
table(clin_merge$pregnant)
clin_merge<- subset(clin_merge,pregnant!="yes" | is.na(pregnant)==T,drop=T)
unique(clin_merge$study_name)

min(as.numeric(clin_merge$age),na.rm=T)



##########################################
##  define the disease variable 
## IBD  
## T2D  
## CRC
## severe obesity  
## CAD (chronic CAD,acute CAD,chronic CAD;HF)
## ACVD
## adenoma
## cirrhosis
## Mets
## IGT
## PD

table(clin_merge$study_name)



nr<-nrow(clin_merge)
CRC_study_name<-unique(subset(clin_merge,disease=="CRC",drop=T)$study_name)
clin_merge_CRC <- clin_merge %>%
  filter(study_name %in% CRC_study_name)
table(clin_merge_CRC$study_name,clin_merge_CRC$disease)


CRC_study <- c("FengQ_2015", "LiuNN_2022", "ThomasAM_2018", 
               "VogtmannE_2016", "WirbelJ_2018", "YachidaS_2019",
               "YuJ_2015", "PiccinnoG_2025_cohort3", "PiccinnoG_2025_cohort4",
               "PiccinnoG_2025_cohort6", "YangJ_2020", "YangY_2021", 
               "ZellerG_2014")

clin_merge$CRC <- NA
idx <- clin_merge$study_name %in% CRC_study & !is.na(clin_merge$disease)
clin_merge$CRC[idx] <- ifelse(grepl("(^|;\\s*)CRC($|;\\s*)", clin_merge$disease[idx]), 1,
                       ifelse(grepl("healthy", clin_merge$disease[idx]), 0, NA))

table(clin_merge$CRC,clin_merge$study_name)



## adenoma

nr<-nrow(clin_merge)
adenoma_study_name<-unique(subset(clin_merge,disease=="adenoma",drop=T)$study_name)
clin_merge_adenoma <- clin_merge %>%
  filter(study_name %in% adenoma_study_name)
table(clin_merge_adenoma$study_name,clin_merge_adenoma$disease)

adenoma_study <- c("FengQ_2015","PiccinnoG_2025_cohort5",
                  "ThomasAM_2018", "YachidaS_2019", "ZellerG_2014")

clin_merge$adenoma <- NA
idx <- clin_merge$study_name %in% adenoma_study & !is.na(clin_merge$disease)
clin_merge$adenoma[idx] <- ifelse(grepl("(^|;\\s*)adenoma($|;\\s*)", clin_merge$disease[idx]), 1,
                       ifelse(grepl("healthy", clin_merge$disease[idx]), 0, NA))

table(clin_merge$adenoma,clin_merge$study_name)




#######################################################################
### IBD 

IBD_study_name<-unique(subset(clin_merge,disease=="IBD",drop=T)$study_name)

clin_merge_IBD <- clin_merge %>%
  filter(study_name %in% IBD_study_name)
table(clin_merge_IBD$study_name,clin_merge_IBD$disease)

IBD_study <-c("FranzosaEA_2018","HallAB_2017","HMP_2019_ibdmdb")

clin_merge$IBD <- NA
idx <- clin_merge$study_name %in% IBD_study & !is.na(clin_merge$disease)
clin_merge$IBD[idx] <- ifelse(grepl("(^|;\\s*)IBD($|;\\s*)", clin_merge$disease[idx]), 1,
                       ifelse(grepl("healthy", clin_merge$disease[idx]), 0, NA))

table(clin_merge$IBD,clin_merge$study_name)



ACVD_study_name<-unique(subset(clin_merge,disease=="ACVD",drop=T)$study_name)

clin_merge_ACVD <- clin_merge %>%
  filter(study_name %in% ACVD_study_name)
table(clin_merge_ACVD$study_name,clin_merge_ACVD$disease)

ACVD_study <- c("JieZ_2017")

clin_merge$ACVD <- NA
idx <- clin_merge$study_name %in% ACVD_study & !is.na(clin_merge$disease)

clin_merge$ACVD[idx] <- ifelse(grepl("(^|;\\s*)ACVD($|;\\s*)", clin_merge$disease[idx]), 1,
                       ifelse(grepl("healthy", clin_merge$disease[idx]), 0, NA))

table(clin_merge$ACVD,clin_merge$study_name)



T2D_study_name<-unique(subset(clin_merge,disease=="T2D",drop=T)$study_name)

clin_merge_T2D <- clin_merge %>%
  filter(study_name %in% T2D_study_name)
table(clin_merge_T2D$study_name,clin_merge_T2D$disease)

T2D_study <- c("KarlssonFH_2013","MetaCardis_2020_a","QinJ_2012","YuJ_2015")

clin_merge$T2D <- NA
idx <- clin_merge$study_name %in% T2D_study & !is.na(clin_merge$disease)
clin_merge$T2D[idx] <- ifelse(grepl("(^|;\\s*)T2D($|;\\s*)", clin_merge$disease[idx]), 1,
                       ifelse(grepl("healthy", clin_merge$disease[idx]), 0, NA))

table(clin_merge$T2D,clin_merge$study_name)



severe_obesity_study_name<-unique(subset(clin_merge,disease=="severe obesity",drop=T)$study_name)

clin_merge_severe_obesity <- clin_merge %>%
  filter(study_name %in% severe_obesity_study_name)
table(clin_merge_severe_obesity$study_name,clin_merge_severe_obesity$disease)

severe_obesity_study <- c("MetaCardis_2020_a")

clin_merge$severe_obesity <- NA
idx <- clin_merge$study_name %in% severe_obesity_study & !is.na(clin_merge$disease)
clin_merge$severe_obesity[idx] <- ifelse(grepl("(^|;\\s*)severe obesity($|;\\s*)", clin_merge$disease[idx]), 1,
                       ifelse(grepl("healthy", clin_merge$disease[idx]), 0, NA))

table(clin_merge$severe_obesity,clin_merge$study_name)



CAD_study_name<-unique(subset(clin_merge,disease=="chronic CAD" | disease=="acute CAD" |
              disease=="chronic CAD;HF",drop=T)$study_name)

clin_merge_CAD <- clin_merge %>%
  filter(study_name %in% CAD_study_name)
table(clin_merge_CAD$study_name,clin_merge_CAD$disease)

CAD_study <- c("MetaCardis_2020_a")

clin_merge$CAD <- NA
idx <- clin_merge$study_name %in% CAD_study & !is.na(clin_merge$disease)
clin_merge$CAD[idx] <- ifelse(grepl("CAD", clin_merge$disease[idx]), 1,
                       ifelse(grepl("healthy", clin_merge$disease[idx]), 0, NA))

table(clin_merge$CAD,clin_merge$study_name)


## cirrhosis
cirrhosis_study_name<-unique(subset(clin_merge,disease=="cirrhosis",drop=T)$study_name)

clin_merge_cirrhosis <- clin_merge %>%
  filter(study_name %in% cirrhosis_study_name)
table(clin_merge_cirrhosis$study_name,clin_merge_cirrhosis$disease)

cirrhosis_study <- c("QinN_2014")

clin_merge$cirrhosis <- NA
idx <- clin_merge$study_name %in% cirrhosis_study & !is.na(clin_merge$disease)
clin_merge$cirrhosis[idx] <- ifelse(grepl("(^|;\\s*)cirrhosis($|;\\s*)", clin_merge$disease[idx]), 1,
                       ifelse(grepl("healthy", clin_merge$disease[idx]), 0, NA))

table(clin_merge$cirrhosis,clin_merge$study_name)


## MetS
MetS_study_name<-unique(subset(clin_merge,disease=="MetS",drop=T)$study_name)

clin_merge_MetS <- clin_merge %>%
  filter(study_name %in% MetS_study_name)
table(clin_merge_MetS$study_name,clin_merge_MetS$disease)

MetS_study <- c("MetaCardis_2020_a")

clin_merge$MetS <- NA
idx <- clin_merge$study_name %in% MetS_study & !is.na(clin_merge$disease)
clin_merge$MetS[idx] <- ifelse(grepl("(^|;\\s*)MetS($|;\\s*)", clin_merge$disease[idx]), 1,
                       ifelse(grepl("healthy", clin_merge$disease[idx]), 0, NA))

table(clin_merge$MetS,clin_merge$study_name)



###########################################################
## PD   

PD_study_name<-unique(subset(clin_merge,disease=="PD",drop=T)$study_name)

clin_merge_PD <- clin_merge %>%
  filter(study_name %in% PD_study_name)
table(clin_merge_PD$study_name,clin_merge_PD$disease)

PD_study <-PD_study_name

clin_merge$PD <- NA
idx <- clin_merge$study_name %in% PD_study & !is.na(clin_merge$disease)
clin_merge$PD[idx] <- ifelse(grepl("PD", clin_merge$disease[idx]), 1,
                       ifelse(grepl("healthy", clin_merge$disease[idx]), 0, NA))

table(clin_merge$PD,clin_merge$study_name)



##################### recode gender 

clin_merge$gender_code<-NA
nr<-nrow(clin_merge)

for(i in 1:nr){
     if(is.na(clin_merge$gender[i])==T){clin_merge$gender_code[i]=NA}
     if(is.na(clin_merge$gender[i])==F){
      if (clin_merge$gender[i]=="female"){clin_merge$gender_code[i]=1}
      else if (clin_merge$gender[i]=="male"){clin_merge$gender_code[i]=0}
      else {clin_merge$gender_code[i]=NA}
    }
    else {clin_merge$gender_code[i]=NA}
 }


##################### define intercontinental (region)
## ISO 3166-1 alpha-3 -> region. One-to-one map, no silent NA fallback.
## Study convention: ISR and TUR are grouped into Europe; state this in Methods.
region_map <- c(
  AUT="Europe", DEU="Europe", DNK="Europe", ESP="Europe", EST="Europe",
  FIN="Europe", FRA="Europe", GBR="Europe", HUN="Europe", IRL="Europe",
  ISL="Europe", ITA="Europe", NLD="Europe", NOR="Europe", SVK="Europe",
  SWE="Europe", CZE="Europe", UK ="Europe", ISR="Europe", TUR="Europe",
  CAN="Americas", PER="Americas", SLV="Americas", USA="Americas",
  CHN="Asia", JPN="Asia", KAZ="Asia", KOR="Asia",
  CMR="Africa", MDG="Africa"
)

country_chr <- as.character(clin_merge$country)
clin_merge$inter <- unname(region_map[country_chr])

## Coverage self-check: report unmapped country codes instead of silently setting NA
unmapped <- setdiff(sort(unique(country_chr[!is.na(country_chr)])), names(region_map))
if (length(unmapped) > 0) {
  warning("Unmapped country code(s): ", paste(unmapped, collapse = ", "))
}

## Map self-check: unique country codes, canonical region names only
stopifnot(!anyDuplicated(names(region_map)),
          all(region_map %in% c("Africa", "Americas", "Asia", "Europe")))

table(clin_merge$country, clin_merge$inter, useNA = "ifany")


table(clin_merge$inter)
table(clin_merge$country)
table(clin_merge$country,clin_merge$study_name)

clin_merge$mean_press<-(clin_merge$systolic_p + 2 * clin_merge$dyastolic_p) / 3


clin_merge <- subset(clin_merge,as.numeric(age)>=18 | age_category=="adult" | age_category=="senior",drop=T)
unique(clin_merge$study_name)

## Age
mean(as.numeric(clin_merge$age), na.rm = T)
se(as.numeric(clin_merge$age))
min(as.numeric(clin_merge$age), na.rm = T)
max(as.numeric(clin_merge$age), na.rm = T)


dir.create("../01.cleanData/phen", recursive = TRUE, showWarnings = FALSE)
write.csv(clin_merge,"clinical/clinical_whole.csv",row.names=F)
write.table(clin_merge, "F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/01.cleanData/phen/Clinical_basic.tsv",sep = '\t')



##########################################################################
### The shapiro.test or Kolmogorov-Smirnov for continues varaibles 

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

clinical<-read.csv("00.rawData/clinical/clinical_whole.csv",header=T)
clinical$study_name <- normalize_study_name(clinical$study_name)

## calculate the sample size and study_name
dataset_variable <- clinical %>%
    group_by(study_name) %>%
    summarise(
        BMI= sum(!is.na(BMI)),
        age = sum(!is.na(age)),
        hba1c = sum(!is.na(hba1c)),
        fasting_glucose = sum(!is.na(fasting_glucose)),
        fasting_insulin = sum(!is.na(fasting_insulin)),
        triglycerides = sum(!is.na(triglycerides)),
        hdl = sum(!is.na(hdl)),
        ldl = sum(!is.na(ldl)),
        cholesterol = sum(!is.na(cholesterol)),
        mean_press = sum(!is.na(mean_press)),
        creatinine = sum(!is.na(creatinine)),
        alt = sum(!is.na(alt)),
        ast = sum(!is.na(ast)),
        total_bilirubin = sum(!is.na(total_bilirubin))
    )

dataset_variable <- dataset_variable %>%
    mutate(across(BMI:total_bilirubin, ~ ifelse(. >= 40, 1, 0)))

print(dataset_variable)

write.csv(dataset_variable,"00.rawData/clinical/dataset_variable.csv",row.names=F)

dataset_variable_norm <- dataset_variable %>%
    rename_with(~paste0(., "_norm"))

write.csv(dataset_variable_norm,"00.rawData/clinical/dataset_variable_norm.csv",row.names=F)



dataset_variable$variable <- apply(dataset_variable, 1, function(x) {paste(colnames(dataset_variable)[which(x==1)], collapse = ",")})
dataset_variable<-dataset_variable[,c("study_name","variable")]

datasets<-dataset_variable$study_name
study_num<-length(datasets)


results<-matrix(nrow=1,ncol=6)
colnames(results)<-c("Variable","Datasets","Statistic_SW","P_SW","Statistic_KS","P_KS")

for (studyname in datasets){
  # studyname<-"MetaCardis_2020_a"
  dataset_clinical<-subset(clinical,study_name==studyname,drop=T)
  covariate<-unlist(strsplit(dataset_variable[dataset_variable$study_name==studyname,]$variable,","))

  for (var in covariate){
     # var<-"BMI"
      data<-dataset_clinical[,var]

      if(length(data)>3 & length(data)<5000){
       ## Shapiro-Wilk
        sw_result<-shapiro.test(data)
        sw_sta<-sw_result$statistic
        sw_p<-sw_result$p.value
        }
       else{
        sw_sta<-NA
        sw_p<-NA
        }

       ## Kolmogorov-Smirnov test 
       ks_result<-ks.test(data, "pnorm", mean = mean(data, na.rm = TRUE), sd = sd(data, na.rm = TRUE))
       ks_sta<-ks_result$statistic
       ks_p<-ks_result$p.value
 
       temp<-c(studyname,var,sw_sta,sw_p,ks_sta,ks_p)
       results<-rbind(results,temp)
    }
 }


write.csv(results,"00.rawData/clinical/normality_test.csv")



##############################################################################
### transform the continues variables into 
##Before association analysis, all continuous variables were standardized to follow a standard normal distribution (N~(0, 1))
## using empirical normal quantile transformation


library(dplyr)

empirical_normal_transform <- function(x) {
  non_na_indices <- !is.na(x)  
  ranks <- rank(x[non_na_indices], ties.method = "average") 
  quantiles <- (ranks - 0.5) / sum(non_na_indices) 
  
  transformed <- rep(NA, length(x)) 
  transformed[non_na_indices] <- qnorm(quantiles, mean = 0, sd = 1)
  return(transformed)
}



setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")
setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

clinical<-read.csv("00.rawData/clinical/clinical_whole.csv",header=T)
clinical$study_name <- normalize_study_name(clinical$study_name)


datasets<-dataset_variable$study_name
study_num<-length(datasets)


transform<-matrix(nrow=1,ncol=14)
colnames(transform)<-c("BMI","age","hba1c","fasting_glucose","fasting_insulin",
"triglycerides","hdl","ldl","cholesterol","mean_press","creatinine",
"alt","ast","total_bilirubin")


for (studyname in datasets){
   # studyname<-"MetaCardis_2020_a"
   dataset_clinical<-subset(clinical,study_name==studyname,drop=T)[,c("match_id","BMI","age",
                                          "hba1c","fasting_glucose","fasting_insulin","triglycerides",
                                          "hdl","ldl","cholesterol","mean_press","creatinine",
                                          "alt","ast","total_bilirubin")]

    row.names(dataset_clinical)<-dataset_clinical$match_id
    dataset_clinical<-dataset_clinical[,-1]
  
    transformed_data <- dataset_clinical %>%
     mutate(across(everything(), empirical_normal_transform))

    transform<-rbind(transform,transformed_data)
  }

transform$match_id<-rownames(transform)

colnames(transform)<-c("BMI_norm","age_norm","hba1c_norm","fasting_glucose_norm","fasting_insulin_norm",
"triglycerides_norm","hdl_norm","ldl_norm","cholesterol_norm","mean_press_norm","creatinine_norm",
"alt_norm","ast_norm","total_bilirubin_norm","match_id")

clinical_whole<-merge(clinical,transform,by.x="match_id",by.y="match_id")


################################################################################
## One sample per donor, among samples with clinical, SV and abundance profiles.
## Select after the existing normalization to preserve the current analysis values.
## Both s02 scripts read the retained sample list exported below.
dsv_sample_ids <- as.character(data.table::fread(
  "00.rawData/SV/Whole_dsgv_20260921.csv", select = 1L)[[1]])
vsv_sample_ids <- as.character(data.table::fread(
  "00.rawData/SV/Whole_vsgv_20260921.csv", select = 1L)[[1]])
abundance_columns <- names(data.table::fread(
  "00.rawData/taxonomic_abundance/All_sample_bracken_frac.txt", nrows = 0L))
abundance_columns <- setdiff(abundance_columns, "taxonomy_lvl")
abundance_sample_ids <- sub("\\._species\\.bracken\\.txt_frac$", "", abundance_columns[-c(1, 2)])
eligible_sample_ids <- Reduce(intersect, list(clinical_whole$match_id,
  dsv_sample_ids, vsv_sample_ids, abundance_sample_ids))
clinical_sub <- clinical_whole[match(eligible_sample_ids, clinical_whole$match_id), , drop = FALSE]
donor_audit_label <- "s00"
## Retain one sample per donor within each study.
## Input: clinical_sub and donor_audit_label; preserve the original row order.
clinical_sub <- clinical_sub[clinical_sub$study_name != "HallAB_2017", , drop = FALSE]
donor_info <- clinical_sub
subject_id <- trimws(as.character(donor_info$subject_id))
missing_id <- is.na(subject_id) | tolower(subject_id) %in% c("", "na", "unknown", "not available")
donor_info$subject_id_missing <- missing_id
donor_info$donor_key <- paste(donor_info$study_name,
  ifelse(missing_id, paste0("sample:", donor_info$match_id), subject_id), sep = "::")

## Collection day, original visit number, then original sample sequence.
collection_day <- suppressWarnings(as.numeric(donor_info$days_from_first_collection))
visit_number <- suppressWarnings(as.numeric(donor_info$visit_number))
sample_order <- rep(NA_real_, nrow(donor_info))
li_samples <- donor_info$study_name == "LiJ_2014" & grepl("_[0-9]+$", donor_info$sample_id)
sample_order[li_samples] <- as.numeric(sub(".*_([0-9]+)$", "\\1", donor_info$sample_id[li_samples]))
costea_samples <- donor_info$study_name == "CosteaPI_2017" &
  grepl("_[0-9]+_[0-9]+_[0-9]+$", donor_info$sample_id)
sample_order[costea_samples] <- as.numeric(sub(".*_[0-9]+_([0-9]+)_[0-9]+$", "\\1",
  donor_info$sample_id[costea_samples]))
german_samples <- donor_info$study_name == "CosteaPI_2017" &
  grepl("^alien_11_[0-9]+_0$", donor_info$sample_id)
sample_order[german_samples] <- as.numeric(sub("^alien_11_([0-9]+)_0$", "\\1",
  donor_info$sample_id[german_samples]))

donor_info$retained <- FALSE
donor_info$selection_rule <- "single_record"
donor_info$chronological_earliest_verified <- NA
donor_info$ordering_proxy <- sample_order
for (donor_rows in split(seq_len(nrow(donor_info)), donor_info$donor_key)) {
  selected_row <- donor_rows[1]
  if (length(donor_rows) > 1L) {
    if (all(is.finite(collection_day[donor_rows]))) {
      selected_row <- donor_rows[order(collection_day[donor_rows], donor_info$match_id[donor_rows])][1]
      selection_rule <- "earliest_collection_day"
    } else if (all(is.finite(visit_number[donor_rows]))) {
      selected_row <- donor_rows[order(visit_number[donor_rows], donor_info$match_id[donor_rows])][1]
      selection_rule <- "first_visit_number"
    } else if (all(is.finite(sample_order[donor_rows]))) {
      selected_row <- donor_rows[order(sample_order[donor_rows], donor_info$match_id[donor_rows])][1]
      selection_rule <- "original_sequence_proxy_time_unverified"
    } else {
      selection_rule <- "original_row_order_proxy_time_unverified"
    }
    donor_info$selection_rule[donor_rows] <- selection_rule
    donor_info$chronological_earliest_verified[donor_rows] <- selection_rule == "earliest_collection_day"
  }
  donor_info$retained[selected_row] <- TRUE
}

## Keep an audit of both retained and excluded records.
stopifnot(!anyDuplicated(donor_info$donor_key[donor_info$retained]),
          all(donor_info$retained[missing_id]))
donor_selection <- donor_info[, c("match_id", "study_name", "sample_id", "subject_id",
  "subject_id_missing", "donor_key", "days_from_first_collection", "visit_number",
  "ordering_proxy", "retained", "selection_rule", "chronological_earliest_verified")]
dir.create("01.cleanData/donor_selection", recursive = TRUE, showWarnings = FALSE)
write.csv(donor_selection, file.path("01.cleanData/donor_selection",
  paste0(donor_audit_label, "_all_samples.csv")), row.names = FALSE)
clinical_sub <- clinical_sub[donor_info$retained, , drop = FALSE]
clin_abun_inter <- clinical_sub$match_id

clinical_whole <- clinical_sub
write.csv(donor_selection[donor_selection$retained, , drop = FALSE],
  "01.cleanData/donor_selection/s00_retained_samples.csv", row.names = FALSE)
dir.create("01.cleanData/phen", recursive = TRUE, showWarnings = FALSE)
################################################################################

write.csv(clinical_whole,"00.rawData/clinical/clinical_whole_add_norm.csv",row.names=F)

write.table(clinical_whole, "F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/01.cleanData/phen/Clinical_basic.tsv",sep = '\t')




################################################################################################
### The shapiro.test or Kolmogorov-Smirnov for continues varaibles  (after normalization)

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

clinical<-read.csv("00.rawData/clinical/clinical_whole_add_norm.csv",header=T)
clinical$study_name <- normalize_study_name(clinical$study_name)

dataset_variable<-read.csv("00.rawData/clinical/dataset_variable_norm.csv",header=T)


dataset_variable$variable <- apply(dataset_variable, 1, function(x) {paste(colnames(dataset_variable)[which(x==1)], collapse = ",")})
dataset_variable<-dataset_variable[,c("study_name_norm","variable")]

datasets<-dataset_variable$study_name
study_num<-length(datasets)

results<-matrix(nrow=1,ncol=6)
colnames(results)<-c("Variable","Datasets","Statistic_SW","P_SW","Statistic_KS","P_KS")

for (studyname in datasets){
  # studyname<-"MetaCardis_2020_a"
  dataset_clinical<-subset(clinical,study_name==studyname,drop=T)
  covariate<-unlist(strsplit(dataset_variable[dataset_variable$study_name==studyname,]$variable,","))

  for (var in covariate){
     # var<-"BMI"
      data<-dataset_clinical[,var]

      if(length(data)>3 & length(data)<5000){
       ## Shapiro-Wilk
        sw_result<-shapiro.test(data)
        sw_sta<-sw_result$statistic
        sw_p<-sw_result$p.value
        }
       else{
        sw_sta<-NA
        sw_p<-NA
        }

       ## Kolmogorov-Smirnov test 
       ks_result<-ks.test(data, "pnorm", mean = mean(data, na.rm = TRUE), sd = sd(data, na.rm = TRUE))
       ks_sta<-ks_result$statistic
       ks_p<-ks_result$p.value
 
       temp<-c(studyname,var,sw_sta,sw_p,ks_sta,ks_p)
       results<-rbind(results,temp)
    }
 }


write.csv(results,"00.rawData/clinical/normality_test_after_norm.csv",row.names=F)


###########################################################################
## the variable that can't be fit norm distribution 

dataset<-subset(clinical,study_name=="KarlssonFH_2013",drop=T)

dataset$age
dataset$age_norm


dataset<-subset(clinical,study_name=="SchirmerM_2016",drop=T)

dataset$age
dataset$age_norm

mean(dataset$age_norm,na.rm=T)
sd(dataset$age_norm,na.rm=T)

dataset<-subset(clinical,study_name=="LiJ_2014",drop=T)

dataset$fasting_insulin
dataset$fasting_insulin_norm





################################################################################
# annotation file 

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/00.rawData/gene_annotation/3_kegg")

kegg<-read.table("user_ko.txt",header=F,sep="\t",stringsAsFactors = FALSE,fill=T)
locas_id<-read.table("locas_id.txt",header=F,sep="\t",stringsAsFactors = FALSE,fill=T)
colnames(locas_id)<-c("locas","genes")

colnames(kegg)<-c("genes","Orthology")
kegg<-subset(kegg,is.na(Orthology)==F,drop=T)

kegg_rename<-merge(locas_id,kegg,by.x="genes",by.y="genes")

kegg_locas<-kegg_rename[,c(2,3)]
write.csv(kegg_locas,"F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/00.rawData/SV/kegg_locas.csv",row.names=F,quote=F)





