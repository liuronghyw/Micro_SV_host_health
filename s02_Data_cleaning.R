### SV data processing
### 2024-09-03
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

source("functions.R")
##install.packages("ggmosaic")
##library("plyr")
##library("dplyr")

library(foreach)
library(doParallel)
library(dplyr)
library(magrittr)
library(car)
library("scales")


#####################################################################
### Read SV files
dsgv<- read.delim("00.rawData/SV/Whole_dsgv_20240801.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")
vsgv<- read.delim("00.rawData/SV/Whole_vsgv_20240801.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")

dsgv_anno<-read.delim("00.rawData/SV/s02.dSVs_anno.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)
vsgv_anno<-read.delim("00.rawData/SV/s03.vSVs_anno.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)

### Read database files
taxa_length <- read.table("00.rawdata/database/Species_genome_size.tsv",
                        sep = "\t", header = T,check.names = F,stringsAsFactors = F)
taxonomy<- read.csv("00.rawdata/database/representatives.genomes.taxonomy.csv",
                   sep = ",", header = T,check.names = F,stringsAsFactors = F)
taxonomy[taxonomy == ""]<-"Unknown"
colnames(taxonomy)[1]<-'X'
ncbi<-read.csv("00.rawData/database/NCBI_accession.txt", sep = "\t",header = T)
tax_relationship<-read.csv("00.rawData/database/progenome1_species_relationship.tsv",sep = "\t",header = F)


############################################################################################
#### 1 Clean SV data
#### Get clean profiles
## overlap with clinical information and abdundance file 

clinical<-unique(read.csv("00.rawData/clinical/clinical_whole.csv",header=T))

rownames(clinical)<-clinical$match_id
SV_abun_s<-read.csv("00.rawData/taxonomic_abundance/All_sample_bracken_frac.csv",header = T)

sample_name<-colnames(SV_abun_s)[-c(1,2)]
colnames(SV_abun_s)<-c("name","taxonomy_id",gsub("._species.bracken.txt_frac","",sample_name))

clin_sv_inter <- intersect(rownames(clinical),unique(rownames(dsgv),rownames(vsgv)))
clin_abun_inter <- intersect(clin_sv_inter,colnames(SV_abun_s)[-c(1,2)])

clinical_sub <- clinical[clin_abun_inter,]
write.table(clinical_sub, "01.cleanData/phen/Clinical_basic_overlap.tsv",sep = '\t')

dsgv<-dsgv[clin_abun_inter,]
vsgv<-vsgv[clin_abun_inter,]

# Change SV names
colnames(dsgv) <- changeSVname(colnames(dsgv))
colnames(vsgv) <- changeSVname(colnames(vsgv))

### the ratio of NAs
dsgv_non_NA_rate<-apply(dsgv[,-1], 2, myfun<-function(x){sum(is.na(x))/length(clin_abun_inter)})
vsgv_non_NA_rate<-apply(vsgv[,-1], 2, myfun<-function(x){sum(is.na(x))/length(clin_abun_inter)})

write.csv(dsgv_non_NA_rate,"01.cleanData/SV/dsgv_non_NA_rate.csv")
write.csv(vsgv_non_NA_rate,"01.cleanData/SV/vsgv_non_NA_rate.csv")

## Outputs
if(!dir.exists("01.cleanData")){dir.create("01.cleanData")}
if(!dir.exists("01.cleanData/SV_all")){dir.create("01.cleanData/SV_all")}

### 2 Get name conversion table
###  Name conversion

sv<-c(colnames(dsgv),colnames(vsgv))

organism<-str_replace_all(sv,"\\:\\d+_\\d+.*","") %>%
  .[!duplicated(.)]

Short_name<- organism %>% 
  str_replace_all('\\[','') %>%
  str_replace_all('\\]', '') %>%
  str_replace_all(' cf\\.','')

Short_name[grep(' sp\\.', organism, invert = F)] <- Short_name[grep(' sp\\.', organism, invert = F)] %>%
  str_replace_all('sp\\..*','sp')

Fst_letter<-Short_name[grep(' sp\\.', organism, invert = T)] %>%
  str_replace_all(' .*','') %>%
  str_sub(start = 1,end = 1)

Spe_name<-Short_name[grep(' sp\\.', organism, invert = T)] %>%
  str_extract_all(' .*') %>%
  str_replace_all('^ ', '') %>%
  str_replace_all(' .*', '')

Short_name[grep(' sp\\.', organism, invert = T)] <-paste(Fst_letter,'.', Spe_name, sep = '')

taxa_name<-data.frame(NCBI_taxonomy_id = taxonomy$X[match(organism,taxonomy$organism)],
                      organism = as.character(organism), 
                      Short_name = as.character(Short_name), stringsAsFactors = F)

taxa_name$Short_name[match('bacterium LF-3',taxa_name$organism)]<-'bacterium LF-3'
taxa_name<-left_join(taxa_name, ncbi, by = "NCBI_taxonomy_id")

if(!dir.exists("01.cleanData/SV_info")){dir.create("01.cleanData/SV_info")}
#write.table(taxa_name, "01.cleanData/SV_info/Species_name.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

#write.csv(taxa_name, "01.cleanData/SV_info/Species_name.csv",row.names=F)
## fjy<-read.csv("01.cleanData/SV_info/article_fjy.csv",header=T)
## overlap<-merge(taxa_name,fjy,by.x="NCBI_taxonomy_id",by.y="NCBI_taxonomy_id")
## write.csv(overlap,"01.cleanData/SV_info/overlap_article_fjy.csv")


### 4.3 Get SV annotation tables
# SV annotation tables
dsgv_info_anno<-data.frame(dsgv_anno,
                           SV_ID=dsgv_anno$SV_id,
                           Taxonomy_Name = taxa_name$organism[match(str_replace_all(dsgv_anno$Taxonomy_id, '\\..*', ''),
                                                                    taxa_name$NCBI_taxonomy_id)],
                           SV_Name = changeSVname(dsgv_anno$SV_id),
                           SV_Name_long=changeSVname_long(dsgv_anno$SV_id),
                           Taxonomy_ID = dsgv_anno$Taxonomy_id,
                           SV_size = calcSVSize(dsgv_anno$SV_id))[,c(10,7,6,8,9,3,11,4,5)]

vsgv_info_anno<-data.frame(vsgv_anno,
                           SV_ID=vsgv_anno$SV_id,
                           Taxonomy_Name = taxa_name$organism[match(str_replace_all(vsgv_anno$Taxonomy_id, '\\..*', ''),
                                                                    taxa_name$NCBI_taxonomy_id)],
                           SV_Name = changeSVname(vsgv_anno$SV_id),
                           SV_Name_long=changeSVname_long(vsgv_anno$SV_id),
                           Taxonomy_ID = vsgv_anno$Taxonomy_id,
                           SV_size = calcSVSize(vsgv_anno$SV_id))[,c(10,7,6,8,9,3,11,4,5)]


write.table(dsgv_info_anno, "01.cleanData/SV_info/dsgv_info_anno.tsv",sep = "\t", quote = F, col.names = T, row.names = F)
write.table(vsgv_info_anno, "01.cleanData/SV_info/vsgv_info_anno.tsv",sep = "\t", quote = F, col.names = T, row.names = F)

write.csv(dsgv_info_anno, "01.cleanData/SV_info/dsgv_info_anno.csv",row.names = F)
write.csv(vsgv_info_anno, "01.cleanData/SV_info/vsgv_info_anno.csv",row.names = F)


###########################################################################################
### 4.4 Get species information table
###  Get SV number per species

species_dsgv_n<-str_replace_all(colnames(dsgv),"\\:\\d+_\\d+.*","") %>%
  table(.) %>%
  as.data.frame(.)
colnames(species_dsgv_n)<-c("Species","Deletion SVs number")

species_vsgv_n<-str_replace_all(colnames(vsgv),"\\:\\d+_\\d+.*","") %>%
  table(.) %>%
  as.data.frame(.)
colnames(species_vsgv_n)<-c("Species","Variable SVs number")

species_sgv_n<-full_join(species_dsgv_n, species_vsgv_n, by = "Species")
species_sgv_n[is.na(species_sgv_n)]<-0

NCBI_taxonomy_id<-species_sgv_n$Species %>%
  match(.,taxonomy$organism) %>%
  taxonomy$X[.]
species_sgv_n<-data.frame(NCBI_taxonomy_id, species_sgv_n)

## Get sample size per species

sv<-cbind(dsgv,vsgv)
sv_name<-c(colnames(dsgv),colnames(vsgv))

sv_infor_sample_n<-str_replace_all(sv_name,"\\:\\d+_\\d+.*","") %>%
  duplicated(.) %>%
  `!`%>%
  sv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(sv_infor_sample_n) <- "Sample_number"
rownames(sv_infor_sample_n) <- rownames(sv_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
sv_infor_sample_n<-data.frame(Species = rownames(sv_infor_sample_n),sv_infor_sample_n)

Taxonomy_name <- match(sv_infor_sample_n$Species,taxa_name$organism) %>%
  taxa_name$Short_name[.]
sample_n<-data.frame(Short_name=Taxonomy_name, sv_infor_sample_n)


### output different intercontinental
clinical <- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")
table(clinical$inter)

Africa<-subset(clinical,inter=="Africa",drop=T)$match_id
Americas<-subset(clinical,inter=="Americas",drop=T)$match_id
Asian<-subset(clinical,inter=="Asian",drop=T)$match_id
Europe<-subset(clinical,inter=="Europe",drop=T)$match_id

Africa_dsgv<- dsgv[rownames(dsgv) %in% Africa,] 
Americas_dsgv<- dsgv[rownames(dsgv) %in% Americas,] 
Asian_dsgv<- dsgv[rownames(dsgv) %in% Asian,] 
Europe_dsgv<- dsgv[rownames(dsgv) %in% Europe,] 

Africa_vsgv<- vsgv[rownames(vsgv) %in% Africa,] 
Americas_vsgv<- vsgv[rownames(vsgv) %in% Americas,] 
Asian_vsgv<- vsgv[rownames(vsgv) %in% Asian,] 
Europe_vsgv<- vsgv[rownames(vsgv) %in% Europe,] 


Africa_sv<-cbind(Africa_dsgv,Africa_vsgv)
Americas_sv<-cbind(Americas_dsgv,Americas_vsgv)
Asian_sv<-cbind(Asian_dsgv,Asian_vsgv)
Europe_sv<-cbind(Europe_dsgv,Europe_vsgv)

Africa_sv_name<-c(colnames(Africa_dsgv),colnames(Africa_vsgv))
Americas_sv_name<-c(colnames(Americas_dsgv),colnames(Americas_vsgv))
Asian_sv_name<-c(colnames(Asian_dsgv),colnames(Asian_vsgv))
Europe_sv_name<-c(colnames(Europe_dsgv),colnames(Europe_vsgv))

## Africa sample size per species
Africa_infor_sample_n<-str_replace_all(Africa_sv_name,"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  Africa_sv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(Africa_infor_sample_n) <- "Africa"
rownames(Africa_infor_sample_n) <- rownames(Africa_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
Africa_infor_sample_n<-data.frame(Species = rownames(Africa_infor_sample_n),Africa_infor_sample_n)


## Americas sample size per species
Americas_infor_sample_n<-str_replace_all(Americas_sv_name,"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  Americas_sv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(Americas_infor_sample_n) <- "Americas"
rownames(Americas_infor_sample_n) <- rownames(Americas_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
Americas_infor_sample_n<-data.frame(Species = rownames(Americas_infor_sample_n),Americas_infor_sample_n)


## Asian sample size per species
Asian_infor_sample_n<-str_replace_all(Asian_sv_name,"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  Asian_sv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(Asian_infor_sample_n) <- "Asian"
rownames(Asian_infor_sample_n) <- rownames(Asian_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
Asian_infor_sample_n<-data.frame(Species = rownames(Asian_infor_sample_n),Asian_infor_sample_n)


## Europe sample size per species
Europe_infor_sample_n<-str_replace_all(Europe_sv_name,"\\:\\d+_\\d+.*","") %>%
 duplicated(.) %>%
  `!`%>%
  Europe_sv[,.] %>%
  is.na(.) %>%
  `!`%>%
  colSums(.) %>%
  as.data.frame(.)

colnames(Europe_infor_sample_n) <- "Europe"
rownames(Europe_infor_sample_n) <- rownames(Europe_infor_sample_n) %>%
  str_replace_all(.,"\\:\\d+_\\d+.*", "")
Europe_infor_sample_n<-data.frame(Species = rownames(Europe_infor_sample_n),Europe_infor_sample_n)

## merge sample size from different dataset 


info_sample<-merge(Africa_infor_sample_n,Americas_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,Asian_infor_sample_n,by.x="Species",by.y="Species")
info_sample<-merge(info_sample,Europe_infor_sample_n,by.x="Species",by.y="Species")

info_sample<-merge(info_sample,sample_n,by.x="Species",by.y="Species")

infor_sample_n <- info_sample[,c(1,6,7,2,3,4,5)]

## Merge sample size and SV number information
species_sample_n<-dplyr::full_join(species_sgv_n,infor_sample_n, by = "Species")
taxa_length$Species<-str_replace_all(taxa_length$Species, '\\..*', '')
species_sample_n$NCBI_taxonomy_id<-as.character(species_sample_n$NCBI_taxonomy_id)
species_sample_n<-dplyr::left_join(species_sample_n, taxa_length, by = c("NCBI_taxonomy_id"="Species"))
species_sample_n<-data.frame(species_sample_n,
                             SVs.number = species_sample_n[,3]+species_sample_n[,4])

## Merge all information
Informative_species_information <- match(species_sample_n$NCBI_taxonomy_id, taxonomy$X)%>%
  taxonomy[.,] %>%
  cbind(.,species_sample_n)

taxa_name_short<-taxa_name[,c(2,4)]

info <- full_join(Informative_species_information[,-11],
                                             taxa_name_short,
                                             by = 'organism')[,c(1:9,10,13,21,11,12,20,15,16,17,18,14,19)]

colnames(info)[c(1,11:21)]<-c("NCBI_taxonomy_id","Short_name","NCBI_bioproject_accession", 
"Deletion_SVs_number", "Variable_SVs_number","SVs_number","Africa_sample_number",
"America_sample_number","Asian_sample_number","Europe_sample_number","Total_samples_number", "Length")

write.table(info, "01.cleanData/SV_info/Informative_species_information.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(info, "01.cleanData/SV_info/Informative_species_information.csv", row.names = F)


##########################################################################
### filter species�� based on abundance and the ratio of sample. sample size is nrow(clinical), in 5% (454)

##info_sub<-subset(info,Total_samples_number>nrow(clinical)*0.05,drop=T)
info_sub<-info


#########################################################################
###  5 Clean species abundance data 

# Read species abundance files
SV_abun_s<-read.csv("00.rawData/taxonomic_abundance/All_sample_bracken_frac.csv",header = T)

sample_name<-colnames(SV_abun_s)[-c(1,2)]
colnames(SV_abun_s)<-c("name","taxonomy_id",gsub("._species.bracken.txt_frac","",sample_name))

# Basic
all_basic <- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

sv_spe_taxid<-info_sub$NCBI_taxonomy_id
sv_spe_taxid[sv_spe_taxid==245018]<-649756
sv_spe_taxid[sv_spe_taxid==1118061]<-2585118

# relative abundance of species detected with SVs
#SV_abun_s[,-c(1:2)]<-apply(SV_abun_s[,-c(1:2)], 2, myfun<-function(x){x/sum(x)})

SV_abun_s <- sv_spe_taxid %>% 
  match(., tax_relationship$V2) %>% 
  tax_relationship$V1[.] %>% 
  match(., SV_abun_s$taxonomy_id) %>% 
  SV_abun_s[., -c(1:2)] %>% 
  t %>%
  as.data.frame

colnames(SV_abun_s)<-info_sub$organism
#rownames(SV_abun_s)<-str_replace_all(rownames(SV_abun_s), ".S.bracken", "")

SV_abun_sub<-SV_abun_s[,-which(colnames(SV_abun_s) %in% c("Alistipes obesi","Anaerostipes hadrus DSM 3319",
"Methanobrevibacter smithii ATCC 35061","Methanosphaera stadtmanae DSM 3091",
"Enterococcus faecium NRRL B-2354"))]

SV_abun_sub<-SV_abun_sub[match(rownames(all_basic),rownames(SV_abun_sub)),]
rownames(SV_abun_sub)<-rownames(all_basic)

spe_abun_mean<-apply(SV_abun_sub,2,mean,na.rm=T)
write.csv(spe_abun_mean,"01.cleanData/spe_abun_mean.csv")

nr<-nrow(SV_abun_sub)
nc<-ncol(SV_abun_sub)

zero_ratio<-c()

for(j in 1:nc){
   count<-0
   for(i in 1:nr){
     if(SV_abun_sub[i,j]==0){count=count+1}
     }
   zero_ratio[j]<-count/nr
  }

zero_ratio<-data.frame(zero_ratio)
zero_ratio$organism<-colnames(SV_abun_sub)
#write.csv(zero_ratio,"01.cleanData/spe_abun_zero_ratio.csv")

zero_ratio<-subset(zero_ratio,zero_ratio<0.4,drop=T)

info_final<-merge(info_sub,zero_ratio,by.x="organism",by.y="organism")[,-22]
head(info_final)

## delect the repliated 
info_final<-subset(info_final,organism!="Oscillibacter sp. KLE 1728" & 
                  organism!="Sutterella wadsworthensis 2_1_59BFAA" & 
                  organism!="Bacteroides massiliensis dnLKV3" &
                  organism!="Lachnospiraceae bacterium 9_1_43BFAA" &
                  organism!="Phascolarctobacterium sp. CAG:266" &
                  organism!="Clostridiales bacterium VE202-14" &
                  organism!="Clostridium sp. KLE 1755" &
                  organism!="Lachnospiraceae bacterium 1_4_56FAA" &
                  organism!="Lachnospiraceae bacterium 3_1_57FAA_CT1",
                  drop=T)

write.table(info_final, "01.cleanData/SV_info/Informative_species_information_final.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(info_final, "01.cleanData/SV_info/Informative_species_information_final.csv", row.names = F)


############################### abundance final

# Read species abundance files
SV_abun_s<-read.csv("00.rawData/taxonomic_abundance/All_sample_bracken_frac.csv",header = T)
sample_name<-colnames(SV_abun_s)[-c(1,2)]
colnames(SV_abun_s)<-c("name","taxonomy_id",gsub("._species.bracken.txt_frac","",sample_name))

# Basic
all_basic <- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

sv_spe_taxid<-info$NCBI_taxonomy_id
sv_spe_taxid[sv_spe_taxid==245018]<-649756
sv_spe_taxid[sv_spe_taxid==1118061]<-2585118

# relative abundance of species detected with SVs
#SV_abun_s[,-c(1:2)]<-apply(SV_abun_s[,-c(1:2)], 2, myfun<-function(x){x/sum(x)})

SV_abun_s <- sv_spe_taxid %>% 
  match(., tax_relationship$V2) %>% 
  tax_relationship$V1[.] %>% 
  match(., SV_abun_s$taxonomy_id) %>% 
  SV_abun_s[., -c(1:2)] %>% 
  t %>%
  as.data.frame

colnames(SV_abun_s)<-info$organism
rownames(SV_abun_s)<-str_replace_all(rownames(SV_abun_s), ".S.bracken", "")

SV_abun_sub<-SV_abun_s[,-which(colnames(SV_abun_s) %in% c("Alistipes obesi","Anaerostipes hadrus DSM 3319",
"Methanobrevibacter smithii ATCC 35061","Methanosphaera stadtmanae DSM 3091",
"Enterococcus faecium NRRL B-2354"))]

SV_abun_sub<-SV_abun_sub[match(rownames(all_basic),rownames(SV_abun_sub)),]
rownames(SV_abun_sub)<-rownames(all_basic)

# outputs
if(!dir.exists("01.cleanData/mbio_all")){dir.create("01.cleanData/mbio_all")}
write.table(SV_abun_sub, "01.cleanData/mbio_all/SV_species_abun.tsv",sep = '\t')

mean(rowSums(SV_abun_sub,na.rm = T),na.rm = T)  # Mean of total relative abundance of species detected with SVs: 0.825154
rowSums(SV_abun_sub,na.rm = T)[rowSums(SV_abun_sub,na.rm = T)!=0] %>% min # Minimum of total relative abundance of species detected with SVs: 
max(rowSums(SV_abun_sub,na.rm = T),na.rm = T)   # Maximum of total relative abundance of species detected with SVs: 


###################################################################################################
###  select dsgv and vsgv information

info_final<-read.delim("01.cleanData/SV_info/Informative_species_information_final.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)
organism<-info_final$organism

select_column<-c()
for(i in 1:ncol(dsgv)){
  name<-str_replace_all(colnames(dsgv[i]),"\\:\\d+_\\d+.*","")
  if(name %in% organism){select_column<-c(select_column,i)}
  }
dsgv_sub<-dsgv[,select_column]

select_column<-c()
for(i in 1:ncol(vsgv)){
  name<-str_replace_all(colnames(vsgv[i]),"\\:\\d+_\\d+.*","")
  if(name %in% organism){select_column<-c(select_column,i)}
  }
vsgv_sub<-vsgv[,select_column]

write.table(dsgv_sub,"01.cleanData/SV_all/dsgv_all.tsv",sep = '\t')
write.table(vsgv_sub,"01.cleanData/SV_all/vsgv_all.tsv",sep = '\t')
save(dsgv_sub, file = "01.cleanData/SV_all/dsgv.RData")
save(vsgv_sub, file = "01.cleanData/SV_all/vsgv.RData")


################################ 4.5 Get distance matrices
#### 4.5.1 All samples

### output different cohort
clinical <- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

study_samples<-as.data.frame(table(clinical$study_name),responseName = "sample_size")
study_name_id<-as.character(study_samples$Var1)
study_num<-nrow(study_samples)

for(j in 1:study_num){
  #j<-1
  samp<-subset(clinical,study_name==study_name_id[j],drop=T)$match_id
  temp_dsgv<- dsgv_sub[rownames(dsgv_sub) %in% samp,] 
  temp_vsgv<- vsgv_sub[rownames(vsgv_sub) %in% samp,] 

  dsgv_tsv<-paste("01.cleanData/SV_all/dsgv_",study_name_id[j],".tsv",sep="")
  vsgv_tsv<-paste("01.cleanData/SV_all/vsgv_",study_name_id[j],".tsv",sep="")

  dsgv_rdata<-paste("01.cleanData/SV_all/dsgv_",study_name_id[j],".RData",sep="")
  vsgv_rdata<-paste("01.cleanData/SV_all/vsgv_",study_name_id[j],".RData",sep="")

  write.table(temp_dsgv,dsgv_tsv,sep = '\t')
  write.table(temp_vsgv,vsgv_tsv,sep = '\t')
  save(temp_dsgv, file = dsgv_rdata)
  save(temp_vsgv, file = vsgv_rdata)

  msv_dist<-paste("01.cleanData/SV_all/distMat/",study_name_id[j],"_msv_dist.RData",sep="")
  msv_std<-paste("01.cleanData/SV_all/distMat/",study_name_id[j],"_msv_dist_std.RData",sep="")

  dataset_msv_dist(temp_vsgv,temp_dsgv,info_final,msv_dist,msv_std)
}


#############################################################################
### all the cohort

#dataset_msv_dist(vsgv_sub,dsgv_sub,info_final,
#"01.cleanData/SV_all/distMat/all_msv_dist.RData",
#"01.cleanData/SV_all/distMat/all_msv_dist_std.RData")


################################ 4.5 Get distance matrices
#### 4.5.1 All samples

## add_matrix
## example
#m1 <- matrix(c(1:9),nrow=3,ncol=3,dimnames=list(c("c1","c2","c3"),c("c1","c2","c3")))
#m2 <- matrix(c(5:13),nrow=3,ncol=3,dimnames=list(c("c1","c4","c3"),c("c1","c4","c3")))
#merge_mat<-average_matrix(m1,m2)

load("01.cleanData/SV_all/dsgv.RData")
load("01.cleanData/SV_all/vsgv.RData")

#test_vsv<-vsgv_sub[seq(1,100),]
#test_dsv<-dsgv_sub[seq(1,100),]
#all_shared_dsv_dis<-shared_sv_dis_jaccard(test_dsv)
#all_shared_vsv_dis<-shared_sv_dis_canberra(test_vsv)

## msv (vsv+dsv) distance
#sgv<-cbind(vsgv_sub, dsgv_sub)

all_shared_vsv_dis<-shared_sv_dis_canberra(vsgv_sub)
all_shared_dsv_dis<-shared_sv_dis_jaccard(dsgv_sub)

all_shared_sv_dis<-average_matrix(all_shared_vsv_dis,all_shared_dsv_dis)


save(all_shared_sv_dis, file = "01.cleanData/SV_all/all_shared_sv_dis.RData")




##start_time <- Sys.time()

##end_time <- Sys.time()
##run_time <- end_time - start_time
##print(run_time)









##########################################################################
## difference between cohorts

### Read SV files
dsgv<- read.delim("00.rawData/SV/France_dsgv.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")
vsgv<- read.delim("00.rawData/SV/France_vsgv.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")

clinical<-read.csv("00.rawData/clinical/clinical_whole_reads.csv",header=T)[,c("id","study")]
rownames(clinical)<-clinical$id

##########################################################
###  for dsgv (chisq test)

dsgv_ID<-colnames(dsgv)
dsgv$id<-rownames(dsgv)

ndsgv<-ncol(dsgv)-1
pvalue<-c()

#i=4

for(i in 1:ndsgv){
   dsgv_ex<-dsgv[,c(ncol(dsgv),i)]
   ID<-dsgv_ex[1,1]
   colnames(dsgv_ex)[2]<-"dsv"

   cor<-merge(clinical,dsgv_ex,by.clinical="id",by.exp="id")
   cor<-subset(cor,is.na(dsv)==F,drop=T)

   mytable <- xtabs(~study+dsv, data=cor)　　　　
   a<-chisq.test(mytable)　　
   pvalue[i]<-a$p.value
  }

pvalue<-as.matrix(pvalue,ndsgv,1)
rownames(pvalue)<-dsgv_ID

write.csv(pvalue,"diff_study_SV/France_dsv_chisq_test.csv")



##############################################################################
## for vsgv (k-s test)
library(rstatix)

vsgv_ID<-colnames(vsgv)
vsgv$id<-rownames(vsgv)

nvsgv<-ncol(vsgv)-1
pvalue<-c()

#i=4
for(i in 1:nvsgv){
   vsgv_ex<-vsgv[,c(ncol(vsgv),i)]
   ID<-vsgv_ex[1,1]
   colnames(vsgv_ex)[2]<-"vsv"

   cor<-merge(clinical,vsgv_ex,by.clinical="id",by.exp="id")
   cor<-subset(cor,is.na(vsv)==F,drop=T)
　　
   a<-kruskal_test(vsv~as.factor(study),data =cor)　　
   pvalue[i]<-a$p
  }

  pvalue<-as.matrix(pvalue,nvsgv,1)
  rownames(pvalue)<-vsgv_ID

write.csv(pvalue,"diff_study_SV/France_vsv_chisq_test.csv")




