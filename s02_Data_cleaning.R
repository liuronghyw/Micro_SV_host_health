### SV data processing
### 2026-09-27
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

## Create output directories before the first write.
dir.create("01.cleanData/phen", recursive = TRUE, showWarnings = FALSE)
dir.create("01.cleanData/SV", recursive = TRUE, showWarnings = FALSE)
dir.create("01.cleanData/SV_info", recursive = TRUE, showWarnings = FALSE)
dir.create("01.cleanData/SV_all/distMat", recursive = TRUE, showWarnings = FALSE)
dir.create("01.cleanData/SV_all/study_pc", recursive = TRUE, showWarnings = FALSE)
dir.create("01.cleanData/mbio_all", recursive = TRUE, showWarnings = FALSE)
## End output directory initialization.


# Match the GB18030 clinical input and downstream cohort filenames.
stopifnot(nzchar(Sys.setlocale("LC_CTYPE", "Chinese_China.936")))

source("functions.R", encoding = "UTF-8")
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
dsgv<- read.delim("00.rawData/SV/Whole_dsgv_20260921.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")
vsgv<- read.delim("00.rawData/SV/Whole_vsgv_20260921.csv",
                 header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")

dsgv_anno<-unique(read.delim("00.rawData/SV/s02.dSVs_anno_20260921.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F))
vsgv_anno<-unique(read.delim("00.rawData/SV/s03.vSVs_anno_20260921.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F))

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

clinical<-unique(read.csv("00.rawData/clinical/clinical_whole_add_norm.csv",header=T,fileEncoding="GB18030"))

rownames(clinical)<-clinical$match_id
SV_abun_s<-read.table("00.rawData/taxonomic_abundance/All_sample_bracken_frac.txt", sep = "\t", quote = "",header = T)
SV_abun_s <- SV_abun_s[, setdiff(colnames(SV_abun_s), "taxonomy_lvl"), drop = FALSE]

sample_name<-colnames(SV_abun_s)[-c(1,2)]
colnames(SV_abun_s)<-c("name","taxonomy_id",sub("\\._species\\.bracken\\.txt_frac$", "", sample_name))

clin_sv_inter <- Reduce(
  intersect,
  list(rownames(clinical), rownames(dsgv), rownames(vsgv))
)

clin_abun_inter <- intersect(clin_sv_inter,colnames(SV_abun_s)[-c(1,2)])

##setdiff(clin_sv_inter, colnames(SV_abun_s)[-c(1,2)])
##setdiff(rownames(clinical), rownames(dsgv))

## Use the donor selection recorded by s00; never substitute another visit.
selected_donors <- read.csv("01.cleanData/donor_selection/s00_retained_samples.csv",
                           fileEncoding = "GB18030", stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(selected_donors$match_id))
clin_abun_inter <- intersect(clin_abun_inter, selected_donors$match_id)
clinical_sub <- clinical[clin_abun_inter, , drop = FALSE]
stopifnot(!anyDuplicated(clinical_sub$match_id),
          !any(clinical_sub$study_name == "HallAB_2017"))

write.table(clinical_sub, "01.cleanData/phen/Clinical_basic_overlap.tsv",sep = '\t')

dsgv<-dsgv[clin_abun_inter, , drop = FALSE]
vsgv<-vsgv[clin_abun_inter, , drop = FALSE]

# Change SV names
colnames(dsgv) <- changeSVname(colnames(dsgv))
colnames(vsgv) <- changeSVname(colnames(vsgv))

### the ratio of NAs
dsgv_NA_rate<-if(ncol(dsgv)>0L) colMeans(is.na(dsgv)) else numeric()
vsgv_NA_rate<-if(ncol(vsgv)>0L) colMeans(is.na(vsgv)) else numeric()

### the ratio of non-delection of the dsv
dsgv_non_deletion_rate <- if(ncol(dsgv)>0L) {
  vapply(dsgv, function(x) {
    n <- sum(!is.na(x))
    if(n==0L) NA_real_ else sum(x==1,na.rm=TRUE)/n
  }, numeric(1))
} else numeric()

write.csv(dsgv_non_deletion_rate,"01.cleanData/SV/dsgv_non_deletion_rate.csv")
write.csv(dsgv_NA_rate,"01.cleanData/SV/dsgv_NA_rate.csv")
write.csv(vsgv_NA_rate,"01.cleanData/SV/vsgv_NA_rate.csv")

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
                           SV_size = calcSVSize(dsgv_anno$SV_id))[,c(11,8,13,9,10,3,12,4,5,6,7)]

vsgv_info_anno<-data.frame(vsgv_anno,
                           SV_ID=vsgv_anno$SV_id,
                           Taxonomy_Name = taxa_name$organism[match(str_replace_all(vsgv_anno$Taxonomy_id, '\\..*', ''),
                                                                    taxa_name$NCBI_taxonomy_id)],
                           SV_Name = changeSVname(vsgv_anno$SV_id),
                           SV_Name_long=changeSVname_long(vsgv_anno$SV_id),
                           Taxonomy_ID = vsgv_anno$Taxonomy_id,
                           SV_size = calcSVSize(vsgv_anno$SV_id))[,c(11,8,13,9,10,3,12,4,5,6,7)]


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

## Get sample size per species. Count a sample once when any SV for that
## species is present; this preserves species found in only one SV matrix.
count_species_samples <- function(sv, sv_names, value_name = "Sample_number") {
  if (is.null(sv) || ncol(sv) == 0L || length(sv_names) == 0L) {
    out <- data.frame(Species = character(), value = numeric(), stringsAsFactors = FALSE)
    names(out)[2] <- value_name
    return(out)
  }
  n <- min(ncol(sv), length(sv_names))
  sv <- sv[, seq_len(n), drop = FALSE]
  species <- str_replace_all(sv_names[seq_len(n)], "\\:\\d+_\\d+.*", "")
  groups <- split(seq_len(n), species)
  counts <- vapply(groups, function(idx) {
    present <- !is.na(sv[, idx, drop = FALSE])
    sum(rowSums(present) > 0L)
  }, numeric(1))
  out <- data.frame(Species = names(counts), value = as.numeric(counts), stringsAsFactors = FALSE)
  names(out)[2] <- value_name
  out
}

sv <- cbind(dsgv, vsgv)
sv_name <- c(colnames(dsgv), colnames(vsgv))
sv_infor_sample_n <- count_species_samples(sv, sv_name, "Sample_number")

Taxonomy_name <- match(sv_infor_sample_n$Species,taxa_name$organism) %>%
  taxa_name$Short_name[.]
sample_n<-data.frame(Short_name=Taxonomy_name, sv_infor_sample_n)


### output different intercontinental
clinical <- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")
table(clinical$inter)

Africa<-subset(clinical,inter=="Africa",drop=T)$match_id
Americas<-subset(clinical,inter=="Americas",drop=T)$match_id
Asia<-subset(clinical,inter=="Asia",drop=T)$match_id
Europe<-subset(clinical,inter=="Europe",drop=T)$match_id

merge_sgv_by_sample <- function(dsv, vsv) {
  sample_ids <- union(rownames(dsv), rownames(vsv))
  dsv_aligned <- dsv[match(sample_ids, rownames(dsv)), , drop = FALSE]
  vsv_aligned <- vsv[match(sample_ids, rownames(vsv)), , drop = FALSE]
  rownames(dsv_aligned) <- sample_ids
  rownames(vsv_aligned) <- sample_ids
  cbind(dsv_aligned, vsv_aligned)
}

Africa_dsgv<- dsgv[rownames(dsgv) %in% Africa, , drop = FALSE]
Americas_dsgv<- dsgv[rownames(dsgv) %in% Americas, , drop = FALSE]
Asia_dsgv<- dsgv[rownames(dsgv) %in% Asia, , drop = FALSE]
Europe_dsgv<- dsgv[rownames(dsgv) %in% Europe, , drop = FALSE]

Africa_vsgv<- vsgv[rownames(vsgv) %in% Africa, , drop = FALSE]
Americas_vsgv<- vsgv[rownames(vsgv) %in% Americas, , drop = FALSE]
Asia_vsgv<- vsgv[rownames(vsgv) %in% Asia, , drop = FALSE]
Europe_vsgv<- vsgv[rownames(vsgv) %in% Europe, , drop = FALSE]


Africa_sv <- merge_sgv_by_sample(Africa_dsgv, Africa_vsgv)
Americas_sv <- merge_sgv_by_sample(Americas_dsgv, Americas_vsgv)
Asia_sv <- merge_sgv_by_sample(Asia_dsgv, Asia_vsgv)
Europe_sv <- merge_sgv_by_sample(Europe_dsgv, Europe_vsgv)

Africa_sv_name<-c(colnames(Africa_dsgv),colnames(Africa_vsgv))
Americas_sv_name<-c(colnames(Americas_dsgv),colnames(Americas_vsgv))
Asia_sv_name<-c(colnames(Asia_dsgv),colnames(Asia_vsgv))
Europe_sv_name<-c(colnames(Europe_dsgv),colnames(Europe_vsgv))


## Africa sample size per species
Africa_infor_sample_n <- count_species_samples(Africa_sv, Africa_sv_name, "Africa")


## Americas sample size per species
Americas_infor_sample_n <- count_species_samples(Americas_sv, Americas_sv_name, "Americas")


## Asia sample size per species
Asia_infor_sample_n <- count_species_samples(Asia_sv, Asia_sv_name, "Asia")


## Europe sample size per species
Europe_infor_sample_n <- count_species_samples(Europe_sv, Europe_sv_name, "Europe")

## merge sample size from different dataset 


info_sample <- Reduce(function(x, y) merge(x, y, by = "Species", all = TRUE),
                      list(Africa_infor_sample_n, Americas_infor_sample_n,
                           Asia_infor_sample_n, Europe_infor_sample_n,
                           sample_n))
for (nm in setdiff(names(info_sample), "Species")) {
  info_sample[[nm]][is.na(info_sample[[nm]])] <- 0
}

infor_sample_n <- info_sample[,c(1,6,7,2,3,4,5)]

## Merge sample size and SV number information
species_sample_n<-dplyr::full_join(species_sgv_n,infor_sample_n, by = "Species")
taxa_length$Species<-str_replace_all(taxa_length$Species, '\\..*', '')
species_sample_n$NCBI_taxonomy_id<-as.character(species_sample_n$NCBI_taxonomy_id)
species_sample_n<-dplyr::left_join(species_sample_n, taxa_length, by = c("NCBI_taxonomy_id"="Species"))
species_sample_n<-data.frame(species_sample_n,
                             SVs.number = species_sample_n[,3]+species_sample_n[,4])

## Merge all information using column names rather than fragile positional
## indexes. Missing columns are created explicitly so an empty/partial
## intermediate result does not block the remainder of the workflow.
taxonomy_base <- taxonomy
colnames(taxonomy_base)[1] <- "NCBI_taxonomy_id"
taxonomy_base$NCBI_taxonomy_id <- as.character(taxonomy_base$NCBI_taxonomy_id)
taxonomy_base$Species <- NULL
taxa_name_short <- taxa_name[, c("organism", "NCBI_bioproject_accession"), drop = FALSE]
info <- dplyr::left_join(species_sample_n, taxonomy_base, by = "NCBI_taxonomy_id")
info <- dplyr::left_join(info, taxa_name_short, by = "organism")
rename_map <- c("Deletion.SVs.number" = "Deletion_SVs_number",
                "Variable.SVs.number" = "Variable_SVs_number",
                "SVs.number" = "SVs_number",
                "Africa" = "Africa_sample_number",
                "Americas" = "America_sample_number",
                "Asia" = "Asia_sample_number",
                "Europe" = "Europe_sample_number",
                "Sample_number" = "Total_samples_number",
                "length" = "Length")
for (old in names(rename_map)) {
  if (old %in% names(info)) names(info)[names(info) == old] <- rename_map[[old]]
}
required_info <- c("NCBI_taxonomy_id", "Superkingdom", "Phylum", "Class",
                   "Order", "Family", "Genus", "Species", "organism",
                   "Short_name", "NCBI_bioproject_accession",
                   "Deletion_SVs_number", "Variable_SVs_number", "SVs_number",
                   "Africa_sample_number", "America_sample_number",
                   "Asia_sample_number", "Europe_sample_number",
                   "Total_samples_number", "Length")
for (nm in setdiff(required_info, names(info))) info[[nm]] <- NA
info <- info[, required_info, drop = FALSE]

write.table(info, "01.cleanData/SV_info/Informative_species_information.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(info, "01.cleanData/SV_info/Informative_species_information.csv", row.names = F)


##########################################################################
### filter species based on the ratio of sample.

##info_sub<-subset(info,Total_samples_number>nrow(clinical)*0.05,drop=T)
info_sub<-info


#########################################################################
###  5 Clean species abundance data 

# Read species abundance files
SV_abun_s<-read.table("00.rawData/taxonomic_abundance/All_sample_bracken_frac.txt", sep = "\t", quote = "",header = T)
SV_abun_s <- SV_abun_s[, setdiff(colnames(SV_abun_s), "taxonomy_lvl"), drop = FALSE]

sample_name<-colnames(SV_abun_s)[-c(1,2)]
colnames(SV_abun_s)<-c("name","taxonomy_id",sub("\\._species\\.bracken\\.txt_frac$", "", sample_name))

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

SV_abun_sub<-SV_abun_s[,!colnames(SV_abun_s) %in% c("Enterococcus faecium NRRL B-2354"),drop=F]

SV_abun_sub<-SV_abun_sub[match(rownames(all_basic),rownames(SV_abun_sub)),]
rownames(SV_abun_sub)<-rownames(all_basic)

spe_abun_mean<-apply(SV_abun_sub,2,mean,na.rm=T)
write.csv(spe_abun_mean,"01.cleanData/spe_abun_mean.csv")

nr<-nrow(SV_abun_sub)
nc<-ncol(SV_abun_sub)

zero_ratio<-c()

for(j in seq_len(nc)){
   count<-0
   for(i in seq_len(nr)){
      if(!is.na(SV_abun_sub[i,j]) && SV_abun_sub[i,j]==0){count=count+1}
     }
   zero_ratio[j]<-count/nr
  }

zero_ratio<-data.frame(zero_ratio)
zero_ratio$organism<-colnames(SV_abun_sub)
#write.csv(zero_ratio,"01.cleanData/spe_abun_zero_ratio.csv")
#zero_ratio<-subset(zero_ratio,zero_ratio<0.4,drop=T)
#info_final<-merge(info_sub,zero_ratio,by.x="organism",by.y="organism")[,-22]
#head(info_sub)


## rename the short name
nr<-nrow(info_sub)

for(i in seq_len(nr)){
  if( info_sub$organism[i]=="Clostridiales bacterium VE202-13"){ info_sub$Short_name[i]="C.bacterium_VE202-13"}
   else if ( info_sub$organism[i]=="Clostridiales bacterium VE202-14"){ info_sub$Short_name[i]="C.bacterium_VE202-14"}
   else if ( info_sub$organism[i]=="Clostridiales bacterium VE202-18"){ info_sub$Short_name[i]="C.bacterium_VE202-18"}
   else if ( info_sub$organism[i]=="Clostridium sp. KLE 1755"){ info_sub$Short_name[i]="Clostridium sp.KLE 1755"}
   else if ( info_sub$organism[i]=="Clostridium sp. L2-50"){ info_sub$Short_name[i]="Clostridium sp_2.L2-50"}
   else if ( info_sub$organism[i]=="Lachnospiraceae bacterium 1_4_56FAA"){ info_sub$Short_name[i]="L.bacterium 1_4_56FAA"}
   else if ( info_sub$organism[i]=="Lachnospiraceae bacterium 3_1_46FAA"){ info_sub$Short_name[i]="L.bacterium 3_1_46FAA"}
   else if ( info_sub$organism[i]=="Lachnospiraceae bacterium 3_1_57FAA_CT1"){ info_sub$Short_name[i]="L.bacterium 1_57FAA_CT1"}
   else if ( info_sub$organism[i]=="Lachnospiraceae bacterium 9_1_43BFAA"){ info_sub$Short_name[i]="L.bacterium 9_1_43BFAA"}
   else if ( info_sub$organism[i]=="Oscillibacter sp. ER4"){ info_sub$Short_name[i]="Oscillibacter sp. ER4"}
   else if ( info_sub$organism[i]=="Oscillibacter sp. KLE 1728"){ info_sub$Short_name[i]="Oscillibacter sp. KLE 1728"}
   else if ( info_sub$organism[i]=="Phascolarctobacterium sp. CAG:207"){ info_sub$Short_name[i]="Phascolarctobacterium sp. CAG:207"}
   else if ( info_sub$organism[i]=="Phascolarctobacterium sp. CAG:266"){ info_sub$Short_name[i]="Phascolarctobacterium sp. CAG:266"}
   else if ( info_sub$organism[i]=="Sutterella wadsworthensis 2_1_59BFAA"){ info_sub$Short_name[i]="S.wadsworthensis 2_1_59BFAA"}
   else if ( info_sub$organism[i]=="Sutterella wadsworthensis HGA0223"){ info_sub$Short_name[i]="S.wadsworthensis HGA0223"}
 }

info_final<-info_sub

write.table(info_final, "01.cleanData/SV_info/Informative_species_information_final.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(info_final, "01.cleanData/SV_info/Informative_species_information_final.csv", row.names = F)


############################### abundance final

# Read species abundance files
SV_abun_s<-read.table("00.rawData/taxonomic_abundance/All_sample_bracken_frac.txt", sep = "\t", quote = "",header = T)
SV_abun_s <- SV_abun_s[, setdiff(colnames(SV_abun_s), "taxonomy_lvl"), drop = FALSE]

sample_name<-colnames(SV_abun_s)[-c(1,2)]
colnames(SV_abun_s)<-c("name","taxonomy_id",sub("\\._species\\.bracken\\.txt_frac$", "", sample_name))

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

#write.csv(SV_abun_s,"SV_abun_s.csv")

SV_abun_sub<-SV_abun_s[,!colnames(SV_abun_s) %in% c("Enterococcus faecium NRRL B-2354"),drop=F]

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


# Raw dSV/vSV may overlap; exclusivity is enforced only in the QC workflow.


select_column<-c()
for(i in seq_len(ncol(dsgv))){
  name<-str_replace_all(colnames(dsgv[i]),"\\:\\d+_\\d+.*","")
  if(name %in% organism){select_column<-c(select_column,i)}
  }
dsgv_sub <- dsgv[, select_column, drop = FALSE]

select_column<-c()
for(i in seq_len(ncol(vsgv))){
  name<-str_replace_all(colnames(vsgv[i]),"\\:\\d+_\\d+.*","")
  if(name %in% organism){select_column<-c(select_column,i)}
  }
vsgv_sub <- vsgv[, select_column, drop = FALSE]

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

for(j in seq_len(study_num)){
  #j<-1
  samp<-subset(clinical,study_name==study_name_id[j],drop=T)$match_id
  temp_dsgv <- dsgv_sub[rownames(dsgv_sub) %in% samp, , drop = FALSE]
  temp_vsgv <- vsgv_sub[rownames(vsgv_sub) %in% samp, , drop = FALSE]
  if (length(samp) == 0L || nrow(temp_dsgv) == 0L || nrow(temp_vsgv) == 0L) {
    warning("Skipping study with no overlapping dSV/vSV samples: ", study_name_id[j])
    next
  }

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



#########################################################################################
##  SV distance matrices of all species


workers <- parallel::detectCores(logical = TRUE)
if (is.na(workers)) workers <- 1L
workers <- min(6L, max(1L, workers - 5L))

cache_dir <- "01.cleanData/SV_all/distMat/species_msv_cache_validated_v2"
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

tasks <- lapply(seq_len(nrow(info)), function(i) {
  spe_name <- str_replace_all(info$organism[i], "\\[", "\\\\\\[") %>%
    str_replace_all("\\]", "\\\\\\]")
  list(
    index = i,
    species = as.character(info$organism[i]),
    cache_dir = cache_dir,
    vsgv_i = vsgv[, grep(spe_name, colnames(vsgv)), drop = FALSE],
    dsgv_i = dsgv[, grep(spe_name, colnames(dsgv)), drop = FALSE]
  )
})

# calc_sv_dist is defined in functions.R, loaded at preparation.

if (length(tasks) == 0L) {
  all_msv_dist <- list()
} else {
  workers <- min(workers, length(tasks))
  message(sprintf("Using %d parallel workers", workers))

  cl <- parallel::makeCluster(workers)

  dist_results <- tryCatch(
    parallel::parLapplyLB(cl, tasks, calc_sv_dist, average_matrix),
    finally = parallel::stopCluster(cl)
  )

  all_msv_dist <- dist_results
  rm(tasks, dist_results)
  gc()
}

names(all_msv_dist) <- paste("msv_", info$organism, sep = "")

all_msv_dist_std <- lapply(
  all_msv_dist,
  function(x) {
    if (!is.matrix(x)) return(x)
    max_value <- max(x, na.rm = TRUE)
    if (!is.finite(max_value) || max_value == 0) return(x)
    x / max_value
  }
)

if (!dir.exists("01.cleanData/SV_all/distMat")) {
  dir.create("01.cleanData/SV_all/distMat")
}

save(all_msv_dist, file = "01.cleanData/SV_all/all_msv_dist.RData", compress = "gzip", compression_level = 1)
save(all_msv_dist_std, file = "01.cleanData/SV_all/all_msv_dist_std.RData", compress = "gzip", compression_level = 1)




##################################################################################
## 4 Species-level PCoA

source("functions.R", encoding = "UTF-8")
clinical<- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")
info<-read.delim("01.cleanData/SV_info/Informative_species_information_final.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)

## PC information for each study
study_samples <- as.data.frame(table(clinical$study_name), responseName = "sample_size")
study_name_id <- as.character(study_samples$Var1)
study_msv_pc_cum0.6_info <- NULL
study_pc_dir <- "01.cleanData/SV_all/study_pc"
dir.create(study_pc_dir, recursive = TRUE, showWarnings = FALSE)

for (j in seq_along(study_name_id)) {
  load(paste0("01.cleanData/SV_all/distMat/", study_name_id[j], "_msv_dist_std.RData"))
  samp_j <- clinical$match_id[clinical$study_name == study_name_id[j]]
  study_msv_pc_cum0.6 <- NULL
  study_pc_info_j <- data.frame(matrix(NA, nrow = nrow(info), ncol = 3))

  for (i in seq_len(nrow(info))) {
    # Keep the first n PCs if they jointly explain at least 60%; otherwise keep at most 10 PCs.
    msv_pc_res_i <- get_PCs(msv_dist_std[[i]], eig.cutoff = 0.6, max_pc = 10L)
    if (is.null(msv_pc_res_i)) next

    msv_pc_i <- msv_pc_res_i$PCoA
    pc_names_i <- paste0(info$organism[i], "_PC", seq_len(ncol(msv_pc_i)))
    pc_full_i <- matrix(
      NA_real_,
      nrow = length(samp_j),
      ncol = ncol(msv_pc_i),
      dimnames = list(samp_j, pc_names_i)
    )
    idx_i <- match(rownames(msv_pc_i), samp_j)
    valid_i <- !is.na(idx_i)
    pc_full_i[idx_i[valid_i], ] <- as.matrix(msv_pc_i)[valid_i, , drop = FALSE]
    study_msv_pc_cum0.6 <- cbind(study_msv_pc_cum0.6, pc_full_i)
    study_pc_info_j[i, ] <- c(info$organism[i], msv_pc_res_i$PC_num, msv_pc_res_i$Total_eig)

    study_msv_pc_cum0.6_info <- rbind(
      study_msv_pc_cum0.6_info,
      data.frame(
        Study = study_name_id[j],
        Species = info$organism[i],
        PC_number = msv_pc_res_i$PC_num,
        Explained_variance_proportion = msv_pc_res_i$Total_eig,
        Sample_number = nrow(msv_pc_res_i$PCoA)
      )
    )
  }

  if (!is.null(study_msv_pc_cum0.6)) {
    study_msv_pc_cum0.6 <- as.data.frame(study_msv_pc_cum0.6, check.names = FALSE)
  } else {
    study_msv_pc_cum0.6 <- data.frame(row.names = samp_j)
  }
  colnames(study_pc_info_j) <- c("Species", "PC_number", "Explained_variance_proportion")
  study_file_id <- if (grepl("^[A-Za-z0-9_.-]+$", study_name_id[j])) {
    study_name_id[j]
  } else {
    tolower(gsub(
      "[^A-Za-z0-9]+", "_",
      iconv(study_name_id[j], from = "GB18030", to = "ASCII//TRANSLIT")
    ))
  }
  save(
    study_msv_pc_cum0.6,
    file = file.path(study_pc_dir, paste0(study_file_id, "_msv_pc_cum0.6.RData"))
  )
  write.table(
    study_pc_info_j,
    file = file.path(study_pc_dir, paste0(study_file_id, "_msv_pc_cum0.6_info.tsv")),
    col.names = TRUE,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

write.table(
  study_msv_pc_cum0.6_info,
  "01.cleanData/SV_all/all_study_msv_pc_cum0.6_info.tsv",
  col.names = TRUE,
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)
save(
  study_msv_pc_cum0.6_info,
  file = "01.cleanData/SV_all/study_msv_pc_cum0.6_info.RData"
)



