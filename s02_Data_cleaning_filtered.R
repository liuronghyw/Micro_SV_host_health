### SV data processing
### 2025-09-25
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

## Create output directories before the first write.
dir.create("01.cleanData/phen_QC", recursive = TRUE, showWarnings = FALSE)
dir.create("01.cleanData/SV_QC", recursive = TRUE, showWarnings = FALSE)
dir.create("01.cleanData/SV_info_QC", recursive = TRUE, showWarnings = FALSE)
dir.create("01.cleanData/SV_all_QC/distMat", recursive = TRUE, showWarnings = FALSE)
dir.create("01.cleanData/SV_all_QC/study_pc", recursive = TRUE, showWarnings = FALSE)
dir.create("01.cleanData/mbio_all_QC", recursive = TRUE, showWarnings = FALSE)


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
dsgv<- read.delim("00.rawData/SV/Whole_dsgv_QC_filter.csv",
                  header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")
vsgv<- read.delim("00.rawData/SV/Whole_vsgv_QC_filter.csv",
                  header = T,sep = ",",check.names = F,stringsAsFactors = F,row.names = "")

dsgv_anno<-read.delim("00.rawData/SV/s02.dSVs_anno_20260921_filter.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)
vsgv_anno<-read.delim("00.rawData/SV/s03.vSVs_anno_20260921_filter.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)

# The dSV/vSV tables are mutually exclusive by construction in s01.  Fail early
# if an old/unfiltered file is accidentally supplied.
sv_overlap <- intersect(colnames(dsgv), colnames(vsgv))
if (length(sv_overlap) > 0L) {
  stop("dSV/vSV overlap after s01 filtering: ", length(sv_overlap),
       " SVs (first: ", paste(head(sv_overlap, 5), collapse = ", "), ")")
}

# Match QC annotations to the exact SV IDs retained by s01.
dsgv_anno<-dsgv_anno[dsgv_anno$SV_id %in% colnames(dsgv),,drop=F]
vsgv_anno<-vsgv_anno[vsgv_anno$SV_id %in% colnames(vsgv),,drop=F]
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

clin_sv_inter <- Reduce(intersect, list(rownames(clinical), rownames(dsgv), rownames(vsgv)))
clin_abun_inter <- intersect(clin_sv_inter,colnames(SV_abun_s)[-c(1,2)])

## Use the donor selection recorded by s00; never substitute another visit.
selected_donors <- read.csv("01.cleanData/donor_selection/s00_retained_samples.csv",
                           fileEncoding = "GB18030", stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(selected_donors$match_id))
clin_abun_inter <- intersect(clin_abun_inter, selected_donors$match_id)
clinical_sub <- clinical[clin_abun_inter, , drop = FALSE]
stopifnot(!anyDuplicated(clinical_sub$match_id),
          !any(clinical_sub$study_name == "HallAB_2017"))

write.table(clinical_sub, "01.cleanData/phen_QC/Clinical_basic_overlap.tsv", sep="\t")
write.table(clinical_sub, "01.cleanData/phen_QC/Clinical_basic_overlap_UTF8.tsv",
            sep="\t", fileEncoding="UTF-8")
writeLines(clinical_sub$match_id, "01.cleanData/phen_QC/retained_sample_ids.txt")
write.csv(as.data.frame(table(clinical_sub$study_name)),
          "01.cleanData/phen_QC/sample_counts_by_dataset.csv", row.names=FALSE)


dsgv <- dsgv[clin_abun_inter, , drop = FALSE]
vsgv <- vsgv[clin_abun_inter, , drop = FALSE]

# Change SV names
colnames(dsgv) <- changeSVname(colnames(dsgv))
colnames(vsgv) <- changeSVname(colnames(vsgv))

### the ratio of NAs
dsgv_non_NA_rate <- if (ncol(dsgv) > 0L) colMeans(!is.na(dsgv)) else numeric()
vsgv_non_NA_rate <- if (ncol(vsgv) > 0L) colMeans(!is.na(vsgv)) else numeric()

write.csv(dsgv_non_NA_rate,"01.cleanData/SV_QC/dsgv_non_NA_rate.csv")
write.csv(vsgv_non_NA_rate,"01.cleanData/SV_QC/vsgv_non_NA_rate.csv")

## Outputs
if(!dir.exists("01.cleanData")){dir.create("01.cleanData")}
if(!dir.exists("01.cleanData/SV_all_QC")){dir.create("01.cleanData/SV_all_QC",recursive=TRUE)}

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

if(!dir.exists("01.cleanData/SV_info_QC")){dir.create("01.cleanData/SV_info_QC")}
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


write.table(dsgv_info_anno, "01.cleanData/SV_info_QC/dsgv_info_anno.tsv",sep = "\t", quote = F, col.names = T, row.names = F)
write.table(vsgv_info_anno, "01.cleanData/SV_info_QC/vsgv_info_anno.tsv",sep = "\t", quote = F, col.names = T, row.names = F)

write.csv(dsgv_info_anno, "01.cleanData/SV_info_QC/dsgv_info_anno.csv",row.names = F)
write.csv(vsgv_info_anno, "01.cleanData/SV_info_QC/vsgv_info_anno.csv",row.names = F)



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
clinical <- read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv")
table(clinical$inter)

Africa<-subset(clinical,inter=="Africa",drop=T)$match_id
Americas<-subset(clinical,inter=="Americas",drop=T)$match_id
Asia<-subset(clinical,inter=="Asia",drop=T)$match_id
Europe<-subset(clinical,inter=="Europe",drop=T)$match_id

Africa_dsgv <- dsgv[rownames(dsgv) %in% Africa, , drop = FALSE] 
Americas_dsgv <- dsgv[rownames(dsgv) %in% Americas, , drop = FALSE] 
Asia_dsgv <- dsgv[rownames(dsgv) %in% Asia, , drop = FALSE] 
Europe_dsgv <- dsgv[rownames(dsgv) %in% Europe, , drop = FALSE] 

Africa_vsgv <- vsgv[rownames(vsgv) %in% Africa, , drop = FALSE] 
Americas_vsgv <- vsgv[rownames(vsgv) %in% Americas, , drop = FALSE] 
Asia_vsgv <- vsgv[rownames(vsgv) %in% Asia, , drop = FALSE] 
Europe_vsgv <- vsgv[rownames(vsgv) %in% Europe, , drop = FALSE] 

Africa_dsgv$id<-row.names(Africa_dsgv)
Africa_vsgv$id<-row.names(Africa_vsgv)
Africa_sv<-merge(Africa_dsgv,Africa_vsgv,by.x="id",by.y="id",all.x=T,all.y=T)
Africa_sv <- Africa_sv[, setdiff(colnames(Africa_sv), "id"), drop = FALSE]

Americas_dsgv$id<-row.names(Americas_dsgv)
Americas_vsgv$id<-row.names(Americas_vsgv)
Americas_sv<-merge(Americas_dsgv,Americas_vsgv,by.x="id",by.y="id",all.x=T,all.y=T)
Americas_sv <- Americas_sv[, setdiff(colnames(Americas_sv), "id"), drop = FALSE]

Asia_dsgv$id<-row.names(Asia_dsgv)
Asia_vsgv$id<-row.names(Asia_vsgv)
Asia_sv<-merge(Asia_dsgv,Asia_vsgv,by.x="id",by.y="id",all.x=T,all.y=T)
Asia_sv <- Asia_sv[, setdiff(colnames(Asia_sv), "id"), drop = FALSE]

Europe_dsgv$id<-row.names(Europe_dsgv)
Europe_vsgv$id<-row.names(Europe_vsgv)
Europe_sv<-merge(Europe_dsgv,Europe_vsgv,by.x="id",by.y="id",all.x=T,all.y=T)
Europe_sv <- Europe_sv[, setdiff(colnames(Europe_sv), "id"), drop = FALSE]

## Count samples in which at least one SV for each species is identified.
## Use the column names after merge(); merge() may append .x/.y to names.
calc_inter_sample_n <- function(sv_mat, inter_name) {
  sv_names <- colnames(sv_mat)
  species_names <- str_replace_all(sv_names, "\\:\\d+_\\d+.*", "")
  species_levels <- unique(species_names)

  sample_n <- vapply(species_levels, function(sp) {
    sp_idx <- which(species_names == sp)
    if (length(sp_idx) == 0L) return(0L)
    as.integer(sum(rowSums(!is.na(sv_mat[, sp_idx, drop = FALSE])) > 0L))
  }, integer(1))

  data.frame(
    Species = species_levels,
    setNames(list(sample_n), inter_name),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

## Africa sample size per species
Africa_infor_sample_n <- calc_inter_sample_n(Africa_sv, "Africa")

## Americas sample size per species
Americas_infor_sample_n <- calc_inter_sample_n(Americas_sv, "Americas")

## Asia sample size per species
Asia_infor_sample_n <- calc_inter_sample_n(Asia_sv, "Asia")

## Europe sample size per species
Europe_infor_sample_n <- calc_inter_sample_n(Europe_sv, "Europe")

## merge sample size from different dataset 

info_sample <- merge(Africa_infor_sample_n, Americas_infor_sample_n, by = "Species", all = TRUE)
info_sample <- merge(info_sample, Asia_infor_sample_n, by = "Species", all = TRUE)
info_sample <- merge(info_sample, Europe_infor_sample_n, by = "Species", all = TRUE)

info_sample <- merge(info_sample, sample_n, by = "Species", all = TRUE)
for (nm in setdiff(names(info_sample), "Species")) info_sample[[nm]][is.na(info_sample[[nm]])] <- 0

infor_sample_n <- info_sample[,c(1,6,7,2,3,4,5)]

## Merge sample size and SV number information
species_sample_n<-dplyr::full_join(species_sgv_n,infor_sample_n, by = "Species")
taxa_length$Species<-str_replace_all(taxa_length$Species, '\\..*', '')
species_sample_n$NCBI_taxonomy_id<-as.character(species_sample_n$NCBI_taxonomy_id)
species_sample_n<-dplyr::left_join(species_sample_n, taxa_length, by = c("NCBI_taxonomy_id"="Species"))
species_sample_n<-data.frame(species_sample_n,
                             SVs.number = species_sample_n[,3]+species_sample_n[,4])

## Merge all information using explicit column names
taxonomy_base <- taxonomy
colnames(taxonomy_base)[1] <- "NCBI_taxonomy_id"
taxonomy_base$NCBI_taxonomy_id <- as.character(taxonomy_base$NCBI_taxonomy_id)
taxonomy_base$Species <- NULL
taxa_name_short <- taxa_name[, c("organism", "NCBI_bioproject_accession"), drop = FALSE]
info <- dplyr::left_join(species_sample_n, taxonomy_base, by = "NCBI_taxonomy_id")
info <- dplyr::left_join(info, taxa_name_short, by = "organism")
rename_map <- c("Deletion.SVs.number" = "Deletion_SVs_number", "Variable.SVs.number" = "Variable_SVs_number", "SVs.number" = "SVs_number", "Africa" = "Africa_sample_number", "Americas" = "America_sample_number", "Asia" = "Asia_sample_number", "Europe" = "Europe_sample_number", "Sample_number" = "Total_samples_number","length"="Length")
for (old in names(rename_map)) if (old %in% names(info)) names(info)[names(info) == old] <- rename_map[[old]]
required_info <- c("NCBI_taxonomy_id", "Superkingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "organism", "Short_name", "NCBI_bioproject_accession", "Deletion_SVs_number", "Variable_SVs_number", "SVs_number", "Africa_sample_number", "America_sample_number", "Asia_sample_number", "Europe_sample_number", "Total_samples_number", "Length")
for (nm in setdiff(required_info, names(info))) info[[nm]] <- NA
info <- info[, required_info, drop = FALSE]

write.table(info, "01.cleanData/SV_info_QC/Informative_species_information.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(info, "01.cleanData/SV_info_QC/Informative_species_information.csv", row.names = F)


##########################################################################
### filter species based on abundance and the ratio of sample. sample size is nrow(clinical), in 5% (454)

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
all_basic <- read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv")

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
write.csv(spe_abun_mean,"01.cleanData/mbio_all_QC/spe_abun_mean.csv")

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
#head(info_final)


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

table(info_final$Short_name)

write.table(info_final, "01.cleanData/SV_info_QC/Informative_species_information_final.tsv", col.names = T, row.names = F, sep = "\t", quote = F)
write.csv(info_final, "01.cleanData/SV_info_QC/Informative_species_information_final.csv", row.names = F)



############################### abundance final

# Read species abundance files
SV_abun_s<-read.table("00.rawData/taxonomic_abundance/All_sample_bracken_frac.txt", sep = "\t", quote = "",header = T)
SV_abun_s <- SV_abun_s[, setdiff(colnames(SV_abun_s), "taxonomy_lvl"), drop = FALSE]

sample_name<-colnames(SV_abun_s)[-c(1,2)]
colnames(SV_abun_s)<-c("name","taxonomy_id",sub("\\._species\\.bracken\\.txt_frac$", "", sample_name))

# Basic
all_basic <- read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv")

sv_spe_taxid<-info$NCBI_taxonomy_id
sv_spe_taxid[sv_spe_taxid==245018]<-649756
sv_spe_taxid[sv_spe_taxid==1118061]<-2585118

# relative abundance of species detected with SVs
#SV_abun_s[,-c(1:2)]<-apply(SV_abun_s[,-c(1:2)], 2, myfun<-function(x){x/sum(x)})

SV_abun_s<- sv_spe_taxid %>% 
  match(., tax_relationship$V2) %>% 
  tax_relationship$V1[.] %>% 
  match(., SV_abun_s$taxonomy_id) %>% 
  SV_abun_s[., -c(1:2)] %>% 
  t %>%
  as.data.frame

colnames(SV_abun_s)<-info$organism
rownames(SV_abun_s)<-str_replace_all(rownames(SV_abun_s), ".S.bracken", "")

SV_abun_sub<-SV_abun_s

SV_abun_sub<-SV_abun_sub[match(rownames(all_basic),rownames(SV_abun_sub)),]
rownames(SV_abun_sub)<-rownames(all_basic)

# Remove species columns that are completely NA after matching QC samples.
# These columns contain no abundance information and can cause downstream
# distance calculations to fail (e.g. undefined columns or empty rows).
all_na_cols <- vapply(
  SV_abun_sub,
  function(x) all(is.na(x)),
  logical(1)
)

if (any(all_na_cols)) {
  removed_na_species <- colnames(SV_abun_sub)[all_na_cols]
  message(
    "Removing ", sum(all_na_cols),
    " all-NA species columns from SV_abun_sub: ",
    paste(removed_na_species, collapse = ", ")
  )
  SV_abun_sub <- SV_abun_sub[, !all_na_cols, drop = FALSE]
}

# outputs
if(!dir.exists("01.cleanData/mbio_all_QC")){dir.create("01.cleanData/mbio_all_QC")}
write.table(SV_abun_sub, "01.cleanData/mbio_all_QC/SV_species_abun.tsv",sep = '\t')

mean(rowSums(SV_abun_sub,na.rm = T),na.rm = T)  # Mean of total relative abundance of species detected with SVs: 0.825154
rowSums(SV_abun_sub,na.rm = T)[rowSums(SV_abun_sub,na.rm = T)!=0] %>% min # Minimum of total relative abundance of species detected with SVs: 
max(rowSums(SV_abun_sub,na.rm = T),na.rm = T)   # Maximum of total relative abundance of species detected with SVs: 




###################################################################################################
###  select dsgv and vsgv information

info_final<-read.delim("01.cleanData/SV_info_QC/Informative_species_information_final.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F)
organism<-info_final$organism

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

if(!dir.exists("01.cleanData/SV_all_QC")){dir.create("01.cleanData/SV_all_QC")}

write.table(dsgv_sub,"01.cleanData/SV_all_QC/dsgv_all.tsv",sep = '\t')
write.table(vsgv_sub,"01.cleanData/SV_all_QC/vsgv_all.tsv",sep = '\t')
save(dsgv_sub, file = "01.cleanData/SV_all_QC/dsgv.RData")
save(vsgv_sub, file = "01.cleanData/SV_all_QC/vsgv.RData")


dsgv_select<-colnames(dsgv_sub)
write.csv(dsgv_select,"01.cleanData/SV_info_QC/dsgv_select.csv",row.names=F)

vsgv_select<-colnames(vsgv_sub)
write.csv(vsgv_select,"01.cleanData/SV_info_QC/vsgv_select.csv",row.names=F)



############################################################################
### table S2

### table A

species_all<-unique(read.csv("01.cleanData/SV_info/Informative_species_information_final.csv",header=T))

species_QC<-unique(read.csv("01.cleanData/SV_info_QC/Informative_species_information_final.csv",header=T))
species_QC<-species_QC[,c("organism","Deletion_SVs_number","Variable_SVs_number","SVs_number")]
species_QC$used<-"Yes"

species_all<-merge(species_all,species_QC,by.x="organism",by.y="organism",all.x=T,all.y=T)

write.csv(species_all,"01.cleanData/SV_info_QC/TableS2a.csv",row.names=F)


### table B (dSV)
dsgv_info_anno<-unique(read.csv("01.cleanData/SV_info/dsgv_info_anno.csv",header=T))

dsgv_select<-data.frame(read.csv("01.cleanData/SV_info_QC/dsgv_select.csv",header=T))
colnames(dsgv_select)<-"SV_Name"
dsgv_select$used<-"Yes"

dsgv_all<-merge(dsgv_info_anno,dsgv_select,by.x="SV_Name",by.y="SV_Name",all.x=T,all.y=T)

write.csv(dsgv_all,"01.cleanData/SV_info_QC/TableS2b.csv",row.names=F)



## table C  (vSV)
vsgv_info_anno<-unique(read.csv("01.cleanData/SV_info/vsgv_info_anno.csv",header=T))

vsgv_select<-data.frame(read.csv("01.cleanData/SV_info_QC/vsgv_select.csv",header=T))
colnames(vsgv_select)<-"SV_Name"
vsgv_select$used<-"Yes"

vsgv_all<-merge(vsgv_info_anno,vsgv_select,by.x="SV_Name",by.y="SV_Name",all.x=T,all.y=T)

write.csv(vsgv_all,"01.cleanData/SV_info_QC/TableS2c.csv",row.names=F)





###########################################################################################
### Generate study-level distance matrices and PCs from QC data


load(file = "01.cleanData/SV_all_QC/dsgv.RData")
load(file = "01.cleanData/SV_all_QC/vsgv.RData")


info_qc <- read.delim(
  "01.cleanData/SV_info_QC/Informative_species_information_final.tsv",
  sep = "\t",
  header = TRUE,
  quote = "",
  stringsAsFactors = FALSE
)
clinical_qc <- read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv")

vsgv_pc <- vsgv_sub
dsgv_pc <- dsgv_sub
#vsgv_pc[is.na(vsgv_pc)] <- 0
# Keep missing dSV values as NA; do not coerce them to zero.

dist_dir_qc <- "01.cleanData/SV_all_QC/distMat"
pc_dir_qc <- "01.cleanData/SV_all_QC/study_pc"
dir.create(dist_dir_qc, recursive = TRUE, showWarnings = FALSE)
dir.create(pc_dir_qc, recursive = TRUE, showWarnings = FALSE)

study_samples_qc <- as.data.frame(table(clinical_qc$study_name))
study_name_qc <- as.character(study_samples_qc$Var1)

for (j in seq_along(study_name_qc)) {
  study_j <- study_name_qc[j]
  study_file_id <- if (grepl("^[A-Za-z0-9_.-]+$", study_j)) {
    study_j
  } else {
    tolower(gsub(
      "[^A-Za-z0-9]+", "_",
      iconv(study_j, from = "GB18030", to = "ASCII//TRANSLIT")
    ))
  }
  samp_j <- clinical_qc$match_id[clinical_qc$study_name == study_j]
  samp_j <- intersect(samp_j, rownames(vsgv_pc))
  samp_j <- Reduce(intersect, list(samp_j, rownames(dsgv_pc)))
  if (length(samp_j) == 0L) {
    warning("Skipping study with no overlapping dSV/vSV samples: ", study_j)
    next
  }

  temp_vsgv <- vsgv_pc[samp_j, , drop = FALSE]
  temp_dsgv <- dsgv_pc[samp_j, , drop = FALSE]

  msv_dist <- calc_study_msv_dist(temp_vsgv, temp_dsgv, info_qc)
  msv_dist_std <- lapply(msv_dist, function(x) {
    if (!is.matrix(x)) return(x)
    max_value <- max(x, na.rm = TRUE)
    if (!is.finite(max_value) || max_value == 0) return(x)
    x / max_value
  })
  save(
    msv_dist,
    file = file.path(dist_dir_qc, paste0(study_file_id, "_msv_dist.RData"))
  )
  save(
    msv_dist_std,
    file = file.path(dist_dir_qc, paste0(study_file_id, "_msv_dist_std.RData"))
  )
  study_msv_pc_cum0.6 <- NULL

  for (i in seq_len(nrow(info_qc))) {
    res_i <- get_PCs(msv_dist_std[[i]], eig.cutoff = 0.6, max_pc = 10L)
    if (is.null(res_i)) next

    pc_i <- res_i$PCoA
    if (!is.data.frame(pc_i) && !is.matrix(pc_i)) pc_i <- as.data.frame(pc_i)
    if (nrow(pc_i) == 0L || ncol(pc_i) == 0L) next
    pc_names_i <- paste0(info_qc$organism[i], "_PC", seq_len(ncol(pc_i)))
    pc_full_i <- matrix(
      NA_real_,
      nrow = length(samp_j),
      ncol = ncol(pc_i),
      dimnames = list(samp_j, pc_names_i)
    )
    idx_i <- match(rownames(pc_i), samp_j)
    valid_i <- !is.na(idx_i)
    if (any(valid_i)) {
      pc_i_mat <- as.matrix(pc_i)
      pc_full_i[idx_i[valid_i], ] <- pc_i_mat[valid_i, , drop = FALSE]
    }
    study_msv_pc_cum0.6 <- cbind(study_msv_pc_cum0.6, pc_full_i)
  }

  if (!is.null(study_msv_pc_cum0.6)) {
    study_msv_pc_cum0.6 <- as.data.frame(study_msv_pc_cum0.6, check.names = FALSE)
  } else {
    study_msv_pc_cum0.6 <- data.frame(row.names = samp_j)
  }
  save(
    study_msv_pc_cum0.6,
    file = file.path(pc_dir_qc, paste0(study_file_id, "_msv_pc_cum0.6.RData"))
  )
}


