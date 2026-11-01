### Microbial GWA: disease and index website tables
### 2026-10-21
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")
rm(list=ls())
options(scipen=200)
library(dplyr)
library(tidyr)
library(data.table)

################################################################################
## 1 Updated species and SV annotations

info<-fread("01.cleanData/SV_info_QC/Informative_species_information_final.tsv",
            data.table=FALSE,colClasses=c(NCBI_taxonomy_id="character"))
info<-info[,c("NCBI_taxonomy_id","Superkingdom","Phylum","Class","Order",
              "Family","Genus","Species","organism","Short_name","NCBI_bioproject_accession")]
vsv_info<-fread("01.cleanData/SV_info_QC/vsgv_info_anno.tsv",data.table=FALSE)
dsv_info<-fread("01.cleanData/SV_info_QC/dsgv_info_anno.tsv",data.table=FALSE)
sv_info<-bind_rows(mutate(vsv_info,sv_type="vsv"),mutate(dsv_info,sv_type="dsv")) %>%
  rename(Product_prokka=product_prokka,Product_kegg=product_kegg) %>%
  mutate(NCBI_taxonomy_id=sub("\\..*","",Taxonomy_ID)) %>%
  left_join(info,by="NCBI_taxonomy_id")
stopifnot(!anyDuplicated(sv_info[,c("SV_Name","sv_type")]),
          !anyDuplicated(info$NCBI_taxonomy_id),!anyNA(sv_info$organism))

################################################################################
## 2 Model inputs and output columns

model_dirs<-c(model1="07.Microbial_GWAS/anno_adjAbun",
              model2="07.Microbial_GWAS/anno_adjAbunPCs")
annotation_cols<-c("NCBI_taxonomy_id","Superkingdom","Phylum","Class","Order","Family",
                   "Genus","Species","organism","Short_name","NCBI_bioproject_accession",
                   "Taxonomy_ID","Taxonomy_Name","SV_ID","SV_Name_long","Segments_num",
                   "Genes","Product_prokka","Product_kegg","Orthology","sv_type")
output_cols<-list(
  disease=c("SV_Name","phenotype","study","OR","SE","p","OR_left","OR_right","N",
            "bonferroni","fdr","I2","het_p","Adjust_Abundance","p_count",annotation_cols,
            "Is_Bonferroni_p_less_than_0.05","Model"),
  index=c("SV_Name","phenotype","study","beta","SE","p","N","I2","het_p",
          "bonferroni","fdr","Adjust_Abundance","p_count",annotation_cols,
          "Is_Bonferroni_p_less_than_0.05","Model"))
dir.create("07.Microbial_GWAS/website",showWarnings=FALSE)
output_files<-c(disease="07.Microbial_GWAS/website/disease_association_website.csv",
                index="07.Microbial_GWAS/website/disease_index_association_website.csv")
first_write<-c(disease=TRUE,index=TRUE)
total_rows<-c(disease=0L,index=0L)

################################################################################
## 3 Cohort results, meta-analysis results and annotations, in model order

for(Model in names(model_dirs)){
  files<-sort(list.files(model_dirs[[Model]],pattern="_sv\\.csv$",full.names=TRUE))
  stopifnot(length(files)>0)
  for(file in files){
    phenotype<-sub("_sv\\.csv$","",basename(file))
    dat<-fread(file,data.table=FALSE,check.names=FALSE)
    stopifnot(all(c("X","type","Adjust_Abundance","pvalue_meta","bonferroni","fdr") %in% names(dat)),
              !anyDuplicated(dat[,c("X","type")]))
    cohort_cols<-grep("\\.(OR|beta|SE|p|OR_left|OR_right|N)$",names(dat),value=TRUE)
    studies<-unique(sub("\\.[^.]+$","",cohort_cols))
    stopifnot(length(studies)>0)
    table_type<-if(any(grepl("\\.OR$",cohort_cols))) "disease" else "index"
    stopifnot(!(table_type=="disease" && any(grepl("\\.beta$",cohort_cols))))
    numeric_cols<-c(cohort_cols,intersect(c("OR_meta","effect_meta","se_meta","pvalue_meta",
                   "cilb","ciub","I2","het_p","p_count","bonferroni","fdr"),names(dat)))
    dat[numeric_cols]<-lapply(dat[numeric_cols],as.numeric)

    long<-dat %>%
      select(SV_Name=X,sv_type=type,Adjust_Abundance,all_of(cohort_cols)) %>%
      pivot_longer(cols=all_of(cohort_cols),names_to=c("study",".value"),
                   names_pattern="^(.*)\\.(OR|beta|SE|p|OR_left|OR_right|N)$")

    if(length(studies)>1){
      meta<-dat %>%
        select(SV_Name=X,sv_type=type,Adjust_Abundance,
               any_of(c("OR_meta","effect_meta","se_meta","pvalue_meta","cilb","ciub",
                        "I2","het_p","p_count","bonferroni","fdr"))) %>%
        rename(any_of(c(OR="OR_meta",beta="effect_meta",SE="se_meta",p="pvalue_meta",
                        OR_left="cilb",OR_right="ciub")))
      meta$study<-"meta_analysis"
      long<-bind_rows(long,meta)
    }else{
      long<-left_join(long,select(dat,SV_Name=X,sv_type=type,I2,het_p,bonferroni,fdr),
                      by=c("SV_Name","sv_type"))
    }

    # Infinite statistics or abs(OR/beta) >= 1000 invalidate statistics; keep annotations.
    numeric_fields<-intersect(c("OR","beta","SE","p","OR_left","OR_right","N",
                                "I2","het_p","bonferroni","fdr","p_count"),names(long))
    effect_fields<-intersect(c("OR","beta"),names(long))
    invalid_effect<-rowSums(abs(as.matrix(long[,effect_fields,drop=FALSE]))>=1000,na.rm=TRUE)>0
    invalid_rows<-invalid_effect | rowSums(is.infinite(as.matrix(long[,numeric_fields,drop=FALSE])))>0
    long[invalid_rows,numeric_fields]<-NA_real_

    # Keep source corrections: multi-cohort rows have no cohort-specific Bonferroni P.
    long$Model<-Model
    long$phenotype<-phenotype
    long<-long %>%
      mutate(Is_Bonferroni_p_less_than_0.05=case_when(
        is.na(bonferroni) ~ NA_character_,
        bonferroni<0.05 ~ "Yes",
        TRUE ~ "No")) %>%
      left_join(sv_info,by=c("SV_Name","sv_type"))
    stopifnot(nrow(long)==nrow(dat)*(length(studies)+(length(studies)>1)),
              !anyNA(long$SV_ID))
    for(column in setdiff(output_cols[[table_type]],names(long))) long[[column]]<-NA_real_
    long<-long[,output_cols[[table_type]]]

    output_tmp<-paste0(output_files[[table_type]],".tmp")
    write.table(long,output_tmp,sep=",",append=!first_write[[table_type]],
                col.names=first_write[[table_type]],row.names=FALSE,na="NA",
                qmethod="double",fileEncoding="UTF-8")
    first_write[[table_type]]<-FALSE
    total_rows[[table_type]]<-total_rows[[table_type]]+nrow(long)
    cat(Model,phenotype,nrow(long),"rows;",sum(invalid_rows),"invalid rows -> NA\n")
  }
}

################################################################################
## 4 Export separate disease (OR) and index (beta) tables

for(table_type in names(output_files)){
  output_tmp<-paste0(output_files[[table_type]],".tmp")
  stopifnot(!first_write[[table_type]],
            file.copy(output_tmp,output_files[[table_type]],overwrite=TRUE))
  unlink(output_tmp)
  cat("WEBSITE_COMPLETE:",table_type,total_rows[[table_type]],"rows ->",
      output_files[[table_type]],"\n")
}
