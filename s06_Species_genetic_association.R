### Species genetic association (all samples for age and BMI, not adjust read count)
### 2026-10-28
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")


## 1 Preparation
### 1.1 Import
##install.packages("tmvnsim")
stopifnot(nzchar(Sys.setlocale("LC_CTYPE", "Chinese_China.936")))
source("functions.R", encoding = "UTF-8")

### 1.2 Inputs
info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)
basic<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")
abun<-read.table("01.cleanData/mbio_all/SV_species_abun.tsv", check.names = F) 

dataset_variable<-read.csv("06.Species_genetic_association/define/dataset_variable.csv", header = T) 
dataset_variable$variable <- apply(dataset_variable, 1, function(x) {paste(colnames(dataset_variable)[which(x==1)], collapse = ",")})
dataset_variable<-dataset_variable[,c("study_name","variable")]
dataset_variable<-dataset_variable[dataset_variable$study_name %in% unique(basic$study_name),,drop=FALSE]

variable_dataset<-read.csv("06.Species_genetic_association/define/variable_dataset.csv", header = T) 
colnames(variable_dataset)[1]<-"variable"
variable_dataset$study_name<- apply(variable_dataset, 1, function(x) {paste(intersect(colnames(variable_dataset)[which(x==1)],dataset_variable$study_name), collapse = ",")})
variable_dataset<-variable_dataset[,c("variable","study_name")]
variable_dataset<-variable_dataset[nzchar(variable_dataset$study_name),,drop=FALSE]


#test<-subset(basic,is.na(total_bilirubin)==F,drop=T)
#table(test$study_name)

######################### 2 Species-level association
### 2.1 Prepare abundance table
if (!dir.exists("06.Species_genetic_association")) {dir.create("06.Species_genetic_association")}
if (!dir.exists("06.Species_genetic_association/RData")) {dir.create("06.Species_genetic_association/RData")}


## IBD covariates are selected within study, before species-specific complete-case filtering.
ibd_studies <- dataset_variable$study_name[vapply(strsplit(dataset_variable$variable, ","),
  function(x) "IBD" %in% x, logical(1))]
ibd_covariate_audit <- do.call(rbind, lapply(ibd_studies, function(study) {
  dat <- basic[basic$study_name==study,,drop=FALSE]
  covariates <- available_species_covariates(dat, c("BMI","age","gender_code"))
  data.frame(study=study, covariates=if(length(covariates)) paste(covariates,collapse=" + ") else "None",
             N_with_IBD=sum(!is.na(dat$IBD)),
             N_complete=sum(complete.cases(dat[,c("IBD",covariates),drop=FALSE])))
}))
write.csv(ibd_covariate_audit, "06.Species_genetic_association/IBD_covariate_audit.csv", row.names=FALSE)

############################################
## calculate the information for all datasets 

datasets<-dataset_variable$study_name
study_num<-length(datasets)

#datasets<-datasets[-seq(1,10)]

for (studyname in datasets){
  set.seed(20261028L + match(studyname, dataset_variable$study_name))
  message("Species associations: ", studyname)
##foreach(j = 1:study_num, .packages = c("stats", "methods", "foreach","magrittr","vegan","microbiome","stringr")) %dopar% {
  sub_adonis_res_BMI<-matrix(nrow=1,ncol=6)
  sub_adonis_res_age<-matrix(nrow=1,ncol=6)
  sub_adonis_res_others<-matrix(nrow=1,ncol=6)

  colnames(sub_adonis_res_BMI)<-c("Species","Prog","R2","P","N","fdr")
  colnames(sub_adonis_res_age)<-c("Species","Prog","R2","P","N","fdr")
  colnames(sub_adonis_res_others)<-c("Species","Prog","R2","P","N","fdr")

  # studyname<-"HallAB_2017"

  sub_basic<-subset(basic,study_name==studyname,drop=T)
  sub_abun<-abun[rownames(sub_basic),]

  dis_dir<-paste("01.cleanData/SV_all/distMat/",studyname,"_msv_dist_std.RData",sep="")
  load(dis_dir)
  sub_msv_dist_std<-msv_dist_std
  stopifnot(identical(names(sub_msv_dist_std), paste0("msv_", info$organism)))
  for (distance in sub_msv_dist_std) validate_species_distance(distance, rownames(sub_basic))

  sub_abun$Others<-1-rowSums(sub_abun)
  sub_abun$Others[sub_abun$Others < 0] <- 0
  sub_abun_clr<-abundances(x=as.data.frame(na.omit(sub_abun)), transform="clr") %>%as.data.frame
  sub_abun_clr <- sub_abun_clr[match(rownames(sub_abun), rownames(sub_abun_clr)),]
  rownames(sub_abun_clr) <- rownames(sub_abun)
  sub_covar<-cbind(sub_basic,sub_abun_clr)

  covariate<-unlist(strsplit(dataset_variable[dataset_variable$study_name==studyname,]$variable,","))

  left_covar <- setdiff(covariate, c("age", "BMI", "age_norm", "BMI_norm", "gender_code", "gender"))
  for (var in covariate){
     if(var=="BMI_norm"){
       sub_BMI<-sub_covar  # Include all disease/control groups; models handle missing values.
       covar_sub<-sub_BMI[,"BMI_norm", drop=FALSE]
       if(nrow(sub_BMI)>=30){
        if("gender" %in% covariate){
          sub_adonis_res_BMI <- my_adonis_terms_noadjAbun(sub_msv_dist_std, covar_sub, sub_BMI,c("gender","gender"),info)$table
         }
         else {
          sub_adonis_res_BMI <- my_adonis_terms_noadjAbun_noadjcov(sub_msv_dist_std, covar_sub,info)$table
         }
        }
       }

     else if(var=="age_norm"){
       sub_age<-sub_covar  # Include all disease/control groups; models handle missing values.
       covar_sub<-sub_age[,"age_norm", drop=FALSE]
       if(nrow(sub_age)>=30){
         if("gender" %in% covariate){
          sub_adonis_res_age <- my_adonis_terms_noadjAbun(sub_msv_dist_std, covar_sub, sub_age,c("gender","gender"),info)$table
          }
         else {
         sub_adonis_res_age <- my_adonis_terms_noadjAbun_noadjcov(sub_msv_dist_std, covar_sub,info)$table
         }
        }
       }

     else if(length(left_covar)>0L && var==left_covar[1L]){
       ## Other phenotypes are fitted once per study, rather than for every variable.
       left_covar<-subset(covariate,covariate!="age" & covariate!="BMI" & covariate!="age_norm" & covariate!="BMI_norm" & covariate!="gender" & covariate!="gender_code")
       ncovar<-length(left_covar)
       covar_sub<-sub_covar[,left_covar,drop=FALSE]
       
       other_traits <- setdiff(left_covar, "IBD")
       sub_adonis_res_others <- data.frame(Species=character(), Prog=character(),
                                          R2=numeric(), P=numeric(), N=numeric(), fdr=numeric())
       if (length(other_traits)) {
         sub_adonis_res_others <- my_adonis_terms_noadjAbun(sub_msv_dist_std,
           sub_covar[,other_traits,drop=FALSE], sub_covar, c("BMI","age","gender"), info)$table
       }
       if ("IBD" %in% left_covar) {
         ibd_covariates <- available_species_covariates(sub_covar, c("BMI","age","gender"))
         ibd_result <- my_adonis_terms_noadjAbun(sub_msv_dist_std,
           sub_covar[,"IBD",drop=FALSE], sub_covar, ibd_covariates, info)$table
         sub_adonis_res_others <- rbind(sub_adonis_res_others, ibd_result)
       }
      }
    sub_adonis_res<-rbind(sub_adonis_res_BMI,sub_adonis_res_age,sub_adonis_res_others)

    out_file<-paste("06.Species_genetic_association/RData/",studyname,"_adonis_res.RData",sep="")
    save(sub_adonis_res, file = out_file)
   }
 }



######################### meta-analysis for different variables 

variable<-variable_dataset$variable

for (covariate in variable){
   # covariate<-"IBD"
   datasets<-unlist(strsplit(variable_dataset[variable_dataset$variable==covariate,]$study_name,","))
   
   study_num<-length(datasets)

   ## the first one 
   if(study_num==1){
     load(paste("06.Species_genetic_association/RData/",datasets[1],"_adonis_res.RData",sep=""))
     tmp_sub_adonis_res<-subset(sub_adonis_res,Prog==covariate,drop=T)
     colnames(tmp_sub_adonis_res)[-c(1,2)]<-c(paste(datasets[1],".",colnames(tmp_sub_adonis_res)[-c(1,2)],sep = ""))
     } 
     else if (study_num>1){
     load(paste("06.Species_genetic_association/RData/",datasets[1],"_adonis_res.RData",sep=""))
     tmp_sub_adonis_res<-subset(sub_adonis_res,Prog==covariate,drop=T)
     colnames(tmp_sub_adonis_res)[-c(1,2)]<-c(paste(datasets[1],".",colnames(tmp_sub_adonis_res)[-c(1,2)],sep = ""))

     for(k in 2:study_num){
       load(paste("06.Species_genetic_association/RData/",datasets[k],"_adonis_res.RData",sep=""))
       tmp_sub_adonis_res_rest<-subset(sub_adonis_res,Prog==covariate,drop=T)
       colnames(tmp_sub_adonis_res_rest)[-c(1,2)]<-c(paste(datasets[k],".",colnames(tmp_sub_adonis_res_rest)[-c(1,2)],sep = ""))
       tmp_sub_adonis_res_rest<-tmp_sub_adonis_res_rest[,-2]
       tmp_sub_adonis_res<-merge(tmp_sub_adonis_res,tmp_sub_adonis_res_rest,by.x="Species",by.y="Species")
     }
    }

   ## my_batch_meta_p   
   covarites_adonis <- my_batch_meta_p(tmp_sub_adonis_res, datasets,seq(5,ncol(tmp_sub_adonis_res),by=4),seq(4,ncol(tmp_sub_adonis_res),by=4))$table
   outfile<-paste("06.Species_genetic_association/",covariate,"_adonis.csv",sep="")
   write.csv(covarites_adonis, outfile,row.names=F)
  }




#############################################################
###  2.3 Association between abundance and prognosis

############################################
## calculate the information for all datasets 

datasets<-dataset_variable$study_name
study_num<-length(datasets)

## Binary variable & logistic regression
binary_var<-c("CRC","IBD","T2D","MetS","CAD","ACVD","cirrhosis","adenoma","PD")

#datasets<-datasets[-seq(1,4)]

for (studyname in datasets){
  set.seed(20261028L + match(studyname, dataset_variable$study_name))
  message("Species associations: ", studyname)
  sub_lm_res_BMI<-matrix(nrow=1,ncol=7)
  sub_lm_res_age<-matrix(nrow=1,ncol=7)
  sub_lg_res_binary<-matrix(nrow=1,ncol=7)
  sub_lm_res_quant<-matrix(nrow=1,ncol=7)

  colnames(sub_lm_res_BMI)<-c("Y","X","Beta","SE","p","N","fdr.p")
  colnames(sub_lm_res_age)<-c("Y","X","Beta","SE","p","N","fdr.p")
  colnames(sub_lg_res_binary)<-c("Y","X","Beta","SE","p","N","fdr.p")
  colnames(sub_lm_res_quant)<-c("Y","X","Beta","SE","p","N","fdr.p")

  # studyname<-"LiJ_2014"

  sub_basic<-subset(basic,study_name==studyname,drop=T)
  sub_abun<-abun[rownames(sub_basic),]

  sub_abun$Others<-1-rowSums(sub_abun)
  sub_abun$Others[sub_abun$Others < 0] <- 0
  sub_abun_clr<-abundances(x=as.data.frame(na.omit(sub_abun)), transform="clr") %>%as.data.frame
  sub_abun_clr <- sub_abun_clr[match(rownames(sub_abun), rownames(sub_abun_clr)),]
  rownames(sub_abun_clr) <- rownames(sub_abun)
  sub_covar<-cbind(sub_basic,sub_abun_clr)

  covariate<-unlist(strsplit(dataset_variable[dataset_variable$study_name==studyname,]$variable,","))

  left_covar <- setdiff(covariate, c("age", "BMI", "age_norm", "BMI_norm", "gender_code", "gender"))
  for (var in covariate){
     if(var=="BMI_norm"){
       covar_sub<-sub_covar[,"BMI_norm", drop=FALSE]
       sub_BMI<-sub_covar  # Include all disease/control groups; models handle missing values.
       covar_BMI_sub<-covar_sub[rownames(sub_BMI),,drop=FALSE]
       sub_abun_BMI_clr<-sub_abun_clr[rownames(sub_BMI),,drop=FALSE]
       if(nrow(sub_BMI)>=30){
        if("gender" %in% covariate){
          sub_lm_res_BMI <- lm_abun_index(covar_BMI_sub,sub_abun_BMI_clr,sub_BMI,c("gender_code","gender_code"))
          }
        else {
         sub_lm_res_BMI <- lm_abun_index_nocov(covar_sub,sub_abun_clr)
         }
        }
       }

     else if(var=="age_norm"){
       covar_sub<-sub_covar[,"age_norm", drop=FALSE]
       sub_age<-sub_covar  # Include all disease/control groups; models handle missing values.
       covar_age_sub<-covar_sub[rownames(sub_age),,drop=FALSE]
       sub_abun_age_clr<-sub_abun_clr[rownames(sub_age),,drop=FALSE]
       if(nrow(sub_age)>=30){
        if("gender" %in% covariate){
          sub_lm_res_age <- lm_abun_index(covar_age_sub,sub_abun_age_clr,sub_age,c("gender_code","gender_code"))
          }
         else {
          sub_lm_res_age <- lm_abun_index_nocov(covar_sub,sub_abun_clr)
          }
         }
        }

    else if(length(left_covar)>0L && var==left_covar[1L]){
       left_covar<-subset(covariate,covariate!="age" & covariate!="BMI" & covariate!="age_norm" & covariate!="BMI_norm" & covariate!="gender_code" & covariate!="gender")
       binary_covar<-subset(left_covar,left_covar %in% binary_var)
       ncovar_binary<-length(binary_covar)
       if(ncovar_binary>0){
          covar_sub_binary<-sub_covar[,binary_covar,drop=FALSE]
          
          other_diseases <- setdiff(binary_covar, "IBD")
          sub_lg_res_binary <- data.frame(Y=character(), X=character(), Beta=numeric(),
                                         SE=numeric(), p=numeric(), N=numeric(), fdr.p=numeric())
          if (length(other_diseases)) {
            sub_lg_res_binary <- lg_abun_index(sub_covar[,other_diseases,drop=FALSE],
              sub_abun_clr, sub_covar, c("BMI","age","gender_code"))
          }
          if ("IBD" %in% binary_covar) {
            ibd_covariates <- available_species_covariates(sub_covar, c("BMI","age","gender_code"))
            ibd_result <- lg_abun_index(sub_covar[,"IBD",drop=FALSE],
              sub_abun_clr, sub_covar, ibd_covariates)
            sub_lg_res_binary <- rbind(sub_lg_res_binary, ibd_result)
          }
       }
       quant_covar<-subset(left_covar,!left_covar %in% binary_var)
       ncovar_quant<-length(quant_covar)
       if(ncovar_quant>0){
          covar_sub_quant<-sub_covar[,quant_covar,drop=FALSE]
          
          sub_lm_res_quant <- lm_abun_index(covar_sub_quant,sub_abun_clr,sub_covar,c("BMI","age","gender_code"))
       }
      }

    lm_res<-rbind(sub_lm_res_BMI,sub_lm_res_age,sub_lg_res_binary,sub_lm_res_quant)
    out_file<-paste("06.Species_genetic_association/RData/",studyname,"_abun_res.RData",sep="")
    save(lm_res, file = out_file)
   }
 }





######################### meta-analysis for different variables 

variable<-variable_dataset$variable

for (covariate in variable){
   # covariate<-"T2D"
   datasets<-unlist(strsplit(variable_dataset[variable_dataset$variable==covariate,]$study_name,","))
   
   study_num<-length(datasets)

   ## the first one 
   if(study_num==1){
     load(paste("06.Species_genetic_association/RData/",datasets[1],"_abun_res.RData",sep=""))
     tmp_sub_lm_res<-subset(lm_res,Y==covariate,drop=T)
     colnames(tmp_sub_lm_res)[-c(1,2)]<-c(paste(datasets[1],".",colnames(tmp_sub_lm_res)[-c(1,2)],sep = ""))
     tmp_sub_lm_res<-tmp_sub_lm_res[,c(2,1,seq(3,ncol(tmp_sub_lm_res)))]

     outfile<-paste("06.Species_genetic_association/",covariate,"_abun.csv",sep="")
     write.csv(tmp_sub_lm_res, outfile,row.names=F)

     } 
     else if (study_num>1){
     load(paste("06.Species_genetic_association/RData/",datasets[1],"_abun_res.RData",sep=""))
     tmp_sub_lm_res<-subset(lm_res,Y==covariate,drop=T)
     colnames(tmp_sub_lm_res)[-c(1,2)]<-c(paste(datasets[1],".",colnames(tmp_sub_lm_res)[-c(1,2)],sep = ""))

     for(k in 2:study_num){
       load(paste("06.Species_genetic_association/RData/",datasets[k],"_abun_res.RData",sep=""))
       tmp_sub_lm_res_rest<-subset(lm_res,Y==covariate,drop=T)
       colnames(tmp_sub_lm_res_rest)[-c(1,2)]<-c(paste(datasets[k],".",colnames(tmp_sub_lm_res_rest)[-c(1,2)],sep = ""))
       tmp_sub_lm_res_rest<-tmp_sub_lm_res_rest[,-1]
       tmp_sub_lm_res<-merge(tmp_sub_lm_res,tmp_sub_lm_res_rest,by.x="X",by.y="X")
     }
     ## my_batch_meta_p   
     covarites_lm <- my_batch_meta_p_lr(tmp_sub_lm_res,seq(3,ncol(tmp_sub_lm_res),by=5),seq(4,ncol(tmp_sub_lm_res),by=5),p_col=seq(5,ncol(tmp_sub_lm_res),by=5),single_study_I2=1)
     outfile<-paste("06.Species_genetic_association/",covariate,"_abun.csv",sep="")
     write.csv(covarites_lm, outfile,row.names=F)
    }
  }





######################################################################################################
### 3 Visualization 

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

### 3.1 Preparation
## Each cell compares species-level SV composition with relative abundance.
## Bonferroni P < 0.05 within each phenotype and analysis; missing results are not negative tests.
info <- read.table("01.cleanData/SV_info_QC/Informative_species_information_final.tsv",
                   sep = "\t", header = TRUE, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(info$Short_name), !anyDuplicated(info$organism))
variable <- as.character(variable_dataset$variable)
stopifnot(!anyDuplicated(variable))
all_adonis.table <- data.frame()
all_abun.table <- data.frame()
bonferroni_column_audit <- data.frame()

for (var in variable) {
  sv_res <- read.csv(paste0("06.Species_genetic_association/", var, "_adonis.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  abun_res <- read.csv(paste0("06.Species_genetic_association/", var, "_abun.csv"),
                       stringsAsFactors = FALSE, check.names = FALSE)
  ## Correct raw P values across all tested species before display filtering.
  sv_p_col <- if ("Meta.p" %in% names(sv_res)) "Meta.p" else
    grep("(^|\\.)P$", names(sv_res), value = TRUE)
  abun_p_col <- if ("pvalue_meta" %in% names(abun_res)) "pvalue_meta" else
    grep("(^|\\.)p$", names(abun_res), value = TRUE)
  if (length(sv_p_col) != 1L || length(abun_p_col) != 1L)
    stop("Cannot identify a unique raw P column for ", var)
  sv_p <- sv_res[[sv_p_col]]
  abun_p <- abun_res[[abun_p_col]]
  for (values in list(sv_p, abun_p)) {
    if (!is.numeric(values) || any(!is.na(values) & (!is.finite(values) | values < 0 | values > 1)))
      stop("Invalid raw P values for ", var)
  }
  sv_bonf <- rep(NA_real_, length(sv_p))
  abun_bonf <- rep(NA_real_, length(abun_p))
  sv_tested <- !is.na(sv_p)
  abun_tested <- !is.na(abun_p)
  sv_bonf[sv_tested] <- p.adjust(sv_p[sv_tested], method = "bonferroni")
  abun_bonf[abun_tested] <- p.adjust(abun_p[abun_tested], method = "bonferroni")
  bonferroni_column_audit <- rbind(bonferroni_column_audit, data.frame(
    phenotype = var, SV_raw_P_column = sv_p_col, SV_tests = sum(sv_tested),
    relative_abundance_raw_P_column = abun_p_col, relative_abundance_tests = sum(abun_tested),
    correction = "bonferroni", adjusted_P_cutoff = 0.05))
  sv_tmp <- data.frame(Species = sv_res$Species, raw.p = sv_p, bonferroni.p = sv_bonf, Prog = var)
  abun_species <- info$Short_name[match(abun_res$X, info$organism)]
  abun_species[abun_res$X %in% info$Short_name] <- abun_res$X[abun_res$X %in% info$Short_name]
  abun_tmp <- data.frame(Species = abun_species, raw.p = abun_p, bonferroni.p = abun_bonf, Prog = var)
  sv_tmp <- sv_tmp[sv_tmp$Species %in% info$Short_name, ]
  abun_tmp <- abun_tmp[!is.na(abun_tmp$Species), ]
  if (anyDuplicated(sv_tmp$Species) || anyDuplicated(abun_tmp$Species))
    stop("Duplicate species-phenotype result for ", var)
  for (values in list(sv_tmp$bonferroni.p, abun_tmp$bonferroni.p)) {
    if (!is.numeric(values) || any(!is.na(values) & (!is.finite(values) | values < 0 | values > 1)))
      stop("Invalid Bonferroni P values for ", var)
  }
  all_adonis.table <- rbind(all_adonis.table, sv_tmp)
  all_abun.table <- rbind(all_abun.table, abun_tmp)
}
write.csv(bonferroni_column_audit, "06.Species_genetic_association/F4A_bonferroni_column_audit.csv", row.names = FALSE)
names(all_abun.table) <- paste0("Abun.", names(all_abun.table))
## Join by species and phenotype, retaining results even if one test is missing.
species.table <- merge(all_adonis.table, all_abun.table,
                       by.x = c("Species", "Prog"), by.y = c("Abun.Species", "Abun.Prog"),
                       all = TRUE, sort = FALSE)
species.table$Abun.Species <- species.table$Species
species.table$Abun.Prog <- species.table$Prog
species.table$sv.MetaSigAssoc <- ifelse(is.na(species.table$bonferroni.p), "Unknown",
                                      ifelse(species.table$bonferroni.p < 0.05, "Yes", "No"))
species.table$abun.MetaSigAssoc <- ifelse(is.na(species.table$Abun.bonferroni.p), "Unknown",
                                        ifelse(species.table$Abun.bonferroni.p < 0.05, "Yes", "No"))
sv_sig <- species.table$sv.MetaSigAssoc == "Yes"
abun_sig <- species.table$abun.MetaSigAssoc == "Yes"
species.sig.table <- species.table[sv_sig | abun_sig, ]
species.sv.sig.table <- species.sig.table[species.sig.table$sv.MetaSigAssoc == "Yes", ]
write.table(species.sig.table, "06.Species_genetic_association/species.sig.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
write.csv(species.sig.table, "06.Species_genetic_association/species.sig.csv", row.names = FALSE)
save(all_adonis.table, file = "06.Species_genetic_association/RData/all_adonis.table.RData")
save(all_abun.table, file = "06.Species_genetic_association/RData/all_abun.table.RData")

### 3.4 Venn diagram
## Count logical sets explicitly; table row indices change when Unknown is absent.
tiff("pics/FS6_species_pvalue_count.venn.tiff", width = 900, height = 900, res = 300)
grid::grid.newpage()
draw.pairwise.venn(sum(sv_sig), sum(abun_sig), sum(sv_sig & abun_sig),
                   category = c("SV composition", "Relative abundance"),
                   lty = "blank", fill = c("#4472c4", "#00b050"), alpha = 0.5,
                   cat.pos = c(0, 0), cat.cex = 0.7, cat.dist = 0.025, scaled = FALSE)
dev.off()

### 3.5 Heatmap: fixed categorical codes, including incomplete test availability
adonis_res.P <- matrix(NA_real_, nrow(info), length(variable),
                       dimnames = list(info$Short_name, variable))
all_abun_res.P <- adonis_res.P
matrix_index <- cbind(match(species.table$Species, rownames(adonis_res.P)),
                      match(species.table$Prog, colnames(adonis_res.P)))
stopifnot(!anyNA(matrix_index))
adonis_res.P[matrix_index] <- species.table$bonferroni.p
all_abun_res.P[matrix_index] <- species.table$Abun.bonferroni.p
sv_available <- !is.na(adonis_res.P)
abun_available <- !is.na(all_abun_res.P)
sv_significant <- sv_available & adonis_res.P < 0.05
abun_significant <- abun_available & all_abun_res.P < 0.05
species.color <- matrix(6L, nrow(adonis_res.P), ncol(adonis_res.P), dimnames = dimnames(adonis_res.P))
species.color[sv_available & abun_available & !sv_significant & !abun_significant] <- 0L
species.color[sv_significant & abun_significant] <- 1L
species.color[sv_available & !sv_significant & abun_significant] <- 2L
species.color[!sv_available & abun_significant] <- 3L
species.color[sv_significant & abun_available & !abun_significant] <- 4L
species.color[sv_significant & !abun_available] <- 5L
species.asso_number <- 1L * sv_significant + 1L * abun_significant
species.asso <- rowSums(species.asso_number)
species_order <- names(sort(species.asso[species.asso > 0], decreasing = TRUE))
## Display all phenotypes with at least one significant association; retain input order.
## Columns with no significant associations are recorded, not silently discarded.
phenotype_order <- c("age_norm", "BMI_norm", "IBD", "CRC", "T2D", "hba1c_norm",
 "fasting_insulin_norm", "CAD", "ACVD", "triglycerides_norm", "hdl_norm",
 "cirrhosis", "PD", "cholesterol_norm", "creatinine_norm", "alt_norm")
phenotype_order <- unique(c(phenotype_order, variable))
species.plot <- phenotype_order[phenotype_order %in% colnames(species.color) &
  phenotype_order %in% names(which(colSums(species.asso_number) > 0))]
if (!length(species_order) || !length(species.plot)) stop("No significant associations to plot")
species.color.plot <- species.color[species_order, species.plot, drop = FALSE]
phenotype_labels <- c(age_norm = "Age", BMI_norm = "BMI", IBD = "IBD", CRC = "CRC",
 T2D = "T2D", hba1c_norm = "HbA1c", fasting_insulin_norm = "Fasting insulin",
 CAD = "CAD", ACVD = "ACVD", triglycerides_norm = "Triglycerides", hdl_norm = "HDL",
 cirrhosis = "Cirrhosis", PD = "PD", cholesterol_norm = "Cholesterol",
 creatinine_norm = "Creatinine", alt_norm = "ALT", MetS = "MetS", adenoma = "Adenoma",
 fasting_glucose_norm = "Fasting glucose", ldl_norm = "LDL", mean_press_norm = "Mean pressure",
 ast_norm = "AST", total_bilirubin_norm = "Total bilirubin")
colnames(species.color.plot) <- unname(phenotype_labels[species.plot])
stopifnot(!anyNA(colnames(species.color.plot)))
category_labels <- c("Neither significant (both tests available)",
 "Both SV and relative abundance significant", "Relative abundance significant; SV not significant",
 "Relative abundance significant; SV unavailable", "SV significant; relative abundance not significant",
 "SV significant; relative abundance unavailable", "No significant result; one or both tests unavailable")
category_colors <- c("#e9e9e9", "#207F4C", "#ECD452", "#EA8958", "#2983BB", "#63BBD0", "#FFFFFF")
category_audit <- species.table
category_audit$category_code <- species.color[matrix_index]
category_audit$category <- category_labels[category_audit$category_code + 1L]
category_audit$displayed <- category_audit$Species %in% species_order & category_audit$Prog %in% species.plot
write.csv(category_audit, "06.Species_genetic_association/F4A_category_audit.csv", row.names = FALSE)
write.csv(data.frame(phenotype = variable, SV_significant = colSums(sv_significant),
                    relative_abundance_significant = colSums(abun_significant),
                    displayed = variable %in% species.plot),
          "06.Species_genetic_association/F4A_phenotype_counts.csv", row.names = FALSE)
write.csv(species.color.plot, "06.Species_genetic_association/F4A_category_matrix.csv")

## Explicit breaks prevent unused categorical codes from changing color mappings.
tiff("pics/F3A_specie_prognosis.heatmap_bonferroni_value.tiff", width = 3600, height = 5200,
     res = 300, compression = "lzw")

heatmap.2(species.color.plot, col = category_colors, breaks = seq(-0.5, 6.5, by = 1),
          scale = "none", trace = "none", Rowv = FALSE, Colv = FALSE, dendrogram = "none",
          density.info = "none", cexCol = 1.5, srtCol = 58, cexRow = 1.5,
          labRow = as.expression(lapply(rownames(species.color.plot), function(a) bquote(italic(.(a))))),
          colsep = seq_len(ncol(species.color.plot) - 1L),
          rowsep = seq_len(nrow(species.color.plot) - 1L),
          sepcolor = "#30303a", sepwidth = c(0.02, 0.02), key = FALSE,
          lmat = rbind(c(4, 4, 3), c(2, 1, 0)), lhei = c(0.05, 3.5), lwid = c(0.05, 8, 3.4),
          key.title = NA, margins = c(16, 8))
par(fig = c(0, 1, 0, 1), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("bottomright", category_labels[c(2:6, 1, 7)], bty = "n",
       title = "Bonferroni-adjusted P < 0.05",
       pch = 22, pt.bg = category_colors[c(2:6, 1, 7)], col = "#606060", cex = 13 / 12, pt.cex = 20 / 12,
       inset = c(0.015, 0.015), xpd = NA)

dev.off()




#######################################################################################################
### 3.6 Species pcoa  (特定的species和预后表型之间关系的展示）

## continues variable 
## F.prausnitzii AsnicarF_2021 age

load("01.cleanData/SV_all/distMat/AsnicarF_2021_msv_dist_std.RData")

clin_AsnicarF_2021<-subset(basic,study_name=="AsnicarF_2021",drop=T)
species_short_name<-"F.prausnitzii"
species_dis<-msv_dist_std[[paste("msv_",info$organism[match(species_short_name,info$Short_name)],sep = "")]]

prog_name<-"age"
prog_vec<-clin_AsnicarF_2021[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

plot_study <- unique(basic[species_prog_inter,"study_name"])
stopifnot(length(plot_study)==1L)
plot_results <- new.env()
load(paste0("06.Species_genetic_association/RData/",plot_study,"_adonis_res.RData"), envir=plot_results)
plot_phenotype <- if(prog_name %in% c("age","BMI")) paste0(prog_name,"_norm") else prog_name
adonis_test <- subset(plot_results$sub_adonis_res, Species==species_short_name & Prog==plot_phenotype)
stopifnot(nrow(adonis_test)==1L)
adonis_test[["Pr(>F)"]] <- adonis_test$P

exp <- substitute(
  atop(
    atop(Age (AsnicarF_2021), italic(F.prausnitzii)),
    atop(adonis~R ^ 2 == r2_value ~italic(P) == p_value, N == n_value)
  ),
  list(r2_value = signif(adonis_test$R2[1L], 3), p_value = adonis_test[1L, "Pr(>F)"], n_value = adonis_test$N[1L])
)

Age<-prog_vec_input

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color =Age))+
  geom_point(size = 2,alpha = 0.8)+
  ggtitle(exp)+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  scale_color_distiller(palette = "Spectral")+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        title = element_text(colour = "black",size=12), 
        legend.position = 'bottom',
        #legend.title = element_blank(),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

tiff(file = "pics/F3C_AsnicarF_2021_age_F.prausnitzii.tiff", width =1000, height =1300, res =300) 
print(p_spe_pcoa)
dev.off()



###### Oscillibacter sp. ER4


load("01.cleanData/SV_all/distMat/AsnicarF_2021_msv_dist_std.RData")

clin_AsnicarF_2021<-subset(basic,study_name=="AsnicarF_2021",drop=T)
species_short_name<-"Oscillibacter sp. ER4"
species_dis<-msv_dist_std[[paste("msv_",info$organism[match(species_short_name,info$Short_name)],sep = "")]]

prog_name<-"age"
prog_vec<-clin_AsnicarF_2021[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

plot_study <- unique(basic[species_prog_inter,"study_name"])
stopifnot(length(plot_study)==1L)
plot_results <- new.env()
load(paste0("06.Species_genetic_association/RData/",plot_study,"_adonis_res.RData"), envir=plot_results)
plot_phenotype <- if(prog_name %in% c("age","BMI")) paste0(prog_name,"_norm") else prog_name
adonis_test <- subset(plot_results$sub_adonis_res, Species==species_short_name & Prog==plot_phenotype)
stopifnot(nrow(adonis_test)==1L)
adonis_test[["Pr(>F)"]] <- adonis_test$P

exp <- substitute(
  atop(
    atop(Age (AsnicarF_2021), italic(D.formicigenerans)),
    atop(adonis~R ^ 2 == r2_value ~italic(P) == p_value, N == n_value)
  ),
  list(r2_value = signif(adonis_test$R2[1L], 3), p_value = adonis_test[1L, "Pr(>F)"], n_value = adonis_test$N[1L])
)

Age<-prog_vec_input

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color =Age))+
  geom_point(size = 2,alpha = 0.8)+
  ggtitle(exp)+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  scale_color_distiller(palette = "Spectral")+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        title = element_text(colour = "black",size=12), 
        legend.position = 'bottom',
        #legend.title = element_blank(),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

tiff(file = "pics/F3C_AsnicarF_2021_age_Oscillibacter sp. ER4.tiff", width =1000, height =1300, res =300) 
print(p_spe_pcoa)
dev.off()



#########################################################################################
## F.prausnitzii AsnicarF_2021 BMI

load("01.cleanData/SV_all/distMat/AsnicarF_2021_msv_dist_std.RData")

clin_AsnicarF_2021<-subset(basic,study_name=="AsnicarF_2021",drop=T)
species_short_name<-"F.prausnitzii"
species_dis<-msv_dist_std[[paste("msv_",info$organism[match(species_short_name,info$Short_name)],sep = "")]]

prog_name<-"BMI"
prog_vec<-clin_AsnicarF_2021[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

plot_study <- unique(basic[species_prog_inter,"study_name"])
stopifnot(length(plot_study)==1L)
plot_results <- new.env()
load(paste0("06.Species_genetic_association/RData/",plot_study,"_adonis_res.RData"), envir=plot_results)
plot_phenotype <- if(prog_name %in% c("age","BMI")) paste0(prog_name,"_norm") else prog_name
adonis_test <- subset(plot_results$sub_adonis_res, Species==species_short_name & Prog==plot_phenotype)
stopifnot(nrow(adonis_test)==1L)
adonis_test[["Pr(>F)"]] <- adonis_test$P

exp <- substitute(
  atop(
    atop(BMI (AsnicarF_2021), italic(F.prausnitzii)),
    atop(adonis~R ^ 2 == r2_value ~italic(P) == p_value, N == n_value)
  ),
  list(r2_value = signif(adonis_test$R2[1L], 3), p_value = adonis_test[1L, "Pr(>F)"], n_value = adonis_test$N[1L])
)

BMI<-prog_vec_input

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color =BMI))+
  geom_point(size = 2,alpha = 0.8)+
  ggtitle(exp)+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  scale_color_distiller(palette = "Spectral")+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        title = element_text(colour = "black",size=12), 
        legend.position = 'bottom',
        #legend.title = element_blank(),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

tiff(file = "pics/F3C_AsnicarF_2021_BMI_F.prausnitzii.tiff", width =1000, height =1300, res =300) 
print(p_spe_pcoa)
dev.off()





#########################################################################################
## Oscillibacter sp. ER4 AsnicarF_2021 BMI

load("01.cleanData/SV_all/distMat/AsnicarF_2021_msv_dist_std.RData")

clin_AsnicarF_2021<-subset(basic,study_name=="AsnicarF_2021",drop=T)
species_short_name<-"Oscillibacter sp. ER4"
species_dis<-msv_dist_std[[paste("msv_",info$organism[match(species_short_name,info$Short_name)],sep = "")]]

prog_name<-"BMI"
prog_vec<-clin_AsnicarF_2021[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

plot_study <- unique(basic[species_prog_inter,"study_name"])
stopifnot(length(plot_study)==1L)
plot_results <- new.env()
load(paste0("06.Species_genetic_association/RData/",plot_study,"_adonis_res.RData"), envir=plot_results)
plot_phenotype <- if(prog_name %in% c("age","BMI")) paste0(prog_name,"_norm") else prog_name
adonis_test <- subset(plot_results$sub_adonis_res, Species==species_short_name & Prog==plot_phenotype)
stopifnot(nrow(adonis_test)==1L)
adonis_test[["Pr(>F)"]] <- adonis_test$P


exp <- substitute(
  atop(
    atop(BMI (AsnicarF_2021), italic(Oscillibacter_sp.ER4)),
    atop(adonis~R ^ 2 == r2_value ~italic(P) == p_value, N == n_value)
  ),
  list(r2_value = signif(adonis_test$R2[1L], 3), p_value = adonis_test[1L, "Pr(>F)"], n_value = adonis_test$N[1L])
)


BMI<-prog_vec_input

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color =BMI))+
  geom_point(size = 2,alpha = 0.8)+
  ggtitle(exp)+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  scale_color_distiller(palette = "Spectral")+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        title = element_text(colour = "black",size=12), 
        legend.position = 'bottom',
        #legend.title = element_blank(),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

tiff(file = "pics/F3C_AsnicarF_2021_BMI_Oscillibacter sp. ER4.tiff", width =1000, height =1300, res =300) 
print(p_spe_pcoa)
dev.off()




##########################################################################################
##  bianary variable 

###########################################################################################
## F.prausnitzii IBD 


load("01.cleanData/SV_all/distMat/FranzosaEA_2018_msv_dist_std.RData")

clin_FranzosaEA_2018<-subset(basic,study_name=="FranzosaEA_2018",drop=T)
species_short_name<-"F.prausnitzii"
species_dis<-msv_dist_std[[paste("msv_",info$organism[match(species_short_name,info$Short_name)],sep = "")]]

prog_name<-"IBD"
prog_vec<-clin_FranzosaEA_2018[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

prog_vec_input[prog_vec_input==0]<-"Healthy"
prog_vec_input[prog_vec_input==1]<-"IBD"

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

plot_study <- unique(basic[species_prog_inter,"study_name"])
stopifnot(length(plot_study)==1L)
plot_results <- new.env()
load(paste0("06.Species_genetic_association/RData/",plot_study,"_adonis_res.RData"), envir=plot_results)
plot_phenotype <- if(prog_name %in% c("age","BMI")) paste0(prog_name,"_norm") else prog_name
adonis_test <- subset(plot_results$sub_adonis_res, Species==species_short_name & Prog==plot_phenotype)
stopifnot(nrow(adonis_test)==1L)
adonis_test[["Pr(>F)"]] <- adonis_test$P

exp <- substitute(
  atop(
    atop("IBD" (FranzosaEA_2018), italic(F.prausnitzii)),
    atop(adonis~R ^ 2 == r2_value ~italic(P) == p_value, N == n_value)
  ),
  list(r2_value = signif(adonis_test$R2[1L], 3), p_value = adonis_test[1L, "Pr(>F)"], n_value = adonis_test$N[1L])
)

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color = as.character(prog_vec_input)))+
  geom_point(size = 2,alpha = 0.8)+ ggtitle(exp)+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  stat_ellipse(aes(group = prog_vec_input, fill = prog_vec_input, 
  color = prog_vec_input) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  scale_color_manual(values=c("#6F99ADFF","#BC3C29FF"))+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        plot.title=element_text(size = 14,family ="sans",color="black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text =  element_text(size = 10,family ="sans",color="black"),
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

#pdf("06.Species_genetic_association/R.gnavus_CA_p.genetics.pdf",width = 4, height = 3)

tiff(file = "pics/F3C_F.prausnitzii_FranzosaEA_2018_IBD.tiff", width =1800, height =2500, res =600) 
print(p_spe_pcoa)
dev.off()



############################################################################################
## E.rectale IBD 


load("01.cleanData/SV_all/distMat/HallAB_2017_msv_dist_std.RData")

clin_HallAB_2017<-subset(basic,study_name=="HallAB_2017",drop=T)
species_short_name<-"E.rectale"
species_dis<-msv_dist_std[[paste("msv_",info$organism[match(species_short_name,info$Short_name)],sep = "")]]

prog_name<-"IBD"
prog_vec<-clin_HallAB_2017[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

prog_vec_input[prog_vec_input==0]<-"Healthy"
prog_vec_input[prog_vec_input==1]<-"IBD"

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

plot_study <- unique(basic[species_prog_inter,"study_name"])
stopifnot(length(plot_study)==1L)
plot_results <- new.env()
load(paste0("06.Species_genetic_association/RData/",plot_study,"_adonis_res.RData"), envir=plot_results)
plot_phenotype <- if(prog_name %in% c("age","BMI")) paste0(prog_name,"_norm") else prog_name
adonis_test <- subset(plot_results$sub_adonis_res, Species==species_short_name & Prog==plot_phenotype)
stopifnot(nrow(adonis_test)==1L)
adonis_test[["Pr(>F)"]] <- adonis_test$P

exp <- substitute(
  atop(
    atop("IBD" (HallAB_2017), italic(E.rectale)),
    atop(adonis~R ^ 2 == r2_value ~italic(P) == p_value, N == n_value)
  ),
  list(r2_value = signif(adonis_test$R2[1L], 3), p_value = adonis_test[1L, "Pr(>F)"], n_value = adonis_test$N[1L])
)

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color = as.character(prog_vec_input)))+
  geom_point(size = 2,alpha = 0.8)+ ggtitle(exp)+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  stat_ellipse(aes(group = prog_vec_input, fill = prog_vec_input, 
  color = prog_vec_input) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  scale_color_manual(values=c("#BC3C29FF","#6F99ADFF"))+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        plot.title=element_text(size = 14,family ="sans",color="black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text =  element_text(size = 10,family ="sans",color="black"),
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

#pdf("06.Species_genetic_association/R.gnavus_CA_p.genetics.pdf",width = 4, height = 3)

tiff(file = "pics/F3C_E.rectale_HallAB_2017_IBD.tiff", width =1800, height =2500, res =600) 
print(p_spe_pcoa)
dev.off()



#################
## CRC ZellerG_2014 E.rectale


load("01.cleanData/SV_all/distMat/ZellerG_2014_msv_dist_std.RData")

clin_ZellerG_2014<-subset(basic,study_name=="ZellerG_2014",drop=T)
species_short_name<-"E.rectale"
species_dis<-msv_dist_std[[paste("msv_",info$organism[match(species_short_name,info$Short_name)],sep = "")]]

prog_name<-"CRC"
prog_vec<-clin_ZellerG_2014[prog_name] %>% na.omit

species_prog_inter<-intersect(rownames(species_dis),rownames(prog_vec))
species_dis_input<-species_dis[match(species_prog_inter,rownames(species_dis)),
                               match(species_prog_inter, colnames(species_dis))]
prog_vec_input<-prog_vec[match(species_prog_inter, rownames(prog_vec)),]

prog_vec_input[prog_vec_input==0]<-"Healthy"
prog_vec_input[prog_vec_input==1]<-"CRC"

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

varExp <- (eigenvals(species_dis_input_mds)/sum(eigenvals(species_dis_input_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

species_dis_input_mds<-cmdscale(species_dis_input, k=5, eig = T)
species_dis_input.pcoa <- data.frame(species_dis_input_mds$points)

plot_study <- unique(basic[species_prog_inter,"study_name"])
stopifnot(length(plot_study)==1L)
plot_results <- new.env()
load(paste0("06.Species_genetic_association/RData/",plot_study,"_adonis_res.RData"), envir=plot_results)
plot_phenotype <- if(prog_name %in% c("age","BMI")) paste0(prog_name,"_norm") else prog_name
adonis_test <- subset(plot_results$sub_adonis_res, Species==species_short_name & Prog==plot_phenotype)
stopifnot(nrow(adonis_test)==1L)
adonis_test[["Pr(>F)"]] <- adonis_test$P

exp <- substitute(
  atop(
    atop("CRC" (ZellerG_2014), italic(E.rectale)),
    atop(adonis~R ^ 2 == r2_value ~italic(P) == p_value, N == n_value)
  ),
  list(r2_value = signif(adonis_test$R2[1L], 3), p_value = adonis_test[1L, "Pr(>F)"], n_value = adonis_test$N[1L])
)

p_spe_pcoa<-ggplot(species_dis_input.pcoa,aes(X1,X2, color = as.character(prog_vec_input)))+
  geom_point(size = 2,alpha = 0.8)+ ggtitle(exp)+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  stat_ellipse(aes(group = prog_vec_input, fill = prog_vec_input, 
  color = prog_vec_input) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  scale_color_manual(values=c("#BC3C29FF","#6F99ADFF"))+
  guides(col = guide_legend(nrow = 2, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        plot.title=element_text(size = 14,family ="sans",color="black"),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text =  element_text(size = 10,family ="sans",color="black"),
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

#pdf("06.Species_genetic_association/R.gnavus_CA_p.genetics.pdf",width = 4, height = 3)

tiff(file = "pics/F3C_E.rectale_ZellerG_2014_CRC.tiff", width =1800, height =2500, res =600) 
print(p_spe_pcoa)
dev.off()



############################################################
## output for Table 3  

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

variable<-variable_dataset$variable

### 1.2 Inputs
info_QC<- read.table("01.cleanData/SV_info_QC/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)


# Derive cohort groups from the same configuration used for the models.
cohort_count <- lengths(strsplit(variable_dataset$study_name, ",", fixed = TRUE))
variable_more <- variable[cohort_count > 1L]
variable_signal <- variable[cohort_count == 1L]
dir.create("06.Species_genetic_association/tableS3", recursive = TRUE, showWarnings = FALSE)
for (phenotype in variable) {
  studies <- strsplit(variable_dataset$study_name[match(phenotype, variable_dataset$variable)], ",", fixed = TRUE)[[1L]]
  sv <- read.csv(paste0("06.Species_genetic_association/", phenotype, "_adonis.csv"), check.names = FALSE)
  abundance <- read.csv(paste0("06.Species_genetic_association/", phenotype, "_abun.csv"), check.names = FALSE)
  result <- build_species_tableS3(sv, abundance, info_QC, studies)
  write.csv(result, paste0("06.Species_genetic_association/tableS3/", phenotype, ".csv"), row.names = FALSE)
}



setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/06.Species_genetic_association")

species.sig<-read.csv("species.sig.csv",header=T)
abun<-subset(species.sig,abun.MetaSigAssoc=="Yes",drop=T)

nrow(abun)
length(unique(abun$Species))

##write.csv(abun,"abun_test.csv")


sv<-subset(species.sig,sv.MetaSigAssoc=="Yes",drop=T)
nrow(sv)
length(unique(sv$Species))

##write.csv(sv,"SV_test.csv")



