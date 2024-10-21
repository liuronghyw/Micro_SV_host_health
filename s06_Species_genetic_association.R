### Species genetic association (healthy age and BMI, not adjust read count)
### 2024-09-4
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

## 1 Preparation
### 1.1 Import
##install.packages("tmvnsim")
source("functions.R")

### 1.2 Inputs
info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)
basic<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")
abun<-read.table("01.cleanData/mbio_all/SV_species_abun.tsv", check.names = F) 

dataset_variable<-read.csv("06.Species_genetic_association/define/dataset_variable.csv", header = T) 
dataset_variable$variable <- apply(dataset_variable, 1, function(x) {paste(colnames(dataset_variable)[which(x==1)], collapse = ",")})
dataset_variable<-dataset_variable[,c("study_name","variable")]

variable_dataset<-read.csv("06.Species_genetic_association/define/variable_dataset.csv", header = T) 
variable_dataset$study_name<- apply(variable_dataset, 1, function(x) {paste(colnames(variable_dataset)[which(x==1)], collapse = ",")})
variable_dataset<-variable_dataset[,c("variable","study_name")]


#test<-subset(basic,is.na(total_bilirubin)==F,drop=T)
#table(test$study_name)

######################### 2 Species-level association
### 2.1 Prepare abundance table
if (!dir.exists("06.Species_genetic_association")) {dir.create("06.Species_genetic_association")}
if (!dir.exists("06.Species_genetic_association/RData")) {dir.create("06.Species_genetic_association/RData")}


############################################
## calculate the information for all datasets 

datasets<-dataset_variable$study_name
study_num<-length(datasets)

#datasets<-datasets[-seq(1,10)]

for (studyname in datasets){
##foreach(j = 1:study_num, .packages = c("stats", "methods", "foreach","magrittr","vegan","microbiome","stringr")) %dopar% {
  sub_adonis_res_BMI<-matrix(nrow=1,ncol=6)
  sub_adonis_res_age<-matrix(nrow=1,ncol=6)
  sub_adonis_res_others<-matrix(nrow=1,ncol=6)

  colnames(sub_adonis_res_BMI)<-c("Species","Prog","R2","P","N","fdr")
  colnames(sub_adonis_res_age)<-c("Species","Prog","R2","P","N","fdr")
  colnames(sub_adonis_res_others)<-c("Species","Prog","R2","P","N","fdr")

  # studyname<-"MetaCardis_2020_a"

  sub_basic<-subset(basic,study_name==studyname,drop=T)
  sub_abun<-abun[rownames(sub_basic),]

  dis_dir<-paste("01.cleanData/SV_all/distMat/",studyname,"_msv_dist_std.RData",sep="")
  load(dis_dir)
  sub_msv_dist_std<-msv_dist_std

  sub_abun$Others<-1-rowSums(sub_abun)
  sub_abun_clr<-abundances(x=as.data.frame(na.omit(sub_abun)), transform="clr") %>%as.data.frame
  sub_abun_clr <- sub_abun_clr[match(rownames(sub_abun), rownames(sub_abun_clr)),]
  rownames(sub_abun_clr) <- rownames(sub_abun)
  sub_covar<-cbind(sub_basic,sub_abun_clr)

  covariate<-unlist(strsplit(dataset_variable[dataset_variable$study_name==studyname,]$variable,","))

  for (var in covariate){
     if(var=="BMI"){
       covar_sub<-sub_covar[,c("BMI","BMI")]
       sub_BMI<-subset(sub_covar,disease=="healthy",drop=T)
       if(nrow(sub_BMI)>=30){
        if("gender" %in% covariate){
          sub_adonis_res_BMI <- my_adonis_terms_noadjAbun(sub_msv_dist_std, covar_sub, sub_BMI,c("gender","gender"),info)$table
         }
         else {
          sub_adonis_res_BMI <- my_adonis_terms_noadjAbun_noadjcov(sub_msv_dist_std, covar_sub,info)$table
         }
        }
       }

     if(var=="age"){
       covar_sub<-sub_covar[,c("age","age")]
       sub_age<-subset(sub_covar,disease=="healthy",drop=T)
       if(nrow(sub_age)>=30){
         if("gender" %in% covariate){
          sub_adonis_res_age <- my_adonis_terms_noadjAbun(sub_msv_dist_std, covar_sub, sub_age,c("gender","gender"),info)$table
          }
         else {
         sub_adonis_res_age <- my_adonis_terms_noadjAbun_noadjcov(sub_msv_dist_std, covar_sub,info)$table
         }
        }
       }

     else{
       left_covar<-subset(covariate,covariate!="age" & covariate!="BMI" & covariate!="gender" & covariate!="gender_code")
       ncovar<-length(left_covar)
       if(ncovar==1){
          covar_sub<-sub_covar[,c(left_covar,left_covar)]
         if(left_covar=="IBD"){
          sub_adonis_res_others <- my_adonis_terms_noadjAbun_noadjcov(sub_msv_dist_std, covar_sub,info)$table
          }
         else{
          sub_adonis_res_others <- my_adonis_terms_noadjAbun(sub_msv_dist_std, covar_sub, sub_covar,c("BMI","age","gender"),info)$table
          }
         }
       if(ncovar>1){
         if(studyname=="LiJ_2014"){
          covar_sub<-sub_covar[,left_covar]
          sub_adonis_res_others <- my_adonis_terms_noadjAbun(sub_msv_dist_std, covar_sub, sub_covar,c("BMI","gender"),info)$table
          }
        else{
          covar_sub<-sub_covar[,left_covar]
          sub_adonis_res_others <- my_adonis_terms_noadjAbun(sub_msv_dist_std, covar_sub, sub_covar,c("BMI","age","gender"),info)$table
          }
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
binary_var<-c("CRC","IBD","T2D","MetS","CAD","ACVD","cirrhosis","adenoma")

#datasets<-datasets[-seq(1,4)]

for (studyname in datasets){
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
  sub_abun_clr<-abundances(x=as.data.frame(na.omit(sub_abun)), transform="clr") %>%as.data.frame
  sub_abun_clr <- sub_abun_clr[match(rownames(sub_abun), rownames(sub_abun_clr)),]
  rownames(sub_abun_clr) <- rownames(sub_abun)
  sub_covar<-cbind(sub_basic,sub_abun_clr)

  covariate<-unlist(strsplit(dataset_variable[dataset_variable$study_name==studyname,]$variable,","))

  for (var in covariate){
     if(var=="BMI"){
       covar_sub<-sub_covar[,c("BMI","BMI")]
       sub_BMI<-subset(sub_covar,disease=="healthy",drop=T)
       covar_BMI_sub<-covar_sub[rownames(sub_BMI),]
       sub_abun_BMI_clr<-sub_abun_clr[rownames(sub_BMI),]
       if(nrow(sub_BMI)>=30){
        if("gender" %in% covariate){
          sub_lm_res_BMI <- lm_abun_index(covar_BMI_sub,sub_abun_BMI_clr,sub_BMI,c("gender_code","gender_code"))
          }
        else {
         sub_lm_res_BMI <- lm_abun_index_nocov(covar_sub,sub_abun_clr)
         }
        }
       }

     if(var=="age"){
       covar_sub<-sub_covar[,c("age","age")]
       sub_age<-subset(sub_covar,disease=="healthy",drop=T)
       covar_age_sub<-covar_sub[rownames(sub_age),]
       sub_abun_age_clr<-sub_abun_clr[rownames(sub_age),]
       if(nrow(sub_age)>=30){
        if("gender" %in% covariate){
          sub_lm_res_age <- lm_abun_index(covar_age_sub,sub_abun_age_clr,sub_age,c("gender_code","gender_code"))
          }
         else {
          sub_lm_res_age <- lm_abun_index_nocov(covar_sub,sub_abun_clr)
          }
         }
        }

    else{
       left_covar<-subset(covariate,covariate!="age" & covariate!="BMI" & covariate!="gender_code" & covariate!="gender")
       binary_covar<-subset(left_covar,left_covar %in% binary_var)
       ncovar_binary<-length(binary_covar)
       if(ncovar_binary==1){
          covar_sub_binary<-sub_covar[,c(binary_covar,binary_covar)]
         if(binary_covar=="IBD"){
          sub_lg_res_binary <- lg_abun_index_nocov(covar_sub_binary,sub_abun_clr)
          }
         else{
          sub_lg_res_binary <- lg_abun_index(covar_sub_binary,sub_abun_clr,sub_covar,c("BMI","age","gender_code"))
          }
         }
       if(ncovar_binary>1){
          covar_sub_binary<-sub_covar[,binary_covar]
          sub_lg_res_binary <- lg_abun_index(covar_sub_binary,sub_abun_clr,sub_covar,c("BMI","age","gender_code"))
        }

       quant_covar<-subset(left_covar,!left_covar %in% binary_var)
       ncovar_quant<-length(quant_covar)
       if(ncovar_quant==1){
          covar_sub_quant<-sub_covar[,c(quant_covar,quant_covar)]
          sub_lm_res_quant <- lm_abun_index(covar_sub_quant,sub_abun_clr,sub_covar,c("BMI","age","gender_code"))
         }
       if(ncovar_quant>1){
         if(studyname=="LiJ_2014"){
          covar_sub_quant<-sub_covar[,quant_covar]
          sub_lm_res_quant <- lm_abun_index(covar_sub_quant,sub_abun_clr,sub_covar,c("BMI","gender_code"))
          }
         else{
          covar_sub_quant<-sub_covar[,quant_covar]
          sub_lm_res_quant <- lm_abun_index(covar_sub_quant,sub_abun_clr,sub_covar,c("BMI","age","gender_code"))
          }
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
     covarites_lm <- my_batch_meta_p_lr(tmp_sub_lm_res,seq(3,ncol(tmp_sub_lm_res),by=5),seq(4,ncol(tmp_sub_lm_res),by=5))
     outfile<-paste("06.Species_genetic_association/",covariate,"_abun.csv",sep="")
     write.csv(covarites_lm, outfile,row.names=F)
    }
  }





######################################################################################################
### 3 Visualization 
### 3.1 Preparation

info<- read.table("01.cleanData/SV_info_QC/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

### SVs association

variable<-variable_dataset$variable
study_num<-length(variable)

## the first variable 
tmp_adonis<-read.csv(paste("06.Species_genetic_association/",variable[1],"_adonis.csv",sep=""),header=T)

### nromal p value 
#tmp_adonis<-tmp_adonis[,c(1,(ncol(tmp_adonis)-4))]

## meta fdr p 
tmp_adonis<-tmp_adonis[,c(1,(ncol(tmp_adonis)-1))]
tmp_adonis$Prog<-variable[1]
colnames(tmp_adonis)<-c("Species","fdr.p","Prog")

  for(k in 2:study_num){
    tmp_adonis_rest<-read.csv(paste("06.Species_genetic_association/",variable[k],"_adonis.csv",sep=""),header=T)
    tmp_adonis_rest<-tmp_adonis_rest[,c(1,(ncol(tmp_adonis_rest)-1))]
    tmp_adonis_rest$Prog<-variable[k]
    colnames(tmp_adonis_rest)<-c("Species","fdr.p","Prog")
    tmp_adonis<-rbind(tmp_adonis,tmp_adonis_rest)
    }

all_adonis.table<-tmp_adonis

### just select the 79 species
all_adonis.table<-subset(all_adonis.table,Species %in% info$Short_name,drop=T)
save(all_adonis.table, file = "06.Species_genetic_association/RData/all_adonis.table.RData")


#### abdundance significant
variable<-variable_dataset$variable

study_num<-length(variable)

## the first variable 
tmp_abun<-read.csv(paste("06.Species_genetic_association/",variable[1],"_abun.csv",sep=""),header=T)

## normal p 
## tmp_abun<-tmp_abun[,c(1,(ncol(tmp_abun)-5))]

## fdr p value 
tmp_abun<-tmp_abun[,c(1,ncol(tmp_abun))]
tmp_abun$Prog<-variable[1]
colnames(tmp_abun)<-c("Species","fdr.p","Prog")

  for(k in 2:study_num){
    tmp_abun_rest<-read.csv(paste("06.Species_genetic_association/",variable[k],"_abun.csv",sep=""),header=T)
    tmp_abun_rest<-tmp_abun_rest[,c(1,ncol(tmp_abun_rest))]
    tmp_abun_rest$Prog<-variable[k]
    colnames(tmp_abun_rest)<-c("Species","fdr.p","Prog")
    tmp_abun<-rbind(tmp_abun,tmp_abun_rest)
    }

all_abun.table<-tmp_abun
all_abun.table$Species<-info$Short_name[match(all_abun.table$Species, info$organism)]

### just select the 79 species
all_abun.table<-subset(all_abun.table,is.na(Species)==F,drop=T)

sv_assoc_id<-paste(all_adonis.table$Species, all_adonis.table$Prog,sep = "_")
abun_assoc_id<-paste(all_abun.table$Species, all_abun.table$Prog, sep = "_")

all_abun.table<-all_abun.table[match(sv_assoc_id,abun_assoc_id),]
colnames(all_abun.table)<-paste("Abun",colnames(all_abun.table),sep = ".")

save(all_abun.table, file = "06.Species_genetic_association/RData/all_abun.table.RData")
### write.csv(all_abun.table,"all_abun.table.csv")


species.table<-cbind(all_adonis.table,all_abun.table)

species.table$sv.MetaSigAssoc<-rep('No', nrow(species.table))
species.table$sv.MetaSigAssoc[species.table$fdr.p< 0.05]<-'Yes'
species.table$sv.MetaSigAssoc[is.na(species.table$fdr.p)==T]<-'Unknown'

species.table$abun.MetaSigAssoc<-rep('No', nrow(species.table))
species.table$abun.MetaSigAssoc[species.table$Abun.fdr.p< 0.05]<-'Yes'

### just include the species which used in the association analysis
species.table<-subset(species.table,Species %in% info$Short_name,drop=T)

species.sig.table<-species.table[,] %>%
  .[.$sv.MetaSigAssoc=="Yes" | .$abun.MetaSigAssoc == "Yes",]

write.table(species.sig.table, "06.Species_genetic_association/species.sig.tsv",sep = "\t", 
   col.names = T, row.names = F, quote = F)

write.csv(species.sig.table, "06.Species_genetic_association/species.sig.csv",row.names = F)

unique(species.sig.table[species.sig.table$sv.MetaSigAssoc=="Yes",]$Species)
unique(species.sig.table[species.sig.table$sv.MetaSigAssoc=="Yes",]$Prog)

unique(species.sig.table[species.sig.table$abun.MetaSigAssoc=="Yes",]$Species)
unique(species.sig.table[species.sig.table$abun.MetaSigAssoc=="Yes",]$Prog)

species.sv.sig.table<-species.sig.table[species.sig.table$sv.MetaSigAssoc=="Yes",]

#### more than two datasets
#cor.test(species_Prog.sv.sig.table$PRJEB22863.R2,species_Prog.sv.sig.table$PRJNA397906.R2)
#plot(species_Prog.sv.sig.table$PRJEB22863.R2,species_Prog.sv.sig.table$PRJNA397906.R2)

### 3.4 Venn diagram
species_count<-table(species.table$sv.MetaSigAssoc,species.table$abun.MetaSigAssoc)

#pdf("06.Species_genetic_association/species_Prog_count.venn.pdf", width = 2, height = 2)
tiff(file = "pics/FS2C_species_pvalue_count.venn.tiff", width =500, height =500, res =300) 

draw.pairwise.venn(species_count[3,1]+species_count[3,2],
                   species_count[1,2]+species_count[2,2]+species_count[3,2],
                   species_count[3,2], 
                   category = c("Genetics", "Abundance"), lty = rep("blank",2), 
                   fill =c("#4472c4", "#00b050"), alpha = rep(0.5, 2), cat.pos = c(0, 0), 
                   cat.cex = c(0.5, 0.5),cat.dist = rep(0.025, 2), scaled = F)

dev.off()
 


##############################################################################
###  3.5 heatmap
load("06.Species_genetic_association/RData/all_adonis.table.Rdata")
load("06.Species_genetic_association/RData/all_abun.table.RData")

## matrix
adonis_res.P<-all_adonis.table[,c("Prog", "Species", "fdr.p")] %>%
  spread("Prog", "fdr.p")
adonis_res.P<-data.frame(adonis_res.P, row.names = "Species")

### just include the species which used in the association analysis
adonis_res.P<-adonis_res.P[info$Short_name,]

all_abun_res.P<-all_abun.table[,c("Abun.Prog", "Abun.Species", "Abun.fdr.p")] %>%
  spread("Abun.Prog", "Abun.fdr.p")

all_abun_res.P<-data.frame(all_abun_res.P, row.names = "Abun.Species")
all_abun_res.P<-all_abun_res.P[match(rownames(adonis_res.P),rownames(all_abun_res.P)),
                                                     match(colnames(adonis_res.P),colnames(all_abun_res.P))]

### just include the species which used in the association analysis
all_abun_res.P<-all_abun_res.P[info$Short_name,]

# color
species.color<-matrix(0, nrow = nrow(adonis_res.P), ncol = ncol(adonis_res.P))
rownames(species.color)<-rownames(adonis_res.P)
colnames(species.color)<-colnames(adonis_res.P)

species.color[adonis_res.P<0.05&all_abun_res.P<0.05]<- 1        #"green", both assoc
species.color[adonis_res.P>=0.05&all_abun_res.P<0.05]<- 2       #"yellow", abundance assoc only 
species.color[is.na(adonis_res.P)==T&all_abun_res.P<0.05&!is.na(all_abun_res.P)]<- 3   #"grey", abundance assoiation and adonis_res.P not avaiable
species.color[adonis_res.P<0.05&all_abun_res.P>=0.05]<- 4       #"blue", genetic assoc only 


### association number 
species.asso_number<-matrix(0, nrow = nrow(adonis_res.P), ncol = ncol(adonis_res.P))
rownames(species.asso_number)<-rownames(adonis_res.P)
colnames(species.asso_number)<-colnames(adonis_res.P)

species.asso_number[adonis_res.P<0.05&all_abun_res.P<0.05]<- 2        #both assoc
species.asso_number[adonis_res.P>=0.05&all_abun_res.P<0.05]<- 1       #abundance assoc
species.asso_number[adonis_res.P<0.05&all_abun_res.P>=0.05]<- 1       #genetic assoc
species.asso_number[is.na(adonis_res.P)==T&all_abun_res.P<0.05]<- 1   #just abundance

species.asso<-apply(species.asso_number,1,sum)
species.asso<-species.asso[which(species.asso>0)]
species.asso<-species.asso[order(species.asso,na.last = TRUE, decreasing = TRUE)]
species_order<-names(species.asso)

## get plot tables
species.sig.table<-species.table[species.table$fdr.p<0.05 | species.table$Abun.fdr.p<0.05,]
species.sig.table<-subset(species.sig.table,is.na(sv.MetaSigAssoc)==F,drop=T)

##write.csv(species.color,"test.csv")
species.plot<- species.sig.table$Prog %>%
  as.character(.) %>%
  na.omit %>%
  .[!duplicated(.)]

species.plot.spe <- species.sig.table$Species %>%
  as.character(.) %>%
  .[!duplicated(.)] %>%
  sort(.,decreasing = F)

species.color.plot    <- species.color %>%
  .[match(species.plot.spe,rownames(.)),match(species.plot,colnames(.))]

## re_order 
species.color.plot<-species.color.plot[,c("age","BMI","IBD","CRC","T2D","hba1c",
"fasting_glucose","fasting_insulin","CAD","ACVD","triglycerides",
"hdl","cirrhosis","alt","ast")]

colnames(species.color.plot)<-c("Age","BMI","IBD","CRC","T2D","HbA1c",
              "Fasting glucose","Fasting insulin","CAD","ACVD","Triglycerides",
              "HDL","Cirrhosis","Alt","Ast")

species.color.plot<-species.color.plot[species_order,]

tiff(file = "pics/F4A_specie_prognosis.heatmap_fdr_value.tiff", width =4200, height =5200, res =300) 

heatmap.2(species.color.plot, 
          col=colorRampPalette(c("#e9e9e9","#207F4C","#ECD452","#EA8958","#2983BB"))(5), 
          trace = "none", Rowv = F, Colv = F, dendrogram = "none",
          density.info="none",
          cexCol =1.5, srtCol = 58, cexRow =1.5,
          labRow=as.expression(lapply(rownames(species.color.plot), function(a) bquote(italic(.(a))))),
          colsep = c(1:(ncol(species.color.plot)-1)), rowsep = c(1:(nrow(species.color.plot)-1)),
          sepcolor="#30303a",sepwidth = c(0.02,0.02),
          key = F,
          lmat=rbind( c(4, 4, 3), c(2, 1, 0 )), lhei = c(0.1, 3.5),lwid=c(3, 8, 4.2 ),key.title = NA,
          margins=c(16,8))

legend('bottomright', c("Both assocation","Abundance association only","Abundance association and SV not avaiable",
"SV association only","No association"), 
box.lty =0,pch=15,
col=c("#207F4C","#ECD452", "#EA8958","#2983BB","#e9e9e9"), cex=1)

dev.off()


#hist(mtcars$mpg,
#     breaks = 15,
#     col ="#EA8958",
#)

# "#6F99ADFF","#E18727FF","#BC3C29FF","#7876B1FF","#0072B5FF",#20854EFF
# c("#ECD452","#BB97C5","#EA8958","#2983BB","#BC3523","#207F4C","#EE4863",
# "#97C24E","#98369E","#63BBD0","#FEBA07","#EA517F","#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF")



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

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

exp <- substitute(
  atop(
    atop(Age (AsnicarF_2021), italic(F.prausnitzii)),
    atop(adonis~R ^ 2 == 0.01   ~italic(P) == 0.001,phantom(0)  )
  ),
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

tiff(file = "pics/F4C_AsnicarF_2021_age_F.prausnitzii.tiff", width =1000, height =1300, res =300) 
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

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

exp <- substitute(
  atop(
    atop(BMI (AsnicarF_2021), italic(F.prausnitzii)),
    atop(adonis~R ^ 2 == 0.006   ~italic(P) == 0.001,phantom(0)  )
  ),
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

tiff(file = "pics/F4C_AsnicarF_2021_BMI_F.prausnitzii.tiff", width =1000, height =1300, res =300) 
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

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

exp <- substitute(
  atop(
    atop("IBD" (FranzosaEA_2018), italic(F.prausnitzii)),
    atop(adonis~R ^ 2 == 0.05   ~italic(P) == 0.001,phantom(0)  )
  ),
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

tiff(file = "pics/F4C_F.prausnitzii_FranzosaEA_2018_IBD.tiff", width =1800, height =2500, res =600) 
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

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

exp <- substitute(
  atop(
    atop("IBD" (HallAB_2017), italic(E.rectale)),
    atop(adonis~R ^ 2 == 0.03   ~italic(P) == 0.001,phantom(0)  )
  ),
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

tiff(file = "pics/F4C_E.rectale_HallAB_2017_IBD.tiff", width =1800, height =2500, res =600) 
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

adonis_test<-adonis2(as.dist(species_dis_input) ~prog_vec_input,  permutations = 999) 
pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

exp <- substitute(
  atop(
    atop("CRC" (ZellerG_2014), italic(E.rectale)),
    atop(adonis~R ^ 2 == 0.02   ~italic(P) == 0.057,phantom(0)  )
  ),
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

tiff(file = "pics/F4C_E.rectale_ZellerG_2014_CRC.tiff", width =1800, height =2500, res =600) 
print(p_spe_pcoa)
dev.off()



############################################################
## output for Table 3  

variable<-variable_dataset$variable

## 去掉单数据集的变量
variable_more<-variable[-c(5,7,8,9,11)]
study_num<-length(variable_more)

for(k in 1:study_num){
  tmp_adonis<-read.csv(paste("06.Species_genetic_association/",variable_more[k],"_adonis.csv",sep=""),header=T)

  tmp_adonis_meta<-tmp_adonis[,c(1,2,seq((ncol(tmp_adonis)-4),ncol(tmp_adonis)))]
  tmp_adonis_meta<-tmp_adonis_meta[,c("Species","Meta.p","Meta.fdr.p")]

  tmp_adonis_dataset<-tmp_adonis[,seq(1,(ncol(tmp_adonis)-5))]
  tmp_adonis_dataset<-tmp_adonis_dataset[,-which(grepl("fdr",colnames(tmp_adonis_dataset))=="TRUE")]
  tmp_adonis_dataset<-tmp_adonis_dataset[,-2]
  tmp_adonis_reorder<-merge(tmp_adonis_meta,tmp_adonis_dataset,by.x="Species",by.y="Species")

  tmp_abun<-read.csv(paste("06.Species_genetic_association/",variable_more[k],"_abun.csv",sep=""),header=T)

  tmp_abun_meta<-tmp_abun[,c(1,2,seq((ncol(tmp_abun)-6),ncol(tmp_abun)))]
  tmp_abun_meta<-tmp_abun_meta[,c("X","effect_meta","pvalue_meta","fdr")]

  tmp_abun_dataset<-tmp_abun[,seq(1,(ncol(tmp_abun)-7))]
  tmp_abun_dataset<-tmp_abun_dataset[,-which(grepl("fdr",colnames(tmp_abun_dataset))=="TRUE")]
  tmp_abun_dataset<-tmp_abun_dataset[,-which(grepl("SE",colnames(tmp_abun_dataset))=="TRUE")]
  tmp_abun_dataset<-tmp_abun_dataset[,-2]

  tmp_abun_reorder<-merge(tmp_abun_meta,tmp_abun_dataset,by.x="X",by.y="X")
  colnames(tmp_abun_reorder)<-paste("Abun",colnames(tmp_abun_reorder),sep = ".")

  tmp_abun_reorder$Species<-info$Short_name[match(tmp_abun_reorder$Abun.X, info$organism)]

  tmp<-merge(tmp_adonis_reorder,tmp_abun_reorder,by.x="Species",by.y="Species")
  tmp<-tmp[,-which(grepl(".X",colnames(tmp))=="TRUE")]

  tmp$sv.MetaSigAssoc<-rep('No', nrow(tmp))
  tmp$sv.MetaSigAssoc[tmp$Meta.fdr.p< 0.05]<-'Yes'
  tmp$sv.MetaSigAssoc[is.na(tmp$Meta.fdr.p)==T]<-'Unknown'

  tmp$abun.MetaSigAssoc<-rep('No', nrow(tmp))
  tmp$abun.MetaSigAssoc[tmp$Abun.fdr< 0.05]<-'Yes'
  tmp$abun.MetaSigAssoc[is.na(tmp$Abun.fdr)==T]<-'Unknown'

  tmp_out<-paste("06.Species_genetic_association/tableS3/",variable_more[k],".csv",sep="")

  write.csv(tmp,tmp_out,row.names=F)
}



####################################################
## signal data

variable_signal<-variable[c(5,7,8,9,11)]
study_num<-length(variable_signal)

for(k in 1:study_num){
  tmp_adonis<-read.csv(paste("06.Species_genetic_association/",variable_signal[k],"_adonis.csv",sep=""),header=T)

  tmp_adonis_dataset<-tmp_adonis[,seq(1,(ncol(tmp_adonis)-5))]
  tmp_adonis_dataset<-tmp_adonis_dataset[,-2]

  tmp_abun<-read.csv(paste("06.Species_genetic_association/",variable_signal[k],"_abun.csv",sep=""),header=T)

  tmp_abun_dataset<-tmp_abun[,-which(grepl("SE",colnames(tmp_abun))=="TRUE")]
  tmp_abun_dataset<-tmp_abun_dataset[,-2]

  colnames(tmp_abun_dataset)<-paste("Abun",colnames(tmp_abun_dataset),sep = ".")

  tmp_abun_dataset$Species<-info$Short_name[match(tmp_abun_dataset$Abun.X, info$organism)]

  tmp<-merge(tmp_adonis_dataset,tmp_abun_dataset,by.x="Species",by.y="Species")

  tmp<-tmp[,-which(grepl(".X",colnames(tmp))=="TRUE")]

  tmp$sv.MetaSigAssoc<-rep('No', nrow(tmp))
  tmp$sv.MetaSigAssoc[tmp[,5]]<-'Yes'
  tmp$sv.MetaSigAssoc[is.na(tmp[,5])==T]<-'Unknown'

  tmp$abun.MetaSigAssoc<-rep('No', nrow(tmp))
  tmp$abun.MetaSigAssoc[tmp[,9]< 0.05]<-'Yes'
  tmp$abun.MetaSigAssoc[is.na(tmp[,9])==T]<-'Unknown'

  tmp_out<-paste("06.Species_genetic_association/tableS3/",variable_signal[k],".csv",sep="")

  write.csv(tmp,tmp_out,row.names=F)
}













#######################################################
##  calculate number  (should be modified)

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



