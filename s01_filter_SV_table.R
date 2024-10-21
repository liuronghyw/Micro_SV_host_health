# Filter the SV table to remove SVs with low call rate and low MAF and remove samples with low call rate.
# Write filtered SV tables per cohort and a file showing for each SV the cohorts that has it
# LiuRong
# 2024-08-01
# R 4.0.5 

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/00.rawData")

args <- commandArgs(trailingOnly = TRUE)

##install.packages("UpSetR")

library(UpSetR)
library(tibble)
library(dplyr)

# Filter the per-cohort SV table by call rate (@param cr), MAF (5%),and number of cases (>10) and controls (>10);
# Draw diagnostic plots, save to @param outpath
# Functions


run_qc_per_dsv <- function(d, cohort_name, outpath, cr = 0.1){
  qc <- as.data.frame(t(apply(d, 2, function(x) table(factor(x, levels = c(0,1,NA)), exclude = NULL))))
  
  qc$num_called <- qc[,1] + qc[,2]
  qc$presence_rate <- qc[,2]/qc$num_called
  qc$call_rate <- qc$num_called / (qc[,1] + qc[,2] + qc[,3])
  
  write.table(qc, file = paste0(outpath, ".qc_per_dSV.txt"),sep = "\t", quote = F, col.names = NA)
  
  # Make plots:
  pdf(paste0(outpath, ".qc_per_dSV.CR",cr, ".pdf"), useDingbats = F)
  par(mfrow=c(2,2))
  hist(qc$num_called, breaks = 100, main = "Number of called dSV", col = "black")
  hist(qc$call_rate, breaks = 100, xlim = c(0,1),  main = "dSV call rate", col = "black")
  abline(v=cr, col = "red")
  plot(qc$presence_rate, qc$call_rate, pch = 16, cex = 0.5, col = "black", main = "Fraction of 1s vs call rate")
  abline(v=0.05, col = "red")
  abline(v=0.95, col = "red")
  abline(h=cr, col = "red")
  hist(qc$presence_rate, breaks = 100, xlim = c(0,1), col = "black", main = "dSV fraction of 1")
  abline(v=0.05, col = "red")
  abline(v=0.95, col = "red")
  dev.off()
  
  # Get the filtered dSVs
  ##qc_05 <- qc[qc[,1] > 10 & qc[,2] > 10 & qc$presence_rate > 0.05 & qc$presence_rate < 0.95 & qc$call_rate > cr,]

  qc_05 <- qc[qc$presence_rate > 0.05 & qc$presence_rate < 0.95 & qc$call_rate > cr,]
  
  cat(cohort_name, ": number dSV before QC:", nrow(qc), "number dSV after filtering:", nrow(qc_05), "\n", sep = " ")
  
  d <- d[colnames(d) %in% row.names(qc_05)]
  return(d)
}


# Filter the per-cohort SV table by call rate (@param cr),call rate 10%;
# Draw diagnostic plots, save to @param outpath
run_qc_per_vsv <- function(d, cohort_name, outpath, cr = 0.1){
  qc <- as.data.frame(sapply(d, function(y) length(which(! is.na(y)))))
  
  n=nrow(d)
  nrow(qc)
  colnames(qc) <- "num_called"
  qc$call_rate <- qc[,1]/n
  
  # Make plots:
  pdf(paste0(outpath, ".qc_per_vSV.CR",cr, ".pdf"), useDingbats = F)
  par(mfrow=c(1,2))
  hist(qc$num_called, breaks = 100, main = "Number of called vSV", col = "black")
  hist(qc$call_rate, breaks = 100, xlim = c(0,1),  main = "vSV call rate", col = "black")
  abline(v=cr, col = "red")
  
  dev.off()
  
  # Get the filtered vSVs
  qc_05 <- qc[qc$call_rate > cr,]
  
  cat(cohort_name, ": number vSV before QC:", nrow(qc), "number vSV after filtering:", nrow(qc_05), "\n", sep = " ")
  d <- d[colnames(d) %in% row.names(qc_05)]
  
  return(d)
}


# Remove the samples with low call rate (@param cr)
run_qc_per_sample <- function(d, cohort_name, outpath, cr = 0.05){
  qc2 <-as.data.frame(apply(d, 1, function(y) length(which(! is.na(y)))))
  n2 <- ncol(d)
  colnames(qc2) <- "num_called"
  
  qc2$call_rate <- qc2$num_called/n2
  
  # Make plots
  pdf(paste0(outpath, ".qc_per_sample.pdf"), useDingbats = F)
  hist(qc2$call_rate, breaks = 100, main = "Call rate per sample", col = "black")
  abline(v=cr, col = "red")
  dev.off()
  
  # Get the filtered samples
  qc2_05 <- qc2[qc2$call_rate > cr,]
  cat(cohort_name, ": number samples before QC:", nrow(qc2), "number samples after filtering:", nrow(qc2_05), "\n", sep = " ")
  
  d <- d[row.names(d) %in% row.names(qc2_05),]
  return(d)
}


#
# Main
#
# Raw SV table:

dsv<-read.csv("SV/Whole_dsgv_20240801.csv",header=T)
vsv<-read.csv("SV/Whole_vsgv_20240801.csv",header=T)

colnames(dsv)[1]<-"sample"
colnames(vsv)[1]<-"sample"

# Table with cohort info for each sample id
clinical<-read.csv("clinical/clinical_whole.csv",header=T)

#for(i in 1:nrow(clinical)){
#  if (clinical$study_name[i]=="HirotsuguS_2022"){
#    clinical$match_id[i]<-gsub("X","",clinical$match_id[i])
#    }
#  }

cohorts<-clinical[,c("match_id","study_name")]
cs<-as.character(data.frame(table(cohorts$study_name))$Var1)

dsv_per_cohort <- matrix(nrow = ncol(dsv), ncol = length(cs))
row.names(dsv_per_cohort) <- colnames(dsv)
colnames(dsv_per_cohort) <- cs

##cs<-cs[c(1,2,3,4,5)]
#### filter dsv 

for (c in cs){
  #per cohort SV table
  # c<-"HirotsuguS_2022"
  d_cohort <- dsv[dsv[,1] %in% cohorts[cohorts$study_name == c, 1],]
  row.names(d_cohort) <- d_cohort[,1]
  d_cohort <- d_cohort[,2:ncol(d_cohort)]
  
  # get samples with covariate data
  clin<-subset(clinical,study_name==c,drop=T)
  cat("Number of samples with dSVs: ", nrow(d_cohort), "\n")
  cat("Number of samples with covariates: ", nrow(clin), "\n")
  sample_overlap <- row.names(d_cohort)[row.names(d_cohort) %in% clin$match_id]
  cat ("Number overlapping samples: ", length(sample_overlap), "\n")
  d_cohort <- d_cohort[row.names(d_cohort) %in% sample_overlap,]
  
  #filter SVs and samples
  d_flt <- run_qc_per_dsv(d_cohort, c, paste0("SV/3_SV_QC/studies/", c, ".dSV_filtering"))
  d_flt <- run_qc_per_sample(d_flt, c, paste0("SV/3_SV_QC/studies/", c, ".dSV_filtering"))
  dsv_per_cohort[row.names(dsv_per_cohort) %in% colnames(d_flt), c] <- 1
  write.table(d_flt, file = paste0("SV/3_SV_QC/dsv/", c, ".","dsv_filtered.txt"), sep = "\t", quote = F) 
}


#####################################################################################
## filter vSV 

vsv_per_cohort <- matrix(nrow = ncol(vsv), ncol = length(cs))
row.names(vsv_per_cohort) <- colnames(vsv)
colnames(vsv_per_cohort) <- cs

for (c in cs){
  #per cohort SV table
  d_cohort <- vsv[vsv[,1] %in% cohorts[cohorts$study_name == c, 1],]
  row.names(d_cohort) <- d_cohort[,1]
  d_cohort <- d_cohort[,2:ncol(d_cohort)]
  
  # get samples with covariate data
  clin<-subset(clinical,study_name==c,drop=T)
  cat("Number of samples with vSVs: ", nrow(d_cohort), "\n")
  cat("Number of samples with covariates: ", nrow(clin), "\n")
  sample_overlap <- row.names(d_cohort)[row.names(d_cohort) %in% clin$match_id]
  cat ("Number overlapping samples: ", length(sample_overlap), "\n")
  d_cohort <- d_cohort[row.names(d_cohort) %in% sample_overlap,]
  
  d_flt <- run_qc_per_vsv(d_cohort, c, paste0("SV/3_SV_QC/studies/", c, ".vSV_filtering"))
  d_flt <- run_qc_per_sample(d_flt, c, paste0("SV/3_SV_QC/studies/", c, ".vSV_filtering"))
    
  vsv_per_cohort[row.names(vsv_per_cohort) %in% colnames(d_flt), c] <- 1
    
  # normalize: apply INT
  d_norm <- apply(d_flt, 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
  d_flt <- as.data.frame(d_norm)
  write.table(d_flt, file = paste0("SV/3_SV_QC/vsv/", c, ".vsv_filtered.txt"), sep = "\t", quote = F) 
}



# plot number of dSVs and their overlap between cohorts
dsv_per_cohort <- as.data.frame(dsv_per_cohort)
dsv_per_cohort[is.na(dsv_per_cohort)] <- 0
pdf(paste0("SV/3_SV_QC/", "dsv_overlap_v3.pdf"),  useDingbats = F)
upset(dsv_per_cohort, order.by = "freq")
dev.off()


# plot number of vSVs and their overlap between cohorts
vsv_per_cohort <- as.data.frame(vsv_per_cohort)
vsv_per_cohort[is.na(vsv_per_cohort)] <- 0
pdf(paste0("SV/3_SV_QC/", "vsv_overlap_v3.pdf"),  useDingbats = F)
upset(vsv_per_cohort, order.by = "freq")
dev.off()



# Write dSV in any number of cohorts
dsv_per_cohort <- dsv_per_cohort[2:nrow(dsv_per_cohort),]
dsv_per_cohort$cohorts <- NA
dsv_per_cohort$cohorts <- apply(dsv_per_cohort, 1, function(x) {paste(colnames(dsv_per_cohort)[which(x==1)], collapse = ",")})
dsv_per_cohort <- dsv_per_cohort[dsv_per_cohort$cohorts != "",]

write.table(dsv_per_cohort, file = paste0("SV/3_SV_QC/", "dsv_1_cohorts.txt"), sep = "\t", quote = F, col.names = NA) 


# Write SVs present in > 2 cohorts
sum <- rowSums(dsv_per_cohort[,1:(ncol(dsv_per_cohort)-1)])
head(sum)
hist(sum)
dsv_per_cohort2 <- dsv_per_cohort[rowSums(dsv_per_cohort[,1:(ncol(dsv_per_cohort)-1)]) >2,]
write.table(dsv_per_cohort2, file = paste0("SV/3_SV_QC/", "dsv_2_cohort.txt"), sep = "\t", quote = F, col.names = NA)


# Write vSV in any number of cohorts
vsv_per_cohort <- vsv_per_cohort[2:nrow(vsv_per_cohort),]
vsv_per_cohort$cohorts <- NA
vsv_per_cohort$cohorts <- apply(vsv_per_cohort, 1, function(x) {paste(colnames(vsv_per_cohort)[which(x==1)], collapse = ",")})
vsv_per_cohort <- vsv_per_cohort[vsv_per_cohort$cohorts != "",]

write.table(vsv_per_cohort, file = paste0("SV/3_SV_QC/", "vsv_1_cohorts.txt"), sep = "\t", quote = F, col.names = NA) 


# Write SVs present in > 2 cohort

sum <- rowSums(vsv_per_cohort[,1:(ncol(dsv_per_cohort)-1)])
head(sum)
hist(sum)

vsv_per_cohort2 <- vsv_per_cohort[rowSums(vsv_per_cohort[,1:(ncol(dsv_per_cohort)-1)]) >2,]
write.table(vsv_per_cohort2, file = paste0("SV/3_SV_QC/", "vsv_2_cohort.txt"), sep = "\t", quote = F, col.names = NA)


### if dsv and vSV both identified, keep vSV
overlap <- intersect(rownames(dsv_per_cohort2),rownames(vsv_per_cohort2))
dsv_per_cohort2_sub<-dsv_per_cohort2[!rownames(dsv_per_cohort2) %in% overlap,]



#####################################################################
### get the dSV file after quality contral

files=list.files("SV/3_SV_QC/dsv") 
files<-data.frame(files)
colnames(files)<-"id"

dir=paste("SV/3_SV_QC/dsv",files$id,sep="/")
n = length(dir)   

#### the first one 
merge.data<-data.frame(t(read.table(dir[1])))
merge.data$id<-rownames(merge.data)

for (j in 2:n){
    data<-data.frame(t(read.table(dir[j])))
    data$id<-rownames(data)
    merge.data<-merge(merge.data,data,by.x="id",by.y="id",all.x=T,all.y=T)
 }


## filter the ones at least in two cohorts 
dsv_id_select<-data.frame(rownames(dsv_per_cohort2_sub))
colnames(dsv_id_select)<-"id"
merge.data<-merge(dsv_id_select,merge.data,by.x="id",by.y="id")

dsv_id<-read.csv("SV/dsgv_id_convert_20240801.csv",header=T)
merge.data$id<-gsub("X","",merge.data$id)

merge.data<-merge(dsv_id,merge.data,by.x="id",by.y="id")[,-1]
t_merge.data<-t(merge.data[,-1])
colnames(t_merge.data)<-merge.data[,1]

write.csv(t_merge.data,file = "SV/Whole_dsgv_QC.csv")  


#####################################################################
### get the vSV file after quality contral

### subset the ones that at least in 2 cohorts

files=list.files("SV/3_SV_QC/vsv") 
files<-data.frame(files)
colnames(files)<-"id"

dir=paste("SV/3_SV_QC/vsv",files$id,sep="/")
n = length(dir)   

#### the first one 
merge.data<-data.frame(t(read.table(dir[1])))
merge.data$id<-rownames(merge.data)

for (j in 2:n){
    data<-data.frame(t(read.table(dir[j])))
    data$id<-rownames(data)
    merge.data<-merge(merge.data,data,by.x="id",by.y="id",all.x=T,all.y=T)
 }

## filter the ones at least in two cohorts 
vsv_id_select<-data.frame(rownames(vsv_per_cohort2))
colnames(vsv_id_select)<-"id"
merge.data<-merge(vsv_id_select,merge.data,by.x="id",by.y="id")

vsv_id<-read.csv("SV/vsgv_id_convert_20240801.csv",header=T)
merge.data$id<-gsub("X","",merge.data$id)
merge.data<-merge(vsv_id,merge.data,by.x="id",by.y="id")[,-1]
t_merge.data<-t(merge.data[,-1])
colnames(t_merge.data)<-merge.data[,1]

write.csv(t_merge.data,file = "SV/Whole_vsgv_QC.csv")  




