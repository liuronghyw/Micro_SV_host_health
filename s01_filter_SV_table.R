# Filter the SV table to remove SVs with low call rate and low MAF
# Write filtered SV tables per cohort and a file showing for each SV the cohorts that has it
# 2026-09-24
# R 4.0.5 

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/00.rawData")

## Create output directories before the first write.
dir.create("SV/3_SV_QC", recursive = TRUE, showWarnings = FALSE)
dir.create("SV/3_SV_QC/studies", recursive = TRUE, showWarnings = FALSE)
dir.create("SV/3_SV_QC/dsv", recursive = TRUE, showWarnings = FALSE)
dir.create("SV/3_SV_QC/vsv", recursive = TRUE, showWarnings = FALSE)
dir.create("../01.cleanData/SV", recursive = TRUE, showWarnings = FALSE)
dir.create("../pics", recursive = TRUE, showWarnings = FALSE)
## End output directory initialization.
# Match the encoding of the clinical cohort names on Windows.
Sys.setlocale("LC_CTYPE","Chinese_China.936")


args <- commandArgs(trailingOnly = TRUE)

##install.packages("UpSetR")

library(UpSetR)
library(tibble)
library(dplyr)

library(ggplot2)


#################################################################################################
## Clinical information: association between read count and sample-level SV call rate

read_count <- read.csv(
  "F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/00.rawData/clinical/clinical_whole.csv",
  header = TRUE
)[, c("match_id", "read_count")]

dsv <- read.csv("SV/Whole_dsgv_20260921.csv", header = TRUE)
vsv <- read.csv("SV/Whole_vsgv_20260921.csv", header = TRUE)

# Combine dSV and vSV calls by sample. Each row is one independent sample and
# each remaining column is one SV.
sv <- merge(dsv, vsv, by = "X")
rownames(sv) <- sv$X
sv <- sv[, setdiff(colnames(sv), "X"), drop = FALSE]

# Sample-level call rate = number of non-NA SV calls / total number of SVs.
total_sv <- ncol(sv)
qc <- data.frame(
  match_id = rownames(sv),
  num_called = rowSums(!is.na(sv)),
  total_sv = total_sv,
  stringsAsFactors = FALSE
)
qc$call_rate <- qc$num_called / qc$total_sv

merge_file <- merge(qc, read_count, by = "match_id")
merge_file$read_count <- suppressWarnings(as.numeric(merge_file$read_count))

# Use only samples with finite read count and call rate for correlation/plotting.
analysis_data <- merge_file[
  is.finite(merge_file$read_count) & is.finite(merge_file$call_rate),,
  drop = FALSE
]


if (nrow(analysis_data) < 3L) {
  stop("Fewer than three samples have both read_count and call_rate.")
}
if (stats::sd(analysis_data$read_count) == 0 || stats::sd(analysis_data$call_rate) == 0) {
  stop("Pearson correlation cannot be calculated because one variable has zero variance.")
}

cor_test <- cor.test(
  analysis_data$read_count,
  analysis_data$call_rate,
  method = "pearson"
)

cor_summary <- data.frame(
  method = "Pearson",
  n = nrow(analysis_data),
  estimate = unname(cor_test$estimate),
  conf_low = cor_test$conf.int[1],
  conf_high = cor_test$conf.int[2],
  p_value = cor_test$p.value,
  p_value_display = format.pval(
    cor_test$p.value, digits = 3, eps = .Machine$double.eps
  ),
  stringsAsFactors = FALSE
)

dir.create("SV/3_SV_QC", recursive = TRUE, showWarnings = FALSE)
write.table(
  analysis_data,
  "SV/3_SV_QC/read_count_call_rate_source_data.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  cor_summary,
  "SV/3_SV_QC/read_count_call_rate_correlation.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

cor_label <- sprintf(
  "r = %.3f\np = %s",
  cor_summary$estimate,
  cor_summary$p_value_display
)

p <- analysis_data
p$read_count_million <- p$read_count / 1e6

p <- ggplot(p, aes(x = read_count_million, y = call_rate)) +
  geom_point(color = "#3182BD", alpha = 0.65, size = 1.6) +
  geom_vline(
    xintercept = 5, colour = "#4D4D4D", linetype = "dashed", size = 0.55
  ) +
  geom_hline(
    yintercept = 0.05, colour = "#4D4D4D", linetype = "dashed", size = 0.55
  ) +
  geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE,
    color = "#D24B40", fill = "#D24B40", alpha = 0.18,
    size = 0.6
  ) +
  annotate(
    "text", x = -Inf, y = Inf, label = cor_label,
    hjust = -0.1, vjust = 1.1, size = 3.0
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0)
  ) +
  labs(
    x = "Read count (million)",
    y = "SV call-rate ratio"
  ) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.line = element_line(size = 0.45, colour = "black"),
    axis.ticks = element_line(size = 0.45, colour = "black"),
    panel.grid.major = element_line(colour = "gray90", linetype = "dashed", size = 0.25),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.margin = margin(6, 7, 5, 6)
  )


# Match the requested raster export: 1200 x 800 pixels at 300 dpi.
tiff(
  file = "F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/pics/0_QC_plot.tiff",
  width = 1200, height = 800, res = 300, compression = "lzw"
)
print(p)
dev.off()


# Distribution of sample-level SV call-rate ratios.
p_call_rate_hist <- ggplot(analysis_data, aes(x = call_rate)) +
  geom_histogram(
    bins = 50, boundary = 0, colour = "white",
    fill = "#3182BD", alpha = 0.85
  ) +
  geom_vline(
    xintercept = 0.05, colour = "#4D4D4D", linetype = "dashed", size = 0.55
  ) +
  scale_x_continuous(
    limits = c(0, 1), breaks = seq(0, 1, by = 0.2), expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "SV call-rate ratio",
    y = "Number of samples"
  ) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.line = element_line(size = 0.45, colour = "black"),
    axis.ticks = element_line(size = 0.45, colour = "black"),
    panel.grid.major = element_line(colour = "gray90", linetype = "dashed", size = 0.25),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.margin = margin(6, 7, 5, 6)
  )

tiff(
  file = "F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/pics/0_call_rate_hist.tiff",
  width = 1200, height = 800, res = 300, compression = "lzw"
)
print(p_call_rate_hist)
dev.off()

print(cor_summary)
print(p)





# Filter each cohort: call rate >=10%, 0.05 < presence rate < 0.95, deletion and presence each >=10, non-NA sample size >=40.
# Draw diagnostic plots, save to @param outpath
# Functions

run_qc_per_dsv <- function(d, cohort_name, outpath, cr = 0.10){
  qc <- as.data.frame(t(apply(d, 2, function(x) table(factor(x, levels = c(0,1,NA)), exclude = NULL))))
  
  qc$num_called <- qc[,1] + qc[,2]
  qc$presence_rate <- qc[,2]/qc$num_called
  qc$call_rate <- qc$num_called / (qc[,1] + qc[,2] + qc[,3])
  qc$num_deletion <- qc[,1]
  qc$num_presence<- qc[,2]
  
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
  qc_05 <- qc[qc$presence_rate > 0.05 & qc$presence_rate < 0.95 & qc$call_rate >= cr &
             qc$num_called >=40 & qc$num_deletion>=10 &  qc$num_presence>=10,]
  
  cat(cohort_name, ": number dSV before QC:", nrow(qc), "number dSV after filtering:", nrow(qc_05), "\n", sep = " ")
  d <- d[colnames(d) %in% row.names(qc_05)]
  return(d)
}


# Filter each cohort: call rate >=10% and non-NA sample size >=40.
# Draw diagnostic plots, save to @param outpath
run_qc_per_vsv <- function(d, cohort_name, outpath, cr = 0.10){
  qc <- as.data.frame(sapply(d, function(y) length(which(! is.na(y)))))
  
  n=nrow(d)
  nrow(qc)
  colnames(qc) <- "num_called"
  qc$call_rate <- qc[,1]/n
  qc$called <- qc[,1]
  
  # Make plots:
  pdf(paste0(outpath, ".qc_per_vSV.CR",cr, ".pdf"), useDingbats = F)
  par(mfrow=c(1,2))
  hist(qc$num_called, breaks = 100, main = "Number of called vSV", col = "black")
  hist(qc$call_rate, breaks = 100, xlim = c(0,1),  main = "vSV call rate", col = "black")
  abline(v=cr, col = "red")
  
  dev.off()
  
  # Get the filtered vSVs
  qc_05 <- qc[qc$call_rate >= cr & qc$num_called>=40,]
  
  cat(cohort_name, ": number vSV before QC:", nrow(qc), "number vSV after filtering:", nrow(qc_05), "\n", sep = " ")
  d <- d[colnames(d) %in% row.names(qc_05)]
  
  return(d)
}


write_pdf_atomic <- function(file_name, plot_fun) {
  out_dir <- dirname(file_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tmp_file <- tempfile(pattern = basename(file_name), tmpdir = out_dir, fileext = ".tmp.pdf")
  device <- NULL
  on.exit(unlink(tmp_file), add = TRUE)
  tryCatch({
    pdf(tmp_file, useDingbats = FALSE)
    device <- dev.cur()
    plot_fun()
  }, finally = {
    if (!is.null(device) && device %in% dev.list()) dev.off(device)
  })
  if (file.exists(file_name)) unlink(file_name)
  if (!file.rename(tmp_file, file_name)) {
    if (!file.copy(tmp_file, file_name, overwrite = TRUE)) stop("Could not write PDF: ", file_name)
    unlink(tmp_file)
  }
}



####################################################################################
# Raw SV table:
dsv<-read.csv("SV/Whole_dsgv_20260921.csv",header=T)
vsv<-read.csv("SV/Whole_vsgv_20260921.csv",header=T)

# Table with cohort info for each sample id
clinical<-read.csv("clinical/clinical_whole_add_norm.csv",header=T)

# Exclude cohorts with no SVs passing QC (non-NA sample size >=40).
exclude_cohorts<-c("HallAB_2017")
clinical<-clinical[!clinical$study_name %in% exclude_cohorts,,drop=F]

# Use the centrally selected donors before cohort-level SV QC and normalization.
selected <- read.csv("../01.cleanData/donor_selection/s00_retained_samples.csv",
                     header = TRUE, fileEncoding = "GB18030")
clinical <- clinical[clinical$match_id %in% selected$match_id,,drop=FALSE]
stopifnot(!anyDuplicated(clinical$match_id))
cohorts<-clinical[,c("match_id","study_name")]

cs<-as.character(data.frame(table(cohorts$study_name))$Var1)

dsv_per_cohort <- matrix(nrow = ncol(dsv), ncol = length(cs))
row.names(dsv_per_cohort) <- colnames(dsv)
colnames(dsv_per_cohort) <- cs

##cs<-cs[c(1,2,3,4,5)]
#### filter dsv 

for (c in cs){
  #per cohort SV table
  # c<-"HMP_2019_ibdmdb"
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
  #d_flt <- run_qc_per_sample(d_flt, c, paste0("SV/3_SV_QC/studies/", c, ".dSV_filtering"), clin)

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
  #d_flt <- run_qc_per_sample(d_flt, c, paste0("SV/3_SV_QC/studies/", c, ".vSV_filtering"), clin)
    
  vsv_per_cohort[row.names(vsv_per_cohort) %in% colnames(d_flt), c] <- 1
    
  # normalize: apply INT
  if(ncol(d_flt)>0){
    d_norm <- apply(d_flt, 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
    d_flt <- as.data.frame(d_norm)
  }
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


# Write dSVs present in >= 2 cohorts
sum <- rowSums(dsv_per_cohort[,1:(ncol(dsv_per_cohort)-1)])
head(sum)
hist(sum)
dsv_per_cohort2 <- dsv_per_cohort[rowSums(dsv_per_cohort[,1:(ncol(dsv_per_cohort)-1)]) >=2,]
write.table(dsv_per_cohort2, file = paste0("SV/3_SV_QC/", "dsv_2_cohort.txt"), sep = "\t", quote = F, col.names = NA)

# Write the dSV cohort counts below, after removing overlap with vSV.


# Write vSV in any number of cohorts
vsv_per_cohort <- vsv_per_cohort[2:nrow(vsv_per_cohort),]
vsv_per_cohort$cohorts <- NA
vsv_per_cohort$cohorts <- apply(vsv_per_cohort, 1, function(x) {paste(colnames(vsv_per_cohort)[which(x==1)], collapse = ",")})
vsv_per_cohort <- vsv_per_cohort[vsv_per_cohort$cohorts != "",]

write.table(vsv_per_cohort, file = paste0("SV/3_SV_QC/", "vsv_1_cohorts.txt"), sep = "\t", quote = F, col.names = NA) 


# Write SVs present in >=2 cohort

sum <- rowSums(vsv_per_cohort[,1:(ncol(vsv_per_cohort)-1)])
head(sum)
hist(sum)

vsv_per_cohort2 <- vsv_per_cohort[rowSums(vsv_per_cohort[,1:(ncol(vsv_per_cohort)-1)]) >=2,]
write.table(vsv_per_cohort2, file = paste0("SV/3_SV_QC/", "vsv_2_cohort.txt"), sep = "\t", quote = F, col.names = NA)

write.csv(sum[sum>=2],file = paste0("SV/3_SV_QC/", "vsv_2_cohort.csv"))


### if dsv and vSV both identified, keep vSV

#vsv_per_cohort2<-read.table("SV/3_SV_QC/vsv_2_cohort.txt",header=T)
#dsv_per_cohort2<-read.table("SV/3_SV_QC/dsv_2_cohort.txt",header=T)

overlap <- intersect(rownames(dsv_per_cohort2),rownames(vsv_per_cohort2))
write.csv(overlap,"SV/3_SV_QC/overlap_sv.csv")

dsv_per_cohort2_sub<-dsv_per_cohort2[!rownames(dsv_per_cohort2) %in% overlap,]
write.table(dsv_per_cohort2_sub, file = paste0("SV/3_SV_QC/", "dsv_2_cohort.txt"), sep = "\t", quote = F, col.names = NA)

# Counts correspond to the final dSV list: >=2 QC-passing cohorts, excluding vSV overlap.
sum <- rowSums(dsv_per_cohort2_sub[,1:(ncol(dsv_per_cohort2_sub)-1)])
write.csv(sum,file = paste0("SV/3_SV_QC/", "dsv_2_cohort.csv"))



#####################################################################
### get the dSV file after quality contral

files=list.files("SV/3_SV_QC/dsv")
files=files[files %in% paste0(cs,".dsv_filtered.txt")]
files<-data.frame(files)
colnames(files)<-"id"

dir=paste("SV/3_SV_QC/dsv",files$id,sep="/")
n = length(dir)   

#### the first one 
merge.data<-data.frame(t(read.table(dir[1],header=T,row.names=1,blank.lines.skip=F)))
merge.data$id<-rownames(merge.data)

for (j in 2:n){
    data<-data.frame(t(read.table(dir[j],header=T,row.names=1,blank.lines.skip=F)))
    data$id<-rownames(data)
    merge.data<-merge(merge.data,data,by.x="id",by.y="id",all.x=T,all.y=T)
 }


## filter the ones at least in 2 cohorts 
dsv_id_select<-data.frame(rownames(dsv_per_cohort2_sub))
colnames(dsv_id_select)<-"id"
merge.data<-merge(dsv_id_select,merge.data,by.x="id",by.y="id")

dsv_id<-read.csv("SV/dsgv_id_convert_20260921.csv",header=T)
merge.data$id<-sub("^X(?=[0-9])","",as.character(merge.data$id),perl=TRUE)

merge.data<-merge(dsv_id,merge.data,by.x="id",by.y="id")[,-1]
t_merge.data<-t(merge.data[,-1])
colnames(t_merge.data)<-merge.data[,1]

write.csv(t_merge.data,file = "SV/Whole_dsgv_QC.csv")  


#####################################################################
### get the vSV file after quality contral

### subset the ones that at least in 2 cohorts

files=list.files("SV/3_SV_QC/vsv")
files=files[files %in% paste0(cs,".vsv_filtered.txt")]
files<-data.frame(files)
colnames(files)<-"id"

dir=paste("SV/3_SV_QC/vsv",files$id,sep="/")
n = length(dir)   

#### the first one 
merge.data<-data.frame(t(read.table(dir[1],header=T,row.names=1,blank.lines.skip=F)))
merge.data$id<-rownames(merge.data)

for (j in 2:n){
    data<-data.frame(t(read.table(dir[j],header=T,row.names=1,blank.lines.skip=F)))
    data$id<-rownames(data)
    merge.data<-merge(merge.data,data,by.x="id",by.y="id",all.x=T,all.y=T)
 }

## filter the ones at least in two cohorts 
vsv_id_select<-data.frame(rownames(vsv_per_cohort2))
colnames(vsv_id_select)<-"id"
merge.data<-merge(vsv_id_select,merge.data,by.x="id",by.y="id")

vsv_id<-read.csv("SV/vsgv_id_convert_20260921.csv",header=T)
merge.data$id<-sub("^X(?=[0-9])","",as.character(merge.data$id),perl=TRUE)
merge.data<-merge(vsv_id,merge.data,by.x="id",by.y="id")[,-1]
t_merge.data<-t(merge.data[,-1])
colnames(t_merge.data)<-merge.data[,1]

write.csv(t_merge.data,file = "SV/Whole_vsgv_QC.csv")  




##################################################################
## get the datasets file 

dsv_id<-read.csv("SV/dsgv_id_convert_20260921.csv",header=T)
dsv_cohort<-read.csv("SV/3_SV_QC/dsv_2_cohort.csv",header=T)
colnames(dsv_cohort)<-c("id","cohort_number")
dsv_cohort$id<-sub("^X(?=[0-9])","",as.character(dsv_cohort$id),perl=TRUE)
merge.data<-merge(dsv_id,dsv_cohort,by.x="id",by.y="id")[,-1]

write.csv(merge.data,"F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/01.cleanData/SV/dsv_cohort_number.csv",row.names=F)



vsv_id<-read.csv("SV/vsgv_id_convert_20260921.csv",header=T)
vsv_cohort<-read.csv("SV/3_SV_QC/vsv_2_cohort.csv",header=T)
colnames(vsv_cohort)<-c("id","cohort_number")
vsv_cohort$id<-sub("^X(?=[0-9])","",as.character(vsv_cohort$id),perl=TRUE)
merge.data<-merge(vsv_id,vsv_cohort,by.x="id",by.y="id")[,-1]

write.csv(merge.data,"F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/01.cleanData/SV/vsv_cohort_number.csv",row.names=F)




##########################################################################################
### replace the ids 

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/00.rawData")


replace_colnames <- function(df, mapping_df, old_col = "id", new_col = "id_final") {
  current_names <- colnames(df)
  replace_vec <- setNames(mapping_df[[new_col]], mapping_df[[old_col]])
  updated_names <- ifelse(current_names %in% names(replace_vec),
                          replace_vec[current_names],
                          current_names)
  colnames(df) <- updated_names
  return(df)
}


# Raw SV table:
dsv<-read.csv("SV/Whole_dsgv_QC.csv",header=T)
vsv<-read.csv("SV/Whole_vsgv_QC.csv",header=T)

colnames(dsv)[1]<-"sample"
colnames(vsv)[1]<-"sample"

vsv_id<-read.csv("SV/vsgv_id_convert_20260921.csv",header=T)
dsv_id<-read.csv("SV/dsgv_id_convert_20260921.csv",header=T)

dsv_id$id<-paste("X",dsv_id$id,sep="")
dsv_output <- replace_colnames(dsv, dsv_id)

vsv_id$id<-paste("X",vsv_id$id,sep="")
vsv_output <- replace_colnames(vsv, vsv_id)

rownames(dsv_output)<-dsv_output$sample
rownames(vsv_output)<-vsv_output$sample

dsv_out<-dsv_output[,-1]
vsv_out<-vsv_output[,-1]

write.csv(dsv_out,"SV/Whole_dsgv_QC_filter.csv")
write.csv(vsv_out,"SV/Whole_vsgv_QC_filter.csv")



# annotation file 
dsgv_anno<-unique(read.delim("SV/s02.dSVs_anno_20260921.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F))
vsgv_anno<-unique(read.delim("SV/s03.vSVs_anno_20260921.tsv",
                      sep = "\t",header = T,quote = '',stringsAsFactors = F))

overlap<-read.csv("SV/3_SV_QC/overlap_sv.csv")[,2]

overlap<-data.frame(overlap)
colnames(overlap)<-"sv"

overlap <- overlap %>%
  mutate(sv = sub("^X(?=[0-9])", "", as.character(sv), perl = TRUE))

dsgv_anno$SV_id_new<-gsub(":", ".", dsgv_anno$SV_id)

rownames(dsgv_anno)<-dsgv_anno$SV_id_new
dsgv_anno_sub<-dsgv_anno[!rownames(dsgv_anno) %in% overlap$sv,]

dsgv_anno_sub <- dsgv_anno_sub[, !(names(dsgv_anno_sub) %in% c("SV_id_new"))]


write.table(dsgv_anno_sub, "SV/s02.dSVs_anno_20260921_filter.tsv",sep = "\t", quote = F, col.names = T, row.names = F)
write.table(vsgv_anno, "SV/s03.vSVs_anno_20260921_filter.tsv",sep = "\t", quote = F, col.names = T, row.names = F)


