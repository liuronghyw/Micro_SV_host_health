### Microbial GWA  £¨plots)
### 2026-10-23
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")


## 1 Preparation
### 1.1 Import
#install.packages("survminer")
rm(list = ls())
#options(scipen = 200)

source("functions.R")
knitr::opts_chunk$set(echo = TRUE)
library(ggthemes)

#library("R.utils")
#library(showtext)
#showtext_auto()

#library("survival")
#library("ggplot2")
#library("survminer")
#library(gridExtra)
#library("grid")
library(reshape2)	  
#library("RColorBrewer")
#library("plyr")


###################### 1.2 Inputs
info<- read.table("01.cleanData/SV_info_QC/Informative_species_information_final.tsv",sep = "\t", header = T, stringsAsFactors = F)

vsv_info<-read.table("01.cleanData/SV_info_QC/vsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")
dsv_info<-read.table("01.cleanData/SV_info_QC/dsgv_info_anno.tsv",sep = "\t",header = T,stringsAsFactors = F, quote = "")

sv_info<-rbind(vsv_info,dsv_info)

all_basic<-read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv")


################################################################################
### 3.3 All-phenotype chord diagram: union of significant model1 and model2 SVs

library(circlize)

## Read significant results from both models, retaining all association statistics.
sig_dirs<-c("07.Microbial_GWAS/analysis_wide_sig/adjAbun",
            "07.Microbial_GWAS/analysis_wide_sig/adjAbunPCs")
sig_files<-unlist(lapply(sig_dirs,list.files,pattern="_sv[.]csv$",full.names=TRUE),use.names=FALSE)
stopifnot(all(dir.exists(sig_dirs)),length(sig_files)>0L)
sv_chord_edges<-data.frame()
for (sig_file in sig_files) {
  sig_table<-read.csv(sig_file,stringsAsFactors=FALSE,check.names=FALSE)
  if (nrow(sig_table)==0L) next
  sig_table$Pheno<-sub("_sv[.]csv$","",basename(sig_file))
  sig_table$SV<-sig_table$X
  sig_table$type<-tolower(sig_table$type)
  sig_table$Model<-ifelse(basename(dirname(sig_file))=="adjAbunPCs","model2","model1")
  sv_chord_edges<-dplyr::bind_rows(sv_chord_edges,sig_table)
}


## Prefer model2 for duplicate phenotype, SV and type combinations.
stopifnot(nrow(sv_chord_edges)>0L)
sv_chord_edges<-sv_chord_edges[order(sv_chord_edges$Model,decreasing=TRUE),,drop=FALSE]
sv_chord_edges<-sv_chord_edges[!duplicated(sv_chord_edges[,c("Pheno","SV","type")]),,drop=FALSE]
stopifnot(!anyNA(sv_chord_edges[,c("Pheno","SV","type","Taxonomy_Name","Model")]),
          all(nzchar(sv_chord_edges$SV)),all(sv_chord_edges$type %in% c("vsv","dsv")))
sv_chord_edges$Species<-info$Short_name[match(sv_chord_edges$Taxonomy_Name,info$organism)]
stopifnot(!anyNA(sv_chord_edges$Species),all(nzchar(sv_chord_edges$Species)))

phenotype_labels<-c(
  age="Age",BMI="BMI",adenoma="Adenoma",cirrhosis="Cirrhosis",hba1c="HbA1c",
  fasting_glucose="Fasting glucose",fasting_insulin="Fasting insulin",
  triglycerides="Triglycerides",hdl="HDL cholesterol",ldl="LDL cholesterol",
  cholesterol="Total cholesterol",mean_press="Mean arterial pressure",
  creatinine="Creatinine",ast="AST",alt="ALT",total_bilirubin="Total bilirubin"
)

sv_chord_edges$Pheno<-ifelse(sv_chord_edges$Pheno %in% names(phenotype_labels),
                           phenotype_labels[sv_chord_edges$Pheno],sv_chord_edges$Pheno)

## Summarize and order by significant SV counts; chord widths represent association counts.
sv_chord_count<-aggregate(list(Count=rep(1L,nrow(sv_chord_edges))),
                          by=sv_chord_edges[,c("Pheno","Species")],FUN=sum)
phenotype_order<-names(sort(tapply(sv_chord_count$Count,sv_chord_count$Pheno,sum),decreasing=TRUE))
species_order<-names(sort(tapply(sv_chord_count$Count,sv_chord_count$Species,sum),decreasing=TRUE))
sv_chord_count<-sv_chord_count[order(match(sv_chord_count$Pheno,phenotype_order),
                                    match(sv_chord_count$Species,species_order)),]
sector_order<-c(rev(phenotype_order),species_order)

## Use predefined colors from functions.R and swap the BMI and Cirrhosis colors.
stopifnot(length(phenotype_order)<=length(mycolor15))
phenotype_colors<-setNames(mycolor15[seq_along(phenotype_order)],phenotype_order)
phenotype_colors[c("BMI","Cirrhosis")]<-c(mycolor2_green_blue[1],mycolor2_blue_yellow[2])
sector_colors<-c(phenotype_colors,setNames(rep(mycolor2_pink_grey[2],length(species_order)),species_order))

write.csv(sv_chord_count,"07.Microbial_GWAS/all_phenotypes_union_sv_prog_count.csv",row.names=FALSE)
dir.create("pics",showWarnings=FALSE)


tiff(file="pics/F4A_all_phenotypes_union_circos.tiff",width=6500,height=5500,res=300,compression="lzw")

par(mar=c(0,0,0,0),xaxs="i",yaxs="i",family="sans")
circos.clear()
circos.par(start.degree=180,clock.wise=TRUE,canvas.ylim=c(-0.95,0.7))
label_track_height<-max(strwidth(sector_order,units="inches",cex=2,font=3))*25.4+4
chordDiagram(
  sv_chord_count,annotationTrack="grid",grid.col=sector_colors,
  col=unname(phenotype_colors[sv_chord_count$Pheno]),transparency=0.55,
  order=sector_order,small.gap=2,big.gap=10,reduce=0,
  preAllocateTracks=list(track.height=uh(label_track_height,"mm"))
)
circos.track(track.index=1,panel.fun=function(x,y) {
  circos.text(CELL_META$xcenter,CELL_META$ylim[1],CELL_META$sector.index,
               cex=2,facing="clockwise",niceFacing=TRUE,adj=c(0,0.5),
               font=ifelse(CELL_META$sector.index %in% phenotype_order,1,3))
},bg.border=NA)
circos.clear()

dev.off()


################################################################################
### 3.4 Oscillibacter sp. ER4: SV associations

library(gplots)

ER4_pheno<-c(BMI="BMI",age="Age")
ER4_sig<-sv_chord_edges[sv_chord_edges$Taxonomy_Name=="Oscillibacter sp. ER4" &
                        sv_chord_edges$Pheno %in% ER4_pheno,,drop=FALSE]
ER4_rows<-sort(unique(paste(ER4_sig$type,ER4_sig$SV,sep=":")))
ER4_columns<-c("Pheno","SV","type","Model","effect_meta","pvalue_meta","se_meta","I2","bonferroni","fdr")
ER4_results<-ER4_sig[,ER4_columns]

## Use the model union for significant cells; fill remaining cells from full results, preferring model2.
for (model in c("adjAbunPCs","adjAbun")) {
  for (pheno in names(ER4_pheno)) {
    ER4_full<-read.csv(file.path("07.Microbial_GWAS",paste0("anno_",model),paste0(pheno,"_sv.csv")),
                      stringsAsFactors=FALSE,check.names=FALSE)
    ER4_full<-ER4_full[ER4_full$Taxonomy_Name=="Oscillibacter sp. ER4" &
                       paste(ER4_full$type,ER4_full$X,sep=":") %in% ER4_rows,,drop=FALSE]
    if (nrow(ER4_full)==0L) next
    ER4_full$Pheno<-ER4_pheno[pheno]
    ER4_full$SV<-ER4_full$X
    ER4_full$Model<-ifelse(model=="adjAbunPCs","model2","model1")
    ER4_results<-rbind(ER4_results,ER4_full[,ER4_columns])
  }
}
ER4_results<-ER4_results[!duplicated(ER4_results[,c("Pheno","SV","type")]),]
ER4_results$Sign<-""
ER4_results$Sign[with(ER4_results,which(pvalue_meta<0.05 & I2<=75))]<-"*"
ER4_results$Sign[with(ER4_results,which(fdr<0.05 & I2<=75))]<-"\u2606"
ER4_results$Sign[with(ER4_results,which(bonferroni<0.05 & I2<=75))]<-"\u2605"




## Raw beta; missing values stay NA; display the strongest significance level
ER4_beta<-matrix(NA_real_,length(ER4_rows),length(ER4_pheno),dimnames=list(ER4_rows,unname(ER4_pheno)))
ER4_sign<-matrix("",nrow(ER4_beta),ncol(ER4_beta),dimnames=dimnames(ER4_beta))
ER4_index<-cbind(match(paste(ER4_results$type,ER4_results$SV,sep=":"),ER4_rows),
                  match(ER4_results$Pheno,colnames(ER4_beta)))
stopifnot(nrow(ER4_beta)>0L,!anyNA(ER4_index))
ER4_beta[ER4_index]<-ER4_results$effect_meta
ER4_sign[ER4_index]<-ER4_results$Sign
ER4_labels<-gsub("sv:","SV: ",gsub("_","-",sub("Oscillibacter sp. ER4:","",ER4_rows,fixed=TRUE)),fixed=TRUE)
ER4_limit<-max(abs(ER4_beta),na.rm=TRUE)

tiff(file="pics/F4B_Oscillibacter_ER4_sv_heatmap.tiff",width=3200,height=4500,res=300,compression="lzw")

par(oma=c(7,0,3,0))
heatmap.2(ER4_beta,scale="none",Rowv=FALSE,Colv=FALSE,dendrogram="none",revC=TRUE,
          col=colorRampPalette(c(mycolor2_blue_red[1],"white",mycolor2_blue_red[2]))(256),
          breaks=seq(-ER4_limit,ER4_limit,length.out=257),na.color=mycolor2_pink_grey[2],
          trace="none",density.info="none",labRow=ER4_labels,labCol=colnames(ER4_beta),
          cexRow=2.2,cexCol=2.1+2/par("ps"),srtCol=45,cellnote=ER4_sign,notecol="black",notecex=3.5,
          add.expr=mtext("Genomic position (kbp)",side=4,line=par("mgp")[2],at=nrow(ER4_beta)+0.9,
                         las=2,adj=0,cex=2.2*par("cex"),font=2,xpd=NA),
          key.title="Beta coefficient",
          key.xlab="",
          key.par=list(mar=c(4,4,5.93,1),cex.axis=2,cex.lab=2,cex.main=2),
          lmat=rbind(c(0,4,3),c(2,1,0)),lhei=c(0.5,3),lwid=c(0.4,4,1),margins=c(11,14))
mtext(expression(italic(Oscillibacter)~"sp. ER4"),side=3,outer=TRUE,line=0.5,cex=2)
mtext(c("* Normal P < 0.05; I\u00b2 \u2264 75%","\u2606 FDR < 0.05; I\u00b2 \u2264 75%","\u2605 Bonferroni P < 0.05; I\u00b2 \u2264 75%"),


      side=1,outer=TRUE,line=c(1,3,5),at=0.58,adj=0,cex=1.6+2/par("ps"))

dev.off()





### Selected vSV associations with BMI and age
setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")
source("functions.R")
library(ggplot2)
library(patchwork)

## Read clinical and vSV data from all participants.
all_basic <- read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv", header = TRUE)
all_basic$study_name <- normalize_study_name(all_basic$study_name)
# Include all health statuses; each plot requires callable vSV and non-missing phenotype.
all_basic$match_id <- row.names(all_basic)
load("01.cleanData/SV_all/vsgv.RData")
vsgv <- vsgv_sub
dir.create("pics/sites", recursive = TRUE, showWarnings = FALSE)

## Oscillibacter sp. ER4:1255_1257 and BMI

vsv_target <- "Oscillibacter sp. ER4:1255_1257"
vsv_pic_er4 <- vsgv[, vsv_target, drop = FALSE]
vsv_pic_er4$match_id <- row.names(vsv_pic_er4)
pic_vsv_er4 <- merge(all_basic, vsv_pic_er4, by = "match_id")
pic_vsv_er4$temp <- pic_vsv_er4[[vsv_target]]
pic_sub <- subset(pic_vsv_er4, !is.na(temp) & !is.na(BMI), drop = TRUE)
study_n <- sort(table(pic_sub$study_name), decreasing = TRUE)
study_order <- names(study_n[study_n > 100L])
if (!length(study_order)) stop("No cohorts have more than 100 complete samples.")
pic_sub <- subset(pic_sub, study_name %in% study_order, drop = TRUE)
study_colors <- setNames(mycolor11[seq_along(study_order)], study_order)
pic_sub$study_name <- factor(pic_sub$study_name, levels = study_order)

study_labels <- setNames(vapply(study_order, function(study) {
  paste(c(strwrap(gsub("_", " ", study, fixed = TRUE), width = 11),
          paste0("(n=", as.integer(study_n[[study]]), ")")),
        collapse = "\n")
}, character(1)), study_order)
plot_list <- lapply(study_order, function(study) {
  plot_data <- subset(pic_sub, study_name == study, drop = TRUE)
  x_limits <- range(plot_data$temp, na.rm = TRUE)
  bmi_limits <- range(plot_data$BMI, na.rm = TRUE)
  study_color <- study_colors[[study]]
  ggplot(plot_data, aes(x = temp, y = BMI)) +
    geom_point(color = study_color, alpha = 0.65, size = 0.8) +
    geom_smooth(method = "lm", color = study_color, fill = study_color,
                linewidth = 0.8) +
    scale_x_continuous(n.breaks = 3, guide = guide_axis(check.overlap = TRUE)) +
    coord_cartesian(xlim = x_limits, ylim = bmi_limits) +
    labs(title = study_labels[[study]], x = "Normalized region coverage",
         y = if (study == study_order[1L]) "BMI" else NULL) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(size = 16, face = "plain", hjust = 0.5,
                                    lineheight = 1.05, margin = margin(b = 8)),
          axis.text = element_text(size = 19), axis.title = element_text(size = 18),
          plot.margin = margin(12, 20, 8, 6))
})
p <- patchwork::wrap_plots(plot_list, nrow = 1) +
  patchwork::plot_layout(axis_titles = "collect_x")

tiff(file = "pics/F4C_Oscillibacter_ER4_BMI_1255_1257.tiff",
     width = 4608, height = 1800, res = 300, compression = "lzw")
print(p)
dev.off()


################################################################################
## Oscillibacter sp. ER4:2001_2003 and BMI

vsv_target <- "Oscillibacter sp. ER4:2001_2003"
vsv_pic_er4 <- vsgv[, vsv_target, drop = FALSE]
vsv_pic_er4$match_id <- row.names(vsv_pic_er4)
pic_vsv_er4 <- merge(all_basic, vsv_pic_er4, by = "match_id")
pic_vsv_er4$temp <- pic_vsv_er4[[vsv_target]]
pic_sub <- subset(pic_vsv_er4, !is.na(temp) & !is.na(BMI), drop = TRUE)
study_n <- sort(table(pic_sub$study_name), decreasing = TRUE)
study_order <- names(study_n[study_n > 100L])
if (!length(study_order)) stop("No cohorts have more than 100 complete samples.")
pic_sub <- subset(pic_sub, study_name %in% study_order, drop = TRUE)
study_colors <- setNames(mycolor11[seq_along(study_order)], study_order)
pic_sub$study_name <- factor(pic_sub$study_name, levels = study_order)

study_labels <- setNames(vapply(study_order, function(study) {
  paste(c(strwrap(gsub("_", " ", study, fixed = TRUE), width = 13),
          paste0("(n=", as.integer(study_n[[study]]), ")")),
        collapse = "\n")
}, character(1)), study_order)
plot_list <- lapply(study_order, function(study) {
  plot_data <- subset(pic_sub, study_name == study, drop = TRUE)
  x_limits <- range(plot_data$temp, na.rm = TRUE)
  bmi_limits <- range(plot_data$BMI, na.rm = TRUE)
  study_color <- study_colors[[study]]
  ggplot(plot_data, aes(x = temp, y = BMI)) +
    geom_point(color = study_color, alpha = 0.65, size = 0.8) +
    geom_smooth(method = "lm", color = study_color, fill = study_color,
                linewidth = 0.8) +
    scale_x_continuous(n.breaks = 3, guide = guide_axis(check.overlap = TRUE)) +
    coord_cartesian(xlim = x_limits, ylim = bmi_limits) +
    labs(title = study_labels[[study]], x = "Normalized region coverage",
         y = if (study == study_order[1L]) "BMI" else NULL) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(size = 18, face = "plain", hjust = 0.5,
                                    lineheight = 0.95),
          axis.text = element_text(size = 19), axis.title = element_text(size = 18),
          plot.margin = margin(5, 10, 5, 10))
})
p <- patchwork::wrap_plots(plot_list, nrow = 1) +
  patchwork::plot_layout(axis_titles = "collect_x")

tiff(file = "pics/F4D_Oscillibacter_ER4_BMI_2001_2003.tiff",
     width = 4608, height = 1800, res = 300, compression = "lzw")
print(p)
dev.off()

## Faecalibacterium cf. prausnitzii KLE1255:1368_1370 and age

vsv_age_target <- "Faecalibacterium cf. prausnitzii KLE1255:1368_1370"
vsv_age_data <- vsgv[, vsv_age_target, drop = FALSE]
vsv_age_data$match_id <- row.names(vsv_age_data)
pic_vsv_age <- merge(all_basic, vsv_age_data, by = "match_id")
pic_vsv_age$temp <- pic_vsv_age[[vsv_age_target]]
pic_sub <- subset(pic_vsv_age, !is.na(temp) & !is.na(age), drop = TRUE)
study_n <- sort(table(pic_sub$study_name), decreasing = TRUE)
study_order <- names(study_n[study_n > 100L])
if (!length(study_order)) stop("No cohorts have more than 100 complete samples.")
pic_sub <- subset(pic_sub, study_name %in% study_order, drop = TRUE)
study_colors <- setNames(mycolor11[seq_along(study_order)], study_order)
pic_sub$study_name <- factor(pic_sub$study_name, levels = study_order)

study_labels <- setNames(vapply(study_order, function(study) {
  paste(c(strwrap(gsub("_", " ", study, fixed = TRUE), width = 11),
          paste0("(n=", as.integer(study_n[[study]]), ")")),
        collapse = "\n")
}, character(1)), study_order)
x_limits <- range(pic_sub$temp, na.rm = TRUE)
plot_list <- lapply(study_order, function(study) {
  plot_data <- subset(pic_sub, study_name == study, drop = TRUE)
  age_limits <- range(plot_data$age, na.rm = TRUE)
  study_color <- study_colors[[study]]
  ggplot(plot_data, aes(x = temp, y = age)) +
    geom_point(color = study_color, alpha = 0.65, size = 0.8) +
    geom_smooth(method = "lm", color = study_color, fill = study_color,
                linewidth = 0.8) +
    scale_x_continuous(n.breaks = 3, guide = guide_axis(check.overlap = TRUE)) +
    coord_cartesian(xlim = x_limits, ylim = age_limits) +
    labs(title = study_labels[[study]], x = "Normalized region coverage",
         y = if (study == study_order[1L]) "Age (years)" else NULL) +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(size = 16, face = "plain", hjust = 0.5,
                                    lineheight = 1.05, margin = margin(b = 8)),
          axis.text = element_text(size = 19), axis.title = element_text(size = 18),
          plot.margin = margin(12, 20, 8, 6))
})
p <- patchwork::wrap_plots(plot_list, nrow = 1) +
  patchwork::plot_layout(axis_titles = "collect_x")
tiff(file = "pics/sites/F4_age_Fprausnitzii_vsv_1368_1370.tiff",
     width = 5760, height = 1800, res = 300, compression = "lzw")
print(p)
dev.off()



