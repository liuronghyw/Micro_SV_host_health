### HUMAnN3 SV gene associations and fixed-effect meta-analysis.
### LiuRong
### 2026-10-18

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

source("functions.R")
library(data.table)
library(foreach)
library(doParallel)
library(dplyr)


##########################################################
## match file 

mapping <- read.delim(file.path("00.rawData/humann3/match_file/SGV_Finder_genes_vs_uniref90_annotation.tsv"))

genes_model1 <- read.csv("07.Microbial_GWAS/models_sig_genes/SV_sig_genes_analysis_wide_model1.csv",header=T)
genes_model2 <- read.csv("07.Microbial_GWAS/models_sig_genes/SV_sig_genes_analysis_wide_model2.csv",header=T)

genes<-rbind(genes_model1,genes_model2)

result <- unique(merge(
  unique(genes[, c("genes", "type")]),
  unique(mapping[, c("gene", "UniRef90_id")]),
  by.x = "genes", by.y = "gene", all.x = TRUE
))

write.csv(result[, c("genes", "type", "UniRef90_id")],
          file.path("00.rawData/humann3/sig_SV/hm3_SV_genes_match_id.csv"),row.names = FALSE)




###################### 1.2 Inputs
all_basic<-read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv")

SVs_humann3<-fread("00.rawData/humann3/sig_SV/SV_no_species_merged_genefamilies_matched.tsv",header=T)
colnames(SVs_humann3) <- gsub("_kneaddata-RPKs", "", colnames(SVs_humann3))

SVs_humann3_t<-t(SVs_humann3[,-1])
colnames(SVs_humann3_t)<-SVs_humann3[,1]$"Gene Family"


#sub_sample<-data.frame(setdiff(all_basic$match_id, colnames(sig_SVs_shortbred)[-1]))
#colnames(sub_sample)<-"match_id"
#study_sample<-all_basic[,c("match_id","study_name")]
#study_sample_sub<-merge(study_sample,sub_sample,by.x="match_id",by.y="match_id")
#write.csv(study_sample_sub,"shortbred_sub.csv")



## 2 Associations between genes and prognosis
### 2.1 Preparation
if (!dir.exists("09.Microbial_humann3_sig_SVs")) {dir.create("09.Microbial_humann3_sig_SVs")}
if (!dir.exists("09.Microbial_humann3_sig_SVs/studies")) {dir.create("09.Microbial_humann3_sig_SVs/studies")}
if (!dir.exists("09.Microbial_humann3_sig_SVs/meta")) {dir.create("09.Microbial_humann3_sig_SVs/meta")}
if (!dir.exists("09.Microbial_humann3_sig_SVs/multi")) {dir.create("09.Microbial_humann3_sig_SVs/multi")}



############################################################################################
#######  index 
#######  2.3 linear regression model 1
## age  BMI fasting_insulin  hba1c  fasting_glucose 
## triglycerides hdl  ldl  cholesterol 
## dyastolic_p systolic_p prothrombin_time inr
## albumine creatinine  alt  ast  total_bilirubin

## Çø·ÖdsvºÍvSV

SV_sig_sites_genes<-read.csv("00.rawData/humann3/sig_SV/hm3_SV_genes_match_id.csv",header=T)
SV_sig_sites_genes_type<-SV_sig_sites_genes[,c("UniRef90_id","type")]

SV_sig_sites_genes_dsv<-subset(SV_sig_sites_genes_type,type=="dsv",drop=T)
SV_sig_sites_genes_vsv<-subset(SV_sig_sites_genes_type,type=="vsv",drop=T)

dsv_genes_to_keep <- unique(SV_sig_sites_genes_dsv$UniRef90_id)
vsv_genes_to_keep <- unique(SV_sig_sites_genes_vsv$UniRef90_id)


############# age all
age_basic<-subset(all_basic,is.na(age)==F ,drop=T)
table(age_basic$study_name)

age_basic<-subset(age_basic,study_name!="LiJ_2014",drop=T)
table(age_basic$study_name)

age_prog<-data.frame(age_basic[,"age_norm"])
row.names(age_prog)<-rownames(age_basic)
colnames(age_prog)<-"Y"

common_rows <- intersect(rownames(age_basic), rownames(SVs_humann3_t))
index_SVs_humann3 <- SVs_humann3_t[common_rows, ]

dsv_cols_keep <- grep(paste(dsv_genes_to_keep, collapse = "|"), 
                  colnames(index_SVs_humann3), 
                  value = TRUE)

vsv_cols_keep <- grep(paste(vsv_genes_to_keep, collapse = "|"), 
                  colnames(index_SVs_humann3), 
                  value = TRUE)

index_dSVs_humann3 <- index_SVs_humann3[, dsv_cols_keep]
index_vSVs_humann3 <- index_SVs_humann3[, vsv_cols_keep]

SVs_humann3_vsv <- linear_shortbred_vsv(age_prog, index_vSVs_humann3, age_basic,c("gender_code"), "age_norm")
write.csv(SVs_humann3_vsv, file = "09.Microbial_humann3_sig_SVs/studies/age_vsv.csv", row.names = FALSE)

SVs_humann3_dsv <- linear_shortbred_dsv(age_prog, index_dSVs_humann3, age_basic,c("gender_code"), "age_norm")
write.csv(SVs_humann3_dsv, file = "09.Microbial_humann3_sig_SVs/studies/age_dsv.csv", row.names = FALSE)



####################################################################################################
############### BMI all
BMI_basic<-subset(all_basic,is.na(BMI)==F,header=T)
table(BMI_basic$study_name)

BMI_basic<-subset(BMI_basic,study_name!="PiccinnoG_2025_cohort6",drop=T)

BMI_prog<-data.frame(BMI_basic[,"BMI_norm"])
row.names(BMI_prog)<-rownames(BMI_basic)
colnames(BMI_prog)<-"Y"

common_rows <- intersect(rownames(BMI_basic), rownames(SVs_humann3_t))
index_SVs_humann3 <- SVs_humann3_t[common_rows, ]

dsv_cols_keep <- grep(paste(dsv_genes_to_keep, collapse = "|"), 
                  colnames(index_SVs_humann3), 
                  value = TRUE)

vsv_cols_keep <- grep(paste(vsv_genes_to_keep, collapse = "|"), 
                  colnames(index_SVs_humann3), 
                  value = TRUE)

index_dSVs_humann3 <- index_SVs_humann3[, dsv_cols_keep]
index_vSVs_humann3 <- index_SVs_humann3[, vsv_cols_keep]

dSVs_humann3 <- linear_shortbred_dsv(BMI_prog, index_dSVs_humann3, BMI_basic,c("gender_code"), "BMI_norm")
write.csv(dSVs_humann3, file = "09.Microbial_humann3_sig_SVs/studies/BMI_dsv.csv", row.names = FALSE)

vSVs_humann3 <- linear_shortbred_vsv(BMI_prog, index_vSVs_humann3, BMI_basic,c("gender_code"), "BMI_norm")
write.csv(vSVs_humann3, file = "09.Microbial_humann3_sig_SVs/studies/BMI_vsv.csv", row.names = FALSE)



############################################################################################################

# Each cohort has beta, SE, P, lower CI, upper CI and N.
# Keep presence (dSV) and abundance (vSV) estimates as separate hypotheses.
# EE inverse-variance pooling matches combine_meta_p_lr in functions.R.
meta_linear <- function(dat) {
  if (names(dat)[1] != "X" || (ncol(dat)-1L) %% 6L != 0L || ncol(dat) < 7L)
    stop("Expected gene ID followed by six-column cohort blocks")
  offsets <- seq.int(2L, ncol(dat), by = 6L)
  for (j in offsets) {
    expected <- paste0(sub("\\.beta$", "", names(dat)[j]),
                       c(".beta", ".SE", ".p", ".beta_left", ".beta_right", ".N"))
    if (!identical(names(dat)[j + 0:5], expected))
      stop("Unexpected cohort columns at: ", names(dat)[j])
  }
  if (anyNA(dat$X) || anyDuplicated(dat$X)) stop("Missing or duplicate gene IDs")
  effects <- as.matrix(dat[, offsets, drop = FALSE])
  ses <- as.matrix(dat[, offsets+1L, drop = FALSE])
  ps <- as.matrix(dat[, offsets+2L, drop = FALSE])
  if (!is.numeric(effects) || !is.numeric(ses) || !is.numeric(ps))
    stop("Non-numeric study statistics")
  valid <- is.finite(effects) & is.finite(ses) & ses > 0
  k <- rowSums(valid)
  weights <- matrix(0, nrow(dat), length(offsets))
  weights[valid] <- 1/ses[valid]^2
  effects[!valid] <- 0
  sum_w <- rowSums(weights)
  beta <- rowSums(weights*effects)/sum_w
  se <- sqrt(1/sum_w)
  p <- 2*pnorm(abs(beta/se), lower.tail = FALSE)
  single <- which(k == 1L)
  for (i in single) {
    study_p <- ps[i, which(valid[i, ])]
    if (is.finite(study_p) && study_p >= 0 && study_p <= 1) p[i] <- study_p
  }
  q <- rowSums(weights*(effects-beta)^2)
  i2 <- rep(NA_real_, nrow(dat))
  multi <- which(k >= 2L)
  i2[multi] <- ifelse(q[multi] > 0,
                     pmax(0, (q[multi]-(k[multi]-1))/q[multi])*100, 0)
  het_p <- rep(NA_real_, nrow(dat))
  het_p[multi] <- pchisq(q[multi], k[multi]-1L, lower.tail = FALSE)
  het_p[single] <- 1
  beta[k == 0L] <- se[k == 0L] <- p[k == 0L] <- NA_real_
  p_count <- rowSums(valid & is.finite(ps) & ps >= 0 & ps <= 0.05, na.rm = TRUE)
  p_count[k == 0L] <- NA_integer_
  cbind(dat, effect_meta = beta, pvalue_meta = p, se_meta = se,
        I2 = i2, het_p = het_p, p_count = p_count, study_count = k,
        cilb = beta-qnorm(0.975)*se, ciub = beta+qnorm(0.975)*se)
}

result_root <- "09.Microbial_humann3_sig_SVs"
for (folder in c("meta", "multi"))
  dir.create(file.path(result_root, folder), recursive = TRUE, showWarnings = FALSE)
outcomes <- c("age", "BMI")
summary_rows <- list()
for (outcome in outcomes) {
  model_results <- list()
  for (sv_type in c("dsv", "vsv")) {
    input <- file.path(result_root, "studies", paste0(outcome, "_", sv_type, ".csv"))
    if (!file.exists(input)) stop("Missing study input: ", input)
    dat <- read.csv(input, check.names = FALSE, stringsAsFactors = FALSE)
    result <- meta_linear(dat)
    result$gene_model <- rep(if (sv_type == "dsv") "presence" else "abundance", nrow(result))
    model_results[[sv_type]] <- result
    write.csv(result, file.path(result_root, "meta", paste0(outcome, "_", sv_type, ".csv")),
              row.names = FALSE)
    summary_rows[[length(summary_rows)+1L]] <- data.frame(
      outcome = outcome, sv_type = sv_type, genes = nrow(result),
      tested = sum(is.finite(result$pvalue_meta)),
      multiple_studies = sum(result$study_count >= 2L),
      single_study = sum(result$study_count == 1L),
      no_valid_study = sum(result$study_count == 0L))
  }
  if (!identical(names(model_results$dsv), names(model_results$vsv)))
    stop("dSV and vSV cohort columns differ: ", outcome)
  result <- rbind(model_results$dsv, model_results$vsv)
  if (anyDuplicated(result[, c("X", "gene_model")])) stop("Duplicate gene-model pairs")
  write.csv(result, file.path(result_root, "meta", paste0(outcome, ".csv")), row.names = FALSE)
  # Correct over all tested gene-model pairs within each phenotype.
  tested <- which(is.finite(result$pvalue_meta))
  result$bonferroni <- result$fdr <- rep(NA_real_, nrow(result))
  result$bonferroni[tested] <- p.adjust(result$pvalue_meta[tested], method = "bonferroni")
  result$fdr[tested] <- p.adjust(result$pvalue_meta[tested], method = "BH")
  write.csv(result, file.path(result_root, "multi", paste0(outcome, ".csv")), row.names = FALSE)
  message("Completed meta-analysis: ", outcome, " (", nrow(result), " gene-model pairs)")
}
write.csv(do.call(rbind, summary_rows), file.path(result_root, "meta_summary.csv"),
          row.names = FALSE)



####################################################################################################
### HUMAnN3 plots: model1/model2 SV union
### LiuRong
# KEGG: significant SVs and gene abundance (P < 0.05 / associated SV genes), with concordant effects.
# At least one significant concordant UniRef90 result qualifies a gene; count once per direction.

library(data.table)
library(ggplot2)

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")
source("functions.R",encoding="UTF-8")
dir.create("pics",showWarnings=FALSE)



## 1 Preparation
### 1.1 SV union: model2 effects and P values take priority
model1<-as.data.frame(fread("07.Microbial_GWAS/models_sig_genes/SV_sig_genes_analysis_wide_model1.csv"))
model2<-as.data.frame(fread("07.Microbial_GWAS/models_sig_genes/SV_sig_genes_analysis_wide_model2.csv"))
model1$source_model<-"model1"
model2$source_model<-"model2"
models<-rbind(model2,model1)
sv_keys<-c("phenotype","SV")

sv_effect<-unique(models[,c(sv_keys,"effect_meta","pvalue_meta","source_model")])
stopifnot(!anyDuplicated(sv_effect[sv_effect$source_model=="model1",sv_keys]),
          !anyDuplicated(sv_effect[sv_effect$source_model=="model2",sv_keys]))
sv_effect<-sv_effect[!duplicated(sv_effect[,sv_keys]),]
sv_genes<-models[!duplicated(models[,c(sv_keys,"genes")]),
                 setdiff(names(models),c("effect_meta","pvalue_meta","source_model"))]
sv_results<-merge(sv_genes,sv_effect,by=sv_keys,all.x=TRUE)
names(sv_results)[names(sv_results)=="X"]<-"SV_short_name"
names(sv_results)[names(sv_results)=="effect_meta"]<-"effect_meta_SV"

### 1.2 Gene annotations and UniRef90 mapping
annotation<-as.data.frame(fread("00.rawData/gene_annotation/representatives.genes.drepped.annotations_renew.csv",
                               select=c(2,7,8)))
names(annotation)[names(annotation)=="gene_name"]<-"genes"
annotation<-unique(annotation[annotation$genes %in% sv_results$genes,])

uniref<-read.table("00.rawData/humann3/sig_SV/SGV_Finder_SV_genes_vs_uniref90.tsv",
                   header=FALSE,stringsAsFactors=FALSE)[,1:2]
names(uniref)<-c("genes","Uniref90_id")
uniref$Uniref90_id<-sub("\\|.*","",uniref$Uniref90_id)
uniref<-unique(uniref)

### 1.3 KEGG BRITE L2 categories, counted once per gene
kegg<-read.csv("00.rawData/KEGG/KEGG_classification.csv",stringsAsFactors=FALSE)
kegg<-unique(kegg[grepl("^091(00|20|30|40) ",kegg$L1),c("KO","L2")])
kegg$category<-sub("^[0-9]+ ","",kegg$L2)

orthology<-ifelse(is.na(annotation$Orthology),"",annotation$Orthology)
ko_list<-regmatches(orthology,gregexpr("K[0-9]{5}",orthology))
gene_ko<-unique(data.frame(genes=rep(annotation$genes,lengths(ko_list)),
                           KO=unlist(ko_list),stringsAsFactors=FALSE))
gene_class<-unique(merge(gene_ko,kegg,by="KO")[,c("genes","category")])
genes<-unique(sv_results$genes[!is.na(sv_results$genes) & nzchar(trimws(sv_results$genes))])
gene_function<-data.frame(genes=genes,category="Unknown / unannotated",stringsAsFactors=FALSE)
gene_function$category[genes %in% gene_ko$genes]<-"KO annotated, no functional class"
for(i in seq_len(nrow(gene_function))){
  categories<-unique(gene_class$category[gene_class$genes==gene_function$genes[i]])
  if(length(categories)==1) gene_function$category[i]<-categories
  if(length(categories)>1) gene_function$category[i]<-"Multiple KEGG subcategories"
}

category_levels<-c("Significant, same direction","Significant, opposite direction","Non-significant")
direction_levels<-c("Same direction with SV","Opposite direction with SV")
category_colors<-mycolor3_green_red_yellow
direction_colors<-setNames(category_colors[1:2],direction_levels)



## 2 BMI
### 2.1 HUMAnN3 gene effects
SV_BMI<-subset(sv_results,phenotype=="BMI_sv")
# Count unique nonempty genes in associated SVs before the HUMAnN3 join or significance filtering.
BMI_genes<-unique(SV_BMI$genes[!is.na(SV_BMI$genes) & nzchar(trimws(SV_BMI$genes))])
BMI_gene_n<-length(BMI_genes)
BMI_gene_p_cutoff<-0.05/BMI_gene_n
hm3_BMI<-as.data.frame(fread("09.Microbial_humann3_sig_SVs/meta/BMI.csv"))
names(hm3_BMI)[1]<-"Uniref90_id"
hm3_BMI<-hm3_BMI[,c("Uniref90_id","effect_meta","pvalue_meta")]
names(hm3_BMI)<-c("Uniref90_id","effect_meta_gene_abun","pvalue_meta_abun")
SV_gene_BMI<-unique(merge(SV_BMI,merge(uniref,hm3_BMI,by="Uniref90_id"),by="genes",all.x=TRUE))
SV_gene_BMI<-merge(SV_gene_BMI,annotation,by="genes",all.x=TRUE)

### 2.2 F5A: mutually exclusive gene P-value tiers
effect_SV<-SV_gene_BMI$effect_meta_SV
effect_gene_abun<-SV_gene_BMI$effect_meta_gene_abun
pvalue<-SV_gene_BMI$pvalue_meta_abun
x_value<-trimws(SV_gene_BMI$SV_short_name)
valid<-is.finite(effect_SV) & is.finite(effect_gene_abun) & is.finite(pvalue) &
       !is.na(x_value) & nzchar(x_value)
n_valid<-sum(valid)
n_excluded<-sum(!valid)
stopifnot(n_valid>0,all(pvalue[valid]>=0 & pvalue[valid]<=1))

# P-value tiers include both effect directions; each record belongs to one tier.
BMI_significance_levels<-c("Strong: P < 1e-7",
                            "Adjusted: 1e-7 <= P < 0.05 / N",
                            "Nominal: 0.05 / N <= P < 0.05",
                            "Non-significant: P >= 0.05")
BMI_significance_colors<-c("#006D2C","#41AB5D","#A1D99B","#ECD452")
BMI_significance_text_colors<-c("white","black","black","black")
category<-cut(pvalue[valid],breaks=c(0,1e-7,BMI_gene_p_cutoff,0.05,Inf),
               labels=BMI_significance_levels,right=FALSE)
all_sv<-unique(x_value[!is.na(x_value) & nzchar(x_value)])
counts<-as.data.frame.matrix(table(factor(x_value[valid],levels=all_sv),
                                    factor(category,levels=BMI_significance_levels)))
counts$Total<-rowSums(counts)
# Order species by their summed SV-gene records, then order SVs within each species.
counts$Species<-SV_gene_BMI$organism[match(rownames(counts),x_value)]
species_totals<-tapply(counts$Total,counts$Species,sum)
BMI_species_totals<-sort(species_totals,decreasing=TRUE)
counts<-counts[order(-species_totals[counts$Species],counts$Species,
                     -counts$Total,rownames(counts)),]
BMI_gene_counts<-counts
no_gene<-is.na(SV_gene_BMI$genes) | !nzchar(trimws(SV_gene_BMI$genes))
compact_labels<-paste0(trimws(SV_gene_BMI$Short_name),":",sub("^.*:","",x_value))
display_labels<-compact_labels[match(rownames(counts),x_value)]
bar_matrix<-t(as.matrix(counts[,BMI_significance_levels]))
totals<-rowSums(bar_matrix)
percentages<-100*totals/sum(totals)
legend_labels<-sprintf("%s\nrecords = %d (%.1f%%)",BMI_significance_levels,as.integer(totals),percentages)
bar_space<-c(0.25,ifelse(head(counts$Species,-1)==tail(counts$Species,-1),0.25,1.6))

figure_width<-max(210,ncol(bar_matrix)*2.5)/25.4
figure_height<-180/25.4


tiff(file = "pics/F5A.tiff", width =5000, height =2500, res =300, pointsize =12)

par(family = "sans", las = 1, xaxs = "i", yaxs = "i")
par(fig = c(0, 1, 0, 1), mai = c(3.3, 1.0, 0.5, 0.12))

bar_midpoints <- barplot(bar_matrix, col = BMI_significance_colors, border = NA, space = bar_space, axes = FALSE, axisnames = FALSE, ylim = c(0, max(counts$Total) * 1.12))

no_testable <- counts$Total == 0L
if (any(no_testable)) points(bar_midpoints[no_testable], rep(0, sum(no_testable)), pch = 4, col = "grey35", cex = 0.8, xpd = NA)

axis(2, las = 1, cex.axis = 0.8+3/12, lwd = 0.6, lwd.ticks = 0.6)
axis(1, at = bar_midpoints, labels = FALSE, tck = -0.008, lwd = 0.6, lwd.ticks = 0.6)
text(bar_midpoints, par("usr")[3] - max(counts$Total) * 0.025, labels = display_labels, srt = 90, adj = 1, xpd = NA, cex = 0.6+3/12)
mtext("Number of SV-gene records", side = 2, line = 2.9, las = 0, cex = 0.9+3/12, xpd = NA)

par(fig = c(0.48, 0.70, 0.63, 0.94), mai = rep(0, 4), new = TRUE)

nonzero <- totals > 0
pie(totals[nonzero], labels = rep("", sum(nonzero)), col = BMI_significance_colors[nonzero], border = "white", lwd = 0.5, clockwise = TRUE, init.angle = 90, radius = 0.85)

angles <- pi / 2 - 2 * pi * (cumsum(totals[nonzero]) - totals[nonzero] / 2) / sum(totals)
show_label <- percentages[nonzero] >= 5
text(0.63*cos(angles[show_label]),0.63*sin(angles[show_label]),
     sprintf("%.1f%%",percentages[nonzero][show_label]),
     col=BMI_significance_text_colors[nonzero][show_label],cex=0.8+3/12)

par(fig = c(0.70, 0.995, 0.60, 0.98), mai = rep(0, 4), new = TRUE)
plot.new()
legend("center", legend = legend_labels, fill = BMI_significance_colors, border = NA, bty = "n", cex = 0.7+3/12, y.intersp = 1.5,
       title = sprintf("Gene meta P-value tiers\nN = %d unique SV genes; 0.05/N = %.2g\nTotal = %d association records",BMI_gene_n,BMI_gene_p_cutoff,sum(totals)))

par(fig = c(0, 1, 0, 1), mai = rep(0, 4), new = TRUE)
plot.new()
#text(0.5, 0.04, paste0(
#  sprintf("BMI | %d SVs retained; %d testable gene-result records; %d unavailable records.\n", nrow(counts), n_valid, n_excluded),
#  sprintf("%d records without gene annotation; cross = no testable gene result.", sum(no_gene))), cex = 0.65)

dev.off()



### 2.3 F5B: all available SV and gene effects (no P-value filter)
# Retain one-to-many UniRef90 matches as separate result records.
plot_data<-subset(SV_gene_BMI,is.finite(effect_meta_SV) & is.finite(effect_meta_gene_abun) &
                   !is.na(genes) & nzchar(trimws(genes)))
stopifnot(nrow(plot_data)>0,all(plot_data$effect_meta_SV!=0 & plot_data$effect_meta_gene_abun!=0))
plot_data$effect_SV<-plot_data$effect_meta_SV
plot_data$effect_gene_abun<-plot_data$effect_meta_gene_abun
same_direction<-sign(plot_data$effect_SV)==sign(plot_data$effect_gene_abun)
plot_data$Direction<-factor(ifelse(same_direction,direction_levels[1],direction_levels[2]),levels=direction_levels)
BMI_scatter<-plot_data

p<-ggplot(plot_data,aes(x=effect_SV,y=effect_gene_abun,color=Direction))+
  geom_hline(yintercept=0,color="#BC3523",linetype="dotted",linewidth=0.5)+
  geom_vline(xintercept=0,color="#BC3523",linetype="dotted",linewidth=0.5)+
  geom_point(size=1.8,alpha=0.85)+
  scale_color_manual(values=direction_colors,breaks=direction_levels,drop=F)+
  scale_x_continuous(expand=expansion(mult=0.08))+
  scale_y_continuous(expand=expansion(mult=0.08))+
  expand_limits(x=0,y=0)+
  labs(x="Effect size of SVs",y="Effect size of bacterial genes",color=NULL,tag="")+
  guides(color=guide_legend(nrow=1,override.aes=list(size=3,alpha=1)))+
  theme_bw(base_size=14,base_family="sans")+
  theme(
    panel.grid.major=element_line(color="grey90",linewidth=0.35),
    panel.grid.minor=element_line(color="grey95",linewidth=0.25),
    panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
    axis.title=element_text(face="bold"),axis.text=element_text(color="black"),
    legend.position="top",legend.text=element_text(size=11.5),
    legend.key.size=grid::unit(3.5,"mm"),legend.key.spacing.x=grid::unit(4,"mm"),
    plot.tag=element_text(face="bold",size=18),plot.tag.position=c(0.01,0.99),
    plot.margin=margin(12,10,8,8),aspect.ratio=1.05
  )

tiff("pics/F5B.tif",width=5.2,height=5.3,units="in",res=300,
     compression="lzw",type="cairo",pointsize=10)
print(p)
dev.off()




### 2.4 F5C: KEGG annotation of positive concordant genes
BMI_concordant<-subset(BMI_scatter,is.finite(pvalue_meta_abun) & pvalue_meta_abun<BMI_gene_p_cutoff &
                         Direction==direction_levels[1])
# Count each selected gene once per effect direction.
foreground_genes<-unique(subset(BMI_concordant,effect_meta_SV>0 & !is.na(genes) & nzchar(trimws(genes)))$genes)
gene_category<-gene_function[gene_function$genes %in% foreground_genes,]
plot_counts<-sort(table(gene_category$category),decreasing=TRUE)
BMI_positive_KEGG<-data.frame(category=names(plot_counts),genes=as.integer(plot_counts),
                                 percentage=100*as.integer(plot_counts)/length(foreground_genes))
stopifnot(length(foreground_genes)>0,sum(plot_counts)==length(foreground_genes))
colors<-mycolor15[seq_along(plot_counts)]

tiff("pics/F5C.tif",width=3.6,height=5,units="in",res=300,
     compression="lzw",type="cairo",pointsize=10)

par(fig=c(0.08,0.92,0.40,0.92),mar=c(0,0,0,0),family="sans")
plot.new()
plot.window(xlim=c(-1,1),ylim=c(-1,1),asp=1)
boundaries<-pi/2-2*pi*c(0,cumsum(plot_counts))/sum(plot_counts)

for(j in seq_along(plot_counts)){
  arc<-seq(boundaries[j],boundaries[j+1],length.out=150)
  polygon(c(0.9*cos(arc),0.37*cos(rev(arc))),c(0.9*sin(arc),0.37*sin(rev(arc))),col=colors[j],border="white",lwd=0.6)
}

text(0,0,paste0(length(foreground_genes)," genes"),cex=1.15,font=2)
angles<-pi/2-2*pi*(cumsum(plot_counts)-plot_counts/2)/sum(plot_counts)
percentages<-100*plot_counts/sum(plot_counts)
label_index<-percentages>=8
text(0.66*cos(angles[label_index]),0.66*sin(angles[label_index]),sprintf("%.1f%%",percentages[label_index]),cex=1)

par(fig=c(0.02,0.98,0.09,0.40),mar=c(0,0,0,0),new=T)
plot.new()
legend_labels<-paste0(names(plot_counts)," (n = ",as.integer(plot_counts),")")

for(i in seq_along(legend_labels)){
  legend_labels[i]<-paste(strwrap(legend_labels[i],width=48),collapse="\n")
}

legend("center",legend=legend_labels,fill=colors,border=NA,bty="n",cex=1,y.intersp=1.1)

par(fig=c(0,1,0,1),mar=c(0,0,0,0),new=T)
plot.new()
text(0.5,0.96,"BMI: SV beta > 0\nKEGG BRITE functional subcategories",font=2,cex=1.1)
dev.off()

### 2.5 F5D: KEGG annotation of negative concordant genes
# Count each selected gene once per effect direction.
foreground_genes<-unique(subset(BMI_concordant,effect_meta_SV<0 & !is.na(genes) & nzchar(trimws(genes)))$genes)
gene_category<-gene_function[gene_function$genes %in% foreground_genes,]
plot_counts<-sort(table(gene_category$category),decreasing=TRUE)
BMI_negative_KEGG<-data.frame(category=names(plot_counts),genes=as.integer(plot_counts),
                                 percentage=100*as.integer(plot_counts)/length(foreground_genes))
stopifnot(length(foreground_genes)>0,sum(plot_counts)==length(foreground_genes))
colors<-mycolor15[seq_along(plot_counts)]

tiff("pics/F5D.tif",width=3.6,height=5,units="in",res=300,
     compression="lzw",type="cairo",pointsize=10)

par(fig=c(0.08,0.92,0.40,0.92),mar=c(0,0,0,0),family="sans")
plot.new()
plot.window(xlim=c(-1,1),ylim=c(-1,1),asp=1)
boundaries<-pi/2-2*pi*c(0,cumsum(plot_counts))/sum(plot_counts)

for(j in seq_along(plot_counts)){
  arc<-seq(boundaries[j],boundaries[j+1],length.out=150)
  polygon(c(0.9*cos(arc),0.37*cos(rev(arc))),c(0.9*sin(arc),0.37*sin(rev(arc))),col=colors[j],border="white",lwd=0.6)
}

text(0,0,paste0(length(foreground_genes)," genes"),cex=1.15,font=2)
angles<-pi/2-2*pi*(cumsum(plot_counts)-plot_counts/2)/sum(plot_counts)
percentages<-100*plot_counts/sum(plot_counts)
label_index<-percentages>=8
text(0.66*cos(angles[label_index]),0.66*sin(angles[label_index]),sprintf("%.1f%%",percentages[label_index]),cex=1)

par(fig=c(0.02,0.98,0.09,0.40),mar=c(0,0,0,0),new=T)
plot.new()
legend_labels<-paste0(names(plot_counts)," (n = ",as.integer(plot_counts),")")

for(i in seq_along(legend_labels)){
  legend_labels[i]<-paste(strwrap(legend_labels[i],width=48),collapse="\n")
}

legend("center",legend=legend_labels,fill=colors,border=NA,bty="n",cex=1,y.intersp=1.1)

par(fig=c(0,1,0,1),mar=c(0,0,0,0),new=T)
plot.new()
text(0.5,0.96,"BMI: SV beta < 0\nKEGG BRITE functional subcategories",font=2,cex=1.1)
dev.off()




## 3 age
### 3.1 HUMAnN3 gene effects
SV_age<-subset(sv_results,phenotype=="age_sv")
# Count unique nonempty genes in associated SVs before the HUMAnN3 join or significance filtering.
age_genes<-unique(SV_age$genes[!is.na(SV_age$genes) & nzchar(trimws(SV_age$genes))])
age_gene_n<-length(age_genes)
age_gene_p_cutoff<-0.05/age_gene_n
hm3_age<-as.data.frame(fread("09.Microbial_humann3_sig_SVs/meta/age.csv"))
names(hm3_age)[1]<-"Uniref90_id"
hm3_age<-hm3_age[,c("Uniref90_id","effect_meta","pvalue_meta")]
names(hm3_age)<-c("Uniref90_id","effect_meta_gene_abun","pvalue_meta_abun")
SV_gene_age<-unique(merge(SV_age,merge(uniref,hm3_age,by="Uniref90_id"),by="genes",all.x=TRUE))
SV_gene_age<-merge(SV_gene_age,annotation,by="genes",all.x=TRUE)

### 3.2 F5_age_A: association categories (gene P < 0.05 / associated SV genes)
effect_SV<-SV_gene_age$effect_meta_SV
effect_gene_abun<-SV_gene_age$effect_meta_gene_abun
pvalue<-SV_gene_age$pvalue_meta_abun
x_value<-trimws(SV_gene_age$SV_short_name)
valid<-is.finite(effect_SV) & is.finite(effect_gene_abun) & is.finite(pvalue) &
       !is.na(x_value) & nzchar(x_value)
n_valid<-sum(valid)
n_excluded<-sum(!valid)
stopifnot(n_valid>0,all(pvalue[valid]>=0 & pvalue[valid]<=1))

p_cutoff<-age_gene_p_cutoff
significant<-pvalue[valid]<p_cutoff
same_direction<-sign(effect_SV[valid])==sign(effect_gene_abun[valid])
stopifnot(all(effect_SV[valid][significant]!=0 & effect_gene_abun[valid][significant]!=0))
category<-ifelse(!significant,category_levels[3],
                 ifelse(same_direction,category_levels[1],category_levels[2]))
all_sv<-unique(x_value[!is.na(x_value) & nzchar(x_value)])
counts<-as.data.frame.matrix(table(factor(x_value[valid],levels=all_sv),
                                    factor(category,levels=category_levels)))
counts$Total<-rowSums(counts)
counts<-counts[order(-counts$Total,rownames(counts)),]
age_gene_counts<-counts
no_gene<-is.na(SV_gene_age$genes) | !nzchar(trimws(SV_gene_age$genes))
compact_labels<-paste0(trimws(SV_gene_age$Short_name),":",sub("^.*:","",x_value))
display_labels<-compact_labels[match(rownames(counts),x_value)]
bar_matrix<-t(as.matrix(counts[,category_levels]))
totals<-rowSums(bar_matrix)
percentages<-100*totals/sum(totals)
legend_labels<-sprintf("%s: n = %d (%.1f%%)",category_levels,as.integer(totals),percentages)
figure_width<-max(210,ncol(bar_matrix)*2.5)/25.4
figure_height<-180/25.4

tiff("pics/F5_age_A.tif",width=figure_width,height=figure_height,units="in",res=300,
     compression="lzw",type="cairo",pointsize=10)

par(family = "sans", las = 1, xaxs = "i", yaxs = "i")
par(fig = c(0, 1, 0, 1), mai = c(2.25, 0.65, 0.5, 0.12))

bar_midpoints <- barplot(bar_matrix, col = category_colors, border = NA, space = 0.25, axes = FALSE, axisnames = FALSE, ylim = c(0, max(counts$Total) * 1.12))

no_testable <- counts$Total == 0L
if (any(no_testable)) points(bar_midpoints[no_testable], rep(0, sum(no_testable)), pch = 4, col = "grey35", cex = 0.8, xpd = NA)

axis(2, las = 1, cex.axis = 0.8, lwd = 0.6, lwd.ticks = 0.6)
axis(1, at = bar_midpoints, labels = FALSE, tck = -0.008, lwd = 0.6, lwd.ticks = 0.6)
text(bar_midpoints, par("usr")[3] - max(counts$Total) * 0.025, labels = display_labels, srt = 90, adj = 1, xpd = NA, cex = 0.6)
mtext("Number of SV-gene records", side = 2, line = 2.9, las = 0, cex = 0.9)
mtext("A", side = 3, adj = 0, line = 1.1, font = 2, cex = 1.2)

par(fig = c(0.48, 0.70, 0.63, 0.94), mai = rep(0, 4), new = TRUE)

nonzero <- totals > 0
pie(totals[nonzero], labels = rep("", sum(nonzero)), col = category_colors[nonzero], border = "white", lwd = 0.5, clockwise = TRUE, init.angle = 90, radius = 0.85)

angles <- pi / 2 - 2 * pi * (cumsum(totals[nonzero]) - totals[nonzero] / 2) / sum(totals)
show_label <- percentages[nonzero] >= 5
text(0.53 * cos(angles[show_label]), 0.53 * sin(angles[show_label]), sprintf("%.1f%%", percentages[nonzero][show_label]), cex = 0.8)

par(fig = c(0.70, 0.995, 0.64, 0.94), mai = rep(0, 4), new = TRUE)
plot.new()
legend("center", legend = legend_labels, fill = category_colors, border = NA, bty = "n", cex = 0.7, y.intersp = 1.5, title = sprintf("Gene meta P < %.2g\n(0.05 / %d SV genes)",p_cutoff,age_gene_n))

par(fig = c(0, 1, 0, 1), mai = rep(0, 4), new = TRUE)
plot.new()
text(0.5, 0.04, paste0(
  sprintf("age | %d SVs retained; %d testable gene-result records; %d unavailable records.\n", nrow(counts), n_valid, n_excluded),
  sprintf("%d records without gene annotation; cross = no testable gene result.", sum(no_gene))), cex = 0.65)

dev.off()

### 3.3 F5_age_B: SV and gene effects (gene P < 0.05 / associated SV genes)
# Retain one-to-many UniRef90 matches as separate result records.
p_cutoff<-age_gene_p_cutoff
plot_data<-SV_gene_age[valid & pvalue<p_cutoff,,drop=FALSE]
stopifnot(nrow(plot_data)>0,all(plot_data$effect_meta_SV!=0 & plot_data$effect_meta_gene_abun!=0))
plot_data$effect_SV<-plot_data$effect_meta_SV
plot_data$effect_gene_abun<-plot_data$effect_meta_gene_abun
same_direction<-sign(plot_data$effect_SV)==sign(plot_data$effect_gene_abun)
plot_data$Direction<-factor(ifelse(same_direction,direction_levels[1],direction_levels[2]),levels=direction_levels)
age_scatter<-plot_data
figure_note<-sprintf("age | %d SV-gene-UniRef records; %d unique genes;\ngene meta P < %.2g (0.05 / %d SV genes)",
                     nrow(plot_data),length(unique(plot_data$genes)),p_cutoff,age_gene_n)

p<-ggplot(plot_data,aes(x=effect_SV,y=effect_gene_abun,color=Direction))+
  geom_hline(yintercept=0,color="#BC3523",linetype="dotted",linewidth=0.5)+
  geom_vline(xintercept=0,color="#BC3523",linetype="dotted",linewidth=0.5)+
  geom_point(size=1.8,alpha=0.85)+
  scale_color_manual(values=direction_colors,breaks=direction_levels,drop=F)+
  scale_x_continuous(expand=expansion(mult=0.08))+
  scale_y_continuous(expand=expansion(mult=0.08))+
  expand_limits(x=0,y=0)+
  labs(x="Effect size of SVs",y="Effect size of bacterial genes",color=NULL,tag="B",caption=figure_note)+
  guides(color=guide_legend(nrow=1,override.aes=list(size=3,alpha=1)))+
  theme_bw(base_size=11,base_family="sans")+
  theme(
    panel.grid.major=element_line(color="grey90",linewidth=0.35),
    panel.grid.minor=element_line(color="grey95",linewidth=0.25),
    panel.border=element_rect(color="black",fill=NA,linewidth=0.5),
    axis.title=element_text(face="bold"),axis.text=element_text(color="black"),
    legend.position="top",legend.text=element_text(size=8.5),
    legend.key.size=grid::unit(3.5,"mm"),legend.key.spacing.x=grid::unit(4,"mm"),
    plot.tag=element_text(face="bold",size=15),plot.tag.position=c(0.01,0.99),
    plot.caption=element_text(size=8,hjust=0),plot.margin=margin(12,10,8,8),aspect.ratio=1.05
  )

tiff("pics/F5_age_B.tif",width=5.2,height=5.3,units="in",res=300,
     compression="lzw",type="cairo",pointsize=10)
print(p)
dev.off()

### 3.4 F5_age_C: KEGG annotation of positive concordant genes
age_concordant<-subset(age_scatter,Direction==direction_levels[1])
# Count each selected gene once per effect direction.
foreground_genes<-unique(subset(age_concordant,effect_meta_SV>0 & !is.na(genes) & nzchar(trimws(genes)))$genes)
gene_category<-gene_function[gene_function$genes %in% foreground_genes,]
plot_counts<-sort(table(gene_category$category),decreasing=TRUE)
age_positive_KEGG<-data.frame(category=names(plot_counts),genes=as.integer(plot_counts),
                                 percentage=100*as.integer(plot_counts)/length(foreground_genes))
stopifnot(sum(plot_counts)==length(foreground_genes))
colors<-mycolor15[seq_along(plot_counts)]

tiff("pics/F5_age_C.tif",width=5,height=5,units="in",res=300,
     compression="lzw",type="cairo",pointsize=10)

if(length(foreground_genes)>0){
  par(fig=c(0.08,0.92,0.40,0.92),mar=c(0,0,0,0),family="sans")
  plot.new()
  plot.window(xlim=c(-1,1),ylim=c(-1,1),asp=1)
  boundaries<-pi/2-2*pi*c(0,cumsum(plot_counts))/sum(plot_counts)

  for(j in seq_along(plot_counts)){
    arc<-seq(boundaries[j],boundaries[j+1],length.out=150)
    polygon(c(0.9*cos(arc),0.37*cos(rev(arc))),c(0.9*sin(arc),0.37*sin(rev(arc))),col=colors[j],border="white",lwd=0.6)
  }

  text(0,0,paste0(length(foreground_genes)," genes"),cex=0.95,font=2)
  angles<-pi/2-2*pi*(cumsum(plot_counts)-plot_counts/2)/sum(plot_counts)
  percentages<-100*plot_counts/sum(plot_counts)
  label_index<-percentages>=8
  text(0.66*cos(angles[label_index]),0.66*sin(angles[label_index]),sprintf("%.1f%%",percentages[label_index]),cex=0.8)

  par(fig=c(0.02,0.98,0.09,0.40),mar=c(0,0,0,0),new=T)
  plot.new()
  legend_labels<-paste0(names(plot_counts)," (n = ",as.integer(plot_counts),")")

  for(i in seq_along(legend_labels)){
    legend_labels[i]<-paste(strwrap(legend_labels[i],width=48),collapse="\n")
  }

  legend("center",legend=legend_labels,fill=colors,border=NA,bty="n",cex=0.8,y.intersp=1.1)

  par(fig=c(0,1,0,1),mar=c(0,0,0,0),new=T)
  plot.new()
  text(0.5,0.96,"age: SV beta > 0\nKEGG BRITE functional subcategories",font=2,cex=0.9)
  text(0.5,0.035,sprintf("Concordant effects; gene P < 0.05/%d (%.2g).",age_gene_n,age_gene_p_cutoff),cex=0.8)
} else {
  par(mar=c(0,0,0,0),family="sans")
  plot.new()
  text(0.5,0.85,"age: SV beta > 0\nKEGG BRITE functional subcategories",font=2,cex=0.9)
  text(0.5,0.5,sprintf("0 genes\nNo concordant positive genes with gene P < %.2g.",age_gene_p_cutoff),cex=0.9)
}
dev.off()

### 3.5 F5_age_D: KEGG annotation of negative concordant genes
# Count each selected gene once per effect direction.
foreground_genes<-unique(subset(age_concordant,effect_meta_SV<0 & !is.na(genes) & nzchar(trimws(genes)))$genes)
gene_category<-gene_function[gene_function$genes %in% foreground_genes,]
plot_counts<-sort(table(gene_category$category),decreasing=TRUE)
age_negative_KEGG<-data.frame(category=names(plot_counts),genes=as.integer(plot_counts),
                                 percentage=100*as.integer(plot_counts)/length(foreground_genes))
stopifnot(length(foreground_genes)>0,sum(plot_counts)==length(foreground_genes))
colors<-mycolor15[seq_along(plot_counts)]

tiff("pics/F5_age_D.tif",width=5,height=5,units="in",res=300,
     compression="lzw",type="cairo",pointsize=10)

par(fig=c(0.08,0.92,0.40,0.92),mar=c(0,0,0,0),family="sans")
plot.new()
plot.window(xlim=c(-1,1),ylim=c(-1,1),asp=1)
boundaries<-pi/2-2*pi*c(0,cumsum(plot_counts))/sum(plot_counts)

for(j in seq_along(plot_counts)){
  arc<-seq(boundaries[j],boundaries[j+1],length.out=150)
  polygon(c(0.9*cos(arc),0.37*cos(rev(arc))),c(0.9*sin(arc),0.37*sin(rev(arc))),col=colors[j],border="white",lwd=0.6)
}

text(0,0,paste0(length(foreground_genes)," genes"),cex=0.95,font=2)
angles<-pi/2-2*pi*(cumsum(plot_counts)-plot_counts/2)/sum(plot_counts)
percentages<-100*plot_counts/sum(plot_counts)
label_index<-percentages>=8
text(0.66*cos(angles[label_index]),0.66*sin(angles[label_index]),sprintf("%.1f%%",percentages[label_index]),cex=0.8)

par(fig=c(0.02,0.98,0.09,0.40),mar=c(0,0,0,0),new=T)
plot.new()
legend_labels<-paste0(names(plot_counts)," (n = ",as.integer(plot_counts),")")

for(i in seq_along(legend_labels)){
  legend_labels[i]<-paste(strwrap(legend_labels[i],width=48),collapse="\n")
}

legend("center",legend=legend_labels,fill=colors,border=NA,bty="n",cex=0.8,y.intersp=1.1)

par(fig=c(0,1,0,1),mar=c(0,0,0,0),new=T)
plot.new()
text(0.5,0.96,"age: SV beta < 0\nKEGG BRITE functional subcategories",font=2,cex=0.9)
text(0.5,0.035,sprintf("Concordant effects; gene P < 0.05/%d (%.2g).",age_gene_n,age_gene_p_cutoff),cex=0.8)
dev.off()




