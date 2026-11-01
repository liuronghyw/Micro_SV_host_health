### SV information
### 2026-10-10
### LiuRong

library("ggplot2")
library("reshape2")	
library("RColorBrewer")
library("dplyr")
library("tidyr")

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

stopifnot(nzchar(Sys.setlocale("LC_CTYPE", "Chinese_China.936")))
source("functions.R", encoding = "UTF-8")

if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
}

technical_dir <- "02.SV_summary_statistics/technical_effects"
dir.create(technical_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("pics", recursive = TRUE, showWarnings = FALSE)

# Keep enlarged study labels within their facet strips.
fs2_study_labels <- function(x) {
  ggplot2::label_wrap_gen(width = 16)(gsub("_", " ", x))
}

fs2_annotation_labels <- function(x) {
  x <- gsub(", p", "\np", x, fixed = TRUE)
  gsub("p = 0e+00", "p < 2.2e-16", x, fixed = TRUE)
}



## 1 Preparation
### 1.1 Import
### 1.2 Inputs
# SV
#all_dsv <- read.table("01.cleanData/SV_all/dsgv_all.tsv",check.names = F)
#all_vsv <- read.table("01.cleanData/SV_all/vsgv_all.tsv",check.names = F)

load("01.cleanData/SV_all/dsgv.RData")
load("01.cleanData/SV_all/vsgv.RData")

all_vsv<-vsgv_sub
all_dsv<-dsgv_sub

info    <- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",
                      sep = "\t", header = T, stringsAsFactors = F)

# abundance
all_abun_sv<-read.table("01.cleanData/mbio_all/SV_species_abun.tsv",check.names = F)

# Basic
all_basic <- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

read_length<-read.table("00.rawData/read_length_count/whole_dataset_read_length_count.txt",header=T)
colnames(read_length)[1]<-"match_id"
read_length$source_row <- seq_len(nrow(read_length))
read_length_fields <- c("Read_count", "Min_read_length", "Max_read_length",
                        "Median_read_length", "Average_read_length")
stopifnot(all(read_length_fields %in% colnames(read_length)))
# Remove only complete repeated header rows introduced by concatenation.
read_length_header <- as.character(read_length$match_id) == "Sample"
for (measure in read_length_fields) {
  read_length_header <- read_length_header & as.character(read_length[[measure]]) == measure
}
read_length_header[is.na(read_length_header)] <- FALSE
write.csv(read_length[read_length_header, , drop = FALSE],
          file.path(technical_dir, "read_length_headers_removed.csv"), row.names = FALSE)
read_length_header_n <- sum(read_length_header)
read_length <- read_length[!read_length_header, , drop = FALSE]
for (measure in read_length_fields) {
  original_value <- read_length[[measure]]
  numeric_value <- as.numeric(as.character(original_value))
  if (any(!is.na(original_value) & !is.finite(numeric_value))) {
    stop("Invalid read-length metadata in ", measure)
  }
  read_length[[measure]] <- numeric_value
}


# A zero placeholder is removed only when the same ID has a nonzero record.
read_length_zero <- rowSums(!is.na(read_length[, read_length_fields]) &
                              read_length[, read_length_fields] == 0) == length(read_length_fields)
read_length_duplicate <- duplicated(read_length$match_id) | duplicated(read_length$match_id, fromLast = TRUE)
read_length_remove <- read_length_zero & read_length_duplicate &
  read_length$match_id %in% read_length$match_id[!read_length_zero]
write.csv(read_length[read_length_remove, , drop = FALSE],
          file.path(technical_dir, "read_length_zero_placeholders_removed.csv"), row.names = FALSE)
read_length_zero_n <- sum(read_length_remove)
read_length <- read_length[!read_length_remove, , drop = FALSE]
if (anyDuplicated(read_length$match_id)) stop("Conflicting read-length records remain for a sample ID")
clinical_ids_before <- all_basic$match_id
stopifnot(!anyDuplicated(clinical_ids_before))
all_basic<-merge(all_basic,read_length,by.x="match_id",by.y="match_id",all.x=T)
stopifnot(nrow(all_basic) == length(clinical_ids_before),
          !anyDuplicated(all_basic$match_id), setequal(all_basic$match_id, clinical_ids_before))
read_length_missing <- all_basic[
  !is.finite(all_basic$Median_read_length) | all_basic$Median_read_length <= 0 |
    !is.finite(all_basic$Average_read_length) | all_basic$Average_read_length <= 0,
  c("match_id", "study_name", "Read_count", "Median_read_length", "Average_read_length"), drop = FALSE]
write.csv(read_length_missing, file.path(technical_dir, "read_length_missing_samples.csv"), row.names = FALSE)
technical_sample_counts <- data.frame(
  metric = c("clinical_rows_before", "clinical_rows_after", "clinical_unique_ids_after",
             "read_length_header_rows_removed", "read_length_zero_placeholder_rows_removed",
             "external_read_length_rows_after", "median_read_length_valid_samples",
             "average_read_length_valid_samples"),
  value = c(length(clinical_ids_before), nrow(all_basic), length(unique(all_basic$match_id)),
            read_length_header_n, read_length_zero_n, nrow(read_length),
            sum(is.finite(all_basic$Median_read_length) & all_basic$Median_read_length > 0),
            sum(is.finite(all_basic$Average_read_length) & all_basic$Average_read_length > 0)))
write.csv(technical_sample_counts, file.path(technical_dir, "technical_sample_counts.csv"), row.names = FALSE)

#test<-subset(all_basic,is.na(Median_read_length)==T,drop=T)
#write.csv(test,"test.csv")




###########################################################################################################
### analyze_technical_effects_on_sv

samp_dsgv_num<-data.frame(apply(all_dsv, 1, myfun<-function(x){sum(!is.na(x))}))
samp_vsgv_num<-data.frame(apply(all_vsv, 1, myfun<-function(x){sum(!is.na(x))}))

samp_dsgv_num$match_id<-rownames(samp_dsgv_num)
samp_vsgv_num$match_id<-rownames(samp_vsgv_num)

samp_sv_num<-merge(samp_dsgv_num,samp_vsgv_num,by.x="match_id",by.y="match_id")
colnames(samp_sv_num)[c(2,3)]<-c("dsv_number","vsv_number")

samp_sv_num$sv_number<-samp_sv_num$dsv_number+samp_sv_num$vsv_number

technique_test<-merge(all_basic,samp_sv_num,by.x="match_id",by.y="match_id")

stopifnot(!anyDuplicated(technique_test$match_id), nrow(technique_test) == length(clinical_ids_before))

technique_test<-technique_test[,c("match_id","study_name","read_count","sequencing_platform",
                                  "DNA_extraction_kit","median_read_length",
                                  "Read_count","Median_read_length","Average_read_length","sv_number","dsv_number","vsv_number")]

#a<-table(technique_test$study_name,technique_test$sequencing_platform)
#write.csv(a,"test.csv")

table(technique_test$sequencing_platform)
b<-table(technique_test$DNA_extraction_kit)
c<-table(technique_test$median_read_length)


##################################################################################################
### difference between sequencing platform

technique_platform<-technique_test

technique_platform$sequence_platform_rename <- technique_platform$sequencing_platform

technique_platform$sequence_platform_rename[technique_platform$sequencing_platform == "BGISEQ-500"] <- "BIGSEG"

technique_platform$sequence_platform_rename[
  technique_platform$sequencing_platform %in% c(
    "illimina HiSeq 2500", "illumina", "Illumina", "Illumina NovaSeq 6000",
    "Illumina NovaSeq6000", "Illumina HiSeq", "IlluminaMiSeq", "IlluminaHiSeq",
    "IlluminaNextSeq", "IlluminaNovaSeq"
  )
] <- "Illumina"


#technique_platform_platform<-subset(technique_platform,is.na(sequence_platform_rename)==F,drop=T)
technique_platform_platform<-subset(technique_platform,sequence_platform_rename!="",drop=T)

# sample size
sample_size = technique_platform_platform %>% group_by(sequence_platform_rename) %>% dplyr::summarize(num=n())

table(technique_platform$sequence_platform_rename)

# dsv_number
aov_dsv <- aov(dsv_number ~ sequence_platform_rename, data = technique_platform)
summary(aov_dsv)

# vsv_number
aov_vsv <- aov(vsv_number ~ sequence_platform_rename, data = technique_platform)
summary(aov_vsv)

# sv_number
aov_sv <- aov(sv_number ~ sequence_platform_rename, data = technique_platform)
summary(aov_sv)


table(technique_platform$sequence_platform_rename)

p_dsv <- summary(aov_dsv)[[1]]$`Pr(>F)`[1]
p_vsv <- summary(aov_vsv)[[1]]$`Pr(>F)`[1]
p_sv  <- summary(aov_sv)[[1]]$`Pr(>F)`[1]

annot_df <- data.frame(
  sv_type = factor(c("Total SV", "dSV", "vSV"),
                   levels = c("Total SV", "dSV", "vSV")),
  label = c(
    paste0("ANOVA p = ", format(p_sv, scientific = TRUE, digits = 2)),
    paste0("ANOVA p = ", format(p_dsv, scientific = TRUE, digits = 2)),
    paste0("ANOVA p = ", format(p_vsv, scientific = TRUE, digits = 2))
  )
)


library(tidyr)

sv_long <- technique_platform %>%
  dplyr::select(match_id, sequence_platform_rename, 
                sv_number, dsv_number, vsv_number) %>%
  pivot_longer(
    cols = c(sv_number, dsv_number, vsv_number),
    names_to = "sv_type",
    values_to = "value"
  ) %>%
  mutate(
    sv_type = factor(sv_type, 
                     levels = c("sv_number", "dsv_number", "vsv_number"),
                     labels = c("Total SV", "dSV", "vSV"))
  )


sample_size2 <- technique_platform %>%
  group_by(sequence_platform_rename) %>%
  dplyr::summarise(num = n())

sv_long <- sv_long %>%
  left_join(sample_size2, by = "sequence_platform_rename") %>%
  mutate(myaxis = paste0(sequence_platform_rename, "\n", "n=", num))

annot_df$label <- fs2_annotation_labels(annot_df$label)
p_multi <- ggplot(sv_long, aes(x = myaxis, y = value, fill = sequence_platform_rename)) +
  geom_violin(width = 0.8, alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.25, color = "black", alpha = 0.3, outlier.shape = NA) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  facet_wrap(~ sv_type, scales = "free_y", ncol = 3) +
  geom_text(data = annot_df,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5,
            size = 4.5 + 3 / ggplot2::.pt,
            inherit.aes = FALSE) +
  labs(x = "", y = "Number of SVs") +
  guides(fill = guide_legend(title = "Platform")) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA),
    strip.text = element_text(face = "bold", size = 17),
    axis.text.x = element_text(colour = "black", size = 17, angle = 60, hjust = 1, vjust = 1),
    axis.text.y = element_text(colour = "black", size = 17),
    axis.title = element_text(colour = "black", size = 19)
  )


ragg::agg_tiff(filename = "pics/FS2_platform_sv_number.tiff", width =3200, height =2000, res =300, compression = "lzw")
print(p_multi)
dev.off()



##############################################################################
## distribution across studies 

# Preserve a stable study order and ensure all studies are shown.
studies <- sort(unique(technique_platform$study_name))
technique_platform$study_name <- factor(technique_platform$study_name, levels = studies)

# One histogram per study, arranged as 5 rows x 7 columns.
p_dis <- ggplot(technique_platform, aes(x = sequence_platform_rename)) +
  geom_bar(fill = "#207F4C", colour = "white", linewidth = 0.1, 
           na.rm = TRUE) +
  facet_wrap(~ study_name, nrow =7, ncol = 6, scales = "free_y", drop = FALSE,
             labeller = labeller(study_name = fs2_study_labels)) +
  labs(x = "Sequencing platform", y = "Number of samples") +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "grey40"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.45, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

ragg::agg_tiff(filename = "pics/FS2_sequence_platform.tiff", width =3200, height =2800, res =300, compression = "lzw")
print(p_dis)
dev.off()




###########################################################################################
##  DNA_extraction_kit


technique_test_kit<-subset(technique_test,DNA_extraction_kit!="" & is.na(DNA_extraction_kit)==F,drop=T)

# sample size
sample_size = technique_test_kit %>% group_by(DNA_extraction_kit) %>% dplyr::summarize(num=n())

table(technique_test_kit$sequence_platform_rename)

# dsv_number
aov_dsv <- aov(dsv_number ~ DNA_extraction_kit, data = technique_test_kit)
summary(aov_dsv)

# vsv_number
aov_vsv <- aov(vsv_number ~ DNA_extraction_kit, data = technique_test_kit)
summary(aov_vsv)

# sv_number
aov_sv <- aov(sv_number ~ DNA_extraction_kit, data = technique_test_kit)
summary(aov_sv)

p_dsv <- summary(aov_dsv)[[1]]$`Pr(>F)`[1]
p_vsv <- summary(aov_vsv)[[1]]$`Pr(>F)`[1]
p_sv  <- summary(aov_sv)[[1]]$`Pr(>F)`[1]

annot_df <- data.frame(
  sv_type = factor(c("Total SV", "dSV", "vSV"),
                   levels = c("Total SV", "dSV", "vSV")),
  label = c(
    paste0("ANOVA p = ", format(p_sv, scientific = TRUE, digits = 2)),
    paste0("ANOVA p = ", format(p_dsv, scientific = TRUE, digits = 2)),
    paste0("ANOVA p = ", format(p_vsv, scientific = TRUE, digits = 2))
  )
)


table(technique_test$DNA_extraction_kit)


sv_long <- technique_test_kit %>%
  dplyr::select(match_id, DNA_extraction_kit, 
                sv_number, dsv_number, vsv_number) %>%
  pivot_longer(
    cols = c(sv_number, dsv_number, vsv_number),
    names_to = "sv_type",
    values_to = "value"
  ) %>%
  mutate(
    sv_type = factor(sv_type, 
                     levels = c("sv_number", "dsv_number", "vsv_number"),
                     labels = c("Total SV", "dSV", "vSV"))
  )


sample_size2 <- technique_test_kit %>%
  group_by(DNA_extraction_kit) %>%
  dplyr::summarise(num = n())

sv_long <- sv_long %>%
  left_join(sample_size2, by = "DNA_extraction_kit") %>%
  mutate(myaxis = paste0(DNA_extraction_kit, "\n", "n=", num))

annot_df$label <- fs2_annotation_labels(annot_df$label)
p_multi <- ggplot(sv_long, aes(x = myaxis, y = value, fill =DNA_extraction_kit)) +
  geom_violin(width = 0.8, alpha = 0.7, trim = TRUE) +
  geom_boxplot(width = 0.25, color = "black", alpha = 0.3, outlier.shape = NA) +
  scale_fill_manual(values = c("#ECD452","#BB97C5","#EA8958","#2983BB","#BC3523","#207F4C","#EE4863","#97C24E")) +
  facet_wrap(~ sv_type, scales = "free_y", nrow = 3) +
  geom_text(data = annot_df,
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5,
            size = 4.5 + 3 / ggplot2::.pt,
            inherit.aes = FALSE) +
  labs(x = "", y = "Number of SVs") +
  guides(fill = guide_legend(title = "Platform")) +
  theme(
   axis.line.x = element_line(),
   axis.line.y = element_line(),
   legend.position = "none",
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
   panel.background = element_rect(fill = NA),
   strip.text = element_text(face = "bold", size = 17),
   axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, colour = "black", size = 17),
   axis.text.y = element_text(colour = "black", size = 17),
   axis.title = element_text(colour = "black", size = 19)
  )


ragg::agg_tiff(filename = "pics/FS2_DNA_kit_sv_number.tiff", width =3200, height =2000, res =300, compression = "lzw")
print(p_multi)
dev.off()



##############################################################################
## distribution across studies 

# Preserve a stable study order and ensure all studies are shown.
studies <- sort(unique(technique_test_kit$study_name))
technique_test_kit$study_name <- factor(technique_test_kit$study_name, levels = studies)

# One histogram per study, arranged as 5 rows x 7 columns.
p_dis <- ggplot(technique_test_kit, aes(x = DNA_extraction_kit)) +
  geom_bar(fill = "#207F4C", colour = "white", linewidth = 0.1, 
           na.rm = TRUE) +
  facet_wrap(~ study_name, nrow = 5, ncol = 4, scales = "free_y", drop = FALSE,
             labeller = labeller(study_name = fs2_study_labels)) +
  labs(x = "DNA extraction kit", y = "Number of samples") +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey90", colour = "grey40"),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.45, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

ragg::agg_tiff(filename = "pics/FS2_DNA_kit.tiff", width =3200, height =2400, res =300, compression = "lzw")
print(p_dis)
dev.off()


table(technique_test_kit$study_name)



################################# the correlation bewteen the read_length and SV number

ThomasAM_2018<-subset(technique_test,study_name=="ThomasAM_2018",drop=T)

sv_long <- ThomasAM_2018 %>%
  dplyr::select(match_id, DNA_extraction_kit, sv_number, dsv_number, vsv_number) %>%
  tidyr::pivot_longer(
    cols = c(sv_number, dsv_number, vsv_number),
    names_to = "sv_type",
    values_to = "value"
  ) %>%
  mutate(
    sv_type = factor(sv_type,
                     levels = c("sv_number", "dsv_number", "vsv_number"),
                     labels = c("Total SV", "dSV", "vSV"))
  )


aov_sv  <- aov(sv_number  ~ DNA_extraction_kit, data = ThomasAM_2018)
aov_dsv <- aov(dsv_number ~ DNA_extraction_kit, data = ThomasAM_2018)
aov_vsv <- aov(vsv_number ~ DNA_extraction_kit, data = ThomasAM_2018)

p_sv  <- summary(aov_sv)[[1]]$`Pr(>F)`[1]
p_dsv <- summary(aov_dsv)[[1]]$`Pr(>F)`[1]
p_vsv <- summary(aov_vsv)[[1]]$`Pr(>F)`[1]

f_sv  <- summary(aov_sv)[[1]]$`F value`[1]

annot_df <- data.frame(
  sv_type = factor(c("Total SV", "dSV", "vSV"),
                   levels = c("Total SV", "dSV", "vSV")),
  label = c(
    paste0("ANOVA F = ", round(f_sv, 2), 
           ", p = ", format(p_sv, scientific = TRUE, digits = 2)),
    paste0("ANOVA p = ", format(p_dsv, scientific = TRUE, digits = 2)),
    paste0("ANOVA p = ", format(p_vsv, scientific = TRUE, digits = 2))
  )
)

annot_df$label <- fs2_annotation_labels(annot_df$label)
p_facet <- ggplot(sv_long, aes(x = DNA_extraction_kit, y = value)) +
  geom_boxplot(aes(fill = DNA_extraction_kit), 
               width = 0.6, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.2, size = 1, color = "#2983BB") +
  geom_text(data = annot_df, aes(x = 1, y = Inf, label = label),
            hjust = 0, vjust = 1.5, size = 4 + 3 / ggplot2::.pt, fontface = "bold",
            inherit.aes = FALSE) +
  facet_wrap(~ sv_type, ncol = 3, scales = "free_y") +
  labs(x = "DNA extraction kit", y = "Number of SVs") +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid.major = element_line(colour = "gray90", linetype = "dashed", linewidth = 0.1),
    panel.grid.minor = element_line(colour = NA),
    panel.background = element_rect(fill = NA),
    strip.text = element_text(face = "bold", size = 17),
    axis.text.x = element_text(angle = 60, hjust = 1, colour = "black", size = 13, vjust = 1),
    axis.text.y = element_text(colour = "black", size = 15),
    axis.title = element_text(colour = "black", size = 17),
    legend.position = "none"
  )

ragg::agg_tiff(filename = "pics/FS2_ThomasAM_2018_DNA_extraction_kit_sv_number.tiff", width =3200, height =2000, res =300, compression = "lzw")
print(p_facet)
dev.off()



####################################################################################
### Median_read_length

technique_read_length<-technique_test

table(technique_read_length$study_name,technique_read_length$Median_read_length)

technique_read_length$study_name <- as.character(technique_read_length$study_name)
technique_read_length$Median_read_length <- suppressWarnings(as.numeric(technique_read_length$Median_read_length))
technique_read_length <- technique_read_length[!is.na(technique_read_length$study_name) & is.finite(technique_read_length$Median_read_length), ]
if (!nrow(technique_read_length)) stop("No valid observations remain after filtering.")

################################# the correlation bewteen the read_length and SV number
MetaCardis_2020_a<-subset(technique_read_length,study_name=="MetaCardis_2020_a",drop=T)
MetaCardis_2020_a$Median_read_length <- suppressWarnings(as.numeric(MetaCardis_2020_a$Median_read_length))

sv_long <- MetaCardis_2020_a %>%
  dplyr::select(match_id, Median_read_length, sv_number, dsv_number, vsv_number) %>%
  tidyr::pivot_longer(
    cols = c(sv_number, dsv_number, vsv_number),
    names_to = "sv_type",
    values_to = "value"
  ) %>%
  mutate(
    sv_type = factor(sv_type,
                     levels = c("sv_number", "dsv_number", "vsv_number"),
                     labels = c("Total SV", "dSV", "vSV"))
  )

cor_sv <- cor.test(MetaCardis_2020_a$sv_number, MetaCardis_2020_a$Median_read_length)
cor_dsv <- cor.test(MetaCardis_2020_a$dsv_number, MetaCardis_2020_a$Median_read_length)
cor_vsv <- cor.test(MetaCardis_2020_a$vsv_number, MetaCardis_2020_a$Median_read_length)

cor_p_labels <- base::format.pval(c(cor_sv$p.value, cor_dsv$p.value, cor_vsv$p.value),
                           digits = 2, eps = .Machine$double.eps)
cor_p_labels <- ifelse(grepl("^<", cor_p_labels), cor_p_labels, paste0("= ", cor_p_labels))

annot_df <- data.frame(
  sv_type = factor(c("Total SV", "dSV", "vSV"),
                   levels = c("Total SV", "dSV", "vSV")),
  label = c(
    paste0("r = ", round(cor_sv$estimate, 3)),
    paste0("r = ", round(cor_dsv$estimate, 3)),
    paste0("r = ", round(cor_vsv$estimate, 3))
  )
)

annot_df$label <- fs2_annotation_labels(annot_df$label)
p_scatter_facet <- ggplot(sv_long, aes(x = Median_read_length, y = value)) +
  geom_point(col = "#2983BB", alpha = 0.4, size = 1) +
  #geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, color = "#BC3523") +
  geom_text(data = annot_df, aes(x = 82, y = Inf, label = label),
            hjust = 0, vjust = 1.5, size = 4 + 3 / ggplot2::.pt, fontface = "bold",
            inherit.aes = FALSE) +
  facet_wrap(~ sv_type, ncol = 3, scales = "free_y") +
  labs(x = "Median read length", y = "Number of SVs",title="MetaCardis_2020_a") +
  coord_cartesian(xlim = c(80, 110)) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid.major = element_line(colour = "gray90", linetype = "dashed", linewidth = 0.1),
    panel.grid.minor = element_line(colour = NA),
    panel.background = element_rect(fill = NA),
    strip.text = element_text(face = "bold", size = 17),
    axis.text = element_text(colour = "black", size = 15),
    axis.title = element_text(colour = "black", size = 17)
  )

ragg::agg_tiff(filename = "pics/FS2_MetaCardis_2020_a_Median_read_length_sv_number.tiff", width =3200, height =2000, res =300, compression = "lzw")
print(p_scatter_facet)
dev.off()


################################# the correlation bewteen the read_length and SV number

table(technique_test$study_name)
CosteaPI_2017<-subset(technique_read_length,study_name=="CosteaPI_2017",drop=T)

sv_long <- CosteaPI_2017 %>%
  dplyr::select(match_id, Median_read_length, sv_number, dsv_number, vsv_number) %>%
  tidyr::pivot_longer(
    cols = c(sv_number, dsv_number, vsv_number),
    names_to = "sv_type",
    values_to = "value"
  ) %>%
  mutate(
    sv_type = factor(sv_type,
                     levels = c("sv_number", "dsv_number", "vsv_number"),
                     labels = c("Total SV", "dSV", "vSV"))
  )

cor_sv <- cor.test(CosteaPI_2017$sv_number, CosteaPI_2017$Median_read_length)
cor_dsv <- cor.test(CosteaPI_2017$dsv_number, CosteaPI_2017$Median_read_length)
cor_vsv <- cor.test(CosteaPI_2017$vsv_number, CosteaPI_2017$Median_read_length)

cor_p_labels <- base::format.pval(c(cor_sv$p.value, cor_dsv$p.value, cor_vsv$p.value),
                           digits = 2, eps = .Machine$double.eps)
cor_p_labels <- ifelse(grepl("^<", cor_p_labels), cor_p_labels, paste0("= ", cor_p_labels))

annot_df <- data.frame(
  sv_type = factor(c("Total SV", "dSV", "vSV"),
                   levels = c("Total SV", "dSV", "vSV")),
  label = c(
    paste0("r = ", round(cor_sv$estimate, 3), ", p ", cor_p_labels[1]),
    paste0("r = ", round(cor_dsv$estimate, 3), ", p ", cor_p_labels[2]),
    paste0("r = ", round(cor_vsv$estimate, 3), ", p ", cor_p_labels[3])
  )
)

annot_df$label <- fs2_annotation_labels(annot_df$label)
p_scatter_facet <- ggplot(sv_long, aes(x = Median_read_length, y = value)) +
  geom_point(col = "#2983BB", alpha = 0.4, size = 1) +
  #geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, color = "#BC3523") +
  geom_text(data = annot_df, aes(x = 80, y = Inf, label = label),
            hjust = 0, vjust = 1.5, size = 4 + 3 / ggplot2::.pt, fontface = "bold",
            inherit.aes = FALSE) +
  facet_wrap(~ sv_type, ncol = 3, scales = "free_y") +
  labs(x = "Median read length", y = "Number of SVs") +
  coord_cartesian(xlim = c(80, 110)) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid.major = element_line(colour = "gray90", linetype = "dashed", linewidth = 0.1),
    panel.grid.minor = element_line(colour = NA),
    panel.background = element_rect(fill = NA),
    strip.text = element_text(face = "bold", size = 17),
    axis.text = element_text(colour = "black", size = 15),
    axis.title = element_text(colour = "black", size = 17)
  )

ragg::agg_tiff(filename = "pics/FS2_CosteaPI_2017_Median_read_length_sv_number.tiff", width =3200, height =2000, res =300, compression = "lzw")
print(p_scatter_facet)
dev.off()



################################# the correlation bewteen the read_length and SV number

table(technique_test$study_name)

ZellerG_2014<-subset(technique_read_length,study_name=="ZellerG_2014",drop=T)

sv_long <- ZellerG_2014 %>%
  dplyr::select(match_id, Median_read_length, sv_number, dsv_number, vsv_number) %>%
  tidyr::pivot_longer(
    cols = c(sv_number, dsv_number, vsv_number),
    names_to = "sv_type",
    values_to = "value"
  ) %>%
  mutate(
    sv_type = factor(sv_type,
                     levels = c("sv_number", "dsv_number", "vsv_number"),
                     labels = c("Total SV", "dSV", "vSV"))
  )

cor_sv <- cor.test(ZellerG_2014$sv_number, ZellerG_2014$Median_read_length)
cor_dsv <- cor.test(ZellerG_2014$dsv_number, ZellerG_2014$Median_read_length)
cor_vsv <- cor.test(ZellerG_2014$vsv_number, ZellerG_2014$Median_read_length)

cor_p_labels <- base::format.pval(c(cor_sv$p.value, cor_dsv$p.value, cor_vsv$p.value),
                           digits = 2, eps = .Machine$double.eps)
cor_p_labels <- ifelse(grepl("^<", cor_p_labels), cor_p_labels, paste0("= ", cor_p_labels))

annot_df <- data.frame(
  sv_type = factor(c("Total SV", "dSV", "vSV"),
                   levels = c("Total SV", "dSV", "vSV")),
  label = c(
    paste0("r = ", round(cor_sv$estimate, 3), ", p ", cor_p_labels[1]),
    paste0("r = ", round(cor_dsv$estimate, 3), ", p ", cor_p_labels[2]),
    paste0("r = ", round(cor_vsv$estimate, 3), ", p ", cor_p_labels[3])
  )
)


annot_df$label <- fs2_annotation_labels(annot_df$label)
p_scatter_facet <- ggplot(sv_long, aes(x = Median_read_length, y = value)) +
  geom_point(col = "#2983BB", alpha = 0.4, size = 1) +
  #geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, color = "#BC3523") +
  geom_text(data = annot_df, aes(x = 72, y = Inf, label = label),
            hjust = 0, vjust = 1.5, size = 4 + 3 / ggplot2::.pt, fontface = "bold",
            inherit.aes = FALSE) +
  facet_wrap(~ sv_type, ncol = 3, scales = "free_y") +
  labs(x = "Median read length", y = "Number of SVs",title = "Study name: ZellerG_2014") +
  coord_cartesian(xlim = c(60, 110)) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid.major = element_line(colour = "gray90", linetype = "dashed", linewidth = 0.1),
    panel.grid.minor = element_line(colour = NA),
    panel.background = element_rect(fill = NA),
    plot.title = element_text(size = 16.2),
    strip.text = element_text(face = "bold", size = 17),
    axis.text = element_text(colour = "black", size = 15),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    axis.title = element_text(colour = "black", size = 17)
  )

ragg::agg_tiff(filename = "pics/FS2_ZellerG_2014_Median_read_length_sv_number.tiff", width =2200, height =1200, res =300, compression = "lzw")
print(p_scatter_facet)
dev.off()



############################################################################
### Median/average read length and SV counts: matching Pearson summaries

read_length_correlations <- data.frame()
for (measure in c("Median_read_length", "Average_read_length", "read_count")) {
  read_length_studies <- if (measure == "read_count") "All samples" else
    c("All samples", "MetaCardis_2020_a", "CosteaPI_2017", "WirbelJ_2018", "ZellerG_2014")
  for (read_length_study in read_length_studies) {
    read_length_data <- if (read_length_study == "All samples") technique_test else
      technique_test[technique_test$study_name == read_length_study, , drop = FALSE]
    for (sv_type in c("sv_number", "dsv_number", "vsv_number")) {
      valid_pair <- is.finite(read_length_data[[measure]]) & read_length_data[[measure]] > 0 &
        is.finite(read_length_data[[sv_type]])
      pair_data <- read_length_data[valid_pair, , drop = FALSE]
      correlation_r <- NA_real_
      correlation_p <- NA_real_
      correlation_status <- "ok"
      if (nrow(pair_data) < 3L) {
        correlation_status <- "insufficient_pairs"
      } else if (length(unique(pair_data[[measure]])) < 2L) {
        correlation_status <- "constant_measure"
      } else if (length(unique(pair_data[[sv_type]])) < 2L) {
        correlation_status <- "constant_sv_count"
      } else {
        correlation_fit <- cor.test(pair_data[[sv_type]], pair_data[[measure]], method = "pearson")
        correlation_r <- unname(correlation_fit$estimate)
        correlation_p <- correlation_fit$p.value
      }
      read_length_correlations <- rbind(read_length_correlations, data.frame(
        study_name = read_length_study, measure = measure, sv_type = sv_type,
        N = nrow(pair_data), r = correlation_r, pvalue = correlation_p,
        status = correlation_status, stringsAsFactors = FALSE))
    }
  }
}
write.csv(read_length_correlations,
          file.path(technical_dir, "read_length_correlations.csv"), row.names = FALSE)


############################################################################
### Median_read_length: distribution across cohorts

technique_median_read_length <- technique_test[
  !is.na(technique_test$study_name) & is.finite(technique_test$Median_read_length) &
    technique_test$Median_read_length > 0, , drop = FALSE]
studies <- sort(unique(technique_median_read_length$study_name))
technique_median_read_length$study_name <- factor(technique_median_read_length$study_name, levels = studies)

p_median_distribution <- ggplot(technique_median_read_length, aes(x = Median_read_length)) +
  geom_histogram(binwidth = 1, boundary = 0, closed = "left",
                 fill = "#207F4C", colour = "white", linewidth = 0.01, na.rm = TRUE) +
  facet_wrap(~ study_name, nrow =5, ncol = 9, scales = "free_y", drop = FALSE,
             labeller = labeller(study_name = fs2_study_labels)) +
  labs(x = "Median read length", y = "Number of samples") +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey90", colour = "grey40"),
        strip.text = element_text(size = 8, face = "bold"),
        axis.title = element_text(size = 9), axis.text = element_text(size = 9),
        panel.grid.minor = element_blank(), panel.spacing = unit(0.45, "lines"),
        plot.margin = margin(5, 5, 5, 5)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

ragg::agg_tiff(filename = "pics/FS2_Median_read_length_sv_number.tiff", width = 3500, height =2000, res = 300, compression = "lzw")
print(p_median_distribution)
dev.off()




################################# the correlation bewteen the sv number and read count

sv_long <- technique_test %>%
  dplyr::select(match_id, read_count, sv_number, dsv_number, vsv_number) %>%
  tidyr::pivot_longer(
    cols = c(sv_number, dsv_number, vsv_number),
    names_to = "sv_type",
    values_to = "value"
  ) %>%
  mutate(
    sv_type = factor(sv_type,
                     levels = c("sv_number", "dsv_number", "vsv_number"),
                     labels = c("Total SV", "dSV", "vSV"))
  )

cor_sv <- cor.test(technique_test$sv_number, technique_test$read_count)
cor_dsv <- cor.test(technique_test$dsv_number, technique_test$read_count)
cor_vsv <- cor.test(technique_test$vsv_number, technique_test$read_count)

cor_p_labels <- base::format.pval(c(cor_sv$p.value, cor_dsv$p.value, cor_vsv$p.value),
                           digits = 2, eps = .Machine$double.eps)
cor_p_labels <- ifelse(grepl("^<", cor_p_labels), cor_p_labels, paste0("= ", cor_p_labels))

annot_df <- data.frame(
  sv_type = factor(c("Total SV", "dSV", "vSV"),
                   levels = c("Total SV", "dSV", "vSV")),
  label = c(
    paste0("r = ", round(cor_sv$estimate, 3), ", p ", cor_p_labels[1]),
    paste0("r = ", round(cor_dsv$estimate, 3), ", p ", cor_p_labels[2]),
    paste0("r = ", round(cor_vsv$estimate, 3), ", p ", cor_p_labels[3])
  )
)

annot_df$label <- fs2_annotation_labels(annot_df$label)
p_scatter_facet <- ggplot(sv_long, aes(x = read_count / 1000000, y = value)) +
  geom_point(col = "#2983BB", alpha = 0.4, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.2, color = "#BC3523") +
  geom_text(data = annot_df, aes(x = 2, y = Inf, label = label),
            hjust = 0, vjust = 1.5, size = 4 + 3 / ggplot2::.pt, fontface = "bold",
            inherit.aes = FALSE) +
  facet_wrap(~ sv_type, ncol = 3, scales = "free_y") +
  labs(x = "Read count (million)", y = "Number of SVs") +
  coord_cartesian(xlim = c(0, 110)) +
  theme(
    axis.line.x = element_line(),
    axis.line.y = element_line(),
    panel.grid.major = element_line(colour = "gray90", linetype = "dashed", linewidth = 0.1),
    panel.grid.minor = element_line(colour = NA),
    panel.background = element_rect(fill = NA),
    strip.text = element_text(face = "bold", size = 17),
    axis.text = element_text(colour = "black", size = 15),
    axis.title = element_text(colour = "black", size = 17)
  )


ragg::agg_tiff(filename = "pics/FS2_sequence_depth_sv_number.tiff", width =3200, height =2000, res =300, compression = "lzw")
print(p_scatter_facet)
dev.off()


###############################################################################################
############# 2.3 SV number factors
##Show relationships between SV numbers and genome size and sample size.

species_sample_n<-info[,c("Length","SVs_number","Total_samples_number","Deletion_SVs_number","Variable_SVs_number")]

length_sv_n_cor<-cor.test(species_sample_n$Length/1000, species_sample_n$SVs_number)
text_r<-paste("r=",round(length_sv_n_cor$estimate,digits = 3),"\np=",format(length_sv_n_cor$p.value,digits = 3),sep = "")

p_scatter_pie<-ggplot()+
  geom_smooth(data=species_sample_n ,aes(x=Length/10000, y=SVs_number),method = "lm", color = "white", alpha = 0.2,linewidth = 0.5)+
  scatterpie::geom_scatterpie(data=species_sample_n ,aes(x=Length/10000, y=SVs_number,r=Total_samples_number/500),
                  cols=c("Deletion_SVs_number","Variable_SVs_number"),color=NA, alpha = 0.75)+
  coord_equal()+
  annotate("text", -Inf, Inf, label = c(text_r),hjust = -0.1, vjust = 1)+
  xlab("Genome size (10 kbp)")+
  ylab("Number of SVs")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name   = NULL,
                      breaks = c("Deletion_SVs_number", "Variable_SVs_number"),
                      labels = c("Deletion SVs              ", "Variable SVs"),
                      values = mycolor2_green_blue)+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position="top",
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = "gray90",linetype = "dashed",linewidth = 0.25),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))


#pdf("pics/Sample_n_SVs_n_scatter_pie.pdf", height = 5, width = 4)
#print(p_scatter_pie)
#dev.off()

tiff(file = "pics/FS3_Sample_n_SVs_n_scatter_pie.tiff", width =1200, height =800, res =300) 
print(p_scatter_pie)
dev.off()



## 2 SV summary
## 2.1 Species number
## Number of species with SVs

if(!dir.exists("02.SV_summary_statistics")){dir.create("02.SV_summary_statistics")}

dsgv_species_n<-sum(info$Deletion_SVs_number>0)
vsgv_species_n<-sum(info$Variable_SVs_number>0)
total_species_n<-sum(info$Variable_SVs_number>0 | info$Deletion_SVs_number>0)

species_n<-data.frame(item = rep("Informative species number", 3),
                      categories = c("SVs","Deletion SVs", "Variable SVs"),
                      value = c(total_species_n,dsgv_species_n, vsgv_species_n))

species_n$categories <- factor(species_n$categories, levels = species_n$categories)
species_n$categories <- factor(species_n$categories,levels(species_n$categories)[c(2,3,1)])

p_species_n<-ggplot(species_n, aes(x=categories, y=value,label = value))+
  geom_bar(aes(fill = categories),stat = 'identity')+
  geom_text(size = ggplot2::GeomText$default_aes$size + 3 / ggplot2::.pt, position = position_stack(vjust = 0.5), color = "white")+
  xlab(NULL)+
  ylab("Number of Informative species")+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name=NULL,
                    breaks=c("Deletion SVs", "Variable SVs", "SVs"),
                    labels=c("Deletion SVs", "Variable SVs", "SVs"),
                    values = mycolor3_red_blue_yellow )+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.key = element_rect(fill = NA), 
        legend.position = "none",
        panel.grid.major = element_line(colour = "gray90",linetype = "dashed",linewidth = 0.25),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 11.8),
        axis.title = element_text(size = 14))

ragg::agg_tiff(filename = "pics/FS2_Informative_species_number.tiff", width =4000, height =3200, res =300, compression = "lzw")
print(p_species_n)
dev.off()


### 2.1 SV number
#### 2.1.1 SV total number

## height<-c(5,10)
## barplot(height,col="#EA517F")
## barplot(height,col="#696969")

dsgv_n<-sum(info$Deletion_SVs_number)
vsgv_n<-sum(info$Variable_SVs_number)
sgv_n<-dsgv_n+vsgv_n

sv_n<-data.frame(items = rep("SVs number", 2),
                      categories = c("Deletion SV", "Variable SV"),
                      value = c(dsgv_n, vsgv_n))

p_pie<-my_pie(sv_n, 'SVs number',mycol =c("#EA517F","#2376B7") )
p_pie_nolabel<-my_pie_nolabel(sv_n, 'SVs number',mycol =c("#EA517F","#2376B7") )

#pdf("pics/Total_SVs_number.pdf",height =5, width = 5)
#print(p_pie)
#dev.off()

tiff(file = "pics/FS3_Total_SVs_number.tiff", width =1000, height =1000, res =300) 
print(p_pie)
dev.off()

tiff(file = "pics/FS3_Total_SVs_number_nolabel.tiff", width =1000, height =1000, res =300) 
print(p_pie_nolabel)
dev.off()



#### 2.1.2 SV number per species

info_svs<-aggregate(info$SVs_number, by=list(Category= info$Short_name), FUN=sum)
info_svs<-data.frame(info_svs)
colnames(info_svs)<-c("Short_name","SVs_number")

species_sgv_n_order<- info_svs$Short_name[order(info_svs$SVs_number, decreasing = T)]

#species_sgv_n_order<- info$Short_name[order(info$SVs_number, decreasing = T)]
species_sgv_n_long<-gather(info[,c(10,12,13)], "type", "number", c(2:3))

p_sv_species<-ggplot(species_sgv_n_long, aes(x=Short_name, y=number, group = type, fill = type))+
  geom_bar(position = "stack",stat="identity")+
  xlab(NULL)+
  ylab("Number of SVs")+
  scale_x_discrete(limits = species_sgv_n_order)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name=NULL,
                      breaks = c("Deletion_SVs_number", "Variable_SVs_number"),
                      labels = c("Deletion SVs              ", "Variable SVs"),
                      values = c("#2376B7","#EA517F"))+ggtitle("")+
  theme(plot.subtitle = element_text(vjust =1), 
        plot.caption = element_text(vjust = 1), 
        plot.margin=unit(c(1,1,1,1),"cm"),
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 65, hjust = 1,size = 12,family ="Arial",color="black",face = "italic"),
        axis.text.y = element_text(size = 15,family ="Arial",color="black"),
        axis.title.y = element_text(size = 15,family ="Arial",color="black"),
        plot.title = element_text(size = 15,family ="Arial",color="black"),
        legend.position="inside",
        legend.position.inside=c(0.2,0.9),
        legend.key = element_rect(fill = NA), 
        legend.text = element_text(size =15,family ="Arial",color="black"),
        panel.grid.major = element_line(colour = NA,linetype = "dashed",linewidth = 0.2),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))


tiff(file = "pics/FS3_SVs_number_species.tiff", width =6500, height =1900, res =300) 
print(p_sv_species)
dev.off()


### 2.2 Sample size per species
### different cohorts (differernt country or wester or non-western)
infor_sample_n_order<- info$Short_name[order(info$Total_samples_number, decreasing = T)]
infor_sample_n <- info[,c("Short_name","Africa_sample_number","America_sample_number","Asia_sample_number","Europe_sample_number")]

infor_sample_n.long <- gather(infor_sample_n,'Cohort', 'Sample_size', c(2:5))

p_sample_n<-ggplot(infor_sample_n.long, aes(x=Short_name, y=Sample_size,group = Cohort))+
  geom_bar(aes(fill = Cohort),position = "stack",stat="identity")+
  xlab(NULL)+
  ylab("Number of samples")+
  scale_x_discrete(limits = infor_sample_n_order)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name=NULL,
                      breaks = c("Africa_sample_number", "America_sample_number","Asia_sample_number","Europe_sample_number"),
                      labels = c("Africa", "Americas","Asia","Europe"),
                      values =c("#2983BB","#EA517F","#ECD452","#207F4C"))+
 theme(plot.subtitle = element_text(vjust =1), 
        plot.caption = element_text(vjust = 1), 
        plot.margin=unit(c(0,1,1,1),"cm"),
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 65, hjust = 1,size = 11,family ="sans",color="black",face = "italic"),
        axis.text.y = element_text(size = 11,family ="sans",color="black"),
        axis.title.y = element_text(size = 15,family ="sans",color="black"),
        legend.position="top",
        legend.key = element_rect(fill = NA), 
        legend.text = element_text(size =15,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA,linetype = "dashed",linewidth = 0.2),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))


#pdf("03.SV_summary_statistics/FS2_Samples_number_species.pdf", height = 6, width = 11) 
#print(p_sample_n)
#dev.off()

tiff(file = "pics/FS3_Samples_number_species.tiff", width =6500, height =1900, res =300) 
print(p_sample_n)
dev.off()



###############################################################################################
############# 2.3 SV number factors
##Show relationships between SV numbers and genome size and sample size.

species_sample_n<-info[,c("Length","SVs_number","Total_samples_number","Deletion_SVs_number","Variable_SVs_number")]

length_sv_n_cor<-cor.test(species_sample_n$Length/1000, species_sample_n$SVs_number)
text_r<-paste("r=",round(length_sv_n_cor$estimate,digits = 3),"\np=",format(length_sv_n_cor$p.value,digits = 3),sep = "")

p_scatter_pie<-ggplot()+
  geom_smooth(data=species_sample_n ,aes(x=Length/10000, y=SVs_number),method = "lm", color = "white", alpha = 0.2,linewidth = 0.5)+
  scatterpie::geom_scatterpie(data=species_sample_n ,aes(x=Length/10000, y=SVs_number,r=Total_samples_number/500),
                  cols=c("Deletion_SVs_number","Variable_SVs_number"),color=NA, alpha = 0.75)+
  coord_equal()+
  annotate("text", -Inf, Inf, label = c(text_r),hjust = -0.1, vjust = 1)+
  xlab("Genome size (10 kbp)")+
  ylab("Number of SVs")+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name   = NULL,
                      breaks = c("Deletion_SVs_number", "Variable_SVs_number"),
                      labels = c("Deletion SVs              ", "Variable SVs"),
                      values = mycolor2_green_blue)+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position="top",
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = "gray90",linetype = "dashed",linewidth = 0.25),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))


#pdf("pics/Sample_n_SVs_n_scatter_pie.pdf", height = 5, width = 4)
#print(p_scatter_pie)
#dev.off()

tiff(file = "pics/FS3_Sample_n_SVs_n_scatter_pie_2.tiff", width =1200, height =800, res =300) 
print(p_scatter_pie)
dev.off()




### 3 Species abundance
sv_species_total_abun<-data.frame(id = rownames(all_abun_sv),
                                  Total_abundance=rowSums(all_abun_sv))

mean(sv_species_total_abun$Total_abundance,na.rm = T)
se(sv_species_total_abun$Total_abundance)
min(sv_species_total_abun$Total_abundance,na.rm = T)
max(sv_species_total_abun$Total_abundance,na.rm = T)

inter<-all_basic[,c("match_id","inter")]
colnames(sv_species_total_abun)[1]<-"match_id"
sv_species_total_abun<-merge(sv_species_total_abun,inter,by.x="match_id",by.y="match_id")

mean_abundance<-c(
mean(subset(sv_species_total_abun,inter=="Europe",drop=T)$Total_abundance,na.rm = T),
mean(subset(sv_species_total_abun,inter=="Americas",drop=T)$Total_abundance,na.rm = T),
mean(subset(sv_species_total_abun,inter=="Asia",drop=T)$Total_abundance,na.rm = T),
mean(subset(sv_species_total_abun,inter=="Africa",drop=T)$Total_abundance,na.rm = T))

p_sv_abun_density<-ggplot(sv_species_total_abun,aes(x=Total_abundance,,color = inter,fill =inter))+
  geom_density(alpha = 0.2)+
  geom_rug(length = unit(0.05, "npc"))+
  geom_vline(xintercept  = mean_abundance, linetype = "dashed",color =c("#207F4C","#EA517F","#ECD452","#2983BB"))+
  ylab('Density')+
  xlab('Abundance')+ggtitle("All")+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(name=NULL,
                    breaks=c("Europe","Americas","Asia","Africa"),
                    labels=c("Europe","Americas","Asia","Africa"),
                    values = c("#207F4C","#EA517F","#ECD452","#2983BB")  )+
  scale_fill_manual(name=NULL,
                    breaks=c("Europe","Americas","Asia","Africa"),
                    labels=c("Europe","Americas","Asia","Africa"),
                    values = c("#207F4C","#EA517F","#ECD452","#2983BB")  )+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(size=10,family ="sans"),
        axis.text.y = element_text(size=10,family ="sans"),
        axis.title.x = element_text(size=10,family ="sans"),
        axis.title.y = element_text(size=10,family ="sans"),
        plot.title= element_text(size=10,family ="sans"),
        legend.position = "right",
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))

#pdf("pics/sv_abun_density.pdf", width = 3, height = 2)
#print(p_sv_abun_density)
#dev.off()

tiff(file = "pics/FS5_sv_abun_density.tiff", width =1500, height =1000, res =300) 
print(p_sv_abun_density)
dev.off()





