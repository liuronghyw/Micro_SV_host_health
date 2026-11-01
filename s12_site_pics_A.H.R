### barplot for sites (BMI and age)
### 2026-10-23
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

## 1 Preparation

library("survival")
library("survminer")
library(gridExtra)
library("grid")

library("patchwork")
source("functions.R")

###################### 1.2 Inputs

all_basic<-read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv",header=T)
all_basic$study_name <- normalize_study_name(all_basic$study_name)

### Include all health statuses 
# Include all eligible participants, consistent with the primary BMI analysis.

load("01.cleanData/SV_all/dsgv.RData")

dsgv<-dsgv_sub

#######################################################################################
### clinical data and dsgv/vsgv
all_basic$match_id<-row.names(all_basic)

dsgv_pic<-dsgv[,c(
"Anaerostipes hadrus SSC/2:1577_1579 and 9 segments"
), drop = FALSE]

### draw pictures
library("ggbeeswarm")

dsgv_pic$match_id<-row.names(dsgv_pic)
pic_dsv<-merge(all_basic,dsgv_pic,by.x="match_id",by.y="match_id")

#########################################################################
##### Anaerostipes hadrus SSC/2:1577_1579 and 9 segments
library("ggplot2")

dsv_target <- "Anaerostipes hadrus SSC/2:1577_1579 and 9 segments"
pic_dsv$temp <- pic_dsv[[dsv_target]]
pic_sub <- subset(pic_dsv, !is.na(temp) & !is.na(BMI), drop = TRUE)
pic_sub$temp <- as.numeric(as.character(pic_sub$temp))
if (!all(pic_sub$temp %in% c(0, 1))) {
  stop("Unexpected dSV coding. Expected 0 = deletion and 1 = no deletion.")
}

f6d_data <- pic_sub
f6d_data$study_name <- sub(
  "Fern[^_]*ndez_Calvet_2024", "Fern\u00e1ndez_Calvet_2024",
  as.character(f6d_data$study_name)
)
# Select every cohort with more than 100 samples in this dSV association.
bmi_dsv_file <- "07.Microbial_GWAS/studies_adjAbun/BMI_dsv.csv"
bmi_dsv_results <- read.csv(bmi_dsv_file, check.names = FALSE)
names(bmi_dsv_results) <- sub("Fern[^_]*ndez_Calvet_2024", "Fern\u00e1ndez_Calvet_2024", names(bmi_dsv_results))
result_row <- which(bmi_dsv_results$X == dsv_target)
if (length(result_row) != 1L) {
  stop("Expected exactly one matching dSV row in ", bmi_dsv_file)
}
sample_columns <- grep("[.]N$", names(bmi_dsv_results), value = TRUE)
study_n <- setNames(as.numeric(bmi_dsv_results[result_row, sample_columns]),
                    sub("[.]N$", "", sample_columns))
study_n <- sort(study_n[is.finite(study_n) & study_n > 100L], decreasing = TRUE)
study_order <- names(study_n)
if (!length(study_order)) stop("No cohorts have more than 100 association samples.")
f6d_data <- subset(f6d_data, study_name %in% study_order, drop = TRUE)
stopifnot(all(as.numeric(table(f6d_data$study_name)[study_order]) == study_n))
study_labels <- setNames(gsub("_", " ", study_order), study_order)
study_labels["PiccinnoG_2025_cohort5"] <- "PiccinnoG 2025\ncohort 5"
study_labels["Fern\u00e1ndez_Calvet_2024"] <- "Fern\u00e1ndez Calvet\n2024"
study_palette <- c("0" = "#7876B1FF", "1" = "#207F4C")

plot_list <- vector("list", length(study_order))
names(plot_list) <- study_order
delta_rows <- vector("list", length(study_order))
names(delta_rows) <- study_order
bmi_limits <- c(min(15, min(f6d_data$BMI, na.rm = TRUE)), max(f6d_data$BMI, na.rm = TRUE))
bmi_breaks <- seq(15, max(45, floor(max(bmi_limits) / 15) * 15), by = 15)

for (study in study_order) {
  plot_data <- subset(f6d_data, study_name == study, drop = TRUE)
  n_del <- sum(plot_data$temp == 0)
  n_no_del <- sum(plot_data$temp == 1)
  if (n_del == 0L || n_no_del == 0L) {
    stop("Both dSV groups are required for study: ", study)
  }

  median_del <- median(plot_data$BMI[plot_data$temp == 0], na.rm = TRUE)
  median_no_del <- median(plot_data$BMI[plot_data$temp == 1], na.rm = TRUE)
  # Delta BMI = median(no deletion, dSV = 1) - median(deletion, dSV = 0).
  delta_median <- median_no_del - median_del
  delta_rows[[study]] <- data.frame(
    study_name = study,
    n_deletion = n_del,
    n_no_deletion = n_no_del,
    median_BMI_deletion = median_del,
    median_BMI_no_deletion = median_no_del,
    delta_median_BMI_no_deletion_minus_deletion = delta_median,
    stringsAsFactors = FALSE
  )

  plot_data$SV_group <- factor(
    as.character(plot_data$temp),
    levels = c("0", "1"),
    labels = c(
      paste0("Deletion\n(n=", n_del, ")"),
      paste0("No deletion\n(n=", n_no_del, ")")
    )
  )
  plot_list[[study]] <- ggplot(
    plot_data,
    aes(x = SV_group, y = BMI, color = as.character(temp))
  ) +
    geom_hline(yintercept = bmi_breaks, linetype = "dashed", color = "grey65", linewidth = 0.8) +
    geom_violin(alpha = 0.05, width = 0.9, linewidth = 1.1, trim = TRUE) +
    ggbeeswarm::geom_quasirandom(width = 0.18, size = 0.7, alpha = 0.55) +
    geom_boxplot(
      alpha = 0.18, width = 0.24, color = "black", linewidth = 0.9,
      outlier.shape = NA
    ) +
    scale_color_manual(values = study_palette, guide = "none") +
    scale_y_continuous(limits = bmi_limits, breaks = bmi_breaks) +
    labs(
      title = paste0(paste(strwrap(gsub("\n", " ", unname(study_labels[study])), width = 14), collapse = "\n"), "\nDelta = ", sprintf("%+.2f", delta_median)),
      x = NULL, y = if (study == study_order[1]) "BMI" else NULL
    ) +
    theme_cowplot() +
    theme(
      panel.grid = element_blank(), legend.position = "none",
      axis.line = element_line(linewidth = 0.9), axis.ticks = element_line(linewidth = 0.8),
      text = element_text(size = 16), axis.text = element_text(size = 16),
      axis.title = element_text(size = 18),
      axis.title.y = if (study == study_order[1]) element_text(size = 18) else element_blank(),
      plot.title = element_text(size = 16, face = "plain", hjust = 0.5, vjust = 0,
                                lineheight = 1.05, margin = margin(t = 12, b = 10)),
      axis.text.y = if (study == study_order[1]) element_text(size = 16) else element_blank(),
      axis.ticks.y = if (study == study_order[1]) element_line() else element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      plot.margin = unit(c(8, 8, 8, 8), "pt")
    )
}

delta_summary <- do.call(rbind, delta_rows)
delta_summary$study_label <- unname(study_labels[delta_summary$study_name])

beta_rows <- lapply(study_order, function(study) {
  beta_col <- paste0(study, ".beta")
  se_col <- paste0(study, ".SE")
  if (!all(c(beta_col, se_col) %in% names(bmi_dsv_results))) {
    stop("Missing beta or SE columns for study: ", study)
  }
  beta <- as.numeric(bmi_dsv_results[result_row, beta_col])
  se <- as.numeric(bmi_dsv_results[result_row, se_col])
  data.frame(
    study_name = study,
    beta = beta,
    SE = se,
    CI_low = beta - 1.96 * se,
    CI_high = beta + 1.96 * se,
    stringsAsFactors = FALSE
  )
})
beta_summary <- do.call(rbind, beta_rows)
beta_summary$study_label <- gsub("\n", " ", unname(study_labels[beta_summary$study_name]))
beta_summary$study_label <- factor(
  beta_summary$study_label,
  levels = rev(gsub("\n", " ", unname(study_labels[study_order])))
)
beta_summary$beta_CI_label <- sprintf(
  "%+.2f (%.2f, %.2f)",
  beta_summary$beta, beta_summary$CI_low, beta_summary$CI_high
)

forest_xlim <- range(c(beta_summary$CI_low, beta_summary$CI_high), na.rm = TRUE)
forest_span <- diff(forest_xlim)
if (!is.finite(forest_span) || forest_span == 0) forest_span <- 1
forest_xlim <- forest_xlim + c(-0.08, 0.08) * forest_span
label_x <- forest_xlim[2] + 0.12 * diff(forest_xlim)
forest_right <- label_x + 2 * diff(forest_xlim)
forest_ticks <- pretty(forest_xlim, n = 2)
forest_ticks <- forest_ticks[forest_ticks >= forest_xlim[1] & forest_ticks <= forest_xlim[2]]

p_forest <- ggplot(beta_summary, aes(y = study_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.8) +
  geom_segment(
    aes(x = CI_low, xend = CI_high, yend = study_label),
    linewidth = 1.1, color = "#3A3A3A"
  ) +
  geom_point(aes(x = beta), size = 3, color = "#2166AC") +
  geom_text(
    aes(label = beta_CI_label), x = label_x,
    hjust = 0, size = 5 + 2 / ggplot2::.pt, lineheight = 0.95, color = "black"
  ) +
  scale_x_continuous(limits = c(forest_xlim[1], forest_right),
                     breaks = forest_ticks) +
  labs(x = "Rank-normalized BMI (beta)", y = NULL) +
  theme_cowplot() +
  theme(
    panel.grid = element_blank(), axis.text = element_text(size = 18),
    axis.line = element_line(linewidth = 0.9), axis.ticks = element_line(linewidth = 0.8),
    axis.title = element_text(size = 20),
    plot.margin = unit(c(5, 12, 5, 5), "pt")
  )

p_whole <- patchwork::wrap_plots(
  c(plot_list, list(p_forest)),
  nrow = 1, widths = c(rep(1, length(plot_list)), 3.2)
)
# Increase the total width by one sixth and allocate more space to the forest.
f6d_width <- 26.6766666667 * 7 / 6
if (!dir.exists("pics/sites")) dir.create("pics/sites", recursive = TRUE)

write.csv(
  merge(delta_summary, beta_summary[, c("study_name", "beta", "SE", "CI_low", "CI_high")],
        by = "study_name", sort = FALSE),
  file = "pics/sites/F4_dSV_BMI_Anaerostipes_hadrus_1577_1579_effect_summary.csv",
  row.names = FALSE
)

tiff(
  file = "pics/F6D_dSV_BMI_Anaerostipes_hadrus_1577_1579_effects.tiff",
  width = f6d_width, height = 11, units = "in", res = 300, compression = "lzw"
)
print(p_whole)
dev.off()

############################################################################
# A. hadrus dSV and metabolic traits: include all QC-available cohorts and samples.
library(ggplot2)
library(patchwork)
library(cowplot)
library(ggbeeswarm)

dsv_target <- "Anaerostipes hadrus SSC/2:1577_1579 and 9 segments"
index_names <- c(fasting_glucose = "Fasting glucose", fasting_insulin = "Fasting insulin",
                 hba1c = "HbA1c", hdl = "HDL", triglycerides = "Triglycerides")
index_units <- c(fasting_glucose = "mmol/L", fasting_insulin = "mU/L",
                 hba1c = "%", hdl = "mmol/L", triglycerides = "mmol/L")
index_dsgv <- new.env()
load("01.cleanData/SV_all_QC/dsgv.RData", envir = index_dsgv)
index_dsgv <- index_dsgv$dsgv_sub
index_clin <- read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv", header = TRUE)
index_clin$match_id <- row.names(index_clin)
index_clin$study_name <- normalize_study_name(index_clin$study_name)
index_clin <- merge(index_clin, data.frame(match_id = row.names(index_dsgv), temp = index_dsgv[[dsv_target]]), by = "match_id")
index_colors <- setNames(mycolor2_green_blue[c(2, 1)], c("0", "1"))
dir.create("pics", recursive = TRUE, showWarnings = FALSE)

for (trait in names(index_names)) {
  # Keep every QC-observed phenotype-dSV pair, regardless of sample size or disease status.
  index_data <- index_clin[!is.na(index_clin$temp) & is.finite(index_clin[[trait]]), ]
  # Exclude LiJ_2014 only from the fasting-insulin figure as requested.
  if (trait == "fasting_insulin") index_data <- subset(index_data, study_name != "LiJ_2014")
  stopifnot(nrow(index_data) > 0L, all(index_data$temp %in% c(0, 1)))
  index_n <- sort(table(index_data$study_name), decreasing = TRUE)
  index_studies <- names(index_n)
  index_labels <- setNames(gsub("_", " ", sub("_([0-9]{4})", "\n\\1", index_studies)), index_studies)
  index_results <- read.csv(paste0("07.Microbial_GWAS/studies_adjAbun/", trait, "_dsv.csv"), check.names = FALSE)
  names(index_results) <- normalize_study_name(names(index_results))
  index_row <- which(index_results$X == dsv_target)
  stopifnot(length(index_row) == 1L)
  index_summary <- data.frame(study = index_studies, n_samples = as.integer(index_n))

  # Use the original adjusted effects; keep unestimated cohorts in the forest plot.
  for (field in c("beta", "SE", "N")) {
    columns <- paste0(index_studies, ".", field)
    available <- columns %in% names(index_results)
    index_summary[[field]] <- NA_real_
    index_summary[[field]][available] <- as.numeric(index_results[index_row, columns[available]])
  }
  index_summary$estimated <- is.finite(index_summary$beta) & is.finite(index_summary$SE) & index_summary$SE > 0
  index_summary$CI_low <- index_summary$beta - 1.96 * index_summary$SE
  index_summary$CI_high <- index_summary$beta + 1.96 * index_summary$SE
  index_breaks <- pretty(range(index_data[[trait]]), n = 5)
  index_plots <- vector("list", length(index_studies))

  for (i in seq_along(index_studies)) {
    study <- index_studies[i]
    x <- subset(index_data, study_name == study)
    n_group <- table(factor(x$temp, levels = c(0, 1)))
    delta <- if (all(n_group > 0)) median(x[[trait]][x$temp == 0]) - median(x[[trait]][x$temp == 1]) else NA_real_
    index_summary$n_deletion[i] <- n_group[1]
    index_summary$n_no_deletion[i] <- n_group[2]
    if (trait %in% c("hba1c", "fasting_insulin", "triglycerides", "fasting_glucose")) delta <- -delta  # No deletion minus deletion, consistent with the forest effect.
    index_summary$delta_median_no_deletion_minus_deletion[i] <- if (trait %in% c("hba1c", "fasting_insulin", "triglycerides", "fasting_glucose")) delta else -delta
    x$SV_group <- factor(x$temp, levels = c(0, 1),
                         labels = paste0(c("Deletion", "No deletion"), "\n(n=", n_group, ")"))
    index_plots[[i]] <- ggplot(x, aes(x = SV_group, y = .data[[trait]], color = as.character(temp))) +
      geom_hline(yintercept = index_breaks, linetype = "dashed", color = "grey65", linewidth = 0.8) +
      geom_violin(alpha = 0.05, width = 0.9, linewidth = 1.1, trim = TRUE) +
      geom_quasirandom(width = 0.18, size = 0.7, alpha = 0.55) +
      geom_boxplot(alpha = 0.18, width = 0.24, color = "black", linewidth = 0.9, outlier.shape = NA) +
      scale_color_manual(values = index_colors, guide = "none") +
      scale_x_discrete(drop = FALSE) +
      scale_y_continuous(limits = range(index_breaks), breaks = index_breaks) +
      labs(title = paste0(if (trait %in% c("hba1c", "fasting_insulin", "triglycerides", "fasting_glucose")) paste(strwrap(gsub("\n", " ", index_labels[study]), width = 14), collapse = "\n") else index_labels[study], "\nDelta = ", if (is.finite(delta)) sprintf("%+.2f", delta) else "NA"),
           x = NULL, y = if (i == 1L) paste0(index_names[trait], " (", index_units[trait], ")") else NULL) +
      theme_cowplot() +
      theme(panel.grid = element_blank(), legend.position = "none",
            axis.line = element_line(linewidth = 0.9), axis.ticks = element_line(linewidth = 0.8),
            text = element_text(size = if (trait %in% c("hba1c", "fasting_insulin", "triglycerides", "fasting_glucose")) 16 else 20), axis.text = element_text(size = if (trait %in% c("hba1c", "fasting_insulin", "triglycerides", "fasting_glucose")) 16 else 19), axis.title = element_text(size = if (trait %in% c("hba1c", "fasting_insulin", "triglycerides", "fasting_glucose")) 18 else 20),
            axis.text.y = if (i == 1L) element_text(size = if (trait %in% c("hba1c", "fasting_insulin", "triglycerides", "fasting_glucose")) 16 else 19) else element_blank(),
            axis.ticks.y = if (i == 1L) element_line() else element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            plot.title = element_text(size = if (trait %in% c("hba1c", "fasting_insulin", "triglycerides", "fasting_glucose")) 16 else 21, face = "plain", hjust = 0.5, vjust = 0,
                                      margin = margin(t = 12, b = 6)),
            plot.margin = margin(8, 8, 8, 8))
  }

  # Model sample sizes can differ from the full observed samples shown in the distributions.
  forest_labels <- paste0(gsub("_", " ", index_studies), "\n(model n=", ifelse(is.finite(index_summary$N), index_summary$N, "NA"), ")")
  index_summary$study_label <- factor(forest_labels, levels = rev(forest_labels))
  index_summary$beta_CI <- ifelse(index_summary$estimated,
                                 sprintf("%+.2f\n(%.2f, %.2f)", index_summary$beta, index_summary$CI_low, index_summary$CI_high),
                                 "Not estimated")
  index_xlim <- range(c(0, index_summary$CI_low[index_summary$estimated], index_summary$CI_high[index_summary$estimated]))
  index_span <- max(diff(index_xlim), 0.1)
  index_xlim <- index_xlim + c(-0.08, 0.08) * index_span
  index_label_x <- index_xlim[2] + 0.12 * diff(index_xlim)
  index_ticks <- if (abs(index_xlim[1]) > abs(index_xlim[2])) ceiling(index_xlim[1] * 10) / 10 else floor(index_xlim[2] * 10) / 10
  index_ticks <- sort(unique(c(0, index_ticks)))
  index_forest <- ggplot(index_summary, aes(y = study_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.8) +
    geom_segment(data = subset(index_summary, estimated), aes(x = CI_low, xend = CI_high, yend = study_label),
                 linewidth = 1.1, color = "#3A3A3A") +
    geom_point(data = subset(index_summary, estimated), aes(x = beta), size = 3, color = "#2166AC") +
    geom_text(aes(label = beta_CI), x = index_label_x, hjust = 0, size = if (trait %in% c("hba1c", "fasting_insulin", "triglycerides", "fasting_glucose")) 5 + 2 / ggplot2::.pt else 7.65, lineheight = 0.95, color = "black") +
    scale_y_discrete(drop = FALSE) +
    scale_x_continuous(limits = c(index_xlim[1], index_label_x + 2 * diff(index_xlim)),
                       breaks = index_ticks, labels = as.character) +
    labs(x = paste0("Beta for rank-normalized\n", index_names[trait]), y = NULL) +
    theme_cowplot() +
    theme(panel.grid = element_blank(), axis.text = element_text(size = if (trait %in% c("hba1c", "fasting_insulin", "triglycerides", "fasting_glucose")) 18 else 21.5),
          axis.line = element_line(linewidth = 0.9), axis.ticks = element_line(linewidth = 0.8),
          axis.title = element_text(size = 20), plot.margin = margin(5, 12, 5, 5))

  index_figure <- wrap_plots(c(index_plots, list(index_forest)), nrow = 1,
                             widths = c(rep(1, length(index_studies)), 1.6))
  # Display the original association P value, not a multiple-testing-adjusted P.
  index_meta <- read.csv(paste0("07.Microbial_GWAS/anno_adjAbun/", trait, "_sv.csv"), check.names = FALSE)
  index_meta <- index_meta[index_meta$X == dsv_target, , drop = FALSE]
  stopifnot(nrow(index_meta) == 1L)
  p_parts <- strsplit(sprintf("%.2e", index_meta$pvalue_meta), "e")[[1]]
  index_title <- bquote(beta[meta] == .(sprintf("%.2f", index_meta$effect_meta)) * "; " ~
    italic(P)[meta] == .(as.numeric(p_parts[1])) %*% 10^.(as.integer(p_parts[2])) ~ "(model 1)")
  index_figure <- index_figure + plot_annotation(title = index_title,
    theme = theme(plot.title = element_text(size = 22, hjust = 0, margin = margin(b = 15))))
  # Preserve cohort plot widths and reduce the forest plot width by half.
  grDevices::pdf(NULL)
  index_fixed_width <- sum(grid::convertWidth(patchwork::patchworkGrob(index_figure)$widths, "in", valueOnly = TRUE))
  grDevices::dev.off()
  index_width <- index_fixed_width + (2.7 * (length(index_studies) + 3.2) - index_fixed_width) *
    (length(index_studies) + 1.6) / (length(index_studies) + 3.2)
  ggsave(paste0("pics/FS7_dSV_", trait, "_Anaerostipes_hadrus_1577_1579_effects.tiff"),
         index_figure, device = grDevices::tiff, width = index_width, height = 8.5, units = "in",
         dpi = 300, compression = "lzw", type = "cairo")
}

############################################################################
# A. hadrus dSV and T2D/IBD: QC-based case/control status proportions.
library(ggplot2)
library(patchwork)
library(cowplot)

disease_dsgv <- new.env()
load("01.cleanData/SV_all_QC/dsgv.RData", envir = disease_dsgv)
disease_dsgv <- disease_dsgv$dsgv_sub
dsv_target <- "Anaerostipes hadrus SSC/2:1577_1579 and 9 segments"
stopifnot(dsv_target %in% names(disease_dsgv))
disease_clin <- read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv", header = TRUE)
disease_clin$match_id <- row.names(disease_clin)
disease_clin$study_name <- normalize_study_name(disease_clin$study_name)
disease_clin <- merge(disease_clin, data.frame(match_id = row.names(disease_dsgv),
                                             temp = disease_dsgv[[dsv_target]]), by = "match_id")
disease_colors <- setNames(mycolor2_green_blue[c(2, 1)], c("Deletion", "No deletion"))
dir.create("pics", recursive = TRUE, showWarnings = FALSE)

for (trait in c("T2D", "IBD", "ACVD")) {
  # Retain every QC-observed disease-dSV pair; do not filter by sample size or P value.
  disease_data <- disease_clin[!is.na(disease_clin$temp) & !is.na(disease_clin[[trait]]), ]
  stopifnot(nrow(disease_data) > 0L, all(disease_data$temp %in% c(0, 1)),
            all(disease_data[[trait]] %in% c(0, 1)))
  disease_n <- sort(table(disease_data$study_name), decreasing = TRUE)
  disease_studies <- names(disease_n)
  disease_labels <- setNames(gsub("_", " ", sub("_([0-9]{4})", "\n\\1", disease_studies)), disease_studies)
  disease_summary <- data.frame(study = disease_studies, n_samples = as.integer(disease_n))
  disease_plots <- vector("list", length(disease_studies))
  disease_counts <- vector("list", length(disease_studies))

  for (i in seq_along(disease_studies)) {
    study <- disease_studies[i]
    x <- subset(disease_data, study_name == study)
    counts <- as.data.frame(table(status = factor(x[[trait]], levels = c(0, 1)),
                                  SV = factor(x$temp, levels = c(0, 1), labels = c("Deletion", "No deletion"))))
    counts$total <- ave(counts$Freq, counts$status, FUN = sum)
    counts$proportion <- ifelse(counts$total > 0, counts$Freq / counts$total, NA_real_)
    counts$label <- ifelse(counts$Freq > 0, sprintf("%.1f%%", 100 * counts$proportion), "")
    counts$study <- study
    disease_counts[[i]] <- counts
    group_n <- as.integer(table(factor(x[[trait]], levels = c(0, 1))))
    disease_summary$n_control[i] <- group_n[1]
    disease_summary$n_case[i] <- group_n[2]
    disease_plots[[i]] <- ggplot(counts, aes(status, proportion, fill = SV)) +
      geom_col(width = 0.85, position = position_stack(reverse = TRUE), color = "black", linewidth = 0.6, na.rm = TRUE) +
      geom_text(aes(label = label), position = position_stack(vjust = 0.5, reverse = TRUE),
                size = 5, color = "white", na.rm = TRUE) +
      scale_fill_manual(values = disease_colors, drop = FALSE) +
      scale_x_discrete(labels = paste0(c("Control", trait), "\n(n=", group_n, ")"), drop = FALSE) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25), expand = c(0, 0)) +
      labs(title = disease_labels[study], x = NULL, y = if (i == 1L) "dSV status (%)" else NULL, fill = NULL) +
      theme_cowplot() +
      theme(text = element_text(family = "sans", size = 20), axis.text = element_text(size = 19),
            axis.title = element_text(size = 20), axis.line = element_line(linewidth = 0.9),
            axis.ticks = element_line(linewidth = 0.8),
            axis.text.y = if (i == 1L) element_text(size = 19) else element_blank(),
            axis.ticks.y = if (i == 1L) element_line() else element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            plot.title = element_text(size = 21, face = "plain", hjust = 0.5, margin = margin(t = 12, b = 6)),
            plot.margin = margin(2, 2, 2, 2))
  }
  disease_counts <- do.call(rbind, disease_counts)

  disease_figure <- wrap_plots(disease_plots, nrow = 1, guides = "collect") &
    theme(legend.position = "bottom", legend.text = element_text(size = 19))
  disease_model <- if (trait == "ACVD") "anno_adjAbunPCs" else "anno_adjAbun"
  disease_result <- read.csv(paste0("07.Microbial_GWAS/", disease_model, "/", trait, "_sv.csv"), check.names = FALSE)
  disease_result <- disease_result[disease_result$X == dsv_target, , drop = FALSE]
  stopifnot(nrow(disease_result) == 1L)
  disease_parts <- strsplit(sprintf("%.2e", disease_result$pvalue_meta), "e")[[1]]
  model_label <- if (trait == "ACVD") "(model 2)" else "(model 1)"
  disease_or_cols <- grep("\\.OR$", names(disease_result), value = TRUE)
  disease_or_values <- as.numeric(disease_result[1, disease_or_cols])
  disease_n_studies <- sum(is.finite(disease_or_values))
  if (disease_n_studies > 1L && "OR_meta" %in% names(disease_result) && is.finite(disease_result$OR_meta)) {
    disease_title <- bquote(OR[meta] == .(sprintf("%.2f", disease_result$OR_meta)) * "; " ~
      italic(P)[meta] == .(as.numeric(disease_parts[1])) %*% 10^.(as.integer(disease_parts[2])) ~ .(model_label))
  } else {
    or_columns <- grep("\\.OR$", names(disease_result), value = TRUE)
    or_values <- as.numeric(disease_result[1, or_columns])
    stopifnot(sum(is.finite(or_values)) == 1L)
    disease_title <- bquote(OR == .(sprintf("%.2f", or_values[is.finite(or_values)])) * "; " ~
      italic(P) == .(as.numeric(disease_parts[1])) %*% 10^.(as.integer(disease_parts[2])) ~ .(model_label))
  }
  disease_figure <- disease_figure + plot_annotation(title = disease_title,
    theme = theme(plot.title = element_text(size = 20, hjust = 0, margin = margin(b = 15))))
  disease_width <- max(if (trait == "ACVD") 7 else 5.4, 3 * length(disease_studies) + 1)
  if (trait == "T2D") disease_width <- disease_width * 2 / 3
  ggsave(paste0("pics/FS7_dSV_", trait, "_Anaerostipes_hadrus_1577_1579_effects.tiff"),
         disease_figure, device = grDevices::tiff, width = disease_width, height = 7.5, units = "in",
         dpi = 300, compression = "lzw", type = "cairo")
}

############################################################################
# A. hadrus dSV and age: include every study with more than 100 valid samples.
library(ggplot2)

dsv_target <- "Anaerostipes hadrus SSC/2:1577_1579 and 9 segments"
age_data <- data.frame(match_id = row.names(dsgv), temp = dsgv[[dsv_target]])
age_data <- merge(all_basic, age_data, by = "match_id")
age_data <- subset(age_data, is.finite(age) & age >= 18 & !is.na(temp))
age_data$study_name <- normalize_study_name(age_data$study_name)
stopifnot(all(age_data$temp %in% c(0, 1)))

# Retain the category codes expected by dsv_age; axis labels show the actual intervals.
age_data$age_cat <- cut(age_data$age, breaks = c(18, 35, 60, Inf), right = FALSE,
                        labels = c("(18,35]", "(35,60]", "[60,+Inf)"))
study_n <- sort(table(age_data$study_name), decreasing = TRUE)
study_n <- study_n[study_n > 100L]
if (!length(study_n)) stop("No studies have more than 100 valid dSV-age samples.")
age_data <- subset(age_data, study_name %in% names(study_n))
ratio_pic <- do.call(rbind, lapply(names(study_n), function(study) {
  dsv_age(age_data[age_data$study_name == study, ], study)
}))
ratio_pic$study <- factor(ratio_pic$study, levels = names(study_n))
study_labels <- setNames(paste0(names(study_n), "\n(n=", as.integer(study_n), ")"), names(study_n))

# Preserve the existing study colors and assign unused predefined colors to new studies.
study_palette <- setNames(mycolor11[c(6, 12, 1, 4, 5, 2, 7, 8)],
                          c("AsnicarF_2021", "Fern\u00e1ndez_Calvet_2024", "HansenLBS_2018",
                            "SchirmerM_2016", "XieH_2016", "YachidaS_2019", "JieZ_2017", "ZellerG_2014"))
extra_studies <- setdiff(names(study_n), names(study_palette))
study_palette <- c(study_palette, setNames(rep(setdiff(mycolor15, study_palette),
                                              length.out = length(extra_studies)), extra_studies))

p_age <- ggplot(ratio_pic, aes(Qt, value, color = study, group = study)) +
  geom_hline(yintercept = c(0.25, 0.5, 0.75), linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 1:3, linetype = "dotted", color = "gray80") +
  geom_line(linewidth = 0.8, na.rm = TRUE) + geom_point(size = 2.2, na.rm = TRUE) +
  scale_x_continuous(breaks = 1:3, labels = c("[18,35)", "[35,60)", "\u226560")) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25), limits = c(0, 1)) +
  scale_color_manual(values = study_palette, breaks = names(study_n), labels = study_labels) +
  labs(x = "Age (years)", y = "Deletion rate (%)", color = "Study") +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(axis.text = element_text(size = 11, color = "black"),
        legend.position = "right", legend.text = element_text(size = 11),
        legend.key.height = grid::unit(0.9, "cm"), plot.margin = margin(6, 8, 6, 6)) +
  guides(color = guide_legend(ncol = 1))

dir.create("pics", recursive = TRUE, showWarnings = FALSE)
ggsave("pics/FS7_dSV_age_Anaerostipes_hadrus_1577_1579.tiff", p_age, device = grDevices::tiff,
       width = 7, height = 4.2, units = "in", dpi = 300, compression = "lzw", type = "cairo")

############################################################################
# A. hadrus

library(data.table)
library(ComplexHeatmap)
library(grid)

mycolor2_blue_yellow <- c("#2983BB", "#ECD452")
input_dir <- "00.rawData/A.H"
dir.create("pics", recursive = TRUE, showWarnings = FALSE)

############################################################################
# hypothetical protein

targets <- as.data.frame(fread(file.path(input_dir, "dSV_gene.csv"), check.names = FALSE))
stopifnot(all(c("order", "gene", "product") %in% names(targets)))

targets$source_row <- seq_len(nrow(targets))
targets$order <- as.numeric(targets$order)
stopifnot(!anyNA(targets$order), !anyNA(targets$product))

targets <- targets[order(targets$order, targets$source_row), , drop = FALSE]
excluded <- tolower(trimws(targets$product)) == "hypothetical protein"
targets <- targets[!excluded, , drop = FALSE]

if (nrow(targets) == 0L) stop("No target genes remain after excluding hypothetical protein.")
if (any(!nzchar(trimws(targets$product)))) stop("Empty target product.")

############################################################################
# 

presence <- fread(file.path(input_dir, "gene_presence_absence.csv"), colClasses = "character", check.names = FALSE)
rtab <- fread(file.path(input_dir, "gene_presence_absence.Rtab"), check.names = FALSE)
stopifnot(all(c("Gene", "Annotation") %in% names(presence)), !anyDuplicated(presence$Gene), !anyDuplicated(rtab[[1]]))

samples <- names(rtab)[-1]
stopifnot(length(samples) >= 2L, !anyDuplicated(samples), all(samples %in% names(presence)))

rtab_rows <- match(presence$Gene, rtab[[1]])
if (anyNA(rtab_rows)) stop("CSV gene clusters are missing from the Rtab file.")

############################################################################
#

terms <- strsplit(ifelse(is.na(presence$Annotation), "", presence$Annotation), ";", fixed = TRUE)
term_table <- data.table(cluster_row = rep(seq_len(nrow(presence)), lengths(terms)), term = unlist(terms, use.names = FALSE))
term_table$key <- tolower(trimws(gsub("[[:space:]]+", " ", gsub(",", "", term_table$term, fixed = TRUE))))
target_key <- tolower(trimws(gsub("[[:space:]]+", " ", gsub(",", "", targets$product, fixed = TRUE))))

matches <- vector("list", nrow(targets))
for (i in seq_len(nrow(targets))) {
  matches[[i]] <- unique(term_table$cluster_row[term_table$key == target_key[i]])
  if (length(matches[[i]]) == 0L) stop("No annotation match for order ", targets$order[i], ": ", targets$product[i])
}

############################################################################
#

selected_rows <- sort(unique(unlist(matches)))
binary_clusters <- matrix(0L, nrow = length(selected_rows), ncol = length(samples), dimnames = list(presence$Gene[selected_rows], samples))
discordance <- data.frame()

for (j in seq_along(samples)) {
  cells <- trimws(presence[[samples[j]]][selected_rows])
  csv_binary <- as.integer(!is.na(cells) & nzchar(cells) & !toupper(cells) %in% c("NA", "N/A", "NULL", "0"))
  rtab_values <- as.numeric(rtab[[samples[j]]][rtab_rows[selected_rows]])
  if (anyNA(rtab_values) || any(!is.finite(rtab_values)) || any(rtab_values < 0)) stop("Invalid Rtab values: ", samples[j])

  rtab_binary <- as.integer(rtab_values > 0)
  bad <- which(csv_binary != rtab_binary)
  if (length(bad) > 0L) {
    discordance <- rbind(discordance, data.frame(cluster_id = presence$Gene[selected_rows[bad]], sample = samples[j], csv = csv_binary[bad], rtab = rtab_binary[bad]))
  }
  binary_clusters[, j] <- csv_binary
}

if (nrow(discordance) > 0L) {
  print(head(discordance, 20))
  stop("CSV and Rtab presence/absence values disagree.")
}

############################################################################
#

presence_matrix <- matrix(0L, nrow = nrow(targets), ncol = length(samples))
for (i in seq_len(nrow(targets))) {
  ix <- match(matches[[i]], selected_rows)
  presence_matrix[i, ] <- as.integer(colSums(binary_clusters[ix, , drop = FALSE]) > 0)
}

row_id <- paste0("order_", targets$order, "_row_", targets$source_row)
rownames(presence_matrix) <- row_id
colnames(presence_matrix) <- samples
stopifnot(all(presence_matrix %in% c(0L, 1L)))

############################################################################
# 

sample_distance <- dist(t(presence_matrix), method = "binary")
stopifnot(all(is.finite(sample_distance)))
sample_cluster <- hclust(sample_distance, method = "average")

############################################################################
#

gene_symbols <- trimws(ifelse(is.na(targets$gene), "", targets$gene))
row_labels <- ifelse(nzchar(gene_symbols), paste(gene_symbols, targets$product, sep = " | "), targets$product)
for (i in seq_along(row_labels)) row_labels[i] <- paste(strwrap(row_labels[i], width = 67), collapse = "\n")

heatmap <- Heatmap(
  presence_matrix, name = "Presence", col = setNames(mycolor2_blue_yellow, c("0", "1")),
  cluster_rows = FALSE, cluster_columns = sample_cluster, column_dend_reorder = FALSE,
  show_row_names = TRUE, row_labels = row_labels, row_names_side = "left", row_names_gp = gpar(fontfamily = "Arial", fontsize = 16, col = "black"), row_names_max_width = unit(170, "mm"),
  show_column_names = FALSE, column_dend_height = unit(22, "mm"), rect_gp = gpar(col = NA),
  use_raster = TRUE, raster_device = "png", raster_quality = 3,
  heatmap_legend_param = list(at = c(0, 1), labels = c("absent", "present"), title_gp = gpar(fontfamily = "Arial", fontsize = 19, fontface = "bold"), labels_gp = gpar(fontfamily = "Arial", fontsize = 18)),
  border = TRUE
)

############################################################################
# 

tiff(filename = "pics/F6C_A.H.heatmap.tiff", width = 20, height = 14 * 2 / 3, units = "in", res = 300, compression = "lzw")

drawn <- draw(heatmap, heatmap_legend_side = "right", padding = unit(c(5, 5, 5, 5), "mm"))
row_order_valid <- identical(as.integer(row_order(drawn)), seq_len(nrow(targets)))
column_order_valid <- identical(as.integer(column_order(drawn)), as.integer(sample_cluster$order))

dev.off()

stopifnot(row_order_valid, column_order_valid)
cat("AH_HEATMAP_VALIDATED\nRows:", nrow(targets), "; MAGs:", length(samples), "; matched clusters:", length(selected_rows), "; CSV/Rtab discrepancies:", nrow(discordance), "\n")

############################################################################
# A. hadrus 

library(ggplot2)

result_file <- "11.gene_function/kegg/results/SV/BMI_dsv_Anaerostipes_hadrus_1577_1610.csv"
annotation_file <- "11.gene_function/kegg/SV_kegg.csv"
symbol_file <- "00.rawData/A.H/dSV_gene.csv"
dir.create("pics", recursive = TRUE, showWarnings = FALSE)

############################################################################
# 1.

enrichment <- read.csv(result_file, stringsAsFactors = FALSE, check.names = FALSE)
sv_annotation <- read.csv(annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
# Standardize the strain name in the archived KEGG annotation.
sv_annotation$SV <- sub("^Anaerostipes hadrus [^:]+:", "Anaerostipes hadrus SSC/2:", sv_annotation$SV)
gene_symbols <- read.csv(symbol_file, stringsAsFactors = FALSE)

stopifnot(all(c("ID", "Description", "pvalue", "p.adjust", "geneID", "Count", "GeneRatio") %in% names(enrichment)))
stopifnot(nrow(enrichment) > 0L, !anyDuplicated(enrichment$ID), all(is.finite(enrichment$p.adjust)), all(enrichment$p.adjust > 0 & enrichment$p.adjust <= 1))
stopifnot(all(c("phenotype", "SV", "genes", "ORTHOLOGY") %in% names(sv_annotation)), all(c("gene", "ko_gene") %in% names(gene_symbols)))

enrichment <- enrichment[order(enrichment$p.adjust, enrichment$pvalue), , drop = FALSE]
enrichment$ratio_denominator <- NA_integer_

############################################################################
# 2.

pathway_KO <- data.frame(ID = character(), KO = character(), stringsAsFactors = FALSE)

for (i in seq_len(nrow(enrichment))) {
  ids <- unique(strsplit(enrichment$geneID[i], "/", fixed = TRUE)[[1]])
  parts <- regmatches(enrichment$GeneRatio[i], gregexpr("[0-9]+", enrichment$GeneRatio[i]))[[1]]

  stopifnot(all(grepl("^K[0-9]{5}$", ids)), length(ids) == enrichment$Count[i])
  stopifnot(length(parts) == 2L, as.integer(parts[1]) == enrichment$Count[i])

  enrichment$ratio_denominator[i] <- as.integer(parts[2])
  pathway_KO <- rbind(pathway_KO, data.frame(ID = enrichment$ID[i], KO = ids))
}

stopifnot(length(unique(enrichment$ratio_denominator)) == 1L)
ko_denominator <- unique(enrichment$ratio_denominator)
selected_KOs <- unique(pathway_KO$KO)

############################################################################
# 3.

sv_genes <- unique(sv_annotation[sv_annotation$phenotype == "BMI_sv" & grepl("^Anaerostipes hadrus SSC/2:1577_", sv_annotation$SV), c("SV", "genes", "ORTHOLOGY")])
stopifnot(nrow(sv_genes) > 0L, length(unique(sv_genes$SV)) == 1L, !anyDuplicated(sv_genes$genes))

gene_KO <- data.frame(genes = character(), KO = character(), stringsAsFactors = FALSE)

for (i in seq_len(nrow(sv_genes))) {
  ids <- regmatches(sv_genes$ORTHOLOGY[i], gregexpr("K[0-9]{5}", sv_genes$ORTHOLOGY[i]))[[1]]
  if (length(ids) > 0L) gene_KO <- rbind(gene_KO, data.frame(genes = sv_genes$genes[i], KO = unique(ids)))
}

stopifnot(all(selected_KOs %in% gene_KO$KO))
links <- unique(merge(pathway_KO, gene_KO, by = "KO"))

############################################################################
#

displayed_genes <- unique(gene_KO$genes[gene_KO$KO %in% selected_KOs])
displayed_genes <- displayed_genes[order(as.integer(sub(".*_", "", displayed_genes)))]

gene_key <- data.frame(genes = displayed_genes, symbol = "", stringsAsFactors = FALSE)

for (i in seq_len(nrow(gene_key))) {
  ids <- unique(gene_KO$KO[gene_KO$genes == gene_key$genes[i]])
  symbol <- unique(gene_symbols$gene[gene_symbols$ko_gene %in% ids & !is.na(gene_symbols$gene) & nzchar(gene_symbols$gene)])

  gene_key$symbol[i] <- if (length(symbol) == 1L) symbol else ""
}

stopifnot(all(nzchar(gene_key$symbol)))

############################################################################
#

membership <- expand.grid(ID = enrichment$ID, genes = gene_key$genes, stringsAsFactors = FALSE)
membership$linked <- paste(membership$ID, membership$genes) %in% paste(links$ID, links$genes)

enrichment$row <- rev(seq_len(nrow(enrichment)))
enrichment$adjusted_P_label <- formatC(enrichment$p.adjust, digits = 3, format = "g")
pathway_labels <- vapply(enrichment$Description, function(x) paste(strwrap(x, width = 14), collapse = "\n"), character(1))
membership$row <- enrichment$row[match(membership$ID, enrichment$ID)]
membership$column <- match(membership$genes, gene_key$genes)

for (i in seq_len(nrow(enrichment))) {
  observed <- unique(links$KO[links$ID == enrichment$ID[i]])
  stopifnot(length(observed) == enrichment$Count[i])
}

############################################################################
# 6. Show enrichment and gene-pathway annotation in one matrix
############################################################################

blue <- "#2983BB"
label_size <- 3.55 + 2 / ggplot2::.pt

figure <- ggplot(membership, aes(x = column, y = row)) +
  geom_tile(fill = "grey96", colour = "white", linewidth = 0.5, width = 0.94, height = 0.82) +
  geom_point(data = membership[membership$linked, , drop = FALSE], colour = blue, size = 2.4) +
  geom_text(data = enrichment, aes(x = -2.1, y = row, label = adjusted_P_label), inherit.aes = FALSE, size = label_size) +
  geom_text(data = enrichment, aes(x = -0.5, y = row, label = Count), inherit.aes = FALSE, size = label_size) +
  annotate("text", x = -2.1, y = nrow(enrichment) + 0.63, label = "Adj. P", size = label_size, fontface = "bold") +
  annotate("text", x = -0.5, y = nrow(enrichment) + 0.63, label = "KOs", size = label_size, fontface = "bold") +
  annotate("text", x = 6.5, y = nrow(enrichment) + 0.63, label = "Gene-pathway links", size = label_size, fontface = "bold") +
  geom_vline(xintercept = 0.25, colour = "grey75", linewidth = 0.35) +
  scale_x_continuous(breaks = seq_along(displayed_genes), labels = gene_key$symbol,
                     limits = c(-2.85, length(displayed_genes) + 0.55), expand = c(0, 0)) +
  scale_y_continuous(breaks = enrichment$row, labels = pathway_labels,
                     limits = c(0.45, nrow(enrichment) + 0.85), expand = c(0, 0)) +
  labs(x = "Gene", y = NULL) +
  theme_classic(base_size = 14, base_family = "sans") +
  theme(axis.text.y = element_text(size = 12.5, colour = "black", lineheight = 0.9),
        axis.text.x = element_text(size = 12, colour = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_blank(), axis.line = element_blank(), axis.title.x = element_text(size = 13.5, margin = margin(t = 6)),
        plot.margin = margin(6, 6, 6, 4))

ggsave("pics/F6_Anaerostipes_hadrus_KEGG_panels.tiff", figure, device = "tiff",
       width = 8.4, height = 4.6, units = "in", dpi = 300, compression = "lzw")

cat("AH_KEGG_PANELS_COMPLETE", nrow(enrichment), "pathways;", length(selected_KOs), "KOs;", length(displayed_genes), "displayed loci from", nrow(sv_genes), "annotated loci\n")





############################################################################
## Gene abundance by A. hadrus dSV status (0=deletion; 1=no deletion).
## Two-sided tests; positive beta supports higher abundance in state 1.
library(data.table)
library(metafor)
ah_out <- "12.AH_dSV_gene_abundance"
dir.create(ah_out, showWarnings=FALSE)
ah_target <- "Anaerostipes hadrus SSC/2:1577_1579 and 9 segments"
ah_clin <- read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv",header=TRUE)
ah_env <- new.env()
load("01.cleanData/SV_all_QC/dsgv.RData",envir=ah_env)
ah_ids <- intersect(rownames(ah_clin),rownames(ah_env$dsgv_sub))
ah_meta <- data.frame(sample=ah_ids,study=ah_clin[ah_ids,"study_name"],
  status=ah_env$dsgv_sub[ah_ids,ah_target])
ah_meta <- ah_meta[!is.na(ah_meta$status),]
stopifnot(!anyDuplicated(ah_meta$sample),all(ah_meta$status %in% c(0,1)))
ah_sv <- fread("01.cleanData/SV_info_QC/dsgv_info_anno.tsv")
ah_sv <- ah_sv[SV_Name == ah_target]
stopifnot(nrow(ah_sv)==1L)
ah_genes <- unique(strsplit(ah_sv$Genes,";",fixed=TRUE)[[1]])
ah_mapping <- fread("00.rawData/humann3/match_file/SGV_Finder_genes_vs_uniref90_annotation.tsv")
ah_mapping <- ah_mapping[gene %in% ah_genes]
fwrite(ah_mapping,file.path(ah_out,"HUMAnN3_gene_mapping.csv"))
ah_short_annotation <- fread("00.rawData/A.H/dSV_gene.csv")
ah_short_ids <- paste0(ah_short_annotation$locas,"_A_hadrus")
ah_sources <- c(ShortBRED="00.rawData/shortBRED/sig_SVs/matched_rows.csv",
  HUMAnN3="00.rawData/humann3/sig_SV/SV_no_species_merged_genefamilies_matched.tsv")
# HUMAnN3 is deferred until conflicting duplicate sample columns are resolved.
ah_sources <- ah_sources["ShortBRED"]
ah_rows <- list()
ah_coverage <- list()
for (method in names(ah_sources)) {
  raw <- fread(ah_sources[[method]])
  names(raw) <- sub("_kneaddata-RPKs$","",names(raw))
  for (dup in unique(names(raw)[duplicated(names(raw))])) {
    cols <- which(names(raw)==dup)
    stopifnot(all(vapply(cols[-1],function(k) identical(raw[[cols[1]]],raw[[k]]),logical(1))))
  }
  raw <- raw[,which(!duplicated(names(raw))),with=FALSE]
  feature_ids <- if (method=="ShortBRED") ah_short_ids else unique(ah_mapping$UniRef90_id)
  raw <- raw[get(names(raw)[1]) %in% feature_ids]
  stopifnot(nrow(raw)>0,!anyDuplicated(raw[[1]]))
  common <- intersect(ah_meta$sample,names(raw)[-1])
  ah_coverage[[method]] <- data.frame(method=method,requested_features=length(feature_ids),
    measured_features=nrow(raw),matched_samples=length(common))
  dat <- ah_meta[match(common,ah_meta$sample),]
  for (j in seq_len(nrow(raw))) {
    dat$abundance <- as.numeric(raw[j,..common])
    for (study_id in unique(dat$study)) {
      d <- dat[dat$study==study_id & is.finite(dat$abundance),]
      n0 <- sum(d$status==0); n1 <- sum(d$status==1)
      row <- data.frame(method=method,feature=raw[[1]][j],study=study_id,n0=n0,n1=n1,
        median0=if(n0) median(d$abundance[d$status==0]) else NA_real_,
        median1=if(n1) median(d$abundance[d$status==1]) else NA_real_,
        beta=NA_real_,SE=NA_real_,P=NA_real_,reason="Insufficient samples")
      if(nrow(d)>=40 && min(n0,n1)>=10 && length(unique(d$abundance))>1) {
        d$abundance_norm <- qnorm((rank(d$abundance,ties.method="average")-0.5)/nrow(d))
        fit <- summary(lm(abundance_norm ~ status,data=d))$coefficients
        row$beta <- fit["status","Estimate"]; row$SE <- fit["status","Std. Error"]
        row$P <- fit["status","Pr(>|t|)"]; row$reason <- "Estimated"
      } else if(length(unique(d$abundance))<=1) row$reason <- "Constant abundance"
      ah_rows[[length(ah_rows)+1L]] <- row
    }
  }
  rm(raw);gc()
}
ah_study <- rbindlist(ah_rows)
fwrite(ah_study,file.path(ah_out,"study_results.csv"))
fwrite(rbindlist(ah_coverage),file.path(ah_out,"input_coverage.csv"))
ah_pool <- list()
for (method in names(ah_sources)) {
  for (feature_id in unique(ah_study$feature[ah_study$method==method])) {
    d <- ah_study[ah_study$method==method & feature==feature_id & is.finite(beta) & is.finite(SE) & SE>0]
    if(!nrow(d)) next
    if(nrow(d)>1) {
      fit <- rma(yi=beta,sei=SE,data=d,method="REML")
      beta <- as.numeric(fit$b); se <- fit$se; p <- fit$pval; i2 <- fit$I2
    } else { beta <- d$beta; se <- d$SE; p <- d$P; i2 <- NA_real_ }
    ah_pool[[length(ah_pool)+1L]] <- data.frame(method=method,feature=feature_id,
      studies=nrow(d),n0=sum(d$n0),n1=sum(d$n1),beta=beta,SE=se,
      CI_low=beta-1.96*se,CI_high=beta+1.96*se,P=p,I2=i2)
  }
}
ah_pooled <- rbindlist(ah_pool)
ah_pooled[,Bonferroni_P:=p.adjust(P,"bonferroni"),by=method]
ah_pooled[,higher_in_no_deletion:=beta>0 & Bonferroni_P<0.05]
fwrite(ah_pooled,file.path(ah_out,"meta_results.csv"))
writeLines(c("All QC samples with callable dSV and measured abundance; no healthy-only restriction.",
  "Each study: >=40 samples and >=10 per dSV group; zeros retained, missing values excluded.",
  "Within-study rank-normalized abundance ~ dSV status; unadjusted two-sided linear model.",
  "Positive beta: no deletion (1) has higher abundance than deletion (0).",
  "Random-effects REML meta-analysis for >=2 studies; single-study effects reported separately.",
  "Bonferroni correction across tested features separately for each profiling method.",
  "HUMAnN3 families can include homologous genes from other taxa/regions; not locus-specific abundance.",
  "Associations may reflect source-species abundance and other confounders."),
  file.path(ah_out,"analysis_notes.txt"))
cat("AH_GENE_ABUNDANCE_COMPLETE\n")

############################################################################
## Forest plot: ShortBRED gene abundance by A. hadrus deletion status.
library(ggplot2)
library(patchwork)
ah_forest <- read.csv("12.AH_dSV_gene_abundance/meta_results.csv")
ah_forest <- ah_forest[ah_forest$method=="ShortBRED",]
ah_labels <- read.csv("00.rawData/A.H/dSV_gene.csv",stringsAsFactors=FALSE)
ah_labels$feature <- paste0(ah_labels$locas,"_A_hadrus")
stopifnot(!anyDuplicated(ah_labels$feature))
ah_match <- match(ah_forest$feature,ah_labels$feature)
stopifnot(!anyNA(ah_match),nrow(ah_forest)==27L)
ah_forest$gene <- ah_labels$gene[ah_match]
ah_forest$product <- ah_labels$product[ah_match]
ah_forest$locus <- ah_labels$locas[ah_match]
ah_forest$order <- ah_labels$order[ah_match]
ah_forest <- ah_forest[order(ah_forest$order),]
ah_forest$label <- ifelse(is.na(ah_forest$gene) | !nzchar(ah_forest$gene),
  ah_forest$locus,paste0(ah_forest$gene," (",ah_forest$locus,")"))
ah_forest$label <- paste0(ah_forest$label," | ",ah_forest$product)
ah_forest$label <- vapply(ah_forest$label,function(x) paste(strwrap(x,width=67),collapse="\n"),character(1))
ah_forest$row <- rev(seq_len(nrow(ah_forest)))
ah_forest$effect_label <- sprintf("%.2f (%.2f, %.2f)",ah_forest$beta,ah_forest$CI_low,ah_forest$CI_high)
ah_forest$p_label <- formatC(ah_forest$P,format="e",digits=2)
ah_ylim <- c(0.4,nrow(ah_forest)+1.6)
ah_table_theme <- theme_void(base_size=16) + theme(plot.margin=margin(5,5,5,5))
ah_name_panel <- ggplot(ah_forest,aes(y=row)) +
  geom_text(aes(x=0,label=label),hjust=0,size=3.5+4/ggplot2::.pt,lineheight=0.95) +
  annotate("text",x=0,y=28.2,label="Gene",hjust=0,fontface="bold",size=4.2+4/ggplot2::.pt) +
  scale_x_continuous(limits=c(0,1)) + coord_cartesian(ylim=ah_ylim,clip="off") + ah_table_theme
ah_effect_panel <- ggplot(ah_forest,aes(y=row,x=beta)) +
  geom_vline(xintercept=0,linetype="dashed",color="grey60") +
  geom_segment(aes(x=CI_low,xend=CI_high,yend=row),color="grey30",linewidth=0.6) +
  geom_point(color="#2166AC",size=2.4) +
  scale_y_continuous(limits=ah_ylim) +
  scale_x_continuous(limits=range(c(0,ah_forest$CI_high))*c(1,1.08),n.breaks=4) +
  labs(x=expression(beta~"(no deletion vs deletion)"),y=NULL) +
  theme_classic(base_size=16) +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
        plot.margin=margin(5,10,5,10))
ah_stats_panel <- ggplot(ah_forest,aes(y=row)) +
  geom_text(aes(x=0,label=effect_label),hjust=0,size=3.7+4/ggplot2::.pt) +
  geom_text(aes(x=0.64,label=p_label),hjust=0,size=3.7+4/ggplot2::.pt) +
  annotate("text",x=0,y=28.2,label="Beta (95% CI)",hjust=0,fontface="bold",size=4.2+4/ggplot2::.pt) +
  annotate("text",x=0.64,y=28.2,label="Meta P",hjust=0,fontface="bold",size=4.2+4/ggplot2::.pt) +
  scale_x_continuous(limits=c(0,1.12)) + coord_cartesian(ylim=ah_ylim,clip="off") + ah_table_theme
ah_forest_plot <- (ah_name_panel + ah_effect_panel + ah_stats_panel) +
  plot_layout(widths=c(5.5,2.5,3.7)) +
  plot_annotation(title=expression(italic("Anaerostipes hadrus")~"dSV: gene abundance by deletion status"),
    theme=theme(plot.title=element_text(size=20,hjust=0)))
ggsave("pics/FS8_AH_dSV_ShortBRED_abundance_forest.tiff",ah_forest_plot,
  width=17,height=13,units="in",dpi=300,compression="lzw",device="tiff")
ggsave("validation/AH_abundance_forest_preview.png",ah_forest_plot,width=17,height=13,dpi=100)
write.csv(ah_forest,"12.AH_dSV_gene_abundance/forest_source.csv",row.names=FALSE)

############################################################################
## carA ShortBRED abundance: dSV groups and study-level forest plot.
library(data.table)
library(ggplot2)
library(patchwork)
ca_id <- "CL2_16170_A_hadrus"
ca_raw <- fread("00.rawData/shortBRED/sig_SVs/matched_rows.csv")
ca_raw <- ca_raw[Family==ca_id]
stopifnot(nrow(ca_raw)==1L,!anyDuplicated(names(ca_raw)))
ca_clin <- read.table("01.cleanData/phen_QC/Clinical_basic_overlap.tsv",header=TRUE)
ca_env <- new.env()
load("01.cleanData/SV_all_QC/dsgv.RData",envir=ca_env)
ca_ids <- Reduce(intersect,list(rownames(ca_clin),rownames(ca_env$dsgv_sub),names(ca_raw)[-1]))
ca_data <- data.frame(sample=ca_ids,study=ca_clin[ca_ids,"study_name"],
  status=ca_env$dsgv_sub[ca_ids,"Anaerostipes hadrus SSC/2:1577_1579 and 9 segments"],
  abundance=as.numeric(ca_raw[,..ca_ids]))
ca_data <- ca_data[!is.na(ca_data$status) & is.finite(ca_data$abundance),]
ca_study <- read.csv("12.AH_dSV_gene_abundance/study_results.csv")
ca_study <- subset(ca_study,method=="ShortBRED" & feature==ca_id & reason=="Estimated")
ca_study <- ca_study[order(-(ca_study$n0+ca_study$n1)),]
ca_meta <- read.csv("12.AH_dSV_gene_abundance/meta_results.csv")
ca_meta <- subset(ca_meta,method=="ShortBRED" & feature==ca_id)
stopifnot(nrow(ca_meta)==1L,nrow(ca_study)>0L)
ca_plots <- list()
for(i in seq_len(nrow(ca_study))) {
  d <- ca_data[ca_data$study==ca_study$study[i],]
  d$group <- factor(d$status,levels=c(0,1),labels=c("Deletion","No deletion"))
  delta <- median(d$abundance[d$status==1])-median(d$abundance[d$status==0])
  ca_plots[[i]] <- ggplot(d,aes(group,abundance,color=group)) +
    geom_violin(trim=TRUE,fill=NA,linewidth=0.5) +
    geom_boxplot(width=0.18,outlier.shape=NA,color="black",fill="white",linewidth=0.4) +
    scale_color_manual(values=c("#8B83BE","#207F4C"),guide="none") +
    scale_x_discrete(labels=c(paste0("Deletion\n(n=",sum(d$status==0),")"),
      paste0("No deletion\n(n=",sum(d$status==1),")"))) +
    labs(title=paste0(gsub("_"," ",ca_study$study[i]),"\nDelta = ",sprintf("%+.2f",delta)),
      x=NULL,y="carA abundance (ShortBRED)") +
    theme_classic(base_size=13) +
    theme(plot.title=element_text(size=13,hjust=0.5),axis.text.x=element_text(angle=45,hjust=1))
}
ca_study$CI_low <- ca_study$beta-1.96*ca_study$SE
ca_study$CI_high <- ca_study$beta+1.96*ca_study$SE
ca_forest_data <- rbind(data.frame(study=gsub("_"," ",ca_study$study),beta=ca_study$beta,
  CI_low=ca_study$CI_low,CI_high=ca_study$CI_high),
  data.frame(study="Meta-analysis",beta=ca_meta$beta,CI_low=ca_meta$CI_low,CI_high=ca_meta$CI_high))
ca_forest_data$study <- factor(ca_forest_data$study,levels=rev(ca_forest_data$study))
ca_forest_data$label <- sprintf("%.2f (%.2f, %.2f)",ca_forest_data$beta,ca_forest_data$CI_low,ca_forest_data$CI_high)
ca_right <- max(ca_forest_data$CI_high)+0.15
ca_forest <- ggplot(ca_forest_data,aes(y=study,x=beta)) +
  geom_vline(xintercept=0,linetype="dashed",color="grey60") +
  geom_segment(aes(x=CI_low,xend=CI_high,yend=study),linewidth=0.6) +
  geom_point(aes(shape=study=="Meta-analysis"),color="#2166AC",size=2.5) +
  scale_shape_manual(values=c(16,18),guide="none") +
  geom_text(aes(x=ca_right,label=label),hjust=0,size=3.5) +
  scale_x_continuous(limits=c(min(0,min(ca_forest_data$CI_low)),ca_right+1.5),breaks=c(0,0.5,1,1.5,2)) +
  labs(x="Beta for rank-normalized abundance\n(no deletion vs deletion)",y=NULL) +
  theme_classic(base_size=13) +
  theme(axis.line.y=element_blank(),axis.ticks.y=element_blank())
ca_figure <- (wrap_plots(ca_plots,ncol=5) | ca_forest) +
  plot_layout(widths=c(5,2.5)) +
  plot_annotation(title=paste0("carA (CL2_16170): abundance by A. hadrus dSV status\n",
    "Meta beta = ",sprintf("%.2f",ca_meta$beta),"; 95% CI ",sprintf("%.2f to %.2f",ca_meta$CI_low,ca_meta$CI_high),
    "; Bonferroni P = ",formatC(ca_meta$Bonferroni_P,format="e",digits=2)),
    theme=theme(plot.title=element_text(size=19)))
ggsave("pics/FS8_carA_dSV_abundance_effects.tiff",ca_figure,width=26,height=16,units="in",dpi=300,compression="lzw")
ggsave("validation/carA_dSV_abundance_preview.png",ca_figure,width=26,height=16,dpi=75)
write.csv(ca_data,"12.AH_dSV_gene_abundance/carA_plot_samples.csv",row.names=FALSE)
