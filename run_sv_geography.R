### Geographic distributions of structural variants before quality filtering
### dSV: 0 = deletion, 1 = no deletion. Original input values are preserved.
### R 4.0.5

options(stringsAsFactors=FALSE, warn=1)
stopifnot(getRversion()=="4.0.5")
library(ggplot2)
library(patchwork)
setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

## 1 Parameters
out_dir <- "02.SV_summary_statistics/SV_geography_preQC_20261027_Europe_Bonferroni"
min_n <- 30L
sensitivity_n <- c(20L, 50L)
mc_B <- 10000L
seed <- 20261026L

## 2 Inputs and statistical helpers
stopifnot(min_n >= 5L, mc_B >= 10000L)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_tab <- function(x, name) write.csv(x, file.path(out_dir, paste0(name, ".csv")),
                                        row.names = FALSE, na = "", fileEncoding="UTF-8")
env <- new.env()
load("01.cleanData/SV_all/dsgv.RData", env)
load("01.cleanData/SV_all/vsgv.RData", env)
d <- as.matrix(env$dsgv_sub); v <- as.matrix(env$vsgv_sub)
stopifnot(is.numeric(d), is.numeric(v), all(is.na(d) | d %in% c(0, 1)))
stopifnot(all(is.na(v) | is.finite(v)))
meta <- read.delim("01.cleanData/phen/Clinical_basic_overlap.tsv",
                   check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(meta$match_id), !anyDuplicated(rownames(d)),
          !anyDuplicated(rownames(v)), !anyDuplicated(colnames(d)),
          !anyDuplicated(colnames(v)))
## Study-defined geographic grouping follows recruitment country.
## Israel and Turkiye are grouped with Europe as specified for this analysis.
## Preserve the source field and record every disagreement explicitly.
continents <- c("Africa", "Americas", "Asia", "Europe")
country_region <- c(CMR="Africa",MDG="Africa",USA="Americas",PER="Americas",SLV="Americas",CHN="Asia",JPN="Asia",KAZ="Asia",
  KOR="Asia",ISR="Europe",TUR="Europe",AUT="Europe",CZE="Europe",DEU="Europe",
  DNK="Europe",ESP="Europe",FRA="Europe",GBR="Europe",IRL="Europe",ITA="Europe",
  NLD="Europe",SWE="Europe")
meta$inter_source <- trimws(meta$inter)
meta$inter <- unname(country_region[trimws(meta$country)])
stopifnot(!anyNA(meta$inter))
write_tab(meta[meta$inter!=meta$inter_source,c("match_id","study_name","country","inter_source","inter")],
          "geography_classification_changes")
ids <- Reduce(intersect,list(rownames(d),rownames(v),meta$match_id))
meta <- meta[match(ids,meta$match_id),,drop=FALSE]
bad_encoding <- !validUTF8(meta$study_name)
meta$study_name[bad_encoding] <- iconv(meta$study_name[bad_encoding],from="GB18030",to="UTF-8")
subject<-trimws(as.character(meta$subject_id));missing<-is.na(subject)|tolower(subject)%in%c("","na","unknown","not available")
meta$subject_id_missing<-missing
meta$donor_key<-paste(meta$study_name,ifelse(missing,paste0("sample:",meta$match_id),subject),sep="::")
day<-suppressWarnings(as.numeric(meta$days_from_first_collection))
visit<-suppressWarnings(as.numeric(meta$visit_number))
## Prefer original visit/timepoint ordering when dates are unavailable.
proxy<-rep(NA_real_,nrow(meta))
li<-meta$study_name=="LiJ_2014"&grepl("_[0-9]+$",meta$sample_id)
proxy[li]<-as.numeric(sub(".*_([0-9]+)$","\\1",meta$sample_id[li]))
costea<-meta$study_name=="CosteaPI_2017"
kaz<-costea&grepl("_[0-9]+_[0-9]+_[0-9]+$",meta$sample_id)
proxy[kaz]<-as.numeric(sub(".*_[0-9]+_([0-9]+)_[0-9]+$","\\1",meta$sample_id[kaz]))
german<-costea&grepl("^alien_11_[0-9]+_0$",meta$sample_id)
proxy[german]<-as.numeric(sub("^alien_11_([0-9]+)_0$","\\1",meta$sample_id[german]))
sort_value<-ifelse(is.finite(day),day,ifelse(is.finite(visit),visit,proxy))
source<-ifelse(is.finite(day),"collection_day",ifelse(is.finite(visit),"visit_number",
        ifelse(is.finite(proxy),"original_sample_sequence_proxy","original_row_order_proxy")))
sort_value[!is.finite(sort_value)]<-seq_len(nrow(meta))[!is.finite(sort_value)]
## Group consistency is assessed explicitly; incomplete known chronology is flagged.
groups<-split(seq_len(nrow(meta)),meta$donor_key)
keep<-integer(length(groups));g<-0L;verified<-rep(TRUE,nrow(meta));choice<-rep("single_record",nrow(meta))
for(ii in groups) {
  g<-g+1L
  if(length(ii)==1L){keep[g]<-ii;next}
  if(all(is.finite(day[ii]))) {k<-ii[order(day[ii],meta$match_id[ii])][1];reason<-"earliest_collection_day"}
  else if(all(is.finite(visit[ii]))) {k<-ii[order(visit[ii],meta$match_id[ii])][1];reason<-"first_visit_number"}
  else if(all(is.finite(proxy[ii]))) {k<-ii[order(proxy[ii],meta$match_id[ii])][1];reason<-"original_sequence_proxy_time_unverified"}
  else {k<-ii[1];reason<-"original_row_order_proxy_time_unverified"}
  keep[g]<-k;choice[ii]<-reason;verified[ii]<-reason=="earliest_collection_day"
}
keep<-sort(keep);meta$retained<-seq_len(nrow(meta))%in%keep
meta$selection_rule<-choice;meta$chronological_earliest_verified<-verified
meta$ordering_proxy<-proxy
selection<-meta[,c("match_id","study_name","sample_id","subject_id","subject_id_missing",
  "donor_key","days_from_first_collection","visit_number","ordering_proxy","retained",
  "selection_rule","chronological_earliest_verified")]

write_tab(selection,"sample_selection_all")
meta <- meta[keep,,drop=FALSE];ids <- meta$match_id
write_tab(meta[,c("match_id","study_name","country","inter","subject_id","subject_id_missing",
                 "selection_rule","chronological_earliest_verified")],"retained_samples")
d <- d[ids,,drop=FALSE];v <- v[ids,,drop=FALSE]
stopifnot(!anyDuplicated(meta$donor_key))
meta$inter <- factor(meta$inter,levels=continents)
write_tab(data.frame(stage=c("input_dSV","input_vSV","aligned","one_per_identified_donor"),
  samples=c(nrow(env$dsgv_sub),nrow(env$vsgv_sub),nrow(selection),length(ids))),"donor_flow")
rm(env)

## Restore pre-INT vSV measurements with the pre-filter matrix's missingness.
an <- read.delim("01.cleanData/SV_info/vsgv_info_anno.tsv",check.names=FALSE,quote="")
an <- an[match(colnames(v),an$SV_Name),]
stopifnot(!anyNA(an$SV_ID),!anyDuplicated(an$SV_ID))
raw <- data.table::fread("00.rawData/SV/Whole_vsgv_20260921.csv",check.names=FALSE,data.table=FALSE)
ri <- match(ids,raw[[1]]);ci <- match(an$SV_ID,names(raw))
if(anyNA(ci)) {
  conversion <- read.csv("00.rawData/SV/vsgv_id_convert_20260921.csv")
  alternative <- conversion$id_final[match(an$SV_ID,conversion$id)]
  ci[is.na(ci)] <- match(alternative[is.na(ci)],names(raw))
}
stopifnot(!anyNA(ri),!anyNA(ci),!anyDuplicated(ci),!anyDuplicated(raw[[1]]))
pre_INT <- as.matrix(raw[ri,ci]);dimnames(pre_INT) <- dimnames(v)
stopifnot(all(is.finite(pre_INT[!is.na(v)])))
pre_INT[is.na(v)] <- NA_real_
stopifnot(identical(is.na(pre_INT),is.na(v)))
write_tab(data.frame(SV=colnames(v),source_id=an$SV_ID,source_column=names(raw)[ci],
  callable_n=colSums(!is.na(v))),"vsv_pre_INT_mapping")
v <- pre_INT;rm(raw,pre_INT)
input_files <- c("01.cleanData/SV_all/dsgv.RData","01.cleanData/SV_all/vsgv.RData",
 "01.cleanData/phen/Clinical_basic_overlap.tsv","01.cleanData/SV_info/vsgv_info_anno.tsv",
 "00.rawData/SV/Whole_vsgv_20260921.csv","00.rawData/SV/vsgv_id_convert_20260921.csv",
 "01.cleanData/mbio_all/SV_species_abun.tsv")
write_tab(data.frame(path=input_files,md5=unname(tools::md5sum(input_files))),"input_md5")
wilson <- function(k, n) {
  z <- qnorm(.975); p <- k / n; den <- 1 + z^2 / n
  mid <- (p + z^2 / (2*n)) / den
  half <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / den
  data.frame(lower = ifelse(n > 0, pmax(0, mid-half), NA_real_),
             upper = ifelse(n > 0, pmin(1, mid+half), NA_real_))
}
adjust_log <- function(x) {
  x$P_Bonferroni <- NA_real_; x$neglog10_P_Bonferroni <- NA_real_
  ix <- which(is.finite(x$log_p) | x$log_p == -Inf)
  if (length(ix)) {
    m <- length(ix)
    adjusted_log_p <- pmin(0, x$log_p[ix] + log(m))
    x$P_Bonferroni[ix] <- exp(adjusted_log_p)
    x$neglog10_P_Bonferroni[ix] <- -adjusted_log_p / log(10)
  }
  x
}
binary_test <- function(k, n) {
  tab <- rbind(deletion = k, no_deletion = n-k)
  fail <- list(statistic = NA_real_, df = NA_real_, p = NA_real_, log_p = NA_real_,
               method = "not_tested", sparse = NA, effect_size = NA_real_, mc_se = NA_real_)
  if (length(n) < 2L || any(rowSums(tab) == 0)) return(fail)
  expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  stat <- sum((tab-expected)^2/expected)
  sparse <- any(expected < 1) || mean(expected < 5) > .2
  df <- ncol(tab)-1L; mc_se <- NA_real_
  if (sparse && ncol(tab) == 2L) {
    p <- fisher.test(tab)$p.value; lp <- log(p); method <- "Fisher exact"
  } else if (sparse) {
    p <- chisq.test(tab, simulate.p.value = TRUE, B = mc_B)$p.value
    lp <- log(p); method <- "Conditional Monte Carlo chi-square"
    mc_se <- sqrt(p*(1-p)/(mc_B+1))
  } else {
    lp <- pchisq(stat, df = df, lower.tail = FALSE, log.p = TRUE)
    p <- exp(lp); method <- "Pearson chi-square"
  }
  list(statistic = stat, df = df, p = p, log_p = lp, method = method,
       sparse = sparse, effect_size = sqrt(stat/sum(tab)), mc_se = mc_se)
}
binary_summary <- function(x, group) {
  do.call(rbind, lapply(continents, function(g) {
    z <- x[which(group == g), , drop = FALSE]
    n <- colSums(!is.na(z)); k <- colSums(z == 0, na.rm = TRUE)
    ci <- wilson(k, n)
    data.frame(SV = colnames(x), continent = g, n = n, deletion_n = k,
               deletion_frequency = ifelse(n > 0, k/n, NA_real_),
               ci_low = ci$lower, ci_high = ci$upper, row.names = NULL)
  }))
}
continuous_summary <- function(x, group) {
  do.call(rbind, lapply(continents, function(g) {
    z <- x[which(group == g), , drop = FALSE]
    quant <- apply(z, 2, quantile, probs = c(.25, .5, .75), na.rm = TRUE, names = FALSE)
    data.frame(SV = colnames(x), continent = g, n = colSums(!is.na(z)),
               q25 = quant[1, ], median = quant[2, ], q75 = quant[3, ], row.names = NULL)
  }))
}
analyse_d <- function(s, threshold) {
  nc <- length(continents); ns <- nrow(s) / nc
  n <- matrix(s$n, ncol = nc); k <- matrix(s$deletion_n, ncol = nc)
  ans <- lapply(seq_len(ns), function(j) {
    ok <- n[j, ] >= threshold; test <- binary_test(k[j, ok], n[j, ok])
    fr <- k[j, ok]/n[j, ok]
    c(list(SV = s$SV[j], n = sum(n[j, ok]), n_continents = sum(ok),
           continents_tested = paste(continents[ok], collapse = ";"),
           status = if (sum(ok)<2) "fewer_than_two_eligible_continents" else
             if (is.na(test$p)) "no_variation" else "tested",
           max_frequency_difference = if (sum(ok)>=2) diff(range(fr)) else NA_real_), test)
  })
  x <- do.call(rbind, lapply(ans, as.data.frame, stringsAsFactors = FALSE))
  x <- adjust_log(x); x$species <- sub(":[0-9]+_[0-9]+.*$", "", x$SV)
  x
}
analyse_v <- function(x, group, s, threshold) {
  ns <- ncol(x); nmat <- matrix(s$n, ncol = length(continents))
  ans <- lapply(seq_len(ns), function(j) {
    ok <- nmat[j, ] >= threshold
    keep <- !is.na(x[, j]) & as.character(group) %in% continents[ok]
    y <- x[keep, j]; g <- droplevels(factor(group[keep]))
    testable <- sum(ok)>=2 && length(unique(y))>=2 && length(y)>sum(ok)
    h <- lp <- eff <- NA_real_
    if (testable) {
      fit <- kruskal.test(y, g); h <- unname(fit$statistic)
      lp <- pchisq(h, df = sum(ok)-1, lower.tail = FALSE, log.p = TRUE)
      eff <- max(0, (h-sum(ok)+1)/(length(y)-sum(ok)))
    }
    data.frame(SV = colnames(x)[j], n = length(y), n_continents = sum(ok),
               continents_tested = paste(continents[ok], collapse = ";"),
               status = if (testable) "tested" else "insufficient_groups_or_variation",
               statistic = h, df = if (testable) sum(ok)-1 else NA_real_,
               p = exp(lp), log_p = lp, effect_size = eff, method = "Kruskal-Wallis")
  })
  adjust_log(do.call(rbind, ans))
}



## 3 Overall geographic comparisons
set.seed(seed)
message("Geography: ", length(ids), " retained samples; one per identified donor; dSV values unchanged")
ds <- binary_summary(d, meta$inter); vs <- continuous_summary(v, meta$inter)
write_tab(ds, "dsv_continent_frequencies"); write_tab(vs, "vsv_continent_summaries")
dr <- analyse_d(ds, min_n); vr <- analyse_v(v, meta$inter, vs, min_n)
write_tab(dr, "dsv_overall"); write_tab(vr, "vsv_overall")
write_tab(dr[which(dr$P_Bonferroni < .05), ], "dsv_significant")
write_tab(vr[which(vr$P_Bonferroni < .05), ], "vsv_significant")
summarise_result <- function(x, type, subset, threshold) {
  data.frame(type = type, subset = subset, min_per_continent = threshold,
             input_SVs = nrow(x), tested = sum(!is.na(x$p)),
             significant = sum(x$P_Bonferroni < .05, na.rm = TRUE),
             min_P_Bonferroni = if (any(!is.na(x$P_Bonferroni))) min(x$P_Bonferroni, na.rm = TRUE) else NA_real_)
}
overview <- rbind(summarise_result(dr, "dSV", "all", min_n),
                  summarise_result(vr, "vSV", "all", min_n))



## 4 Sample-size and healthy-donor sensitivity analyses
for (th in sensitivity_n) {
  message("Minimum sample-size sensitivity: ", th)
  td <- analyse_d(ds, th); tv <- analyse_v(v, meta$inter, vs, th)
  write_tab(td, paste0("dsv_min_n", th)); write_tab(tv, paste0("vsv_min_n", th))
  overview <- rbind(overview, summarise_result(td,"dSV","all",th),
                    summarise_result(tv,"vSV","all",th))
}
healthy <- which(tolower(trimws(meta$disease)) == "healthy")
message("Healthy-only sensitivity: ", length(healthy), " donors")
hd <- analyse_d(binary_summary(d[healthy, , drop=FALSE],meta$inter[healthy]), min_n)
hv <- analyse_v(v[healthy, , drop=FALSE], meta$inter[healthy],
                continuous_summary(v[healthy, , drop=FALSE],meta$inter[healthy]), min_n)
write_tab(hd,"dsv_healthy"); write_tab(hv,"vsv_healthy")
overview <- rbind(overview, summarise_result(hd,"dSV","healthy",min_n),
                  summarise_result(hv,"vSV","healthy",min_n))
write_tab(overview, "analysis_summary")


## 5 Pairwise continent comparisons
# Test all eligible dSV continent pairs; Bonferroni covers the full pairwise family,
# including pairs whose overall test is not significant (no selected-only Bonferroni).
message("All eligible dSV pairwise comparisons")
pairs <- combn(continents, 2, simplify = FALSE)
pw <- do.call(rbind, lapply(pairs, function(pair) {
  a <- ds[ds$continent == pair[1], ]; b <- ds[ds$continent == pair[2], ]
  stopifnot(identical(a$SV,b$SV))
  jj <- which(a$n>=min_n & b$n>=min_n)
  if (!length(jj)) return(NULL)
  do.call(rbind,lapply(jj,function(j) {
    fit <- binary_test(c(a$deletion_n[j],b$deletion_n[j]),c(a$n[j],b$n[j]))
    data.frame(SV=a$SV[j],continent1=pair[1],continent2=pair[2],
               n1=a$n[j],n2=b$n[j],frequency1=a$deletion_frequency[j],
               frequency2=b$deletion_frequency[j],
               frequency_difference=a$deletion_frequency[j]-b$deletion_frequency[j],
               p=fit$p,log_p=fit$log_p,method=fit$method)
  }))
}))
if (!is.null(pw)) { pw <- adjust_log(pw); write_tab(pw,"dsv_pairwise_all") }



## 6 Representative SVs and cohort sensitivity
# Display selection: P_Bonferroni<0.05, >=3 eligible continents, >=10 percentage-point
# range; largest frequency range first, at most one SV per species, at most 20.
candidates <- dr[which(dr$P_Bonferroni < .05 & dr$n_continents>=3 & dr$max_frequency_difference>=.10), ]
candidates <- candidates[order(-candidates$max_frequency_difference,candidates$log_p), ]
top <- head(candidates[!duplicated(candidates$species), ],20)
top$display_id <- sprintf("D%02d",seq_len(nrow(top)))
write_tab(top,"representative_dsv_key")
study_tab <- as.data.frame(table(meta$study_name,meta$inter),stringsAsFactors=FALSE)
names(study_tab) <- c("study","continent","donors")
study_tab <- study_tab[study_tab$donors>0, ]
write_tab(study_tab,"study_continent_counts")
sam <- data.frame(continent=continents,donors=as.integer(table(meta$inter)),
                  studies=vapply(continents,function(g) length(unique(meta$study_name[meta$inter==g])),integer(1)))
write_tab(sam,"sampling_summary")
callability <- do.call(rbind,lapply(c("dSV","vSV"),function(type) {
  z <- if(type=="dSV") ds else vs
  do.call(rbind,lapply(continents,function(g) {
    n <- z$n[z$continent==g]
    data.frame(type=type,continent=g,any_call=sum(n>0),eligible=sum(n>=min_n),
               median_n=median(n),q25=unname(quantile(n,.25)),q75=unname(quantile(n,.75)))
  }))
}))
write_tab(callability,"callability_summary")
# Cohort-specific deletion frequencies and fixed high-vs-low continent contrasts.
# Selection is exploratory. LOCO does not reselect the compared continents.
cohort <- loco <- within <- list()
studies <- unique(meta$study_name)
for (j in seq_len(nrow(top))) {
  sv <- top$SV[j]; vals <- d[,match(sv,colnames(d))]
  s <- ds[ds$SV==sv & ds$n>=min_n, ]
  high <- s$continent[which.max(s$deletion_frequency)]
  low <- s$continent[which.min(s$deletion_frequency)]
  top$high_continent[j] <- high; top$low_continent[j] <- low
  for (h in seq_len(nrow(study_tab))) {
    ii <- which(meta$study_name==study_tab$study[h] & meta$inter==study_tab$continent[h])
    n <- sum(!is.na(vals[ii])); k <- sum(vals[ii]==0,na.rm=TRUE); ci<-wilson(k,n)
    cohort[[length(cohort)+1L]] <- data.frame(SV=sv,display_id=top$display_id[j],
        study=study_tab$study[h],continent=study_tab$continent[h],n=n,deletion_n=k,
        frequency=if(n) k/n else NA_real_,ci_low=ci$lower,ci_high=ci$upper)
  }
  for (st in studies) {
    get_counts <- function(g,retain) {
      z <- vals[which(meta$inter==g & retain)]
      c(k=sum(z==0,na.rm=TRUE),n=sum(!is.na(z)))
    }
    a<-get_counts(high,meta$study_name!=st); b<-get_counts(low,meta$study_name!=st)
    ok<-min(a["n"],b["n"])>=min_n
    fit<-if(ok) binary_test(c(a["k"],b["k"]),c(a["n"],b["n"])) else list(p=NA_real_,log_p=NA_real_)
    loco[[length(loco)+1L]]<-data.frame(SV=sv,display_id=top$display_id[j],excluded_study=st,
        high_continent=high,low_continent=low,n_high=a["n"],n_low=b["n"],eligible=ok,
        difference=if(ok) a["k"]/a["n"]-b["k"]/b["n"] else NA_real_,p=fit$p,log_p=fit$log_p,
        row.names=NULL)
    a<-get_counts(high,meta$study_name==st); b<-get_counts(low,meta$study_name==st)
    ok<-min(a["n"],b["n"])>=min_n
    if(ok) {
      fit<-binary_test(c(a["k"],b["k"]),c(a["n"],b["n"]))
      within[[length(within)+1L]]<-data.frame(SV=sv,study=st,high_continent=high,low_continent=low,
           n_high=a["n"],n_low=b["n"],difference=a["k"]/a["n"]-b["k"]/b["n"],
           p=fit$p,log_p=fit$log_p,row.names=NULL)
    }
  }
}
write_tab(top,"representative_dsv_key")
if(length(cohort)) { cohort<-do.call(rbind,cohort); write_tab(cohort,"representative_cohort_frequencies") }
if(length(loco)) { loco<-adjust_log(do.call(rbind,loco)); write_tab(loco,"representative_leave_one_study_out") }
if(length(within)) {
  write_tab(adjust_log(do.call(rbind,within)),"representative_within_study_contrasts")
} else {
  write_tab(data.frame(status="No representative had both its selected high- and low-frequency continents represented by >= min_n callable donors within the same study."),
            "representative_within_study_status")
}
saveRDS(list(meta=meta,overview=overview,dr=dr,vr=vr,ds=ds,vs=vs,top=top,
             cohort=cohort,loco=loco,sampling=sam,callability=callability),
        file.path(out_dir,"geography_results.rds"))


#################################################################################################
## 7 Representative SV covariate models
r <- readRDS(file.path(out_dir,"geography_results.rds"))
top <- r$top
if (nrow(top) > 0L) {
  e <- new.env(); load("01.cleanData/SV_all/dsgv.RData",e)
  d <- as.matrix(e$dsgv_sub[r$meta$match_id,,drop=FALSE]); m <- r$meta
  ab <- read.table("01.cleanData/mbio_all/SV_species_abun.tsv",header=TRUE,check.names=FALSE)
  ab <- ab[match(m$match_id,rownames(ab)),,drop=FALSE]
  results <- list()
  for(j in seq_len(nrow(top))) {
    sp <- top$species[j]
    if(!sp %in% colnames(ab)) stop("Abundance missing for species: ",sp)
    keep <- m$inter %in% c(top$high_continent[j],top$low_continent[j])
    health <- tolower(trimws(m$disease))
    x <- data.frame(dsv_state=d[,match(top$SV[j],colnames(d))],
       continent=factor(ifelse(m$inter==top$high_continent[j],"high","low"),levels=c("low","high")),
       age=suppressWarnings(as.numeric(m$age)),
       gender=factor(ifelse(tolower(m$gender) %in% c("male","female"),tolower(m$gender),NA)),
       health=factor(ifelse(is.na(health) | health=="",NA,ifelse(health=="healthy","healthy","nonhealthy"))),
       log_reads=log1p(suppressWarnings(as.numeric(m$read_count))),
       species_abundance=as.numeric(ab[,sp]))
    for(nm in c("age","log_reads","species_abundance")) x[[nm]][!is.finite(x[[nm]])] <- NA_real_
    xx <- droplevels(x[keep & complete.cases(x),,drop=FALSE])
    terms <- c("age","gender","health","log_reads","species_abundance")
    terms <- terms[vapply(xx[terms],function(z) length(unique(z))>=2,logical(1))]
    for(model in c("same_donors_unadjusted","covariate_adjusted")) {
      warnings <- character(); coef <- se <- p <- NA_real_; status <- "insufficient_data"
      form <- reformulate(c("continent",if(model=="covariate_adjusted") terms),response="dsv_state")
      if(nrow(xx)>=30 && length(unique(xx$continent))==2 &&
         length(unique(xx$dsv_state))==2 && min(table(xx$continent))>=30) {
        fit <- tryCatch(withCallingHandlers(glm(form,data=xx,family=binomial()),
           warning=function(w) { warnings<<-c(warnings,conditionMessage(w)); invokeRestart("muffleWarning") }),
           error=function(e) NULL)
        if(!is.null(fit)) {
          cf <- coef(summary(fit))
          if("continenthigh" %in% rownames(cf) && fit$converged &&
             is.finite(cf["continenthigh",2]) && cf["continenthigh",2]<10 &&
             !any(grepl("0 or 1|converge",warnings))) {
            coef <- cf["continenthigh",1]; se <- cf["continenthigh",2]; p <- cf["continenthigh",4]; status <- "estimated"
          } else status <- "unstable_or_separated_fit"
        } else status <- "model_error"
      }
      results[[length(results)+1L]] <- data.frame(SV=top$SV[j],display_id=top$display_id[j],
          high_continent=top$high_continent[j],low_continent=top$low_continent[j],model=model,
          N_before_complete_cases=sum(keep & !is.na(x$dsv_state)),N_complete=nrow(xx),
          n_high=sum(xx$continent=="high"),n_low=sum(xx$continent=="low"),
          outcome="1=no deletion; 0=deletion",coefficient=coef,OR_no_deletion=exp(coef),
          CI_low=exp(coef-1.96*se),CI_high=exp(coef+1.96*se),p=p,
          adjusted_terms=if(model=="covariate_adjusted") paste(terms,collapse=";") else "",
          status=status,warnings=paste(unique(warnings),collapse=";"))
    }
  }
  ans <- do.call(rbind,results)
  ans$exploratory_P_Bonferroni <- NA_real_
  for(model in unique(ans$model)) {
    ix <- which(ans$model==model & is.finite(ans$p)); ans$exploratory_P_Bonferroni[ix] <- p.adjust(ans$p[ix],"bonferroni")
  }
  write.csv(ans,file.path(out_dir,"representative_covariate_models.csv"),row.names=FALSE,na="")
}


#################################################################################################
capture.output(sessionInfo(),file=file.path(out_dir,"sessionInfo.txt"))
print(overview)
write_tab(data.frame(type=c("dSV","vSV"),
  eligible_tests=c(sum(!is.na(dr$p)),sum(!is.na(vr$p))),
  raw_p_cutoff=.05/c(sum(!is.na(dr$p)),sum(!is.na(vr$p)))),"primary_Bonferroni_families")
cat("GEOGRAPHY_ANALYSIS_COMPLETE",out_dir,"\n")

## 8 One supplementary figure. Main tests use pre-INT vSV values only.
## Two scatterplots show continental differences in dSVs and vSVs.
r <- readRDS(file.path(out_dir,"geography_results.rds"))
theme_set(theme_classic(base_size=12,base_family="Arial") +
  theme(axis.text=element_text(size=11,colour="black"),axis.title=element_text(size=12),
    plot.title=element_text(size=14,face="bold"),plot.subtitle=element_text(size=11),
    legend.text=element_text(size=11),legend.title=element_text(size=11),
    plot.margin=margin(10,12,10,10)))
sig_cols <- c("Bonf. P >= 0.05"="#B8B8B8","Bonf. P < 0.05"="#0072B2")
dr <- r$dr[!is.na(r$dr$P_Bonferroni),];vr <- r$vr[!is.na(r$vr$P_Bonferroni),]
dr$significance <- factor(ifelse(dr$P_Bonferroni<.05,"Bonf. P < 0.05","Bonf. P >= 0.05"),levels=names(sig_cols))
vr$significance <- factor(ifelse(vr$P_Bonferroni<.05,"Bonf. P < 0.05","Bonf. P >= 0.05"),levels=names(sig_cols))
dr <- dr[order(dr$P_Bonferroni,decreasing=TRUE),];vr <- vr[order(vr$P_Bonferroni,decreasing=TRUE),]
p_b <- ggplot(dr,aes(max_frequency_difference,neglog10_P_Bonferroni,colour=significance)) +
  geom_hline(yintercept=-log10(.05),linetype=2,colour="grey50",linewidth=.4) +
  geom_point(size=.65,alpha=.55) +
  scale_colour_manual(values=sig_cols,drop=FALSE,name=NULL) +
  scale_x_continuous(limits=c(0,1),breaks=c(0,.25,.5,.75,1)) +
  labs(x="Maximum deletion-frequency difference",y="-log10(Bonferroni P)") +
  theme(legend.position="bottom")
p_c <- ggplot(vr,aes(effect_size,neglog10_P_Bonferroni,colour=significance)) +
  geom_hline(yintercept=-log10(.05),linetype=2,colour="grey50",linewidth=.4) +
  geom_point(size=.65,alpha=.55) +
  scale_colour_manual(values=sig_cols,drop=FALSE,name=NULL) +
  labs(x="Rank effect size",y="-log10(Bonferroni P)") +
  theme(legend.position="bottom")
panel_ids <- c("dsv","vsv")
figure <- wrap_plots(p_b,p_c,ncol=2,nrow=1) +
  plot_layout(widths=c(1,1))
figure_dir <- file.path(out_dir,"supplementary_figure")
dir.create(figure_dir,showWarnings=FALSE)
width_in <- 11;height_in <- 4.5
## Measure the drawn plot areas; null gtable units are not physical rectangles.
grDevices::cairo_pdf(file.path(figure_dir,"geography_QA.pdf"),width=width_in,height=height_in,family="Arial")
grob <- patchwork::patchworkGrob(figure)
grid::grid.newpage();grid::grid.draw(grob);grid::grid.force()
panels <- grob$layout[grepl("^panel(?:-[0-9]+)?$",grob$layout$name,perl=TRUE),]
panels <- panels[order(panels$t,panels$l),]
rectangles <- matrix(NA_real_,nrow(panels),4);layout <- vector("list",nrow(panels))
for(k in seq_len(nrow(panels))) {
  q <- panels[k,]
  grid::seekViewport(paste0(q$name,".",q$t,"-",q$l,"-",q$b,"-",q$r))
  xy <- grid::deviceLoc(grid::unit(c(0,1),"npc"),grid::unit(c(0,1),"npc"),valueOnly=TRUE)
  rectangles[k,] <- c(xy$x[1],xy$y[1],xy$x[2],xy$y[2])*72
  layout[[k]] <- list(id=panel_ids[k],bbox_pt=rectangles[k,],grid_id="geography-grid",
    row_start=q$t-1L,row_stop=q$b,col_start=q$l-1L,col_stop=q$r)
  grid::upViewport(0)
}
dev.off()
stopifnot(nrow(rectangles)==2L,all(rectangles[,3]>rectangles[,1]),all(rectangles[,4]>rectangles[,2]),
 diff(range(rectangles[,3]-rectangles[,1]))<=1.5,diff(range(rectangles[,4]-rectangles[,2]))<=1.5,
 diff(range(rectangles[,2]))<=1.5,diff(range(rectangles[,4]))<=1.5)
jsonlite::write_json(list(schema_version=1L,backend="r-patchwork-drawn-viewports",
 figure=list(width_pt=width_in*72,height_pt=height_in*72),panels=layout,
 exemptions=list()),
 file.path(figure_dir,"panel_layout.json"),auto_unbox=TRUE,pretty=TRUE,digits=10)
ragg::agg_tiff(file.path(figure_dir,"Supplementary_SV_geographic_distributions.tiff"),
 width=width_in,height=height_in,units="in",res=300,compression="lzw")
grid::grid.newpage();grid::grid.draw(grob);dev.off()
ragg::agg_png(file.path(figure_dir,"geography_preview.png"),
 width=width_in,height=height_in,units="in",res=140)
grid::grid.newpage();grid::grid.draw(grob);dev.off()
file.copy(file.path(figure_dir,"Supplementary_SV_geographic_distributions.tiff"),
 "pics/Supplementary_SV_geographic_distributions_Europe_Bonferroni_axis_fixed.tiff",overwrite=TRUE)
file.copy(file.path(figure_dir,"Supplementary_SV_geographic_distributions.tiff"),
 "pics/Supplementary_SV_geographic_distributions_preQC_Europe_Bonferroni.tiff",overwrite=TRUE)
write.csv(dr,file.path(figure_dir,"panel_D_source.csv"),row.names=FALSE)
write.csv(vr,file.path(figure_dir,"panel_E_source.csv"),row.names=FALSE)
cat("GEOGRAPHY_FIGURE_COMPLETE\n")

## 9 Reporting notes and reproducibility
writeLines(c(
 "Input SVs and samples are from the combined dataset before quality filtering (SV_all, phen and SV_info).",
 "Primary vSV analysis uses pre-INT values only, preserving the pre-filter matrix's missingness; no post-QC mask is applied.",
 "Original dSV values are unchanged. Deletion frequency = count(x==0)/count(!is.na(x)).",
 paste0("One record per identified donor within each study; ",sum(r$meta$subject_id_missing)," records without donor IDs retained."),
 paste0(sum(!r$meta$chronological_earliest_verified)," repeated-donor selections use original sample sequence with chronology unverified."),
 "Study-defined recruitment-country grouping assigns ISR and TUR to Europe, as specified by the author. Input metadata are unchanged.",
 "This grouping is not described as the UN M49 continental classification.",
 "Each SV requires >=30 callable samples in each of >=2 eligible continents. Other continents are excluded for that SV.",
 "dSV: Pearson chi-square; sparse tables use Fisher exact (2 groups) or Monte Carlo chi-square (B=10000).",
 paste0(sum(r$dr$method=="Pearson chi-square")," primary dSV tests in this data use the Pearson chi-square branch."),
 "vSV: Kruskal-Wallis. Rank effect = max(0,(H-k+1)/(N-k)). Pre-INT input scores are not absolute copy numbers.",
 "Bonferroni correction is separate for all eligible dSV tests and all eligible vSV tests in each analysis.",
 "All eligible dSV pairwise continent comparisons form a separate full Bonferroni family, without screening by the overall P value.",
 "Threshold sensitivity uses 20 and 50 callable samples per continent. Healthy-only results are secondary.",
 "Representatives are selected by P_Bonferroni<0.05, >=3 eligible continents, >=0.10 frequency range, descending range, one locus per species.",
 "20 representative dSVs are checked in sensitivity analyses. The figure shows all eligible dSVs and vSVs in two scatterplots.",
 "Leave-one-study-out contrasts fix the original high and low continents; they do not reselect them.",
 "Representative logistic models preserve outcome 1=no deletion; an OR<1 means more deletion in the high-deletion continent.",
 "Adjusted models use age, sex, healthy/nonhealthy status, log1p(read count), and source-species relative abundance.",
 "Only estimable models can support adjustment claims. Selected-model Bonferroni results are exploratory.",
 "Pooled continental comparisons do not separate geography from study effects. Study representation differs across continents.",
 "Small P values alone do not establish large or biologically important differences; effect sizes are shown.",
 "TIFF is the submission figure at 300 dpi. PNG and PDF are internal viewing/QA files.",
 "Rendered panel geometry is measured in R and checked to 1.5 pt. No automatic PDF collision-audit pass is claimed when PyMuPDF is unavailable."
),file.path(out_dir,"analysis_methods_and_limits.txt"))
stopifnot(identical(unname(tools::md5sum(input_files)),read.csv(file.path(out_dir,"input_md5.csv"))$md5))
file.copy("run_sv_geography.R",file.path(out_dir,"run_sv_geography.source.R"),overwrite=TRUE)
