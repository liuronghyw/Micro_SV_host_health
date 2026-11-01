### Supplementary Tables 3-7: retain all Table 3 results; Normal P < 0.05 for Tables 4, 5, 7; Normal P < 1e-7 for Table 6
### Formatting follows supplementary_files/Supplementary Table 3.xlsx and Supplementary Table 4.xlsx
### LiuRong
### 2026-10-23


setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")
Sys.setlocale("LC_CTYPE","Chinese_China.936")
library(data.table)
library(openxlsx)
output_dir<-"F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R/Supplementary_tables"
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)

phenotype_labels<-c(age="Age",BMI="BMI",IBD="IBD",T2D="T2D",CRC="CRC",MetS="MetS",
  CAD="CAD",ACVD="ACVD",adenoma="Adenoma",cirrhosis="Cirrhosis",PD="PD",hba1c="HbA1c",
  fasting_glucose="Fasting glucose",fasting_insulin="Fasting insulin",triglycerides="Triglycerides",
  hdl="HDL",ldl="LDL",cholesterol="Cholesterol",mean_press="Mean pressure",
  creatinine="Creatinine",alt="ALT",ast="AST",total_bilirubin="Total bilirubin")
binary_phenotypes<-c("CRC","IBD","T2D","MetS","CAD","ACVD","cirrhosis","adenoma","PD")

header_style<-createStyle(textDecoration="bold",fgFill="#D9D9D9",halign="center",
                          valign="center",wrapText=TRUE,border="bottom",borderColour="#A6A6A6")
body_style<-createStyle(valign="center")
identifier_style<-createStyle(wrapText=TRUE,valign="center")
number_style<-createStyle(numFmt="0.0000E+00",halign="right")
pvalue_style<-createStyle(numFmt="0.00E+00",halign="right")
integer_style<-createStyle(numFmt="0",halign="right")
I2_style<-createStyle(numFmt="0.00",halign="right")
caption_style<-createStyle(fontSize=10,fontColour="#595959",valign="center")
rename_columns<-c(X="SV_name_short",SV_Name_long="SV_name",type="SV_type",
  effect_meta="meta.beta",OR_meta="meta.OR",pvalue_meta="meta.p",se_meta="meta.SE",
  cilb="meta.OR_left",ciub="meta.OR_right",I2="meta.I2 (%)",het_p="meta.het_p",
  p_count="meta.p_count",fdr="meta.fdr",bonferroni="meta.bonferroni",
  product_prokka="Products",product_kegg="Product_kegg")
normalize_calvet <- function(x) {
  gsub("[^\\t\\r\\n;,|()]*_Calvet_2024", "Fern\u00e1ndez_Calvet_2024", as.character(x), perl=TRUE)
}
normalize_calvet_table <- function(x) {
  names(x) <- normalize_calvet(names(x))
  for(j in seq_along(x)) if(is.character(x[[j]]) || is.factor(x[[j]]))
    x[[j]] <- normalize_calvet(x[[j]])
  x
}

## Supplementary Table 3: species-level SV composition and relative abundance
## Follow the reference table's sheet order; include every phenotype and every species.
table3_files<-list.files("06.Species_genetic_association/tableS3",pattern="\\.csv$",full.names=TRUE)
table3_phenotypes<-sub("\\.csv$","",basename(table3_files))
table3_labels<-c(age_norm="Age",BMI_norm="BMI",IBD="IBD",CRC="CRC",adenoma="Adenoma",
  T2D="T2D",MetS="MetS",hba1c_norm="HbA1c",fasting_glucose_norm="Fasting glucose",
  fasting_insulin_norm="Fasting insulin",PD="PD",CAD="CAD",ACVD="ACVD",
  mean_press_norm="Mean pressure",cholesterol_norm="Cholesterol",triglycerides_norm="Triglycerides",
  hdl_norm="HDL",ldl_norm="LDL",cirrhosis="Cirrhosis",alt_norm="ALT",ast_norm="AST",
  total_bilirubin_norm="Total bilirubin",creatinine_norm="Creatinine")
table3_order<-unique(c(names(table3_labels),table3_phenotypes))
table3_phenotypes<-table3_order[table3_order %in% table3_phenotypes]
table3_files<-table3_files[match(table3_phenotypes,sub("\\.csv$","",basename(table3_files)))]
stopifnot(length(table3_files)>0,!anyDuplicated(table3_phenotypes),length(table3_files)<=length(LETTERS))
wb3<-createWorkbook(creator="LiuRong",title="Supplementary Table 3",
                    subject="Species-level associations of SV composition and relative abundance")
modifyBaseFont(wb3,fontName="Calibri",fontSize=11,fontColour="#000000")

for(i in seq_along(table3_files)){
  phenotype3<-table3_phenotypes[i]
  label3<-if(phenotype3 %in% names(table3_labels)) table3_labels[[phenotype3]] else phenotype3
  result3<-read.csv(table3_files[i],fileEncoding="GB18030",check.names=FALSE,stringsAsFactors=FALSE)
  result3<-normalize_calvet_table(result3)
  stopifnot(all(c("Short_name","sv.MetaSigAssoc","abun.MetaSigAssoc") %in% names(result3)),
            !anyNA(result3$Short_name),!anyDuplicated(result3$Short_name))
  names(result3)[names(result3)=="Short_name"]<-"Species"
  flag_columns3<-match(c("sv.MetaSigAssoc","abun.MetaSigAssoc"),names(result3))
  stopifnot(all(unlist(result3[,flag_columns3]) %in% c("Yes","No","Unknown")))
  numeric_columns3<-setdiff(seq_len(ncol(result3)),c(1L,flag_columns3))
  for(j in numeric_columns3) result3[[j]]<-as.numeric(result3[[j]])
  infinite_rows3<-which(Reduce(`|`,lapply(result3[,numeric_columns3,drop=FALSE],is.infinite)))
  if(length(infinite_rows3)>0L) result3[infinite_rows3,numeric_columns3]<-NA_real_
  stopifnot(!any(vapply(result3[,numeric_columns3,drop=FALSE],function(x) any(is.infinite(x)),logical(1))))

  sheet3<-paste0("ST3",LETTERS[i],"_",label3)
  method3<-if(sub("_norm$","",phenotype3) %in% binary_phenotypes) "logistic regression" else "linear regression"
  addWorksheet(wb3,sheet3,gridLines=TRUE,zoom=90)
  writeData(wb3,sheet3,paste0("Table S3",LETTERS[i],". Species-level associations of ",label3,
    " with SV composition (PERMANOVA) and relative abundance (",method3,")."),startRow=1,colNames=FALSE)
  writeData(wb3,sheet3,paste0("All ",nrow(result3)," species retained. Significance: FDR < 0.05. ",
    "Unknown denotes unavailable FDR; NA denotes an unavailable numeric result."),startRow=2,colNames=FALSE)
  writeData(wb3,sheet3,result3,startRow=3,headerStyle=header_style,keepNA=TRUE,na.string="NA",withFilter=TRUE)
  addStyle(wb3,sheet3,caption_style,rows=2,cols=1)
  setColWidths(wb3,sheet3,cols=seq_len(ncol(result3)),widths=26)
  setColWidths(wb3,sheet3,cols=1,widths=28)
  setRowHeights(wb3,sheet3,rows=1:3,heights=c(21,19,55))
  freezePane(wb3,sheet3,firstActiveRow=4,firstActiveCol=2)
  rows3<-seq_len(nrow(result3))+3L
  addStyle(wb3,sheet3,body_style,rows=rows3,cols=seq_len(ncol(result3)),gridExpand=TRUE,stack=TRUE)
  addStyle(wb3,sheet3,number_style,rows=rows3,cols=numeric_columns3,gridExpand=TRUE,stack=TRUE)
  p_columns3<-grep("(^|\\.)(p|pvalue_meta|fdr)(\\.p)?$",names(result3),ignore.case=TRUE)
  addStyle(wb3,sheet3,pvalue_style,rows=rows3,cols=p_columns3,gridExpand=TRUE,stack=TRUE)
  count_columns3<-grep("\\.N$",names(result3))
  addStyle(wb3,sheet3,integer_style,rows=rows3,cols=count_columns3,gridExpand=TRUE,stack=TRUE)
  setRowHeights(wb3,sheet3,rows=rows3,heights=17)
  message(label3,": ",nrow(result3)," species; ",ncol(result3)," columns")
}

for(i in seq_along(wb3$worksheets_rels)){
  wb3$worksheets_rels[[i]]<-grep("/drawing|/vmlDrawing",wb3$worksheets_rels[[i]],invert=TRUE,value=TRUE)
}
wb3$Content_Types<-grep("/xl/drawings/",wb3$Content_Types,invert=TRUE,value=TRUE)
saveWorkbook(wb3,file.path(output_dir,"Supplementary Table 3.xlsx"),overwrite=TRUE)
rm(wb3,result3)
invisible(gc())

## Supplementary Tables 4 and 5: species abundance; species abundance and SV PCs
sv_input_dirs<-c("07.Microbial_GWAS/anno_adjAbun","07.Microbial_GWAS/anno_adjAbunPCs")
for(table_number in 4:5){
  input_dir<-sv_input_dirs[table_number-3L]
  adjustment<-if(table_number==4L) "species abundance" else "species abundance and SV PCs"
  files<-list.files(input_dir,pattern="_sv\\.csv$",full.names=TRUE)
  phenotypes<-sub("_sv\\.csv$","",basename(files))
  phenotype_order<-unique(c(names(phenotype_labels),phenotypes))
  phenotypes<-phenotype_order[phenotype_order %in% phenotypes]
  files<-files[match(phenotypes,sub("_sv\\.csv$","",basename(files)))]
  stopifnot(length(files)>0,!anyDuplicated(phenotypes),length(files)<=length(LETTERS))
  wb<-createWorkbook(creator="LiuRong",title=paste0("Supplementary Table ",table_number),
                     subject=paste0("Microbial SV associations with Normal P < 0.05 adjusted for ",adjustment))
  modifyBaseFont(wb,fontName="Calibri",fontSize=11,fontColour="#000000")
  for(i in seq_along(files)){
    phenotype<-phenotypes[i]
    label<-if(phenotype %in% names(phenotype_labels)) phenotype_labels[[phenotype]] else phenotype
    result<-as.data.frame(fread(files[i],na.strings=c("NA","")))
    result<-normalize_calvet_table(result)
    ## Explicit numeric types also cover extreme scientific notation read as text by fread.
    numeric_columns<-grep("\\.(beta|OR|SE|p|OR_left|OR_right|N)$|^(effect_meta|OR_meta|pvalue_meta|se_meta|cilb|ciub|I2|het_p|p_count|bonferroni|fdr|Segments_num)$",names(result))
    for(j in numeric_columns){
      result[[j]][result[[j]] %in% c("NA","")]<-NA
      result[[j]]<-as.numeric(result[[j]])
    }
    numeric_columns<-union(numeric_columns,which(vapply(result,is.numeric,logical(1))))
    infinite_rows<-which(Reduce(`|`,lapply(result[,numeric_columns,drop=FALSE],is.infinite)))
    if(length(infinite_rows)>0L) result[infinite_rows,numeric_columns]<-NA_real_
    required<-c("X","type","SV_Name_long","SV_ID","Segments_num","Genes","product_prokka",
                "product_kegg","Orthology","pvalue_meta","I2","bonferroni","fdr")
    stopifnot(all(required %in% names(result)),!anyNA(result$SV_Name_long),!anyNA(result$SV_ID),
              !anyDuplicated(result[,c("SV_ID","type")]),all(result$type %in% c("vsv","dsv")))
    cohort_columns<-grep("\\.(beta|OR|SE|p|OR_left|OR_right|N)$",names(result),value=TRUE)
    cohort_effect<-grep(if(phenotype %in% binary_phenotypes) "\\.OR$" else "\\.beta$",
                        cohort_columns,value=TRUE)
    stopifnot(length(cohort_effect)>0)

    ## Single-cohort files store their estimate and confidence interval in cohort columns.
    if(!any(c("effect_meta","OR_meta") %in% names(result))){
      stopifnot(length(cohort_effect)==1L)
      effect_name<-if(phenotype %in% binary_phenotypes) "OR_meta" else "effect_meta"
      result[[effect_name]]<-result[[cohort_effect]]
      if(phenotype %in% binary_phenotypes){
        result$cilb<-result[[sub("\\.OR$",".OR_left",cohort_effect)]]
        result$ciub<-result[[sub("\\.OR$",".OR_right",cohort_effect)]]
      }
    }

    ## SV coordinates are 1-kbp bins; size is the sum of end minus start across segments.
    regions<-strsplit(sub("^.*:","",result$SV_Name_long),";",fixed=TRUE)
    result$SV_size<-vapply(regions,function(x){
      bounds<-do.call(rbind,strsplit(x,"_",fixed=TRUE))
      stopifnot(ncol(bounds)==2L)
      sum(as.numeric(bounds[,2])-as.numeric(bounds[,1]))
    },numeric(1))
    result$SV_size[is.na(result$Segments_num)]<-NA_real_
    stopifnot(all(is.na(result$SV_size) | is.finite(result$SV_size) & result$SV_size>0),
              all(is.na(result$Segments_num) | lengths(regions)==result$Segments_num))
    result$FDR_less_or_equal_0.05<-ifelse(is.na(result$fdr),NA_character_,ifelse(result$fdr<=0.05,"Yes","No"))
    result$Is_Bonferroni_p_less_than_0.05<-ifelse(is.na(result$bonferroni),NA_character_,
                                                ifelse(result$bonferroni<0.05,"Yes","No"))
    result<-result[!is.na(result$pvalue_meta) & result$pvalue_meta<0.05,,drop=FALSE]
    result<-result[order(match(result$type,c("vsv","dsv")),result$pvalue_meta,result$SV_Name_long),,drop=FALSE]
    result$type<-c(vsv="vSV",dsv="dSV")[result$type]
    rename_index<-match(names(result),names(rename_columns))
    names(result)[!is.na(rename_index)]<-unname(rename_columns[rename_index[!is.na(rename_index)]])
    names(result)[names(result)=="SV_size"]<-"SV_size (kbp)"
    first_columns<-c("SV_name","SV_name_short","SV_type","meta.beta","meta.OR","meta.p",
      "meta.I2 (%)","meta.fdr","meta.bonferroni","Segments_num","SV_size (kbp)","Genes",
      "Products","Orthology","FDR_less_or_equal_0.05","Is_Bonferroni_p_less_than_0.05")
    first_columns<-intersect(first_columns,names(result))
    result<-result[,c(first_columns,setdiff(names(result),c(first_columns,cohort_columns)),cohort_columns),drop=FALSE]

    sheet<-paste0("ST",table_number,LETTERS[i],"_",label)
    addWorksheet(wb,sheet,gridLines=TRUE,zoom=90)
    writeData(wb,sheet,paste0("Table S",table_number,LETTERS[i],". Associations between microbial SVs and ",label,
                             " adjusted for ",adjustment,"."),startRow=1,colNames=FALSE)
    caption<-paste0(nrow(result)," associations with Normal P < 0.05. Flags: FDR <= 0.05; Bonferroni P < 0.05. NA denotes an unavailable result.")
    if(nrow(result)>0L && (length(cohort_effect)==1L || "meta.p_count" %in% names(result) && any(result$meta.p_count==1,na.rm=TRUE)))
      caption<-paste0(caption," Single-cohort I2 is not estimated (source value: 1).")
    stopifnot(!any(vapply(result,function(x) is.numeric(x) && any(is.infinite(x)),logical(1))))
    writeData(wb,sheet,caption,startRow=2,colNames=FALSE)
    writeData(wb,sheet,result,startRow=3,headerStyle=header_style,keepNA=TRUE,na.string="NA",withFilter=nrow(result)>0)
    addStyle(wb,sheet,caption_style,rows=2,cols=1)
    setColWidths(wb,sheet,cols=seq_len(ncol(result)),widths=16)
    setColWidths(wb,sheet,cols=match(c("SV_name","SV_name_short"),names(result)),widths=c(62,48))
    setColWidths(wb,sheet,cols=match("SV_type",names(result)),widths=10)
    setColWidths(wb,sheet,cols=match(c("Genes","Products","Product_kegg","Orthology"),names(result)),widths=38)
    setColWidths(wb,sheet,cols=match(c("FDR_less_or_equal_0.05","Is_Bonferroni_p_less_than_0.05"),names(result)),widths=c(23,29))
    setColWidths(wb,sheet,cols=match(cohort_columns,names(result)),widths=22)
    setRowHeights(wb,sheet,rows=1:3,heights=c(21,19,46))
    freezePane(wb,sheet,firstActiveRow=4,firstActiveCol=2)

    if(nrow(result)>0L){
      rows<-seq_len(nrow(result))+3L
      addStyle(wb,sheet,body_style,rows=rows,cols=seq_len(ncol(result)),gridExpand=TRUE,stack=TRUE)
      addStyle(wb,sheet,identifier_style,rows=rows,cols=match(c("SV_name","SV_name_short"),names(result)),gridExpand=TRUE,stack=TRUE)
      numeric_columns<-which(vapply(result,is.numeric,logical(1)))
      addStyle(wb,sheet,number_style,rows=rows,cols=numeric_columns,gridExpand=TRUE,stack=TRUE)
      p_columns<-grep("^meta\\.(p|fdr|bonferroni|het_p)$|\\.p$",names(result))
      addStyle(wb,sheet,pvalue_style,rows=rows,cols=p_columns,gridExpand=TRUE,stack=TRUE)
      count_columns<-grep("^Segments_num$|^SV_size|^meta\\.p_count$|\\.N$",names(result))
      addStyle(wb,sheet,integer_style,rows=rows,cols=count_columns,gridExpand=TRUE,stack=TRUE)
      addStyle(wb,sheet,I2_style,rows=rows,cols=match("meta.I2 (%)",names(result)),gridExpand=TRUE,stack=TRUE)
      row_heights<-17*pmax(1,ceiling(nchar(result$SV_name)/60),ceiling(nchar(result$SV_name_short)/46))
      setRowHeights(wb,sheet,rows=rows,heights=pmin(409,row_heights))
    }
    message("Table ",table_number," ",label,": ",nrow(result)," associations (",sum(result$SV_type=="vSV")," vSV; ",sum(result$SV_type=="dSV")," dSV)")
  }

  # Remove unused drawing references emitted by openxlsx 4.2.5.2 for table-only sheets.
  for(i in seq_along(wb$worksheets_rels)){
    wb$worksheets_rels[[i]]<-grep("/drawing|/vmlDrawing",wb$worksheets_rels[[i]],invert=TRUE,value=TRUE)
  }
  wb$Content_Types<-grep("/xl/drawings/",wb$Content_Types,invert=TRUE,value=TRUE)
  saveWorkbook(wb,file.path(output_dir,paste0("Supplementary Table ",table_number,".xlsx")),overwrite=TRUE)
  rm(wb,result)
  invisible(gc())
}

if (!isTRUE(getOption("sv.skip_humann3", FALSE))) {
## Supplementary Table 6: HUMAnN3 UniRef90 abundance associations with Normal P < 1e-7
table6_phenotypes<-names(phenotype_labels)
table6_files<-file.path("10.Microbial_humann3/meta",paste0(table6_phenotypes,".csv"))
if(any(!file.exists(table6_files))) stop("Table 6 missing phenotype inputs: ",paste(table6_phenotypes[!file.exists(table6_files)],collapse=", "))
table6_p_cutoff<-1e-7
table6_annotation<-fread("00.rawData/humann3/match_file/SGV_Finder_genes_vs_uniref90_annotation.tsv",
                         select=c("UniRef90_id","gene","ko","name","SV_location","dest"),na.strings=c("NA",""),encoding="UTF-8")
stopifnot(!anyNA(table6_annotation$UniRef90_id))
table6_annotation_audit<-data.frame()

table6_rename<-c(gene="Genes",dest="Dest",effect_meta="meta.beta",OR_meta="meta.OR",
  pvalue_meta="meta.p",se_meta="meta.SE",I2="meta.I2 (%)",het_p="meta.het_p",p_count="meta.p_count",fdr="meta.fdr",bonferroni="meta.bonferroni")
table6_numeric_pattern<-"\\.(beta|OR|SE|p|beta_left|beta_right|OR_left|OR_right|N)$|^(effect_meta|OR_meta|pvalue_meta|se_meta|cilb|ciub|I2|het_p|p_count|study_count|bonferroni|fdr)$"
table6_headers<-vector("list",length(table6_files))
wb6<-createWorkbook(creator="LiuRong",title="Supplementary Table 6",subject="HUMAnN3 UniRef90 abundance associations with Normal P < 1e-7")
modifyBaseFont(wb6,fontName="Calibri",fontSize=11,fontColour="#000000")

## Prepare a small formatting template; write the full data to worksheet XML in batches below.
for(i in seq_along(table6_files)){
  phenotype6<-table6_phenotypes[i]
  label6<-if(phenotype6 %in% names(phenotype_labels)) phenotype_labels[[phenotype6]] else phenotype6
  result6<-as.data.frame(fread(table6_files[i],nrows=1,na.strings=c("NA",""),encoding="UTF-8"))
  result6<-normalize_calvet_table(result6)
  if("X" %in% names(result6)) names(result6)[names(result6)=="X"]<-"UniRef90_id"
  stopifnot(all(c("UniRef90_id","pvalue_meta","I2") %in% names(result6)))
  for(column6 in c("gene","ko","name","SV_location","dest")) result6[[column6]]<-NA_character_
  result6$bonferroni<-result6$fdr<-NA_real_
  numeric_columns6<-grep(table6_numeric_pattern,names(result6))
  for(j in numeric_columns6) result6[[j]]<-as.numeric(result6[[j]])
  result6$FDR_less_or_equal_0.05<-ifelse(is.na(result6$fdr),NA_character_,ifelse(result6$fdr<=0.05,"Yes","No"))
  result6$Is_Bonferroni_p_less_than_0.05<-ifelse(is.na(result6$bonferroni),NA_character_,ifelse(result6$bonferroni<0.05,"Yes","No"))
  rename6<-c(table6_rename,cilb=if(phenotype6 %in% binary_phenotypes) "meta.OR_left" else "meta.beta_left",
    ciub=if(phenotype6 %in% binary_phenotypes) "meta.OR_right" else "meta.beta_right")
  rename_index6<-match(names(result6),names(rename6))
  names(result6)[!is.na(rename_index6)]<-unname(rename6[rename_index6[!is.na(rename_index6)]])
  cohort_columns6<-grep("\\.(beta|OR|SE|p|beta_left|beta_right|OR_left|OR_right|N)$",setdiff(names(result6),unname(rename6)),value=TRUE)
  first_columns6<-intersect(c("UniRef90_id","Genes","ko","name","meta.beta","meta.OR","meta.p","meta.I2 (%)","meta.fdr","meta.bonferroni",
    "SV_location","Dest","FDR_less_or_equal_0.05","Is_Bonferroni_p_less_than_0.05"),names(result6))
  result6<-result6[,c(first_columns6,setdiff(names(result6),c(first_columns6,cohort_columns6)),cohort_columns6),drop=FALSE]
  stopifnot(identical(names(result6)[1:4],c("UniRef90_id","Genes","ko","name")))
  table6_headers[[i]]<-names(result6)
  sheet6<-paste0("ST6",LETTERS[i],"_",label6)
  addWorksheet(wb6,sheet6,gridLines=TRUE,zoom=90)
  writeData(wb6,sheet6,paste0("Table S6",LETTERS[i],". Associations between HUMAnN3 UniRef90 abundance and ",label6,"."),startRow=1,colNames=FALSE)
  writeData(wb6,sheet6,"UniRef90 results with Normal P < 1e-7. ko and name matched by UniRef90_id; missing annotations shown as NA. Flags: FDR <= 0.05; Bonferroni P < 0.05. Only raw meta P determines inclusion. Nonfinite estimates are shown as NA.",startRow=2,colNames=FALSE)
  writeData(wb6,sheet6,result6,startRow=3,headerStyle=header_style,keepNA=TRUE,na.string="NA",withFilter=TRUE)
  addStyle(wb6,sheet6,caption_style,rows=2,cols=1)
  setColWidths(wb6,sheet6,cols=seq_len(ncol(result6)),widths=22)
  setColWidths(wb6,sheet6,cols=match(c("UniRef90_id","Genes","ko","name","Dest"),names(result6)),widths=c(26,48,28,48,32))
  setColWidths(wb6,sheet6,cols=match(c("FDR_less_or_equal_0.05","Is_Bonferroni_p_less_than_0.05"),names(result6)),widths=c(23,29))
  setRowHeights(wb6,sheet6,rows=1:4,heights=c(21,19,55,17))
  freezePane(wb6,sheet6,firstActiveRow=4,firstActiveCol=3)
  addStyle(wb6,sheet6,body_style,rows=4,cols=seq_len(ncol(result6)),gridExpand=TRUE,stack=TRUE)
  addStyle(wb6,sheet6,identifier_style,rows=4,cols=match(c("UniRef90_id","name"),names(result6)),gridExpand=TRUE,stack=TRUE)
  addStyle(wb6,sheet6,number_style,rows=4,cols=which(vapply(result6,is.numeric,logical(1))),gridExpand=TRUE,stack=TRUE)
  addStyle(wb6,sheet6,pvalue_style,rows=4,cols=grep("^meta\\.(p|fdr|bonferroni|het_p)$|\\.p$",names(result6)),gridExpand=TRUE,stack=TRUE)
  addStyle(wb6,sheet6,integer_style,rows=4,cols=grep("^study_count$|^meta\\.p_count$|\\.N$",names(result6)),gridExpand=TRUE,stack=TRUE)
  addStyle(wb6,sheet6,I2_style,rows=4,cols=match("meta.I2 (%)",names(result6)),gridExpand=TRUE,stack=TRUE)
}
for(i in seq_along(wb6$worksheets_rels)) wb6$worksheets_rels[[i]]<-grep("/drawing|/vmlDrawing",wb6$worksheets_rels[[i]],invert=TRUE,value=TRUE)
wb6$Content_Types<-grep("/xl/drawings/",wb6$Content_Types,invert=TRUE,value=TRUE)
table6_dir<-tempfile("Table6_",tmpdir=normalizePath(output_dir,winslash="/"))
dir.create(table6_dir)
table6_xml_dir<-file.path(table6_dir,"xlsx")
saveWorkbook(wb6,file.path(table6_dir,"template.xlsx"),overwrite=TRUE)
unzip(file.path(table6_dir,"template.xlsx"),exdir=table6_xml_dir)
rm(wb6,result6)
invisible(gc())

for(i in seq_along(table6_files)){
  phenotype6<-table6_phenotypes[i]
  meta6<-fread(table6_files[i],na.strings=c("NA",""),encoding="UTF-8")
  if("X" %in% names(meta6)) setnames(meta6,"X","UniRef90_id")
  stopifnot(!anyNA(meta6$UniRef90_id),!anyDuplicated(meta6$UniRef90_id))
  valid6<-which(is.finite(meta6$pvalue_meta))
  stopifnot(all(meta6$pvalue_meta[valid6]>=0 & meta6$pvalue_meta[valid6]<=1))
  meta6[,c("bonferroni","fdr"):=list(NA_real_,NA_real_)]
  meta6[valid6,bonferroni:=p.adjust(pvalue_meta,"bonferroni")]
  meta6[valid6,fdr:=p.adjust(pvalue_meta,"BH")]
  input_rows6<-nrow(meta6)
  keep6<-which(is.finite(meta6$pvalue_meta) & meta6$pvalue_meta<table6_p_cutoff)
  result6<-as.data.frame(meta6[keep6])
  rm(meta6)
  annotation6<-table6_annotation[UniRef90_id %in% result6$UniRef90_id,lapply(.SD,function(x){
    x<-unique(trimws(x[!is.na(x)])); x<-x[nzchar(x)]
    if(length(x)) paste(x,collapse="; ") else NA_character_
  }),by=UniRef90_id]
  annotation_index6<-match(result6$UniRef90_id,annotation6$UniRef90_id)
  for(column6 in c("gene","ko","name","SV_location","dest")) result6[[column6]]<-annotation6[[column6]][annotation_index6]
  result6<-normalize_calvet_table(result6)
  for(j in grep(table6_numeric_pattern,names(result6))){
    result6[[j]]<-as.numeric(result6[[j]])
    result6[[j]][!is.finite(result6[[j]])]<-NA_real_
  }
  result6$FDR_less_or_equal_0.05<-ifelse(result6$fdr<=0.05,"Yes","No")
  result6$Is_Bonferroni_p_less_than_0.05<-ifelse(result6$bonferroni<0.05,"Yes","No")
  result6<-result6[order(result6$pvalue_meta,result6$UniRef90_id),,drop=FALSE]
  stopifnot(nrow(result6)+3L<=1048576L)
  table6_annotation_audit<-rbind(table6_annotation_audit,data.frame(
    phenotype=phenotype6,input_rows=input_rows6,tested=length(valid6),exported_rows=nrow(result6),
    exported_matched=sum(!is.na(annotation_index6)),exported_ko_annotated=sum(!is.na(result6$ko)),exported_name_annotated=sum(!is.na(result6$name))))
  rename6<-c(table6_rename,cilb=if(phenotype6 %in% binary_phenotypes) "meta.OR_left" else "meta.beta_left",
    ciub=if(phenotype6 %in% binary_phenotypes) "meta.OR_right" else "meta.beta_right")
  rename_index6<-match(names(result6),names(rename6))
  names(result6)[!is.na(rename_index6)]<-unname(rename6[rename_index6[!is.na(rename_index6)]])
  result6<-result6[,table6_headers[[i]],drop=FALSE]
  numeric_columns6<-which(vapply(result6,is.numeric,logical(1)))
  stopifnot(!any(vapply(result6[,numeric_columns6,drop=FALSE],function(x) any(is.infinite(x)),logical(1))))

  xml_file6<-file.path(table6_xml_dir,"xl","worksheets",paste0("sheet",i,".xml"))
  xml6<-paste0(readLines(xml_file6,warn=FALSE,encoding="UTF-8"),collapse="")
  style_nodes6<-xml2::xml_find_all(xml2::read_xml(xml6),"//*[local-name()='row'][@r='4']/*[local-name()='c']")
  style6<-setNames(xml2::xml_attr(style_nodes6,"s"),sub("[0-9]+$","",xml2::xml_attr(style_nodes6,"r")))
  columns6<-openxlsx::int2col(seq_len(ncol(result6)))
  stopifnot(all(columns6 %in% names(style6)),!anyNA(style6[columns6]))
  prefix6<-sub('<row r="4".*$',"",xml6)
  suffix6<-sub("^.*</sheetData>","</sheetData>",xml6)
  last_cell6<-paste0(tail(columns6,1),nrow(result6)+3L)
  prefix6<-sub('<dimension ref="[^"]*"/>',paste0('<dimension ref="A1:',last_cell6,'"/>'),prefix6)
  suffix6<-sub('<autoFilter ref="[^"]*"',paste0('<autoFilter ref="A3:',last_cell6,'"'),suffix6)
  connection6<-file(xml_file6,open="wb")
  writeLines(prefix6,connection6,useBytes=TRUE)
  for(first6 in if(nrow(result6)>0L) seq.int(1L,nrow(result6),by=5000L) else integer(0)){
    rows6<-seq.int(first6,min(first6+4999L,nrow(result6)))
    xml_columns6<-vector("list",ncol(result6))
    for(j in seq_len(ncol(result6))){
      values6<-result6[[j]][rows6]
      references6<-paste0('<c r="',columns6[j],rows6+3L,'" s="',style6[[columns6[j]]],'"')
      if(j %in% numeric_columns6){
        xml_columns6[[j]]<-ifelse(is.na(values6),paste0(references6,' t="inlineStr"><is><t>NA</t></is></c>'),
          paste0(references6,'><v>',format(values6,digits=17,scientific=TRUE,trim=TRUE),'</v></c>'))
      }else{
        values6<-ifelse(is.na(values6),"NA",enc2utf8(as.character(values6)))
        stopifnot(all(nchar(values6)<=32767L))
        values6<-gsub(">","&gt;",gsub("<","&lt;",gsub("&","&amp;",values6,fixed=TRUE),fixed=TRUE),fixed=TRUE)
        xml_columns6[[j]]<-paste0(references6,' t="inlineStr"><is><t xml:space="preserve">',values6,'</t></is></c>')
      }
    }
    heights6<-17*pmax(1,ceiling(nchar(ifelse(is.na(result6$name[rows6]),"",result6$name[rows6]))/46))
    writeLines(paste0('<row r="',rows6+3L,'" ht="',pmin(409,heights6),'" customHeight="1">',do.call(paste0,xml_columns6),'</row>'),connection6,useBytes=TRUE)
    if(first6==1L || tail(rows6,1)==nrow(result6) || first6 %% 100000L==1L) message("Table 6 ",phenotype6,": ",tail(rows6,1),"/",nrow(result6)," UniRef90 results written")
  }
  writeLines(suffix6,connection6,useBytes=TRUE)
  close(connection6)
  message("Table 6 ",phenotype6,": ",nrow(result6)," of ",length(valid6)," tested families exported")
  rm(result6)
  invisible(gc())
}
shared_file6<-file.path(table6_xml_dir,"xl","sharedStrings.xml")
shared_xml6<-readLines(shared_file6,warn=FALSE,encoding="UTF-8")
shared_xml6<-sub('count="[0-9]+"',paste0('count="',sum(lengths(table6_headers))+2L*length(table6_headers),'"'),shared_xml6)
writeLines(shared_xml6,shared_file6,useBytes=TRUE)
table6_staged<-file.path(table6_dir,"Supplementary Table 6.xlsx")
zip::zipr(table6_staged,files=list.files(table6_xml_dir,recursive=TRUE,all.files=TRUE,no..=TRUE),
  root=table6_xml_dir,mode="mirror",include_directories=FALSE,compression_level=1)
stopifnot(identical(openxlsx::getSheetNames(table6_staged),paste0("ST6",LETTERS[seq_along(table6_phenotypes)],"_",unname(phenotype_labels[table6_phenotypes]))))
if(!file.copy(table6_staged,file.path(output_dir,"Supplementary Table 6.xlsx"),overwrite=TRUE)) stop("Could not replace Table 6; verified workbook at ",table6_staged)
stopifnot(dirname(normalizePath(table6_dir,winslash="/"))==normalizePath(output_dir,winslash="/"))
unlink(table6_dir,recursive=TRUE)
message("Supplementary Table 6 exported")
write.csv(table6_annotation_audit,file.path(output_dir,"Supplementary_Table_6_annotation_audit.csv"),row.names=FALSE)


} else {
  message("Table 6 not updated: awaiting new server HUMAnN3 results")
}

## Supplementary Table 7: shortBRED MAG gene abundance associations with Normal P < 0.05
## Match the requested MAG annotations by New_Gene_Name without expanding rows.
table7_annotation_file<-"08.Microbial_shortbred/MAGs/HRGM2_select_gene_info_0.15.csv"
table7_annotation<-read.csv(table7_annotation_file,fileEncoding="UTF-8",check.names=FALSE,stringsAsFactors=FALSE)
table7_annotation_columns<-c("New_Gene_Name","Gene","Non-unique Gene name","Annotation","Nonzero_Ratio")
stopifnot(all(table7_annotation_columns %in% names(table7_annotation)))
table7_annotation_header<-!is.na(table7_annotation$New_Gene_Name) &
  table7_annotation$New_Gene_Name=="New_Gene_Name" & table7_annotation$Gene=="Gene" &
  table7_annotation$Annotation=="Annotation"
table7_annotation_header[is.na(table7_annotation_header)]<-FALSE
message("MAG annotation: removing ",sum(table7_annotation_header)," repeated header rows")
table7_annotation<-table7_annotation[!table7_annotation_header,table7_annotation_columns,drop=FALSE]
table7_annotation$New_Gene_Name<-trimws(table7_annotation$New_Gene_Name)
stopifnot(!anyNA(table7_annotation$New_Gene_Name),all(nzchar(table7_annotation$New_Gene_Name)))
table7_annotation<-unique(table7_annotation)
stopifnot(!anyDuplicated(table7_annotation$New_Gene_Name))
table7_annotation$Nonzero_Ratio<-as.numeric(table7_annotation$Nonzero_Ratio)
stopifnot(all(is.na(table7_annotation$Nonzero_Ratio) |
  (is.finite(table7_annotation$Nonzero_Ratio) & table7_annotation$Nonzero_Ratio>=0 & table7_annotation$Nonzero_Ratio<=1)))
table7_annotation_audit<-data.frame()
table7_files<-list.files("08.Microbial_shortbred/MAGs/anno",pattern="\\.csv$",full.names=TRUE)
table7_phenotypes<-sub("\\.csv$","",basename(table7_files))
table7_order<-unique(c(names(phenotype_labels),table7_phenotypes))
table7_phenotypes<-table7_order[table7_order %in% table7_phenotypes]
table7_files<-table7_files[match(table7_phenotypes,sub("\\.csv$","",basename(table7_files)))]
stopifnot(length(table7_files)>0,!anyDuplicated(table7_phenotypes),length(table7_files)<length(LETTERS))
table7_rename<-c(X="Gene_ID",effect_meta="meta.beta",OR_meta="meta.OR",
  pvalue_meta="meta.p",se_meta="meta.SE",cilb="meta.OR_left",ciub="meta.OR_right",
  I2="meta.I2 (%)",het_p="meta.het_p",p_count="meta.p_count",fdr="meta.fdr",bonferroni="meta.bonferroni",product="Product",
  Non.unique.Gene.name="Non_unique_Gene_name")
wb7<-createWorkbook(creator="LiuRong",title="Supplementary Table 7",
                    subject="shortBRED MAG gene abundance associations with Normal P < 0.05")
modifyBaseFont(wb7,fontName="Calibri",fontSize=11,fontColour="#000000")

for(i in seq_along(table7_files)){
  phenotype7<-table7_phenotypes[i]
  label7<-if(phenotype7 %in% names(phenotype_labels)) phenotype_labels[[phenotype7]] else phenotype7
  result7<-read.csv(table7_files[i],fileEncoding="GB18030",check.names=FALSE,stringsAsFactors=FALSE)
  result7<-normalize_calvet_table(result7)
  stopifnot(all(c("X","New_Gene_Name","pvalue_meta","I2","bonferroni","fdr") %in% names(result7)),
            !anyNA(result7$X),!anyDuplicated(result7$X))
  table7_original_ids<-result7$X
  annotation_index7<-match(trimws(result7$New_Gene_Name),table7_annotation$New_Gene_Name)
  result7$MAG_annotation_matched<-!is.na(annotation_index7)
  ## Refresh the requested source fields; retain existing Orthology and product.
  annotation_target7<-c(Gene="Gene",Non.unique.Gene.name="Non-unique Gene name",
                         Annotation="Annotation",Nonzero_Ratio="Nonzero_Ratio")
  for(annotation_column7 in names(annotation_target7)){
    if(!annotation_column7 %in% names(result7)) result7[[annotation_column7]]<-
      if(annotation_column7=="Nonzero_Ratio") NA_real_ else NA_character_
    result7[[annotation_column7]][result7$MAG_annotation_matched]<-
      table7_annotation[[annotation_target7[[annotation_column7]]]][annotation_index7[result7$MAG_annotation_matched]]
  }
  stopifnot(identical(result7$X,table7_original_ids))
  numeric_columns7<-union(which(vapply(result7,is.numeric,logical(1))),
    grep("\\.(beta|OR|SE|p|beta_left|beta_right|OR_left|OR_right|N)$|^(effect_meta|OR_meta|pvalue_meta|se_meta|cilb|ciub|I2|het_p|p_count|study_count|bonferroni|fdr)$",names(result7)))
  for(j in numeric_columns7) result7[[j]]<-as.numeric(result7[[j]])
  infinite_rows7<-which(Reduce(`|`,lapply(result7[,numeric_columns7,drop=FALSE],is.infinite)))
  if(length(infinite_rows7)>0L) result7[infinite_rows7,numeric_columns7]<-NA_real_
  stopifnot(!any(vapply(result7[,numeric_columns7,drop=FALSE],function(x) any(is.infinite(x)),logical(1))))
  result7$FDR_less_or_equal_0.05<-ifelse(is.na(result7$fdr),NA_character_,ifelse(result7$fdr<=0.05,"Yes","No"))
  result7$Is_Bonferroni_p_less_than_0.05<-ifelse(is.na(result7$bonferroni),NA_character_,ifelse(result7$bonferroni<0.05,"Yes","No"))
  result7<-result7[!is.na(result7$pvalue_meta) & result7$pvalue_meta<0.05,,drop=FALSE]
  result7<-result7[order(result7$pvalue_meta,result7$X),,drop=FALSE]
  table7_annotation_audit<-rbind(table7_annotation_audit,data.frame(
    phenotype=phenotype7,input_rows=length(table7_original_ids),
    input_matched=sum(!is.na(annotation_index7)),exported_rows=nrow(result7),
    exported_matched=sum(result7$MAG_annotation_matched),
    exported_unmatched=sum(!result7$MAG_annotation_matched)))
  if(any(!result7$MAG_annotation_matched)) warning(label7,": ",sum(!result7$MAG_annotation_matched),
    " exported genes lack the requested MAG annotation")
  ## Matching fields and model identifiers are not included in the workbook.
  result7<-result7[,setdiff(names(result7),c("New_Gene_Name","Gene_name","gene_model",
    "Nonzero_Ratio","MAG_annotation_matched")),drop=FALSE]
  cohort_columns7<-grep("\\.(beta|OR|SE|p|beta_left|beta_right|OR_left|OR_right|N)$",names(result7),value=TRUE)
  rename_index7<-match(names(result7),names(table7_rename))
  names(result7)[!is.na(rename_index7)]<-unname(table7_rename[rename_index7[!is.na(rename_index7)]])
  first_columns7<-intersect(c("Gene_ID","Annotation","meta.beta","meta.OR","meta.p","meta.I2 (%)",
    "meta.fdr","meta.bonferroni","Gene","Non_unique_Gene_name","Product","Orthology",
    "FDR_less_or_equal_0.05","Is_Bonferroni_p_less_than_0.05"),names(result7))
  result7<-result7[,c(first_columns7,setdiff(names(result7),c(first_columns7,cohort_columns7)),cohort_columns7),drop=FALSE]
  stopifnot(identical(names(result7)[1:2],c("Gene_ID","Annotation")))

  sheet7<-paste0("ST7",LETTERS[i+1L],"_",label7)
  addWorksheet(wb7,sheet7,gridLines=TRUE,zoom=90)
  writeData(wb7,sheet7,paste0("Table S7",LETTERS[i+1L],". Associations between shortBRED MAG gene abundance and ",label7,"."),startRow=1,colNames=FALSE)
  writeData(wb7,sheet7,paste0(nrow(result7)," genes with Normal P < 0.05. Flags: FDR <= 0.05; Bonferroni P < 0.05. ",
    "Inf/-Inf rows cleared before filtering: ",length(infinite_rows7),". MAG annotations matched by New_Gene_Name: ",
    table7_annotation_audit$exported_matched[i],"/",nrow(result7)," genes."),startRow=2,colNames=FALSE)
  writeData(wb7,sheet7,result7,startRow=3,headerStyle=header_style,keepNA=TRUE,na.string="NA",withFilter=nrow(result7)>0L)
  addStyle(wb7,sheet7,caption_style,rows=2,cols=1)
  setColWidths(wb7,sheet7,cols=seq_len(ncol(result7)),widths=22)
  setColWidths(wb7,sheet7,cols=match("Gene_ID",names(result7)),widths=30)
  setColWidths(wb7,sheet7,cols=match(c("Gene","Annotation","Product","Orthology"),names(result7)),widths=c(22,48,48,28))
  setColWidths(wb7,sheet7,cols=match("Non_unique_Gene_name",names(result7)),widths=30)
  setColWidths(wb7,sheet7,cols=match(c("FDR_less_or_equal_0.05","Is_Bonferroni_p_less_than_0.05"),names(result7)),widths=c(23,29))
  setRowHeights(wb7,sheet7,rows=1:3,heights=c(21,19,55))
  freezePane(wb7,sheet7,firstActiveRow=4,firstActiveCol=3)
  if(nrow(result7)>0L){
    rows7<-seq_len(nrow(result7))+3L
    addStyle(wb7,sheet7,body_style,rows=rows7,cols=seq_len(ncol(result7)),gridExpand=TRUE,stack=TRUE)
    addStyle(wb7,sheet7,identifier_style,rows=rows7,cols=match(c("Gene_ID","Gene","Non_unique_Gene_name","Annotation","Product","Orthology"),names(result7)),gridExpand=TRUE,stack=TRUE)
    addStyle(wb7,sheet7,number_style,rows=rows7,cols=which(vapply(result7,is.numeric,logical(1))),gridExpand=TRUE,stack=TRUE)
    addStyle(wb7,sheet7,pvalue_style,rows=rows7,cols=grep("^meta\\.(p|fdr|bonferroni|het_p)$|\\.p$",names(result7)),gridExpand=TRUE,stack=TRUE)
    addStyle(wb7,sheet7,integer_style,rows=rows7,cols=grep("^study_count$|^meta\\.p_count$|\\.N$",names(result7)),gridExpand=TRUE,stack=TRUE)
    addStyle(wb7,sheet7,I2_style,rows=rows7,cols=match("meta.I2 (%)",names(result7)),gridExpand=TRUE,stack=TRUE)
    setRowHeights(wb7,sheet7,rows=rows7,heights=17*pmax(1,ceiling(nchar(ifelse(is.na(result7$Annotation),"",result7$Annotation))/46),
      ceiling(nchar(ifelse(is.na(result7$Product),"",result7$Product))/46)))
  }
  message(label7,": ",nrow(result7)," genes; ",length(infinite_rows7)," Inf/-Inf rows cleared")
}

for(i in seq_along(wb7$worksheets_rels)){
  wb7$worksheets_rels[[i]]<-grep("/drawing|/vmlDrawing",wb7$worksheets_rels[[i]],invert=TRUE,value=TRUE)
}
wb7$Content_Types<-grep("/xl/drawings/",wb7$Content_Types,invert=TRUE,value=TRUE)
saveWorkbook(wb7,file.path(output_dir,"Supplementary Table 7.xlsx"),overwrite=TRUE)
write.csv(table7_annotation_audit,file.path(output_dir,"Supplementary_Table_7_MAG_annotation_audit.csv"),row.names=FALSE)

