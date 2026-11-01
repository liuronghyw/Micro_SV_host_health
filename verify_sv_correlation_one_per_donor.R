options(stringsAsFactors=FALSE,warn=1)
base_out<-"02.SV_summary_statistics/SV_correlation_one_per_donor"
out<-readLines(file.path(base_out,"latest_complete_result.txt"),warn=FALSE)[1]
stopifnot(readLines(file.path(out,"status.txt"))[1]=="COMPLETE")
## Load function definitions only, never execute the production analysis again.
code<-parse("run_sv_correlation_one_per_donor.R")
for(x in code)if(is.call(x)&&identical(x[[1]],as.name("<-"))&&
                 as.character(x[[2]])%in%c("sv_fisher_vector","sv_rank_matrix","sv_group_fit","sv_exact_tile"))eval(x)
a<-read.csv(file.path(out,"SV_annotation.csv"))
m<-read.csv(file.path(out,"sample_selection_all.csv"))
retained<-m[m$retained,,drop=FALSE]
stopifnot(!anyDuplicated(retained$donor_key),all(m$retained[m$subject_id_missing]))
e<-new.env();load("01.cleanData/SV_all_QC/dsgv.RData",envir=e);load("01.cleanData/SV_all_QC/vsgv.RData",envir=e)
X<-cbind(as.matrix(e$dsgv_sub),as.matrix(e$vsgv_sub))
colnames(X)<-c(paste0("dSV:",colnames(e$dsgv_sub)),paste0("vSV:",colnames(e$vsgv_sub)))
rm(e);X<-X[match(retained$match_id,rownames(X)),a$SV,drop=FALSE]
stopifnot(identical(rownames(X),retained$match_id))
manifest<-read.csv(file.path(out,"edge_blocks_manifest.csv"))
## Pair samples are chosen throughout the full tiling, including missing data.
set.seed(20261026)
records<-list()
for(file in manifest$file[unique(round(seq(1,nrow(manifest),length.out=45)))]) {
  z<-readRDS(file);ii<-which(z$eligible)
  if(!length(ii))next
  for(k in sample(ii,min(6,length(ii)))) {
    x<-X[,z$i[k]];y<-X[,z$j[k]];obs<-is.finite(x)&is.finite(y);x<-x[obs];y<-y[obs]
    d1<-a$type[z$i[k]]=="dSV";d2<-a$type[z$j[k]]=="dSV"
    rho<-cor(x,y,method="spearman")
    if(d1&&d2)p<-fisher.test(table(factor(x,levels=0:1),factor(y,levels=0:1)),conf.int=FALSE)$p.value
    else if(d1||d2) {state<-if(d1)x else y;value<-if(d1)y else x
      p<-wilcox.test(value[state==0],value[state==1],exact=FALSE)$p.value
    } else p<-cor.test(x,y,method="spearman",exact=FALSE)$p.value
    log_error<-if(p>1e-280&&z$P[k]>1e-280)abs(log(p)-log(z$P[k]))else NA_real_
    stopifnot(length(x)==z$N[k],abs(rho-z$R[k])<1e-10,
              abs(p-z$P[k])<1e-9,is.na(log_error)||log_error<1e-6)
    records[[length(records)+1L]]<-data.frame(i=z$i[k],j=z$j[k],N=z$N[k],
      pair_type=if(d1&&d2)"dSV-dSV"else if(d1||d2)"dSV-vSV"else"vSV-vSV",
      stored_rho=z$R[k],reference_rho=rho,stored_P=z$P[k],reference_P=p,
      abs_rho_error=abs(rho-z$R[k]),abs_P_error=abs(p-z$P[k]),log_P_error=log_error)
  }
}
qachecks<-do.call(rbind,records)
write.csv(qachecks,file.path(out,"independent_pair_recalculation.csv"),row.names=FALSE)
cat("Reference pair recalculation passed:",nrow(qachecks),"pairs\n")
rm(X);gc()
## Recheck global BH independently from the saved final edge blocks.
pp<-list();qq<-list()
for(k in seq_len(nrow(manifest))) {
  z<-readRDS(manifest$file[k]);stopifnot(all(z$i<z$j),!anyDuplicated(paste(z$i,z$j)))
  pp[[k]]<-z$P[z$eligible];qq[[k]]<-z$FDR[z$eligible]
}
p<-unlist(pp,use.names=FALSE);q<-unlist(qq,use.names=FALSE);rm(pp,qq);gc()
qref<-p.adjust(p,"BH")
stopifnot(identical(q,qref),length(p)==sum(manifest$tested))
write.csv(data.frame(eligible_unique_pairs=length(p),global_BH_identical=TRUE,
                    reference_pairs=nrow(qachecks),max_abs_rho_error=max(qachecks$abs_rho_error),
                    max_abs_P_error=max(qachecks$abs_P_error)),
          file.path(out,"independent_validation_summary.csv"),row.names=FALSE)
cat("Independent global BH passed:",length(p),"tests\n")
