### SV information
### 2024-9-3
### LiuRong

library("ggplot2")
library("reshape2")	
library("RColorBrewer")

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")
source("functions.R")


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

### correlation between read count and number of SVs in a sample 
samp_dsgv_num<-apply(all_dsv, 1, myfun<-function(x){sum(!is.na(x))})
samp_vsgv_num<-apply(all_vsv, 1, myfun<-function(x){sum(!is.na(x))})

samp_sv_num<-data.frame(samp_dsgv_num+samp_vsgv_num)
samp_sv_num$match_id<-rownames(samp_sv_num)
colnames(samp_sv_num)[1]<-"sv_number"

whole<-merge(all_basic,samp_sv_num,by.x="match_id",by.y="match_id")
cor.test(whole$sv_number,whole$read_count)

#whole_sub<-subset(whole,read_count<10000000,drop=T)
#cor.test(whole_sub$sv_number,whole_sub$read_count)

p_scatter<-ggplot(data=whole,aes(x=read_count/1000000, y=sv_number),color = "white", alpha = 0.2,size = 0.5)+
  #geom_smooth(data=whole ,aes(x=read_count/1000000, y=sv_number),method = "lm", color = "white", alpha = 0.2,size = 0.5)+
  geom_point(col="#b5182b")+
  xlab("Read count (million)")+
  ylab("Number of SVs")+xlim(0,50)+
  #scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position="top",
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = "gray90",linetype = "dashed",size = 0.25),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))

tiff(file = "pics/FS2_read_count_sv_number.tiff", width =1000, height =800, res =300) 
print(p_scatter)
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
  geom_text(position = position_stack(vjust = 0.5), color = "white")+
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
        panel.grid.major = element_line(colour = "gray90",linetype = "dashed",size = 0.25),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))

tiff(file = "pics/FS1_Informative_species_number.tiff", width =4000, height =3200, res =300) 
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

p_pie<-my_pie(sv_n, 'SVs number',mycol =c("#2376B7","#EA517F") )
p_pie_nolabel<-my_pie_nolabel(sv_n, 'SVs number',mycol =c("#2376B7","#EA517F") )

#pdf("pics/Total_SVs_number.pdf",height =5, width = 5)
#print(p_pie)
#dev.off()

tiff(file = "pics/FS2_Total_SVs_number.tiff", width =1000, height =1000, res =300) 
print(p_pie)
dev.off()

tiff(file = "pics/FS2_Total_SVs_number_nolabel.tiff", width =1000, height =1000, res =300) 
print(p_pie_nolabel)
dev.off()



#### 2.1.2 SV number per species

info_svs<-aggregate(info$SVs_number, by=list(Category= info$Short_name), FUN=sum)
info_svs<-data.frame(info_svs)
colnames(info_svs)<-c("Short_name","SVs_number")

species_sgv_n_order<- info_svs$Short_name[order(info_svs$SVs_number, decreasing = T)]

#species_sgv_n_order<- info$Short_name[order(info$SVs_number, decreasing = T)]
species_sgv_n_long<-gather(info[,c(11,13,14)], "type", "number", c(2:3))

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
        axis.text.x = element_text(angle = 65, hjust = 1,size = 13,family ="Arial",color="black",face = "italic"),
        axis.text.y = element_text(size = 15,family ="Arial",color="black"),
        axis.title.y = element_text(size = 15,family ="Arial",color="black"),
        plot.title = element_text(size = 15,family ="Arial",color="black"),
        legend.position=c(0.2,0.9),
        legend.key = element_rect(fill = NA), 
        legend.text = element_text(size =15,family ="Arial",color="black"),
        panel.grid.major = element_line(colour = NA,linetype = "dashed",size = 0.2),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))


tiff(file = "pics/FS2_SVs_number_species.tiff", width =5800, height =2500, res =300) 
print(p_sv_species)
dev.off()


### 2.2 Sample size per species
### different cohorts (differernt country or wester or non-western)
infor_sample_n_order<- info$Short_name[order(info$Total_samples_number, decreasing = T)]
infor_sample_n <- info[,c(11,16,17,18,19)]

infor_sample_n.long <- gather(infor_sample_n,'Cohort', 'Sample_size', c(2:5))

p_sample_n<-ggplot(infor_sample_n.long, aes(x=Short_name, y=Sample_size,group = Cohort))+
  geom_bar(aes(fill = Cohort),position = "stack",stat="identity")+
  xlab(NULL)+
  ylab("Number of samples")+
  scale_x_discrete(limits = infor_sample_n_order)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name=NULL,
                      breaks = c("Africa_sample_number", "America_sample_number","Asian_sample_number","Europe_sample_number"),
                      labels = c("Africa", "Americas","Asian","Europe"),
                      values =c("#2983BB","#EA517F","#ECD452","#207F4C"))+
 theme(plot.subtitle = element_text(vjust =1), 
        plot.caption = element_text(vjust = 1), 
        plot.margin=unit(c(0,1,1,1),"cm"),
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        axis.text.x = element_text(angle = 65, hjust = 1,size = 11,family ="sans",color="black",face = "italic"),
        axis.text.y = element_text(size = 11,family ="sans",color="black"),
        axis.title.y = element_text(size = 20,family ="sans",color="black"),
        legend.position="top",
        legend.key = element_rect(fill = NA), 
        legend.text = element_text(size =20,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA,linetype = "dashed",size = 0.2),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))


#pdf("03.SV_summary_statistics/FS2_Samples_number_species.pdf", height = 6, width = 11) 
#print(p_sample_n)
#dev.off()

tiff(file = "pics/FS1_Samples_number_species.tiff", width =5500, height =1900, res =300) 
print(p_sample_n)
dev.off()


### 2.3 SV number factors
##Show relationships between SV numbers and genome size and sample size.
species_sample_n<-info[,c("Length","SVs_number","Total_samples_number","Deletion_SVs_number","Variable_SVs_number")]

length_sv_n_cor<-cor.test(species_sample_n$Length/1000, species_sample_n$SVs_number)
text_r<-paste("r=",round(length_sv_n_cor$estimate,digits = 3),"\np=",format(length_sv_n_cor$p.value,digits = 3),sep = "")

p_scatter_pie<-ggplot()+
  geom_smooth(data=species_sample_n ,aes(x=Length/10000, y=SVs_number),method = "lm", color = "white", alpha = 0.2,size = 0.5)+
  geom_scatterpie(data=species_sample_n ,aes(x=Length/10000, y=SVs_number,r=Total_samples_number/500),
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
        panel.grid.major = element_line(colour = "gray90",linetype = "dashed",size = 0.25),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))


#pdf("pics/Sample_n_SVs_n_scatter_pie.pdf", height = 5, width = 4)
#print(p_scatter_pie)
#dev.off()

tiff(file = "pics/FS1_Sample_n_SVs_n_scatter_pie.tiff", width =1200, height =800, res =300) 
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
mean(subset(sv_species_total_abun,inter=="Asian",drop=T)$Total_abundance,na.rm = T),
mean(subset(sv_species_total_abun,inter=="Africa",drop=T)$Total_abundance,na.rm = T))

p_sv_abun_density<-ggplot(sv_species_total_abun,aes(x=Total_abundance,,color = inter,fill =inter))+
  geom_density(alpha = 0.2)+
  geom_rug(length = unit(0.05, "npc"))+
  geom_vline(xintercept  = mean_abundance, linetype = "dashed",color =c("#207F4C","#EA517F","#ECD452","#2983BB"))+
  ylab('Density')+
  xlab('Abundance')+ggtitle("All")+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(name=NULL,
                    breaks=c("Europe","Americas","Asian","Africa"),
                    labels=c("Europe","Americas","Asian","Africa"),
                    values = c("#207F4C","#EA517F","#ECD452","#2983BB")  )+
  scale_fill_manual(name=NULL,
                    breaks=c("Europe","Americas","Asian","Africa"),
                    labels=c("Europe","Americas","Asian","Africa"),
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

tiff(file = "pics/FS1_sv_abun_density.tiff", width =1500, height =1000, res =300) 
print(p_sv_abun_density)
dev.off()


## 4 SV correlation

vsv_corr <- rcorr(as.matrix(all_vsv))
vsv_corr.r<-vsv_corr$r
vsv_corr.p<-vsv_corr$P
vsv_corr.n<-vsv_corr$n
vsv_corr.r.edge<-melt(vsv_corr.r)
vsv_corr.p.edge<-melt(vsv_corr.p)
vsv_corr.n.edge<-melt(vsv_corr.n)

vsv_corr.edge<-cbind(vsv_corr.r.edge, vsv_corr.p.edge[,-c(1:2)], vsv_corr.n.edge[,-c(1:2)])
colnames(vsv_corr.edge)<-c("SV1", "SV2", "R", "P", "N")

save(vsv_corr.edge, file = "04.SV_summary_statistics/vsv_corr.edge.RData")

dsv_corr <- rcorr(as.matrix(all_dsv))
dsv_corr.r<-dsv_corr$r
dsv_corr.p<-dsv_corr$P
dsv_corr.n<-dsv_corr$n
dsv_corr.r.edge<-melt(dsv_corr.r)
dsv_corr.p.edge<-melt(dsv_corr.p)
dsv_corr.n.edge<-melt(dsv_corr.n)

dsv_corr.edge<-cbind(dsv_corr.r.edge, dsv_corr.p.edge[,-c(1:2)], dsv_corr.n.edge[,-c(1:2)])
colnames(dsv_corr.edge)<-c("SV1", "SV2", "R", "P", "N")

save(dsv_corr.edge, file = "04.SV_summary_statistics/dsv_corr.edge.RData")


all_vsv_M_plot.list<-list()
all_dsv_M_plot.list<-list()

for (i in c(1:nrow(info))){
#  i<-1
  
  file_name<-str_replace_all(info$organism[i], "\\/", "\\.")
  spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
  
  
  vsv_i<-all_vsv[,grep(spe_name,colnames(all_vsv))]
  dsv_i<-all_dsv[,grep(spe_name,colnames(all_dsv))]
  
  if(info$organism[i]=="Bifidobacterium adolescentis"){
    p_vsv_M_i<-ggplot() + theme_void()
    
  }else{
      p_vsv_M_i<- ggcorr(vsv_i,method = c("pairwise"), alpha = 0)+
    ggtitle(info$Short_name[i])+
    theme(legend.position = "none",
          plot.title = element_text(size = 4))
    
  }

  p_dsv_M_i<- ggcorr(dsv_i,method = c("pairwise"), alpha = 0)+
    ggtitle(info$Short_name[i])+
    theme(legend.position = "none",
          plot.title = element_text(size = 4))
  
  all_vsv_M_plot.list[[i]] <- p_vsv_M_i
  all_dsv_M_plot.list[[i]] <- p_dsv_M_i
}

#pdf("pics/vsv_correlation.pdf")
#plot_grid(plotlist=all_vsv_M_plot.list)
#dev.off()

#pdf("pics/dsv_correlation.pdf")
#plot_grid(plotlist=all_dsv_M_plot.list)
#dev.off()

tiff(file = "pics/FS1_vsv_correlation.tiff", width =2000, height =2000, res =300) 
plot_grid(plotlist=all_vsv_M_plot.list)
dev.off()

tiff(file = "pics/FS1_dsv_correlation.tiff", width =2000, height =2000, res =300) 
plot_grid(plotlist=all_dsv_M_plot.list)
dev.off()




