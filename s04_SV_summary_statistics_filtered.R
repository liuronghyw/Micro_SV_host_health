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

load("01.cleanData/SV_all_QC/dsgv.RData")
load("01.cleanData/SV_all_QC/vsgv.RData")

all_vsv<-vsgv_sub
all_dsv<-dsgv_sub

info    <- read.table("01.cleanData/SV_info_QC/Informative_species_information_final.tsv",
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

tiff(file = "pics/F2_Total_SVs_number.tiff", width =1000, height =1000, res =300) 
print(p_pie)
dev.off()

tiff(file = "pics/F2_Total_SVs_number_nolabel.tiff", width =1000, height =1000, res =300) 
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


tiff(file = "pics/F2_SVs_number_species.tiff", width =5200, height =2500, res =300) 
print(p_sv_species)
dev.off()

