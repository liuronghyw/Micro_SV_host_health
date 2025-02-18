### Microbiome population structure
### 2024-09-03
### LiuRong

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")

## 1 Preparation
### 1.1 Import
##install.packages("promises")

library("promises")
source("functions.R")

### 1.2 Inputs
load("01.cleanData/SV_all/all_shared_sv_dis.RData")
#load("01.cleanData/SV_all/distMat/all_msv_dist_std.RData")

all_abun_sv <- read.table("01.cleanData/mbio_all/SV_species_abun.tsv",check.names = F)
all_basic <- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

info<- read.table("01.cleanData/SV_info/Informative_species_information_final.tsv",
                      sep = "\t", header = T, stringsAsFactors = F)

## 2 Microbiome PCoA
### 2.1 SV PCoA
## based on shared sv

all_msv_dist_std_avg<-all_shared_sv_dis
all_msv_dist_std_avg_rmna<-dist_rmna(all_msv_dist_std_avg)
all_basic_sv_rmna <- all_basic[match(rownames(all_msv_dist_std_avg_rmna), rownames(all_basic)),]

all_sv_avg_dist_mds<-cmdscale(all_msv_dist_std_avg_rmna, k=5, eig = T)
all_sv_avg_dist_pcoa <- data.frame(all_sv_avg_dist_mds$points)

varExp <- (eigenvals(all_sv_avg_dist_mds)/sum(eigenvals(all_sv_avg_dist_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

#adonis_test<-adonis2(as.dist(all_msv_dist_std_avg_rmna) ~all_basic_sv_rmna$study,  permutations = 999) 
#pic_adonis <- paste0("adonis R2: ",round(adonis_test$R2,2), "; P-value: ", adonis_test$'Pr(>F)')

p_all_sv_pcoa<-ggplot(all_sv_avg_dist_pcoa,aes(X1,X2, color = all_basic_sv_rmna$inter))+
  geom_point(size = 0.8,alpha = 0.4)+ ##ggtitle("France")+
  stat_ellipse(aes(group = all_basic_sv_rmna$dataset, fill = all_basic_sv_rmna$inter, color = all_basic_sv_rmna$inter) ,type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  scale_color_manual(name=NULL, 
                    breaks = c("Europe","Americas","Asian","Africa"),
                    labels = c("Europe","Americas","Asian","Africa"),
                    values = c("#207F4C","#EA517F","#ECD452","#2983BB"))+
  #guides(col = guide_legend(nrow = 3, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size =12,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

p_all_sv_pcoa<-ggExtra::ggMarginal(p_all_sv_pcoa, type = "histogram", groupColour = F, groupFill = TRUE,
                                 xparams = list(bins = 60, alpha = 0.4,position = 'identity', color = 'white'),
                                 yparams = list(bins = 60, alpha = 0.4,position = 'identity', color = 'white'))

if(!dir.exists("05.Microbiome_population_structure")){dir.create("05.Microbiome_population_structure")}

#pdf("05.Microbiome_population_structure/all_SV_PCoA.pdf", width = 4.5, height = 5)
#print(p_all_sv_pcoa)
#dev.off()

tiff(file = "pics/F3C_all_SV_PCoA.tiff", width =3000, height =3000, res =600) 
print(p_all_sv_pcoa)
dev.off()


### 方差分析
library(gplots)
library(multcomp)

all_sv_avg<-cbind(all_basic_sv_rmna,all_sv_avg_dist_pcoa)

aggregate(all_sv_avg$X1,by=list(all_sv_avg$inter),FUN=mean)
fit <- aov(X1 ~inter,data=all_sv_avg)
summary(fit) 
TukeyHSD(fit)

plotmeans(all_sv_avg$X1 ~ all_sv_avg$inter,
          xlab = "inter",
          ylab = "PC1",
          main = "Mean Plot \nwith 95% CI")

aggregate(all_sv_avg$X2,by=list(all_sv_avg$inter),FUN=mean)
fit <- aov(X2 ~inter,data=all_sv_avg)
summary(fit) 
TukeyHSD(fit)

plotmeans(all_sv_avg$X2 ~ all_sv_avg$inter,
          xlab = "Study",
          ylab = "PC2",
          main = "Mean Plot \nwith 95% CI")


################################# 2.2 Abundance PCoA

all_abun_sv<-all_abun_sv[rownames(all_basic),]

all_abun_dist <- as.matrix(vegdist(na.omit(as.data.frame(all_abun_sv)),method = "canberra"))
all_basic_abun_rmna <- all_basic[na.omit(match(rownames(all_abun_dist), rownames(all_basic))),]

all_abun_dist_mds<-cmdscale(all_abun_dist, k=5, eig = T)
all_abun_dist_pcoa <- data.frame(all_abun_dist_mds$points)

varExp <- (eigenvals(all_abun_dist_mds)/sum(eigenvals(all_abun_dist_mds)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)

p_all_abun_pcoa<-ggplot(all_abun_dist_pcoa,aes(X1,X2, color = all_basic_abun_rmna$inter))+
  geom_point(size = 0.8,alpha = 0.5)+
  stat_ellipse(aes(group = all_basic_abun_rmna$inter,
                   fill = all_basic_abun_rmna$inter,
                   color = all_basic_abun_rmna$inter),
               type = "norm",linetype = 2,geom = "polygon",alpha = 0.05,show.legend = F)+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  scale_color_manual(name=NULL, 
                    breaks = c("Europe","Americas","Asian","Africa"),
                    labels = c("Europe","Americas","Asian","Africa"),
                    values = c("#207F4C","#EA517F","#ECD452","#2983BB"))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size =10,family ="sans",color="black"),
        legend.key = element_rect(fill = NA), 
        axis.text.x = element_text(size = 10,family ="sans",color="black"),
        axis.text.y = element_text(size = 10,family ="sans",color="black"),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))

p_all_abun_pcoa<-ggExtra::ggMarginal(p_all_abun_pcoa, type = "histogram", groupColour = F, groupFill = TRUE,
                                 xparams = list(bins = 60, alpha = 0.5,position = 'identity', color = 'white'),
                                 yparams = list(bins = 60, alpha = 0.5,position = 'identity', color = 'white'))

#pdf("05.Microbiome_population_structure/all_abun_sv_PCoA.pdf", width = 4.5, height = 5)
#print(p_all_abun_pcoa)
#dev.off()

tiff(file = "pics/FigS2_all_abun_PCoA.tiff", width =3000, height =3000, res =600) 
print(p_all_abun_pcoa)
dev.off()



### 2.3 Effect of abundance on population structure of SV

all_sv_abun_inter<-intersect(rownames(all_basic_sv_rmna), rownames(all_basic_abun_rmna))

all_inter_basic_rmna<-all_basic[match(all_sv_abun_inter, rownames(all_basic)),]
all_inter_sv_pc   <- all_sv_avg_dist_pcoa[match(all_sv_abun_inter,rownames(all_sv_avg_dist_pcoa)),]
all_inter_abun_pc <- all_sv_avg_dist_pcoa[match(all_sv_abun_inter,rownames(all_abun_dist_pcoa)),]
all_inter_sv_abun_pc<-cbind(all_inter_sv_pc, all_inter_abun_pc)
colnames(all_inter_sv_abun_pc)<-c("SV_PC1", "SV_PC2", "SV_PC3", "SV_PC4", "SV_PC5",
                                  "Abun_PC1", "Abun_PC2", "Abun_PC3", "Abun_PC4", "Abun_PC5")

## Linear model: Remove confounding effect of abundance
# SV PC1
all_sv_pc1_abun_pc_lm.res<-lm(SV_PC1~Abun_PC1+Abun_PC2+Abun_PC3+Abun_PC4+Abun_PC5, data = all_inter_sv_abun_pc)
Abun_PC1_resi<-all_sv_pc1_abun_pc_lm.res$residuals

bartlett.test(Abun_PC1_resi~all_inter_basic_rmna$inter[match(names(all_sv_pc1_abun_pc_lm.res$residuals), rownames(all_inter_basic_rmna))]
,data=all_inter_sv_abun_pc)

aov_test<-aov(Abun_PC1_resi~all_inter_basic_rmna$inter[match(names(all_sv_pc1_abun_pc_lm.res$residuals), rownames(all_inter_basic_rmna))]
,data=all_inter_sv_abun_pc)

summary(aov_test)

##wilcox.test(Abun_PC1_resi~all_inter_basic_rmna$inter[match(names(all_sv_pc1_abun_pc_lm.res$residuals), rownames(all_inter_basic_rmna))]) # p-value = 0.01357

# SV PC2
all_sv_pc2_abun_pc_lm.res<-lm(SV_PC2~Abun_PC1+Abun_PC2+Abun_PC3+Abun_PC4+Abun_PC5, data = all_inter_sv_abun_pc)
Abun_PC2_resi<-all_sv_pc2_abun_pc_lm.res$residuals

bartlett.test(Abun_PC2_resi~all_inter_basic_rmna$inter[match(names(all_sv_pc2_abun_pc_lm.res$residuals), rownames(all_inter_basic_rmna))]
,data=all_inter_sv_abun_pc)

aov_test<-aov(Abun_PC2_resi~all_inter_basic_rmna$inter[match(names(all_sv_pc2_abun_pc_lm.res$residuals), rownames(all_inter_basic_rmna))]
,data=all_inter_sv_abun_pc)

summary(aov_test)

##wilcox.test(Abun_PC2_resi~all_inter_basic_rmna$inter[match(names(all_sv_pc2_abun_pc_lm.res$residuals), rownames(all_inter_basic_rmna))]) # p-value = 1.826e-05

## PERMANOVA
all_sv_dist_abun_inter<-intersect(rownames(all_msv_dist_std_avg_rmna), rownames(all_abun_dist_pcoa))
all_inter_sv_dist<-all_msv_dist_std_avg_rmna[match(all_sv_dist_abun_inter,rownames(all_msv_dist_std_avg_rmna)),
                                             match(all_sv_dist_abun_inter,colnames(all_msv_dist_std_avg_rmna))]

all_inter_abun_pc<-all_abun_dist_pcoa[match(all_sv_dist_abun_inter, rownames(all_abun_dist_pcoa)),]
all_inter_basic<-all_basic[match(all_sv_dist_abun_inter, rownames(all_basic)),]

## 自变量中不允许出现NA
all_inter_basic<-subset(all_inter_basic,is.na(gender_code)!=T & is.na(age)!=T & is.na(inter)!=T,drop=T)
all_inter_abun_pc<-all_inter_abun_pc[rownames(all_inter_basic),]

all_inter_abun_pc_basic<-cbind(all_inter_abun_pc,all_inter_basic[,c("gender_code","age","inter")])
all_inter_sv_dist<-all_inter_sv_dist[rownames(all_inter_basic),rownames(all_inter_basic)]

all_sv_dist.adonis<-adonis(as.dist(all_inter_sv_dist)~X1+X2+X3+X4+X5, data = all_inter_abun_pc_basic)
save(all_sv_dist.adonis, file = "05.Microbiome_population_structure/all_sv_dist.adonis.RData")


############ 3 Effect of basic phenotypes
### 3.1 PERMANOVA 

all_basic_sv_p_adonis.df<-as.data.frame(matrix(NA, nrow = 8, ncol = 5))
demo_factors<-c("inter","X1","X2","X3","X5","X4","age","gender_code")

for (i in 1:8) {
  #i<-1
  all_basic_sv_p.adonis.i<-adonis(
    as.dist(all_inter_sv_dist)~., data = as.data.frame(all_inter_abun_pc_basic[,demo_factors[i]])) 
     all_basic_sv_p_adonis.df[i,c(1:3)]<-c(demo_factors[i],
     all_basic_sv_p.adonis.i$aov.tab$R2[1], 
     all_basic_sv_p.adonis.i$aov.tab$`Pr(>F)`[1])
}


all_basic_sv_p.adonis.comb<-adonis(as.dist(all_inter_sv_dist)~inter+X1+X2+X3+X5+X4+age+gender_code,data = all_inter_abun_pc_basic)
all_basic_sv_p.adonis.comb.mat<-all_basic_sv_p.adonis.comb$aov.tab[seq(1,8),] %>%
  .[match(demo_factors, rownames(.)),]


#save.image()

all_basic_sv_p_adonis.df[,4]<-all_basic_sv_p.adonis.comb.mat$R2
all_basic_sv_p_adonis.df[,5]<-all_basic_sv_p.adonis.comb.mat$`Pr(>F)`
colnames(all_basic_sv_p_adonis.df)<-c("Factor", 
                                      "Individual_R2", "Individual_p",
                                      "Combined_R2", "Combined_p")

all_basic_sv_p_adonis.df$Cumulative_R2<-cumsum(all_basic_sv_p_adonis.df$Combined_R2)
all_basic_sv_p_adonis.df$Factor[1:8]<-c( "Continent","abun_PC1","abun_PC2","abun_PC3","abun_PC5","abun_PC4","Age","Gender")
write.table(all_basic_sv_p_adonis.df, "05.Microbiome_population_structure/all_basic_sv_p_adonis.tsv",
            sep = "\t", col.names = T, row.names = F, quote = F)

save(all_basic_sv_p_adonis.df, file = "05.Microbiome_population_structure/all_basic_sv_p_adonis.df.RData")



### 3.2 Visulization
demo_factors<-c("Continent","abun_PC1", "abun_PC2","abun_PC3","abun_PC5","abun_PC4","Age","Gender")

## ProPortion
load("05.Microbiome_population_structure/all_basic_sv_p_adonis.df.RData")
knitr::kable(all_basic_sv_p_adonis.df)

all_basic_sv_p_adonis.df.long<-gather(all_basic_sv_p_adonis.df[, c(1,2,6)], 'R2', 'Value', -1)
all_basic_sv_p_adonis.df.long$Value<-as.numeric(all_basic_sv_p_adonis.df.long$Value)
all_basic_sv_p_adonis.df.long$R2 <- factor(all_basic_sv_p_adonis.df.long$R2, 
                                              levels = c("Individual_R2", "Cumulative_R2"))

p_all_basic_sv_p_adonis<-ggplot(all_basic_sv_p_adonis.df.long, aes(Factor, Value,group = R2))+
  geom_bar(aes(fill = R2),position =position_dodge(),stat="identity")+
  xlab(NULL)+
  ylab(bquote("R"^2))+
  ggtitle("All")+
  scale_x_discrete(limits = demo_factors)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(name=NULL,
                      breaks = c("Individual_R2", "Cumulative_R2"),
                      labels = c(bquote("Univariate R"^2),bquote("Cumulative R"^2) ),
                      values = mycolor2_green_blue)+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        plot.margin=unit(c(2,1,1,2),"cm"),plot.title = element_text(size = 20),
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        axis.title.y= element_text(size=20),
        axis.text.x = element_text(angle = 45, hjust = 1,size=20),
        axis.text.y = element_text(size=20),
        legend.position="top",
        legend.key = element_rect(fill = NA), 
        legend.text = element_text(size = 20),
        panel.grid.major = element_line(colour = NA,linetype = "dashed",size = 0.25),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA))

#pdf("05.Microbiome_population_structure/basic_sv_adonis_bar.pdf", width = 5.5, height = 4)
#print(p_all_basic_sv_p_adonis)
#dev.off()

tiff(file = "pics/FS1_basic_sv_adonis_bar.tiff", width =2200, height =2200, res =300) 
print(p_all_basic_sv_p_adonis)
dev.off()


