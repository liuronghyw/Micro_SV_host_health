###################################################################################################
####### draw pictures (Figure 1)
####### draw starked barplot for response  
####### Liurong
####### 2024-09-03

##install.packages("car")

library("ggplot2")
library("reshape2")	
library("RColorBrewer")
library("survminer")

setwd("F:/analysis_lab/pharmacomicrobiomics/3_metagenome_health/3_R")
source("functions.R")

# Basic
clin<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

table(clin$gender)
table(clin$age_category)
a<-data.frame(table(clin$country))


table(clin$inter)

temp<-table(clin$inter,clin$study_name)
write.csv(temp,"01.cleanData/phen/inter_study_name.csv")


Europe<-subset(clin,inter=="Europe",drop=T)
Europe_temp<-table(Europe$study_name)
write.csv(Europe_temp,"01.cleanData/phen/Europe_study.csv")

Asian<-subset(clin,inter=="Asian",drop=T)
Asian_temp<-table(Asian$study_name)
write.csv(Asian_temp,"01.cleanData/phen/Asian_study.csv")

Africa<-subset(clin,inter=="Africa",drop=T)
Africa_temp<-table(Africa$study_name)
write.csv(Africa_temp,"01.cleanData/phen/Africa_study.csv")

Americas<-subset(clin,inter=="Americas",drop=T)
Americas_temp<-table(Americas$study_name)
write.csv(Americas_temp,"01.cleanData/phen/Americas_study.csv")



Europe<-subset(clin,inter=="Europe",drop=T)
Europe_temp<-table(Europe$country,Europe$study_name)
write.csv(Europe_temp,"01.cleanData/phen/Europe_country_study.csv")

Asian<-subset(clin,inter=="Asian",drop=T)
Asian_temp<-table(Asian$country,Asian$study_name)
write.csv(Asian_temp,"01.cleanData/phen/Asian_country_study.csv")

Africa<-subset(clin,inter=="Africa",drop=T)
Africa_temp<-table(Africa$country,Africa$study_name)
write.csv(Africa_temp,"01.cleanData/phen/Africa_country_study.csv")

Americas<-subset(clin,inter=="Americas",drop=T)
Americas_temp<-table(Americas$country,Americas$study_name)
write.csv(Americas_temp,"01.cleanData/phen/Americas_country_study.csv")


# clin<-read.csv("00.rawData/clinical/clinical_whole.csv",header=T)
# table(clin$cancer_type)

### diseaase
### IBD T2D CRC severe obesity CAD
### ACVD adenoma cirrhosis

CRC<-subset(clin,is.na(CRC)==F,drop=T)
table(CRC$study_name)

CRC<-subset(clin,is.na(CRC)==F & study_name!="FengQ_2015",drop=T)
nrow(CRC)
table(CRC$study_name,CRC$CRC)
table(CRC$CRC)

adenoma<-subset(clin,is.na(adenoma)==F,drop=T)
table(adenoma$study_name,adenoma$adenoma)

adenoma<-subset(clin,is.na(adenoma)==F & study_name!="FengQ_2015",drop=T)
nrow(adenoma)
table(adenoma$study_name)
table(adenoma$study_name,adenoma$adenoma)
table(adenoma$adenoma)

IBD<-subset(clin,is.na(IBD)==F,drop=T)
table(IBD$study_name)
table(IBD$study_name,IBD$IBD)

IBD<-subset(clin,is.na(IBD)==F & study_name!="LiJ_2014",drop=T)
nrow(IBD)
table(IBD$study_name,IBD$IBD)
table(IBD$study_name)
table(IBD$IBD)

T2D<-subset(clin,is.na(T2D)==F,drop=T)
table(T2D$study_name,T2D$T2D)

T2D<-subset(clin,is.na(T2D)==F & study_name!="HMP_2019_t2d" & study_name!="LiJ_2014",drop=T)
nrow(T2D)
table(T2D$study_name)
table(T2D$study_name,T2D$T2D)
table(T2D$T2D)

MetS<-subset(clin,is.na(MetS)==F,drop=T)
table(MetS$study_name,MetS$MetS)

CAD<-subset(clin,is.na(CAD)==F,drop=T)
nrow(CAD)
table(CAD$CAD)

ACVD<-subset(clin,is.na(ACVD)==F,drop=T)
nrow(ACVD)
table(ACVD$ACVD)

cirrhosis<-subset(clin,is.na(cirrhosis)==F,drop=T)
nrow(cirrhosis)
table(cirrhosis$cirrhosis)



temp<-subset(clin,is.na(age)==F & disease=="healthy",drop=T)
nrow(temp)
table(temp$study_name)

## FranzosaEA_2018 with no gender 
temp<-subset(clin,is.na(age)==F & disease=="healthy" & study_name!="FengQ_2015" 
     & study_name!="HMP_2019_t2d" &  study_name!="LiJ_2014" & study_name!="HMP_2019_ibdmdb" 
& study_name!="FranzosaEA_2018",drop=T)
nrow(temp)
a<-table(temp$study_name)
write.csv(a,"age.csv")


temp<-subset(clin,is.na(BMI)==F & disease=="healthy",drop=T)
nrow(temp)
table(temp$study_name)

temp<-subset(clin,is.na(BMI)==F & disease=="healthy" & study_name!="FengQ_2015" 
     & study_name!="HMP_2019_t2d" &  study_name!="LiJ_2014" & study_name!="HMP_2019_ibdmdb",drop=T)
nrow(temp)
table(temp$study_name)

a<-table(temp$study_name)
write.csv(a,"BMI.csv")


temp<-subset(clin,is.na(hba1c)==F,drop=T)
nrow(temp)
table(temp$study_name)

a<-table(temp$study_name)
write.csv(a,"temp.csv")

temp<-subset(clin,is.na(fasting_glucose)==F,drop=T)
nrow(temp)
table(temp$study_name)

a<-table(temp$study_name)
write.csv(a,"temp.csv")

temp<-subset(clin,is.na(fasting_insulin)==F,drop=T)
nrow(temp)
table(temp$study_name)

a<-table(temp$study_name)
write.csv(a,"temp.csv")


temp<-subset(clin,is.na(cholesterol)==F,drop=T)
nrow(temp)
table(temp$study_name)

a<-table(temp$study_name)
write.csv(a,"temp.csv")


temp<-subset(clin,is.na(triglycerides)==F,drop=T)
nrow(temp)
table(temp$study_name)

a<-table(temp$study_name)
write.csv(a,"temp.csv")

temp<-subset(clin,is.na(hdl)==F,drop=T)
nrow(temp)
table(temp$study_name)

a<-table(temp$study_name)
write.csv(a,"temp.csv")

temp<-subset(clin,is.na(ldl)==F,drop=T)
nrow(temp)
table(temp$study_name)

a<-table(temp$study_name)
write.csv(a,"temp.csv")


temp<-subset(clin,is.na(mean_press)==F,drop=T)
nrow(temp)
table(temp$study_name)
a<-table(temp$study_name)
write.csv(a,"temp.csv")

temp<-subset(clin,is.na(ast)==F,drop=T)
nrow(temp)
table(temp$study_name)
a<-table(temp$study_name)
write.csv(a,"temp.csv")

temp<-subset(clin,is.na(alt)==F,drop=T)
nrow(temp)
table(temp$study_name)
a<-table(temp$study_name)
write.csv(a,"temp.csv")

temp<-subset(clin,is.na(creatinine)==F,drop=T)
nrow(temp)
table(temp$study_name)
a<-table(temp$study_name)
write.csv(a,"temp.csv")


temp<-subset(clin,is.na(total_bilirubin)==F,drop=T)
nrow(temp)
table(temp$study_name)
a<-table(temp$study_name)
write.csv(a,"temp.csv")




###################################################################
###  map (R 4.2.3)
#install.packages("rnaturalearth")
#install.packages("rnaturalearthdata")
#install.packages("rgeos")
#install.packages("ggplot2")
#install.packages("viridis")

library("ggplot2")
library("rnaturalearth")
library("rnaturalearthdata")
library(dplyr) 
library(RColorBrewer)
library("viridis")

# Basic
clin<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

country<-data.frame(table(clin$country))
#write.csv(country,"country.csv",row.names=F)
#country<-read.csv("country.csv",header=T)
colnames(country)<-c("gu_a3","sample_size")

world <- ne_countries(scale = "medium",returnclass = "sf")

world_cou<-data.frame(world$gu_a3)
colnames(world_cou)<-"gu_a3"

country_sample<-merge(country,world_cou,by.x="gu_a3",by.y="gu_a3",all.x=T,all.y=T)
##country_sample$sample_size[is.na(country_sample$sample_size)==T]<-0

world<-merge(world,country_sample,by.x="gu_a3",by.y="gu_a3")

##### world map

dt<-ggplot(data = world) +
  geom_sf(aes(geometry = geometry,fill =sample_size)) +
  #scale_fill_viridis_c(option = "B",direction = -1) +
  scale_fill_viridis_c(option = "mako", begin=0,direction = -1, alpha = 0.8,na.value = "white") + 
  #scale_fill_gradient(low="white",high="#BC3523")+
  labs(x = "Longitude", y = "Latitude") +
  #ggtitle("Sample distribution over the world",
  #         subtitle = paste0("(", 
  #                           length(unique(world$admin)),
  #                           " countries)"))+
  theme_bw()+
  theme(panel.border = element_blank(),
        legend.position="bottom",
        legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.2, "cm"),
        plot.margin = unit(c(1, 1, 1, 1), "cm"),
        plot.background = element_rect(fill = "white"),
        panel.grid = element_blank())

tiff(file = "pics/Fig2_world.tiff", width =2000, height =1300, res =300) 
dt
dev.off()



##########################################################################
##  clinical characteristic 
## Gender

all_gender_tbl <- table(clin$inter, clin$gender)
chisq.test(all_gender_tbl) 

clin_gender<-subset(clin,is.na(clin$gender)==F,drop=T)

nr<-nrow(clin_gender)
for(i in 1:nr){
  #if(is.na(clin_gender$gender[i])==T) {clin_gender$Gender[i]="unknown"}
   if(is.na(clin_gender$gender[i])==F){
    if (clin_gender$gender[i]=="female"){clin_gender$Gender[i]="female"}
    else if (clin_gender$gender[i]=="male"){clin_gender$Gender[i]="male"}
  }
 }

clin_gender$Gender<-as.factor(clin_gender$Gender)
##par(family = "sans")

p_gender<-ggplot(data=clin_gender)+
  geom_mosaic(aes(x = product(inter), fill=Gender))+
  ylab("Gender")+
  xlab("Intercontinental")+
  scale_fill_manual(values=c("#EE4863","#696969")) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        axis.ticks.length = unit(0, "cm"), 
        axis.text.x = element_text(angle = 60, hjust = 1,colour = "black",size=15),
        axis.text.y = element_text(colour = "black",size=15), 
        axis.title = element_text(colour = "black",size=18), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_gender.tiff", width =1600, height =1000, res =300) 
p_gender
dev.off()



# test color 
## height<-c(5,10)
## barplot(height,col="#EA517F")


## Age
mean(as.numeric(clin$age), na.rm = T)
se(as.numeric(clin$age))
min(as.numeric(clin$age), na.rm = T)
max(as.numeric(clin$age), na.rm = T)

clin$age<-as.numeric(clin$age)

## annova test 
##wilcox.test(clin$age~clin$study)

clin_age<-subset(clin,age!="NA",drop=T)

p_age_cont<-ggplot(clin_age,aes(x=age, color = inter, fill =inter))+
  geom_density(alpha = 0.2)+
  geom_rug(alpha = 0.5,length = unit(0.05, "npc"))+
  ylab('Density')+
  xlab('Age (year)')+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(name=NULL,
                    breaks=c("Europe","Americas","Asian","Africa"),
                    labels=c("Europe","Americas","Asian","Africa"),
                    values = c("#207F4C","#EA517F","#ECD452","#2983BB")  )+
  scale_fill_manual(name=NULL,
                    breaks=c("Europe","Americas","Asian","Africa"),
                    labels=c("Europe","Americas","Asian","Africa"),
                    values = c("#207F4C","#EA517F","#ECD452","#2983BB")  )+
  guides(col = guide_legend(nrow =1, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "bottom",
        legend.text=element_text(colour = "black",size=20),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=20),
        axis.text.y = element_text(colour = "black",size=20), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_age_num.tiff", width =2000, height =1400, res =300) 
p_age_cont
dev.off()


## Age age_category
table(clin$age_category)

## annova test 
##wilcox.test(clin$age~clin$study)
clin_age<-subset(clin,age_category!="NA" & inter!="NA",drop=T)

p_age_category<-ggplot(data=clin_age)+
  geom_mosaic(aes(x = product(inter), fill=age_category))+
  ylab("age_category")+
  xlab("inter")+
  scale_fill_manual(values=mycolor4) +
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'),
        axis.ticks.length = unit(0, "cm"), 
        axis.text.x = element_text(angle = 60, hjust = 1,colour = "black",size=15),
        axis.text.y = element_text(colour = "black",size=15), 
        axis.title = element_text(colour = "black",size=18), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_age_bins.tiff", width =2300, height =1800, res =300) 
p_age_category
dev.off()



########################################################################################
#### BMI
mean(as.numeric(clin$BMI), na.rm = T)
se(as.numeric(clin$BMI))
min(as.numeric(clin$BMI), na.rm = T)
max(as.numeric(clin$BMI), na.rm = T)

clin$BMI<-as.numeric(clin$BMI)
##wilcox.test(clin$BMI~clin$study)

clin_BMI<-subset(clin,BMI!="NA",drop=T)
table(clin_BMI$inter)

p_bmi<-ggplot(clin_BMI,aes(x=BMI, color = inter, fill =inter))+
  geom_density(alpha = 0.2)+
  geom_rug(alpha = 0.5,length = unit(0.05, "npc"))+
  ylab('Density')+
  xlab('Body mass index')+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(name=NULL,
                    breaks=c("Europe","Americas","Asian","Africa"),
                    labels=c("Europe","Americas","Asian","Africa"),
                     values =mycolor4)+
  scale_fill_manual(name=NULL,
                    breaks=c("Europe","Americas","Asian","Africa"),
                    labels=c("Europe","Americas","Asian","Africa"),
                    values = mycolor4)+
  guides(col = guide_legend(nrow =1, byrow = TRUE))+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "bottom",
        legend.text=element_text(colour = "black",size=20),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=20),
        axis.text.y = element_text(colour = "black",size=20), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=18),
        text=element_text(family="sans"))

tiff(file = "pics/F1_BMI.tiff",width =2000, height =1400, res =300) 
p_bmi
dev.off()

#mycolor6
#hist(mtcars$mpg,
#     breaks = 15,
#     col ="#20854EFF",
#)
#"#6F99ADFF","#E18727FF","#BC3C29FF","#7876B1FF","#0072B5FF",#20854EFF


##################################################################################
### disease   
### IBD T2D CRC severe obesity MetS CAD
### ACVD adenoma cirrhosis 

library(ggplot2)
library(dplyr)
library(tidyverse)

clin<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv",header=T)

clin_T2D<-subset(clin,is.na(T2D)==F,drop=T)
df_T2D<-data.frame(table(clin_T2D$T2D,clin_T2D$study_name))
df_T2D$disease<-"T2D"
df_T2D<-subset(df_T2D,Var2!="HMP_2019_t2d" & Var2!="LiJ_2014",drop=T)
sum(df_T2D$Freq)

clin_CRC<-subset(clin,is.na(CRC)==F,drop=T)
df_CRC<-data.frame(table(clin_CRC$CRC,clin_CRC$study_name))
df_CRC<-subset(df_CRC,Var2!="FengQ_2015",drop=T)
df_CRC$disease<-"CRC"
sum(df_CRC$Freq)

clin_adenoma<-subset(clin,is.na(adenoma)==F,drop=T)
df_adenoma<-data.frame(table(clin_adenoma$adenoma,clin_adenoma$study_name))
df_adenoma<-subset(df_adenoma,Var2!="FengQ_2015",drop=T)
df_adenoma$disease<-"Adenoma"
sum(df_adenoma$Freq)

clin_IBD<-subset(clin,is.na(IBD)==F,drop=T)
df_IBD<-data.frame(table(clin_IBD$IBD,clin_IBD$study_name))
df_IBD$disease<-"IBD"
df_IBD<-subset(df_IBD,Var2!="LiJ_2014",drop=T)
sum(df_IBD$Freq)

clin_MetS<-subset(clin,is.na(MetS)==F,drop=T)
df_MetS<-data.frame(table(clin_MetS$MetS,clin_MetS$study_name))
df_MetS$disease<-"MetS"
sum(df_MetS$Freq)


clin_CAD<-subset(clin,is.na(CAD)==F,drop=T)
df_CAD<-data.frame(table(clin_CAD$CAD,clin_CAD$study_name))
df_CAD$disease<-"CAD"
sum(df_CAD$Freq)

clin_ACVD<-subset(clin,is.na(ACVD)==F,drop=T)
df_ACVD<-data.frame(table(clin_ACVD$ACVD,clin_ACVD$study_name))
df_ACVD$disease<-"ACVD"
sum(df_ACVD$Freq)

clin_cirrhosis<-subset(clin,is.na(cirrhosis)==F,drop=T)
df_cirrhosis<-data.frame(table(clin_cirrhosis$cirrhosis,clin_cirrhosis$study_name))
df_cirrhosis$disease<-"Cirrhosis"
sum(df_cirrhosis$Freq)


data<-rbind(df_CRC,df_adenoma,df_IBD,df_T2D,df_MetS,df_CAD,df_ACVD,df_cirrhosis)
colnames(data)<-c("disease_state","study_name","number","disease_name")


# Set a number of 'empty bar' to add at the end of each group
empty_bar <-2
nObsType <- nlevels(as.factor(data$disease_state))
to_add <- data.frame( matrix(NA, empty_bar*nlevels(as.factor(data$disease_name))*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$disease_name <- rep(levels(as.factor(data$disease_name)), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(disease_name, study_name)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)
 
# Get the name and the y position of each label
label_data <- data %>% group_by(id, study_name) %>% summarize(tot=sum(number))

number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- data %>% 
  group_by(disease_name) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]


# Make the plot
p <- ggplot(data) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=number, fill=disease_state), stat="identity", alpha=0.8) +
  #scale_fill_viridis(discrete=TRUE) +
  scale_fill_manual(breaks=c("0","1"),labels=c("healthy","disease"),values =c("#207F4C","#98369E"))+
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 200, xend = start, yend = 200), colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 400, xend = start, yend = 400), colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 600, xend = start, yend = 600), colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 800, xend = start, yend = 800), colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 1000, xend = start, yend = 1000), colour = "grey", alpha=1, linewidth=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each lines
  ggplot2::annotate("text", x = rep(max(data$id),6), y = c(0, 200, 400, 600, 800,1000), 
            label = c("0", "200", "400", "600", "800","1000") , color="black", size=3, angle=0, fontface="bold", hjust=0.5) +
  
  ylim(-300,max(label_data$tot, na.rm=T)) +
  theme_minimal() +
  theme(
    legend.position = "left",legend.title = element_blank(), 
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=tot+10, label=study_name, hjust=hjust), color="black",alpha=1, size=3, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE)+
  geom_text(data=base_data, aes(x = title, y = -18, label=disease_name),hjust=c(1,1,1,0,0,0,0,0), angle=c(90,50,0,110,60,30,-10,-50),
                          colour = "black", alpha=1, size=3, inherit.aes = FALSE)

tiff(file = "pics/F1_disease_dis.tiff", width =2500, height =1800, res =300) 
p 
dev.off()




############################################################################
## indexes

#if (!require(devtools)) {
#    install.packages('devtools')
#}
#devtools::install_github('erocoar/gghalves')

library(ggplot2)
library(gghalves)
library(dplyr)
library(hrbrthemes)

clin<- read.table("01.cleanData/phen/Clinical_basic_overlap.tsv",header=T)

### index 
## hba1c  triglycerides hdl  ldl  cholesterol
## mean_press creatinine  alt  ast  total_bilirubin


# test color 
height<-c(5,10)
barplot(height,col="#2983BB")

## c("#207F4C","#EA517F","#ECD452","#2983BB")  

## Europe #207F4C
## American #EA517F
## Asian #ECD452
## African #2983BB


##############################################################################
### hba1c

clin_hba1c<-subset(clin,is.na(hba1c)==F,drop=T)
head(clin_hba1c) 
clin_hba1c<-clin_hba1c[,c("inter","hba1c")]

# sample size
sample_size = clin_hba1c %>% group_by(inter) %>% summarize(num=n())

table(clin_hba1c$inter)

# Plot
hba1c_pic<-clin_hba1c %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(inter, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=hba1c, fill=inter)) +
    geom_violin(width=0.9) +
    geom_boxplot(width=0.3, color="black", alpha=0.2) +
    #scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(breaks=c("Americas","Asian","Europe"),labels=c("Americas","Asian","Europe"),values =c("#2983BB","#ECD452","#207F4C"))+
    ggtitle("") +
    ylab("HbA1c (%)")+
    xlab("")+
    guides(fill = guide_legend(title = NULL))+
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "",
        legend.text=element_text(colour = "black",size=18),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_hba1c.tiff", width =1500, height =1300, res =300) 
hba1c_pic
dev.off()



##############################################################################
### fasting_insulin

clin_fasting_insulin<-subset(clin,is.na(fasting_insulin)==F,drop=T)
head(clin_fasting_insulin) 
clin_fasting_insulin<-clin_fasting_insulin[,c("inter","fasting_insulin")]

# sample size
sample_size = clin_fasting_insulin %>% group_by(inter) %>% summarize(num=n())

table(clin_fasting_insulin$inter)

# Plot
fasting_insulin_pic<-clin_fasting_insulin %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(inter, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=fasting_insulin, fill=inter)) +
    geom_violin(width=0.9) +
    geom_boxplot(width=0.3, color="black", alpha=0.2) +
    #scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(breaks=c("Asian","Europe"),labels=c("Asian","Europe"),values =c("#ECD452","#207F4C"))+
    ggtitle("") +
    ylab("Fasting insulin (%)")+
    xlab("")+
    guides(fill = guide_legend(title = NULL))+
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "",
        legend.text=element_text(colour = "black",size=18),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_fasting_insulin.tiff", width =1200, height =1300, res =300) 
fasting_insulin_pic
dev.off()



##############################################################################
### fasting_glucose

clin_fasting_glucose<-subset(clin,is.na(fasting_glucose)==F,drop=T)
head(clin_fasting_glucose) 
clin_fasting_glucose<-clin_fasting_glucose[,c("inter","fasting_glucose")]

# sample size
sample_size = clin_fasting_glucose %>% group_by(inter) %>% summarize(num=n())

table(clin_fasting_glucose$inter)

# Plot
fasting_glucose_pic<-clin_fasting_glucose %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(inter, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=fasting_glucose, fill=inter)) +
    geom_violin(width=0.9) +
    geom_boxplot(width=0.3, color="black", alpha=0.2) +
    #scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(breaks=c("Asian","Europe"),labels=c("Asian","Europe"),values =c("#ECD452","#207F4C"))+
    ggtitle("") +
    ylab("Fasting glucose (%)")+
    xlab("")+
    guides(fill = guide_legend(title = NULL))+
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "",
        legend.text=element_text(colour = "black",size=18),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_fasting_glucose.tiff", width =1000, height =1300, res =300) 
fasting_glucose_pic
dev.off()




##############################################################################
### triglycerides

clin_triglycerides<-subset(clin,is.na(triglycerides)==F,drop=T)
head(clin_triglycerides) 
clin_triglycerides<-clin_triglycerides[,c("inter","triglycerides")]

# sample size
sample_size = clin_triglycerides %>% group_by(inter) %>% summarize(num=n())

table(clin_triglycerides$inter)

# Plot
triglycerides_pic<-clin_triglycerides %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(inter, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=triglycerides, fill=inter)) +
    geom_violin(width=0.9) +
    geom_boxplot(width=0.3, color="black", alpha=0.2) +
    #scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(breaks=c("Asian","Europe"),labels=c("Asian","Europe"),values =c("#ECD452","#207F4C"))+
    ggtitle("") +
    ylab("Triglycerides (mmol/L)")+
    xlab("")+
    guides(fill = guide_legend(title = NULL))+
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "",
        legend.text=element_text(colour = "black",size=18),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_triglycerides.tiff", width =1000, height =1300, res =300) 
triglycerides_pic
dev.off()



##############################################################################
### ldl

clin_ldl<-subset(clin,is.na(ldl)==F,drop=T)
head(clin_ldl) 
clin_ldl<-clin_ldl[,c("inter","ldl")]

# sample size
sample_size = clin_ldl %>% group_by(inter) %>% summarize(num=n())

# Plot
ldl_pic<-clin_ldl %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(inter, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=ldl, fill=inter)) +
    geom_violin(width=0.9) +
    geom_boxplot(width=0.3, color="black", alpha=0.2) +
    #scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(breaks=c("Asian","Europe"),labels=c("Asian","Europe"),values =c("#ECD452","#207F4C"))+
    ggtitle("") +
    ylab("LDL (mmol/L)")+
    xlab("")+
    guides(fill = guide_legend(title = NULL))+
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "",
        legend.text=element_text(colour = "black",size=18),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_ldl.tiff", width =1000, height =1300, res =300) 
ldl_pic
dev.off()


##############################################################################
### hdl

clin_hdl<-subset(clin,is.na(hdl)==F,drop=T)
head(clin_hdl) 
clin_hdl<-clin_hdl[,c("inter","hdl")]

# sample size
sample_size = clin_hdl %>% group_by(inter) %>% summarize(num=n())

table(clin_hdl$inter)

# Plot
hdl_pic<-clin_hdl %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(inter, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=hdl, fill=inter)) +
    geom_violin(width=0.9) +
    geom_boxplot(width=0.3, color="black", alpha=0.2) +
    #scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(breaks=c("Asian","Europe"),labels=c("Asian","Europe"),values =c("#ECD452","#207F4C"))+
    ggtitle("") +
    ylab("HDL (mmol/L)")+
    xlab("")+
    guides(fill = guide_legend(title = NULL))+
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "",
        legend.text=element_text(colour = "black",size=18),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_hdl.tiff", width =1000, height =1300, res =300) 
hdl_pic
dev.off()


##############################################################################
### cholesterol

clin_cholesterol<-subset(clin,is.na(cholesterol)==F,drop=T)
head(clin_cholesterol) 
clin_cholesterol<-clin_cholesterol[,c("inter","cholesterol")]

# sample size
sample_size = clin_cholesterol %>% group_by(inter) %>% summarize(num=n())

table(clin_cholesterol$inter)

# Plot
cholesterol_pic<-clin_cholesterol %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(inter, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=cholesterol, fill=inter)) +
    geom_violin(width=0.9) +
    geom_boxplot(width=0.3, color="black", alpha=0.2) +
    #scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(breaks=c("Asian","Europe"),labels=c("Asian","Europe"),values =c("#ECD452","#207F4C"))+
    ggtitle("") +
    ylab("Total cholesterol (mmol/L)")+
    xlab("")+
    guides(fill = guide_legend(title = NULL))+
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "",
        legend.text=element_text(colour = "black",size=18),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_cholesterol.tiff", width =1000, height =1300, res =300) 
cholesterol_pic
dev.off()


##############################################################################
### mean_press

clin_mean_press<-subset(clin,is.na(mean_press)==F,drop=T)
head(clin_mean_press) 
clin_mean_press<-clin_mean_press[,c("inter","mean_press")]

# sample size
sample_size = clin_mean_press %>% group_by(inter) %>% summarize(num=n())

table(clin_mean_press$inter)

# Plot
mean_press_pic<-clin_mean_press %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(inter, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=mean_press, fill=inter)) +
    geom_violin(width=0.9) +
    geom_boxplot(width=0.3, color="black", alpha=0.2) +
    #scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(breaks=c("Asian","Europe"),labels=c("Asian","Europe"),values =c("#ECD452","#207F4C"))+
    ggtitle("") +
    ylab("Mean arterial press \n(mm/Hg)")+
    xlab("")+
    guides(fill = guide_legend(title = NULL))+
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "",
        legend.text=element_text(colour = "black",size=18),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_mean_press.tiff", width =1100, height =1300, res =300) 
mean_press_pic
dev.off()



##############################################################################
### creatinine

clin_creatinine<-subset(clin,is.na(creatinine)==F,drop=T)
head(clin_creatinine) 
clin_creatinine<-clin_creatinine[,c("inter","creatinine")]

# sample size
sample_size = clin_creatinine %>% group_by(inter) %>% summarize(num=n())

table(clin_creatinine$inter)

# Plot
creatinine_pic<-clin_creatinine %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(inter, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=creatinine, fill=inter)) +
    geom_violin(width=0.9) +
    geom_boxplot(width=0.3, color="black", alpha=0.2) +
    #scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(breaks=c("Americas","Asian","Europe"),labels=c("Americas","Asian","Europe"),values =c("#2983BB","#ECD452","#207F4C"))+
    ggtitle("") +
    ylab("Creatinine (mmol/L)")+
    xlab("")+
    guides(fill = guide_legend(title = NULL))+
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "",
        legend.text=element_text(colour = "black",size=18),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_creatinine.tiff", width =1400, height =1300, res =300) 
creatinine_pic
dev.off()



##############################################################################
### alt

clin_alt<-subset(clin,is.na(alt)==F,drop=T)
head(clin_alt) 
clin_alt<-clin_alt[,c("inter","alt")]

# sample size
sample_size = clin_alt %>% group_by(inter) %>% summarize(num=n())

table(clin_alt$inter)

# Plot
alt_pic<-clin_alt %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(inter, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=alt, fill=inter)) +
    geom_violin(width=0.9) +
    geom_boxplot(width=0.3, color="black", alpha=0.2) +
    #scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(breaks=c("Americas","Asian","Europe"),labels=c("Americas","Asian","Europe"),values =c("#EA517F","#ECD452","#207F4C"))+
    ggtitle("") +
    ylab("ALT (U/L)")+
    xlab("")+
    guides(fill = guide_legend(title = NULL))+
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "",
        legend.text=element_text(colour = "black",size=18),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_alt.tiff", width =1400, height =1300, res =300) 
alt_pic
dev.off()


##############################################################################
### ast

clin_ast<-subset(clin,is.na(ast)==F,drop=T)
head(clin_ast) 
clin_ast<-clin_ast[,c("inter","ast")]

# sample size
sample_size = clin_ast %>% group_by(inter) %>% summarize(num=n())

table(clin_ast$inter)

# Plot
ast_pic<-clin_ast %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(inter, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=ast, fill=inter)) +
    geom_violin(width=0.9) +
    geom_boxplot(width=0.3, color="black", alpha=0.2) +
    #scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(breaks=c("Americas","Asian","Europe"),labels=c("Americas","Asian","Europe"),values =c("#EA517F","#ECD452","#207F4C"))+
    ggtitle("") +
    ylab("AST (U/L)")+
    xlab("")+
    guides(fill = guide_legend(title = NULL))+
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "",
        legend.text=element_text(colour = "black",size=18),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_ast.tiff", width =1400, height =1300, res =300) 
ast_pic
dev.off()



##############################################################################
### total_bilirubin

clin_total_bilirubin<-subset(clin,is.na(total_bilirubin)==F,drop=T)
head(clin_total_bilirubin) 
clin_total_bilirubin<-clin_total_bilirubin[,c("inter","total_bilirubin")]

# sample size
sample_size = clin_total_bilirubin %>% group_by(inter) %>% summarize(num=n())

table(clin_total_bilirubin$inter)

# Plot
total_bilirubin_pic<-clin_total_bilirubin %>%
  left_join(sample_size) %>%
  mutate(myaxis = paste0(inter, "\n", "n=", num)) %>%
  ggplot( aes(x=myaxis, y=total_bilirubin, fill=inter)) +
    geom_violin(width=0.9) +
    geom_boxplot(width=0.3, color="black", alpha=0.2) +
    #scale_fill_viridis(discrete = TRUE) +
    scale_fill_manual(breaks=c("Asian"),labels=c("Asian"),values =c("#ECD452"))+
    ggtitle("") +
    ylab("Total bilirubin (umol/L)")+
    xlab("")+
    guides(fill = guide_legend(title = NULL))+
    theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        legend.position = "",
        legend.text=element_text(colour = "black",size=18),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA), 
        panel.background = element_rect(fill = NA),
        axis.text.x = element_text(colour = "black",size=18),
        axis.text.y = element_text(colour = "black",size=18), 
        axis.title = element_text(colour = "black",size=20), 
        plot.title = element_text(hjust = 0.5,size=15),
        text=element_text(family="sans"))

tiff(file = "pics/F1_total_bilirubin.tiff", width =800, height =1300, res =300) 
total_bilirubin_pic
dev.off()





#########################################################################
### output information for table 1

table_S1<-unique(read.csv("00.rawData/table_S1_check.csv",header=T))


# Basic
clin<-read.table("01.cleanData/phen/Clinical_basic_overlap.tsv")

clin_sub<-clin[,c("study_name","country","sequencing_platform")]

sample<-data.frame(table(clin_sub$study_name))
colnames(sample)[1]<-"study_name"

countray<-data.frame(table(clin_sub$study_name,clin_sub$country))
countray<-subset(countray,Freq!=0,drop=T)
colnames(countray)<-c("study_name","country","country_num")

sequencing_platform<-data.frame(table(clin_sub$study_name,clin_sub$sequencing_platform))
sequencing_platform<-subset(sequencing_platform,Freq!=0,drop=T)
colnames(sequencing_platform)<-c("study_name","platform","platform_num")

table_S1_final<-merge(table_S1,sample,by.x="study_name",by.y="study_name")
table_S1_final<-merge(table_S1_final,countray,by.x="study_name",by.y="study_name")
table_S1_final<-merge(table_S1_final,sequencing_platform,by.x="study_name",by.y="study_name")

write.csv(table_S1_final,"00.rawData/table_S1_check_final.csv")



