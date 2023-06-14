library(readxl)
library(ggplot2)
library(vegan)
library(grid)
library(dplyr)
library(MASS)
library(tidyr)
library(RColorBrewer)
library(stats)

setwd("~/Documents/Documents/Project Exudates (my copies)/16S:18S data on 2019 microcosms experiment")

##################################################
##NMDS for 16S data
##################################################
#features as columns, samples as rows
features = read_excel("~/Documents/Documents/Project Exudates (my copies)/16S:18S data on 2019 microcosms experiment/16S_feature_table_template from BM.xlsx", sheet="feature_table_abundances", col_names = TRUE)

#make a metadata frame
factors_16S = features[,1]
factors_16S$time = as.numeric(c(1,10,14,20,3,5,7,1,10,14,20,3,5,7,1,10,14,20,3,5,7,1,10,14,20,3,5,7,1,10,14,20,5,7,1,10,14,3,5,7,1,10,14,20,3,5,7,1,10,14,20,3,5,7,0,0,0))
factors_16S$treatment = c("ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","inoc","inoc","inoc")
factors_16S$category = paste(factors_16S$treatment, factors_16S$time, sep = "_")

features = features[,-1]

features = as.data.frame(sapply(features, as.numeric))
features = features*100
log_features = log(features+1)

########
##NMDS for log+1 transformed data
########
set.seed(3)
bray_dist = metaMDSdist(log_features, distance = "bray", 
                        autotransform = FALSE, noshare = 0.2, trace = 1)

NMDS_Bray <-metaMDS(log_features, distance = "bray",
                    autotransform = FALSE, maxit=800, noshare = 0.2, trace = 1)
#stress = 0.09243046
ord.scrs = as.data.frame(scores(NMDS_Bray), display="sites")

##plot with ggplot
plot_16S = ggplot(ord.scrs, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=factors_16S$treatment, size=factors_16S$time)) +
  geom_text(aes(label=factors_16S$time),color="black")+
  theme_classic()
plot_16S
ggsave("../raw_figures/16S_NMDS.pdf")

##anosim and mrpp
anosim <- anosim(bray_dist, factors_16S$treatment, permutations = 999)
anosim
#ANOSIM statistic R: 0.7503 
#Significance: 0.001 
mrpp_data = mrpp(bray_dist, factors_16S$category)
mrpp_data
#Class means and counts:
#       ct     inoc   x     
#delta 0.6666 0.4655 0.4891
#n     28     3      26    
#Chance corrected within-group agreement A: 0.2024 
#Based on observed delta 0.5751 and expected delta 0.721 
#Significance of delta: 0.001 
betadisp=betadisper(bray_dist, factors_16S$category, type="centroid")
TukeyHSD(betadisp)

##################################################
##ASV diversity
##################################################
features = read_excel("../data/16S_feature_table.xlsx", sheet="feature_table_abundances", col_names = TRUE)

factors = features[,1]
factors$time = c(1,10,14,20,3,5,7,1,10,14,20,3,5,7,1,10,14,20,3,5,7,1,10,14,20,3,5,7,1,10,14,20,5,7,1,10,14,3,5,7,1,10,14,20,3,5,7,1,10,14,20,3,5,7,0,0,0)
factors$treatment = c("ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","ct","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","inoc","inoc","inoc")
factors$category = paste(factors$treatment, factors$time, sep = "_")

features = features[,-1]

features = as.data.frame(sapply(features, as.numeric))
features = features*100

######
##diversity
######
shan = diversity(features, index="shannon")
simp = diversity(features, index="simpson")
invsimp = diversity(features, index="invsimp")
#total_rich = rowSums(features)
#convert counts into presence/absence for ASV richness
ASV_rich = features
ASV_rich = 1*(ASV_rich>0)
ASV_rich = rowSums(ASV_rich)

metrics = factors
metrics$h = shan
metrics$simp = simp
metrics$invsimp = invsimp
#metrics$total_rich = total_rich
metrics$ASV_rich = ASV_rich
metrics$time = as.character(metrics$time)
metrics = metrics %>%
  unite(group, time, treatment, remove=FALSE)
all_metrics = metrics %>%
  gather("h","invsimp","simp", "ASV_rich", key="index", value="value")
all_metrics$time=as.numeric(all_metrics$time)

write.table(all_metrics, file="ASV_diversity_metrics.txt", sep = "\t")

color=c("#b300b3","#b300b3","#b300b3","#b300b3","#b300b3","#b300b3","#b300b3","#1a001a","#ff9900","#ff9900","#ff9900","#ff9900","#ff9900","#ff9900","#ff9900","#ff9900")

diversity_plot = ggplot(all_metrics, aes(x=time, y=value, color=category)) +
  facet_grid(index~., scales = "free_y")+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.75))+
  scale_color_manual(values = color)+
  theme_bw()
diversity_plot
ggsave("../raw_figures/ASV_diversity_metrics.pdf")

inv_simp = ggplot(metrics, aes(x=time, y=invsimp, color=group))+
  geom_boxplot()+
  geom_point(position=position_dodge(width=0.75))+
  theme_bw()
inv_simp

##diversity stats
ASV_rich_aov=aov(ASV_rich~group, metrics)
summary(ASV_rich_aov)
TukeyHSD(ASV_rich_aov)

h_aov=aov(h~group, metrics)
summary(h_aov)
TukeyHSD(h_aov)

invsimp_aov=aov(invsimp~group, metrics)
summary(invsimp_aov)
TukeyHSD(invsimp_aov)

simp_aov=aov(simp~group, metrics)
summary(simp_aov)
TukeyHSD(simp_aov)



######
##ASV rank abundance
######

rank_data = read_excel("../data/16S_feature_table.xlsx", sheet="combined_rank")

rank_data = rank_data %>%
  mutate(ave_abun = ave_abun*100) %>%
  mutate(std = std*100)

rank_top50 = rank_data %>%
  filter(rank <= 100)

rank_top100 = rank_data %>%
  filter(rank <= 50)

##plot barplot with x-axis rank, y-axis average abundance
rank = ggplot(rank_top50, aes(x=rank, y=ave_abun)) +
  geom_bar(aes(fill=phylum), stat="identity", position="dodge")+
  geom_errorbar(aes(x=rank, ymin=ave_abun-std, ymax=ave_abun+std, color=phylum), position = "dodge")+
  facet_grid(Day~treatment)+
  theme_classic()
rank
