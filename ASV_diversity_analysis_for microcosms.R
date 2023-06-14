library(readxl)
library(ggplot2)
library(vegan)
library(grid)
library(dplyr)
library(MASS)
library(tidyr)
library(RColorBrewer)
library(stats)
library(forcats)
library(ggthemes)

setwd("/Users/valerielindstrom/Documents/Documents - Valerie’s MacBook Air/Documents/Project Exudates (my copies)/16S_microcosms_2/Merged Data")

##################################################
##NMDS for 16S data
##################################################
#features as columns, samples as rows
features = read_excel("/Users/valerielindstrom/Documents/Documents - Valerie’s MacBook Air/Documents/Project Exudates (my copies)/16S_microcosms_2/Merged Data/feature-table-2-final.xlsx", sheet="AT_noKC", col_names = TRUE)

#make a metadata frame
factors_16S = features[,1]
factors_16S$time = as.numeric(c(0,	0,	0,	0,	0,	0,	0,	0, 0, 1,	1,	1,	1,	1,	1,	1,	1, 1, 2, 2,	2,	2,	2,	2,	2,	2,	2,	3,	3,	3,	3,	3,	3,	3,	3,	3,	4,	4,	4,	4,	4,	4,	4,	4,	4,	5, 5, 5, 5, 5, 5, 5, 5,	7,	7,	7,	7,	7,	7,	7,	7,	7,	11,	11,	11,	11,	11,	11,	11,	11,	11,20,20,20,20,20,20,20,20,20))
factors_16S$treatment = c("E1",	"E1",	"E1",	"E2",	"E2",	"E2","ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",		"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct")
factors_16S$category = paste(factors_16S$treatment, factors_16S$time, sep = "_")

features = features[,-1]

features = as.data.frame(sapply(features, as.numeric))
features = features*100
log_features = log(features+1)

########
##NMDS for log+1 transformed data
########
set.seed(3)
bray_dist = metaMDSdist(log_features, distance = "jaccard", 
                        autotransform = FALSE, noshare = 0.2, trace = 1)

NMDS_Bray <-metaMDS(log_features, distance = "jaccard",
                    autotransform = FALSE, maxit=800, noshare = 0.2, trace = 1)
#stress = 0.09243046
ord.scrs = as.data.frame(scores(NMDS_Bray), display="sites")

##plot with ggplot
plot_16S = ggplot(ord.scrs, aes(x=NMDS1, y=NMDS2)) +
  geom_point(aes(color=factors_16S$treatment, size=factors_16S$time)) +
  geom_text(aes(label=factors_16S$time),color="black")+
  theme_classic() +
  labs(color = "Treatment", size = "Time")
plot_16S
ggsave("../16S_NMDS_noKC.pdf")

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
#Chance corrected within-group agreement A: 0.273 
#Based on observed delta 0.6185 and expected delta 0.8508
#Significance of delta: 0.001 
betadisp=betadisper(bray_dist, factors_16S$category, type="centroid")
TukeyHSD(betadisp)

##################################################
##ASV diversity
##################################################
features = read_excel("/Users/valerielindstrom/Documents/Documents - Valerie’s MacBook Air/Documents/Project Exudates (my copies)/16S_microcosms_2/Merged Data/feature-table-2-final.xlsx", sheet="AT_noKC", col_names = TRUE)

##create the dataframe for features.
factors = features[,1]
factors_16S$time = as.numeric(c(0,	0,	0,	0,	0,	0,	0,	0, 0, 1,	1,	1,	1,	1,	1,	1,	1, 1, 2, 2,	2,	2,	2,	2,	2,	2,	2,	3,	3,	3,	3,	3,	3,	3,	3,	3,	4,	4,	4,	4,	4,	4,	4,	4,	4,	5, 5, 5, 5, 5, 5, 5, 5,	7,	7,	7,	7,	7,	7,	7,	7,	7,	11,	11,	11,	11,	11,	11,	11,	11,	11,20,20,20,20,20,20,20,20,20))
factors_16S$treatment = c("E1",	"E1",	"E1",	"E2",	"E2",	"E2","ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",		"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct")
factors_16S$category = paste(factors_16S$treatment, factors_16S$time, sep = "_")

features = features[,-1]

features = as.data.frame(sapply(features, as.numeric))
features = features*100

######
##diversity
######

#create the diversity metrics you want
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

#Clean up the data using dplyr to make a nice data frame for inputting into diversity metrics
metrics <- metrics %>%
  mutate(time = c(0,	0,	0,	0,	0,	0,	0,	0, 0, 1,	1,	1,	1,	1,	1,	1,	1, 1, 2, 2,	2,	2,	2,	2,	2,	2,	2,	3,	3,	3,	3,	3,	3,	3,	3,	3,	4,	4,	4,	4,	4,	4,	4,	4,	4,	5, 5, 5, 5, 5, 5, 5, 5,	7,	7,	7,	7,	7,	7,	7,	7,	7,	11,	11,	11,	11,	11,	11,	11,	11,	11,20,20,20,20,20,20,20,20,20))%>%
  mutate(treatment = c("E1",	"E1",	"E1",	"E2",	"E2",	"E2","ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",		"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct",	"E1",	"E1",	"E1",	"E2",	"E2",	"E2",	"ct",	"ct",	"ct"))
metrics$time = as.character(metrics$time)
metrics = metrics %>%
  unite(group, time, treatment, remove=FALSE)
all_metrics = metrics %>%
  gather("h","invsimp","simp", "ASV_rich", key="index", value="value")
all_metrics$time=as.numeric(all_metrics$time)

write.table(all_metrics, file="ASV_diversity_metrics.txt", sep = "\t")

color=c("#b300b3","#b300b3","#b300b3","#b300b3","#b300b3","#b300b3","#b300b3","#1a001a","#ff9900","#ff9900","#ff9900","#ff9900","#ff9900","#ff9900","#ff9900","#ff9900")


#I changed this from her original code to be able to have a grid with all time points in separate boxes
#and the x axis equal to time
#removed the number from the x axis becuase it was too busy and the numbers at the top tell a better story.
diversity_plot = ggplot(all_metrics, aes(x=" ", y=value, color=treatment)) +
  facet_grid(index~time, scales = "free_y")+
  geom_boxplot(outlier.shape=NA)+
  theme_bw() +
  theme_few() + 
  labs(x="Time")
diversity_plot
ggsave("../ASV_diversity_metrics.pdf")

metrics$time <- as.character(metrics$time)
metrics$time <- factor(metrics$time, levels=unique(metrics$time))

#make a single shannon's plot
shannon = ggplot(metrics, aes(x=time, y= h, color = treatment))+ 
  geom_boxplot()+
  theme_bw() + 
  theme_few()
  labs(x= "Time", y="Shannons") 
shannon
ggsave("../shannons_R.pdf")

#make an inverse simpson plot
inv_simp = ggplot(metrics, aes(x=time, y=invsimp, color=treatment))+
  geom_boxplot()+
  theme_bw()+
  theme_few()
inv_simp
ggsave("../inverse_simpsons_R.pdf")


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

