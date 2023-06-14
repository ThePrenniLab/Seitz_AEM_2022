library(readxl)
library(ggplot2)
library(vegan)
library(grid)
library(dplyr)
library(MASS)
library(tidyr)
library(RColorBrewer)
library(stats)

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

##export the scores to put into SIMCA
data.scores = as.data.frame(scores(ord.scrs))
data.scores
write.csv(data.scores, "/Users/valerielindstrom/Documents/Documents - Valerie’s MacBook Air/Documents/Project Exudates (my copies)/16S_microcosms_2/Merged Data//data.scores.nmds.csv")
