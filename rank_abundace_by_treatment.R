library(readxl)
library(ggplot2)
library(vegan)
library(grid)
library(dplyr)
library(MASS)
library(tidyr)
library(RColorBrewer)
library(stats)
library(ggthemes)

######
##ASV rank abundance
######
setwd("/Users/valerielindstrom/Documents/Documents - Valerie’s MacBook Air/Documents/Project Exudates (my copies)/16S_microcosms_2/Merged Data/Rank Abundance")

rank_data = read_excel("/Users/valerielindstrom/Documents/Documents - Valerie’s MacBook Air/Documents/Project Exudates (my copies)/16S_microcosms_2/Merged Data/Rank Abundance/rank_abundance.xlsx", 
                       sheet="E1-rank-abund")

rank_data = rank_data %>%
  mutate(ave_abun = ave_abun*100) %>%
  mutate(std = std*100)
rank_data

rank_top40 = rank_data %>%
  filter(rank <= 40)

rank_top100 = rank_data %>%
  filter(rank <= 100)

##plot barplot with x-axis rank, y-axis average abundance
##use fill=treatment to fill in bar with treatment level.
##if doing more than one treatment, fill have to include fill = "darkseagreen4" outside of the aes 
##parenthesis on line 36 in order to chagne the color. 
##decided to run two separate plots for each treatment for clarity since the % abundances were 
##so different.

rank = ggplot(rank_top40, aes(x=rank, y=ave_abun)) +
  geom_bar(aes(fill=class), stat="identity", position="dodge")+
  geom_errorbar(aes(x=rank, ymin=ave_abun, ymax=ave_abun), position = "dodge")+
  labs(x= "Rank", y= "Average % Abundance") +
  theme_classic()
rank
ggsave("/Users/valerielindstrom/Documents/Documents - Valerie’s MacBook Air/Documents/Project Exudates (my copies)/16S_microcosms_2/Merged Data/Rank Abundance/rank_abundance_E1_class.pdf")

