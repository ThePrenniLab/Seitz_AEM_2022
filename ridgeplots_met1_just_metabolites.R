---
  title: "Ridge plots for metabolites"
author: "Valerie Lindstrom"
output: "ridge plots"
---
  
  library(ggplot2)
#library(gplots)
library(RColorBrewer)
library(dplyr)
library(tidyr)
#library(clusterSim)
#library(vegan)
#library(WGCNA)
library(tidyverse)
#library(limma)
#library(XLConnect)
#library(gg3D)
#library(plotly)
#library(colorRamps)
library(readxl)
#library(pheatmap)
library(BBmisc)
library(ggridges)
library(viridis)

##############################
#Notes 11.9.2020
##############################
# Dr. Anderson's youtube video ("ggplot extensions" at min 31) goes over some 
# really nice ways to clean up the ridge plots and make them look way better)
# tried experimenting with it, couldnt get it to actuallly use the viridis package
# becuase we are using "group" as our fill which is a character vector and not a numeric
# vector. numeric vectors need to be passed through the scale_fill_viridis() function. 

# I think one way to get around this would be to separate the data into 3 graphs
# one for each treatment, then I might be able to use this. 


##############################
#Import Data 
######################a########
test = read_excel("ridgeplots_metabolite_spec_abund.xlsx", sheet="only_metabolites_met1", col_names = TRUE)


##processing
test = as.data.frame(test)
row.names(test) = test$id
test = test[,-1]
#scale peak heigts 0-100
norm_test = normalize(test, method = "range", range = c(0,100), margin=2)
range(norm_test[,3])
range(norm_test[,4])
norm_test$day=as.numeric(c('0', '0', '0', '1', '1', '1', '2', '2', '2', '3', '3','3', '4', '4', '4', '5', '5', '5', '7', '7', '7', '11', '11', '11', '20', '20', '20'))
tidy_norm_test = norm_test %>%
  gather(-day,-group, key='compound', value= 'value')


#reorder the compounds to be more informative
tidy_norm_test <- tidy_norm_test %>%
  mutate(compound = fct_relevel(compound, "Succinic Acid", "Pyruvic Acid", "Fumaric Acid", "Lactic Acid", "Malonic Acid", "Phosphoric Acid", "Oxalic Acid", "3-Phenyllactic Acid", "3-Methylsalicyclic Acid"))


##plotting
test_ridge = ggplot(tidy_norm_test, aes(x=day, y=compound, height=value)) +
  geom_density_ridges_gradient(aes(fill=group), stat="identity", scale=0.75) +
  theme_ridges(grid = FALSE) 
test_ridge
ggsave("~/Documents/Documents/Project Exudates (my copies)/ridge_plots_microcosms_averages.pdf")
