library(readxl)
library(tidyverse)
library(viridis)
library(gplots)

# read in data
data <- read.csv("z-scores-for-heatmap.csv", header = TRUE, check.names = FALSE)

# rename weird column
data <- data %>%
  rename(Name = "﻿")

# create matrix and set row names to be metabolites
data1 <- data.matrix(data)
row.names(data1) <- data$Name

# remove the weird name column
data1 <- data1[, colnames(data1) != "Name"]


#one types of heatmap - don't like the colors though
#heatmap(data1, Colv = NA, Rowv = NA, scale = "column", col = viridis::viridis_pal())

#heatmap from the gplots pkg, able to use viridis here
heatmap.2(data1, # data frame a matrix
          margins = c(8,15), # Adds margins below and to the right
          density.info = "none", # Remove density legend lines
          trace = "none", # Remove the blue trace lines from heatmap
          Rowv = FALSE, # Do not reorder the rows
          Colv = FALSE, #do no reorder the columns
          dendrogram = "none", # Only plot column dendrogram
          colsep=1:nrow(data1), # Add vertical grid lines
          rowsep=1:nrow(data1), # Add horizontal grid lines
          sepcolor = "black", # Color gridlines black
          key.title = "Raw Z-Score", #gets rid of key title
          col = viridis::viridis_pal()) # Make colors viridis


# More Info
https://ryjohnson09.netlify.app/post/how-to-make-a-heatmap-in-r/ 
  
