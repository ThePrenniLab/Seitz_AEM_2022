library(readxl)
library(ggplot2)
library(tidyverse)
library(BBmisc)
library(ggthemes)




##############################
#Import Data 
##############################
test = read_excel("/Volumes/jprenni/Projects/VL_Cover Crops/Microcosms Data/exudate_turnover_consumption_and_production.xlsx", 
                  sheet="test", col_names = TRUE)


#transform

test = as.data.frame(test)
row.names(test) = test$id

norm_test = normalize(test, method = "range", range = c(0,100), margin=2)
range(norm_test[,3])
range(norm_test[,4])
norm_test$day=as.numeric(c('1', '1', '1', '0', '0', '0', '1', '1', '1', '0', '0', '0'))
tidy_norm_test = norm_test %>%
  gather(-day,-trt, -rep, -group, key='compound', value= 'value')

#tidy
#test2 <- test %>%
  #pivot_longer(cols = -group:-day, names_to = "compound", values_to = "value") %>%
  #group_by(compound) %>%
  #arrange(desc(compound)) 


#plot
# graph shows boxplots for each treatment, faceted by compound. 
ggplot(tidy_norm_test, aes(x = group, y = value, color = trt)) +
  facet_wrap(~compound) +
  geom_boxplot() +
  labs(title = "Change in synthetic root exudates from d0 to d1") +
  theme_few()

