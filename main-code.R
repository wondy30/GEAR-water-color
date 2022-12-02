library(tidyverse)
library(lubridate)

dataset1 <- read_csv("dataset/site9991799.csv")

#import the second dataset and combine it here first


metaData <- read_csv("dataset/locations_metadata.txt")

#Add lake name to the reluctance data
##Lakes with no names in the metadata, name was assigned with Pour_lat value
metaData <- metaData %>% 
  mutate(Lake_name = coalesce(Lake_name, as.character(Pour_lat)))

fidLakeName <- metaData %>% 
  select(FID, Lake_name) 

dataset1 <- left_join(dataset1, fidLakeName, by=c("ID"="FID"))

#aggregate to annual average
dataset1$Date <- mdy(dataset1$Date)
dataset1$Year <- year(dataset1$Date)

#Decided to calculate color DW first and then aggregate per year annual


